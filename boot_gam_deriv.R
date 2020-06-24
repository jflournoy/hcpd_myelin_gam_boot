library(mgcv) 
library(gratia) 
library(boot)
PLOT <- TRUE
MAX_NREG <- 400

boot_out_dir <- 'boot_out_dir'
png_out_dir <- 'png'
if(!dir.exists(png_out_dir)){
  message('Creating directory: ', file.path(png_out_dir))
  dir.create(png_out_dir)
}
if(!dir.exists(boot_out_dir)){
  message('Creating directory: ', file.path(boot_out_dir))
  dir.create(boot_out_dir)
}
if(Sys.getenv('HOME') != '/users/jflournoy' || Sys.getenv('RSTUDIO') == '1'){
  cpus_per_task <- 2
  task_id <- 0
  message('Not running as SLURM batch')
  R <- 1e4 #100
} else {
  cpus_per_task <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
  task_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  if(task_id > MAX_NREG) stop(paste0('Task number, ', task_id, ', excedes number of regions: ', MAX_NREG))
  R <- 1e4 #10,000
  message('Running as SLURM batch')
  message('CPUs: ', cpus_per_task)
  message('task: ', task_id)
  message('Iterations: ', R)
}
outfile <- file.path(boot_out_dir, paste0('bootstrap_', sprintf('%dk_', R/1e3), task_id, '.rds'))
outfile.ci <- file.path(boot_out_dir, paste0('bootstrap_ci_', sprintf('%dk_', R/1e3), task_id, '.rds'))
outfile.dt <- file.path(boot_out_dir, paste0('bootstrap_dt_', sprintf('%dk_', R/1e3), task_id, '.rds'))
outfile.plot <- file.path(png_out_dir, paste0('first_deriv_boot_', sprintf('%dk_', R/1e3), task_id, '.png'))

gam_deriv_boot <- function(d, i, form){
  library(mgcv)
  library(gratia)
  d <- d[i,]
  fit <- mgcv::gam(form, data=d, method="REML")
  deriv1 <- gratia::derivatives(fit, n = 500)
  out <- deriv1[['derivative']]
  ages <- unlist(c(deriv1[which.max(deriv1$derivative), c('data')],
                   deriv1[which.min(abs(deriv1$derivative)), c('data')]))
  return(c(out, ages))
}

hcpd_data <- readr::read_csv("hcpd_n762_myelin_cog_data.csv")

right_side <- '~ s(Age) + Sex + Scanner + numNavs_sum'

if(task_id == 0){
  left_side <- 'mean_wholebrain_T1wT2w'
} else {
  left_side <- paste0('myelin_schaefer_v', task_id)
}

mf_form <- as.formula(paste0(left_side, gsub('s\\((.*)\\)', '\\1', right_side)))
form <- as.formula(paste0(left_side, right_side))

d <- model.frame(formula = mf_form, data = hcpd_data)


RUN_BOOT <- TRUE
RUN_BOOT.CI <- TRUE
RUN_DATA_MUNGE <- TRUE
if(file.exists(outfile)) {
  message('File ', outfile, ' exists. Skipping bootstrap.')
  RUN_BOOT <- FALSE 
}
if(file.exists(outfile.ci)){
  message('File ', outfile.ci, ' exists. Skipping CI computation.')
  RUN_BOOT.CI <- FALSE 
}
if(file.exists(outfile.dt)){
  message('Files ', outfile.dt, ', and ', outfile.model, ' exist. Skipping data munging.')
  RUN_DATA_MUNGE <- FALSE 
}

if(any(c(RUN_BOOT, RUN_BOOT.CI, RUN_DATA_MUNGE))){
  if(RUN_BOOT | RUN_BOOT.CI){
    cl <- parallel::makePSOCKcluster(cpus_per_task - 1)  
  }
  
  if(RUN_BOOT){
    message('Bootstrapping...')
    system.time({boot_fit <- boot(hcpd_data, 
                                  statistic = gam_deriv_boot,
                                  R = R,
                                  sim = 'ordinary',
                                  stype = 'i',
                                  form = form,
                                  parallel = 'snow',
                                  cl = cl,
                                  ncpus = cpus_per_task - 1)})
    
    saveRDS(object = boot_fit, file = outfile)
  } else {
    message('Loading ', outfile, '...')
    boot_fit <- readRDS(file = outfile)
  }
  
  if(RUN_BOOT.CI){
    parallel::clusterExport(cl = cl, 'boot_fit')
    message('Computing CIs...')
    ci <- parallel::parLapply(cl = cl, split(1:dim(boot_fit$t)[[2]], 1:cpus_per_task), function(chunk){
      lapply(chunk, function(index){
        library(boot)
        ci_out <- lapply(c('basic', 'perc', 'bca'), function(t){
          t_name <- c('basic' = 'basic', 'perc' = 'percent', 'bca' = 'bca')[[t]]
          ci <- try(boot.ci(boot_fit, index = index, type = t))
          if(!inherits(ci, 'try-error')){
            message(t)
            r <- try(data.frame(index = index, type = t,
                                l = ci[[t_name]][4], u = ci[[t_name]][5]))
            if(inherits(r, 'try-error')) r <- data.frame()
            return(r)
          } else {
            return(data.frame())
          }
        })
        adf <- data.frame(data.table::rbindlist(ci_out))
        return(adf)
      })
    })
    saveRDS(object = ci, file = outfile.ci)
  } else {
    message('Loading ', outfile.ci, '...')
    ci <- readRDS(file = outfile.ci)
  }
  
  if(RUN_DATA_MUNGE){
    library(data.table)
    setDTthreads(cpus_per_task - 1)
    ci_dt <- data.table::rbindlist(unlist(ci, recursive = FALSE))
    
    ci_dt_age <- ci_dt[index %in% 501:502]
    ci_dt_age[, curve := factor(index, levels = 501:502, labels = c('steepest', 'plateau'))][, index := NULL]
    ci_dt_age_l <- melt(ci_dt_age, id.vars = c('type', 'curve'))
    ci_dt_age_l[, variable := paste(curve, variable, sep = '_')][, curve := NULL]
    ci_dt_age <- dcast(ci_dt_age_l, type ~ variable)
    ci_dt <- ci_dt[!index %in% 501:502]
    ci_dt <- merge(ci_dt, ci_dt_age, by = 'type')
    fit <- mgcv::gam(form, data = d, method = "REML")
    deriv1 <- gratia::derivatives(fit, n = 500)
    est_dt <- data.table(index = 1:500,
                         est = boot_fit$t0[-(501:502)],
                         boot_mean = apply(boot_fit$t[, -(501:502)], 2L, mean, na.rm = TRUE),
                         age = deriv1$data)
    ci_dt <- merge(ci_dt, est_dt, by = 'index')
    ci_dt[, bc_est := 2*est - boot_mean]
    ci_dt[, curve := dplyr::case_when(
      age < steepest_u & age > steepest_l ~ 'steepest',
      age < plateau_u & age > plateau_l ~ 'plateau',
      TRUE ~ as.character(index))] #setting as index helps with plotting

    saveRDS(file = outfile.dt, object = ci_dt)
    
    if(PLOT){
      library(ggplot2)
      library(patchwork)
      p <- lapply(unique(ci_dt$type), function(t){
        ggplot(ci_dt[type == t], aes(x = age)) +
          geom_ribbon(aes(ymin = l, ymax = u), alpha = .1) +
          geom_ribbon(aes(ymin = l, ymax = u, fill = curve), alpha = .8) +
          geom_line(aes(y = est, color = 'black')) +
          geom_line(aes(y = bc_est, color = 'red')) +
          scale_color_manual(name = 'Estimate', values = c('black', 'red'), labels = c('Model', 'Bias corrected')) + 
          scale_fill_discrete(limits = c('steepest', 'plateau')) + 
          theme_minimal() + 
          labs(title = t)
      })
      pp <- patchwork::wrap_plots(p, guides = 'collect')
      try(ggsave(filename = outfile.plot,
                 plot = pp,
                 width = 8, height = 5, units = 'in', dpi = 300))
      
    }
  } #else nothing -- if this already exists we don't need to do anything else... yet.
  
  message('Done!')
} else {
  message('Nothing to do, exiting.')
}

