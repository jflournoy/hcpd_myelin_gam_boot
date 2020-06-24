library(mgcv) 
library(gratia) 
library(boot)

MAX_NREG <- 400

boot_out_dir <- 'boot_out_dir'
if(!dir.exists(boot_out_dir)){
  message('Creating directory: ', file.path(boot_out_dir))
  dir.create(boot_out_dir)
}
if(Sys.getenv('HOME') != '/users/jflournoy' || Sys.getenv('RSTUDIO') == '1'){
  cpus_per_task <- 4
  task_id <- 0
  message('Not running as SLURM batch')
  R <- 1e3 #100
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

cl <- parallel::makePSOCKcluster(cpus_per_task - 1)

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

parallel::clusterExport(cl = cl, 'boot_fit')

message('Computing CIs...')
ci <- parallel::parLapply(cl = cl, split(1:dim(boot_fit$t)[[2]], 1:cpus_per_task), function(chunk){
  lapply(chunk, function(index){
    library(boot)
    ci <- boot.ci(boot_fit, index = index, type = c('basic', 'perc', 'bca'))
    ci$basic
    adf <- data.frame(index = index, 
                      type = c('basic', 'perc', 'bca'), 
                      l = c(ci$basic[[4]], ci$perc[[4]], ci$bca[[4]]),
                      u = c(ci$basic[[5]], ci$perc[[5]], ci$bca[[5]]))
    return(adf)
  })
})

saveRDS(object = ci, file = outfile.ci)
message('Done!')

