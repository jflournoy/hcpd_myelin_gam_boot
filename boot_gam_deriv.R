library(mgcv) 
library(gratia) 
library(boot)

MAX_NREG <- 400

boot_out_dir <- 'boot_out_dir'
if(!dir.exists(boot_out_dir)){
  message('Creating directory: ', file.path(boot_out_dir))
  dir.create(boot_out_dir)
}
if(Sys.getenv('HOME') != '/users/jflournoy'){
  cpus_per_task <- 4
  task_id <- 0
  message('Not running on SLURM system')
  R <- 1e2 #100
} else {
  cpus_per_task <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
  task_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  if(task_id > MAX_NREG) stop(paste0('Task number, ', task_id, ', excedes number of regions: ', MAX_NREG))
  R <- 1e4 #10,000
  message('Running on SLURM system')
}

gam_deriv_boot <- function(d, i, form){
  d <- d[i,]
  fit <- gam(form, data=d, method="REML")
  deriv1 <- derivatives(fit, n = 500)
  out <- unlist(c(deriv1[which.max(deriv1$derivative), c('data', 'derivative')],
                  deriv1[which.min(abs(deriv1$derivative)), c('data', 'derivative')]))
  names(out) <- c('steep_age', 'steep_deriv', 'plat_age', 'plat_deriv')
  return(out)
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

system.time({boot_fit <- boot(hcpd_data, 
                              statistic = gam_deriv_boot,
                              R = 1000,
                              sim = 'ordinary',
                              stype = 'i',
                              form = formula(mean_wholebrain_T1wT2w ~ s(Age) + Sex + Scanner + numNavs_sum),
                              parallel = 'multicore',
                              ncpus = 4)})
