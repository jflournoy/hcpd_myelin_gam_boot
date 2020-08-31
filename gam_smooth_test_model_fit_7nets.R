#+ include=FALSE
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)

#'
#' # Model Comparison
#'

#+ include=TRUE
NCFSBATCH <- TRUE
ncpus <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
if(is.na(ncpus)){
  ncpus <- 4
  NCFSBATCH <- FALSE
}

#Read in data
hcpd_data <- readr::read_csv("hcpd_n762_myelin_cog_data.csv")
#Column names of myelin measurements from 7 networks
seven_nets <- c('mean_wholebrain_T1wT2w', paste0('Yeo', 1:7, '_myelin')) #the first element is whole-brain

#+results='asis'
NADA <- lapply(1:length(seven_nets), function(i){
  right_side_s <- '~ s(Age) + Sex + Scanner + numNavs_sum'
  right_side_l <- '~ Age + Sex + Scanner + numNavs_sum'
  left_side <- seven_nets[[i]]
  
  mf_form_s <- as.formula(paste0(left_side, gsub('s\\((.*)\\)', '\\1', right_side_s)))
  form_s <- as.formula(paste0(left_side, right_side_s))
  form_l <- as.formula(paste0(left_side, right_side_l))
  
  d <- model.frame(formula = mf_form_s, data = hcpd_data)
  agerange <- range(d$Age)
  
  cat(sprintf("\n\n# Model comparison for %s\n\n", left_side))
  
  model_name_s <- ifelse(i == 1, 'fit_brm', paste0('fit_brm_', left_side))
  model_name_l <- ifelse(i == 1, 'fit_brm_l', paste0('fit_brm_l_', left_side))
  if(file.exists(paste0(model_name_s, '.rds'))){
    fit_brm_s <- readRDS(paste0(model_name_s, '.rds'))
  } else {
    stop(sprintf('Could not find previously estimated GAM: %s.rds', model_name_s))
  }
  if(file.exists(paste0(model_name_l, '.rds'))){
    fit_brm_l <- readRDS(paste0(model_name_l, '.rds'))
  } else {
    message('Sampling linear age model for comparison to GAM.\nModel forumla: ', mf_form_l)
    fit_brm_l <- brms::brm(brms::bf(form_l), 
                           data=d,
                           chains = ncpus, cores = ncpus,
                           iter = 4500, warmup = 2000,
                           control = list(adapt_delta = .999, max_treedepth = 20),
                           file = model_name_l, silent = FALSE) 
  }
  system.time({fit_brm_s <- brms::add_criterion(fit_brm_s, criterion = c("loo", "waic"))})
  system.time({fit_brm_l <- brms::add_criterion(fit_brm_l, criterion = c("loo", "waic"))})
  comparison_ls <- brms::loo_compare(fit_brm_s, fit_brm_l, criterion = 'loo', model_names = c('GAM', 'linear'))
  print(knitr::kable(comparison_ls, caption = 'Model comparison: first model is best fit.', digits = 1))
  
  return(NULL)
})
  