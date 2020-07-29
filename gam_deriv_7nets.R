library(brms)
library(data.table)
source('conditional_smooth_sample.R') #Loads functions for extracting posterior samples

get_smooth_posterior <- function(fit_brm, d, res, eps){
  c_smooths_1 <- conditional_smooth_sample(fit_brm, smooths = 's(Age)', 
                                           smooth_values = list('Age' = seq(min(d$Age), max(d$Age), length.out = res)),
                                           probs = c(0.025, 0.975),
                                           spaghetti = TRUE)
  c_smooths_2 <- conditional_smooth_sample(fit_brm, smooths = 's(Age)', 
                                           smooth_values = list('Age' = seq(min(d$Age), max(d$Age), length.out = res) + eps),
                                           probs = c(0.025, 0.975),
                                           spaghetti = TRUE)
  
  c_smooths_1_sample_data <- data.table(attr(c_smooths_1$`mu: s(Age)`, 'spaghetti'))
  c_smooths_2_sample_data <- data.table(attr(c_smooths_2$`mu: s(Age)`, 'spaghetti'))
  
  c_smooths_1_sample_data[, estimate__eps := c_smooths_2_sample_data[,estimate__]]
  c_smooths_1_sample_data[, deriv1 := (estimate__eps - estimate__) / eps]
}

NCFSBATCH <- TRUE
ncpus <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
if(is.na(ncpus)){
  ncpus <- 1
  message('Not running as sbatch job; will not complete a full run.')
  NCFSBATCH <- FALSE
}

#Read in data
hcpd_data <- readr::read_csv("hcpd_n762_myelin_cog_data.csv")
#Column names of myelin measurements from 7 networks
seven_nets <- c('mean_wholebrain_T1wT2w', paste0('Yeo', 1:7, '_myelin')) #the first element is whole-brain
#Age difference for finite differences
eps = 1e-07
#How many ages to get marginal posterior at -- doesn't need to be a lot, necessarily
res = 50

NADA <- lapply(1:length(seven_nets), function(i){
  right_side <- '~ s(Age) + Sex + Scanner + numNavs_sum'
  left_side <- seven_nets[[i]]
  
  mf_form <- as.formula(paste0(left_side, gsub('s\\((.*)\\)', '\\1', right_side)))
  form <- as.formula(paste0(left_side, right_side))
  
  d <- model.frame(formula = mf_form, data = hcpd_data)
  agerange <- range(d$Age)
  
  message(sprintf('Estimating derivative for %s', left_side))
  message(sprintf('Setting resolution for derivative across ages %.0f-%.0f to %d increments of %.1f years each.', agerange[[1]], agerange[[2]], res, diff(agerange)/res))
  
  model_name <- ifelse(i == 1, 'fit_brm', paste0('fit_brm_', left_side))
  if(!file.exists(paste0(model_name, '.rds')) & !NCFSBATCH){
    stop('Not running as sbatch job, will not run Bayesian estimation')
  }
  fit_brm <- brms::brm(brms::bf(form), 
                       data=d,
                       chains = 4, cores = 4,
                       iter = 4500, warmup = 2000,
                       control = list(adapt_delta = .999, max_treedepth = 20),
                       file = model_name, silent = FALSE)
  
  posterior_samples <- get_smooth_posterior(fit_brm = fit_brm, d = d, res = res, eps = eps)
  
  smooth_post_summary <- posterior_samples[, list(est = median(estimate__),
                                                  'lower' = quantile(estimate__, probs = .025), 
                                                  'upper' = quantile(estimate__, probs = .975)),
                                           by = 'Age']
  
  deriv1_post_summary <- posterior_samples[, list('deriv1' = median(deriv1), 
                                                  'lower' = quantile(deriv1, probs = .025), 
                                                  'upper' = quantile(deriv1, probs = .975)),
                                           by = 'Age']
  
  deriv1_age_at_max_post <- posterior_samples[, list('max_age' = Age[which(deriv1 == max(deriv1))]),
                                              by = 'sample__']
  age_at_max_med_deriv <- deriv1_post_summary[deriv1 == max(deriv1), Age]
  deriv1_diff_from_max_post <- posterior_samples[, diff_from_max := deriv1 - deriv1[Age == age_at_max_med_deriv],
                                                 by = 'sample__']
  deriv1_age_at_max_post_summary <- deriv1_age_at_max_post[, list('est' = median(max_age),
                                                                  'lower' = quantile(max_age, probs = .025), 
                                                                  'upper' = quantile(max_age, probs = .975))]
  
  deriv1_diff_from_max_post_summary <- deriv1_diff_from_max_post[, list('diff_from_max' = median(diff_from_max), 
                                                                        'lower' = quantile(diff_from_max, probs = .025), 
                                                                        'upper' = quantile(diff_from_max, probs = .975)),
                                                                 by = 'Age']
  deriv1_diff_from_max_post_summary[, compared_to_max := dplyr::case_when(upper < 0 ~ 'less_steep',
                                                                          lower > 0 ~ 'steeper',
                                                                          TRUE ~ 'as_steep')]
  
  deriv1_post_summary <- deriv1_post_summary[deriv1_diff_from_max_post_summary[, .(Age, compared_to_max)], on = 'Age']
  
  age_summary_fn <- paste0(model_name, '_age_summary.csv')
  deriv_summary_fn <- paste0(model_name, '_deriv_summary.csv')
  readr::write_csv(x = deriv1_age_at_max_post_summary, path = age_summary_fn)
  readr::write_csv(x = deriv1_diff_from_max_post_summary, path = deriv_summary_fn)
})