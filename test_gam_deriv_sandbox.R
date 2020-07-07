#'---
#'title: "Finding region of steepest slope"
#'author: "John Flournoy"
#'date: "`r Sys.Date()`"
#'output:
#'  html_document:
#'    toc: true
#'    toc_float: true
#'    code_folding: hide
#'---
#'

#+ setup, include = FALSE
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)
cpal <- paste0('#', c('003E51', '007DBA', 'FFFFFF', 'D6D2C4', 'B78B20'))

#+ load packages
library(mgcv) 
library(gratia) 
library(boot)
library(dplyr)
library(ggplot2)

#' Loading and preparing data, including fitting GAM models with smooth terms
#' using both maximum-likelihood and Bayesian methods.
#'

#+ load data
hcpd_data <- readr::read_csv("hcpd_n762_myelin_cog_data.csv")
eps = 1e-07
res = 50

right_side <- '~ s(Age) + Sex + Scanner + numNavs_sum'
left_side <- 'mean_wholebrain_T1wT2w'

mf_form <- as.formula(paste0(left_side, gsub('s\\((.*)\\)', '\\1', right_side)))
form <- as.formula(paste0(left_side, right_side))

d <- model.frame(formula = mf_form, data = hcpd_data)
agerange <- range(d$Age)

#+results='asis'
cat(sprintf('Setting resolution for derivative across ages %.0f-%.0f to %d increments of %.1f years each.', agerange[[1]], agerange[[2]], res, diff(agerange)/res))

fit <- mgcv::gam(form, data=d, method="REML")
deriv1 <- gratia::derivatives(fit, n = res, eps = eps, level = .95, interval = 'simultaneous', unconditional = TRUE, n_sim = 10000)
deriv1aa <- gratia::derivatives(fit, n = res, eps = eps, level = .95, interval = 'simultaneous', unconditional = FALSE, n_sim = 10000)
deriv1a <- gratia::derivatives(fit, n = res, eps = eps, level = .95, interval = 'confidence', unconditional = TRUE)
deriv1b <- gratia::derivatives(fit, n = res, eps = eps, level = .95, interval = 'confidence', unconditional = FALSE)

library(brms)
library(data.table)

fit_brm <- brms::brm(brms::bf(form), 
                     data=d,
                     chains = 4, cores = 4,
                     iter = 4500, warmup = 2000,
                     control = list(adapt_delta = .99, max_treedepth = 20),
                     file = 'fit_brm')


source('conditional_smooth_sample.R')
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

smooth_post_summary <- c_smooths_1_sample_data[, list(est = median(estimate__),
                                                      'lower' = quantile(estimate__, probs = .025), 
                                                      'upper' = quantile(estimate__, probs = .975)),
                                               by = 'Age']

deriv1_post_summary <- c_smooths_1_sample_data[, list('deriv1' = median(deriv1), 
                                                      'lower' = quantile(deriv1, probs = .025), 
                                                      'upper' = quantile(deriv1, probs = .975)),
                                               by = 'Age']

deriv1_age_at_max_post <- c_smooths_1_sample_data[, list('max_age' = Age[which(deriv1 == max(deriv1))]),
                                                  by = 'sample__']
age_at_max_med_deriv <- deriv1_post_summary[deriv1 == max(deriv1), Age]
deriv1_diff_from_max_post <- c_smooths_1_sample_data[, diff_from_max := deriv1 - deriv1[Age == age_at_max_med_deriv],
                                                     by = 'sample__']
deriv1_diff_from_max_post_summary <- deriv1_diff_from_max_post[, list('diff_from_max' = median(diff_from_max), 
                                                                      'lower' = quantile(diff_from_max, probs = .025), 
                                                                      'upper' = quantile(diff_from_max, probs = .975)),
                                                               by = 'Age']
deriv1_diff_from_max_post_summary[, compared_to_max := case_when(upper < 0 ~ 'less_steep',
                                                                 lower > 0 ~ 'steeper',
                                                                 TRUE ~ 'as_steep')]

deriv1_post_summary <- deriv1_post_summary[deriv1_diff_from_max_post_summary[, .(Age, compared_to_max)], on = 'Age']

#'
#' # Summary of ML and Bayesian fits {.tabset}
#'
#' ## ML
#' 

#+ summaryml
summary(fit)

#'
#' ## Bayesian
#' 

#+ summarybayes
summary(fit_brm)

#' # Sources of uncertainty
#'
#' I'll start off by considering the important issue of simultaneous confidence
#' intervals described in this [blog
#' post](https://fromthebottomoftheheap.net/2014/06/16/simultaneous-confidence-intervals-for-derivatives/)
#' by the author of the `gratia` package. Essentially, if we are to compare,
#' point by point, the confidence interval, we need to correct for multiple
#' comparisons. See below for how this changes the width of the interval around
#' the first derivative plot. Also note how making the interval unconditional on
#' the smooth terms slightly widens it. This is to add in the uncertainty about
#' the exact shape of the smooth.
#'
#' The bayesian method does not face this particular multiple comparison
#' problem. I use the posterior of the derivative to compute the posterior for
#' the difference between the derivative at a given age to the derivative at the
#' age where the median posterior derivative is at its maximum (i.e., at `r sprintf('%0.1f', age_at_max_med_deriv)` 
#' years). I then compute the 95% credible interval for that statistic at each age.
#' If at each age 95% of the posterior falls in a certain range, then we can also be sure that 
#' 95% of the posterior across all ages falls in this range, so we are not inflating the probability
#' that we make a mistake in deciding that this range includes the true region of steepest slope.
#' 
#' In the plots below, you will notice that the more sources of uncertainty we account for 
#' in the ML models (that is, the non-Bayesian models), the wider the confidence regions
#' are. The widest interval is when we account for the multiple tests (i.e., a test at each age
#' we're interested in examining the derivative), and when we account for uncertainty in the 
#' exact shape of the spline. "Simul" refers to the use of simultaneous confidence intervals 
#' that account for multiple testing, whereas "Conf" refers to traditional confidence intervals.
#' "Uncon" refers to intervals that are unconditional on the choice of shape of the spline and "Cond" 
#' refers to intervals that are conditional on the exact trajectory as estimated.
#' 
#' The derivative of the Bayesian model falls somewhere in the middle, and has a slightly different shape.
#' 
#' ## Plots {.tabset}
#' 
#' ### ML Intervals
#' 
ggplot(deriv1, aes(x = data, y = derivative)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper, color = 'Simul/Uncond', fill = 'Simul/Uncond'), alpha = .25) + 
  geom_ribbon(data = deriv1aa, aes(ymin = lower, ymax = upper, color = 'Simul/Cond', fill = 'Simul/Cond'), alpha = .25) + 
  geom_ribbon(data = deriv1a, aes(ymin = lower, ymax = upper, color = 'Conf/Uncond', fill = 'Conf/Uncond'), alpha = .25) + 
  geom_ribbon(data = deriv1b, aes(ymin = lower, ymax = upper, color = 'Conf/Cond', fill = 'Conf/Cond'), alpha = .25) + 
  scale_color_manual(breaks = c('Simul/Uncond', 'Simul/Cond', 'Conf/Uncond', 'Conf/Cond', 'Bayes'),
                     values = c(cpal[c(1:2, 4:5)], '#000000'),
                     name = 'Type of interval') + 
  scale_fill_manual(breaks = c('Simul/Uncond', 'Simul/Cond', 'Conf/Uncond', 'Conf/Cond', 'Bayes'),
                     values = c(cpal[c(1:2, 4:5)], '#000000'),
                     name = 'Type of interval') + 
  geom_line() + 
  theme_minimal()
#'
#' ### Bayesian interval (w/ML)
#'
ggplot(deriv1, aes(x = data, y = derivative)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper, color = 'Simul/Uncond', fill = 'Simul/Uncond'), alpha = .1) + 
  geom_ribbon(data = deriv1aa, aes(ymin = lower, ymax = upper, color = 'Simul/Cond', fill = 'Simul/Cond'), alpha = .1) + 
  geom_ribbon(data = deriv1a, aes(ymin = lower, ymax = upper, color = 'Conf/Uncond', fill = 'Conf/Uncond'), alpha = .1) + 
  geom_ribbon(data = deriv1b, aes(ymin = lower, ymax = upper, color = 'Conf/Cond', fill = 'Conf/Cond'), alpha = .1) + 
  geom_line(alpha = .5) + 
  geom_ribbon(data = deriv1_post_summary, aes(x = Age, y = deriv1, ymin = lower, ymax = upper, color = 'Bayes', fill = 'Bayes'), 
              alpha = .7) + 
  geom_line(data = deriv1_post_summary, aes (x = Age, y = deriv1)) + 
  scale_color_manual(breaks = c('Simul/Uncond', 'Simul/Cond', 'Conf/Uncond', 'Conf/Cond', 'Bayes'),
                     values = c(cpal[c(1:2, 4:5)], '#aaaaaa'),
                     name = 'Type of interval') + 
  scale_fill_manual(breaks = c('Simul/Uncond', 'Simul/Cond', 'Conf/Uncond', 'Conf/Cond', 'Bayes'),
                    values = c(cpal[c(1:2, 4:5)], '#ffffff'),
                    name = 'Type of interval') + 
  theme_minimal()


#'
#' # Testing for regions of steepest slope
#'

fit_eval_ignore <- gratia:::confint.gam(fit, parm = 's(Age)', n = 100, unconditional = FALSE, type = 'confidence')
fit_eval <- gratia:::confint.gam(fit, parm = 's(Age)', n = 100, unconditional = TRUE, type = 'simultaneous', nsim = 10000)
#' In order to test where the GAM smooth derivative becomes different than the peak, 
#' we can explore a few different options:
#' 
#' 1. Where does the the 83% CI for the maximum derivative no longer overlap the 83% CI for the derivative at other points?
#' 2. Where is the maximum derivative different from the derivative at other points using pooled standard errors?
#' 3. Using the posterior of the estimated smooth, compute the posterior of other quantities of interest, including
#'     1. The posterior for the derivative of the smooth
#'     2. The age where the median derivative is at its maximum
#'     2. The age where the derivative is maximized (across all posterior samples)
#'     3. The difference between the derivative, and the deriviative at the age computed in (2)
#'
#' For the ML smooths, we will look at what happens when we take into account the two sources of 
#' uncertainty discussed above, and when we ignore them.
#'
#' ## Difference in 83% CIs {.tabset}
#'

deriv1_83 <- gratia::derivatives(fit, n = 100, level = .83, interval = 'simultaneous', unconditional = TRUE, n_sim = 10000)
i_max <- which(deriv1_83$derivative == max(deriv1_83$derivative))
max_lower <- deriv1_83$lower[i_max]

deriv1_83 <- deriv1_83 %>%
  mutate(compared_to_max = case_when(
    upper < max_lower ~ 'less_steep',
    TRUE ~ 'as_steep'
  ))


#'
#' ### Derivative
#'
ggplot(deriv1_83, aes(x = data, y = derivative)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .5, fill = cpal[4]) + 
  geom_line() + 
  geom_line(aes(color = compared_to_max, group = 1), size = 2, alpha = .8) + 
  scale_color_manual(breaks = c('less_steep', 'as_steep'), values = cpal[c(5,2)], 
                     labels = c('Less steep', 'As steep'),
                     name = 'Region, as compared \nto maximum derivative, \nis...') + 
  theme_minimal()

#'
#' ### Smooth estimate
#'
ggplot(cbind(fit_eval, select(deriv1_83, compared_to_max)), aes(x = Age, y = est)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .5, fill = cpal[4]) +  
  geom_line() + 
  geom_line(aes(color = compared_to_max, group = 1), size = 2, alpha = .8) + 
  scale_color_manual(breaks = c('less_steep', 'as_steep'), values = cpal[c(5,2)], 
                     labels = c('Less steep', 'As steep'),
                     name = 'Region, as compared \nto maximum derivative, \nis...') + 
  theme_minimal()

#'
#' ### Derivative (ignoring uncertainty)
#'

deriv1b_83 <- gratia::derivatives(fit, n = 100, level = .83, interval = 'confidence', unconditional = FALSE)
i_max <- which(deriv1b_83$derivative == max(deriv1b_83$derivative))
max_lower <- deriv1b_83$lower[i_max]

deriv1b_83 <- deriv1b_83 %>%
  mutate(compared_to_max = case_when(
    upper < max_lower ~ 'less_steep',
    TRUE ~ 'as_steep'
  ))

ggplot(deriv1b_83, aes(x = data, y = derivative)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .5, fill = cpal[4]) +  
  geom_line() + 
  geom_line(aes(color = compared_to_max, group = 1), size = 2, alpha = .8) + 
  scale_color_manual(breaks = c('less_steep', 'as_steep'), values = cpal[c(5,2)], 
                     labels = c('Less steep', 'As steep'),
                     name = 'Region, as compared \nto maximum derivative, \nis...') +
  theme_minimal()

#'
#' ### Smooth estimate (ignoring uncertainty)
#'

ggplot(cbind(fit_eval_ignore, select(deriv1b_83, compared_to_max)), aes(x = Age, y = est)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .5, fill = cpal[4]) +  
  geom_line() + 
  geom_line(aes(color = compared_to_max, group = 1), size = 2, alpha = .8) + 
  scale_color_manual(breaks = c('less_steep', 'as_steep'), values = cpal[c(5,2)], 
                     labels = c('Less steep', 'As steep'),
                     name = 'Region, as compared \nto maximum derivative, \nis...') +
  theme_minimal()

#'
#' ##  Difference in point estimates using pooled standard errors {.tabset}
#'

deriv1_se <- gratia::derivatives(fit, n = 100, level = .95, interval = 'simultaneous', unconditional = TRUE, n_sim = 10000)
i_max <- which(deriv1_se$derivative == max(deriv1_se$derivative))
max_est_se <- deriv1_se[i_max, c('derivative', 'se')]

deriv1_se <- deriv1_se %>%
  mutate(diff = derivative - max_est_se$derivative,
         se_diff = sqrt(se^2 + max_est_se$se^2),
         Z = diff/se_diff,
         compared_to_max = case_when(
           Z < -crit ~ 'less_steep',
           TRUE ~ 'as_steep'))

#'
#' ### Derivative
#'

ggplot(deriv1_se, aes(x = data, y = derivative)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .5, fill = cpal[4]) +  
  geom_line() + 
  geom_line(aes(color = compared_to_max, group = 1), size = 2, alpha = .8) + 
  scale_color_manual(breaks = c('less_steep', 'as_steep'), values = cpal[c(5,2)], 
                     labels = c('Less steep', 'As steep'),
                     name = 'Region, as compared \nto maximum derivative, \nis...') +
  theme_minimal()
#'
#' ### Difference in derivative (from maximum)
#'
ggplot(deriv1_se, aes(x = data, y = diff)) + 
  geom_ribbon(aes(ymin = diff - crit*se_diff, ymax = diff + crit*se_diff), alpha = .5, fill = cpal[4]) +  
  geom_line() + 
  geom_line(aes(color = compared_to_max, group = 1), size = 2, alpha = .8) + 
  scale_color_manual(breaks = c('less_steep', 'as_steep'), values = cpal[c(5,2)], 
                     labels = c('Less steep', 'As steep'),
                     name = 'Region, as compared \nto maximum derivative, \nis...') +
  theme_minimal()
#'
#' ### Smooth estimate
#'
ggplot(cbind(fit_eval, select(deriv1_se, compared_to_max)), aes(x = Age, y = est)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .5, fill = cpal[4]) +  
  geom_line() + 
  geom_line(aes(color = compared_to_max, group = 1), size = 2, alpha = .8) + 
  scale_color_manual(breaks = c('less_steep', 'as_steep'), values = cpal[c(5,2)], 
                     labels = c('Less steep', 'As steep'),
                     name = 'Region, as compared \nto maximum derivative, \nis...') +
  theme_minimal()

#'
#' ### Derivative (ignoring uncertainty)
#'
deriv1b_se <- gratia::derivatives(fit, n = 100, level = .95, interval = 'confidence', unconditional = FALSE)
i_max <- which(deriv1b_se$derivative == max(deriv1b_se$derivative))
max_est_se <- deriv1b_se[i_max, c('derivative', 'se')]

deriv1b_se <- deriv1b_se %>%
  mutate(diff = derivative - max_est_se$derivative,
         se_diff = sqrt(se^2 + max_est_se$se^2),
         Z = diff/se_diff,
         compared_to_max = case_when(
           Z < -crit ~ 'less_steep',
           TRUE ~ 'as_steep'))

ggplot(deriv1b_se, aes(x = data, y = derivative)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .5, fill = cpal[4]) +  
  geom_line() + 
  geom_line(aes(color = compared_to_max, group = 1), size = 2, alpha = .8) + 
  scale_color_manual(breaks = c('less_steep', 'as_steep'), values = cpal[c(5,2)], 
                     labels = c('Less steep', 'As steep'),
                     name = 'Region, as compared \nto maximum derivative, \nis...') +
  theme_minimal()
#'
#' ### Difference in derivative (from maximum, ignoring uncertainty)
#'
ggplot(deriv1b_se, aes(x = data, y = diff)) + 
  geom_ribbon(aes(ymin = diff - crit*se_diff, ymax = diff + crit*se_diff), alpha = .5, fill = cpal[4]) +  
  geom_line() + 
  geom_line(aes(color = compared_to_max, group = 1), size = 2, alpha = .8) + 
  scale_color_manual(breaks = c('less_steep', 'as_steep'), values = cpal[c(5,2)], 
                     labels = c('Less steep', 'As steep'),
                     name = 'Region, as compared \nto maximum derivative, \nis...') +
  theme_minimal()
#'
#' ### Smooth estimate (ignoring uncertainty)
#'
ggplot(cbind(fit_eval_ignore, select(deriv1b_se, compared_to_max)), aes(x = Age, y = est)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .5, fill = cpal[4]) +  
  geom_line() + 
  geom_line(aes(color = compared_to_max, group = 1), size = 2, alpha = .8) + 
  scale_color_manual(breaks = c('less_steep', 'as_steep'), values = cpal[c(5,2)], 
                     labels = c('Less steep', 'As steep'),
                     name = 'Region, as compared \nto maximum derivative, \nis...') +
  theme_minimal()

#'
#' ## Use posterior of Bayesian smooth {.tabset}
#'
# This may be helpful: https://github.com/paul-buerkner/brms/issues/738
#'
#' ### Derivative
#'
#' Age where the median posterior derivative is at its max is marked.
#'
ggplot(deriv1_post_summary, aes(x = Age, y = deriv1)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .5, fill = cpal[4]) +  
  geom_segment(x = age_at_max_med_deriv, 
               xend = age_at_max_med_deriv, 
               y = deriv1_post_summary[Age == age_at_max_med_deriv, lower],
               yend = deriv1_post_summary[Age == age_at_max_med_deriv, upper]) +
  geom_line() + 
  geom_line(aes(color = compared_to_max, group = 1), size = 2, alpha = .8) + 
  geom_point(x = age_at_max_med_deriv, y = deriv1_post_summary[Age == age_at_max_med_deriv, deriv1], color = cpal[1], size = 2) +
  scale_color_manual(breaks = c('less_steep', 'as_steep'), values = cpal[c(5,2)], 
                     labels = c('Less steep', 'As steep'),
                     name = 'Region, as compared \nto maximum derivative, \nis...') +
  theme_minimal()

#'
#' ### Age at maximum derivative
#'
bayesplot::color_scheme_set('purple')
bayesplot::mcmc_areas(data.frame(max_age = deriv1_age_at_max_post$max_age, Chain = rep(1:4, each = 1000)),
                      prob = 0.95,
                      prob_outer = 0.99)

#'
#' ### Difference from derivative at steepest age
#'
ggplot(deriv1_diff_from_max_post_summary, aes(x = Age, y = diff_from_max)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .5, fill = cpal[4]) +  
  geom_line() + 
  geom_line(aes(color = compared_to_max, group = 1), size = 2, alpha = .8) + 
  scale_color_manual(breaks = c('less_steep', 'as_steep'), values = cpal[c(5,2)], 
                     labels = c('Less steep', 'As steep'),
                     name = 'Region, as compared \nto maximum derivative, \nis...') + 
  theme_minimal()
#'
#' ### Smooth estimate
#'
ggplot(smooth_post_summary[deriv1_diff_from_max_post_summary[, .(Age, compared_to_max)], on = 'Age'],
       aes(x = Age, y = est)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .5, fill = cpal[4]) +  
  geom_line() + 
  geom_line(aes(color = compared_to_max, group = 1), size = 2, alpha = .8) + 
  scale_color_manual(breaks = c('less_steep', 'as_steep'), values = cpal[c(5,2)], 
                     labels = c('Less steep', 'As steep'),
                     name = 'Region, as compared \nto maximum derivative, \nis...') + 
  theme_minimal()
