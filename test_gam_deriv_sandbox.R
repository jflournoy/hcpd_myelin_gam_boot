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

library(mgcv) 
library(gratia) 
library(boot)
library(dplyr)

hcpd_data <- readr::read_csv("hcpd_n762_myelin_cog_data.csv")

right_side <- '~ s(Age) + Sex + Scanner + numNavs_sum'
left_side <- 'mean_wholebrain_T1wT2w'

mf_form <- as.formula(paste0(left_side, gsub('s\\((.*)\\)', '\\1', right_side)))
form <- as.formula(paste0(left_side, right_side))

d <- model.frame(formula = mf_form, data = hcpd_data)

fit <- mgcv::gam(form, data=d, method="REML")
deriv1 <- gratia::derivatives(fit, n = 100, level = .95, interval = 'simultaneous', unconditional = TRUE, n_sim = 10000)
deriv1aa <- gratia::derivatives(fit, n = 100, level = .95, interval = 'simultaneous', unconditional = FALSE, n_sim = 10000)
deriv1a <- gratia::derivatives(fit, n = 100, level = .95, interval = 'confidence', unconditional = TRUE)
deriv1b <- gratia::derivatives(fit, n = 100, level = .95, interval = 'confidence', unconditional = FALSE)


ggplot(deriv1, aes(x = data, y = derivative)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25, color = 'black') + 
  geom_ribbon(data = deriv1aa, aes(ymin = lower, ymax = upper), alpha = .25, color = 'yellow') + 
  geom_ribbon(data = deriv1a, aes(ymin = lower, ymax = upper), alpha = .25, color = 'blue') + 
  geom_ribbon(data = deriv1b, aes(ymin = lower, ymax = upper), alpha = .25, color = 'red') + 
  geom_line()

fit_eval <- gratia:::confint.gam(fit, parm = 's(Age)', n = 100, unconditional = TRUE, type = 'simultaneous', nsim = 10000)

#' Test where the GAM smooth derivative becomes different than the peak. We'll
#' explore a few different options:
#' 
#' 1. difference in 83% CIs
#' 2. compute difference in point estimates using pooled standard errors
#' 3. use brms to derive a quantity from the posterior

#'
#' # difference in 83% CIs
#'

deriv1_83 <- gratia::derivatives(fit, n = 100, level = .83, interval = 'simultaneous', unconditional = TRUE, n_sim = 10000)
i_max <- which(deriv1_83$derivative == max(deriv1_83$derivative))
max_lower <- deriv1_83$lower[i_max]

deriv1_83 <- deriv1_83 %>%
  mutate(compared_to_max = case_when(
    upper < max_lower ~ 'less_steep',
    TRUE ~ 'as_steep'
  ))

ggplot(deriv1_83, aes(x = data, y = derivative)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) + 
  geom_line() + 
  geom_line(aes(color = compared_to_max), size = 2, alpha = .8)

ggplot(cbind(fit_eval, select(deriv1_83, compared_to_max)), aes(x = Age, y = est)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .15) + 
  geom_line() + 
  geom_line(aes(color = compared_to_max), size = 2, alpha = .8) + 
  theme_minimal()

any(deriv1_83$compared_to_max == 'less_steep')

#'
#' What if we ignore all the uncertainty and multiple comparisons?
#'

deriv1b_83 <- gratia::derivatives(fit, n = 100, level = .83, interval = 'confidence', unconditional = FALSE)
i_max <- which(deriv1b_83$derivative == max(deriv1b_83$derivative))
max_lower <- deriv1b_83$lower[i_max]

deriv1b_83 <- deriv1b_83 %>%
  mutate(compared_to_max = case_when(
    upper < max_lower ~ 'less_steep',
    TRUE ~ 'as_steep'
  ))

any(deriv1b_83$compared_to_max == 'less_steep')

ggplot(deriv1b_83, aes(x = data, y = derivative)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) + 
  geom_line() + 
  geom_line(aes(color = compared_to_max), size = 2, alpha = .8)

ggplot(cbind(fit_eval, select(deriv1b_83, compared_to_max)), aes(x = Age, y = est)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .15) + 
  geom_line() + 
  geom_line(aes(color = compared_to_max), size = 2, alpha = .8) + 
  theme_minimal()

#'
#' #  compute difference in point estimates using pooled standard errors
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

ggplot(deriv1_se, aes(x = data, y = derivative)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) + 
  geom_line() + 
  geom_line(aes(color = compared_to_max), size = 2, alpha = .8)

ggplot(deriv1_se, aes(x = data, y = diff)) + 
  geom_ribbon(aes(ymin = diff - crit*se_diff, ymax = diff + crit*se_diff), alpha = .25) + 
  geom_line() + 
  geom_line(aes(color = compared_to_max), size = 2, alpha = .8)

ggplot(cbind(fit_eval, select(deriv1_se, compared_to_max)), aes(x = Age, y = est)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .15) + 
  geom_line() + 
  geom_line(aes(color = compared_to_max), size = 2, alpha = .8) + 
  theme_minimal()

#'
#' What if we ignore all the uncertainty and multiple comparisons?
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
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) + 
  geom_line() + 
  geom_line(aes(color = compared_to_max), size = 2, alpha = .8)

ggplot(deriv1b_se, aes(x = data, y = diff)) + 
  geom_ribbon(aes(ymin = diff - crit*se_diff, ymax = diff + crit*se_diff), alpha = .25) + 
  geom_line() + 
  geom_line(aes(color = compared_to_max), size = 2, alpha = .8)

ggplot(cbind(fit_eval, select(deriv1b_se, compared_to_max)), aes(x = Age, y = est)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .15) + 
  geom_line() + 
  geom_line(aes(color = compared_to_max), size = 2, alpha = .8) + 
  theme_minimal()

#'
#' # use brms to derive a quantity from the posterior
#'
#' This may be helpful: https://github.com/paul-buerkner/brms/issues/738
#'
