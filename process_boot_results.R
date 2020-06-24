####READ DATA, THEN,

ci_dt <- data.table::rbindlist(unlist(ci, recursive = FALSE))
ci_dt_l <- data.table::melt(ci_dt, id.vars = c('index', 'type'))
ci_dt_l[, variable := paste(type, variable, sep = '_')][, type := NULL]
ci_dt <- data.table::dcast(ci_dt_l, index ~ variable)
ci_dt$est <- boot_fit$t0
ci_dt$boot_mean <- apply(boot_fit$t, 2L, mean, na.rm = TRUE)
ci_dt[, bc_est := 2*est - boot_mean]
ci_dt_age <- ci_dt[501:502]
ci_dt_age[, curve := factor(index, levels = 501:502, labels = c('steepest', 'plateau'))]
ci_dt <- ci_dt[-(501:502)]
fit <- mgcv::gam(form, data = d, method = "REML")
deriv1 <- gratia::derivatives(fit, n = 500)
ci_dt$age <- deriv1$data

ci_dt[, curve := dplyr::case_when(
  age < ci_dt_age[1,bca_u] & age > ci_dt_age[1,bca_l] ~ 'steepest',
  age < ci_dt_age[2,bca_u] & age > ci_dt_age[2,bca_l] ~ 'plateau',
  TRUE ~ as.character(index))]

library(ggplot2)
ggplot(ci_dt, aes(x = age)) +
  geom_ribbon(aes(ymin = perc_l, ymax = perc_u), alpha = .1) +
  geom_ribbon(aes(ymin = perc_l, ymax = perc_u, fill = curve), alpha = .5) +
  geom_line(aes(y = est)) +
  geom_line(aes(y = bc_est), color = 'red') +
  scale_fill_discrete(limits = c('plateau', 'steepest'))