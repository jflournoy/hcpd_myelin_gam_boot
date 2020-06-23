##############################
# Load relevant Libraries ----
##############################
library(ggplot2) 
library(mgcv) 
library(stringr)  
library(stats) 
library(Formula) 
library(lattice) 
library(gratia) 
library(fda.usc) 
library(colorspace) 
library(data.table)

######################################
# Load imaging and cognitive data ----
######################################
hcpd_data <- readr::read_csv("hcpd_n762_myelin_cog_data.csv")

# Define number of brain regions in parcellation
nreg <- 400

#######################################################
# Estimate Age Curve for Wholebrain T1w/T2w Myelin ----
#######################################################
mean_wholebrain_myelin_gam <- gam(mean_wholebrain_T1wT2w ~ s(Age) + Sex + Scanner + numNavs_sum, data=hcpd_data, method="REML")
summary(mean_wholebrain_myelin_gam)
## Gratia package for derivative visualization and analysis
deriv2_gam <- derivatives(mean_wholebrain_myelin_gam, order = 2, n = 100)
draw(mean_wholebrain_myelin_gam)
draw(deriv_gam)
draw(deriv2_gam)

system.time(blah1 <- replicate(100, {
  deriv_gam <- derivatives(mean_wholebrain_myelin_gam, n = 500)
  deriv_gam_dt <- data.table(deriv_gam)
  deriv_gam_dt[, nearest := derivative]
  setkey(deriv_gam_dt, nearest) # sorts the data
  out <- deriv_gam_dt[J(c(max(derivative), 0)), roll = "nearest"][, .(data, derivative)]
}))/100*1e5/60/32

system.time(blah2 <- replicate(100, {
  deriv_gam <- derivatives(mean_wholebrain_myelin_gam, n = 500)
  out <- rbind(deriv_gam[which.max(deriv_gam$derivative), c('data', 'derivative')],
               deriv_gam[which.min(abs(deriv_gam$derivative)), c('data', 'derivative')])
}))/100*1e5/60/32

## visreg visualization of marginal effects
visreg(mean_wholebrain_myelin_gam, "Age", type="conditional") 

##################################
# Estimate Age of Peak Growth ----
##################################
peak_age <- as.vector(matrix(NA, nrow=nreg, ncol=1))

for (i in 1:nreg){
  m <- as.formula(paste("myelin_schaefer_v", i, "~ s(Age) + Sex + Scanner + numNavs_sum", sep=""))
  
  # print(m)
  m.gam <- eval(substitute(gam(m, data=hcpd_data, method="REML")))
  
  deriv_gam <- derivatives(m.gam)
  deriv_gam$data[which.max(deriv_gam$derivative)]
  peak <- deriv_gam$data[which.max(deriv_gam$derivative)]
  peak_age[i] <- peak
}

hist(peak_age, col="slategray3")


###
# Questions to be answered by bootstrapping:
##
#
# - 95% CI of peak age for each region
# - 95% CI of peak age across regions?
# - Differences in peak age between regions (e.g., bootstrap of peak_age_R1 - peak_age_R2)
# - Differences in slopes at peaks and plateaus, perhaps

library(boot)
gam_deriv_boot <- function(d, i, form){
  d <- d[i,]
  fit <- gam(form, data=d, method="REML")
  deriv1 <- derivatives(fit, n = 500)
  out <- unlist(c(deriv1[which.max(deriv1$derivative), c('data', 'derivative')],
                  deriv1[which.min(abs(deriv1$derivative)), c('data', 'derivative')]))
  names(out) <- c('steep_age', 'steep_deriv', 'plat_age', 'plat_deriv')
  return(out)
}
gam_deriv_boot(hcpd_data,
               1:dim(hcpd_data)[1],
               form = as.formula(mean_wholebrain_T1wT2w ~ s(Age) + Sex + Scanner + numNavs_sum))
system.time({boot_fit <- boot(hcpd_data, 
                 statistic = gam_deriv_boot,
                 R = 1000,
                 sim = 'ordinary',
                 stype = 'i',
                 form = formula(mean_wholebrain_T1wT2w ~ s(Age) + Sex + Scanner + numNavs_sum),
                 parallel = 'multicore',
                 ncpus = 4)})
print(boot_fit)
plot(boot_fit, index = 1)
plot(boot_fit, index = 2)
plot(boot_fit, index = 3)
plot(boot_fit, index = 4)


boot.ci(boot_fit, type = c('basic', 'perc', 'bca'), index = 1)
boot.ci(boot_fit, type = c('basic', 'perc', 'bca'), index = 2)
boot.ci(boot_fit, type = c('basic', 'perc', 'bca'), index = 3)
boot.ci(boot_fit, type = c('basic', 'perc', 'bca'), index = 4)
