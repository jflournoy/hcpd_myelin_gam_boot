##############################
# Load relevant Libraries ----
##############################
rm(list = ls())

if(!require("ggplot2")) {install.packages("ggplot2");require("ggplot2")} 
if(!require("mgcv")) {install.packages("mgcv");require("mgcv")} 
if(!require("stringr")) {install.packages("stringr");require("stringr")}  
if(!require("stats")) {install.packages("stats");require("stats")} 
if(!require("Formula")) {install.packages("Formula");require("Formula")} 
if(!require("lattice")) {install.packages("lattice");require("lattice")} 
if(!require("gratia")) {install.packages("gratia");require("gratia")} 
if(!require("fda.usc")) {install.packages("fda.usc");require("fda.usc")} 
if(!require("colorspace")) {install.packages("colorspace");require("colorspace")} 

######################################
# Load imaging and cognitive data ----
######################################
hcpd_data <- read.csv("~/Downloads/hcpd_n762_myelin_cog_data.csv")

# Define number of brain regions in parcellation
nreg <- 400

#######################################################
# Estimate Age Curve for Wholebrain T1w/T2w Myelin ----
#######################################################
mean_wholebrain_myelin_gam <- gam(mean_wholebrain_T1wT2w ~ s(Age) + Sex + Scanner + numNavs_sum, data=hcpd_data, method="REML")

## Gratia package for derivative visualization and analysis
deriv_gam <- derivatives(mean_wholebrain_myelin_gam)
draw(deriv_gam)

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
