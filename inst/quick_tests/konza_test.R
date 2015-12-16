##
##  R Script to Test Patches Functions on Well-Known Dataset
##
##  Author: Andrew Tredennick
##  Email: atredenn@gmail.com
##  Date created: 12-16-2015
##

## Clear the workspace...
rm(list=ls())



####
####  Load Libraries
####
library(patches)
library(plyr)
library(reshape2)
library(dplyr)



####
####  Load Konza-RaMPS Dataset
####
setwd("/Users/atredenn/Dropbox/Patches/Data_Nov2015/Data/")
konza_data <- read.csv("abs abundance_sp as cols_incomplete plots or yrs rmvd - best_Nov2015.csv")%>%
  filter(site_proj_comm=="KNZ_RaMPs_0")
konza_data <- konza_data[,2:ncol(konza_data)]

knzm <- melt(konza_data, id.vars = c("site_proj_comm", "calendar_year", "block",
                                     "treatment", "plot_id", "community_type"))
knzag <- ddply(knzm, .(calendar_year, plot_id), summarise,
               plot_cover = sum(value))

plot_area <- 100*100 #1m^2 plot area in cm
knzm$value <- with(knzm, (value/100)*plot_area)
konza_data_area <- dcast(knzm, site_proj_comm+calendar_year+block+treatment+plot_id+community_type~variable, value.var = "value")

species_synchrony(konza_data)
species_synchrony(konza_data_area)
ind_spp_var(konza_data)
ind_spp_var(konza_data_area)
gamma_var(konza_data)
gamma_var(konza_data_area)
community_var(konza_data)
avg_spp_var(konza_data)
avg_spp_var(konza_data_area)
alpha_var(konza_data)
alpha_var(konza_data_area)
