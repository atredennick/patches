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

species_synchrony(konza_data)



