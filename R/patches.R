####
##  Alpha variability
####
#' Calculate alpha variability from community matrix
#' @author Andrew Tredennick
#' @param D Dataframe with time series of species 
#'          level measurements of abundance across multiple patches.
#' @param id_var Column name within dataframe that identifies unique patches.
#' @param measurement_var Column name within dataframe that identifies the 
#'                        measurement of abundance used.
#' @param time_var Column name within dataframe for the temporal identifier.

# Weighted average of coefficients of temporal variation across local patches
alpha_var <- function(D, id_var, measurement_var, time_var){
  patch_cvs <- ddply(D, .(id_var), summarise,
                     cv = sd(measurement_var)/mean(measurement_var),
                     sds = sd(measurement_var))
  total_metacomm <- ddply(patch_anpp, .(time_var), summarise,
                          total_biomass = sum(comm_anpp))
  mean_metacomm <- mean(total_metacomm$total_biomass)
  weighted_local_cv <- sum(patch_cvs$sds) / mean_metacomm
  return(weighted_local_cv)
}