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
#' @param species_id Character flag for columns indicating species id. Defaults
#'                   to "sp".

# Weighted average of coefficients of temporal variation across local patches
patches_var <- function(D, id_var, measurement_var, time_var, species_id="sp"){
  # Within in patch synchrony (species synchrony)
  num_plots <- length(unique(D$id_var))
  plot_ids <- unique(D$id_var)
  species_columns <- colnames(D)[grep(species_id, colnames(D))]
  species_columns <- species_columns[-grep("species", species_columns)]
  species_synchrony <- matrix(ncol=2, nrow=num_plots)
  for(plot_now in 1:num_plots){
    tmp_data <- subset(D, plot_id==plot_ids[plot_now])
    comm_matrix <- tmp_data[,species_columns]
    species_synchrony[plot_now,1] <- plot_ids[plot_now]
    species_synchrony[plot_now,2] <- as.numeric(community.sync(comm_matrix)[1])
  }
  
  patch_cvs <- ddply(D, .(id_var), summarise,
                     cv = sd(measurement_var)/mean(measurement_var),
                     sds = sd(measurement_var))
  total_metacomm <- ddply(patch_anpp, .(time_var), summarise,
                          total_biomass = sum(comm_anpp))
  mean_metacomm <- mean(total_metacomm$total_biomass)
  weighted_local_cv <- sum(patch_cvs$sds) / mean_metacomm
  return(weighted_local_cv)
}

