
#' Calculate species synchrony within each plot
#' @author Andrew Tredennick
#' @param D Dataframe with time series of species level measurements of abundance/biomass across multiple patches. The dataframe can contain many additional columns that will be ignored.
#' @param id_var Column name within dataframe that identifies unique patches. Defaults to "plot_id".
#' @param time_var Column name within dataframe for the temporal identifier. Defaults to "calendar_year".
#' @param species_id Character flag for columns indicating species id. Defaults to "sp".
#' @return species_synchrony Synchrony of species' cover through time in each observed plot.
species_synchrony <- function(D,
                              id_var="plot_id", 
                              time_var="calendar_year", 
                              species_id="sp"){
  require(plyr)
  require(dplyr)
  
  ####
  ####  Housekeeping for column names
  ####
  colnames(D)[which(colnames(D)==id_var)] <- "id_var"
  colnames(D)[which(colnames(D)==time_var)] <- "time_var"
  num_plots <- length(unique(D$id_var))
  plot_ids <- unique(D$id_var)
  
  ####
  ####  Test to make sure plots occur at all time points
  ####
  yrs <- unique(D$time_var)
  yr_plots <- list()
  for(y in 1:length(yrs)){
    tmp <- subset(D, time_var==yrs[y])
    yr_plots <- unique(tmp$id_var)
    if(all(plot_ids %in% yr_plots)==FALSE){
      stop("plot numbers are uneven through years; \nconsider breaking up the analysis by years")
    } # end T/F
  } # end years loop
  
  ####
  ####  Within patch species synchrony
  ####
  species_columns <- colnames(D)[grep(species_id, colnames(D))]
  species_columns <- species_columns[-grep("species", species_columns)]
  species_synchrony <- matrix(ncol=2, nrow=num_plots)
  for(plot_now in 1:num_plots){
    tmp_data <- subset(D, id_var==plot_ids[plot_now])
    comm_matrix <- tmp_data[,species_columns]
    # Remove species columns that only include 0s
    keepers <- which(colSums(comm_matrix)!=0)
    comm_matrix <- comm_matrix[,keepers]
    
    cov_matrix <- cov(comm_matrix)
    intra_cov_vector <- diag(cov_matrix)
    # diag(cov_matrix) <- 0
    summed_inter <- sum(cov_matrix)
    summed_sq_intra <- sum(sqrt(intra_cov_vector))
    
    # Save output
    species_synchrony[plot_now,1] <- plot_ids[plot_now]
    species_synchrony[plot_now,2] <- summed_inter / summed_sq_intra^2
  } # end plot loop
  
  return(species_synchrony)
  
} # end species synchrony function



#' Calculate species variability within each plot, for each species
#' @author Andrew Tredennick
#' @param D Dataframe with time series of species level measurements of abundance/biomass across multiple patches. The dataframe can contain many additional columns that will be ignored.
#' @param id_var Column name within dataframe that identifies unique patches. Defaults to "plot_id".
#' @param time_var Column name within dataframe for the temporal identifier. Defaults to "calendar_year".
#' @param species_id Character flag for columns indicating species id. Defaults to "sp".
#' @return output_df Dataframe of temporal coefficient of variation for each species in each plot.
ind_spp_var <- function(D, 
                        id_var="plot_id", 
                        time_var="calendar_year", 
                        species_id="sp"){
  require(plyr)
  require(dplyr)
  
  ####
  ####  Housekeeping for column names
  ####
  colnames(D)[which(colnames(D)==id_var)] <- "id_var"
  colnames(D)[which(colnames(D)==time_var)] <- "time_var"
  num_plots <- length(unique(D$id_var))
  plot_ids <- unique(D$id_var)
  
  ####
  ####  Test to make sure plots occur at all time points
  ####
  yrs <- unique(D$time_var)
  yr_plots <- list()
  for(y in 1:length(yrs)){
    tmp <- subset(D, time_var==yrs[y])
    yr_plots <- unique(tmp$id_var)
    if(all(plot_ids %in% yr_plots)==FALSE){
      stop("plot numbers are uneven through years; \nconsider breaking up the analysis by years")
    } # end T/F
  } # end years loop
  
  # Get index of species columns
  species_columns <- colnames(D)[grep(species_id, colnames(D))]
  species_columns <- species_columns[-grep("species", species_columns)]
  
  ####
  ####  Variability of j-th species in patch i (CV_j(i)^2)
  ####
  spp_within_plot_cv <- list()
  for(i in 1:num_plots){
    tmp_data <- subset(D, id_var==plot_ids[i])
    comm_data <- tmp_data[,species_columns]
    # Remove species columns that only include 0s
    keepers <- which(colSums(comm_data)!=0)
    comm_data <- comm_data[,keepers]
    cv_ji <- apply(comm_data, MARGIN = 2, sd)/apply(comm_data, MARGIN = 2, mean)
    cv_df <- data.frame(species=colnames(comm_data),
                        cv = cv_ji)
    spp_within_plot_cv[[as.character(plot_ids[i])]] <- cv_df
  } # end plots loop
  
  output_df <- melt(data = spp_within_plot_cv, id_vars="species")
  colnames(output_df) <- c("species_id", "variable", "value", "plot_id")
  return(output_df)
  
} # end function



#' Calculate alpha, beta, and gamma variability metrics from community matrix
#' @author Andrew Tredennick
#' @param D Dataframe with time series of species level measurements of abundance/biomass across multiple patches. The dataframe can contain many additional columns that will be ignored.
#' @param id_var Column name within dataframe that identifies unique patches. Defaults to "plot_id".
#' @param time_var Column name within dataframe for the temporal identifier. Defaults to "calendar_year".
#' @param species_id Character flag for columns indicating species id. Defaults to "sp".

patches_var <- function(D, 
                        id_var="plot_id", 
                        time_var="calendar_year", 
                        species_id="sp"){
  require(plyr)
  require(dplyr)
  require(synchrony)
  
  ####
  ##  Housekeeping for column names
  ####
  colnames(D)[which(colnames(D)==id_var)] <- "id_var"
  colnames(D)[which(colnames(D)==time_var)] <- "time_var"
  num_plots <- length(unique(D$id_var))
  plot_ids <- unique(D$id_var)
  
  ####
  ##  Test to make sure plots occur at all time points
  ####
  yrs <- unique(D$time_var)
  yr_plots <- list()
  for(y in 1:length(yrs)){
    tmp <- subset(D, time_var==yrs[y])
    yr_plots <- unique(tmp$id_var)
    if(all(plot_ids %in% yr_plots)==FALSE){
      stop("plot numbers are uneven through years; \nconsider breaking up the analysis by years")
    }
  }
  
  ####
  ##  Within patch species synchrony
  ####
  species_columns <- colnames(D)[grep(species_id, colnames(D))]
  species_columns <- species_columns[-grep("species", species_columns)]
  species_synchrony <- matrix(ncol=2, nrow=num_plots)
  for(plot_now in 1:num_plots){
    tmp_data <- subset(D, id_var==plot_ids[plot_now])
    comm_matrix <- tmp_data[,species_columns]
    species_synchrony[plot_now,1] <- plot_ids[plot_now]
    species_synchrony[plot_now,2] <- as.numeric(community.sync(comm_matrix)[1])
  }
  
  
  ####
  ##  Species level variability within local patches: squared coefficient of 
  ##  temporal variation of species biomass
  ####
  Dtmp <- D[,c("id_var", "time_var", species_columns)]
  Dmelt <- melt(Dtmp, id.vars = c("id_var", "time_var"))
  species_cvs <- ddply(Dmelt, .(id_var, variable), summarise,
                       cv_square = (sd(value)/mean(value))^2)
  species_cvs <- species_cvs[complete.cases(species_cvs),]
  ##TODO: calculate weighted average species CV across patches???
  
  
  ####
  ##  Patch level richness and evenness, averaged over time
  ####
  Dtmp <- D[,c("id_var", "time_var", species_columns)]
  Dmelt <- melt(Dtmp, id.vars = c("id_var", "time_var"))
  species_abundance <- ddply(Dmelt, .(id_var, time_var), summarise,
                             abundance = length(which(value>0)))
  avg_spp_rich <- ddply(species_abundance, .(id_var), summarise,
                        avg_richness = mean(abundance))
  ##TODO: add in species evenness metric
  
  
  ####
  ##  Patch synchrony and CVs
  ####
  spp_long <- melt(D, measure.vars = species_columns)
  patch_abund <- ddply(spp_long, .(time_var, id_var), summarise,
                     comm_anpp = sum(value))
  patch_cvs <- ddply(patch_abund , .(id_var), summarise,
                     cv = sd(comm_anpp)/mean(comm_anpp),
                     sds = sd(comm_anpp))
  total_metacomm <- ddply(spp_long, .(time_var), summarise,
                          total_biomass = sum(value))
  mean_metacomm <- mean(total_metacomm$total_biomass)
  weighted_local_cv <- sum(patch_cvs$sds) / mean_metacomm
  
  comm_matrix <- dcast(patch_abund, time_var~id_var, value.var="comm_anpp")
  w <- cov(comm_matrix[,2:ncol(comm_matrix)])
  metacomm_cv <- sqrt(sum(w))/mean_metacomm
  patch_synchrony <- metacomm_cv^2 / weighted_local_cv^2
  
  
  ####
  ##  Gamma (regional) variability
  ####
  gamma_var <- sd(total_metacomm$total_biomass)/mean(total_metacomm$total_biomass)
  gamma_test <- weighted_local_cv^2 / (1/patch_synchrony)
  if(round(gamma_var^2,2) != round(gamma_test,2)){
    warning("alpha and beta variability not producing gamma variability as expected")
  }
  
  ####
  ##  Collate output
  ####
  out_list <- list(species_synchrony = species_synchrony,
                   species_var = species_cvs,
                   plot_species_richness = avg_spp_rich,
                   weighted_patch_cv = weighted_local_cv,
                   alpha_var = weighted_local_cv^2,
                   patch_synchrony = patch_synchrony,
                   beta_var = 1/patch_synchrony,
                   number_of_patches = num_plots,
                   gamma_var = gamma_var^2)
  return(out_list)
}

