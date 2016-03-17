
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
  tmprms <- grep("species", species_columns) # check for extra species col
  if(length(tmprms)>0){
    species_columns <- species_columns[-tmprms]
  }
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
  
  species_synchrony <- as.data.frame(species_synchrony)
  colnames(species_synchrony) <- c("plot", "spp_synch")
  return(species_synchrony)
  
} # end species synchrony function



#' Calculate species variability within each plot, for each species
#' @author Andrew Tredennick
#' @param D Dataframe with time series of species level measurements of abundance/biomass across multiple patches. The dataframe can contain many additional columns that will be ignored.
#' @param id_var Column name within dataframe that identifies unique patches. Defaults to "plot_id".
#' @param time_var Column name within dataframe for the temporal identifier. Defaults to "calendar_year".
#' @param species_id Character flag for columns indicating species id. Defaults to "sp".
#' @return Dout Dataframe of temporal coefficient of variation for each species in each plot.
ind_spp_var <- function(D, 
                        id_var="plot_id", 
                        time_var="calendar_year", 
                        species_id="sp"){
  require(plyr)
  require(dplyr)
  require(reshape2)
  
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
  tmprms <- grep("species", species_columns) # check for extra species col
  if(length(tmprms)>0){
    species_columns <- species_columns[-tmprms]
  }
  
  ####
  ####  Variability of j-th species in patch i (CV_j(i)^2)
  ####
  Dtmp <- D[,c("id_var", "time_var", species_columns)]
  Dmelt <- melt(Dtmp, id.vars = c("id_var", "time_var"))
  Dout <- ddply(Dmelt, .(id_var, variable), summarise,
                avg = mean(value),
                sdev = sd(value),
                cv = sd(value)/mean(value))
  colnames(Dout) <- c("plot", "species", "avg", "sdev", "cv")
  Dout <- Dout[which(Dout$avg>0),]
  return(Dout)
  
} # end function





#' Calculate averaged species variability within each plot
#' @author Andrew Tredennick
#' @param D Dataframe with time series of species level measurements of abundance/biomass across multiple patches. The dataframe can contain many additional columns that will be ignored.
#' @param id_var Column name within dataframe that identifies unique patches. Defaults to "plot_id".
#' @param time_var Column name within dataframe for the temporal identifier. Defaults to "calendar_year".
#' @param species_id Character flag for columns indicating species id. Defaults to "sp".
#' @return output_df Dataframe of averaged temporal coefficient of variation for each plot.
avg_spp_var <- function(D, 
                        id_var="plot_id", 
                        time_var="calendar_year", 
                        species_id="sp"){
  require(plyr)
  require(dplyr)
  require(reshape2)
  
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
  tmprms <- grep("species", species_columns) # check for extra species col
  if(length(tmprms)>0){
    species_columns <- species_columns[-tmprms]
  }
  
  
  ####
  ####  Variability of j-th species in patch i (CV_j(i)^2)
  ####
  Dtmp <- D[,c("id_var", "time_var", species_columns)]
  Dmelt <- melt(Dtmp, id.vars = c("id_var", "time_var"))
  Dout <- ddply(Dmelt, .(id_var, variable), summarise,
                avg = mean(value),
                sdev = sd(value),
                cv = sd(value)/mean(value))
  colnames(Dout) <- c("plot", "species", "avg", "sdev", "cv")
  Dout <- Dout[which(Dout$avg>0),] #remove spp absent throughout time series
  
  # Get mean community biomass per plot
  Dplot <- Dmelt[which(Dmelt$value>0),]
  plot_biomass <- ddply(Dplot, .(id_var, time_var), summarise,
                        tot_plot_biom = sum(value))
  avg_plot_biomass <- ddply(plot_biomass, .(id_var), summarise,
                             avg_biomass = mean(tot_plot_biom))
  names(avg_plot_biomass) <- c("plot", "avg_plot_biomass")
  
  Dout <- merge(Dout, avg_plot_biomass) 
  
  # Make spp weights
  Dout$sppweights <- with(Dout, avg/avg_plot_biomass)
  Dout$var <- with(Dout, sppweights*cv)
  
  output_df <- ddply(Dout, .(plot), summarise,
                     sum(var))
  names(output_df) <- c("plot", "avgd_spp_var")
  return(output_df)
  
} # end function







#' Variability of community plot i
#' @author Andrew Tredennick
#' @param D Dataframe with time series of species level measurements of abundance/biomass across multiple patches. The dataframe can contain many additional columns that will be ignored.
#' @param id_var Column name within dataframe that identifies unique patches. Defaults to "plot_id".
#' @param time_var Column name within dataframe for the temporal identifier. Defaults to "calendar_year".
#' @param species_id Character flag for columns indicating species id. Defaults to "sp".
#' @return output_df Dataframe of community variability for each plot; CV^2.
community_var <- function(D, 
                        id_var="plot_id", 
                        time_var="calendar_year", 
                        species_id="sp"){
  require(plyr)
  require(dplyr)
  require(reshape2)
  
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
  tmprms <- grep("species", species_columns) # check for extra species col
  if(length(tmprms)>0){
    species_columns <- species_columns[-tmprms]
  }
  
  ####
  ####  Within patch species synchrony
  ####
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
  
  ####
  ####  Variability of j-th species in patch i (CV_j(i)^2)
  ####
  Dtmp <- D[,c("id_var", "time_var", species_columns)]
  Dmelt <- melt(Dtmp, id.vars = c("id_var", "time_var"))
  Dout <- ddply(Dmelt, .(id_var, variable), summarise,
                avg = mean(value),
                sdev = sd(value),
                cv = sd(value)/mean(value))
  colnames(Dout) <- c("plot", "species", "avg", "sdev", "cv")
  Dout <- Dout[which(Dout$avg>0),] #remove spp absent throughout time series
  
  # Get mean community biomass per plot
  Dplot <- Dmelt[which(Dmelt$value>0),]
  plot_biomass <- ddply(Dplot, .(id_var, time_var), summarise,
                        tot_plot_biom = sum(value))
  avg_plot_biomass <- ddply(plot_biomass, .(id_var), summarise,
                            avg_biomass = mean(tot_plot_biom))
  names(avg_plot_biomass) <- c("plot", "avg_plot_biomass")
  
  Dout <- merge(Dout, avg_plot_biomass) 
  
  # Make spp weights
  Dout$sppweights <- with(Dout, avg/avg_plot_biomass)
  Dout$var <- with(Dout, sppweights*cv)
  var_df <- ddply(Dout, .(plot), summarise,
                     sum(var))
  avgd_spp_var <- var_df[,2]
  
  var_comm <- (avgd_spp_var^2)*species_synchrony[,2]
  output_df <- data.frame(plot=plot_ids, var_comm_sq=var_comm)
  return(output_df)
  
} # end function






#' Alpha variability
#' @author Andrew Tredennick
#' @param D Dataframe with time series of species level measurements of abundance/biomass across multiple patches. The dataframe can contain many additional columns that will be ignored.
#' @param id_var Column name within dataframe that identifies unique patches. Defaults to "plot_id".
#' @param time_var Column name within dataframe for the temporal identifier. Defaults to "calendar_year".
#' @param species_id Character flag for columns indicating species id. Defaults to "sp".
#' @return output Numeric scalar for alpha variability (CV) at local scale.
alpha_var <- function(D, 
                        id_var="plot_id", 
                        time_var="calendar_year", 
                        species_id="sp"){
  require(plyr)
  require(dplyr)
  require(reshape2)
  
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
  tmprms <- grep("species", species_columns) # check for extra species col
  if(length(tmprms)>0){
    species_columns <- species_columns[-tmprms]
  }
  
  # Get mean community biomass, cv per plot
  spp_long <- melt(D, measure.vars = species_columns)
  patch_abund <- ddply(spp_long, .(time_var, id_var), summarise,
                       comm_anpp = sum(value))
  patch_cvs <- ddply(patch_abund , .(id_var), summarise,
                     cv = sd(comm_anpp)/mean(comm_anpp),
                     sds = sd(comm_anpp))
  total_metacomm <- ddply(patch_abund, .(time_var), summarise,
                          total_biomass = sum(comm_anpp))
#   if(var_type=="cover"){
#     total_metacomm <- ddply(patch_abund, .(time_var), summarise,
#                             total_biomass = mean(comm_anpp))
#   }
#   if(var_type=="biomass"){
#     total_metacomm <- ddply(patch_abund, .(time_var), summarise,
#                             total_biomass = sum(comm_anpp))
#   }
  
  mean_metacomm <- mean(total_metacomm$total_biomass)
  weighted_local_cv <- sum(patch_cvs$sds) / mean_metacomm
  
  alpha_cv_sq <- weighted_local_cv^2
  
  return(data.frame(alpha_var = alpha_cv_sq,
                    weighted_local_cv = weighted_local_cv))
  
} # end function






#' Calculate beta variability; 1/spatial synchrony
#' @author Andrew Tredennick
#' @param D Dataframe with time series of species level measurements of abundance/biomass across multiple patches. The dataframe can contain many additional columns that will be ignored.
#' @param id_var Column name within dataframe that identifies unique patches. Defaults to "plot_id".
#' @param time_var Column name within dataframe for the temporal identifier. Defaults to "calendar_year".
#' @param species_id Character flag for columns indicating species id. Defaults to "sp".
#' @return Numeric scalar of beta variability
beta_var <- function(D,
                     id_var="plot_id", 
                     time_var="calendar_year", 
                     species_id="sp"){
  require(plyr)
  require(dplyr)
  
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

  
  species_columns <- colnames(D)[grep(species_id, colnames(D))]
  tmprms <- grep("species", species_columns) # check for extra species col
  if(length(tmprms)>0){
    species_columns <- species_columns[-tmprms]
  }
  
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
  
  return(data.frame(beta_variability = 1/patch_synchrony,
                    patch_synchrony = patch_synchrony))
} # end beta var function





#' Calculate gamma variability
#' @author Andrew Tredennick
#' @param D Dataframe with time series of species level measurements of abundance/biomass across multiple patches. The dataframe can contain many additional columns that will be ignored.
#' @param id_var Column name within dataframe that identifies unique patches. Defaults to "plot_id".
#' @param time_var Column name within dataframe for the temporal identifier. Defaults to "calendar_year".
#' @param species_id Character flag for columns indicating species id. Defaults to "sp".
#' @return Numeric scalar of gamma variability 
gamma_var <- function(D, 
                      id_var="plot_id", 
                      time_var="calendar_year", 
                      species_id="sp"){
  require(plyr)
  require(dplyr)
  
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
  
  species_columns <- colnames(D)[grep(species_id, colnames(D))]
  tmprms <- grep("species", species_columns) # check for extra species col
  if(length(tmprms)>0){
    species_columns <- species_columns[-tmprms]
  }
  
  spp_long <- melt(D, measure.vars = species_columns)
  patch_abund <- ddply(spp_long, .(time_var, id_var), summarise,
                       comm_anpp = sum(value))
  patch_cvs <- ddply(patch_abund , .(id_var), summarise,
                     cv = sd(comm_anpp)/mean(comm_anpp),
                     sds = sd(comm_anpp))
  total_metacomm <- ddply(spp_long, .(time_var), summarise,
                          total_biomass = sum(value))
  
  gamma_cv <- sd(total_metacomm$total_biomass)/mean(total_metacomm$total_biomass)
  alpha_variability <- alpha_var(D = D)$alpha_var
  beta_variability <- beta_var(D = D)$beta_variability
  gamma_test <- alpha_variability / beta_variability
  if(round(gamma_cv^2,2) != round(gamma_test,2)){
    warning("alpha and beta variability not producing gamma variability as expected")
  }
  
  return(gamma_cv^2)

} # end gamma variability function
