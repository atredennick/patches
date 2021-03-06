% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/patches_funcs.R
\name{beta_var}
\alias{beta_var}
\title{Calculate beta variability; 1/spatial synchrony}
\usage{
beta_var(D, id_var = "plot_id", time_var = "calendar_year",
  species_id = "sp")
}
\arguments{
\item{D}{Dataframe with time series of species level measurements of abundance/biomass across multiple patches. The dataframe can contain many additional columns that will be ignored.}

\item{id_var}{Column name within dataframe that identifies unique patches. Defaults to "plot_id".}

\item{time_var}{Column name within dataframe for the temporal identifier. Defaults to "calendar_year".}

\item{species_id}{Character flag for columns indicating species id. Defaults to "sp".}
}
\value{
Numeric scalar of beta variability
}
\description{
Calculate beta variability; 1/spatial synchrony
}
\author{
Andrew Tredennick
}

