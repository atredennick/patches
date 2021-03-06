% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/patches_funcs.R
\name{community_var}
\alias{community_var}
\title{Variability of community plot i}
\usage{
community_var(D, id_var = "plot_id", time_var = "calendar_year",
  species_id = "sp")
}
\arguments{
\item{D}{Dataframe with time series of species level measurements of abundance/biomass across multiple patches. The dataframe can contain many additional columns that will be ignored.}

\item{id_var}{Column name within dataframe that identifies unique patches. Defaults to "plot_id".}

\item{time_var}{Column name within dataframe for the temporal identifier. Defaults to "calendar_year".}

\item{species_id}{Character flag for columns indicating species id. Defaults to "sp".}
}
\value{
output_df Dataframe of community variability for each plot; CV^2.
}
\description{
Variability of community plot i
}
\author{
Andrew Tredennick
}

