# patches
`patches` is an R package to calculate metrics of alhpa, beta, and gamma variability from species level time series data where the data is collected over several patches (plots) in a landscape or region. All metrics are based on equations derived in [Wang and Loreau 2014, Ecology Letters](http://onlinelibrary.wiley.com/doi/10.1111/ele.12292/abstract). 

The function `patches_var()` returns a list of metrics, currently including (note that all descriptions of the metrics described below are *ver batim* from [Wang and Loreau 2014](http://onlinelibrary.wiley.com/doi/10.1111/ele.12292/abstract)):

| Output name | Description |
| ----------- | ----------- |
| `species_synchrony` | Species synchrony within local patches |
| `species_var` | Species-level variability within local patches: squared coefficient of temporal variation of species biomass |
| `plot_species_richness` | Number of species observed in a focal plot, averaged over years |
| `weighted_patch_cv` | Weighted average of coefficients of temporal variation across local patches |
| `alpha_var` | *Alpha* variability: temporal variability at local scale |
| `patch_synchrony` | Spatial synchrony among local patches |
| `beta_var` | Multiplicative *beta* variability: spatial asynchrony or the reciprocal of spatial synchrony |
| `number_of_patches` | Number of unique patches (plots) in the dataset |
| `gamma_var` | *Gamma* variability: temporal variability at the metacommunity scale |

See Table 1 in [Wang and Loreau 2014](http://onlinelibrary.wiley.com/doi/10.1111/ele.12292/abstract) for mathematical details.




