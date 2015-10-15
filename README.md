# patches
[![Build Status](https://travis-ci.org/atredennick/patches.svg?branch=master)](https://travis-ci.org/atredennick/patches)

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

##License
The MIT License (MIT)

Copyright (c) 2015 Andrew Tredennick

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


