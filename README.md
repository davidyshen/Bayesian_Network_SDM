# Bayesian Network SDM function
A function at applies a bayesian network to SDMs. Based on Staniczenko et al. (2017), [doi:10.1111/ele.12770](https://doi.org/10.1111/ele.12770)
### Calling function
`bnSDM(in_dir, out_dir, focal, direction, method)`

Outputs a tif in the output directory named *focal.species*_bnSDM.tif

### Function inputs:
* `in_dir`    = character string of directory of folder ending in "/" containing .tif input raster SDMs of all species e.g. `"folder/"`
* `out_dir`   = character string ending in "/" for output directory, defaults to `"BN_out/"`
* `focal`     = character string file name of focal species inside in_dir folder including file extension e.g. `"Gonimbrasia_belina.tif"`
* `direction` = vector of 1s (positive interaction) and -1s (negative interactions) that denote the direction of interaction of interacting species on focal species in alphabetical order e.g. `c(1, -1, 1)`
* `method`    = character string `"or"` (default) or `"and"` specifies state table filling method OR vs AND (box 2 of Staniczenko et al)
* `ncores`    = integer value or character string `auto`, specifies number of cores/threads to run in multithreaded mode, defaults to `auto` using all threads on all cores

SDM rasters are loaded and used alphabetically, so the `direction` vector must be in the same order

### Calling
`devtools::source_url("https://github.com/davidyshen/bayesian_network_sdm/blob/main/bnSDM_rewrite.R?raw=TRUE")`  
or download `bnSDM_rewrite.R`

### Function dependencies
* raster
* ~~gRain~~
* dplyr
* parallel

~~Some dependencies of gRain (RBGL) are no-longer hosted in CRAN, however can be obtained through Bioconductor~~ No longer requires gRain package

