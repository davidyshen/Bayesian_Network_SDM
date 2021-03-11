# Bayesian Network SDM function
A function at applies a bayesian network to SDMs. Based on Staniczenko et al. (2017), [doi:10.1111/ele.12770](https://doi.org/10.1111/ele.12770)

----
For optimal performance this should be run in parallel  
More Cores = More Better  
Currently uses multisession because I'm writing this on Windows but will change to forking later

----

### Calling function
`bnSDM(in_dir, out_dir = "BN_out/", focal, direction, method = "or", ncores = "auto")`  
Outputs a tif in the output directory named *focal.species*_bnSDM.tif

### Function inputs:
* `in_dir`    = (required) character string of directory of folder ending in "/" containing .tif input raster SDMs of all species e.g. `"folder/"`
* `out_dir`   = character string ending in "/" for output directory, defaults to `"BN_out/"`
* `focal`     = (required) character string file name of focal species inside in_dir folder including file extension e.g. `"Gonimbrasia_belina.tif"`
* `direction` = (required) vector of 1s (positive interaction) and -1s (negative interactions) that denote the direction of interaction of interacting species on focal species in alphabetical order e.g. `c(1, -1, 1)`
* `method`    = character string `"or"` (default) or `"and"` specifies state table filling method OR vs AND (box 2 of Staniczenko et al)
* `ncores`    = integer value or character string `auto`, specifies number of cores/threads to run in multithreaded mode, defaults to `auto` using all threads on all cores

SDM rasters are loaded and used alphabetically, so the `direction` vector must be in the same order

### Calling
`devtools::source_url("https://github.com/davidyshen/bayesian_network_sdm/blob/main/bnSDM_rewrite.R?raw=TRUE")`  
or download `bnSDM_rewrite.R`

### Known issues
Networks with more than ~12 dependencies are not computationally possible (I tried running it over 30 cores and its still not possible 😢 - size of matrix increases by 2^n)

### Function dependencies
* parallel
* raster
* ~~gRain~~
* ~~dplyr~~

~~Some dependencies of gRain (RBGL) are no-longer hosted in CRAN, however can be obtained through Bioconductor~~ No longer requires gRain package

Future steps: Reinclude gRain package to draw interaction network plot
