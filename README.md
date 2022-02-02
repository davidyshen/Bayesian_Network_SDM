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
`devtools::source_url("https://github.com/davidyshen/bayesian_network_sdm/blob/main/bnSDM.R?raw=TRUE")`  
or download `bnSDM.R`

### Known issues
Networks with more than ~12 dependencies are not computationally possible (I tried running it over 30 cores and its still not possible 😢 - size of matrix increases by 2^n)

### Function dependencies
* parallel
* raster
* stringr

----
## Other functions
`smush` function compresses many species into a single raster (or more depending on how many you want). Takes the highest probability value in each cell over all input rasters, and selects `n` of them. Used to reduce many species down into a functional group and reduce to a single or few interactors.

### Calling function:
`smush(in_dir, out_dir, out_name, n, ncores = "auto")`
Outputs `n` tifs in the output directory with name *out_name*.tif

### Function inputs:
* `in_dir`    = (required) character string of path to folder containing the input data
* `out_dir`   = (required) character string of output directory
* `out_name`  = (required) character string of output tif filename
* `n`         = integer of how many files to output. If more than 1, takes `n` of the highest probabilities in each cell over all the input rasters. Defaults to `1`
* `ncores`    = number of cores to use. Defaults to `auto` and uses all cores available
