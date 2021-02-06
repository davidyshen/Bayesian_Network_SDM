# Bayesian Network SDM function
A function at applies a bayesian network to SDMs. Based on Staniczenko et al. (2017), [doi:10.1111/ele.12770](https://doi.org/10.1111/ele.12770)
# Function inputs:
* in_dir    = character string of directory of folder ending in "/" containing .tif input raster SDMs of all species e.g. "folder/"
* out_dir   = character string ending in "/" for output directory e.g. "BN_out/"
* focal     = character string file name of focal species inside in_dir folder including file extension e.g. "Gonimbrasia_belina.tif"
* direction = a vector of 1s and -1s that denote the direction of interaction of interacting species on focal species in alphabetical order e.g. c(1, -1, 1)
* method    = character string "or" (default) or "and" specifies state table filling method OR vs AND (box 2 of Staniczenko et al)

SDM rasters are loaded and used alphabetically, so the `direction` vector must be in the same order

### Source
`devtools::source_url("https://github.com/psijure/bayesian_network_sdm/blob/main/bnSDM_solve_func.R?raw=TRUE")`
