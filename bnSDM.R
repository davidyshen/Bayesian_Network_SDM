# David Shen
# bnSDM_rewrite.R
# A function that applies post-hoc bayesian networks to a species distribution
# model
# Based on Staniczenko et al. 2017    doi:10.111/ele.12770
# ----------------------------------------------------------------------------------- #
# Help ----
# Function inputs:
#   in_dir    = character string of directory of folder ending in "/" containing .tif input 
#               raster SDMs of all species 
#               e.g. "folder/"
#   out_dir   = character string ending in "/" for output directory 
#               e.g. "BN_out/"
#   focal     = character string file name of focal species inside in_dir folder including
#               file extension 
#               e.g. "Gonimbrasia_belina.tif"
#   direction = a vector of 1s and -1s that denote the direction of interaction of interacting
#               species on focal species in alphabetical order 
#               e.g. c(1, -1, 1)
#   method    = character string "or" (default) or "and"
#               specifies state table filling method OR vs AND (box 2 of Staniczenko et al)
# 
# SDM rasters will be loaded and used alphabetically

# Function ----
bnSDM <- function(in_dir, 
                  out_dir = "BN_out/", 
                  focal, 
                  direction, 
                  method = "or", 
                  ncores = "auto")
{
  # # Export variables to global env because the cluster doesn't like it normally???
  # direction <- direction
  # method <- method
  
  # If output directory doesn't exist, make it
  if(!dir.exists(out_dir)){
    dir.create(out_dir, recursive = T)
  }
  files <- list.files(in_dir)
  # interactors are the files excluding the focal species
  # Updated to allow partial string matching instead of verbatim filenames
  interactors <- files[!stringr::str_detect(files,focal)]
  
  # Multicore processing
  if(ncores == "auto"){
    cores <- parallel::detectCores()
  } else if(is.numeric(ncores)) {
    if(ncores > parallel::detectCores()){stop("More cores specified than available: ", parallel::detectCores(), " usable")}
    cores <- ncores
  }
  
  
  # A stack of rasters of species with focal species last
  stack <- raster::stack(file.path(in_dir, interactors), file.path(in_dir, focal))
  cat("Extracting values... \n")
  # Extract the values of the stack, focal species last column
  value <- raster::values(stack)
  
  # ## New stuff ##
  # # Set values below threshold to NA 
  # value <- matrix(nrow = nrow(tempvalue), ncol = ncol(tempvalue))
  # for (c in 1:length(threshold)){
  #   value[,c] <- unlist(sapply(tempvalue[,c], function(q, thr){
  #     if(is.na(q) == T){
  #       NA
  #     } else if(q < thr){
  #       q <- NA
  #     } else {q}}, 
  #     thr = threshold[c]))
  # }
  # 
  # sp_pres <- apply(value, 1, function(x) {length(na.omit(x))})
  # hist(sp_pres)
  
  cat("Done \n")
  
  # Make state table once and reuse for each cell rather than making a new table every cell
  st <- .stateTable(direction)
  
  # Make empty raster for posterior occurrence of focal species
  out <- raster::raster(nrows = nrow(stack), ncols = ncol(stack), ext = raster::extent(stack), crs = raster::crs(stack))
  
  cat("Calculating posterior values for each cell... \n")
  # Working cell by cell of raster
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl, c("value", "direction", "method", "st", ".interP", ".focalP", ".orFunc", ".andFunc"), envir = environment())
  outvals <- parallel::parRapply(cl, value, .interP, direction, method, st)
  parallel::stopCluster(cl)
  
  cat("\n Done \n")
  # Write posterior values to output raster
  raster::values(out) <- outvals
  cat("Writing output to: ", paste0("~", out_dir, names(out), "_bnSDM.tif"), "... \n")
  suppressWarnings(raster::writeRaster(out, file.path(out_dir, paste0(names(out), "_bnSDM.tif")), format = "GTiff", overwrite = T))
}

# Dependent Functions ----
## Function that fills state table for interacting species ----
# inter = vector of occurrence probabilities at a single point
# st = blank state table with nrow = number of species
.interP <- function(inter, direction, method, st) {
  # Generate a state table for the number of interacting species; rows = 2^n_species, cols = n_species
  t0 <- st
  # Multiply the each column of state table by the occurrence probability of each species
  t0 <- t(t(t0)*inter[-length(inter)])
  t1 <- st
  # Make state table where if species is absent, value = 1-(probability of occurrence)
  t1 <- abs(t(t(1-t1) * inter[-length(inter)])-1)
  t2 <- abs(st-1)
  # Combine state tables: 1|>probability of occurence, 0|>1-probability of occurrence
  t <- t0+t1*t2
  # Bind final column with conditional probability of occurrence of focal species
  t <- cbind(t, .focalP(fp = inter[length(inter)], direction, method, st))
  # Return the product of each row all summed as the posterior occurrence of the focal species
  return(sum(apply(t, 1, prod)))
}

## Function that solves state table for focal species ----
.focalP <- function(fp, direction, method, st) {
  # Make a state table
  m <- st
  # Multiply each column by the direction of interaction of each interacting species
  m <- t(t(m)*direction)
  # Calculate the cumulative impact
  d <- apply(m, 1, sum)
  if(method == "or") {
    # If using OR, calculate probability of occurrence of focal species based on OR method
    return(sapply(d, .orFunc, fp))
  } else if(method == "and") {
    # If using AND, calculate probability of occurrence of focal species based on AND method
    return(sapply(d, .andFunc, fp, direction))
  }
}

## Function that makes state tables from direction vector ----
.stateTable <- function(direction) {
  m1 <- matrix(nrow = 2^length(direction), ncol = length(direction))
  for (k in 1:nrow(m1)) {m1[k,] <- as.numeric(intToBits(k-1))[1:length(direction)]}
  m1 <- m1[nrow(m1):1,]
  return(m1)
}

## Function that evaluates state using OR ----
.orFunc <- function(x, fp) {
  if(x == 0) {
    # If cumulative direction of interaction is 0, no change
    return(fp)
  } else if(x > 0) {
    # If cumulative direction of interaction > 0, raise occurrence prob (P) by P+min(P, 1-P)
    return(fp + min(fp, 1-fp))
  } else if(x < 0) {
    # If cumulative direction of interaction < 0, lower occurrence prob by P-min(P, 1-P)
    return(fp - min(fp, 1-fp))
  }
}

## Function that evaluates state using AND ----
.andFunc <- function(x, fp, direction) {
  dt <- table(direction)
  if (x == dt[2]) {
    # If only positive interactors are present, increase prob by P+min(P, 1-P)
    return(fp + min(fp, 1-fp))
  } else if (x == dt[1]*-1) {
    # If only negative interactors are present, decrease prob by P-min(P,1-P)
    return(fp - min(fp, 1-fp))
    # Otherwise, no change
  } else {return(fp)}
}

# Other functions ----
## Function to combine multiple rasters into a single raster ----
# in_dir    = directory containing rasters to be smushed into a single/reduced number of rasters, must end in /
# out_dir   = directory to output in, must end in /
# out_name  = name of output file(s)
# n         = number of be reduced to, defaults to 1
smush <- function(in_dir, out_dir, out_name, n = 1, ncores = "auto") {
  files <- list.files(in_dir)
  
  if(!dir.exists(out_dir)){
    dir.create(out_dir, recursive = T)
  }
  
  stack <- raster::stack(as.list(file.path(in_dir, files)))
  outRas <- raster::raster(nrows = nrow(stack), 
                           ncols = ncol(stack), 
                           ext = raster::extent(stack), 
                           crs = raster::crs(stack))
  
  # Multicore processing
  if(ncores == "auto"){
    cores <- parallel::detectCores()
  } else if(is.numeric(ncores)) {
    if(ncores > parallel::detectCores()){stop("More cores specified than available: ", parallel::detectCores(), " usable")}
    cores <- ncores
  }
  
  values <- raster::values(stack)
  cl <- parallel::makeCluster(cores)
  maxVals <- parallel::parApply(cl, values, 1, function(x, n){
    sort(x, decreasing = T)[1:n]
  }, n)
  parallel::stopCluster(cl)
  
  if(n > 1) {
    for(i in 1:n){
      raster::values(outRas) <- maxVals[i,]
      raster::writeRaster(outRas, paste0(file.path(out_dir,out_name), "_", i, ".tif"), format = "GTiff", overwrite = T)
    }
  } else {
    raster::values(outRas) <- maxVals
    raster::writeRaster(outRas, paste0(file.path(out_dir,out_name),".tif"), format = "GTiff", overwrite = T)
  }
}
