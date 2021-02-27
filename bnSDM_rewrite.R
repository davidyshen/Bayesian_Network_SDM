# David Shen
# 27/02/2021
# bnSDM_rewrite.R v5
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
  if(!dir.exists(out_dir)){
    dir.create(out_dir, recursive = T)
  }
  files <- list.files(in_dir)
  interactors <- files[-which(files == focal)]
  
  # Multicore processing
  if(ncores == "auto"){
    cores <- parallel::detectCores()
  } else if(is.numeric(ncores)) {
    if(ncores > parallel::detectCores()){stop("More cores specified than available: ", parallel::detectCores(), " usable")}
    cores <- ncores
  }
  
  
  # A stack of rasters of species with focal species last
  stack <- raster::stack(paste0(in_dir, interactors), paste(in_dir, focal, sep = "/"))
  cat("Extracting values... \n")
  values <- raster::values(stack)
  
  print(head(values))
  
  cat("Done \n")
  
  out <- raster::raster(nrows = nrow(stack), ncols = ncol(stack), ext = raster::extent(stack), crs = raster::crs(stack))
  
  cat("Calculating posterior values for each cell... \n")
  # Working cell by cell of raster
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl, c("values", "direction", "method", "stateTable", "focalP", "orFunc", "andFunc"))
  outvals <- parallel::parRapply(cl, values, interP, direction, method)
  parallel::stopCluster(cl)
  
  cat("\n Done \n")
  # Write posterior values to output raster
  raster::values(out) <- outvals
  cat("Writing output to: ", paste0("~", out_dir, names(out), "_bnSDM.tif"), "... \n")
  suppressWarnings(raster::writeRaster(out, paste0(out_dir, names(out), "_bnSDM.tif"), format = "GTiff", overwrite = T))
}

# Dependent Functions ----
## Function that fills state table for interacting species ----
interP <- function(inter, direction, method) {
  t0 <- stateTable(direction)
  t0 <- t(t(t0)*inter[-length(inter)])
  t1 <- stateTable(direction)
  t1 <- abs(t(t(1-t1) * inter[-length(inter)])-1)
  t2 <- abs(stateTable(direction)-1)
  t <- t0+t1*t2
  t <- cbind(t, focalP(fp = inter[length(inter)], direction, method))
  return(sum(apply(t, 1, prod)))
}

## Function that solves state table for focal species ----
focalP <- function(fp, direction, method) {
  m <- stateTable(direction)
  m <- t(t(m)*direction)
  d <- apply(m, 1, sum)
  if(method == "or") {
    return(sapply(d, orFunc, fp))
  } else if(method == "and") {
    return(sapply(d, andFunc, fp, direction))
  }
}

## Function that makes state tables from direction vector ----
stateTable <- function(direction) {
  m1 <- matrix(nrow = 2^length(direction), ncol = length(direction))
  for (k in 1:nrow(m1)) {m1[k,] <- as.numeric(intToBits(k-1))[1:length(direction)]}
  m1 <- m1[nrow(m1):1,]
  return(m1)
}

## Function that evaluates state using OR ----
orFunc <- function(x, fp) {
  if(x == 0) {
    return(fp)
  } else if(x > 0) {
    return(fp + min(fp, 1-fp))
  } else if(x < 0) {
    return(fp - min(fp, 1-fp))
  }
}

## Function that evaluates state using AND ----
andFunc <- function(x, fp, direction) {
  dt <- table(direction)
  if (x == dt[2]) {
    return(fp + min(fp, 1-fp))
  } else if (x == dt[1]*-1) {
    return(fp - min(fp, 1-fp))
  } else {return(fp)}
}
