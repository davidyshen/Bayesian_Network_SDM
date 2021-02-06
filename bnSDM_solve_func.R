# David Shen
# 06/02/2021
# bnSDM_solve_func.R
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
#   direction = vector of 1s (positive interaction) and -1s (negative interactions) that 
#               denote the direction of interaction of interacting species on focal species 
#               in alphabetical order 
#               e.g. c(1, -1, 1)
#   method    = character string "or" (default) or "and"
#               specifies state table filling method OR vs AND (box 2 of Staniczenko et al)
# 
# SDM rasters are loaded and used alphabetically so direction vector must be in the same 
# order

# Required packages ----
require(gRain)
# Some dependencies of gRain (RBGL) are no-longer hosted in CRAN, however can be 
# obtained through Bioconductor
require(raster)
require(dplyr)

# Function ----
bnSDM <- function(in_dir, out_dir = "BN_out/", focal, direction, method = "or")
{
  if(!dir.exists(out_dir)){
    dir.create(out_dir, recursive = T)
  }
  pb <- txtProgressBar(min = 1, max = ncell(stack), style = 3)
  levels = c("pres", "abs")
  files <- list.files(in_dir)
  interactors <- files[-which(files == focal)]
  
  # A stack of rasters of the non-focal species
  stack <- raster::stack(paste0(in_dir, interactors))
  cat("Extracting values from interacting species... \n")
  stackValues <- raster::values(stack)
  cat("Done \n")
  
  # Raster of focal species
  focalSp <- raster::raster(paste(in_dir, focal, sep = "/"))
  cat("Extracting values from focal species... \n")
  focalValues <- raster::values(focalSp)
  cat("Done \n")
  
  # Names of non-focal species
  names <- names(stack)
  # Name of focal species
  focalname <- names(focalSp)
  # Output focal species raster
  out <- focalSp
  outvals <- c()
  
  cat("Calculating posterior values for each cell... \n")
  # Working cell by cell of raster
  for (i in 1:ncell(stack))
  {
    # Run if focal species presence != 0
    if(focalValues[i] != 0)
    {
      # Make table for each interacting species
      for (j in 1:length(names))
      {
        assign(names[j], gRain::cptable(
          as.formula(paste0("~", names[j])),
          levels = levels,
          # Get probability values at cell
          values = c(stackValues[i, j], 1 - stackValues[i, j])
        ))
      }
      # Calculate focal species probability
      f <- matrix(nrow = 2^length(names), ncol = length(names))
      for (k in 1:nrow(f)) {f[k,] <- as.numeric(intToBits(k-1))[1:length(names)]}
      f <- as.data.frame(f[nrow(f):1,])
      colnames(f) <- names
      for (l in 1:ncol(f)) {f[,l] <- f[,l] * direction[l]}
      f$focal <- apply(f, 1, sum)
      
      if(method == "or")
      {
        for(m in 1:nrow(f))
        {
          # If net = 0, no change
          if(f$focal[m] == 0)
          {
            f$p[m] <- focalValues[i]
          }
          # If net > 0, increased probability
          else if(f$focal[m] > 0)
          {
            f$p[m] <- focalValues[i] + min(focalValues[i], 1-focalValues[i])
          }
          # If net < 0, decreased probability
          else if(f$focal[m] < 0)
          {
            f$p[m] <- focalValues[i] - min(focalValues[i], 1-focalValues[i])
          }
        }
      }
      if(method == "and")
      {
        for(m in 1:nrow(f))
        {
          # If all positive interactors are present
          if(f$focal[m] == table(direction)[2])
          {
            f$p[m] <- focalValues[i] + min(focalValues[i], 1-focalValues[i])
          }
          # If all negative interactors are present
          else if(f$focal[m] == table(direction)[1]*-1)
          {
            f$p[m] <- focalValues[i] - min(focalValues[i], 1-focalValues[i])
          }
          # Else no chance
          else {f$p[m] <- focalValues[i]}
        }
      }
      f$a <- 1-f$p
      focalpa <- c(t(f %>% dplyr::select(p, a)))
      
      # Complete table for focal species
      assign(focalname, gRain::cptable(
        as.formula(paste0("~", focalname, "|", paste(names, collapse = "+"))),
        levels = levels,
        values = focalpa)
      )
      
      # Solve bayesian network
      comp <- mget(c(names, focalname))
      plist <- compileCPT(comp)
      BN <- grain(plist)
      posteriors <- querygrain(BN)
      outvals[i] <-  eval(parse(text = paste0("posteriors$", focalname)))[1]
    }
    # Skip above and focal species presence = 0 if previously 0
    if(focalValues[i] == 0) {outvals[i] <- 0}
    
    setTxtProgressBar(pb, i)
  }
  cat("\n Done \n")
  # Write posterior values to output raster
  raster::values(out) <- outvals
  cat("Writing output to: ", paste0("~", out_dir, focalname, "_bnSDM.tif"), "... \n")
  suppressWarnings(raster::writeRaster(out, paste0(out_dir, focalname, "_bnSDM.tif"), format = "GTiff", overwrite = T))
}
