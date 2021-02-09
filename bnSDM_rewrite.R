bnSDM2 <- function(in_dir, out_dir = "BN_out/", focal, direction, method = "or")
{
  if(!dir.exists(out_dir)){
    dir.create(out_dir, recursive = T)
  }
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
  
  values <- cbind(stackValues, focalValues)
  
  out <- focalSp
  
  cat("Calculating posterior values for each cell... \n")
  # Working cell by cell of raster
  #### TO BE REWRITTEN ####
  outvals <- apply(values, 1, interP, direction, method)
  
  cat("\n Done \n")
  # Write posterior values to output raster
  raster::values(out) <- outvals
  cat("Writing output to: ", paste0("~", out_dir, focalname, "_bnSDM.tif"), "... \n")
  suppressWarnings(raster::writeRaster(out, paste0(out_dir, focalname, "_bnSDM.tif"), format = "GTiff", overwrite = T))
}

# Function that fills state table for interacting species ----
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

testfunc <- function(inter) {
  
}

# Function that solves state table for focal species ----
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

# Function that makes state tables from direction vector ----
stateTable <- function(direction) {
  m1 <- matrix(nrow = 2^length(direction), ncol = length(direction))
  for (k in 1:nrow(m1)) {m1[k,] <- as.numeric(intToBits(k-1))[1:length(direction)]}
  m1 <- m1[nrow(m1):1,]
  return(m1)
}

# Function that evaluates state using OR ----
orFunc <- function(x, fp) {
  if(x == 0) {
    return(fp)
  } else if(x > 0) {
    return(fp + min(fp, 1-fp))
  } else if(x < 0) {
    return(fp - min(fp, 1-fp))
  }
}

# Function that evaluates state using AND ----
andFunc <- function(x, fp, direction) {
  dt <- table(direction)
  if (x == dt[2]) {
    return(fp + min(fp, 1-fp))
  } else if (x == dt[1]*-1) {
    return(fp - min(fp, 1-fp))
  } else {return(fp)}
}
