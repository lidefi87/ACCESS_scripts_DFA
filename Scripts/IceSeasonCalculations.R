#Script taken from 'aceecostats' package developed by AAD


# Load libraries ----------------------------------------------------------
library(raster)


# Define calc_ice_season function --------------------------------------
#Inputs:
#yfile takes a filepath to a .grd file
#threshval refers to the min ice concentration value that will be considered for ice seasonality calculations
calc_ice_season <- function(yfile, threshval = 15) {
  #Extract what hemisphere the data represents
  hemi <- substr(basename(yfile[1]), 1, 5)
  threshold.value <- threshval
  #Minimum number of days a pixel has 15% or more ice concentration
  ndays <- 5
  
  ## north_ice_1979_02_15_1980_02_14.grd"
  #Creates a multilayer raster file using the files in the file path
  ice <- brick(yfile)
  #Number of layers in raster brick
  year_n <- nlayers(ice)
  #Create a new raster layer with 0 using the already existing grid
  template <- ice[[1]] * 0
  #Extracting all values of the raster into a matrix. Each column contains all values for one layer
  icemat <- values(ice)
  #Replace any values above 100 with 0
  icemat[icemat > 100] <- 0
  
  ## here we need to get the next day for the interpolation . . .
  ##icemat[is.na(icemat)] <- 0
  #Creating vectors with 0 using the same number of rows in icemat for advance and retreat of ice
  adv <- numeric(nrow(icemat))
  ret <- numeric(nrow(icemat))
  #Create a boolean matrix. True will include all values more or equal to the threshold value
  threshold <- icemat >= threshold.value
  
  #Adding up all values per row
  rsum <- .rowSums(threshold, nrow(threshold), ncol(threshold))
  ## If rsum is 0, it means that that pixel did not have ice concentrations >15% during the entire year
  alllt <- rsum == 0
  ## If rsum is the same at the number of columns in the matrix, it means that ice concentration of 15% or above was recorded all year
  allgt <- rsum == ncol(threshold)
  ## values missing
  ##miss <- icemat[,1] > 100
  ## Identifies cells with ice concentration above 15% at least once
  visit <- which(!alllt & !allgt)
  for (ii in seq_along(visit)) {
    #Identifies unique values (values) and how many times these values were repeated in a row (lengths)
    rl <- rle(threshold[visit[ii], ])
    
    ##    annual day of advance is the time when the ice
    ##    concentration in a given pixel first exceeds 15% (taken to
    ##    approximate the ice edge) for at least 5 days
    for (ir in seq_along(rl$lengths)) {
      #If value is TRUE and its length is at least 5
      if (rl$values[ir] & rl$lengths[ir] >= ndays) {
        #Get the location of 
        adv[visit[ii]] <- if (ir == 1) 1L else sum(head(rl$lengths, ir - 1))
        break;
      }
    }
    #If value is 0, then change to NA
    if (adv[visit[ii]] == 0) adv[visit[ii]] <- NA
    
    ##while day of retreat
    ## ## is the time when concentration remains below 15% until the end
    ## ## of the given sea ice year
    revlengths <- rev(rl$lengths)
    revvals <- rev(rl$values)
    for (ri in seq_along(revlengths)) {
      if (revvals[ri]) {
        ret[visit[ii]] <- if (ri == 1) length(year_n) else sum(revlengths[ri:length(revlengths)])
        break;
      }
    }
    if (ret[visit[ii]] == 0) ret[visit[ii]] <- NA
  }
  
  
  #Pixels with no ice concentration above 15% set to NA
  adv[alllt] <- NA
  ## adv[miss] <- NA
  #Pixels where ice was present all year long set to 1
  adv[allgt] <- 1
  #Pixels with no ice concentration above 15% set to NA
  ret[alllt] <- NA
  #Pixels where ice was present all year long set to 1
  ret[allgt] <- length(year_n)
  
  #Assign values to template
  list(adv = setValues(template, adv), ret = setValues(template, ret))
}

yfile <- "../Downloads/south_ice_2014_02_15_2015_02_15.grd"

x <- calc_ice_season(yfile)
