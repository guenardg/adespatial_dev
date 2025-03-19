
# rm(list=ls())
if(FALSE) {
  try(dyn.unload("src/geoddist.so"))
  system("R CMD SHLIB src/*.f")
  dyn.load("src/geoddist.so")
  
  source("R/geoDist.R")
  source("R/geoXY.R")
}

if(FALSE) {
  load("gpsObject1.rda")
  data.frame(
    latitude = gpsObject1@latitude,
    longitude = gpsObject1@longitude,
    elevation = gpsObject1@elevation,
    time = gpsObject1@time
  ) -> gpsObject1
  ## class(gpsObject1$time)
  save(gpsObject1, file="inst/extdata/gpsObject1.rda")
}
## load("inst/extdata/gpsObject1.rda")

library(adespatial)

load(system.file("extdata/gpsObject1.rda", package = "adespatial"))

xy <- geoXY(gpsObject1$latitude, gpsObject1$longitude, unit = 1000)

plot(xy[,1], xy[,2], asp = 1)
