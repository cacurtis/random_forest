#load packages and set working directory
library(raster)
library(sp)
library(maptools)
library(rgdal)
setwd("E:/artr_obl_spp")

#load raster layers (in this case, 19 bioclim layers) and turn them into a single raster stack.
files <-list.files(("pheno_metrics/TINDVI_2001-2014/TINDVI_2001-2014_clipped_TIFs/final_rasters2/tifs"), pattern = 'tif$', full.names=TRUE)
bioclm <- stack(files)
plot(bioclm)

#import raster mask
r_mask<-raster("sb_outlineNEW/sbrush_mask.tif")
plot(r_mask)

#apply the mask to the raster layers.
#set maskvalue = to whatever value you defined as the NA value for the tif
x <- mask(bioclm, r_mask, maskvalue=0)
plot(x)


if (require(rgdal)) {
  rf <- writeRaster(x, filename="pheno_metrics/TINDVI_2001-2014/masked.tif", format="GTiff", overwrite=TRUE)
#   bf <- writeRaster(b, filename="multi.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
}