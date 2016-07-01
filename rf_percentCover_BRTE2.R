#############################################################################
# This script reads training data from the CSV file created using the "percentCoverResample.R 
# script.  The script then uses the X and Y coordinates from the training data file to select
# the pixel values (predictor values) for each sample point in the input image. The predictor 
# values and the percent cover data from the training data file (response variable) are 
# combined and used as input to the random forests model. After the model is created percent 
# cover predictions are made on the input image to create an output image with percent cover 
# values ranging from 0 to 1. 
# 
# Set the variables below in the "SET VARIABLES HERE" section of the script. 
#
# This script was written by Ned Horning [horning@amnh.org]
# Support for writing and maintaining this script comes from The John D. and 
# Catherine T. MacArthur Foundation.
#
# This script is free software; you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software Foundation
# either version 2 of the License, or ( at your option ) any later version.                                   
#
#############################################################################
#Load libraries
require(maptools)
require(sp)
require(randomForest)
require(raster)
require(rgdal)

#
#############################   SET VARIABLES HERE  ###################################
#############################   point data     #######################################
# setwd("C:/Users/Bethany/Documents/My Dropbox/Bethany/Research/Invasive_mapping/Cheatgrass_RS_2015")
setwd("E:/artr_obl_spp")
# The CSV file containing X, Y, and percent cover point data created by the percentCoverResample.R script.
artr_dat <- read.csv('artr_pts/artr_clip2gb.csv', header=TRUE)
long <-artr_dat$long_proj
lat <-artr_dat$lat_proj
pctcov <-artr_dat$pctcov

#create pointdata from csv
pointData <- cbind(long, lat, pctcov)




#############################   raster data - choose one option  #########################
#############################   option 1: bioclm data with mask   #######################
#load raster layers (in this case, 19 bioclim layers) and turn them into a single raster stack.
files <- list.files(("_pvs"), pattern = 'tif$', full.names=TRUE)
# files <-list.files(("clim/Daymet_monthly/dymt_ppt20inc_bioclm_rmvcorrelated/TIFS"), pattern = 'tif$', full.names=TRUE)

bioclm <- stack(files)
plot(bioclm)

#import raster mask
r_mask<-raster("sb_outlineNEW/sbrush_mask.tif")
plot(r_mask)

#apply the mask to the raster layers.
#set maskvalue = to whatever value you defined as the NA value for the tif
x <- mask(bioclm, r_mask, maskvalue=0)
plot(x)
#write the masked file to save it
DaymetBioclim <- writeRaster(x, "Daymet_monthly/dymt_curr_bioclm/dymt_bioclm_tifs/19bclm_masked.tif", format='GTiff')

############################## option 2: bioclim data no mask   ##########################
# create raster stack of PV layers
files <-list.files(("Daymet_monthly/dymt_curr_bioclm/dymt_bioclm_tifs"), pattern = 'tif$', full.names=TRUE)
predictors <- stack(files)
# predictors <- dropLayer(x=predictors, i=1, 3)  #remove any unwated layers
predictors
s <- stack(predictors)


############################ option 3: single raster or satellite layer   ################
# Name and path for the input satellite image 
inImage <-'clim/stacks/annual_gb.bil'


############################    define output   ###########################################
# Name and path of the output GeoTiff image
outImage <- 'rf_output/artr_23PVs2.tif'

# No data value for satellite image or raster stack
# nd <- -9999
nd <- -1.7e+308

######################################################################################
#
# Start processing
print("Set variables and start processing")
startTime <- Sys.time()
cat("Start time", format(startTime),"\n")

# pointTable <- read.csv(pointData, header=TRUE)
xy <- SpatialPoints(pointData[,1:2])
response <- as.numeric(pointData[,3])

# Load the moderate resolution image
satImage <- stack(x)

#don't use if using stacked raster ('s')
# for (b in 1:nlayers(satImage)) { NAvalue(satImage@layers[[b]]) <- nd }

# Get pixel DNs from the input image for each sample point
print("Getting the pixel values under each point")
trainvals <- cbind(response, extract(satImage, xy)) 

# Remove NA values from trainvals
trainvals_no_na <- na.omit(trainvals)



#####################################################################

##see http://www.statistik.uni-dortmund.de/useR-2008/slides/Strobl+Zeileis.pdf for other options
## also, varSelRF for backward elimination



library (varSelRF)
rf.vs1 <- varSelRF(x, cl, ntree = 200, ntreeIterat = 100,
                   vars.drop.frac = 0.2)
rf <- randomForest(x, cl, ntree = 200, importance = TRUE)
rf.rvi <- randomVarImpsRF(x, cl, 
                          rf, 
                          numrandom = 20, 
                          usingCluster = FALSE) 

randomVarImpsRFplot(rf.rvi, rf)










#OR











##add rfcv and tuneRF code here - figure out which PVs to drop
head(trainvals)
alldat <- read.csv("_pvs/pvs_extractto_artrGB.csv", head=TRUE)
pvs <- cbind(alldat[10:32])
response <- cbind(alldat[5])

result <-rfcv(pvs, alldat$pctcov, cv.fold=5)
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))


#############################################################################
# Run Random Forest
print("Starting to calculate random forest object")
randfor <- randomForest(response ~. , data=trainvals_no_na, importance=TRUE)

#save randfor model
# save(randfor, file = "rf_output/randfor_model_curr.rda") #only needed if projecting onto different layers than what the model was trained on

#if projecting onto "future" climate, change satImage to future layers now. 
# Start predictions
print("Starting predictions")
predict(satImage, randfor, filename=outImage, progress='text', format='GTiff', datatype='FLT4S', type='response', overwrite=TRUE)
#
# Calculate processing time
timeDiff <- Sys.time() - startTime
cat("Processing time", format(timeDiff), "\n")

#Variable importance plotted in decreasing order (most important at bottom)
varImpPlot(randfor, sort=TRUE)
importance(randfor) #increase in node impurity = 'residual sum of squares'

#Plot error rates vs. number of trees
plot(randfor)

#Plot response curves of individual predictor variables to regression
partialPlot(randfor,trainvals_no_na, bio1_mskd) #Change BRTE_predictors.X accordingly


