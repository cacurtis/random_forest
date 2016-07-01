################################# GETTING SET UP ##########################################
#Load libraries
require(maptools)
require(sp)
require(randomForest)
require(raster)
require(rgdal)
#set your working directory
setwd("F:/artr_obl_spp")




################################ INPUTS AND OUTPUTS ########################################

#load the CSV file containing X, Y, and percent cover point data.  
#Then, assign your lat, long, and response variable (e.g. percent cover) and bind them into a new dataset
artr_dat <- read.csv('artr_pts/artr_clip2gb.csv', header=TRUE)
long <-artr_dat$long_proj
lat <-artr_dat$lat_proj
pctcov <-artr_dat$pctcov
pointData <- cbind(long, lat, pctcov)
xy <- SpatialPoints(pointData[,1:2])
response <- as.numeric(pointData[,3])

# Name and path of the  GeoTiff image output from the Random Forest
outImage <- 'rf_output/artr_clim_9pvs.tif'




#############################   RASTER DATA - CHOOSE ONE OPTION  #########################



#OPTION ONE: MASKED RASTER STACK 
#replace "_pvs" with the name of the folder holding your rasters
files <- list.files(("_pvs/selected"), pattern = 'tif$', full.names=TRUE)
rasterstack <- stack(files)
#visualize the raster files you loaded.  this will take a while for large/many rasters, so feel free to skip it.
plot(rasterstack)
#import raster mask
r_mask<-raster("sb_outlineNEW/sbrush_mask.tif")
#visualize the mask (if you want to).
plot(r_mask)
#apply the mask to the raster layers.
#set maskvalue = to whatever  you defined as the NA value for the mask 
  #(i.e., keep the areas that are "1" and maks out areas that are "0")
x <- mask(rasterstack, r_mask, maskvalue=0)
#visualize your new, masked raster layers (if you want to)
plot(x)
#write the masked file to save it.  change the directory to the loction you want to save the file to. 
DaymetBioclim <- writeRaster(x, "Daymet_monthly/dymt_curr_bioclm/dymt_bioclm_tifs/19bclm_masked.tif", format='GTiff')
#define the NoData value from your raster stack or image.
#If you have predictor variables from multiple sources, you might need to change their NoData values to be consistent.
nd <- -1.7e+308


#OPTION TWO:  RASTER STACK - NO MASK
files <-list.files(("folder_holding_rasters"), pattern = 'tif$', full.names=TRUE)
predictors <- stack(files)
predictors <- dropLayer(x=predictors, i=1, 3)  #remove any unwated layers
predictors
s <- stack(predictors)
#define the NoData value from your raster stack or image.
#If you have predictor variables from multiple sources, you might need to change their NoData values to be consistent.
nd <- -1.7e+308


#OPTION THREE: SINGLE RASTER
# Name and path for the input image 
inImage <-'directory/file_name.bil'
#define the NoData value from your raster stack or image.
#If you have predictor variables from multiple sources, you might need to change their NoData values to be consistent.
nd <- -1.7e+308



################################# PREP AND RUN THE RF #######################################
# Get pixel DNs from the input image for each sample point
satImage <- stack(x)
trainvals <- cbind(response, extract(satImage, xy)) 
trainvals_no_na <- na.omit(trainvals)

# Run Random Forest
randfor <- randomForest(response ~. , data=trainvals_no_na, importance=TRUE)

#if projecting onto "future" climate (i.e., different layers than were used to train the model), use the following 
#NOTE: your "future" cliamte layers will need to have exactly the same names as those used to train the models
files <- list.files(("folder_holding_future_layers"), pattern = 'tif$', full.names=TRUE)
rasterstack <- stack(files)
x <- mask(rasterstack, r_mask, maskvalue=0)
satImage <- stack(x)


#use the RF model defined above to predict land cover for current or future layers
predict(satImage, randfor, filename=outImage, progress='text', format='GTiff', datatype='FLT4S', type='response', overwrite=TRUE)

#Variable importance plotted in decreasing order (most important at bottom)
varImpPlot(randfor, sort=TRUE)
importance(randfor) #increase in node impurity = 'residual sum of squares'
#Plot error rates vs. number of trees
plot(randfor)
#Plot response curves of individual predictor variables to regression
partialPlot(randfor,trainvals_no_na, name_of_predictorvariable) 





