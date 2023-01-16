#install.packages("raster")
library(raster)


#descargamos los rasters
Srtm_1_arc <- raster("DEMs _DF_RM/SRTM_s34_w071_1arc_v3.tif")
Srtm_void_filled_3_arc <-  raster("DEMs _DF_RM/SRTM_void_filled_s34_w071_3arc_v2.tif")
Alos_lon_lat <- raster("DEMs _DF_RM/Alos_lon_lat.tif")
nasadem_1_arc <- raster("DEMs _DF_RM/s34w071.hgt")
tanDEM <-  raster("DEMs _DF_RM/TDM1_DEM__30_S34W071_DEM.tif")
Merit_DEM <-  raster("DEMs _DF_RM/s35w075_dem.tif")

