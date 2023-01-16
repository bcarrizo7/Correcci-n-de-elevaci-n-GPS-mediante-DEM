#install.packages(c("rgl", "rglwidget", "rayshader"))
library(rayshader)
library(rgl)
library(rglwidget)
library(sp)


three_plot_seg <- function (raster, segment , name){
  # Check if the raster is in the list of rasters with low resolution
  if (name %in% c("SRTM1arc", "NASADEM1arc")) {
    zscalenum <- 10
  } 
  
  else if(name %in% c("Alos Palsar 2")) {
    zscalenum <- 4
  } 
  
  else if(name %in% c( "TEG DCT" ,  "TEG SVD")) {
    zscalenum <- 2
  } 
  
  else {
    zscalenum <- 25
  }
  rgl.open()
  ll_prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  
  puntos_1 = sp::SpatialPointsDataFrame(coords = data.frame(cbind(segment$lon ,segment$lat )),
                                        data = data.frame(cbind(segment$lon ,segment$lat )),
                                        proj4string = sp::CRS(ll_prj))
  
  #veamos si calzan los features con las rows
  
  raster_extraction  = raster::extract(raster, puntos_1)
  
  e3 <- extent(min(segment$lon) - 0.01, # xmin
               max(segment$lon) + 0.01, # xmax
               min(segment$lat) - 0.01, # ymin
               max(segment$lat) + 0.01) # ymax
  
  raster_crop <- crop(raster, e3)
  
  #And convert it to a matrix:
  elmat = raster_to_matrix(raster_crop)
  
  #rgl.open()
  elmat %>%
    sphere_shade(sunangle = 60,texture = "imhof2") %>%
    plot_3d(elmat, zscale = zscalenum, fov = 0, theta = 135, zoom = 0.75, phi = 45) 
  
  render_points(extent =  attr(raster_crop,"extent"), 
                lat = segment$lat, long = segment$lon, 
                altitude = raster_extraction , zscale=zscalenum, color="blue",size = 10)
  
}



