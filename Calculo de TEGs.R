source("DCT & SVD Functions.R")
source("Rayshader_3D_plot.R")

library(elevatr)
library(dplyr)
library(sp)

TEG_DEMs <- function(segment, kx = 0.3, vl=0.5 ) {
  
  
  ll_prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  
  puntos = sp::SpatialPointsDataFrame(coords = data.frame(cbind(segment$lon ,segment$lat )),
                                      data = data.frame(cbind(segment$lon ,segment$lat )),
                                      proj4string = sp::CRS(ll_prj))
  
  #Define the extent of the crop by clicking on the plot
  cropbox1 <- extent(min(segment$lon) - 0.001, # xmin
                     max(segment$lon) + 0.01, # xmax
                     min(segment$lat) - 0.001, # ymin
                     max(segment$lat) + 0.001) # ymax
  #crop the raster, then plot the new cropped raster
  raster_elevatr_aws <- elevatr::get_elev_raster(locations = puntos, units = "meters",src="aws",z=14, override_size_check = TRUE) # src = c("aws", "gl3", "gl1", "alos", "srtm15plus") , se puede cambiar el dem 
  raster_elevatr_aws
  #plot all raster
  DEMcrop1 <- crop(raster_elevatr_aws, cropbox1)
  
  
  m <- as.matrix(DEMcrop1)
  head(m)
  
  rotate_clockwise  <- function(x) { t(apply(x, 2, rev))}
  m_90 <- rotate_clockwise(m) 
  
  volcano <- m_90 %>%
    reshape2::melt(value.name = "original") %>%
    transform(noisy = original + 1.5*rnorm(length(original))) %>% 
    transform(smooth_dct = smooth_dct(Var2, Var1, noisy, kx = kx)) %>%
    transform(smooth_svd = smooth_svd(Var2, Var1, noisy, variance_lost = vl)) %>%
    reshape2::melt(id.vars = c("Var1", "Var2"))
  

  
  C_l_original <- volcano %>% dplyr::filter(variable == "original") %>% dplyr::select(Var1,Var2,value)
  C_l_noizy <- volcano %>% dplyr::filter(variable =="noisy") %>% dplyr::select(Var1,Var2,value)
  C_l_smooth_dct <- volcano %>% dplyr::filter(variable == "smooth_dct") %>% dplyr::select(Var1,Var2,value)
  C_l_smooth_svd <- volcano %>% dplyr::filter(variable == "smooth_svd") %>% dplyr::select(Var1,Var2,value)
  
  
  C_l_original_raster <- rasterFromXYZ(C_l_original, crs =ll_prj)  
  C_l_noizy_raster <- rasterFromXYZ(C_l_noizy, crs =ll_prj)  
  C_l_smooth_dct_raster <- rasterFromXYZ(C_l_smooth_dct, crs =ll_prj)  
  C_l_smooth_svd_raster <- rasterFromXYZ(C_l_smooth_svd, crs =ll_prj)
  
  extent(C_l_original_raster) <- extent(DEMcrop1)
  extent(C_l_noizy_raster) <- extent(DEMcrop1)
  extent(C_l_smooth_dct_raster) <- extent(DEMcrop1)
  extent(C_l_smooth_svd_raster) <- extent(DEMcrop1)
  
  #three_plot_seg(C_l_smooth_dct_raster, segment , "C_l_smooth_dct_raster")
  #three_plot_seg(C_l_smooth_svd_raster, segment , "C_l_smooth_dct_raster")
  
  #extreaemos puntos de los dif DEMs
  Srtm_1_arc_extraction  <-  raster::extract(Srtm_1_arc, puntos)
  C_l_noizy_raster_extraction  <- raster::extract(C_l_noizy_raster, puntos)
  C_l_smooth_dct_raster_extraction   <-  raster::extract(C_l_smooth_dct_raster, puntos)
  C_l_smooth_svd_raster_extraction   <-  raster::extract(C_l_smooth_svd_raster, puntos)
  Srtm_void_filled_3_arc_extraction = raster::extract(Srtm_void_filled_3_arc, puntos)
  nasadem_1_arc_extraction   <-  raster::extract(nasadem_1_arc, puntos)
  Alos_lon_lat_extraction <-  raster::extract(Alos_lon_lat, puntos)
  tanDEM_extraction   <-  raster::extract(tanDEM, puntos)
  Merit_Dem_extraction <-  raster::extract(Merit_DEM, puntos)
  correct_elevation <- segment %>% ele_correction(replace = TRUE, z = 14)

  
  stream_raw_clean_dem_correct_ele_dct <- segment %>% mutate(ele = C_l_smooth_dct_raster_extraction  )
  stream_raw_clean_dem_correct_ele_svd  <- segment %>% mutate(ele = C_l_smooth_svd_raster_extraction )

  
  stream_clean_dem_correct_ele_smooth_srtm <- correct_elevation %>% smooth_stream(interpolate = FALSE, alpha = 0.07, replace = FALSE)
  stream_clean_dem_correct_ele_smooth_dct <- stream_raw_clean_dem_correct_ele_dct %>% smooth_stream(interpolate = FALSE, alpha = 0.07, replace = FALSE)
  stream_clean_dem_correct_ele_smooth_svd <- stream_raw_clean_dem_correct_ele_svd %>% smooth_stream(interpolate = FALSE, alpha = 0.07, replace = FALSE)

  
  real_dem_smooth <-  data.frame(Gps_data=segment$ele, Srtm1arc = Srtm_1_arc_extraction  , C_l_noizy_raster = C_l_noizy_raster_extraction, C_l_smooth_dct_raster = C_l_smooth_dct_raster_extraction,  C_l_smooth_svd_raster = C_l_smooth_svd_raster_extraction ,Gpstream_correction = correct_elevation$ele, Srtm3arcvoidfilled = Srtm_void_filled_3_arc_extraction, 
                                 Alos_palsar = Alos_lon_lat_extraction ,Nasadem1arc = nasadem_1_arc_extraction , TanDEM=tanDEM_extraction , MeritDEM = Merit_Dem_extraction 
                              )
  
  #Calculamos TEG
  
  real_dem_smooth$delta_ele_Gps_data<- c(0,diff(real_dem_smooth$Gps_data,lag=1))
  real_dem_smooth$dplus_Gps_data <- ifelse(real_dem_smooth$delta_ele_Gps_data>0, real_dem_smooth$delta_ele_Gps_data, 0)
  TEG_Gps_data<- sum(real_dem_smooth$dplus_Gps_data)
  TEG_Gps_data
  
  real_dem_smooth$delta_ele_Srtm1arc<- c(0,diff(real_dem_smooth$Srtm1arc,lag=1))
  real_dem_smooth$dplus_Srtm1arc <- ifelse(real_dem_smooth$delta_ele_Srtm1arc>0, real_dem_smooth$delta_ele_Srtm1arc, 0)
  TEG_Srtm1arc <- sum(real_dem_smooth$dplus_Srtm1arc)
  TEG_Srtm1arc
  
  real_dem_smooth$delta_ele_C_l_noizy_raster<- c(0,diff(real_dem_smooth$C_l_noizy_raster,lag=1))
  real_dem_smooth$dplus_C_l_noizy_raster <- ifelse(real_dem_smooth$delta_ele_C_l_noizy_raster>0, real_dem_smooth$delta_ele_C_l_noizy_raster, 0)
  TEG_C_l_noizy_raster <- sum(real_dem_smooth$dplus_C_l_noizy_raster)
  TEG_C_l_noizy_raster
  
  real_dem_smooth$delta_ele_C_l_smooth_dct_raster<- c(0,diff(real_dem_smooth$C_l_smooth_dct_raster,lag=1))
  real_dem_smooth$dplus_C_l_smooth_dct_raster <- ifelse(real_dem_smooth$delta_ele_C_l_smooth_dct_raster>0, real_dem_smooth$delta_ele_C_l_smooth_dct_raster, 0)
  TEG_C_l_smooth_dct_raster <- sum(real_dem_smooth$dplus_C_l_smooth_dct_raster)
  TEG_C_l_smooth_dct_raster
  
  real_dem_smooth$delta_ele_C_l_smooth_svd_raster<- c(0,diff(real_dem_smooth$C_l_smooth_svd_raster,lag=1))
  real_dem_smooth$dplus_C_l_smooth_svd_raster <- ifelse(real_dem_smooth$delta_ele_C_l_smooth_svd_raster>0, real_dem_smooth$delta_ele_C_l_smooth_svd_raster, 0)
  TEG_C_l_smooth_svd_raster<- sum(real_dem_smooth$dplus_C_l_smooth_svd_raster)
  TEG_C_l_smooth_svd_raster
  
  real_dem_smooth$delta_ele_Gpstream_correction <- c(0,diff(real_dem_smooth$Gpstream_correction,lag=1))
  real_dem_smooth$dplus_Gpstream_correction <- ifelse(real_dem_smooth$delta_ele_Gpstream_correction>0, real_dem_smooth$delta_ele_Gpstream_correction, 0)
  TEG_Gpstream_correction <- sum(real_dem_smooth$dplus_Gpstream_correction)
  TEG_Gpstream_correction
  
  
  real_dem_smooth$delta_ele_Srtm3arcvoidfilled <- c(0,diff(real_dem_smooth$Srtm3arcvoidfilled,lag=1))
  real_dem_smooth$dplus_Srtm3arcvoidfilled <- ifelse(real_dem_smooth$delta_ele_Srtm3arcvoidfilled >0, real_dem_smooth$delta_ele_Srtm3arcvoidfilled, 0)
  TEG_Srtm3voidfilled <- sum(real_dem_smooth$dplus_Srtm3arcvoidfilled)
  TEG_Srtm3voidfilled
  
  real_dem_smooth$delta_ele_Alos_palsar<- c(0,diff(real_dem_smooth$Alos_palsar,lag=1))
  real_dem_smooth$dplus_Alos_palsar <- ifelse(real_dem_smooth$delta_ele_Alos_palsar>0, real_dem_smooth$delta_ele_Alos_palsar, 0)
  TEG_Alos_palsar <- sum(real_dem_smooth$dplus_Alos_palsar)
  TEG_Alos_palsar
  
  real_dem_smooth$delta_ele_Nasadem1arc<- c(0,diff(real_dem_smooth$Nasadem1arc,lag=1))
  real_dem_smooth$dplus_Nasadem1arc <- ifelse(real_dem_smooth$delta_ele_Nasadem1arc>0, real_dem_smooth$delta_ele_Nasadem1arc, 0)
  TEG_Nasadem1arc <- sum(real_dem_smooth$dplus_Nasadem1arc)
  TEG_Nasadem1arc
  
  real_dem_smooth$delta_ele_TanDEM<- c(0,diff(real_dem_smooth$TanDEM,lag=1))
  real_dem_smooth$dplus_TanDEM <- ifelse(real_dem_smooth$delta_ele_TanDEM>0, real_dem_smooth$delta_ele_TanDEM, 0)
  TEG_TanDEM <- sum(real_dem_smooth$dplus_TanDEM)
  TEG_TanDEM
  
  real_dem_smooth$delta_ele_MeritDEM<- c(0,diff(real_dem_smooth$MeritDEM,lag=1))
  real_dem_smooth$dplus_MeritDEM <- ifelse(real_dem_smooth$delta_ele_MeritDEM>0, real_dem_smooth$delta_ele_MeritDEM, 0)
  TEG_MeritDEM <- sum(real_dem_smooth$dplus_MeritDEM)
  TEG_MeritDEM
  
  stream_raw_clean_dem_correct_ele_dct <- segment %>% mutate(ele = real_dem_smooth$C_l_smooth_dct_raster  )
  stream_raw_clean_dem_correct_ele_svd  <- segment %>% mutate(ele = real_dem_smooth$C_l_smooth_svd_raster )
  
  stream_clean_dem_correct_ele_smooth_srtm <- correct_elevation %>% smooth_stream(interpolate = FALSE, alpha = 0.07, replace = TRUE)
  stream_clean_dem_correct_ele_smooth_dct <- stream_raw_clean_dem_correct_ele_dct %>% smooth_stream(interpolate = FALSE, alpha = 0.07, replace = TRUE)
  stream_clean_dem_correct_ele_smooth_svd <- stream_raw_clean_dem_correct_ele_svd %>% smooth_stream(interpolate = FALSE, alpha = 0.07, replace = TRUE)
  
  dif_dem_correct_ele_smooth_srtm <- stream_clean_dem_correct_ele_smooth_srtm %>% differential_stream()
  dif_dem_correct_ele_smooth_dct <- stream_clean_dem_correct_ele_smooth_dct %>% differential_stream()
  dif_dem_correct_ele_smooth_svd <- stream_clean_dem_correct_ele_smooth_svd  %>% differential_stream()
  
  TEG_dem_correct_ele_smooth_srtm <- sum(dif_dem_correct_ele_smooth_srtm$dplus)
  TEG_dem_correct_ele_smooth_dct <- sum(dif_dem_correct_ele_smooth_dct$dplus)
  TEG_dem_correct_ele_smooth_svd <-sum(dif_dem_correct_ele_smooth_svd$dplus)
  
  
  total_distance <- sum(segment$delta_distance)
  
  df_teg_seg_estu <- data.frame(TEG_Gps_data, TEG_Srtm1arc,TEG_C_l_noizy_raster,TEG_C_l_smooth_dct_raster, TEG_C_l_smooth_svd_raster,
                                TEG_Gpstream_correction, TEG_Srtm3voidfilled, TEG_Alos_palsar ,TEG_Nasadem1arc,TEG_TanDEM , TEG_MeritDEM ,TEG_dem_correct_ele_smooth_srtm,
                                TEG_dem_correct_ele_smooth_dct , TEG_dem_correct_ele_smooth_svd, total_distance)

  colnames(df_teg_seg_estu) <- c("TEG GPS", "TEG Srtm1arc","TEG Srtm con ruido","TEG DCT","TEG SVD","TEG SRTM- Elevatr","TEG Srtm3arc","TEG Alos Palsar 2","TEG NASADEM","TEG TanDEM","TEG MeritDEM","TEG SRTM-Elevatr+ LOESS","TEG DCT + LOESS","TEG SVD + LOESS","total_distance")
  

  return(df_teg_seg_estu)
  
  
}