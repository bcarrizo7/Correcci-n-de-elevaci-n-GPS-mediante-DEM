#https://gist.github.com/eliocamp/5850c629a032fbeae39931d236e82fc7
#' Smooths a 2D field
#' 
#' @param x,y Vector of x and y coordinates
#' @param value Vector of values
#' @param kx,ky Proportion of components to keep in the x and 
#' y direction respectively. Lower values increased the smoothness.

library(elevatr)
library(dplyr)

smooth_dct <- function (x, y, value, kx = 0.1, ky = kx) {
  force(ky)  # Evaluate ky now, because kx will change later.
  m <- mean(value)
  value <- value - m
  
  o <- order(y, x)
  # Transform into a matrix
  matrix <- matrix(value[o], nrow = length(unique(x)),  byrow = FALSE)
  
  # Compute the Discreete Cosine Transform
  # We could also use fft, but it suffers from edge effects
  dct <- gsignal::dct2(matrix)
  
  # Define the components to keep 
  # Remmeber that the FFT is symmetrical
  kx <- c(0, seq_len(floor(nrow(dct)/2 * kx)))
  kx <- c(kx + 1, nrow(dct) - kx[kx != 0] + 1)
  
  ky <- c(0, seq_len(floor(ncol(dct)/2 * ky)))
  ky <- c(ky + 1, ncol(dct) - ky[ky != 0] + 1)
  
  # Replace with zeros and invert
  dct[, -ky] <- 0
  dct[-kx, ] <- 0
  c(gsignal::idct2(dct))[order(o)] + m
}

#' @param  variance_lost Maximum percentage of variance lost after smoothing.
smooth_svd <- function(x, y, value, variance_lost = 0.03) {
  m <- mean(value)
  value <- value - m
  o <- order(y, x)
  # Transform into a matrix
  matrix <- matrix(value[o], nrow = length(unique(x)),  byrow = FALSE)
  total_variance <- norm(abs(matrix), type = "F")
  
  svd <- svd(matrix)
  
  variance_accumulated <- cumsum(svd$d^2/total_variance^2)
  variance_kept <- 1 - variance_lost
  
  keep <- seq_len(which(variance_accumulated - variance_kept > 0)[1])
  
  smooth <- svd$u[, keep] %*% diag(svd$d[keep], nrow = length(keep)) %*% t(svd$v[, keep])
  smooth <- smooth + m
  c(smooth)[order(o)]
}

dct_or_dct <-  function(segment, dct =TRUE ,kx = 0.1, vl=0.05){
  
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
  
  if (dct==TRUE){
    return(C_l_smooth_dct_raster)
  }
  
  else {return(C_l_smooth_svd_raster)}
}


