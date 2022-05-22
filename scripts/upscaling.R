library(here)
library(terra)
library(dplyr)
library(APUpscale)

data_path <- here('data')
fname <- file.path(data_path, 'landcover.tif')

#################### scales #####################
# resolution: 30m, 100m, 250m, 500m, 1000m
# scales: fine-scale, medium-scale, medium-scale, 
# medium-scale, coarse-scale
#################################################
# Upscale to 1000m
landscape <- rast(fname)
landscape <- upscale(landscape, cellsize = 1000, verbose = TRUE)
writeRaster(landscape, file.path(data_path, 'landcover_1000m.tif'),
            overwrite = T,
            wopt = list(datatype = 'INT1U',
                        gdal=c("COMPRESS=LZW")))

# Upscale to 500m
landscape <- rast(fname)
landscape <- upscale(landscape, cellsize = 500)
writeRaster(landscape, file.path(data_path, 'landcover_500m.tif'),
            overwrite = T,
            wopt = list(datatype = 'INT1U',
                        gdal=c("COMPRESS=LZW")))

# Upscale to 250m
landscape <- rast(fname)
landscape <- upscale(landscape, cellsize = 250, verbose = TRUE)
writeRaster(landscape, file.path(data_path, 'landcover_250m.tif'),
            overwrite = T,
            wopt = list(datatype = 'INT1U',
                        gdal=c("COMPRESS=LZW")))

# Upscale to 100m
landscape <- rast(fname)
landscape <- upscale(landscape, cellsize = 100, verbose = TRUE)
writeRaster(landscape, file.path(data_path, 'landcover_100m.tif'),
            overwrite = T,
            wopt = list(datatype = 'INT1U',
                        gdal=c("COMPRESS=LZW")))

# Upscale to 30m
landscape <- rast(fname)
landscape <- upscale(landscape, cellsize = 30, verbose = TRUE)
writeRaster(landscape, file.path(data_path, 'landcover_30m.tif'),
            overwrite = T,
            wopt = list(datatype = 'INT1U',
                        gdal=c("COMPRESS=LZW")))
