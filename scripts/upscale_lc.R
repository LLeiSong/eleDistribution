library(terra)
library(APUpscale)

data_path <- here('data')
fname <- file.path(data_path, 'landcover.tif')

#################### scales #####################
# resolution: 30m, 100m, 250m, 500m, 1000m
# scales: fine-scale, medium-scale, medium-scale, 
# medium-scale, coarse-scale
#################################################
upscale_lc <- function(fname, cellsize, verbose, dst_fname) {
    landscape <- rast(fname)
    # the factor of upscaling to e.g. 30m is too small. Almost all pixels would 
    # have the majority class, so here we just use majority method for upscaling 
    # for speed.
    if (cellsize < 100) {
        output <- rast(ext(landscape), crs = crs(landscape),
                       resolution = cellsize,
                       vals = NA)
        landscape <- resample(landscape, output, method = "mode")
        writeRaster(landscape, dst_fname,
                    overwrite = T,
                    wopt = list(datatype = 'INT1U',
                                gdal=c("COMPRESS=LZW")))
    # Use area preserving method for 
        # coarse resolution of 100m, 250m, 500m, 1000m
    } else {
        landscape <- upscale(landscape, 
                             cellsize = cellsize, 
                             verbose = verbose)
        writeRaster(landscape, file.path(dst_fname),
                    overwrite = T,
                    wopt = list(datatype = 'INT1U',
                                gdal=c("COMPRESS=LZW")))
    }
}
