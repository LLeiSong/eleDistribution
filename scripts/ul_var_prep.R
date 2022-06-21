ul_var_prep <- function(data_path, scale = "1km") {
    checkmate::check_choice(
        scale, choices = c("1km", "2km", "4km", "6km", "8km", "10km"))
    
    # Boundary
    bry <- read_sf(file.path(data_path, "tanzania.geojson"))
    bry <- vect(bry["name"])
    
    # resources
    ndvi <- rast(file.path(data_path, sprintf("ndvi_harmonic_%s.tif", scale)))
    ndvi_mean <- rast(file.path(data_path, sprintf("ndvi_mean_%s.tif", scale)))
    names(ndvi_mean) <- "ndvi_mean"
    names(ndvi) <- paste0("ndvi_", names(ndvi))
    ndvi <- c(ndvi, ndvi_mean); rm(ndvi_mean)
    
    # climate
    wc <- rast(file.path(data_path, "worldclim2_30s.tif"))
    wc <- resample(wc, ndvi[[1]], method = "bilinear")
    nms <- as.integer(str_extract(names(wc), "[0-9]+$"))
    wc <- subset(wc, order(nms)); rm(nms)
    names(wc) <- paste0("bio_", 1:19)
    
    # elevation
    dsm <- rast(file.path(data_path, "dsm_tanzania.tif"))
    dsm <- resample(dsm, ndvi[[1]], method = "bilinear")
    names(dsm) <- "DSM"
    terrains <- terrain(dsm, v = c("slope", "roughness", "TRI"), neighbors = 8)
    
    # Concatenate them
    vars <- c(ndvi, wc, dsm, terrains)
    vars <- mask(crop(vars, bry), bry)
    rm(bry, ndvi, wc, dsm, terrains)
    
    # Write out
    writeRaster(vars, 
                file.path(data_path, sprintf("ul_variables_%s.tif", scale)))
}