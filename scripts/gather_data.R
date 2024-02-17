# Title: Functions to prepare data for modeling
# Created by: Lei Song
# Created on: 09/09/23
# Note      :
# - gather_data_cs: prepare data for coarse-scale modeling
# - gather_data_fs: prepare data for fine-scale modeling

gather_data_cs <- function(bry, data_path, scale) {
    # Load datasets
    bios <- rast(file.path(data_path, sprintf("bios_%s.tif", scale)))
    ndvis <- rast(file.path(data_path, sprintf("ndvi_%s.tif", scale)))
    roads_density <- rast(file.path(data_path, 
                               sprintf("rivers_density_%s.tif", scale)))
    rivers_density <- rast(file.path(data_path, 
                                     sprintf("roads_density_%s.tif", scale)))
    stm_density <- rast(file.path(data_path, 
                                  sprintf("settlement_density_%s.tif", scale)))
    densities <- c(roads_density, rivers_density, stm_density)
    names(densities) <- c("roads_density", "rivers_density", 
                          "settlement_density")
    rm(roads_density, rivers_density, stm_density)
    terrains <- rast(file.path(data_path, 
                               sprintf("terrain_roughness_%s.tif", scale)))
    ls_l_metrics <- rast(file.path(data_path, 
                                   sprintf("lsp_l_metrics_%s.tif", scale)))
    
    # Need to fill some NAs because not each cell have all types
    ls_c_metrics <- rast(file.path(data_path, 
                                   sprintf("lsp_c_metrics_%s.tif", scale)))
    # Remove clumpy because it is hard to interpret a value for cell without
    # a class, which is common in this case.
    # And clumpy can be described by other variables.
    # All missing values are filled with zero according to their
    # ecological meanings.
    ls_c_metrics <- subset(
        ls_c_metrics, setdiff(names(ls_c_metrics), c("1_clumpy", "2_clumpy")))
    ls_c_metrics <- classify(ls_c_metrics, cbind(NA, 0))
    # Rename
    nms <- gsub("1", "cropland", names(ls_c_metrics))
    nms <- gsub("2", "tree", nms)
    nms <- gsub("3", "savanna", nms)
    nms <- gsub("4", "water", nms)
    names(ls_c_metrics) <- nms
    
    # Put them together, mask and save out
    vars <- c(bios, ndvis, densities, terrains, ls_l_metrics, ls_c_metrics)
    rm(bios, ndvis, densities, terrains, ls_l_metrics, ls_c_metrics)
    vars <- mask(vars, bry, touches = TRUE)
    writeRaster(vars, file.path(data_path, sprintf("variables_%s.tif", scale)))
}

gather_data_fs <- function(bry, data_path) {
    library(purrr)
    library(terra)
    
    # Load datasets
    bios <- rast(file.path(data_path, "bio4.tif"))
    ndvis <- rast(file.path(data_path, "ndvi.tif"))
    dists <- rast(file.path(data_path, "dists_1km.tif"))
    names(dists) <- gsub("_1km", "", names(dists))
    terrains <- rast(file.path(data_path, "vrms.tif"))
    names(terrains) <- paste0("vrm_", c(3, 5, 7))
    ls_c_metrics <- rast(
        file.path(data_path, 
                  sprintf("variables_1km_%sbuf.tif", c(3, 5, 7))))
    names(ls_c_metrics) <- cross(
        list(names(ls_c_metrics)[1:(nlyr(ls_c_metrics) %/% 3)], 
             c(3, 5, 7))) %>% 
        map_chr(paste, sep = "_", collapse = "_")
    
    # # Read land cover and project
    # lc <- rast(file.path(data_path, "landcover_1km.tif"))
    # names(lc) <- "landcover"
    # lc <- project(lc, crs(bios), method = "near")
    # lc <- resample(lc, bios, method = "near")
    # 
    # # Put them together, mask and save out
    # vars <- c(lc, bios, ndvis, dists, terrains, ls_c_metrics)
    # rm(lc, bios, ndvis, dists, terrains, ls_c_metrics)
    
    # Put them together, mask and save out
    vars <- c(bios, ndvis, dists, terrains, ls_c_metrics)
    rm(bios, ndvis, dists, terrains, ls_c_metrics)
    
    vars <- mask(vars, bry, touches = TRUE)
    writeRaster(vars, file.path(data_path, "variables.tif"))
}
