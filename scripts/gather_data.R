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

# gather data for Maxent
gather_data_maxent <- function(census_block, occ, scale, data_path, dst_path) {
    library(purrr)
    library(raster)
    
    # Write asc for Maxent
    vars <- stack(
        file.path(data_path, sprintf("variables_final_%s.tif", scale)))
    
    asc_path <- file.path(dst_path, scale)
    dir.create(asc_path)
    walk(1:nlayers(vars), function(n) {
        layer <- vars[[n]]
        writeRaster(
            layer, file.path(asc_path, sprintf('%s.asc', names(layer))))
    })
    
    # generate pseudo samples
    # Make a non-NA mask
    msk <- lapply(1:nlayers(vars), function(n) !is.na(vars[[n]]))
    msk <- sum(stack(msk)) == nlayers(vars)
    msk[msk == 0] <- NA
    vars <- vars * msk
    
    # Add group id for census_block
    census_block <- census_block %>% 
        mutate(weight = estimate / sum(.$estimate)) %>% 
        mutate(group = 1:nrow(.))
    
    # Get a template for pseudo records
    template <- vars[[1]]
    values(template) <- NA
    census_block <- stack(
        rasterize(census_block, template, field = "group"),
        rasterize(census_block, template, field = "weight"),
        rasterize(census_block, template, field = "density"))
    names(census_block) <- c("group", "weight", "density")
    occ <- rasterize(occ, template, fun = "count")
    rm(template, msk)
    
    num_cv <- 10
    num_sample <- floor(sum(!is.na(values(census_block[[1]]))) * 0.2)
    # Thin the samples to [300, 2000] to reduce overfit
    num_sample <- min(max(300, num_sample), 2000)
    
    # Convert to points with all info
    census_block <- rasterToPoints(
        census_block, na.rm = TRUE, sp = TRUE) %>% 
        st_as_sf() %>% mutate(id = 1:nrow(.)) %>% 
        mutate(num = ceiling(num_sample * weight)) %>% 
        select(id, group, num, density)
    
    pseudo_occ <- do.call(rbind, lapply(1:num_cv, function(n_cv) {
        # Dynamically make pseudo samples
        set.seed(123 + n_cv)
        do.call(rbind, lapply(unique(census_block$group), 
                              function(grp_id) {
                set.seed(grp_id)
                census_block %>% filter(group == grp_id) %>% 
                    sample_n(size = min(unique(.$num), nrow(.))) %>% 
                    select(density)})) %>% 
            st_coordinates() %>% data.frame() %>% 
            mutate(Long = X, Lat = Y, Species = paste0("Elephant_", n_cv)) %>% 
            select(Species, Long, Lat)
    }))
    
    # Save out
    smp_path <- file.path(dst_path, "pseudo_samples")
    if (!dir.exists(smp_path)) dir.create(smp_path)
    write.csv(pseudo_occ, row.names = FALSE,
              file.path(smp_path, sprintf("maxent_pseudos_%s.csv", scale)))
    
    ## Use real occurrence to make independent test dataset
    occ <- rasterToPoints(occ) %>% data.frame()
    occ <- data.frame(
        species = paste0("Elephant_", sort(rep(1:10, nrow(occ)))),
        Long = rep(occ$x, 10), Lat = rep(occ$y, 10))
    write.csv(occ, row.names = FALSE,
              file.path(smp_path, sprintf("maxent_test_samples_%s.csv", scale)))
}
