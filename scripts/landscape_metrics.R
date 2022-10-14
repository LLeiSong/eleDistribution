######## Landscape level #############
# Aggregation metric:
## Aggregation index (AI)
## Patch density (PD)
# Complexity metric:
## Marginal entropy (ENT)
## Conditional entropy
## Joint entropy
## Mutual information
# Area and edge metric:
## Edge density (ED)
# Diversity metric:
## Patch richness density (PRD)
## Simpson’s diversity index (SIDI)
## Shannon’s diversity index (SHDI)
landscape_l_metrics <- function(landscape, template, dst_path) {
    # loading packages
    library(landscapemetrics)
    library(parallel)
    library(pbmcapply)
    library(raster)
    library(dplyr)
    library(terra)
    library(tidyr)
    
    # landscape level
    # Calculate the metrics in parallel
    # Define metrics to calculate
    metrics <- c("lsm_l_ai", "lsm_l_pd", "lsm_l_ent", "lsm_l_condent",
                 "lsm_l_joinent", "lsm_l_mutinf", "lsm_l_ed", 
                 "lsm_l_prd", "lsm_l_sidi", "lsm_l_shdi")
    dst_rst <- template
    # template <- project(template, crs(landscape), method = 'near')
    
    # Convert to raster for parallel
    landscape <- raster(landscape)
    template <- raster(template)
    
    # Calculate the metrics for each cell
    message("Calculate landscape metrics for each BIO cell.")
    vals <- do.call(rbind, pbmclapply(1:ncell(template), function(id) {
        blk <- classify(rast(template), cbind(id, 1), others = NA)
        if (isTRUE(global(blk, fun = "isNA") == ncell(blk))) {
            return(NA)
        } else {
            blk <- trim(blk)
            blk <- project(blk, crs(rast(landscape)))
            
            # Do calculation for each window size
            if (relate(as.polygons(ext(landscape)), 
                       as.polygons(ext(blk)), "intersects")) {
                # Crop the landscape
                lc_blk <- crop(rast(landscape), blk)
                
                # Calculate the metrics
                if (isTRUE(global(lc_blk, fun = "isNA") != ncell(lc_blk))) {
                    calculate_lsm(lc_blk, what = metrics) %>% 
                        dplyr::select(metric, value) %>% 
                        pivot_wider(names_from = metric)
                } else {
                    return(NA)
                }
            } else return(NA)
        }
    }, mc.cores = detectCores() - 1))
    
    # Burn values into raster stack
    message("Burn the values into raster stack.")
    save(vals, file = gsub(".tif", ".rda", dst_path))
    spatial_metrics <- rep(dst_rst, ncol(vals))
    values(spatial_metrics) <- as.matrix(vals)
    names(spatial_metrics) <- names(vals)
    writeRaster(spatial_metrics, dst_path)
}

########## Class level ###############
# All class
## class ratio
## dominant class
# Dense-tree and cropland
## Clumpiness index (CLUMPY)
## Edge density (ED)
# Savanna (shrub + grass)
## Largest patch index (LPI)
## Patch density (PD)
## Contiguity index distribution (mean) (CONTIG_MN)
landscape_c_metrics <- function(landscape, template, dst_path) {
    # loading packages
    library(landscapemetrics)
    library(parallel)
    library(pbmcapply)
    library(raster)
    library(dplyr)
    library(terra)
    library(tidyr)
    
    # landscape level
    # Calculate the metrics in parallel
    # Define metrics to calculate
    metrics <- c("lsm_c_clumpy", "lsm_c_contig_mn", 
                 "lsm_c_lpi", "lsm_c_pd", "lsm_c_ed")
    dst_rst <- template
    # template <- project(template, crs(landscape), method = 'near')
    
    # Convert to raster for parallel
    landscape <- raster(landscape)
    template <- raster(template)
    
    # Calculate the metrics for each cell
    message("Calculate landscape metrics for each BIO cell.")
    vals <- do.call(rbind, pbmclapply(1:ncell(template), function(id) {
        blk <- classify(rast(template), cbind(id, 1), others = NA)
        if (isTRUE(global(blk, fun = "isNA") == ncell(blk))) {
            return(NA)
        } else {
            blk <- trim(blk)
            blk <- project(blk, crs(rast(landscape)))
            
            # Read ocean mask
            ocean <- vect("data/vectors/ocean_mask.geojson")
            ocean <- project(ocean, crs(blk))
            
            # Do calculation for each window size
            if (relate(as.polygons(ext(landscape)), 
                       as.polygons(ext(blk)), "intersects")) {
                # Crop the landscape
                lc_blk <- crop(rast(landscape), blk)
                
                if (relate(as.polygons(ext(blk)), 
                           ocean, "intersects")) {
                    lc_blk <- terra::mask(lc_blk, ocean, inverse = TRUE)
                }
                
                # Calculate the metrics
                if (isTRUE(global(lc_blk, fun = "isNA") != ncell(lc_blk))) {
                    values <- calculate_lsm(lc_blk, what = metrics)
                    
                    # Add NA for missing class
                    ## some hardcoded values here
                    classes <- unique(values$class)
                    missing_class <- setdiff(1:3, classes)
                    if (length(missing_class) > 0) {
                        nms <- unique(values$metric)
                        values <- rbind(
                            values,
                            data.frame(
                                layer = 1,
                                level = "class",
                                class = sort(rep(missing_class, length(nms))),
                                id = NA,
                                metric = rep(nms, length(missing_class)),
                                value = NA)) %>% 
                            arrange(metric, class)}
                    
                    # Group metrics
                    crop_tree <- values %>% 
                        filter(class %in% c(1, 2) & 
                                   metric %in% c("clumpy", "ed")) %>%
                        mutate(metric = paste(class, metric, sep = "_")) %>%
                        dplyr::select(metric, value) %>% 
                        pivot_wider(names_from = metric)
                    savanna <- values %>% 
                        filter(class %in% c(3) & 
                                   metric %in% c("lpi", "pd", "contig_mn")) %>%
                        mutate(metric = paste(class, metric, sep = "_")) %>%
                        dplyr::select(metric, value) %>% 
                        pivot_wider(names_from = metric)
                    
                    # ratio
                    ratios <- matrix(freq(lc_blk)[, c("value", 'count')], 
                                     ncol = 2, 
                                     dimnames = list(NULL, c("value", "count")))
                    ratios <- right_join(data.frame(ratios), 
                                         data.frame(value = 1:4, 
                                                    metric = paste(1:4, "ratio", 
                                                                   sep = "_")),
                                         by = "value") %>% arrange(value) %>% 
                        mutate(value = ifelse(is.na(count), 0, 
                                              count / ncell(lc_blk))) %>% 
                        dplyr::select(value, metric) %>% 
                        pivot_wider(names_from = metric)
                    
                    cbind(crop_tree, savanna, ratios)
                } else {
                    return(NA)
                }
            } else return(NA)
        }
    }, mc.cores = detectCores() - 1))
    
    # Burn values into raster stack
    message("Burn the values into raster stack.")
    save(vals, file = gsub(".tif", ".rda", dst_path))
    spatial_metrics <- rep(dst_rst, ncol(vals))
    values(spatial_metrics) <- as.matrix(vals)
    names(spatial_metrics) <- names(vals)
    writeRaster(spatial_metrics, dst_path)
}

landscape_c_metrics_mv <- function(data_path, buf) {
    # loading packages
    library(landscapemetrics)
    library(parallel)
    library(pbmcapply)
    library(raster)
    library(dplyr)
    library(terra)
    library(tidyr)
    
    temp <- rast(file.path(data_path, "worldclim/bios_0.5.tif"), lyrs = 1)
    values(temp) <- 1:ncell(temp)
    lc_class <- raster(
        file.path(data_path, "landcover/crop_tree_savanna_water.tif"))
    
    template <- raster(temp)
    n_col <- ncol(template); crs_to <- crs(template)
    xre <- res(template)[1]; yre <- res(template)[2]
    xre_buf <- xre * buf; yre_buf <- yre * buf
    xmin <- extent(template)[1]; ymax <- extent(template)[4]
    
    dst_path <- file.path(data_path, "vars_patch", 
                          sprintf("variables_1km_%sbuf.tif", buf))
    csv_path <- gsub(".tif", ".csv", dst_path)
    
    ################ Fresh new process ###################
    file.create(csv_path)
    write.table(setNames(data.frame(matrix(ncol = 7, nrow = 0)), 
                         c("id", 
                           "cropland_ed", "savanna_area_mn", "savanna_contig_mn",
                           "savanna_pd", "cropland_ratio", "savanna_ratio")), 
                file = csv_path, append = TRUE, 
                sep = ',', row.names = FALSE)
    
    ########### Recover the interrupted process ############
    done <- read.csv(csv_path, stringsAsFactors = FALSE)
    ids <- setdiff(1:ncell(template), done$id)
    
    sf::sf_use_s2(FALSE)
    vals <- do.call(rbind, pbmclapply(ids, function(id) {
        x_id <- id %% n_col; y_id <- id %/% n_col
        
        # Pixel
        x_min <- xmin + xre * (x_id - 1); x_max <- x_min + xre
        y_max <- ymax - yre * y_id; y_min <- y_max - yre
        
        # centroid
        x_ctd <- x_min + (x_max - x_min) / 2
        y_ctd <- y_min + (y_max - y_min) / 2
        
        # Buffered area
        x_min <- x_ctd - xre_buf / 2; x_max <- x_ctd + xre_buf / 2
        y_min <- y_ctd - yre_buf / 2; y_max <- y_ctd + yre_buf / 2
        
        # Get the buffered area
        blk <- rast(raster(matrix(1),
                           xmn = x_min, xmx = x_max,
                           ymn = y_min, ymx = y_max, 
                           crs = crs_to))
        
        # Start to calculate the variables
        # Landscape
        metrics_c <- c("lsm_c_contig_mn", "lsm_c_area_mn", "lsm_c_pd", "lsm_c_ed")
        
        blk_lc <- project(blk, crs(rast(lc_class)))
        # Read ocean mask
        ocean <- vect("data/vectors/ocean_mask.geojson")
        ocean <- project(ocean, crs(blk_lc))
        
        # Do calculation for each window size
        if (relate(as.polygons(ext(lc_class)), 
                   as.polygons(ext(blk_lc)), "intersects")) {
            # Crop the landscape
            lc_blk <- crop(rast(lc_class), blk_lc)
            
            if (relate(as.polygons(ext(blk_lc)), 
                       ocean, "intersects")) {
                lc_blk <- terra::mask(lc_blk, ocean, inverse = TRUE)
            }
            
            # Calculate the metrics
            if (isTRUE(global(lc_blk, fun = "isNA") != ncell(lc_blk))) {
                values <- calculate_lsm(lc_blk, what = metrics_c)
                
                # Add NA for missing class
                ## some hardcoded values here
                classes <- unique(values$class)
                missing_class <- setdiff(1:3, classes)
                if (length(missing_class) > 0) {
                    nms <- unique(values$metric)
                    values <- rbind(
                        values,
                        data.frame(
                            layer = 1,
                            level = "class",
                            class = sort(rep(missing_class, length(nms))),
                            id = NA,
                            metric = rep(nms, length(missing_class)),
                            value = 0)) %>% 
                        arrange(metric, class)}
                
                # Group metrics
                crop <- values %>% 
                    filter(class == 1 & metric  == "ed") %>%
                    mutate(metric = paste(class, metric, sep = "_")) %>%
                    dplyr::select(metric, value) %>% 
                    pivot_wider(names_from = metric)
                savanna <- values %>% 
                    filter(class %in% c(3) & metric != "ed") %>%
                    mutate(metric = paste(class, metric, sep = "_")) %>%
                    dplyr::select(metric, value) %>% 
                    pivot_wider(names_from = metric)
                
                # ratio
                ratios <- data.frame(freq(lc_blk)[, c("value", 'count')])
                ratios <- right_join(ratios, 
                                     data.frame(value = c(1, 3), 
                                                metric = paste(c(1, 3), "ratio", 
                                                               sep = "_")),
                                     by = "value") %>% arrange(value) %>% 
                    mutate(value = ifelse(is.na(count), 0, 
                                          count / ncell(lc_blk))) %>% 
                    dplyr::select(value, metric) %>% 
                    pivot_wider(names_from = metric)
                lc_metrics_c <- as_tibble(cbind(crop, savanna, ratios))
            } else {
                lc_metrics_c <- tibble(
                    "1_ed" = 0, "3_area_mn" = 0, "3_contig_mn" = 0,
                    "3_pd" = 0, "1_ratio" = 0, "3_ratio" = 0)
            }
        } else lc_metrics_c <- tibble(
            "1_ed" = 0, "3_area_mn" = 0, "3_contig_mn" = 0,
            "3_pd" = 0, "1_ratio" = 0, "3_ratio" = 0)
        
        nms <- gsub("1", "cropland", names(lc_metrics_c))
        nms <- gsub("2", "tree", nms)
        nms <- gsub("3", "savanna", nms)
        nms <- gsub("4", "water", nms)
        names(lc_metrics_c) <- nms
        
        # Collect all of them
        write.table(cbind(data.frame(id = id), lc_metrics_c), 
                    file = csv_path, append = TRUE, 
                    sep = ',', row.names = FALSE, 
                    col.names = FALSE)
        lc_metrics_c
    }, mc.cores = detectCores(), ignore.interactive = TRUE))
    
    # Burn values into raster stack
    message("Burn the values into raster stack.")
    save(vals, file = gsub(".tif", ".rda", dst_path))
    vars <- rep(temp, ncol(vals))
    values(vars) <- as.matrix(vals)
    names(vars) <- names(vals)
    writeRaster(vars, dst_path)
}