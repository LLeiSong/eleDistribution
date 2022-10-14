param_tuning_cs <- function(census_block, scale, data_path, 
                            dst_path, sample_weight = TRUE){
    # Load packages
    library(pbmcapply)
    library(caret, quietly = TRUE)
    library(itsdm, quietly = TRUE)
    library(parallel, quietly = TRUE)
    library(stars)
    library(sf)
    
    # Read variables
    vars <- stack(
        file.path(data_path, sprintf("variables_final_%s.tif", scale)))
    
    # Make a non-NA mask
    msk <- lapply(1:nlayers(vars), function(n) !is.na(vars[[n]]))
    msk <- sum(stack(msk)) == nlayers(vars)
    msk[msk == 0] <- NA
    
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
    census_block <- census_block * msk
    names(census_block) <- c("group", "weight", "density")
    rm(template, msk)
    
    vars <- st_as_stars(vars)
    vars <- split(vars, "band")
    
    # Monte Carlo simulation
    ## Important parameters: 
    ## prob_pick_pooled_gain (prefer to lower value here) - 0, 0.1, 0.2
    ## ntry (only if use prob_pick_pooled_gain) - 1, 10, 15, 20
    ## sample_size (prefer to higher value here) - 0.8, 0.9, 1.0
    ## max_depth (prefer to slightly higher value) - 15, 20, 30
    ## ndim (recommended to use small value) - 2, 3, 4
    ## scoring_metric: depth, adj_depth
    ## Since the objective is to rank the environmental conditions, not to
    ## really isolate the outliers, so not use prob_pick_pooled_gain and ntry.
    ## max_depth is rather be higher. 
    ## Not use density based methods because the density of elephants in Tanzania
    ## is not linear related to the environmental conditions.
    
    message("Start tuning.")
    num_sample <- floor(sum(!is.na(values(census_block[[1]]))) * 0.2)
    # Thin the samples to [300, 2000] to reduce overfit
    num_sample <- min(max(300, num_sample), 2000)
    num_sample <- num_sample / 0.9 # add extra for CV
    
    # Convert to points with all info
    census_block <- rasterToPoints(
        census_block, na.rm = TRUE, sp = TRUE) %>% 
        st_as_sf() %>% mutate(id = 1:nrow(.)) %>% 
        mutate(num = ceiling(num_sample * weight)) %>% 
        select(id, group, num, density)
    
    # Define the sets of preferred parameters
    params <- expand.grid(
        sample_size = c(0.8, 0.85, 0.9, 0.95),
        max_depth = seq(30, min(floor(num_sample / 2), 180), 30),
        ndim = c(2, 3, 4),
        scoring_metric = c("depth", "adj_depth"))
    
    tuning_cv <- do.call(rbind, pbmclapply(1:nrow(params), function(n) {
        # Parameters
        params_run <- params[n, ]
        
        # Dynamically make pseudo samples,
        ## make sure set the same seed to reproduce
        pseudo_occ <- do.call(
            rbind, lapply(unique(census_block$group), function(grp_id) {
                set.seed(grp_id)
                census_block %>% filter(group == grp_id) %>% 
                    sample_n(size = min(unique(.$num), nrow(.))) %>% 
                    select(density)}))
        
        # Split occurrence with 5 folds
        set.seed(123)
        flds <- createFolds(1:nrow(pseudo_occ), 
                            k = 10, list = TRUE, 
                            returnTrain = FALSE)
        
        # Cross validation
        it_sdms <- lapply(flds, function(ids) {
            # Split
            occ_sf <- pseudo_occ[setdiff(1:nrow(pseudo_occ), ids), ]
            occ_test_sf <- pseudo_occ[ids, 'geometry']
            
            # sample weight or not
            if (sample_weight) {
                weights <- occ_sf$density
                if (length(weights) != nrow(occ_sf)) {
                    warning("Failed to extract weights!")
                }
            } else weights <- NULL
            
            # Do modeling
            isotree_po(
                occ = occ_sf[, "geometry"],
                occ_test = occ_test_sf,
                variables = vars,
                sample_size = params_run$sample_size,
                ndim = params_run$ndim, 
                max_depth = params_run$max_depth,
                scoring_metric = params_run$scoring_metric,
                seed = 10L,
                # set sample weights as sampling probability
                # so the samples are made population-weighted
                # and the model is affected more by areas with high density
                weights_as_sample_prob = TRUE,
                sample_weights = weights,
                response = FALSE,
                spatial_response = FALSE,
                check_variable = FALSE,
                visualize = FALSE)})
        
        # Collect results
        eval_mean <- do.call(rbind, lapply(it_sdms, function(run) {
            eval_test <- run$eval_test
            tibble("cvi25" = eval_test$po_evaluation$cvi$`cvi with 0.25`,
                   "cvi50" = eval_test$po_evaluation$cvi$`cvi with 0.5`,
                   "cvi75" = eval_test$po_evaluation$cvi$`cvi with 0.75`,
                   "cbi" = eval_test$po_evaluation$boyce$cor,
                   "auc_ratio" = eval_test$po_evaluation$roc_ratio$auc_ratio,
                   "sensitivity" = eval_test$pb_evaluation$sensitivity,
                   "specificity" = eval_test$pb_evaluation$specificity,
                   "TSS" = eval_test$pb_evaluation$TSS$`Optimal TSS`,
                   "auc" = eval_test$pb_evaluation$roc$auc,
                   `Jaccard's similarity index` = eval_test$pb_evaluation$`Jaccard's similarity index`,
                   "f-measure" = eval_test$pb_evaluation$`f-measure`,
                   `Overprediction rate` = eval_test$pb_evaluation$`Overprediction rate`,
                   `Underprediction rate` = eval_test$pb_evaluation$`Underprediction rate`)
        })) %>% summarise(across(everything(), mean))
        
        # Stack the results
        cbind(params_run, eval_mean)
    }, mc.cores = 5))
    
    save(tuning_cv, 
         file = file.path(dst_path, sprintf("tuning_cv_%s.rda", scale)))
}

modeling_cs <- function(params, census_block, occ, 
                        scale, data_path, dst_path,
                        sample_weight = TRUE) {
    # Load packages
    library(pbmcapply, quietly = TRUE)
    library(caret, quietly = TRUE)
    library(itsdm, quietly = TRUE)
    library(parallel, quietly = TRUE)
    library(stars)
    library(sf)
    
    # Read variables
    vars <- stack(
        file.path(data_path, sprintf("variables_final_%s.tif", scale)))
    
    # Make a non-NA mask
    msk <- lapply(1:nlayers(vars), function(n) !is.na(vars[[n]]))
    msk <- sum(stack(msk)) == nlayers(vars)
    msk[msk == 0] <- NA
    
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
    census_block <- census_block * msk
    names(census_block) <- c("group", "weight", "density")
    occ <- rasterize(occ, template, fun = "count")
    occ <- occ * msk
    rm(template, msk)
    
    vars <- st_as_stars(vars)
    vars <- split(vars, "band")
    
    # Start to run the simulations
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
    
    # Print
    message(sprintf("Params: sample_size: %s, max_depth: %s, ndim: %s",
                    params$sample_size, params$max_depth, params$ndim))
    
    # Run 10 fold runs for each set of parameters
    runs <- pbmclapply(1:num_cv, function(n_cv) {
        # Dynamically make pseudo samples
        ## Unlike tuning, make sure use dynamic seed each run
        pseudo_occ <- do.call(
            rbind, lapply(unique(census_block$group), function(grp_id) {
                set.seed(123 + n_cv * grp_id)
                census_block %>% filter(group == grp_id) %>% 
                    sample_n(size = min(unique(.$num), nrow(.))) %>% 
                    select(density)}))
        
        # Get real occurrence for testing
        occ <- rasterToPoints(occ, spatial = TRUE) %>% 
            st_as_sf()
        
        # Do modeling
        if (sample_weight) {
            weights <- pseudo_occ$density
            if (length(weights) != nrow(pseudo_occ)) {
                warning("Failed to extract weights!")
            }
        } else weights <- NULL
        
        isotree_po(occ = pseudo_occ[, "geometry"],
                   occ_test = occ[, "geometry"],
                   variables = vars,
                   ndim = params$ndim,
                   sample_size = params$sample_size,
                   max_depth = params$max_depth,
                   scoring_metric = params$scoring_metric,
                   seed = 10L,
                   # set sample weights as distribution density
                   weights_as_sample_prob = TRUE,
                   sample_weights = weights,
                   # switch spatial response off
                   spatial_response = FALSE)
    }, mc.cores = 5)
    names(runs) <- paste(
        sprintf("%s_%s_%s_%s", 
                params$sample_size, params$max_depth,
                params$ndim, params$scoring_metric), 
        1:num_cv, sep = "_")
    
    save(runs, file = file.path(dst_path, sprintf("runs_%s.rda", scale)))
    
    # Return
    runs
}

# Define a function to gather modeling metrics
gather_metrics <- function(eval_list) {
    po <- eval_list$po_evaluation
    pb <- eval_list$pb_evaluation
    tibble(
        cvi_25 = po$cvi$`cvi with 0.25`,
        cvi_50 = po$cvi$`cvi with 0.5`,
        cvi_75 = po$cvi$`cvi with 0.75`,
        boyce_index = po$boyce$cor,
        roc_ratio = po$roc_ratio$auc_ratio, # a typo here
        sensitivity = pb$sensitivity,
        specificity = pb$specificity,
        tss = pb$TSS$`Optimal TSS`,
        cutoff = pb$TSS$`Recommended threshold`,
        roc = pb$roc$auc, # a typo here
        `Jaccard's similarity index` = pb$`Jaccard's similarity index`,
        `f-measure` = pb$`f-measure`,
        `Overprediction rate` = pb$`Overprediction rate`,
        `Underprediction rate` = pb$`Underprediction rate`)
}

univariate_model <- function(occ, variable) {
    # Load packages
    library(caret, quietly = TRUE)
    library(itsdm, quietly = TRUE)
    library(stars)
    library(sf)
    library(raster)
    
    # Get a template for pseudo records
    template <- as(variable, "Raster")
    values(template) <- NA
    occ <- rasterize(occ, template, fun = "count")
    occ <- rasterToPoints(occ, spatial = TRUE) %>% 
        st_as_sf() %>% select()
    
    # Split occurrence with 10 folds
    set.seed(123)
    flds <- createFolds(1:nrow(occ), 
                        k = 10, list = TRUE, 
                        returnTrain = FALSE)
    
    # Cross validation
    it_sdms <- lapply(flds, function(ids) {
        # Split
        occ_sf <- occ[setdiff(1:nrow(occ), ids), ]
        occ_test_sf <- occ[ids, 'geometry']
        
        # Do modeling
        isotree_po(
            occ = occ_sf[, "geometry"],
            occ_test = occ_test_sf,
            variables = variable,
            sample_size = 0.8,
            ndim = 1, 
            seed = 10L,
            response = FALSE,
            spatial_response = FALSE,
            check_variable = FALSE,
            visualize = FALSE)})
    
    # Collect results
    do.call(rbind, lapply(it_sdms, function(run) {
        gather_metrics(run$eval_test)
    })) %>% mutate(fold = paste("fold", 1:10, sep = "_"))
}

param_tuning_fs <- function(occ, data_path, dst_path){
    # Load packages
    library(pbmcapply)
    library(caret, quietly = TRUE)
    library(itsdm, quietly = TRUE)
    library(parallel, quietly = TRUE)
    library(stars)
    library(raster)
    library(sf)
    
    # Read variables
    vars <- stack(file.path(data_path, "variables_final.tif"))
    
    # Get a template for pseudo records
    template <- vars[[1]]
    values(template) <- NA
    occ <- rasterize(occ, template, fun = "count")
    
    vars <- st_as_stars(vars)
    vars <- split(vars, "band") %>% 
        mutate(landcover = factor(landcover))
    
    # Convert to points with all info
    occ <- rasterToPoints(occ, spatial = TRUE) %>% 
        st_as_sf() %>% select()
    
    # Define the sets of preferred parameters
    params <- expand.grid(
        sample_size = c(0.85, 0.9, 0.95, 0.1),
        max_depth = seq(30, min(floor(nrow(occ) / 2), 180), 30),
        ndim = c(2, 3, 4),
        scoring_metric = c("depth", "adj_depth"))
    
    tuning_cv <- do.call(rbind, pbmclapply(1:nrow(params), function(n) {
        # Parameters
        params_run <- params[n, ]
        
        # Split occurrence with 10 folds
        set.seed(123)
        flds <- createFolds(1:nrow(occ), 
                            k = 10, list = TRUE, 
                            returnTrain = FALSE)
        
        # Cross validation
        it_sdms <- lapply(flds, function(ids) {
            # Split
            occ_sf <- occ[setdiff(1:nrow(occ), ids), ]
            occ_test_sf <- occ[ids, 'geometry']
            
            # Do modeling
            isotree_po(
                occ = occ_sf[, "geometry"],
                occ_test = occ_test_sf,
                variables = vars,
                categ_vars = "landcover",
                sample_size = params_run$sample_size,
                ndim = params_run$ndim, 
                max_depth = params_run$max_depth,
                scoring_metric = params_run$scoring_metric,
                seed = 10L,
                response = FALSE,
                spatial_response = FALSE,
                check_variable = FALSE,
                visualize = FALSE)})
        
        # Collect results
        eval_mean <- do.call(rbind, lapply(it_sdms, function(run) {
            gather_metrics(run$eval_test)
        })) %>% summarise(across(everything(), mean))
        
        # Stack the results
        cbind(params_run, eval_mean)
    }, mc.cores = detectCores() - 1))
    
    save(tuning_cv, file = file.path(dst_path, "tuning_cv_fs.rda"))
}

modeling_fs <- function(params, occ, data_path, dst_path, seed = 10) {
    # Load packages
    library(caret, quietly = TRUE)
    library(itsdm, quietly = TRUE)
    library(parallel, quietly = TRUE)
    library(stars, quietly = TRUE)
    library(raster, quietly = TRUE)
    library(sf, quietly = TRUE)
    library(spThin, quietly = TRUE)
    
    # Read variables
    vars <- stack(file.path(data_path, "variables_final.tif"))
    
    # Get a template for pseudo records
    template <- vars[[1]]
    values(template) <- NA
    occ <- rasterize(occ, template, fun = "count")
    
    vars <- st_as_stars(vars)
    vars <- split(vars, "band") %>% 
        mutate(landcover = factor(landcover))
    
    # Run 10 fold runs with the best parameters
    # Convert to points with all info
    occ <- rasterToPoints(occ, spatial = TRUE) %>% 
        st_as_sf() %>% select()
    
    # Split occurrence with 10 folds
    set.seed(seed)
    flds <- createFolds(1:nrow(occ), 
                        k = 10, list = TRUE, 
                        returnTrain = FALSE)
    
    # Print
    message(sprintf("Params: sample_size: %s, max_depth: %s, ndim: %s",
                    params$sample_size, params$max_depth, params$ndim))
    
    # Cross validation
    runs <- lapply(1:length(flds), function(n) {
        message(sprintf("Run %s", n))
        
        ids <- flds[[n]]
        # Split
        occ_sf <- occ[setdiff(1:nrow(occ), ids), ]
        occ_test_sf <- occ[ids, 'geometry']
        
        # Do modeling
        isotree_po(
            occ = occ_sf[, "geometry"],
            occ_test = occ_test_sf,
            variables = vars,
            categ_vars = "landcover",
            ndim = params$ndim,
            sample_size = params$sample_size,
            max_depth = params$max_depth,
            scoring_metric = params$scoring_metric,
            seed = 10L,
            spatial_response = FALSE)
    })
    
    # Gather things together
    names(runs) <- paste(
        sprintf("%s_%s_%s_%s", 
                params$sample_size, params$max_depth,
                params$ndim, params$scoring_metric), 
        1:length(flds), sep = "_")
    
    save(runs, file = file.path(dst_path, "runs_fs.rda"))
    
    # Return
    runs
}


modeling_fs_thin <- function(params, occ, cluster, 
                             data_path, dst_path, 
                             seed = 10) {
    # Load packages
    library(caret, quietly = TRUE)
    library(itsdm, quietly = TRUE)
    library(parallel, quietly = TRUE)
    library(stars, quietly = TRUE)
    library(raster, quietly = TRUE)
    library(sf, quietly = TRUE)
    library(spThin, quietly = TRUE)
    
    # Read variables
    vars <- stack(file.path(data_path, "variables_final.tif"))
    
    # Get a template for pseudo records
    template <- vars[[1]]
    values(template) <- NA
    occ <- rasterize(occ, template, fun = "count")
    
    vars <- st_as_stars(vars)
    vars <- split(vars, "band") %>% 
        mutate(landcover = factor(landcover))
    vars <- vars %>% 
        select(-c(landcover, savanna_contig_mn_3, dist_to_big_roads))
    
    # Run 10 fold runs with the best parameters
    # Convert to points with all info
    occ <- rasterToPoints(occ, spatial = TRUE) %>% 
        st_as_sf() %>% select()
    
    # Run the model 10 folds with thinned occurrence
    num_to_run <- 10
    # Print
    message(sprintf("Params: sample_size: %s, max_depth: %s, ndim: %s",
                    params$sample_size, params$max_depth, params$ndim))
    
    # Cross validation
    ## Split data
    occ_keep_id <- st_intersects(cluster, occ) %>% unlist() %>% unique()
    occ_keep <- occ %>% slice(setdiff(1:nrow(occ), occ_keep_id)) %>% 
        mutate(id = 1:nrow(.))
    occ_to_thin <- occ %>% slice(occ_keep_id) %>% st_coordinates() %>% 
        data.frame() %>% mutate(SPEC = "SPEC")
    
    ## Run
    runs <- lapply(1:num_to_run, function(n) {
        message(sprintf("Run %s", n))
        
        # Get the train and test dataset
        set.seed(seed + n)
        occ_keep_train <- occ_keep %>% sample_frac(size = 0.7)
        occ_keep_test <- occ_keep %>% filter(!id %in% occ_keep_train$id)
        
        # Thin the samples in northeast
        set.seed(seed * 10 + n)
        occ_thin <- thin(
            loc.data = occ_to_thin, 
            lat.col = "Y", long.col = "X", 
            spec.col = "SPEC", 
            thin.par = 20, reps = 100,
            locs.thinned.list.return = TRUE, 
            write.files = FALSE, 
            write.log.file = FALSE) %>% 
            `[[`(100)
        
        occ_sf <- occ_thin %>% 
            st_as_sf(coords = c("Longitude", "Latitude"),
                     crs = 4326) %>% 
            rbind(occ_keep_train %>% select(-id))
        occ_test_sf <- occ_to_thin %>% select(-SPEC) %>% 
            slice(setdiff(1:nrow(occ_to_thin), 
                          as.integer(row.names(occ_thin)))) %>% 
            st_as_sf(coords = c("X", "Y"),
                     crs = 4326) %>% 
            rbind(occ_keep_test %>% select(-id))
        
        rm(occ_keep_train, occ_keep_test, occ_thin)
        
        # Do modeling
        isotree_po(
            occ = occ_sf[, "geometry"],
            occ_test = occ_test_sf,
            variables = vars,
            # categ_vars = "landcover",
            ndim = params$ndim,
            sample_size = params$sample_size,
            max_depth = params$max_depth,
            scoring_metric = params$scoring_metric,
            seed = 10L,
            spatial_response = FALSE)
    })
    
    save(runs, file = file.path(dst_path, "runs_fs_thin_noname.rda"))
    
    # Gather things together
    names(runs) <- paste(
        sprintf("%s_%s_%s_%s", 
                params$sample_size, params$max_depth,
                params$ndim, params$scoring_metric), 
        1:num_to_run, sep = "_")
    
    save(runs, file = file.path(dst_path, "runs_fs_thin.rda"))
    
    # Return
    runs
}
