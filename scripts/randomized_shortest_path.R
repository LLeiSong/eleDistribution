generate_nodes <- function(nodes, value) {
    source_nodes <- nodes
    source_nodes <- source_nodes == value
    source_nodes[is.na(source_nodes)] <- FALSE
    values(source_nodes) <- as.logical(values(source_nodes))
    return(source_nodes)
}

rand_shortest_path <- function(
        habitat_clusters, pas, 
        cluster_id, type_index, 
        cnt_path,
        cnt_exp_path,
        n_iteration = 20,
        nthreads = 6) { 
    # Load packages
    library(raster)
    library(sf)
    library(terra)
    library(dplyr)
    library(gdistance)
    library(pbmcapply)
    
    # Types
    types <- c("mnt_block", "human_block", "human_block_complete")
    type <- types[type_index]
    thetas <- c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5)
    
    # Subset the whole region
    cluster_msk <- habitat_clusters %>% 
        filter(cluster == cluster_id) %>% 
        vect() %>% ext()
    
    # Subset suitability
    suit_map <- file.path(
        cnt_path, sprintf("suitability_%s.tif", type)) %>% 
        rast() %>% crop(cluster_msk) %>% raster()
    
    # Get the transition matrix
    tr <- suit_map
    tr[tr == 0] <- 0.0001
    tr <- transition(
        tr, transitionFunction = mean, 
        directions = 16, symm = FALSE)
    tr <- geoCorrection(tr, type = "c", scl = TRUE)
    tr <- geoCorrection(tr, type = "r", scl = TRUE)
    
    message("Preparation is finished. Start the iteration.")
    
    # Iterations
    cum_raster <- suit_map
    values(cum_raster) <- 0
    thetas_select <- c()
    for (n in 1:n_iteration) {
        message(sprintf("Iteration -- %s", n))
        
        # Get the nodes
        set.seed(cluster_id * 100 + 10 * type_index + n)
        nodes <- pas %>% filter(cluster == cluster_id) %>% 
            st_sample(size = rep(1, nrow(.))) %>%
            st_as_sf() %>% mutate(id = 1:nrow(.)) %>% 
            vect() %>% buffer(width = 5000) %>% st_as_sf() %>% 
            rasterize(suit_map, field = "id") %>%
            rasterToPoints() %>% as.data.frame() %>% 
            rename(id = layer)
        
        # Do experiments in each node pairs
        ids <- unique(nodes$id)
        for (s in ids) {
            source_node <- nodes %>% filter(id == s) %>% 
                dplyr::select(-id) %>% as.matrix()
             
            cum_raster <- pbmclapply(setdiff(ids, s), function(d) {
                # Pick a theta
                theta_select <- sample(thetas, size = 1)
                
                dst_node <- nodes %>% filter(id == d) %>% 
                    dplyr::select(-id) %>% as.matrix()
                
                flow <- passage(
                    tr, source_node, dst_node, 
                    theta = theta_select, totalNet = "net")
                flow[c(cellFromXY(suit_map, source_node), 
                       cellFromXY(suit_map, dst_node))] <- 0
                
                list("cum_raster" = cum_raster, "theta_select" = theta_select)
            }, mc.cores = min(length(ids), nthreads))
            
            cum_rasters <- do.call(
                stack, lapply(cum_raster, function(each) each$cum_raster))
            thetas_select <- do.call(
                c, lapply(cum_rasters, function(each) each$thetas_select))
            
            # Accumulate them
            cum_raster <-  cum_raster + flow
            thetas_select <- c(thetas_select, theta_select)
        }
    }
    
    # Write out the final cumulative raster
    writeRaster(cum_raster, file.path(
        cnt_exp_path, sprintf("cum_rsp_%s_%s.tif", type, cluster_id)))
    
    # Write out the selected thetas
    save(thetas_select, file = file.path(
        cnt_exp_path, sprintf("thetas_rsp_%s_%s.rda", type, cluster_id)))
}


# Start the experiments
library(here)
library(terra)
library(dplyr)
library(sf)
library(optparse)

data_path <- here("data")
result_path <- here("results")
cnt_path <- file.path(data_path, "connectivity")
cnt_exp_path <- file.path(result_path, "connectivity")

# Read the fixed files
habitat_clusters <- read_sf(
    file.path(data_path, "observations",
              "habitat_clusters_aggre.geojson"))
pas <- read_sf(file.path(cnt_path, "pas.geojson"))

# Parse inline arguments
option_list <- list(
    make_option(c("-c", "--cluster_id"), 
                action = "store", default = 1, type = 'integer', 
                help = "the index of cluster to process from [1, 2, 3]"),
    make_option(c("-t", "--type_index"), 
                action = "store", default = 1, type = 'integer',
                help = paste0(
                    "the condition to process from",
                    " [mnt_block, human_block, human_block_complete]")),
    make_option(c("-n", "--n_iteration"), 
                action = "store", default = 10, type = 'integer',
                help = "No. of time to iteration [default %default]"),
    make_option(c("-t", "--nthreads"), 
                action = "store", default = 10, type = 'integer',
                help = "No. of threads to use [default %default]"))
opt <- parse_args(OptionParser(option_list = option_list))

# Run it
rand_shortest_path(habitat_clusters, pas, opt$cluster_id, 
                   opt$type_index, cnt_path, cnt_exp_path, 
                   opt$n_iteration, opt$nthreads)
