library(here)
library(terra)
library(ggplot2)
library(dplyr)
library(stars)
library(ini)
library(purrr)
library(JuliaCall)
library(parallel)
library(optparse)

# Set directories
data_path <- here("data")
result_path <- here("results")
cnt_path <- file.path(data_path, "connectivity")
cnt_exp_path <- file.path(result_path, "connectivity")

# Read PAs
pas <- read_sf(file.path(cnt_path, "pas.geojson"))

# Parse inline arguments
option_list <- list(
    make_option(c("-c", "--cluster_id"), 
                action = "store", default = 1, type = 'integer', 
                help = "the index of cluster to process from [1, 2, 3]"),
    make_option(c("-n", "--n_iteration"), 
                action = "store", default = 10, type = 'integer',
                help = "No. of time to iteration [default %default]"))
opt <- parse_args(OptionParser(option_list = option_list))
n_iteration <- opt$n_iteration
cluster_id <- opt$cluster_id

# set up Julia
# julia_setup(installJulia = TRUE)
julia_install_package_if_needed("Circuitscape")
julia_library("Circuitscape")
config <- read.ini(file.path(cnt_path, "circuitscape_setting_template.ini"))

# Make a folder for each cluster
# Make a folder for each cluster
clt_dir <- file.path(cnt_path, sprintf("cluster_%s", cluster_id))
dir.create(clt_dir)
clt_out_dir <- file.path(cnt_exp_path, sprintf("cluster_%s", cluster_id))
dir.create(clt_out_dir)

cluster_msk <- pas %>% filter(cluster == cluster_id) %>% 
    st_union() %>% vect() %>% ext()

# Different conditions
types <- c("mnt_block", "human_block", "human_block_complete")
walk(1:length(types), function(type_index) {
    type <- types[type_index]
    
    # Make folder for each condition
    sub_dir <- file.path(clt_dir, type)
    dir.create(sub_dir)
    out_dir <- file.path(clt_out_dir, type)
    dir.create(out_dir)
    
    # Subset suitability
    suit_map <- file.path(
        cnt_path, sprintf("suitability_%s_spath.asc", type)) %>% 
        rast() %>% 
        crop(cluster_msk)
    suit_map[is.na(suit_map)] <- -9999
    suit_fname <- file.path(
        sub_dir, sprintf("suit_%s_cluster%s.asc", type, cluster_id))
    writeRaster(suit_map, 
                suit_fname,
                wopt = list(NAflag = -9999))
    
    # Iterations
    nodes <- pas %>% filter(cluster == cluster_id) %>% 
        mutate(id = 1:nrow(.)) %>% select(id)
    
    # Get the possible combination of non-touched PAs
    remove_neighbor <- function(orig_node, plys){
        ids_ex <- unique(
            unlist(st_intersects(
                orig_node %>% st_buffer(0.01), plys)))
        plys %>% slice(-ids_ex) %>% rbind(orig_node)}
    
    unique_index <- function(lst) {
        sapply(1:length(lst), function(id) {
            sapply(setdiff(1:length(lst), id), function(id_2) {
                identical(lst[[id]], lst[[id_2]])
            }) %>% which() %>% c(id) %>% min()
        }) %>% unique()}
    
    sf_use_s2(FALSE)
    set.seed(cluster_id * 100 + 10 * type_index)
    comb_nodes <- mclapply(rep(nodes$id, n_iteration), function(index) {
        ids_ex <- unique(
            unlist(st_intersects(
                nodes %>% filter(id == index) %>% st_buffer(0.01), nodes)))
        nodes_select <- nodes %>% slice(-ids_ex) %>% 
            rbind(nodes %>% filter(id == index))
        
        while (sum(st_overlaps(nodes_select %>% st_buffer(0.01), 
                               sparse = FALSE) * 1) > 0) {
            start_node <- nodes_select %>% sample_n(size = 1)
            nodes_select <- remove_neighbor(start_node, nodes_select)
        }
        
        # Return
        nodes_select
    }, mc.cores = detectCores())
    
    ids <- unique_index(comb_nodes)
    comb_nodes <- comb_nodes[ids]
    
    cum_rasters <- lapply(1:5, function(n) {
        nodes <- comb_nodes[[n]] %>% rasterize(suit_map, field = "id")
        nodes[nodes == 0] <- -9999
        
        # Name it
        nm <- sprintf('cluster%s_%s_%s', cluster_id, type, n)
        message(sprintf("Run experiment for the case: %s.", nm))
        
        # Save out nodes
        node_name <- file.path(sub_dir, sprintf("nodes_%s.asc", nm))
        writeRaster(nodes, node_name, wopt = list(NAflag = -9999))
        
        # Revise the config parameters
        config$`Output options`$write_cur_maps <- 0
        config$`Output options`$write_volt_maps <- 0
        config$`Output options`$write_cum_cur_map_only <- "True"
        config$`Output options`$write_max_cur_maps <- "False"
        config$`Output options`$output_file <-
            file.path(out_dir, nm)
        config$`Logging Options`$log_file <-
            file.path(out_dir, sprintf("%s.log", nm))
        config$`Options for pairwise and one-to-all and all-to-one modes`$point_file <-
            node_name
        config$`Habitat raster or graph`$habitat_file <- suit_fname
        fname <- file.path(sub_dir, sprintf("%s.ini", nm))
        write.ini(config, fname)
        
        # Run the experiment
        if (file.exists(fname)) {
            julia_command(sprintf('compute("%s")', fname))
        } else {
            warning("Failed to write out config file.")
        }
        
        # Read the files
        cum_raster <- rast(
            file.path(out_dir, sprintf("%s_cum_curmap.asc", nm)))
        mask(cum_raster, nodes, inverse = TRUE, updatevalue = NA)
    })
    
    # Gather results
    cum_raster <- do.call(c, cum_rasters) %>% 
        sum(na.rm = TRUE)
    writeRaster(
        cum_raster, 
        file.path(out_dir, 
                  sprintf("cum_curmap_%s_cluster%s.tif", type, cluster_id)))
})
