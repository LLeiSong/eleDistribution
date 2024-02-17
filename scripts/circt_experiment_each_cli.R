# Title     : Calculate landscape connectivity for each habitat cluster
# Objective : The main function to calculate the species level habitat 
#             connectivity. Use the mean suitability values as the conductance. 
# Created by: Lei Song
# Created on: 09/09/23
# Note      : Use cg+amg solver to avoid RAM crisis. The user can modify the
#             loop part to an independent script to run on HPC to speed up.
#             The users should install and set up Circuitscape before use it.

# Load libraries
library(sf)
library(here)
library(terra)
library(dplyr)
library(stringr)
library(ini)
library(purrr)
library(JuliaCall)
library(parallel)
library(optparse)
sf_use_s2(FALSE)

# Set directories
data_path <- here("data")
result_path <- here("results")
vct_path <- file.path(data_path, "vectors")
cnt_path <- file.path(data_path, "connectivity")
cnt_exp_path <- file.path(result_path, "connectivity")
if (!dir.exists(cnt_path)) dir.create(cnt_path)
if (!dir.exists(cnt_exp_path)) dir.create(cnt_exp_path)

# Parse inline arguments
option_list <- list(
    make_option(c("-i", "--initial"), 
                action = "store_true", default = FALSE,
                help = "Initialize inputs [default]"),
    make_option(c("-c", "--cluster_id"), 
                action = "store", default = 1, type = 'integer', 
                help = sprintf("The index of cluster to process ",
                               "from [1, 2, 3]. [default %default]")),
    make_option(c("-n", "--n_iteration"), 
                action = "store", default = 10, type = 'integer',
                help = "No. of time to iteration [default %default]"))
opt <- parse_args(OptionParser(option_list = option_list))
initial <- opt$initial
cluster_id <- opt$cluster_id
n_iteration <- opt$n_iteration

if (initial){
    # Prepare inputs for circuitscape
    ## Filter the protected areas
    ## Remove culture and forest based PAs by assuming dense forest reserve highly 
    ## impossibly serve as a long-term home range as they are mainly savanna animals.
    ## These areas can still be used by animals if they are important, 
    ## just not be used for nodes in circuitscape.
    pas <- st_read(file.path(cnt_path, "pas_selected.geojson"))
    habitat_clusters <- st_read(file.path(cnt_path, "habitat_clusters.geojson"))
    pas <- st_join(pas, habitat_clusters) %>% filter(cluster == cluster_id)
    st_write(pas, file.path(cnt_path, sprintf("pas_%s.geojson", cluster_id)))
    
    ## Write out asc file for suitability
    suitability <- rast(
        file.path(result_path, "landscape_utility_1km_integrated.tif"),
        lyrs = 1) %>% crop(st_union(pas) %>% vect() %>% ext())
    suitability[is.na(suitability)] <- -9999
    fname <- file.path(cnt_path, sprintf('suitability_%s.asc', cluster_id))
    writeRaster(suitability, fname, wopt = list(NAflag = -9999))
}

# Run simulations
## Reload the data
pas <- st_read(file.path(cnt_path, sprintf("pas_%s.geojson", cluster_id))) %>% 
    select()
suit_fname <- file.path(cnt_path, sprintf('suitability_%s.asc', cluster_id))

# set up Julia
# julia_setup(installJulia = TRUE)
julia_install_package_if_needed("Circuitscape")
julia_library("Circuitscape")
# Copy the file to the directory
config <- read.ini(
    file.path(cnt_path, "circuitscape_setting_template.ini"))

# Make folder to save nodes and simulations
src_dir <- file.path(cnt_path, sprintf("nodes_%s", cluster_id))
dst_dir <- file.path(cnt_path, sprintf("simulations_%s", cluster_id))
for (dir_to in c(src_dir, dst_dir)){
    if (!dir.exists(dir_to)) dir.create(dir_to)}

# Iterations
suit_map <- rast(suit_fname)
cum_rasters <- lapply(1:n_iteration, function(n) {
    set.seed(100 + n)
    
    # Get the PAs
    ## Use half PAs and randomly generate mini habitat patches as nodes
    nodes <- pas %>% 
        st_sample(size = rep(1, nrow(.))) %>% 
        st_as_sf() %>% mutate(id = 1:nrow(.)) %>% vect() %>% 
        buffer(width = 5000) %>% 
        rasterize(suit_map, field = "id")
    nodes[nodes == 0] <- -9999
    nodes <- mask(nodes, suit_map)
    
    # Name it
    nm <- sprintf('iter_%s', n)
    message(sprintf("Run experiment for iteration No.%s.", n))
    
    # Save out nodes
    node_name <- file.path(src_dir, sprintf("nodes_%s.asc", nm))
    writeRaster(nodes, node_name, wopt = list(NAflag = -9999), overwrite = TRUE)
    
    # Revise the config parameters
    config$`Output options`$write_cur_maps <- 0
    config$`Output options`$write_volt_maps <- 0
    config$`Output options`$write_cum_cur_map_only <- "True"
    config$`Output options`$write_max_cur_maps <- "False"
    config$`Output options`$output_file <-
        file.path(dst_dir, nm)
    config$`Logging Options`$log_file <-
        file.path(dst_dir, sprintf("%s.log", nm))
    config$`Options for pairwise and one-to-all and all-to-one modes`$point_file <-
        node_name
    config$`Habitat raster or graph`$habitat_file <- suit_fname
    config$`Calculation options`$max_parallel <- parallel::detectCores() - 2
    config$`Calculation options`$solver <- "cg+amg"
    config$`Calculation options`$use_64bit_indexing = 'True'
    config$`Calculation options`$print_timings <- 'True'
    config$`Calculation options`$print_rusages <- 'True'
    config$`Calculation options`$preemptive_memory_release <- 'True'
    fname <- file.path(src_dir, sprintf("%s.ini", nm))
    write.ini(config, fname)
    
    # Run the experiment
    if (file.exists(fname)) {
        julia_command(sprintf('compute("%s")', fname))
    } else {
        warning("Failed to write out config file.")
    }
    
    # Read the files
    cum_raster <- rast(
        file.path(dst_dir, sprintf("%s_cum_curmap.asc", nm)))
    mask(cum_raster, nodes, inverse = TRUE, updatevalue = NA)
})

# Gather results
mean_raster <- do.call(c, cum_rasters) %>% 
    mean(na.rm = TRUE)
mean_fname <- file.path(cnt_exp_path, sprintf("mean_curmap_%s.tif", cluster_id))
writeRaster(mean_raster, mean_fname)
