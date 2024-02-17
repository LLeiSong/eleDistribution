# Title     : Extract possible shortest paths from current density map
# Objective : To get some example paths from the current density map generated
#             by circuitscape. 
# Note      : The script uses GRASS GIS, so the users need to install and set
#             up GRASS GIS before use it.
# Created by: Lei Song
# Created on: 09/09/23

# Load libraries
library(sf)
library(terra)
library(gdistance)
library(pbmcapply)
library(rgrass)
library(sfnetworks)
library(rmapshaper)

# Define the function to convert raster to polylines
raster_to_lines <- function(src_path, 
                            gisBase = "/opt/local/lib/grass82",
                            gisDbase = NULL){
    if (is.null(gisDbase)) gisDbase <- tempdir()
    template <- rast(src_path)
    initGRASS(gisBase = gisBase,
              home = gisDbase,
              gisDbase = gisDbase,
              SG = template,
              mapset = "PERMANENT",
              location = "r_to_l",
              override = TRUE)
    # Import the lcps
    execGRASS("r.in.gdal",
              flags = c("o", "overwrite"),
              input = src_path,
              band = 1,
              output = "lcps")
    execGRASS("g.region", raster = "lcps")
    # Thin the raster
    execGRASS('r.thin', flags = c("overwrite"),
              parameters = list(input = 'lcps', 
                                output = 'lcps_thin',
                                iterations = 1000))
    # Convert to lines
    execGRASS('r.to.vect', flags = c("overwrite"),
              parameters = list(input = 'lcps_thin', 
                                output = 'lcps',
                                type = 'line'))
    # Return the result as a sf vector
    read_VECT("lcps") %>% st_as_sf()
}

# Set directories
pat_path <- file.path("data/vars_patch")
cnt_path <- file.path("data/connectivity")
cnt_exp_path <- file.path("results/connectivity")

# Read the "cost" to calculate the LSP
habitat_clusters <- st_read(file.path(cnt_path, "habitat_clusters.geojson"))
curmap <- rast(file.path(cnt_path, "suitability.asc"))
curmap[is.na(curmap)] <- 0
tr <- transition(raster(curmap), transitionFunction = mean, directions = 16)

borders <- st_cast(habitat_clusters,"LINESTRING")
n_iteration <- 1000
lcps <- do.call(rbind, pbmclapply(1:n_iteration, function(n){
    seed <- 100 + n
    # Get the random coordinates along the borders
    set.seed(seed)
    link_points <- st_sample(borders, size = rep(1, 3))
    coords <- st_coordinates(link_points)
    
    # Loop over to calculate the distance
    pairs <- list(c(1, 2), c(1, 3), c(2, 3))
    names(pairs) <- c("1to2", "1to3", "2to3") # Not consider asymmetric
    do.call(rbind, lapply(1:length(pairs), function(ind){
        pair <- pairs[[ind]]; nm <- names(pairs)[ind]
        coords_pair <- as.matrix(coords[pair, c("X", "Y")])
        tryCatch({
            lcp <- shortestPath(
                tr, coords_pair[1, ], coords_pair[2, ], 
                output = "SpatialLines")
            st_as_sf(lcp) %>% mutate(path = nm) %>% st_set_crs(4326)},
            error = function(e) NULL)
    }))
}, mc.cores = 10))

# Save out the raw file
st_write(lcps, file.path(cnt_exp_path, "least_cost_paths.geojson"))

# Clean up
vars <- rast(file.path(pat_path, "variables.tif"))
cropland <- subset(vars, "cropland_ratio_3")
cropland[cropland < 0.97] <- NA
lcp_accum <- rasterize(lcps, curmap, fun = 'count') %>% 
    mask(habitat_clusters, inverse = TRUE, updatevalue = NA)
lcps_all <- lcp_accum > 50
lcps_all[isFALSE(lcps_all)] <- NA
lcps_all <- mask(lcps_all, cropland, inverse = TRUE, updatevalue = NA)
writeRaster(lcps_all, file.path(cnt_exp_path, "lcps_all.tif"), 
            datatype = "INT1U", overwrite = TRUE)
lcps_common <- lcp_accum > 200
lcps_common[isFALSE(lcps_common)] <- NA
lcps_common <- mask(lcps_common, cropland, inverse = TRUE, updatevalue = NA)
writeRaster(lcps_common, file.path(cnt_exp_path, "lcps_common.tif"), 
            datatype = "INT1U", overwrite = TRUE)

# Convert to lines
gisBase <- "/opt/local/lib/grass82"
gisDbase <- here("data/grass")
src_path <- file.path(cnt_exp_path, "lcps_all.tif")
src_path_common <- file.path(cnt_exp_path, "lcps_common.tif")

# Get the polylines
lcps_all <- raster_to_lines(src_path, gisBase, gisDbase)
lcps_common <- raster_to_lines(src_path_common, gisBase, gisDbase)

# Merge the LINESTRINGs for all
com <- components(graph.adjlist(st_touches(lcps_all)))
lcps_all <- st_sf(geometry = lcps_all)
lcps_all$group <- com$membership
lcps_all <- lcps_all |> group_by(group) |>
    summarise() %>% st_cast("MULTILINESTRING")

# Merge the LINESTRINGs for common
com <- components(graph.adjlist(st_touches(lcps_common)))
lcps_common <- st_sf(geometry = lcps_common)
lcps_common$group <- com$membership
lcps_common <- lcps_common |> group_by(group) |>
    summarise() %>% st_cast("MULTILINESTRING") %>% 
    mutate(length = st_length(.)) %>% 
    filter(length > 10000 %>% units::set_units("m"))

# Select the lcps_all using lcps_common, the objective is to select the whole lines
lcps <- lcps_all %>% slice(unique(unlist(st_intersects(lcps_common, .)))) %>% 
    ms_simplify(0.05)

# Delete temporary files
file.remove(src_path, src_path_common)

# Now prune the selected lines
borders_buf <- habitat_clusters %>% 
    st_buffer(res(curmap)[1] * 1.5) %>% st_cast("LINESTRING")
lcps_simplified <- do.call(rbind, lapply(lcps$group, function(ind){
    lcp <- lcps %>% filter(group == ind)
    net <- as_sfnetwork(lcp %>% st_cast("LINESTRING"), directed = FALSE) %>% 
        activate("edges")
    
    # Get points
    pts <- st_intersection(borders_buf, lcp) %>% 
        st_cast("MULTIPOINT") %>% st_cast("POINT")
    
    clusters <- unique(pts$cluster)
    pts_cluster1 <- pts %>% filter(cluster == clusters[1])
    pts_cluster2 <- pts %>% filter(cluster == clusters[2])
    
    do.call(
        rbind, lapply(1:nrow(pts_cluster1), function(m){
        do.call(rbind, lapply(1:nrow(pts_cluster2), function(n){
            paths <- st_network_paths(
                net, from = pts_cluster1[m, ], to = pts_cluster2[n, ])
            ids <- paths %>% slice(1) %>% pull(edge_paths) %>% unlist()
            st_as_sf(slice(activate(net, "edges"), ids), "edges") %>% 
                st_cast("MULTILINESTRING") |> group_by(group) |>
                summarise() %>% st_cast("MULTILINESTRING")
        }))
    })) %>% mutate(length = st_length(.)) %>% 
        filter(length == min(length)) %>% select(-length)
}))

st_write(lcps_simplified, file.path(cnt_exp_path, "least_cost_examples.geojson"))
