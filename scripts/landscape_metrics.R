library(landscapemetrics)
library(raster)

# List all available metrics
metrics_avail <- list_lsm()

# 1000m
data_path <- here('data')
landcover <- raster(file.path(data_path, 'landcover_1000m.tif'))
check_landscape(landcover)
metrics_1000m <- calculate_lsm(landcover)
write.csv(metrics_1000m, 
          file.path(data_path, 'landscape_metrics_1000m.csv'),
          row.names = FALSE)
save(metrics_1000m, file = file.path(data_path, 'metrics_1000m.rda'))
rm(landcover, metrics_1000m)

# 500m
landcover <- raster(file.path(data_path, 'landcover_500m.tif'))
check_landscape(landcover)
metrics_500m <- calculate_lsm(landcover)
write.csv(metrics_500m, 
          file.path(data_path, 'landscape_metrics_500m.csv'),
          row.names = FALSE)
save(metrics_500m, file = file.path(data_path, 'metrics_500m.rda'))
rm(landcover, metrics_500m)

# 250m
landcover <- raster(file.path(data_path, 'landcover_250m.tif'))
check_landscape(landcover)
metrics_250m <- calculate_lsm(landcover)
write.csv(metrics_250m, 
          file.path(data_path, 'landscape_metrics_250m.csv'),
          row.names = FALSE)
save(metrics_250m, file = file.path(data_path, 'metrics_250m.rda'))
rm(landcover, metrics_250m)

# 100m
landcover <- raster(file.path(data_path, 'landcover_100m.tif'))
check_landscape(landcover)
metrics_100m <- calculate_lsm(landcover)
write.csv(metrics_100m, 
          file.path(data_path, 'landscape_metrics_100m.csv'),
          row.names = FALSE)
save(metrics_100m, file = file.path(data_path, 'metrics_100m.rda'))
rm(landcover, metrics_100m)
