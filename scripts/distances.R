# Title: Function to calculate distance-based features
# Created by: Lei Song
# Created on: 09/09/23
# Note      : The script uses GRASS GIS, so the users need to install and set up
#             GRASS GIS before use it.

distances <- function(data_path, vct_path, cs_path, scale) {
    ## Use GRASS GIS for speed
    library(sf)
    library(terra)
    library(rgrass7)
    
    ## set up
    gisBase <- '/Applications/GRASS-8.2.app/Contents/Resources'
    crs_mer <- crs(rast(file.path(cs_path, sprintf("bios_%s.tif", scale))), 
                   proj = T)
    initGRASS(gisBase = gisBase,
              home = tempdir(),
              gisDbase = tempdir(),  
              mapset = 'PERMANENT', 
              location = 'osm', 
              override = TRUE)
    execGRASS("g.proj", flags = "c", 
              proj4 = crs_mer)
    execGRASS('r.in.gdal', flags = c("o", "overwrite"),
              input = file.path(cs_path, sprintf("bios_%s.tif", scale)),
              band = 1,
              output = "template")
    execGRASS("g.region", raster = "template")
    
    # Big roads
    big_roads <- read_sf(file.path(vct_path, "big_roads.geojson"))
    writeVECT(big_roads, 'big_roads', v.in.ogr_flags = 'overwrite')
    rm(big_roads)
    execGRASS('v.to.rast', flags = c("overwrite"),
              parameters = list(input = 'big_roads', 
                                output = 'big_roads',
                                use = 'val',
                                value = 1))
    execGRASS('r.grow.distance', flags = c("overwrite"),
              parameters = list(input = 'big_roads', 
                                distance = 'big_roads'))
    execGRASS('r.out.gdal', flags = c("m", "overwrite"),
              output = file.path(cs_path, 
                                 sprintf("dist_to_big_roads_%s.tif", 
                                         scale)),
              input = "big_roads")
    
    # Waterbodies
    waterbodies <- read_sf(file.path(vct_path, "waterbodies.geojson"))
    writeVECT(waterbodies, 'waterbodies', v.in.ogr_flags = 'overwrite')
    rm(waterbodies)
    execGRASS('v.to.rast', flags = c("overwrite"),
              parameters = list(input = 'waterbodies', 
                                output = 'waterbodies',
                                use = 'val',
                                value = 1))
    execGRASS('r.grow.distance', flags = c("overwrite"),
              parameters = list(input = 'waterbodies', 
                                distance = 'waterbodies'))
    execGRASS('r.out.gdal', flags = c("m", "overwrite"),
              output = file.path(cs_path, 
                                 sprintf("dist_to_waterbodies_%s.tif",
                                         scale)),
              input = "waterbodies")
    
    # Rivers
    rivers <- read_sf(file.path(vct_path, "rivers.geojson"))
    writeVECT(rivers, 'rivers', v.in.ogr_flags = 'overwrite')
    rm(rivers)
    execGRASS('v.to.rast', flags = c("overwrite"),
              parameters = list(input = 'rivers', 
                                output = 'rivers',
                                use = 'val',
                                value = 1))
    execGRASS('r.grow.distance', flags = c("overwrite"),
              parameters = list(input = 'rivers', 
                                distance = 'rivers'))
    execGRASS('r.out.gdal', flags = c("m", "overwrite"),
              output = file.path(cs_path, 
                                 sprintf("dist_to_rivers_%s.tif",
                                         scale)),
              input = "rivers")
    
    # Streams
    streams <- read_sf(file.path(vct_path, "streams.geojson"))
    writeVECT(streams, 'streams', v.in.ogr_flags = 'overwrite')
    rm(streams)
    execGRASS('v.to.rast', flags = c("overwrite"),
              parameters = list(input = 'streams', 
                                output = 'streams',
                                use = 'val',
                                value = 1))
    execGRASS('r.grow.distance', flags = c("overwrite"),
              parameters = list(input = 'streams', 
                                distance = 'streams'))
    execGRASS('r.out.gdal', flags = c("m", "overwrite"),
              output = file.path(cs_path, 
                                 sprintf("dist_to_streams_%s.tif",
                                         scale)),
              input = "streams")
    
    # Cropland
    cropland <- rast(file.path(data_path, "landcover/landcover_1000m.tif"))
    cropland[cropland != 1] <- NA
    fname <- tempfile(fileext = ".tif")
    cropland <- resample(
        cropland, rast(file.path(cs_path, sprintf("bios_%s.tif", scale))), 
        method = "near", filename = fname)
    execGRASS('r.in.gdal', flags = c("o", "overwrite"),
              input = fname,
              band = 1,
              output = "cropland")
    execGRASS('r.grow.distance', flags = c("overwrite"),
              parameters = list(input = 'cropland', 
                                distance = 'cropland'))
    execGRASS('r.out.gdal', flags = c("m", "overwrite"),
              output = file.path(cs_path, 
                                 sprintf("dist_to_crop_%s.tif",
                                         scale)),
              input = "cropland")
}
