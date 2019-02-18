# Take large dataset on many ecoregions
# pull out the larger ecoregions (n = 8)
# and smaller vegetative classes, make them into their own polygons
# save them

wwfecoregions <- readOGR(dsn = "./SpatialFiles/", "tnc_terr_ecoregions")

wwfecoregions$WWF_REALM2

biomes <- unionSpatialPolygons(ecoregions, ecoregions@data$WWF_REALM2)
biomes2 <- biomes
pid <- sapply(slot(biomes2, "polygons"), function(x) slot(x, "ID"))
# Create dataframe with correct rownames
p.df <- data.frame( ID=1:length(biomes2),row.names = pid)  
biom.vec <- vector(length = 8)
for (i in 1:nrow(p.df)){
   biom.vec[i] <- biomes@polygons[[p.df$ID[i]]]@ID
}
p.df$biome <- biom.vec
biomes2 <- gSimplify(biomes, .1)
biomes2 <- SpatialPolygonsDataFrame(biomes2, p.df)
writeOGR(obj = biomes2 , dsn = "./SpatialFiles/", layer = "Lg_Ecoregions", driver = "ESRI Shapefile",
         overwrite = TRUE)

regions <- unionSpatialPolygons(ecoregions, ecoregions@data$RealmMHT)
over(regions2, biomes2)
plot(regions)
rid <- sapply(slot(regions, "polygons"), function(x) slot(x, "ID"))
# Create dataframe with correct rownames
r.df <- data.frame( ID=1:length(regions),row.names = rid)  
region.vec <- vector()
for (i in 1:nrow(r.df)){
   region.vec[i] <- regions@polygons[[r.df$ID[i]]]@ID
}
r.df$region <- region.vec

regions2 <- gSimplify(regions, .1)
plot(regions2)
regions.spdf <- SpatialPolygonsDataFrame(regions, r.df)
writeOGR(obj = regions.spdf , dsn = "./SpatialFiles/", layer = "Sm_Ecoregions", driver = "ESRI Shapefile",
         overwrite = TRUE)

# allregions[] <- lapply(ecoregion$region, as.character)
# ecoregion$biome <- substr(unlist(allregions, use.names = FALSE), 1, 2)

# AA AN AT IM NA NT OC PA
# Afrotropic Antarctic Australasia Indo-Malay Nearctic Neotropic Oceania Palearctic

reference.df <- data.frame("Cont.Abbr" = c("AA", "AN", "AT", "IM", "NA", "NT", "OC", "PA"),
                           "Cont." = c("Afrotropic", "Antarctic", "Australasia", "Indo-Malay", "Nearctic",
                                       "Neotropic", "Oceania", "Palearctic"))

require(raster)
require(rgdal)
require(maptools)
# Making the regions a raster
ecoregion <- readOGR(dsn = "./SpatialFiles/", "Sm_Ecoregions")
blank.raster <- raster(ext = extent(ecoregion), crs = crs(ecoregion), resolution = res(clim))
ecoregion.raster <- rasterize(as(ecoregion, "SpatialLines"), blank.raster, progress = "text")
ecoregion.raster <- rasterize(ecoregion, ecoregion.raster, progress =  "text",
                              update = TRUE)
ecoregion.raster.rat <- ratify(ecoregion.raster)
rat <- levels(ecoregion.raster.rat)[[1]]
rat$region.names <- ecoregion$region
rat$biome.names <- ecoregion$biome
levels(ecoregion.raster.rat) <- rat
writeRaster(ecoregion.raster.rat, "ecoregionraster10m", overwrite = TRUE)
levels(ecoregion.raster.rat)
