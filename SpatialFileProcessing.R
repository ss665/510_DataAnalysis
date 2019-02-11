# Take large dataset on many ecoregions
# pull out the larger ecoregions (n = 8)
# and smaller vegetative classes, make them into their own polygons
# save them

ecoregions <- readOGR(dsn = "./SpatialFiles/", "tnc_terr_ecoregions")

ecoregions$RealmMHT

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
