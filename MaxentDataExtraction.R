# Code to extract climate data we need from species points
# Also generate random points


require(rgdal)
require(rgeos)
require(raster)
require(enmSdm)
require(dismo)
require(rJava)
require(dplyr)

species.df.full <- read.csv("Species_Table_021619.csv", row.names = 1, stringsAsFactors = FALSE)
simple.species <- read.csv("simple.species.csv", row.names = 1, stringsAsFactors = FALSE)
gbif.prj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
clim <- brick("clim10m")

species.df <- simple.species[simple.species$nObs<100000,]

# Sample DF
sample.df <- species.df[c(154, 376, 393, 19, 409, 361, 131, 26, 281, 133),]
# Read get species points, buffer them based around continents and ecoregions.
# This loop will also create the maxent data frame, so we need climate as well.
ecoregion <- readOGR(dsn = "./SpatialFiles/", "Sm_Ecoregions")
ecoregion.raster <- raster("ecoregionraster10m.gri")
ecoregions.df <- levels(ecoregion.raster)[[1]]

clim <- brick("clim10m")


for (i in 1:nrow(species.df)){
   species <- species.df[i, 1] #ith species
   print(species)
   print(i)
   # spatial points file that I made earlier
   species.pts.orig <- readOGR(dsn = paste0("./Species_Pts/", gsub(" ", ".", species), "/"), 
                               layer = paste0(gsub(" ", ".", species),"sp.pts"))
   # reference raster for thinning
   print("Thinning...")
   ref.raster <- raster(crs = crs(species.pts.orig), ext = extent(clim), 
                        resolution = c(0.00833333, 0.00833333))
   species.pts.thinned <- elimCellDups(coordinates(species.pts.orig), ref.raster)
   species.pts <- over(SpatialPoints(species.pts.thinned, proj4string=gbif.prj), species.pts.orig)
   species.pts.spdf <- SpatialPointsDataFrame(coords = data.frame(species.pts$dcmlLng, species.pts$dcmlLtt),
                                              proj4string=gbif.prj, data=species.pts)
   
   # get information about ecoregions that each point lies within, creating a buffer of
   # 10000 to avoid getting NAs because of coarse borders and islands
   print("Extracting...")
   sp.over.df <-as.data.frame(extract(y = SpatialPoints(coordinates(species.pts.spdf), proj4string=gbif.prj), 
                                      buffer = 20000, x = ecoregion.raster, small = TRUE, df = TRUE, 
                                      fun = function(x){max(x)}))
   sp.over <- merge(x = sp.over.df, y = ecoregions.df[,1:3], by.x = "layer", by.y = "ID", all.x = TRUE)     
   print("Extract 1 Finished")
   # need to do this stuff to keep the data frame nice
   species.pts.spdf$ecoregion <- sp.over$region
   species.pts.spdf$biome <- sp.over$biome
   # which regions/biomes do the species reside in?
   regions <- unique(species.pts.spdf$ecoregion)
   biomes <- unique(species.pts.spdf$biome)
   # clipping the ecoregions layer so we only have the regions/biomes
   # where the species has been found
   these.regions.spdf <- ecoregion[as.vector(!is.na(match(ecoregion$region,regions))),]
   these.biomes.spdf <- ecoregion[as.vector(!is.na(match(ecoregion$biome,biomes))),]
   # generate random points, same number as observations
   print("Generating random points...")
   random.pts.region <- spsample(these.regions.spdf, n=10000, type="random")
   random.pts.biome <- spsample(these.biomes.spdf, n=10000, type = "random")
   
   # I want to save this info now
   species.pts.df <- as.data.frame(species.pts.spdf@data)
   write.csv(species.pts.df, paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                   gsub(" ", ".", species), ".editeddata.csv", sep = ""))
   print("Extracting climate data...")
   
   # making a two simple data frames for maxent; region and biome
   simpledata.df <- cbind(rep(1, nrow(species.pts.df)),species.pts.df[,1:2])
   colnames(simpledata.df) <- c("resp", "x", "y")
   
   random.df.region <- data.frame("x"=random.pts.region@coords[,1], "y"=random.pts.region@coords[,2],"resp"=0)
   for.maxent.pts.region <- rbind(simpledata.df, random.df.region)
   
   random.df.biome <- data.frame("x"=random.pts.biome@coords[,1], "y"=random.pts.biome@coords[,2],"resp"=0)
   for.maxent.pts.biome <- rbind(simpledata.df, random.df.biome)
   
   # Extract climate data
   
   df.clim.region <- as.data.frame(extract(clim, cbind(for.maxent.pts.region$x, for.maxent.pts.region$y)))
   df.clim.biome <- as.data.frame(extract(clim, cbind(for.maxent.pts.biome$x, for.maxent.pts.biome$y)))
   print(paste(as.character(species), "climate extraction completed"))
   # gotta rename the climate variable names for the columns. 
   bionames <- c("annualmeantemp", "meandiurnalrange", "isothermality", "tempseasonality",
                 "maxtempwarmestmonth", "mintempcoldestmonth", "tempannualrange",
                 "meantempwettestquart", "meantempdriestquart", 
                 "meantempwarmestquart", "meantempcoldestquart",
                 "annualprecip", "precipwettestmo", "precipdriestmo", "precipseasonality",
                 "precipwettestquart", "precipdriestquart", "precipwarmestquart", 
                 "precipcoldestquart")
   colnames(df.clim.region) <- bionames
   colnames(df.clim.biome) <- bionames
   # Binds climate data and long lat and pres together, in a way that the
   # maxent function likes
   for.maxent.full.region <- cbind(for.maxent.pts.region, df.clim.region)
   for.maxent.full.biome <- cbind(for.maxent.pts.biome, df.clim.biome)
   # Save 'em!
   write.csv(for.maxent.full.region, paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                           gsub(" ", ".", species), ".maxentdatafile.region.csv", sep = ""))
   
   write.csv(for.maxent.full.biome, paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                          gsub(" ", ".", species), ".maxentdatafile.biome.csv", sep = ""))
   write.csv(for.maxent.full.biome, paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                          gsub(" ", ".", species), ".maxentdatafile.csv", sep = ""))
}
