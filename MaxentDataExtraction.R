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
i = 1

for (i in 1:nrow(sample.df)){
   species <- sample.df[i, 1] #ith species
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
   #species.pts <- over(SpatialPoints(species.pts.thinned, proj4string=gbif.prj), species.pts.orig)
   species.pts.spdf <- SpatialPointsDataFrame(coords = data.frame(species.pts.thinned[,1], species.pts.thinned[,2]),
                                              proj4string=gbif.prj, data= data.frame("x" = species.pts.thinned[,1],
                                                                                     "y" = species.pts.thinned[,2],
                                                                                     "ID" = 1:nrow(species.pts.thinned)))
   
   # get information about ecoregions that each point lies within, creating a buffer of
   # 10000 to avoid getting NAs because of coarse borders and islands
   print("Extracting...")
   sp.over.df <-as.data.frame(extract(y = species.pts.spdf@data[,1:2], proj4string=gbif.prj, 
                                      buffer = 20000, x = ecoregion.raster, small = TRUE, df = TRUE, 
                                      fun = function(x){max(x)}))
   sp.over <- merge(x = sp.over.df, y = ecoregions.df, by.x = "layer", by.y = "ID", all.x = TRUE)     
   
   
   print("Extract 1 Finished")

   # Removing points from biomes that have less tat 5 pts AND less than
   # 1 percent of observations from data
   species.pts.spdf@data <- base::merge(y = sp.over, x = species.pts.spdf@data, all.x = TRUE)
   species.pts.spdf <- species.pts.spdf[!is.na(species.pts.spdf$region.names),]
   remove.biomes <- factor(levels=(levels(species.pts.spdf$biome.names)))
   j = 1
   for (i in 1:length(biomes)){
      biome <- unique(species.pts.spdf$biome.names)[i]
      if (nrow(species.pts.spdf[species.pts.spdf$biome.names==biome,]) < 5 & 
          nrow(species.pts.spdf[species.pts.spdf$biome.names==biome,])/
          nrow(species.pts.spdf[species.pts.spdf$biome.names,]) < 0.01){
         remove.biomes[j] <- biome
         j = j + 1
      }
   }
   if (length(remove.biomes) > 0){
      species.pts.spdf.rm <- species.pts.spdf[!species.pts.spdf$biome.names %in% remove.biomes,]
      
   }
   
   regions <- unique(species.pts.spdf.rm$region.names)

   # clipping the ecoregions layer so we only have the regions/biomes
   # where the species has been found
   these.regions.spdf <- ecoregion[as.vector(!is.na(match(ecoregion$region,regions))),]
   # generate random points, same number as observations
   print("Generating random points...")
   random.pts <- spsample(these.regions.spdf, n=10000, type="random")
   # I want to save this info now
   species.pts.df <- as.data.frame(species.pts.spdf@data)
   write.csv(species.pts.df, paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                   gsub(" ", ".", species), ".editeddata.csv", sep = ""))
   print("Extracting climate data...")
   
   # making a two simple data frames for maxent; region and biome
   simpledata.df <- cbind(rep(1, nrow(species.pts.df)),species.pts.df[,1:2])
   colnames(simpledata.df) <- c("resp", "x", "y")
   
   random.df <- data.frame("x"=random.pts@coords[,1], "y"=random.pts@coords[,2],"resp"=0)
   for.maxent.pts <- rbind(simpledata.df, random.df)

   # Extract climate data
   
   df.clim <- as.data.frame(extract(clim, cbind(for.maxent.pts$x, for.maxent.pts$y)))
   print(paste(as.character(species), "climate extraction completed"))
   # gotta rename the climate variable names for the columns. 
   bionames <- c("annualmeantemp", "meandiurnalrange", "isothermality", "tempseasonality",
                 "maxtempwarmestmonth", "mintempcoldestmonth", "tempannualrange",
                 "meantempwettestquart", "meantempdriestquart", 
                 "meantempwarmestquart", "meantempcoldestquart",
                 "annualprecip", "precipwettestmo", "precipdriestmo", "precipseasonality",
                 "precipwettestquart", "precipdriestquart", "precipwarmestquart", 
                 "precipcoldestquart")
   colnames(df.clim) <- bionames
   # Binds climate data and long lat and pres together, in a way that the
   # maxent function likes
   for.maxent.full <- cbind(for.maxent.pts, df.clim)
   # Save 'em!
   write.csv(for.maxent.full, paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                           gsub(" ", ".", species), ".maxentdatafile.csv", sep = ""))
   
}
