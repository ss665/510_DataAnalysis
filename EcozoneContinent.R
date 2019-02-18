# Sarah Schooler
# 2/16/19
# Continent vs. Ecozone comparisons
# 10 species
require(rgdal)
require(rgeos)
require(raster)
require(enmSdm)
require(dismo)
species.df.full <- read.csv("Species_Table_021619.csv", row.names = 1, stringsAsFactors = FALSE)
species.df <- read.csv("simple.species.csv", row.names = 1, stringsAsFactors = FALSE)

# Sample DF
sample.df <- species.df[c(154, 376, 393, 19, 409),]
sample.df

# Read get species points, buffer them based around continents and ecoregions.
# This loop will also create the maxent data frame, so we need climate as well.
ecoregion <- readOGR(dsn = "./SpatialFiles/", "Sm_Ecoregions")
ecoregion.raster <- raster("ecoregionraster10m.gri")
ecoregions.df <- levels(ecoregion.raster)[[1]]

clim <- brick("clim10m")
require(raster)
require(dplyr)
for (i in 1:nrow(sample.df)){
   species <- sample.df[i, 1] #ith species
   print(species)
   print(i)
   # spatial points file that I made earlier
   species.pts.spdf <- readOGR(dsn = paste0("./Species_Pts/", gsub(" ", ".", species), "/"), 
                                         layer = paste0(gsub(" ", ".", species),"sp.pts"))
   #species.pts.thinned <- elimCellDups(SpatialPoints(coordinates(species.pts.orig), 
   #                                               proj4string = crs(clim)), clim[[5]])
   #species.pts <- over(species.pts.thinned, species.pts.orig)
   #species.pts.spdf <- SpatialPointsDataFrame(coords = data.frame(species.pts$dcmlLng, species.pts$dcmlLtt),
   #                                          proj4string=gbif.prj, data=species.pts)
   # get information about ecoregions that each point lies within, creating a buffer of
   # 10000 to avoid getting NAs because of coarse borders and islands
   
   #sp.buffer <- buffer(species.pts.spdf, width = 10000, dissolve = FALSE, progress= TRUE)
   #print("Buffered, overlapping")
   sp.over.df <- as.data.frame(extract(y = SpatialPoints(coordinates(species.pts.spdf), proj4string=gbif.prj), 
                      buffer = 200000, x = ecoregion.raster, small = TRUE, df = TRUE,
                      progress = TRUE))
   
   sp.over.sort <- sp.over.df[order(sp.over.df$layer),]
   sp.over <- sp.over.sort[!duplicated(sp.over.sort$ID),]
   sp.over <- merge(sp.over, ecoregions.df[,2:3], all.x = TRUE)
   #sp.over <- over(sp.buffer, ecoregion[,c("region", "biome")])
   print("Extract 1 Finished")
   # need to do this stuff to keep the data frame nice
   species.pts.spdf$ecoregion <- sp.over$region.names
   species.pts.spdf$biome <- sp.over$biome.names
   # which regions/biomes do the species reside in?
   regions <- unique(species.pts.spdf$ecoregion)
   biomes <- unique(species.pts.spdf$biome)
   # clipping the ecoregions layer so we only have the regions/biomes
   # where the species has been found
   these.regions.spdf <- ecoregion[as.vector(!is.na(match(ecoregion$region,regions))),]
   these.biomes.spdf <- ecoregion[as.vector(!is.na(match(ecoregion$biome,biomes))),]
   # generate random points, same number as observations
   print("Generating random points")
   random.pts.region <- spsample(these.regions.spdf, n=nrow(species.pts.spdf), type="random")
   random.pts.biome <- spsample(these.biomes.spdf, n=nrow(species.pts.spdf), type = "random")
   
   # I want to save this info now
   species.pts.df <- as.data.frame(species.pts.spdf@data)
   write.csv(species.pts.df, paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                 gsub(" ", ".", species), ".editeddata.csv", sep = ""))
   print("Extracting...")
   
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

}


# Want to compare temperatures at occupied and random locations
par(mfrow = c(3, 2))
for (i in 1:nrow(sample.df)){
   species <- sample.df[i, 1]
   for.maxent.full.region <- read.csv(paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                            gsub(" ", ".", species), ".maxentdatafile.region.csv", sep = ""), row.names = 1)
   for.maxent.full.biomes <- read.csv(paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                            gsub(" ", ".", species), ".maxentdatafile.biome.csv", sep = ""), row.names = 1)
   compare <- data.frame("Region.Maxtemp" = for.maxent.full.region$maxtempwarmestmonth,
                         "Biome.Maxtemp" = for.maxent.full.biomes$maxtempwarmestmonth,
                         "Region.Mintemp" = for.maxent.full.region$mintempcoldestmonth,
                         "Biome.Mintemp" = for.maxent.full.biomes$mintempcoldestmonth)
   boxplot.matrix(as.matrix(compare[,1:2]), main = paste0(sample.df[i, 4], ", ", species), 
                  sub = sample.df[i, 3])
   
}

par(mfrow = c(3, 2))
for (i in 1:nrow(sample.df)){
   species <- sample.df[i, 1]
   for.maxent.full.region <- read.csv(paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                            gsub(" ", ".", species), ".maxentdatafile.region.csv", sep = ""), row.names = 1)
   for.maxent.full.biomes <- read.csv(paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                            gsub(" ", ".", species), ".maxentdatafile.biome.csv", sep = ""), row.names = 1)
   compare <- data.frame("Region.Maxtemp" = for.maxent.full.region$maxtempwarmestmonth,
                         "Biome.Maxtemp" = for.maxent.full.biomes$maxtempwarmestmonth,
                         "Region.Mintemp" = for.maxent.full.region$mintempcoldestmonth,
                         "Biome.Mintemp" = for.maxent.full.biomes$mintempcoldestmonth)
   boxplot.matrix(as.matrix(compare[,3:4]), main = paste0(sample.df[i, 4], ", ", species), 
                  sub = sample.df[i, 3])
   
}


species.df.full <- read.csv("Species_Table_021619.csv", row.names = 1, stringsAsFactors = FALSE)
species.df <- read.csv("simple.species.csv", row.names = 1, stringsAsFactors = FALSE)
compare <- data.frame("species" = character(), "thresholdRegion" = numeric(), "SuitRegion" = numeric(), 
                      "AreaRegion" = numeric(), "thresholdBiome" = numeric(), "SuitBiome" = numeric(), 
                      "AreaBiome" = numeric())
gbif.prj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


species <- species.df[i, 1]

raster_area <- function(area_raster){
   area_raster[area_raster==0] <- NA
   cell_size <- area(area_raster, na.rm = TRUE, weights = FALSE)
   cell_size <- cell_size[!is.na(cell_size)]
   raster_area <- length(cell_size) * median(cell_size)
   if (length(raster_area) == 0) {return(0)}
   return(raster_area)
}

for (i in 1:sample.df){
   species <- sample.df[i, 1]
   # Run maxent
   # Region first
   for.maxent.full.region <- read.csv(paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                    gsub(" ", ".", species), ".maxentdatafile.region.csv", sep = ""), row.names = 1)
   predsnums <- c(1, 2, 5, 6, 7)
   preds.region <- for.maxent.full.region[,predsnums+3]
   resp.region <- for.maxent.full.region[,"resp"] 
   thisRegMult <- 1
   params <- c(paste("betamultiplier=", thisRegMult, sep=""), "jackknife=false")
   require(rJava)
   maxent.model.region <- maxent(x = preds.region, p = resp.region, path = "./MaxentModels/ScratchDir/",
                       args = params)
   curPredsRaster <- clim[[predsnums]]
   names(curPredsRaster) <- colnames(preds.region)
   dir.create(paste0("./MaxentModels/",gsub(" ", ".", species)))
   sp.raster.region <- predict(model = maxent.model.region, object = curPredsRaster, progress = "text",
                        filename = paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species),
                                       ".full.region.tif"), overwrite = TRUE)
   thresh.region <- maxent.model.region@results[38]
   sp.raster.thresh.region <- sp.raster.region
   sp.raster.thresh.region[sp.raster.thresh.region<thresh.region] = 0
   avg.suit.region <- mean(sp.raster.region[sp.raster.region>0])
   thresh.area.region <- raster_area(sp.raster.thresh.region)
   plot(sp.raster.region)
   plot(sp.raster.thresh.region)
   writeRaster(sp.raster.thresh.region, filename = paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species),
                                                ".region.thresh.tif"), overwrite = TRUE)
   
# Biome
   for.maxent.full.biome <- read.csv(paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                  gsub(" ", ".", species), ".maxentdatafile.biome.csv", sep = ""), row.names = 1)
   preds.biome <- for.maxent.full.biome[,predsnums+3]
   resp.biome <- for.maxent.full.biome[,"resp"] 
   maxent.model.biome <- maxent(x = preds.biome, p = resp.biome, path = "./MaxentModels/ScratchDir/",
                       args = params)
   sp.raster.biome <- predict(model = maxent.model.biome, object = curPredsRaster, progress = "text",
                     filename = paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species),
                                       ".full.biome.tif"), overwrite = TRUE)

   thresh.biome <- maxent.model.biome@results[38]
   sp.raster.thresh.biome <- sp.raster.biome
   sp.raster.thresh.biome[sp.raster.thresh.biome<thresh.biome] = 0
   avg.suit.biome <- mean(sp.raster.biome[sp.raster.biome>0])
   thresh.area.biome <- raster_area(sp.raster.thresh.biome)
   writeRaster(sp.raster.thresh.biome, filename = paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species),
                                                      ".biome.thresh.tif"), overwrite = TRUE)
   this.row <- data.frame("species" = species, "thresholdRegion" = thresh.region, "SuitRegion" = avg.suit.region, 
                       "AreaRegion" = thresh.area.region, "thresholdBiome" = thresh.biome, 
                       "SuitBiome" = avg.suit.biome, "AreaBiome" = thresh.area.biome)

   compare <- rbind(this.row, compare)
   
}



### Comparison method

overlap <- function(raster1, raster2){
   raster1[raster1>0] <- 1
   raster2[raster2>0] <- 1
   add <- raster1+raster2
   add[add<2] <- 0
   return(add)
}
raster_area <- function(area_raster){
   area_raster[area_raster==0] <- NA
   cell_size <- area(area_raster, na.rm = TRUE, weights = FALSE)
   cell_size <- cell_size[!is.na(cell_size)]
   raster_area <- length(cell_size) * median(cell_size)
   if (length(raster_area) == 0) {return(0)}
   return(raster_area)
}



clim <- brick("clim10m")


overlap.df <- data.frame("species" =character(), "Area.B" = double(),
                         "Area.Overlap.Upper.B" = double(), "Area.Overlap.Lower.B" = double(), 
                         "Area.R" = double(), "Area.Overlap.Upper.R" = double(), 
                         "Area.Overlap.Lower.R" = double(), "Area.Upper" = double(),
                         "Area.Lower" = double())
species.df.full <- read.csv("Species_Table_021619.csv", row.names = 1, stringsAsFactors = FALSE)
starttime <-Sys.time()

for (i in 1:nrow(sample.df)){
   print(paste(i, "/", nrow(sample.df)))
   species <- sample.df[i, 1]
   print(species)
   
   upper.crit <- species.df.full$UCT...C.[species.df.full$Species==species]
   lower.crit <-species.df.full$LCT...C.[species.df.full$Species==species]
   
   upper.crit.raster <- clim[[5]] <= upper.crit
   lower.crit.raster <- clim[[6]] >= lower.crit
   data.thresh.region <- raster(paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species),
                                       ".region.thresh.tif"))
   overlap.upper.region <- overlap(data.thresh.region, upper.crit.raster)
   overlap.lower.region <- overlap(data.thresh.region, lower.crit.raster)
   
   data.thresh.biome <- raster(paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species),
                                      ".biome.thresh.tif"))
   overlap.upper.biome <- overlap(data.thresh.biome, upper.crit.raster)
   overlap.lower.biome <- overlap(data.thresh.biome, lower.crit.raster)
   this.row <- data.frame("species" =species, "Area.b" = raster_area(data.thresh.biome)/1000, 
                          "Area.Overlap.Upper.B" = raster_area(overlap.upper.biome)/1000, 
                          "Area.Overlap.Lower.B" = raster_area(overlap.lower.biome)/1000,
                          "Area.R" = raster_area(data.thresh.region)/1000, 
                          "Area.Overlap.Upper.R" = raster_area(overlap.upper.region)/1000, 
                          "Area.Overlap.Lower.R" = raster_area(overlap.lower.region)/1000,
                          "Area.Upper" = raster_area(upper.crit.raster)/1000, 
                          "Area.Lower" = raster_area(lower.crit.raster)/1000)
   print(this.row)
   overlap.df <- rbind(this.row, overlap.df)
   write.csv(overlap.df, "./SummaryTables/areaoverlap.csv")
   #cur.runtime <- Sys.time()-starttime
}


#cur.runtime <- Sys.time()-starttime




### misc code
#region.spdf <- SpatialPointsDataFrame(coords = data.frame(for.maxent.full.region$x, for.maxent.full.region$y),
#                                   proj4string=gbif.prj, data=for.maxent.full.region)

plot(region.spdf, add = TRUE)
