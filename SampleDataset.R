## Sample datasets
## Maybe will split up later
require(rgdal)
require(rgeos)
require(raster)
require(enmSdm)
require(dismo)
require(rJava)
species.df.full <- read.csv("Species_Table_021619.csv", row.names = 1, stringsAsFactors = FALSE)
species.df <- read.csv("simple.species.csv", row.names = 1, stringsAsFactors = FALSE)
gbif.prj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
clim <- brick("clim10m")
# Sample DF
sample.df <- species.df[c(154, 376, 393, 19, 409, 361, 131, 26, 281, 133),]
# Read get species points, buffer them based around continents and ecoregions.
# This loop will also create the maxent data frame, so we need climate as well.
ecoregion <- readOGR(dsn = "./SpatialFiles/", "Sm_Ecoregions")
ecoregion.raster <- raster("ecoregionraster10m.gri")
ecoregions.df <- levels(ecoregion.raster)[[1]]

clim <- brick("clim10m")
require(raster)
require(dplyr)
require(enmSdm)

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

####################################################################
require(DescTools)

species.df <- read.csv("simple.species.csv", row.names = 1, stringsAsFactors = FALSE)
species.df.full <- read.csv("Species_Table_021619.csv", row.names = 1, stringsAsFactors = FALSE)
df <- data.frame("species" = character(), "full.area" = numeric(), "above.lower.crit" = numeric(), "below.lower.crit" = numeric(),
           "above.upper.crit" = numeric(), "below.upper.crit" = numeric(), 
           "area.between.crits" = numeric())
pts.climate.df <- data.frame("species" = character(), "below.maxtemp" = numeric(), 
                             "above.maxtemp" = numeric(), "below.mintemp" = numeric(), 
                             "above.mintemp" = numeric(), "maxtemp.loc" = numeric(), 
                             "mintemp.loc" = numeric(), "maxcrit" = numeric(), "mincrit" = numeric(),
                             "nobs" = integer())
sample.df <- species.df[c(154, 376, 393, 19, 409, 361, 131, 26, 281, 133),]
# type in biome or region here to get appropriate points
bor = "biome"
for (i in 1:nrow(sample.df)){
   species = sample.df[i, 1]
   print(species)
   for.maxent.full <- read.csv(paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                      gsub(" ", ".", species), ".maxentdatafile.",bor,".csv", sep = ""), row.names = 1)
   upper.crit <- species.df.full$UCT...C.[species.df.full$Species==species]
   lower.crit <-species.df.full$LCT...C.[species.df.full$Species==species]
   predsnums <- c(1, 2, 5, 6, 7)
   preds <- for.maxent.full[,predsnums+3]
   resp <- for.maxent.full[,"resp"] 
   thisRegMult <- 1
   params <- c(paste("betamultiplier=", thisRegMult, sep=""), "jackknife=false")
   dir.create(paste0("./MaxentModels/",gsub(" ", ".", species)))
   maxent.model <- maxent(x = preds, p = resp, path = paste0("./MaxentModels/",gsub(" ", ".", species)),
                          args = params, silent = FALSE)
   
   areas.row <- AreaUnderResponseCurves(x = maxent.model, var = c("annualmeantemp", "maxtempwarmestmonth", "mintempcoldestmonth"), at = median,
                                        expand = 10, data = NULL, fun = predict, upper.crit = upper.crit, lower.crit = lower.crit)
   df <- rbind(cbind("species" = species, areas.row), df)
   write.csv(df, paste0("./SummaryTables/response.suitability.",bor,".csv"))
   
   
   #Evaluating climate at occurence points
   species.pts.climate <- dplyr::filter(for.maxent.full, for.maxent.full$resp == 1)
   species.pts.climate <- species.pts.climate[!is.na(species.pts.climate$maxtempwarmestmonth),]
   below.maxtemp <- sum(species.pts.climate$maxtempwarmestmonth <= upper.crit)
   above.mintemp <- sum(species.pts.climate$mintempcoldestmonth >= lower.crit)
   above.maxtemp <- sum(species.pts.climate$maxtempwarmestmonth > upper.crit)
   below.mintemp <- sum(species.pts.climate$mintempcoldestmonth < lower.crit)
   
   
   maxtemp.loc <- max(species.pts.climate$maxtempwarmestmonth)
   mintemp.loc <- min(species.pts.climate$mintempcoldestmonth)
   meantemp.loc <- mean(species.pts.climate$annualmeantemp)
   nobs <- nrow(species.pts.climate)
   
   pts.climate.df <- rbind(data.frame("species" = species, "below.maxtemp" = below.maxtemp, 
                                      "above.maxtemp" = above.maxtemp, "below.mintemp" = below.mintemp, 
                                      "above.mintemp" = above.mintemp, "maxtemp.loc" = maxtemp.loc, 
                                      "mintemp.loc" = mintemp.loc, "meantemp.loc" = meantemp.loc,
                                      "nobs" = nobs, "maxcrit" = upper.crit, "mincrit" = lower.crit), pts.climate.df)
   write.csv(pts.climate.df, paste0("./SummaryTables/ptssuitability.",bor,".csv"))
}
######################################################################
# species.df <- read.csv("simple.species.csv", row.names = 1, stringsAsFactors = FALSE)
# species.df.full <- read.csv("Species_Table_021619.csv", row.names = 1, stringsAsFactors = FALSE)
# 
# sample.df <- species.df[c(154, 376, 393, 19, 409),]
# pts.climate.df <- data.frame("species" = character(), "below.maxtemp" = numeric(), 
#                              "above.maxtemp" = numeric(), "below.mintemp" = numeric(), 
#                              "above.mintemp" = numeric(), "maxtemp.loc" = numeric(), 
#                              "mintemp.loc" = numeric(), "nobs" = integer())
# for (i in 1:nrow(sample.df)){
#    species = sample.df[i, 1]
#    
#    # Just getting maxtemp and mintemp where the species is found and comparing it
#    
#    for.maxent.full <- read.csv(paste("./Species_Pts/", gsub(" ", ".", species), "/", 
#                                      gsub(" ", ".", species), ".maxentdatafile.csv", sep = ""), row.names = 1)
#    species.pts.climate <- dplyr::filter(for.maxent.full, for.maxent.full$resp == 1)
#    species.pts.climate <- species.pts.climate[!is.na(species.pts.climate$maxtempwarmestmonth),]
#    below.maxtemp <- sum(species.pts.climate$maxtempwarmestmonth <= upper.crit)
#    above.mintemp <- sum(species.pts.climate$mintempcoldestmonth >= lower.crit)
#    above.maxtemp <- sum(species.pts.climate$maxtempwarmestmonth > upper.crit)
#    below.mintemp <- sum(species.pts.climate$mintempcoldestmonth < lower.crit)
#    
#    
#    maxtemp.loc <- max(species.pts.climate$maxtempwarmestmonth)
#    mintemp.loc <- min(species.pts.climate$mintempcoldestmonth)
#    meantemp.loc <- mean(species.pts.climate$annualmeantemp)
#    nobs <- nrow(species.pts.climate)
#    
#    pts.climate.df <- rbind(data.frame("species" = species, "below.maxtemp" = below.maxtemp, 
#                                       "above.maxtemp" = above.maxtemp, "below.mintemp" = below.mintemp, 
#                                       "above.mintemp" = above.mintemp, "maxtemp.loc" = maxtemp.loc, 
#                                       "mintemp.loc" = mintemp.loc, "meantemp.loc" = meantemp.loc,
#                                       "nobs" = nobs), pts.climate.df)
#    write.csv(pts.climate.df, "./SummaryTables/ptssuitability")
# }
