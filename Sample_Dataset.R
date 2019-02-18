# This is a sample code for sample dataset
# Annas Hummingbird - Temperate
# Tufted Coquet - Tropical
# Steps:
# 1. Pull GBIF Data
# 2. Get biogeographic realm shapefile
# 3. Figure out which biogeographic realms points are in; gen random points in those realms
# 4. Load in climate data
# 5. Generate available points within shapefile
# 5. Run maxent
# 6. ???



## Bioclim variables: BIO1=Ann Mean Temp; BIO2=Mean Diurnal Range (Mean of monthly 
# (max temp - min temp)); BIO3=Isothermality (BIO2/BIO7) (* 100); BIO4=Temp Seasonality 
# (SD *100);BIO5 = Max Temp Warmest Mon;BIO6=Min Temp Coldest Mon; BIO7=Temp Ann Range;
# BIO8=Mean Temp Wettest Q;BIO9=Mean Temp Driest Q; BIO10=Mean Temp Warmest Q;
# BIO11=Mean Temp Coldest Q; BIO12 = Annual Precip; BIO13=Precip Wettest Month;
# BIO14=Precip Driest Month; BIO15=Precip Seasonality (CV); BIO16=Precip Wettest Q;
# BIO17=Precip Driest Q; BIO18=Precip Warmest Q; BIO19 = Precip Coldest Q
require(dismo)
require(rgbif)
require(rgdal)
require(rgeos)
require(maptools)
require(sp)
require(raster)

# carpodacus mexicanus
# estrilda melpoda

gbif_user <- "ss665"
gbif_pwd <- "9eKUSuZfMY3r"
gbif_email <- "ss665@humboldt.edu"

species.df.full <- read.csv("Species_Table_021019.csv", row.names = 1, stringsAsFactors = FALSE)
species.df <- read.csv("simple.species.csv", row.names = 1, stringsAsFactors = FALSE)
sum(species.df$nObs>200000)
species.try <- c("Estrilda melpoda", "Colius striatus", "")
species <- "Colius striatus"
key <- species.df[species.df$Species == species, 2]

species.pts <- occ_search(taxonKey = key, hasCoordinate = TRUE, hasGeospatialIssue = FALSE,
                          eventDate = '1970, 2019', return = "data", 
                          fields = c('decimalLatitude','decimalLongitude','eventDate'
                                     ,'locality', 'geodeticDatum'), limit = 200000, curlopts=list(noprogress=FALSE))
dir.create( paste0("./Species_Pts/",gsub(" ", ".", species)))
write.csv(species.pts, paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species), ".gbifdata.csv"), 
          overwrite = TRUE)

gbif.prj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
speciespts.spdf <- SpatialPointsDataFrame(coords = data.frame(species.pts$decimalLongitude, species.pts$decimalLatitude),
                                    proj4string=gbif.prj, data=species.pts)

# ecozone <- readOGR(dsn = "./SpatialFiles/", "Lg_Ecoregions")
ecoregion <- readOGR(dsn = "./SpatialFiles/", "Sm_Ecoregions")
plot(gSimplify(ecoregion, .1))
speciespts.spdf$ecoregion <- as.character(over(speciespts.spdf, ecoregion[,"region"])[,1])
#nas$ecozone <- as.character(over(buffer(nas, width = 10000, dissolve = FALSE), ecozone[,"biome"])[,1])
nas <- speciespts.spdf[as.vector(is.na(speciespts.spdf$ecoregion)),]
nas$ecoregion <- as.character(over(buffer(nas, width = 10000, dissolve = FALSE), ecoregion[,"region"])[,1])

notnas <- speciespts.spdf[as.vector(!is.na(speciespts.spdf$ecoregion)),]


#plot(gSimplify(ecozone, .1), add = TRUE)\
speciespts.spdf <- union(nas, notnas)

regions <- unique(speciespts.spdf$ecoregion)

these.regions <- ecoregion[as.vector(!is.na(match(ecoregion$region,regions))),]
random.pts <- spsample(these.regions, n=nrow(speciespts.spdf), type="random")

speciespts.df <- as.data.frame(speciespts.spdf@data)
write.csv(speciespts.df, paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                               gsub(" ", ".", species), ".editeddata.csv", sep = ""))
simpledata.df <- speciespts.df[, c(7,1:2)]
colnames(simpledata.df) <- c("resp", "x", "y")

random.df <- data.frame("x"=random.pts@coords[,1], "y"=random.pts@coords[,2],"resp"=0)
for.maxent.pts <- rbind(simpledata.df, random.df)

# # Extract climate data
# clim <- stack(file.names)
# ?brick
# clim <- brick(clim)
# writeRaster(clim, "./10m.bioclim.data/clim10m")
# file.names <- list.files("./30s.bioclim.data")
# clim30 <- stack(file.names)
# clim30 <- brick(clim30)
# writeRaster(clim30, "clim30s")
# 
# clim <- raster("clim30s")
clim <- brick("clim10m")
df.clim <- as.data.frame(extract(clim, cbind(for.maxent.pts$x, for.maxent.pts$y)))
print(paste(as.character(species), "climate extraction completed"))
bionames <- c("annualmeantemp", "meandiurnalrange", "isothermality", "tempseasonality",
              "maxtempwarmestmonth", "mintempcoldestmonth", "tempannualrange",
              "meantempwettestquart", "meantempdriestquart", 
              "meantempwarmestquart", "meantempcoldestquart",
              "annualprecip", "precipwettestmo", "precipdriestmo", "precipseasonality",
              "precipwettestquart", "precipdriestquart", "precipwarmestquart", 
              "precipcoldestquart")
colnames(df.clim) <- bionames
for.maxent.full <- cbind(for.maxent.pts, df.clim)

write.csv(for.maxent.full, paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                               gsub(" ", ".", species), ".maxentdatafile.csv", sep = ""))

# run the maxent model. I've decided to use the following predictors:
# annualmeantemp [1], meandiurnalrange [2], maxtempwarmestmonth [5],
# mintempcoldestmonth [6], annualrange [7]
predsnums <- c(1, 2, 5, 6, 7)
preds <- for.maxent.full[,predsnums+3]
resp <- for.maxent.full[,"resp"] 
thisRegMult <- 1
params <- c(paste("betamultiplier=", thisRegMult, sep=""), "jackknife=false")
require(rJava)
maxent.model <- maxent(x = preds, p = resp, path = "./MaxentModels/ScratchDir/",
                       args = params)
curPredsRaster <- clim[[predsnums]]
names(curPredsRaster) <- colnames(preds)
#dir.create(paste0("./MaxentModels/",gsub(" ", ".", species)))
sp.raster <- predict(model = maxent.model, object = curPredsRaster, progress = "text",
        filename = paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species),
                         ".full.tif"), overwrite = TRUE)
thresh <- maxent.model@results[38]
sp.raster.thresh <- sp.raster
sp.raster.thresh[sp.raster.thresh<thresh] = 0
plot(sp.raster)
plot(sp.raster.thresh)
writeRaster(sp.raster.thresh, filename = paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species),
                                                ".thresh.tif"), overwrite = TRUE)

##########################################################
# Trying spatial overlap method

upper.crit <- species.df.full$UCT...C.[species.df.full$Species==species]
lower.crit <-species.df.full$LCT...C.[species.df.full$Species==species]
upper.crit.raster <- clim[[5]] <= upper.crit
lower.crit.raster <- clim[[6]] >= lower.crit
plot(lower.crit.raster)

# Simple function to calculate area from a raster
raster_area <- function(area_raster){
   area_raster[area_raster==0] <- NA
   cell_size <- area(area_raster, na.rm = TRUE, weights = FALSE)
   cell_size <- cell_size[!is.na(cell_size)]
   raster_area <- length(cell_size) * median(cell_size)
   if (length(raster_area) == 0) {return(0)}
   return(raster_area)
}

# Writing a simple function to get a raster of overlapping areas of > 0 suitability.
overlap <- function(raster1, raster2){
   raster1[raster1>0] <- 1
   raster2[raster2>0] <- 1
   add <- raster1+raster2
   add[add<2] <- 0
   return(add)
}

overlap.df <- data.frame("species" =character(), "Area" = double(),
                         "Area.Overlap.Upper" = double(), "Area.Overlap.Lower" = double())
starttime <-Sys.time()
#for (i in 1:nrow(bird.file)){
   # print(paste(i, "/ 77"))
   # cur.bird <- bird.file$Scientific.Name[i]
   # print(paste("Current Bird:", cur.bird))
   upper.crit.raster <- clim10[[5]] <= upper.crit
   lower.crit.raster <- clim10[[6]] >= lower.crit
   data.thresh <- raster(paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species),
                                ".thresh.tif"))
   overlap.upper <- overlap(data.thresh, upper.crit.raster)
   overlap.lower <- overlap(data.thresh, lower.crit.raster)
   this.row <- data.frame("species" = cur.bird, "Area" = raster_area(data.thresh)/1000, 
                          "Area.Overlap.Upper" = raster_area(overlap.upper)/1000, 
                          "Area.Overlap.Lower" = raster_area(overlap.lower)/1000,
                          "Area.Upper" = raster_area(upper.crit.raster)/1000, 
                          "Area.Lower" = raster_area(lower.crit.raster)/1000)
   print(this.row)
   overlap.df <- rbind(this.row, overlap.df)
   #(overlap.df, "areaoverlap.csv")
   cur.runtime <- Sys.time()-starttime
   #print(cur.runtime)
   #print(paste("Time Per Iteration:",round((cur.runtime/i),3)))
   #print(paste("Estimated Time Remaining:", round(((cur.runtime/i)*(nrow(bird.file)-i)),3)))
#}

require(MIAmaxent)
require(maxnet)
# Trying Response curve method


for.maxent.full <- read.csv(paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                  gsub(" ", ".", species), ".maxentdatafile.csv", sep = ""), row.names = 1)
upper.crit <- species.df.full$UCT...C.[species.df.full$Species==species]
lower.crit <-species.df.full$LCT...C.[species.df.full$Species==species]
predsnums <- c(1, 2, 5, 6, 7)
preds <- for.maxent.full[,predsnums+3]
resp <- for.maxent.full[,"resp"] 
thisRegMult <- 1
params <- c(paste("betamultiplier=", thisRegMult, sep=""), "jackknife=false", "writeplotdata=true")
require(rJava)
maxent.model <- maxent(x = preds, p = resp, path = "./MaxentModels/ScratchDir/",
                       args = params)
# maxnet.model <- maxnet( p = resp, data = preds, maxnet.formula(resp, preds, classes="lq"))

savePlot(type = "jpeg", filename = paste0("./Species_Pts/", gsub(" ", ".", species), "/", 
           gsub(" ", ".", species), ".jpeg"))
response(maxent.model, "maxtempwarmestmonth")
abline(v = c(upper.crit, lower.crit))
dev.off()
#########################################
# Trying value extraction method

upper.crit <- species.df.full$UCT...C.[species.df.full$Species==species]
lower.crit <-species.df.full$LCT...C.[species.df.full$Species==species]
require(dplyr)
data.thresh <- raster(paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species),
                             ".full.tif"))
maxent.vals <- getValues(data.thresh)[!is.na(getValues(data.thresh))]
max.temp.vals <- getValues(clim[[5]])[!is.na(getValues(clim[[5]]))]
min.temp.vals <- getValues(clim[[6]])[!is.na(getValues(clim[[6]]))]

vals.df <- data.frame(maxent.vals, max.temp.vals, min.temp.vals)
vals.df <- dplyr::filter(vals.df, max.temp.vals != 0 & min.temp.vals != 0)
# Take mean of maxent values for values that are inside range and outside range.
maxent.below.maxcrit <- mean(dplyr::filter(vals.df, max.temp.vals <= upper.crit)$maxent.vals)
maxent.above.maxcrit <- mean(dplyr::filter(vals.df, max.temp.vals > upper.crit)$maxent.vals)
maxent.below.mincrit <- mean(dplyr::filter(vals.df, min.temp.vals < lower.crit)$maxent.vals)
maxent.above.mincrit <- mean(dplyr::filter(vals.df, min.temp.vals >= lower.crit)$maxent.vals)
mean.maxent <- mean(maxent.vals)
####################################################

# Just getting maxtemp and mintemp where the species is found and comparing it
for.maxent.full <- read.csv(paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                  gsub(" ", ".", species), ".maxentdatafile.csv", sep = ""), row.names = 1)
species.pts.climate <- dplyr::filter(for.maxent.full, for.maxent.full$resp == 1)
inside.maxtemp <- sum(species.pts.climate$maxtempwarmestmonth <= upper.crit)
inside.mintemp <- sum(species.pts.climate$mintempcoldestmonth >= lower.crit)
outside.maxtemp <- sum(species.pts.climate$maxtempwarmestmonth > upper.crit)
outside.mintemp <- sum(species.pts.climate$mintempcoldestmonth < lower.crit)

maxtemp.loc <- max(species.pts.climate$maxtempwarmestmonth)
mintemp.loc <- min(species.pts.climate$mintempcoldestmonth)
nobs <- length(species.pts.climate)



# Characteristics: Diurnal, Non-Migratory, Non-Hibernating, Non-Microhabitat



   
# test enmSDM train maxent function for getting best beta parameter