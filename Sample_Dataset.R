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
# 6. 



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

species <- "Estrilda melpoda"
key <- species.df[species.df$Species == species, 2]

species.pts <- occ_search(taxonKey = key, hasCoordinate = TRUE, hasGeospatialIssue = FALSE,
                          eventDate = '1970, 2019', return = "data", 
                          fields = c('decimalLatitude','decimalLongitude','eventDate'
                                     ,'locality', 'geodeticDatum'), limit = 200000, curlopts=list(noprogress=FALSE))
dir.create( "./Species_Pts/Estrilda.melpoda/")
write.csv(species.pts, "./Species_Pts/Estrilda.melpoda/Estrilda.melpoda.gbifdata.csv")

gbif.prj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
speciespts.spdf <- SpatialPointsDataFrame(coords = data.frame(species.pts$decimalLongitude, species.pts$decimalLatitude),
                                    proj4string=gbif.prj, data=species.pts)

# ecozone <- readOGR(dsn = "./SpatialFiles/", "Lg_Ecoregions")
ecoregion <- readOGR(dsn = "./SpatialFiles/", "Sm_Ecoregions")

speciespts.spdf$ecoregion <- as.character(over(speciespts.spdf, ecoregion[,"region"])[,1])
#nas$ecozone <- as.character(over(buffer(nas, width = 10000, dissolve = FALSE), ecozone[,"biome"])[,1])
nas <- speciespts.spdf[as.vector(is.na(speciespts.spdf$ecoregion)),]
nas$ecoregion <- as.character(over(buffer(nas, width = 10000, dissolve = FALSE), ecoregion[,"region"])[,1])

notnas <- speciespts.spdf[as.vector(!is.na(speciespts.spdf$ecoregion)),]


#plot(gSimplify(ecozone, .1), add = TRUE)\
speciespts.spdf <- union(nas, notnas)

regions <- unique(speciespts.spdf$ecoregion)

these.regions <- ecoregion[as.vector(!is.na(match(ecoregion$region,regions))),]
start <- Sys.time()
random.pts <- spsample(these.regions, n=nrow(speciespts.spdf), type="random")

speciespts.df <- as.data.frame(speciespts.spdf@data)
write.csv(speciespts.df, paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                               gsub(" ", ".", species), ".editeddata.csv", sep = ""))
simpledata.df <- speciespts.df[, c(7,1:2)]
colnames(simpledata.df) <- c("resp", "x", "y")

random.df <- data.frame("x"=random.pts@coords[,1], "y"=random.pts@coords[,2],"resp"=0)
for.maxent.pts <- rbind(simpledata.df, random.df)

# Extract climate data
clim <- stack(file.names)
?brick
clim <- brick(clim)
writeRaster(clim, "./10m.bioclim.data/clim10m")

clim <- raster("./10m.bioclim.data/clim10m")
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
