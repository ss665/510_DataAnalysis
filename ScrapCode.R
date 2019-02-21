require(curl)
require(rgbif)
## Old code for later reference
gbif_user <- "ss665"
gbif_pwd <- "9eKUSuZfMY3r"
gbif_email <- "ss665@humboldt.edu"

key <- as.character(paste("taxonKey =",name_suggest(q=species, rank='species')$key[[1]]))
occ.dl <- occ_download(key, "hasCoordinate = TRUE", "hasGeospatialIssue = FALSE",
                       "year >= 1970", "year <= 2019",
                       user = gbif_user, pwd = gbif_pwd, email = gbif_email,
                       curlopts = list())
occ.dl <- occ_download(body = query1, user = gbif_user, pwd = gbif_pwd, email = gbif_email)
get.output <- occ_download_get(key = occ.dl[[1]], "./Species_Pts/", overwrite = TRUE)
df.orig <- occ_download_import(get.output, path = "./Species_Pts")
data.pts <- read.csv("./Species_Pts/Calypte.anna/occurrence.csv")
unzip("./Species_Pts/0041203-181108115102211.zip", exdir = paste0("./Species_Pts/",
                                                                  gsub(" ", ".", species)))
file.remove(paste0("./Species_Pts/", "0041203-181108115102211", ".zip"))
key <- "0041203-181108115102211"
data.orig <- read.table(paste0("./Species_Pts/",
                gsub(" ", ".", species), "/", key, ".csv"), row.names = NULL, header = TRUE,
                fill = TRUE)
data.orig <- read.csv(paste0("./Species_Pts/",
                               gsub(" ", ".", species), "/", key, ".csv"), row.names = NULL, header = TRUE)

query1 <- '{"creator":"ss665",
   "notification_address":["ss665@gmail.com"],
   "format": "SIMPLE_CSV",
   "predicate": {
      "type":"and","predicates":[
         {
            "type":"equals",
            "key":"TAXON_KEY",
            "value":"5229467"
         },
         {
            "type":"equals",
            "key":"HAS_COORDINATE",
            "value":"TRUE"
         },
         {
            "type":"equals",
            "key":"HAS_GEOSPATIAL_ISSUE",
            "value":"FALSE"
         },
         {
            "type":"greaterThanOrEquals",
            "key":"YEAR",
            "value":"1970"
         },


      ]
   }
}'


unzip("./Species_Pts/0038393-181108115102211.zip", exdir = "./Species_Pts/Calypte.anna")
read.delim("./Species_Pts/Calypte.anna/occurrence.txt")
dforig <- gbif(genus = "Calypte", species = "anna", geo = TRUE, removeZeros = TRUE, start = 1,
               end = 30000)
dforig <- occ_download(scientificName = species, hasCoordinate = TRUE,
                       eventDate = '1970, 2017', return = "data", 
                       fields = c('decimalLatitude','decimalLongitude','eventDate'
                                  ,'locality'), limit = 200000, curlopts=list(noprogress=FALSE),
                       download = FALSE)

dforig <- occ_search(scientificName = species, hasCoordinate = TRUE,
                     eventDate = '1970, 2017', return = "data", 
                     fields = c('decimalLatitude','decimalLongitude','eventDate'
                                ,'locality'), limit = 200000, curlopts=list(noprogress=FALSE))

write.csv(dforig, "./Species_Pts/Lophornis.ornatus/gbif.file.csv")
key <- as.character(paste("taxonKey =",name_suggest(q='Lophornis ornatus', rank='species')$key[[1]]))
rdforig <- occ_download(key, "hasCoordinate = TRUE", "hasGeospatialIssue = FALSE",
                        "year >= 1970", "year <= 2017",
                        user = gbif_user, pwd = gbif_pwd, email = gbif_email)
get.output <- occ_download_get(key = dforig[[1]], "./Species_Pts/", overwrite = TRUE)
df.orig <- occ_download_import(as.download("./Species_Pts/0038395-181108115102211.zip"), path = "./Species_Pts")
data.pts <- read.csv("./Species_Pts/Calypte.anna/occurrence.csv")


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


overlap.df <- data.frame("species" =character(), "Area" = double(),
                         "Area.Overlap.Upper" = double(), "Area.Overlap.Lower" = double(), 
                         "Area.Upper" = double(), "Area.Lower" = double())
species.df.full <- read.csv("Species_Table_021619.csv", row.names = 1, stringsAsFactors = FALSE)

for (i in 1:nrow(sample.df)){
   print(paste(i, "/", nrow(sample.df)))
   species <- sample.df[i, 1]
   print(species)
   
   upper.crit <- species.df.full$UCT...C.[species.df.full$Species==species]
   lower.crit <-species.df.full$LCT...C.[species.df.full$Species==species]
   
   upper.crit.raster <- clim[[5]] <= upper.crit
   lower.crit.raster <- clim[[6]] >= lower.crit
   data.thresh <- raster(paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species),
                                ".thresh.tif"))
   overlap.upper <- overlap(data.thresh, upper.crit.raster)
   overlap.lower <- overlap(data.thresh, lower.crit.raster)
   this.row <- data.frame("species" =species, "Area" = raster_area(data.thresh)/1000, 
                          "Area.Overlap.Upper" = raster_area(overlap.upper)/1000, 
                          "Area.Overlap.Lower" = raster_area(overlap.lower)/1000,
                          "Area.Upper" = raster_area(upper.crit.raster)/1000, 
                          "Area.Lower" = raster_area(lower.crit.raster)/1000)
   print(this.row)
   overlap.df <- rbind(this.row, overlap.df)
   write.csv(overlap.df, "./SummaryTables/areaoverlap.csv")
}
raster_area <- function(area_raster){
   area_raster[area_raster==0] <- NA
   cell_size <- area(area_raster, na.rm = TRUE, weights = FALSE)
   cell_size <- cell_size[!is.na(cell_size)]
   raster_area <- length(cell_size) * median(cell_size)
   if (length(raster_area) == 0) {return(0)}
   return(raster_area)
}
info.df <- data.frame("species" = character(), "threshold" = numeric(), "Suit" = numeric(), 
                      "Area" = numeric())

for (i in 1:nrow(sample.df)){
   species <- sample.df[i, 1]
   print(species)
   print(i)
   # Run maxent
   for.maxent.full <- read.csv(paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                     gsub(" ", ".", species), ".maxentdatafile.region.csv", sep = ""), row.names = 1)
   predsnums <- c(1, 2, 5, 6, 7)
   preds <- for.maxent.full[,predsnums+3]
   resp <- for.maxent.full[,"resp"] 
   thisRegMult <- 1
   params <- c(paste("betamultiplier=", thisRegMult, sep=""), "jackknife=false")
   
   maxent.model <- maxent(x = preds, p = resp, path = "./MaxentModels/ScratchDir/",
                          args = params)
   curPredsRaster <- clim[[predsnums]]
   names(curPredsRaster) <- colnames(preds)
   dir.create(paste0("./MaxentModels/",gsub(" ", ".", species)))
   sp.raster <- predict(model = maxent.model, object = curPredsRaster, progress = "text",
                        filename = paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species),
                                          ".full.tif"), overwrite = TRUE)
   thresh <- maxent.model@results[38]
   sp.raster.thresh <- sp.raster
   sp.raster.thresh[sp.raster.thresh<thresh] = 0
   avg.suit <- mean(sp.raster[sp.raster>0])
   thresh.area <- raster_area(sp.raster.thresh)
   plot(sp.raster)
   plot(sp.raster.thresh)
   writeRaster(sp.raster.thresh, filename = paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species),
                                                   ".thresh.tif"), overwrite = TRUE)
   this.row <- data.frame("species" = species, "threshold" = thresh, "Suit" = avg.suit, 
                          "Area" = thresh.area)
   info.df <- rbind(this.row, info.df)
}
write.csv(info.df, "./SummaryTables/maxent.info.csv")

###############################################################################################

require(dplyr)
species.df.full <- read.csv("Species_Table_021619.csv", row.names = 1, stringsAsFactors = FALSE)

summary.df <- data.frame("Species" = character(), "BelowMaxCrit" = numeric(), "AboveMaxCrit" = numeric(),
                         "BelowMinCrit" = numeric(), "AboveMinCrit" = numeric(),
                         "BtwnCrit" = numeric(), "MeanMaxent" = numeric())

for (i in 1:nrow(sample.df)){
   species <- sample.df[i, 1]
   print(species)
   upper.crit <- species.df.full$UCT...C.[species.df.full$Species==species]
   lower.crit <-species.df.full$LCT...C.[species.df.full$Species==species]
   
   full.model <- raster(paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species),
                               ".full.tif"))
   maxent.vals <- getValues(full.model)[!is.na(getValues(full.model))]
   max.temp.vals <- getValues(clim[[5]])[!is.na(getValues(clim[[5]]))]
   min.temp.vals <- getValues(clim[[6]])[!is.na(getValues(clim[[6]]))]
   
   vals.merge <- data.frame(maxent.vals, max.temp.vals, min.temp.vals)
   vals.df <- dplyr::filter(vals.merge, max.temp.vals != 0 & min.temp.vals != 0 &
                               !is.na(max.temp.vals) & !is.na(min.temp.vals) &
                               !is.na(maxent.vals))
   
   # Take mean of maxent values for values that are inside range and outside range.
   maxent.below.maxcrit <- mean(dplyr::filter(vals.df, max.temp.vals <= upper.crit)$maxent.vals)
   maxent.above.maxcrit <- mean(dplyr::filter(vals.df, max.temp.vals > upper.crit)$maxent.vals)
   maxent.below.mincrit <- mean(dplyr::filter(vals.df, min.temp.vals < lower.crit)$maxent.vals)
   maxent.above.mincrit <- mean(dplyr::filter(vals.df, min.temp.vals >= lower.crit)$maxent.vals)
   maxent.btwn.crit <- mean(dplyr::filter(vals.df, min.temp.vals >= lower.crit & max.temp.vals <= upper.crit)$maxent.vals)
   mean.maxent <- mean(maxent.vals)
   
   this.row <- data.frame("Species" = species, "BelowMaxCrit" = maxent.below.maxcrit, "AboveMaxCrit" = maxent.above.maxcrit,
                          "BelowMinCrit" = maxent.below.mincrit, "AboveMinCrit" = maxent.above.mincrit,
                          "BtwnCrit" = maxent.btwn.crit, "MeanMaxent" = mean.maxent)
   print(this.row)
   summary.df <- rbind(this.row, summary.df)
   write.csv(summary.df, "./SummaryTables/valueextraction.csv")
}


raster_area <- function(area_raster){
   area_raster[area_raster==0] <- NA
   cell_size <- area(area_raster, na.rm = TRUE, weights = FALSE)
   cell_size <- cell_size[!is.na(cell_size)]
   raster_area <- length(cell_size) * median(cell_size)
   if (length(raster_area) == 0) {return(0)}
   return(raster_area)
}
info.df <- data.frame("species" = character(), "threshold" = numeric(), "Suit" = numeric(), 
                      "Area" = numeric())

for (i in 1:nrow(sample.df)){
   species <- sample.df[i, 1]
   print(species)
   print(i)
   # Run maxent
   for.maxent.full <- read.csv(paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                     gsub(" ", ".", species), ".maxentdatafile.region.csv", sep = ""), row.names = 1)
   predsnums <- c(1, 2, 5, 6, 7)
   preds <- for.maxent.full[,predsnums+3]
   resp <- for.maxent.full[,"resp"] 
   thisRegMult <- 1
   params <- c(paste("betamultiplier=", thisRegMult, sep=""), "jackknife=false")
   
   maxent.model <- maxent(x = preds, p = resp, path = "./MaxentModels/ScratchDir/",
                          args = params)
   curPredsRaster <- clim[[predsnums]]
   names(curPredsRaster) <- colnames(preds)
   dir.create(paste0("./MaxentModels/",gsub(" ", ".", species)))
   sp.raster <- predict(model = maxent.model, object = curPredsRaster, progress = "text",
                        filename = paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species),
                                          ".full.tif"), overwrite = TRUE)
   thresh <- maxent.model@results[38]
   sp.raster.thresh <- sp.raster
   sp.raster.thresh[sp.raster.thresh<thresh] = 0
   avg.suit <- mean(sp.raster[sp.raster>0])
   thresh.area <- raster_area(sp.raster.thresh)
   plot(sp.raster)
   plot(sp.raster.thresh)
   writeRaster(sp.raster.thresh, filename = paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species),
                                                   ".thresh.tif"), overwrite = TRUE)
   this.row <- data.frame("species" = species, "threshold" = thresh, "Suit" = avg.suit, 
                          "Area" = thresh.area)
   info.df <- rbind(this.row, info.df)
}
write.csv(info.df, "./SummaryTables/maxent.info.csv")

# Compare
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


overlap.df <- data.frame("species" =character(), "Area" = double(),
                         "Area.Overlap.Upper" = double(), "Area.Overlap.Lower" = double(), 
                         "Area.Upper" = double(), "Area.Lower" = double())
species.df.full <- read.csv("Species_Table_021619.csv", row.names = 1, stringsAsFactors = FALSE)

for (i in 1:nrow(sample.df)){
   print(paste(i, "/", nrow(sample.df)))
   species <- sample.df[i, 1]
   print(species)
   
   upper.crit <- species.df.full$UCT...C.[species.df.full$Species==species]
   lower.crit <-species.df.full$LCT...C.[species.df.full$Species==species]
   
   upper.crit.raster <- clim[[5]] <= upper.crit
   lower.crit.raster <- clim[[6]] >= lower.crit
   data.thresh <- raster(paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species),
                                ".thresh.tif"))
   overlap.upper <- overlap(data.thresh, upper.crit.raster)
   overlap.lower <- overlap(data.thresh, lower.crit.raster)
   this.row <- data.frame("species" =species, "Area" = raster_area(data.thresh)/1000, 
                          "Area.Overlap.Upper" = raster_area(overlap.upper)/1000, 
                          "Area.Overlap.Lower" = raster_area(overlap.lower)/1000,
                          "Area.Upper" = raster_area(upper.crit.raster)/1000, 
                          "Area.Lower" = raster_area(lower.crit.raster)/1000)
   print(this.row)
   overlap.df <- rbind(this.row, overlap.df)
   write.csv(overlap.df, "./SummaryTables/areaoverlap.csv")
}

#######################################################
species = species.df[i, 1]
species.df <- read.csv("simple.species.csv", row.names = 1, stringsAsFactors = FALSE)
species.df.full <- read.csv("Species_Table_021619.csv", row.names = 1, stringsAsFactors = FALSE)

for.maxent.full <- read.csv(paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                  gsub(" ", ".", species), ".maxentdatafile.region.csv", sep = ""), row.names = 1)
upper.crit <- species.df.full$UCT...C.[species.df.full$Species==species]
lower.crit <-species.df.full$LCT...C.[species.df.full$Species==species]
predsnums <- c(1, 2, 5, 6, 7)
preds <- for.maxent.full[,predsnums+3]
resp <- for.maxent.full[,"resp"] 
thisRegMult <- 1
params <- c(paste("betamultiplier=", thisRegMult, sep=""), "jackknife=false", "writeplotdata=true")
maxent.model <- maxent(x = preds, p = resp, path = "./MaxentModels/ScratchDir/",
                       args = params)
areas.row <- AreaUnderResponseCurves(x = maxent.model, var = c("annualmeantemp", "maxtempwarmestmonth", "mintempcoldestmonth"), at = median,
                                     expand = 10, data = NULL, fun = predict, upper.crit = upper.crit, lower.crit = lower.crit)




require(dplyr)
species.df.full <- read.csv("Species_Table_021619.csv", row.names = 1, stringsAsFactors = FALSE)

summary.df <- data.frame("Species" = character(), "BelowMaxCrit" = numeric(), "AboveMaxCrit" = numeric(),
                         "BelowMinCrit" = numeric(), "AboveMinCrit" = numeric(),
                         "BtwnCrit" = numeric(), "MeanMaxent" = numeric())

for (i in 1:nrow(sample.df)){
   species <- sample.df[i, 1]
   print(species)
   upper.crit <- species.df.full$UCT...C.[species.df.full$Species==species]
   lower.crit <-species.df.full$LCT...C.[species.df.full$Species==species]
   
   full.model <- raster(paste0("./MaxentModels/",gsub(" ", ".", species),"/", gsub(" ", ".", species),
                               ".full.tif"))
   maxent.vals <- getValues(full.model)[!is.na(getValues(full.model))]
   max.temp.vals <- getValues(clim[[5]])[!is.na(getValues(clim[[5]]))]
   min.temp.vals <- getValues(clim[[6]])[!is.na(getValues(clim[[6]]))]
   
   vals.merge <- data.frame(maxent.vals, max.temp.vals, min.temp.vals)
   vals.df <- dplyr::filter(vals.merge, max.temp.vals != 0 & min.temp.vals != 0 &
                               !is.na(max.temp.vals) & !is.na(min.temp.vals) &
                               !is.na(maxent.vals))
   
   
   
   
   
# Working on the code for the raster read-in
   require(rgdal)
   
   species.df.full <- read.csv("Species_Table_021619.csv", row.names = 1, stringsAsFactors = FALSE)
   species.df <- read.csv("simple.species.csv", row.names = 1, stringsAsFactors = FALSE)
   gbif.prj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
   clim <- brick("clim10m")
   # Sample DF
   sample.df <- species.df[c(154, 376, 393, 19, 409),]
   
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
      species.pts.thinned <- elimCellDups(SpatialPoints(coordinates(species.pts.orig), 
                                                        proj4string = crs(clim)), clim[[5]])
      #species.pts <- over(species.pts.thinned, species.pts.orig)
      species.pts.spdf <- SpatialPointsDataFrame(coords = data.frame(species.pts$dcmlLng, species.pts$dcmlLtt),
                                                 proj4string=gbif.prj, data=species.pts)
      
      # get information about ecoregions that each point lies within, creating a buffer of
      # 10000 to avoid getting NAs because of coarse borders and islands
      
      #sp.buffer <- buffer(species.pts.spdf, width = 10000, dissolve = FALSE, progress= TRUE)
      #print("Buffered, overlapping")
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
      write.csv(for.maxent.full.biome, paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                             gsub(" ", ".", species), ".maxentdatafile.csv", sep = ""))
      
   }
   
   for (i in 1:nrow(sample.df)){
      species <- sample.df[i, 1] #ith species
      print(species)
      print(i)
      # spatial points file that I made earlier
      species.pts.orig <- readOGR(dsn = paste0("./Species_Pts/", gsub(" ", ".", species), "/"), 
                                  layer = paste0(gsub(" ", ".", species),"sp.pts"))
      species.pts.thinned <- elimCellDups(SpatialPoints(coordinates(species.pts.orig), 
                                                        proj4string = crs(clim)), clim[[5]])
      species.pts <- over(species.pts.thinned, species.pts.orig)
      species.pts.spdf <- SpatialPointsDataFrame(coords = data.frame(species.pts$dcmlLng, species.pts$dcmlLtt),
                                                 proj4string=gbif.prj, data=species.pts)
      # get information about ecoregions that each point lies within, creating a buffer of
      # 10000 to avoid getting NAs because of coarse borders and islands
      
      sp.buffer <- buffer(species.pts.spdf, width = 10000, dissolve = FALSE, progress= TRUE)
      print("Buffered, overlapping")
      # sp.over.df <- as.data.frame(extract(y = SpatialPoints(coordinates(species.pts.spdf), proj4string=gbif.prj), 
      #                                     buffer = 200000, x = ecoregion.raster, small = TRUE, df = TRUE,
      #                                     progress = TRUE))
      # 
      # sp.over.sort <- sp.over.df[order(sp.over.df$layer),]
      # sp.over <- sp.over.sort[!duplicated(sp.over.sort$ID),]
      # sp.over <- merge(sp.over, ecoregions.df[,2:3], all.x = TRUE)
      sp.over <- over(sp.buffer, ecoregion[,c("region", "biome")])
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
      write.csv(for.maxent.full.biome, paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                             gsub(" ", ".", species), ".maxentdatafile.csv", sep = ""))
      
   }
   
   
   
   
   species.df <- read.csv("simple.species.csv", row.names = 1, stringsAsFactors = FALSE)
   species.df.full <- read.csv("Species_Table_021619.csv", row.names = 1, stringsAsFactors = FALSE)
   df <- data.frame("species" = character(), "full.area" = numeric(), "above.lower.crit" = numeric(), "below.lower.crit" = numeric(),
                    "above.upper.crit" = numeric(), "below.upper.crit" = numeric(), 
                    "area.between.crits" = numeric())
   sample.df <- species.df[c(154, 376, 393, 19, 409, 361, 131, 26, 281, 133),]
   
   for (i in 1:nrow(sample.df)){
      species = sample.df[i, 1]
      print(species)
      for.maxent.full <- read.csv(paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                        gsub(" ", ".", species), ".maxentdatafile.region.csv", sep = ""), row.names = 1)
      upper.crit <- species.df.full$UCT...C.[species.df.full$Species==species]
      lower.crit <-species.df.full$LCT...C.[species.df.full$Species==species]
      predsnums <- c(1, 2, 5, 6, 7)
      preds <- for.maxent.full[,predsnums+3]
      resp <- for.maxent.full[,"resp"] 
      thisRegMult <- c(0.5, 1, 2, 5, 10)
      params <- c(paste("betamultiplier=", thisRegMult, sep=""), "jackknife=false", "writeplotdata=true")
      maxent.model <- maxent(x = preds, p = resp, path = "./MaxentModels/ScratchDir/",
                             args = params)
      areas.row <- AreaUnderResponseCurves(x = maxent.model, var = c("annualmeantemp", "maxtempwarmestmonth", "mintempcoldestmonth"), at = median,
                                           expand = 10, data = NULL, fun = predict, upper.crit = upper.crit, lower.crit = lower.crit)
      this.row <- cbind("species" = species, areas.row)
      df <- rbind(this.row, df)
      write.csv(df, "./SummaryTables/response.suitability.csv")
   }
   
   
   species.df <- read.csv("simple.species.csv", row.names = 1, stringsAsFactors = FALSE)
   species.df.full <- read.csv("Species_Table_021619.csv", row.names = 1, stringsAsFactors = FALSE)
   df <- data.frame("species" = character(), "full.area" = numeric(), "above.lower.crit" = numeric(), "below.lower.crit" = numeric(),
                    "above.upper.crit" = numeric(), "below.upper.crit" = numeric(), 
                    "area.between.crits" = numeric())
   sample.df <- species.df[c(154, 376, 393, 19, 409, 361, 131, 26, 281, 133),]
   
   for (i in 1:nrow(sample.df)){
      species = sample.df[i, 1]
      print(species)
      for.maxent.full <- read.csv(paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                        gsub(" ", ".", species), ".maxentdatafile.region.csv", sep = ""), row.names = 1)
      upper.crit <- species.df.full$UCT...C.[species.df.full$Species==species]
      lower.crit <-species.df.full$LCT...C.[species.df.full$Species==species]
      predsnums <- c(1, 2, 5, 6, 7)
      preds <- for.maxent.full[,predsnums+3]
      resp <- for.maxent.full[,"resp"] 
      thisRegMult <- (1)
      #params <- c(paste("betamultiplier=", thisRegMult, sep=""), "jackknife=false", "writeplotdata=true")
      maxent.model <- trainMaxEnt(data = for.maxent.full, scratchDir = "./MaxentModels/ScratchDir/",
                                  regMult = thisRegMult, jacknife = FALSE, verbose = TRUE, out = c("model", "tuning"),
                                  testClasses = FALSE, silent = FALSE, classes = "lpqht")
      areas.row <- AreaUnderResponseCurves(x = maxent.model[[1]], var = c("annualmeantemp", "maxtempwarmestmonth", "mintempcoldestmonth"), at = median,
                                           expand = 10, data = NULL, fun = predict, upper.crit = upper.crit, lower.crit = lower.crit)
      this.row <- cbind("species" = species, areas.row)
      df <- rbind(this.row, df)
      write.csv(df, "./SummaryTables/response.suitability.region.csv")
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
   
   
   
   
   
   #cur.runtime <- Sys.time()-starttime
   
   
   

