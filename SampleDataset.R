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
sample.df <- species.df[c(154, 376, 393, 19, 409),]

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
######################################################################
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
require(dplyr)

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
