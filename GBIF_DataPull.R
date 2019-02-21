## 2/10/18
## WLDF 510 - Physiology
## Script to get GBIF information for species in the species list and save it in
## /Species_Pts/Species.Name/Species.Name_GBIFPull


## First step: Determine which species have too many observations
## to pull using regular occ_search (limit 200,000)
## 1. Load in CSV
## 2. Get taxa key for all species, put in DF (data frame has species name, key)
## 3. Get observation numbers for all species, put in DF
require(dismo)
require(rgbif)
require(rgeos)
require(raster)
require(rgdal)
require(curl)
species.table.full <- read.csv("Species_Table_021619.csv", row.names = 1, as.is = TRUE,
                               stringsAsFactors = FALSE)

species.table <- species.table.full[,1:2]
species.df <- data.frame("Species" = character(), "taxonKey" = integer(),
                         "nObs" = integer(), stringsAsFactors = FALSE)
for (i in i:nrow(species.table)){
   species <- species.table[i, 2]
   print(paste("Status =", i, "/", nrow(species.table)))
   key.df <- name_suggest(q=species, rank='species')$key
   nobs <- 0
   for (j in 1:length(key.df)){
      this.nobs <- occ_count(taxonKey = key.df[j], georeferenced = TRUE)
      if (this.nobs > nobs){ 
         nobs = this.nobs
         key = key.df[j] 
      } 
   }
   this.row <- data.frame("Species" = species, "taxonKey" = key, "nObs" = nobs)
   print(this.row)
   species.df <- rbind(this.row, species.df, deparse.level = 0)
}

write.csv(simple.species, "simple.species.csv")
require(taxize)
simple.species <- read.csv("simple.species.csv", row.names = 1, stringsAsFactors = FALSE)
simple.species$common.names <- "NA"
for (i in 428:nrow(simple.species)){
   species <- simple.species[i, 1]
   print(paste("Status =", i, "/", nrow(simple.species)))
   common.df <- sci2comm(scinames = species, db = "ncbi")
   if (length(common.df[[1]])==0){
      simple.species$common.names[i] <- NA
   } else {
      simple.species$common.names[i] <- as.character(common.df[[1]])
   }
}

na.species <- na.species[is.na(na.species$common.names),]
these.species <- na.species[!is.na(na.species$common.names),]
some.species <- rbind(na.species, these.species)
some.species$rownames <- rownames(some.species)
rownames(some.species) <- some.species$rownames

test.species <- simple.species
testytest<- cbind(test.species, some.species[, "common.names"][match(test.species$Species, some.species$Species)])
colnames(testytest)[5] <- "na.species.name"
testytest$common.names[is.na(testytest$common.names)] <- as.character(testytest$na.species.name[is.na(testytest$common.names)])
simple.species <- testytest[,1:4]

i = 0
require(rredlist)
for (i in (i+1):nrow(na.species)){
   species <- na.species[i, 1]
   print(paste("Status =", i, "/", nrow(na.species)))
   common.df <- sci2comm(get_tsn(species), db = "itis", simplify = TRUE, verbose = TRUE)
   if (length(common.df[[1]])==0){
      na.species$common.names[i] <- NA
   } else {
      na.species$common.names[i] <- as.character(common.df[[1]][1])
   }
}

na.species$common.names <- simple.species$common.names[1:96]
na.species[rownames(na.species)==442, 4] <- "Uniform swiftlet"

# Pull occ. records and save them for species that have < 100,000 observations
require(dplyr)
require(rgbif)
simple.species <- read.csv("simple.species.csv", row.names = 1, stringsAsFactors = FALSE)
low.obs.species <- filter(simple.species, simple.species$nObs<100000)
gbif.prj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

low.obs.species[359,1]
for (i in 418:nrow(low.obs.species)){
   key <- low.obs.species$taxonKey[i]
   species <- low.obs.species$Species[i]
   species.pts <- occ_search(taxonKey = key, hasCoordinate = TRUE, hasGeospatialIssue = FALSE,
                          eventDate = '1970, 2019', return = "data", 
                          fields = c('decimalLatitude','decimalLongitude','eventDate'
                                     ,'locality', 'geodeticDatum'), limit = 200000, curlopts=list(noprogress=FALSE))
   if (!dir.exists(paste0("./Species_Pts/",gsub(" ", ".", species)))){
      dir.create(paste0("./Species_Pts/",gsub(" ", ".", species)))
   }
   write.csv(species.pts, paste0("./Species_Pts/",gsub(" ", ".", species),"/", gsub(" ", ".", species), ".gbifdata.csv"))
   speciespts.spdf <- SpatialPointsDataFrame(coords = data.frame(species.pts$decimalLongitude, 
                                                                 species.pts$decimalLatitude),
                                             proj4string=gbif.prj, data=species.pts)
   writeOGR(obj = speciespts.spdf, dsn = "./Species_Pts/", layer = paste0(gsub(" ", ".", species), "/",
                                                                          layer = paste0(gsub(" ", ".", species),"sp.pts")), 
            driver = "ESRI Shapefile", overwrite = TRUE)
}


for (i in 1:nrow(simple.species)){
   species <- simple.species[i, 1]
   if (!dir.exists(paste0("./Species_Pts/",gsub(" ", ".", species)))){
      print(species)
   }
}

folders <- list.files("./Species_Pts/")
test <- cbind(sort(folders), sort(simple.species$Species))



