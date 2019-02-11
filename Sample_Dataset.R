# This is a sample code for sample dataset
# Annas Hummingbird - Temperate
# Tufted Coquet - Tropical
# Steps:
# 1. Pull GBIF Data
# 2. Get biogeographic realm shapefile
# 3. Load in climate data and clip it by shapefile
# 4. Generate available points within shapefile
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


gbif_user <- "ss665"
gbif_pwd <- "9eKUSuZfMY3r"
gbif_email <- "ss665@humboldt.edu"
species = "Calypte anna"
species = "Lophornis ornatus"
# Lophornis ornatus

key <- as.character(paste("taxonKey =",name_suggest(q='Calypte anna', rank='species')$key[[1]]))
dforig <- occ_download(key, "hasCoordinate = TRUE", "hasGeospatialIssue = FALSE",
                       "year >= 1970", "year <= 2017", "fields = c('decimalLatitude', 'decimalLongitude', 'eventDate', 'continent')",
                     user = gbif_user, pwd = gbif_pwd, email = gbif_email)
get.output <- occ_download_get(key = dforig[[1]], "./Species_Pts/", overwrite = TRUE)
df.orig <- occ_download_import(as.download("./Species_Pts/0038395-181108115102211.zip"), path = "./Species_Pts")
data.pts <- read.csv("./Species_Pts/Calypte.anna/occurrence.csv")

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

