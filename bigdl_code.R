require(rgbif)
require(rgdal)
require(dplyr)
gbif_user <- "ss665"
gbif_pwd <- "9eKUSuZfMY3r"
gbif_email <- "ss665@humboldt.edu"

simple.species <- read.csv("simple.species.csv", row.names = 1, stringsAsFactors = FALSE)
high.obs.species <- filter(simple.species, simple.species$nObs>100000) 
gbif.prj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#species <- "Peromyscus maniculatus"
#key <- as.character(paste("taxonKey =",name_suggest(q=species, rank='species')$key[[1]]))
# Using function
for (i in 24:nrow(high.obs.species)){
   key <- as.character(paste("taxonKey =",high.obs.species$taxonKey[i]))
   species <- high.obs.species$Species[i]
   nobs <- high.obs.species$nObs[i]
   print(species)
   print(nobs)
   print(paste(i, "/",nrow(high.obs.species)))
   occ.dl <- occ_download(key, "hasCoordinate = TRUE", "hasGeospatialIssue = FALSE",
                                   "year >= 1970", "year <= 2019",
                                   user = gbif_user, pwd = gbif_pwd, email = gbif_email)
   meta <- occ_download_meta(occ.dl)
   print(meta)
   while(meta[[8]]!="SUCCEEDED"){
      Sys.sleep(120)
      meta <- occ_download_meta(occ.dl)
      print(meta)
   }
   print(meta[[8]])
   get.output <- occ_download_get(key = occ.dl[[1]], "./Species_Pts/", overwrite = TRUE)
   species.pts <- occ_download_import(get.output, path = "./Species_Pts/", fill=FALSE, select = c('decimalLatitude','decimalLongitude','eventDate'
                                                                              ,'locality')) 
   file.remove(paste0("./Species_Pts/", occ.dl[[1]], ".zip"))
   species.pts$geodeticDatum = "NA"
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







# Using JSON
key <- as.character(paste("taxonKey =",name_suggest(q=species, rank='species')$key[[1]]))
query1 <- '{"creator":"ss665",
   "notification_address":["ss665@gmail.com"],
"format": "SIMPLE_CSV",
"predicate": {
"type":"and","predicates":[
{
   "type":"equals",
   "key":"TAXON_KEY",
   "value":"2436248"
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

occ.dl <- occ_download(body = query1, user = gbif_user, pwd = gbif_pwd, email = gbif_email)
get.output <- occ_download_get(key = occ.dl[[1]], ".", overwrite = TRUE)
df.orig <- occ_download_import(get.output, path = ".") # This causes my r-session to terminate/abort


# Trying to unzip the CSV and import it (from the JSON query)
unzip("PUT YOUR FILE NAME HERE")
data.orig <- data.orig <- fread(paste0("./Species_Pts/",
                                            gsub(" ", ".", species), "/", key, ".csv"))