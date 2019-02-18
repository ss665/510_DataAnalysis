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

