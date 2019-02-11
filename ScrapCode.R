## Old code for later reference


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

