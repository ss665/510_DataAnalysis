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
species.table.full <- read.csv("Species_Table.csv", row.names = 1, as.is = TRUE,
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

write.csv(species.df, "simple.species.csv")
