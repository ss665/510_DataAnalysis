# Code to run maxent model for all species (loading in the species spatial dataframe) made
# in previous code, then also calculate area under the curve (run AreaUnderResponseCurves 
# function first), and pull out temperature information for species. Creates 
# 3 summary csvs - output of best model (for beta selection), mean temperature values for sp. points,
# and area under curves.

species.df <- read.csv("simple.species.csv", row.names = 1, stringsAsFactors = FALSE)
species.df.full <- read.csv("Species_Table_021619.csv", row.names = 1, stringsAsFactors = FALSE)
df <- data.frame("species" = character(), "full.area" = numeric(), "above.lower.crit" = numeric(), "below.lower.crit" = numeric(),
                 "above.upper.crit" = numeric(), "below.upper.crit" = numeric(), 
                 "area.between.crits" = numeric())
pts.climate.df <- data.frame("species" = character(), "below.maxtemp" = numeric(), 
                             "above.maxtemp" = numeric(), "below.mintemp" = numeric(), 
                             "above.mintemp" = numeric(), "maxtemp.loc" = numeric(), 
                             "mintemp.loc" = numeric(), "maxcrit" = numeric(), "mincrit" = numeric(),
                             "nobs" = integer())
best.reg.mults <- data.frame("regMult" = numeric(), "logLik" = numeric(),
                             "K" = numeric(), "AICc" = numeric(), 
                             "deltaAICc" = numeric(), "aicWeight" = numeric())
sample.df <- species.df[c(154, 376, 393, 19, 409, 361, 131, 26, 281, 133),]
# type in biome or region here to get appropriate points
bor = "biome"
for (i in 1:nrow(sample.df)){
   species = sample.df[i, 1]
   print(species)
   for.maxent.full <- read.csv(paste("./Species_Pts/", gsub(" ", ".", species), "/", 
                                     gsub(" ", ".", species), ".maxentdatafile.",bor,".csv", sep = ""), row.names = 1)
   upper.crit <- species.df.full$UCT...C.[species.df.full$Species==species]
   lower.crit <-species.df.full$LCT...C.[species.df.full$Species==species]
   predsnums <- c(1, 2, 5, 6, 7)
   preds <- for.maxent.full[,predsnums+3]
   resp <- for.maxent.full[,"resp"] 
   maxent.data <- cbind(resp, preds)
   thisRegMult <- c(0.5, 1, 2, 5)
   # params <- c(paste("betamultiplier=", thisRegMult, sep=""), "jackknife=false")
   dir.create(paste0("./MaxentModels/",gsub(" ", ".", species)))
   maxent.model <- trainMaxEnt(data = maxent.data , 
                               regMult = thisRegMult, testClasses = FALSE, out = c("model", "tuning"),
                               scratchDir = paste0("./MaxentModels/",gsub(" ", ".", species)),
                               jackknife = FALSE, verbose = TRUE, anyway = FALSE, classes = "lpqht")
   areas.row <- AreaUnderResponseCurves(x = maxent.model[[2]], var = c("annualmeantemp", "maxtempwarmestmonth", "mintempcoldestmonth"), at = median,
                                        expand = 10, data = NULL, fun = predict, upper.crit = upper.crit, lower.crit = lower.crit)
   tuning.df <- maxent.model[[1]]
   best.reg.mults <- rbind(best.reg.mults, data.frame("regMult" = tuning.df[1,1], "logLik" = tuning.df[1, 9],
                                                      "K" = tuning.df[1, 10], "AICc" =  tuning.df[1, 11], 
                                                      "deltaAICc" = tuning.df[2, 12], "aicWeight" = tuning.df[1, 14]))
   write.csv(best.reg.mults, paste0("./SummaryTables/reg.mults.", bor, ".csv"))
   df <- rbind(cbind("species" = species, areas.row), df)
   write.csv(df, paste0("./SummaryTables/response.suitability.",bor,".csv"))
   
   
   #Evaluating climate at occurence points
   species.pts.climate <- dplyr::filter(for.maxent.full, for.maxent.full$resp == 1)
   species.pts.climate <- species.pts.climate[!is.na(species.pts.climate$maxtempwarmestmonth),]
   below.maxtemp <- sum(species.pts.climate$maxtempwarmestmonth <= upper.crit)
   above.mintemp <- sum(species.pts.climate$mintempcoldestmonth >= lower.crit)
   above.maxtemp <- sum(species.pts.climate$maxtempwarmestmonth > upper.crit)
   below.mintemp <- sum(species.pts.climate$mintempcoldestmonth < lower.crit)
   
   
   maxtemp.loc <- max(species.pts.climate$maxtempwarmestmonth)
   mintemp.loc <- min(species.pts.climate$mintempcoldestmonth)
   meantemp.loc <- mean(species.pts.climate$annualmeantemp)
   nobs <- nrow(species.pts.climate)
   
   pts.climate.df <- rbind(data.frame("species" = species, "below.maxtemp" = below.maxtemp, 
                                      "above.maxtemp" = above.maxtemp, "below.mintemp" = below.mintemp, 
                                      "above.mintemp" = above.mintemp, "maxtemp.loc" = maxtemp.loc, 
                                      "mintemp.loc" = mintemp.loc, "meantemp.loc" = meantemp.loc,
                                      "nobs" = nobs, "maxcrit" = upper.crit, "mincrit" = lower.crit), pts.climate.df)
   write.csv(pts.climate.df, paste0("./SummaryTables/ptssuitability.",bor,".csv"))
}