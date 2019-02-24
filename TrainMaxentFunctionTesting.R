trainMaxEnt(preds = preds, resp = resp, 
            regMult = thisRegMult, testClasses = FALSE, out = c("model, tuning"),
            scratchDir = paste0("./MaxentModels/",gsub(" ", ".", species)),
            jackknife = FALSE, verbose = TRUE)


function (data, resp = names(data)[1], preds = names(data)[2:ncol(data)], 
          regMult = c(seq(0.5, 5, by = 0.5), 6:10, 12.5, 15, 17.5, 
                      20), classes = "default", testClasses = TRUE, forceLinear = TRUE, 
          jackknife = TRUE, args = "", out = "model", anyway = TRUE, 
          scratchDir = NULL, verbose = FALSE, ...) {
   if (class(resp) %in% c("integer", "numeric")){ 
      resp <- names(data)[resp]}
   if (class(preds) %in% c("integer", "numeric")) {
      preds <- names(data)[preds]}
   presentBg <- data[, resp]
   data <- data[, preds, drop = FALSE]
   if (classes == "default") {
      classesToTest <- c("l", "p", "q", "h")
   } else {
      classesToTest <- rep(NA, nchar(classes))
      for (i in 1:nchar(classes)) classesToTest[i] <- substr(classes, 
                                                             i, i)
   }
   if (any("p" %in% classesToTest) & ncol(data) == 1) {
      product <- FALSE
      warning("Data has only one variable so forcing product features to FALSE.")
   }
   scratchDir <- if (is.null(scratchDir)) {
      base::tempfile(pattern = "/_maxentTempFiles/")
   }
   else {
      base::tempfile(pattern = "/_maxentTempFiles/", tmpdir = scratchDir)
   }
   dirCreate(scratchDir)
   dirCreate(scratchDir, "/plots")
   allPres <- data[presentBg == 1, , drop = FALSE]
   allBg <- data[presentBg == 0, , drop = FALSE]
   if (testClasses) {
      classGrid <- expand.grid(rep(list(c(1, 0)), length(classesToTest)))
      classGrid <- classGrid[-which(rowSums(classGrid) == 0), 
                             ]
   } else {
      classGrid <- data.frame(matrix(rep(1, length(classesToTest)), 
                                     nrow = 1))
   }
   names(classGrid) <- classesToTest
   if (forceLinear & any(classGrid$l == 0)){ 
      classGrid <- classGrid[-which(classGrid$l == 0), ]}
   if (verbose) {
      say("Testing models with regularization multiplier:", 
          post = 0)}
   for (thisRegMult in regMult) {
      if (verbose) {
         omnibus::say(thisRegMult, post = 0)}
      for (countParam in 1:nrow(classGrid)) {
         params <- c(paste0("betamultiplier=", thisRegMult), 
                     paste0("linear=", ifelse("l" %in% names(classGrid) && 
                                                 classGrid$l[countParam] == 1, "true", "false")), 
                     paste0("product=", ifelse("p" %in% names(classGrid) && 
                                                  classGrid$p[countParam] == 1, "true", "false")), 
                     paste0("quadratic=", ifelse("q" %in% names(classGrid) && 
                                                    classGrid$q[countParam] == 1, "true", "false")), 
                     paste0("hinge=", ifelse("h" %in% names(classGrid) && 
                                                classGrid$h[countParam] == 1, "true", "false")), 
                     paste0("threshold=", ifelse("t" %in% names(classGrid) && 
                                                    classGrid$t[countParam] == 1, "true", "false")), 
                     "jackknife=false")
         if (args != "") {
            params <- c(params, args)}
         trialModel <- dismo::maxent(x = data, p = as.vector(presentBg), 
                                     path = scratchDir, args = params)
         predPres <- dismo::predict(object = trialModel, x = allPres, 
                                    na.rm = TRUE, args = "outputformat=raw")
         predBg <- dismo::predict(object = trialModel, x = allBg, 
                                  na.rm = TRUE, args = "outputformat=raw")
         
         bgSum <- sum(predBg)
         ll <- sum(log(predPres/bgSum), na.rm = TRUE)
         K <- 0
         for (thisLambda in trialModel@lambdas) {
            if (!grepl(thisLambda, pattern = "linearPredictorNormalizer") & 
                !grepl(thisLambda, pattern = "densityNormalizer") & 
                !grepl(thisLambda, pattern = "numBackgroundPoints") & 
                !grepl(thisLambda, pattern = "entropy")) {
               split <- strsplit(thisLambda, ", ")
               paramValue <- as.numeric(split[[1]][2])
               if (paramValue != 0) 
                  K <- K + 1
            }
         }
         AICc <- -2 * ll + 2 * K + (2 * K * (K + 1))/(sum(presentBg) - 
                                                         K - 1)
         thisAicFrame <- data.frame(regMult = thisRegMult, 
                                    linear = as.logical(ifelse("l" %in% names(classGrid), 
                                                               classGrid$l[countParam], FALSE)), quadratic = as.logical(ifelse("q" %in% 
                                                                                                                                  names(classGrid), classGrid$q[countParam], 
                                                                                                                               FALSE)), product = as.logical(ifelse("p" %in% 
                                                                                                                                                                       names(classGrid), classGrid$p[countParam], 
                                                                                                                                                                    FALSE)), hinge = as.logical(ifelse("h" %in% 
                                                                                                                                                                                                          names(classGrid), classGrid$h[countParam], 
                                                                                                                                                                                                       FALSE)), threshold = as.logical(ifelse("t" %in% 
                                                                                                                                                                                                                                                 names(classGrid), classGrid$t[countParam], 
                                                                                                                                                                                                                                              FALSE)), numFeats = sum(classGrid[countParam, 
                                                                                                                                                                                                                                                                                ]), n = sum(presentBg), logLik = ll, K = K, 
                                    AICc = AICc)
         tuning <- if (exists("tuning", inherits = FALSE)) {
            rbind(tuning, thisAicFrame)
         }
         else {
            thisAicFrame
         }
      }
   }
   if (!anyway) 
      tuning <- tuning[which(tuning$n >= tuning$K), ]
   if (nrow(tuning) > 0) {
      tuning <- tuning[order(tuning$K, decreasing = FALSE), 
                       ]
      tuning <- tuning[order(tuning$linear, tuning$quadratic, 
                             tuning$threshold, tuning$hinge, tuning$product, decreasing = TRUE), 
                       ]
      tuning <- tuning[order(tuning$numFeats, decreasing = FALSE), 
                       ]
      tuning <- tuning[order(tuning$regMult, decreasing = TRUE), 
                       ]
      tuning <- tuning[order(tuning$AICc, decreasing = FALSE), 
                       ]
      tuning$deltaAICc <- tuning$AICc - min(tuning$AICc)
      tuning$relLike <- exp(-0.5 * tuning$deltaAICc)
      tuning$aicWeight <- tuning$relLike/sum(tuning$relLike)
   }
   if (verbose) {
      omnibus::say("")
      print(tuning)
      omnibus::say("")
   }
   if ("model" %in% out) {
      if (nrow(tuning) == 0 & !anyway) {
         warning("No models had fewer coefficients than predictors. No model returned.", 
                 immediate. = TRUE)
         model <- "No MAXENT model had number of parameters < number of training presences."
      }
      else {
         if (nrow(tuning) > 0) {
            params <- c(paste0("betamultiplier=", tuning$regMult[1]), 
                        paste0("linear=", tolower(tuning$linear[1])), 
                        paste0("quadratic=", tolower(tuning$quadratic[1])), 
                        paste0("product=", tolower(tuning$product[1])), 
                        paste0("hinge=", tolower(tuning$hinge[1])), 
                        paste0("threshold=", tolower(tuning$threshold[1])), 
                        paste0("jackknife=", tolower(jackknife)))
         }
         else if (anyway) {
            warning("Returning model with multipler = 1 even though no model had fewer coefficients than predictors.", 
                    immediate. = TRUE)
            params <- c(paste0("betamultiplier=1"), paste0("linear=", 
                                                           tolower(grepl(pattern = "l", classes) | classses == 
                                                                      "default")), paste0("quadratic=", tolower(grepl(pattern = "q", 
                                                                                                                      classes) | classses == "default")), paste0("product=", 
                                                                                                                                                                 tolower(grepl(pattern = "p", classes) | classses == 
                                                                                                                                                                            "default")), paste0("hinge=", tolower(grepl(pattern = "h", 
                                                                                                                                                                                                                        classes) | classses == "default")), paste0("threshold=", 
                                                                                                                                                                                                                                                                   tolower(grepl(pattern = "t", classes))), paste0("jackknife=", 
                                                                                                                                                                                                                                                                                                                   tolower(jackknife)))
         }
         if (args != "") 
            params <- c(params, args)
         model <- dismo::maxent(x = data, p = as.vector(presentBg), 
                                removeDuplicates = FALSE, path = scratchDir, 
                                args = params)
      }
   }
   write.csv(NULL, paste0(scratchDir, "/species.lambdas"))
   if (file.exists(paste0(scratchDir, "/presences"))) 
      write.csv(NULL, paste0(scratchDir, "/presences"))
   if (file.exists(paste0(scratchDir, "/absences"))) 
      write.csv(NULL, paste0(scratchDir, "/absences"))
   unlink(paste0(scratchDir, "/plots"), recursive = TRUE, force = TRUE)
   unlink(scratchDir, recursive = TRUE, force = TRUE)
   if ("model" %in% out & !("tuning" %in% out)) {
      model
   }
   else if (!("model" %in% out) & "tuning" %in% out) {
      tuning
   }
   else {
      list(tuning = tuning, model = model)
   }
}
