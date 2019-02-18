{
   .local <- function (x, var = NULL, at = median, range = "pa", 
                       expand = 10, rug = FALSE, data = NULL, fun = predict, 
                       ylim = c(0, 1), col = "red", lwd = 2, add = FALSE, ...) 
   {
      stopifnot(range %in% c("p", "pa")) #range is presence or presence absence
      if (is.null(d)) { # Combines maxent feedout to presences and absences
         d <- x@presence
         if (range == "pa" & x@hasabsence) {
            data <- rbind(d, x@absence)
         }
      }
      cn <- colnames(d) # cn is colnames of data
      # if (is.null(var)) { # if var not provided, makes var all of column names
      #    var <- cn
      # }
      # if (is.numeric(var)) {
      #    var <- cn[var]
      # }
      # if (length(var) == 1) {
      # }
      # var <- var[var %in% cn]
      # if (length(var) == 0) {
      #    stop("var not found")
      # }
      .doResponse(x, var, at, data, cn, expand, rug, ylim, 
                  col, lwd, add, fun, ...)
   }
   .local(x, ...)
}

function (x, var, at, d, cn, expand, rug, ylim, col, lwd, add, 
          fun, ...) 
{
   # if (length(var) > 1 & !add) { # setting parameters if var is > 1
   #    old.par <- graphics::par(no.readonly = TRUE)
   #    on.exit(graphics::par(old.par))
   #    xs <- floor(sqrt(length(var)))
   #    ys <- ceiling(length(var)/xs)
   #    graphics::par(mfrow = c(xs, ys))
   # }
   f <- sapply(d, is.factor) # d is data - presence ans absence values. checks if any are factors
   notf <- !f # data that is not factors
   m <- matrix(nrow = 1, ncol = ncol(d)) # matrix with number of columns for data columns
   # if (is.null(at)) { # at is a function to indicate what level other variables should be fixed at
   #    m <- d # putting data in the matrix 
   #    nrm <- nrow(m)
   # }
   # else 
   if (is.function(at)) { # Dealing with functions in case you put one in for at
      if (sum(notf) > 0) {
         m[notf] <- as.numeric(apply(d[, notf, drop = FALSE],
                                     2, at))
      }
      # if (sum(f) > 0) {
      #    m[f] <- as.numeric(apply(d[, f, drop = FALSE], 2,
      #                             modal))
      # }
      m <- matrix(m, nrow = 100, ncol = length(m), byrow = TRUE)
      colnames(m) <- cn
   }
   # else {
   #    at <- at[cn]
   #    m <- as.vector(at)
   #    m <- matrix(m, nrow = 100, ncol = length(m), byrow = TRUE)
   #    colnames(m) <- names(at)
   # }
   m <- data.frame(m)
   for (vr in var) {
      i <- which(cn == vr)
      # if (is.null(at)) {
      #    nr <- ifelse(length(var) == 1, 25, 10)
      #    v <- d[, i]
      #    if (is.factor(v)) {
      #       v <- as.numeric(levels(v))
      #       fact <- TRUE
      #    }
      #    else {
      #       fact <- FALSE
      #       r <- range(v)
      #       expand <- round(abs(expand))
      #       v <- (r[1] - expand) + 0:(nr - 1) * (r[2] - r[1] + 
      #                                               2 * expand)/(nr - 1)
      #    }
      #    mm <- m[rep(1:nrm, length(v)), ]
      #    mm[, vr] <- rep(v, each = nrm)
      #    p <- fun(x, mm)
      #    pd <- cbind(v, colMeans(matrix(p, nrow = nrm), na.rm = TRUE))
      # }   else {
         nr <- 100
         v <- d[, i]
         # if (is.factor(v)) {
         #    v <- as.numeric(levels(v))
         #    v <- rep(v, ceiling(100/length(v)))
         #    v <- v[1:100]
         #    fact <- TRUE
         # }
         # else {
            fact <- FALSE
            r <- range(v)
            expand <- round(abs(expand))
            v <- (r[1] - expand) + 0:(nr - 1) * (r[2] - r[1] + 
                                                    2 * expand)/(nr - 1)
         # }
         mm <- m
         mm[, vr] <- v
         p <- fun(x, mm)
         pd <- cbind(mm[, vr], p)
      #}
      # if (add) {
      #    if (fact) {
      #       points(pd, col = col, lwd = lwd, ...)
      #    }
      #    else {
      #       points(pd, col = col, lwd = lwd, type = "l", 
      #              ...)
      #    }
      # }
      #else {
         # if (fact) {
         #    plot(pd, xlab = vr, ylab = "predicted value",
         #         col = col, lwd = lwd, ylim = ylim, ...)
         # }
         # else {
            plot(pd, xlab = vr, ylab = "predicted value",
                 col = col, lwd = lwd, ylim = ylim, type = "l")
            AUC(x = pd[,1]>=lower.crit, y = pd[,2])
         #}
         # if (rug) {
         #    if (!is.factor(d[, i])) {
         #       rug(quantile(d[, i], probs = seq(0, 1, 0.1)),
         #           col = "blue")
         #    }
         # }
      }
   }
   if (length(var) == 1) {
      return(invisible(pd))
   }
}