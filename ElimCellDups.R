elimCellDups_mine <- function (x, r, longLat = NULL, priority = NULL) 
{
   if (class(x) %in% c("SpatialPoints", "SpatialPointsDataFrame")) {
      if (is.na(raster::projection(r))) {
         warning("Raster will be assumed to have same coordinate reference system as points.", 
                 .immediate = TRUE)
         raster::projection(r) <- raster::projection(x)
      }
      else if (raster::projection(r) != raster::projection(x)) {
         stop("Raster and points do not have the same coordinate reference system.")
      }
   }
   xy <- xToCoords(x, longLat, sp = FALSE)
   cellNum <- raster::cellFromXY(r, xy)
   index <- 1:nrow(xy)
   if (is.null(priority)) 
      priority <- 1:nrow(xy)
   removeThese <- integer()
   for (thisCell in unique(cellNum)) {
      if (sum(cellNum %in% thisCell) > 1) {
         thisRow <- index[thisCell == cellNum]
         thisPriority <- priority[thisCell == cellNum]
         thisRow <- thisRow[order(thisPriority)]
         removeThese <- c(removeThese, thisRow[2:length(thisRow)])
      }
   }
   x <- if (class(x) == "data.frame" | class(x) == "matrix") {
      x[-removeThese, ]
   }
   else {
      x[-removeThese]
   }
   x
}