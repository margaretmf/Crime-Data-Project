

## SOURCE CODE ##
## AUTHOR: MATT GERBER, UVA SYSTEMS ENGINEERING 

library(ks)
library(rgdal)
library(rgeos)
library(maptools)
library(lubridate)
library(zoo)

#' Gets a crime sample from a Socrata CSV file.
#' @param crime.csv.path Path to Socrata CSV file.
#' @param sample.size Size of sample. If -1 or larger than the number of actual records, all records will be returned.
#' @param start.month Starting month of sample (1-12).
#' @param end.month Ending month of sample (1-12).
#' @return Data frame containing crime sample with columns for \code{x}-location, \code{y}-location, \code{timestamp}, crime \code{hour}, crime \code{day.of.week}, and crime \code{month}.
sample.crime = function(crime.csv.path, sample.size, start.month = 1, end.month = 12)
{
  # read crime data and take sample
  crimes = read.csv(crime.csv.path, header=TRUE)
  
  # reproject lon/lat points to meters
  crimes.locations.lonlat = cbind(crimes$Longitude, crimes$Latitude)
  crimes.locations.meters = project(crimes.locations.lonlat, proj="+init=epsg:26971")
  
  # convert strings to dates
  crimes.dates = as.data.frame(strptime(crimes[,"Date"],"%m/%d/%Y %I:%M:%S %p"))
  
  # reassemble more friendly crimes.sample matrix
  crimes = as.data.frame(cbind(crimes.locations.meters[,1],  # x value
                               crimes.locations.meters[,2],  # y value
                               crimes.dates,                 # crime date/time
                               hour(crimes.dates[[1]]),      # crime hour
                               wday(crimes.dates[[1]]),      # crime day of week
                               month(crimes.dates[[1]])))    # crime month
  
  # set column names for convenience
  names(crimes) = c("x","y","timestamp","hour","day.of.week","month")
  
  # remove NA values
  crimes = crimes[!(is.na(crimes$x) | is.na(crimes$y) | is.na(crimes$timestamp)),]
  
  # filter by month
  crimes = crimes[crimes$month >= start.month & crimes$month <= end.month,]
  
  # filter by size
  if(sample.size == -1 | sample.size > nrow(crimes))
    sample.size = nrow(crimes)
  
  crimes.sample.rows = sample(nrow(crimes), size=sample.size)
  crimes.sample = crimes[crimes.sample.rows,]
  
  return (crimes.sample)
}

#' Reads an ESRI shapefile.
#' @param path Filesystem path (relative or absolute) to shapefile.
#' @param type Type of shapefile:  "poly" or "points"
#' @param sourceProj Source projection of shapefile in proj4 format (e.g., "+init=epsg:3435").
#' @param targetProj Target projection of shapefile in proj4 format (e.g., "+init=epsg:26971").
#' @return A \code{SpatialPolygonDataFrame} projected as desired.
read.shapefile = function(path, type, sourceProj, targetProj)
{
  if(type == "poly")
    shapefile = readShapePoly(path)
  else if(type == "points")
    shapefile = readShapePoints(path)
  else
    message(paste("Unrecognized shapefile type:", type))
  
  # reproject shapefile
  proj4string(shapefile) = sourceProj
  shapefile = spTransform(shapefile, CRS(targetProj))
  
  return (shapefile)
}

#' Gets the minimum distance from all points in one collection to all points in another.
#' @param points.1 An n-by-2 matrix of points.
#' @param points.2 An n-by-2 matrix of points.
#' @return Vector of distances from each point in \code{points.1} to points in \code{points.2}.
get.min.distances = function(points.1, points.2)
{
  return (apply(points.1, 
                1,
                function(point.1)
                {
                  distances = apply(points.2,
                                    1,
                                    function(point.2, point.1)
                                    {
                                      return (get.euc.distance(point.1, point.2))
                                    },
                                    point.1)
                  
                  return (min(distances))
                }))
}

#' Gets Euclidean distance between two points.
#' @param point.1 First point.
#' @param point.2 Second point.
#' @return Euclidean distance.
get.euc.distance = function(point.1, point.2)
{
  return (sqrt(sum((point.1 - point.2)^2)))
}

#' Gets the coordinate range of a shapefile
#' @param shapefile Shapefile to get range for.
#' @return The coordinate range of \code{shapefile} in the following format:  \code{c(xl, xu, yl, yu)}.
get.shapefile.range = function(shapefile)
{
  bbox = bbox(shapefile)
  return (c(range(bbox["x",]), range(bbox["y",])))
}

#' Gets the number of prediction points in north-south and west-east directions.
#' @param prediction.boundary A \code{shapefile} to analyze.
#' @param north.south Whether to analyze in the north-south (TRUE) or west-east (FALSE) direction.
#' @param prediction.resolution.meters Resolution (spacing) of prediction points in meters.
#' @return Number of points that will achieve the given resolution.
get.num.prediction.points = function(prediction.boundary, north.south=TRUE, prediction.resolution.meters)
{
  city.bbox = bbox(prediction.boundary)
  
  range = c()
  if(north.south)
  {
    range = range(city.bbox["y",])
  }
  else
  {
    range = range(city.bbox["x",])
  }
  
  return (trunc((range[2] - range[1]) / prediction.resolution.meters))
}

#' Gets grid of points covering the bounding box of a shapefile.
#' @param shapefile Shapefile to analyze.
#' @param resolution.meters Resolution of grid in meters.
#' @param filter Whether or not to filter out grid points residing outside the given shapefile.
#' @return An n-by-2 data frame of evenly spaced points points.
get.grid.points = function(shapefile, resolution.meters, filter=TRUE)
{
  shapefile.range = get.shapefile.range(shapefile)
  x.values = seq(shapefile.range[1], shapefile.range[2], resolution.meters)
  y.values = seq(shapefile.range[3], shapefile.range[4], resolution.meters)
  grid.points = expand.grid(x.values, y.values)
  names(grid.points) = c("x", "y")
  
  if(filter)
  {
    grid.points = filter.points.to.boundary(grid.points, proj4string(shapefile), shapefile)
  }
  
  return (grid.points)
}

#' Filters out points that are not within a boundary.
#' @param points Matrix (n-by-2) of points in the projection given by \code{pointsProj4String}.
#' @param pointsProj4String Projection string for points.
#' @param boundary Shapefile boundary to use.
#' @return Filtered version of \code{points}.
filter.points.to.boundary = function(points, pointsProj4String, boundary)
{
  sp.points = SpatialPoints(points, proj4string=CRS(pointsProj4String))
  points.in.boundary = !is.na(over(sp.points, boundary)$OBJECTID)
  points = points[points.in.boundary,]
  return (points)
}

#' Plots a surveillance curve for a threat prediction and set of actual crime locations.
#' @param threat.prediction Dataframe with columns \code{x}, \code{y}, and \code{threat}, which give the x-y location of a prediction and its real-valued threat prediction.
#' @param eval.crime.points Dataframe with columns \code{x} and \code{y}, which give the location of actual crimes that occurred during the evaluation period.
#' @param prediction.resolution.meters Resolution of prediction points in meters.
#' @param add Whether or not to add the plot to the current plot.
plot.surveillance.curve = function(threat.prediction, eval.crime.points, prediction.resolution.meters, add=FALSE)
{
  surv.plot.points = get.surveillance.plot.points(threat.prediction, eval.crime.points, prediction.resolution.meters)
  
  if(add)
  {
    lines(surv.plot.points, lty=2, col="blue")
  }
  else
  {
    plot(surv.plot.points, type="l", xlab="% Area Surveilled", ylab = "% Crime Captured", main="Surveillance Plot")
    
    # show random guess
    abline(0,1,lty="dashed")
  }   
  
  # add AUC 
  x = surv.plot.points[,1]
  y = surv.plot.points[,2]
  idx = order(x)
  auc = round(sum(diff(x[idx])*rollmean(y[idx],2)),digits=2)
  
  auc.x.location = 0.6
  auc.y.location = 0.4
  if(add)
  {
    auc.y.location = 0.3
  }
  
  text(auc.x.location, auc.y.location, paste("AUC=",auc,sep=""))
}

#' Gets the plot points for a surveillance curve for a threat prediction and set of actual crime locations.
#' @param threat.prediction Dataframe with columns \code{x}, \code{y}, and \code{threat}, which give the x-y location of a prediction and its real-valued threat prediction.
#' @param eval.crime.points Dataframe with columns \code{x} and \code{y}, which give the location of actual crimes that occurred during the evaluation period.
#' @param prediction.resolution.meters Resolution of prediction points in meters.
#' @return Plot points for surveillance curve.
get.surveillance.plot.points = function (threat.prediction, eval.crime.points, prediction.resolution.meters)
{
  # sort prediction rows by threat level
  threat.prediction = threat.prediction[order(threat.prediction$threat, decreasing=TRUE),]
  
  # augment prediction with columns for IDs, prediction squares, and captured crime counts
  threat.prediction = as.data.frame(cbind(seq(1,nrow(threat.prediction)),
                                          threat.prediction$threat,
                                          threat.prediction$x - prediction.resolution.meters / 2,
                                          threat.prediction$y - prediction.resolution.meters / 2,
                                          threat.prediction$x + prediction.resolution.meters / 2,
                                          threat.prediction$y + prediction.resolution.meters / 2,
                                          rep(0,nrow(threat.prediction))))
  
  names(threat.prediction) = c("id", "threat", "llx","lly","urx","ury", "captured.crimes")
  
  # get the square for each actual crime point
  eval.crime.points.squares = apply(eval.crime.points, 1, square.index.of.crime, threat.prediction)
  
  # get the number of actual crimes that fell into each square
  for(square in eval.crime.points.squares)
  {
    if(!is.na(square))
    {
      threat.prediction$captured.crimes[square] = threat.prediction$captured.crimes[square] + 1
    }
  }  
  
  threat.prediction$captured.crimes = cumsum(threat.prediction$captured.crimes)
  
  plot.x.points = threat.prediction$id / nrow(threat.prediction)
  plot.y.points = threat.prediction$captured.crimes / nrow(eval.crime.points)
  
  return (cbind(plot.x.points, plot.y.points))
}

#' Gets the index of the predicted square that an actual crime location falls into.
#' @param crime.point Location vector with \code{x} and \code{y} columns giving the location of an actual crime.
#' @param threat.prediction Dataframe with columns \code{x}, \code{y}, and \code{threat}, which give the x-y location of a prediction and its real-valued threat prediction.
square.index.of.crime = function(crime.point, threat.prediction)
{
  indices = threat.prediction[!(threat.prediction$llx > crime.point["x"] |
                                  threat.prediction$urx < crime.point["x"] |
                                  threat.prediction$lly > crime.point["y"] |
                                  threat.prediction$ury < crime.point["y"]), "id"]
  
  if(length(indices) == 1)
    return (indices[1])
  else
  {
    message(paste("Crime point fell into", length(indices), "squares"))
    return (NA)
  }
}

#' Runs a two-dimensional (spatial) KDE.
#' @param sample.points Sample points.
#' @param est.points Where to estimate the KDE.
#' @param max.sample.size Maximum size of the data sample.
#' @return Vector of estimates at \code{est.points}.
run.spatial.kde = function(sample.points, est.points, max.sample.size = -1)
{
  if(max.sample.size == -1 | max.sample.size > nrow(sample.points))
  {
    max.sample.size = nrow(sample.points)
  }
  
  sample.points = sample.points[sample(nrow(sample.points), size=max.sample.size),]
  
  # compute optimal KDE bandwidth
  h = Hpi(sample.points, pilot="dscalar")
  
  # run KDE
  est = kde(sample.points, H=h, eval.points=est.points)$estimate
  
  return (est)
}

plot.spatial.kde = function(kde.est, kde.est.points)
{
  image.x.values = sort(unique(kde.est.points[,1]))
  image.y.values = sort(unique(kde.est.points[,2]))
  image(x = image.x.values, 
        y = image.y.values,
        z = matrix(kde.est, ncol=length(image.y.values), byrow=FALSE),
        col = colorRampPalette(rev(brewer.pal(11,'Spectral')))(32),
        xlab="West-to-East (M)", ylab="South-to-North (M)",
        asp=1)
}




