#########################################

## EXPLORATORY CRIME DATA ANALYSIS ##
## AUTHOR: MARGARET FURR ##

#########################################

library(ks); library(rgdal)
library(rgeos); library(RColorBrewer)
library(ggplot2); library(ggmap)
library(zoo); attach(mtcars)
library(car); library(leaps)

# source code
source("CrimeUtil.R")

# KDE
kde.resolution.meters = 200

# read chicago boundary
city.boundary = read.shapefile("City_Boundary/City_Boundary", "poly", "+init=epsg:3435", "+init=epsg:26971")

#########################################

# 2014 Theft Crime By Month

# 2014 theft data
theft = sample.crime("2014_THEFT.csv", -1, 2, 12)
theft.prior.month <- c()
theft.post.month <- c()
for (i in 2:10) {
  theft.prior.month[[i]] <- theft[theft$month==i,]
}
for (i in 3:11) {
  theft.post.month[[i]] <- theft[theft$month==i,]
}

# get KDE sample data from months 2-10
theft.kde.sample.points.month <- c()
for (i in 2:10) {
  theft.kde.sample.points.month[[i]] = theft.prior.month[[i]][,c("x","y")]
}

# get estimation points of KDE
theft.kde.est.points.month = get.grid.points(city.boundary, kde.resolution.meters, FALSE) 

# run KDE
theft.kde.est.month <- c()
for (i in 2:10) {
  theft.kde.est.month[[i]] = run.spatial.kde(theft.kde.sample.points.month[[i]], theft.kde.est.points.month, 1000) 
}

# multiple plots - 2 rows, 5 columns
par(mfrow=c(2,5))
# plot KDE 
for (i in 2:10) {
  plot.spatial.kde(theft.kde.est.month[[i]], theft.kde.est.points.month)  
}

# for evaluation, remove kde estimation points that are not within the city boundary
theft.sp.kde.est.points.month = SpatialPoints(theft.kde.est.points.month, proj4string=CRS("+init=epsg:26971"))
theft.kde.est.points.in.boundary.month = !is.na(over(theft.sp.kde.est.points.month, city.boundary)$OBJECTID)
theft.kde.est.points.month = theft.kde.est.points.month[theft.kde.est.points.in.boundary.month,]
for (i in 2:10) {
  theft.kde.est.month[[i]] = theft.kde.est.month[[i]][theft.kde.est.points.in.boundary.month]
}

# build prediction matrix for evaluation -- requires x, y, and threat columns
theft.prediction.month <- c()
for (i in 2:10) {
  theft.prediction.month[[i]] = as.data.frame(cbind(theft.kde.est.points.month, theft.kde.est.month[[i]]))
  names(theft.prediction.month[[i]]) = c("x","y","threat")
}

# build crime point matrix for evaluation -- evaluate on months 3-11
theft.eval.crime.points.month <- c()
for (i in 3:11) {
  if (nrow(theft.post.month[[i]]) > 0) {
    theft.eval.crime.points.month[[i]] = filter.points.to.boundary(theft.post.month[[i]][,c("x","y")], "+init=epsg:26971", city.boundary)
    theft.eval.crime.points.month[[i]] = as.data.frame(theft.eval.crime.points.month[[i]])
    names(theft.eval.crime.points.month[[i]]) = c("x","y")
  }
  else {
    theft.eval.crime.points.month[[i]] = "NA"
  }
}

# multiple plots - 2 rows, 5 columns
par(mfrow=c(2,5))
for (i in 2:10) {
  plot.surveillance.curve(theft.prediction.month[[i]], theft.eval.crime.points.month[[i+1]], kde.resolution.meters)
}
