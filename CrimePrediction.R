#########################################

## PREDICTIVE MODELS ##
## AUTHOR: MARGARET FURR ##

#########################################

source("CrimeUtil.R") 

# read chicago boundary
city.boundary = read.shapefile("City_Boundary/City_Boundary", "poly", "+init=epsg:3435", "+init=epsg:26971")

# set prediction resolution
prediction.resolution.meters = 200

# get crime samples
theft.month.to.month = sample.crime("2014_THEFT.csv", -1, 2, 12)
theft.prior.month <- c()
for (i in 2:9) {
  theft.prior.month[[i]] = theft.month.to.month[theft.month.to.month$month==i,]
}
theft.month <- c()
for (i in 2:9) {
  theft.month[[i]] = theft.month.to.month[theft.month.to.month$month==i+1,]
}
theft.post.month <- c()
for (i in 2:9) {
  theft.post.month[[i]] = theft.month.to.month[theft.month.to.month$month==i+2,]
}

##### train model on responses from month2 data, using predictors from month1. #####

# get negative observations within chicago
non.crime.points = cbind(0, get.grid.points(city.boundary, prediction.resolution.meters, TRUE))
names(non.crime.points)[1] = "response"

# get positive observations from month (months 3-10) within chicago
training.crime.points <- c()
for (i in 2:9) {
  training.crime.points[[i]] = cbind(1, theft.month[[i]][,c("x","y")])
  names(training.crime.points[[i]])[1] = "response"
}

# combine positive and negative points
training.data <- c()
for (i in 2:9) {
  training.data[[i]] = rbind(non.crime.points, training.crime.points[[i]])
}

# calculate crime density for each training observation based on prior.month (months 2-9) records
theft.density <- c()
for (i in 2:9) {
  theft.density[[i]] = run.spatial.kde(theft.prior.month[[i]][,c("x","y")],training.data[[i]][,c("x","y")], max.sample.size=1000)
}

# theft - calculate distances to CTA l'rail stations
lrail.stations = read.shapefile("CTA_Stations", "points", "+init=epsg:3435", "+init=epsg:26971")@coords  # 290
lrail.stations.min.distance <- c()
for (i in 2:9) {
  lrail.stations.min.distance[[i]] = get.min.distances(training.data[[i]][,c("x","y")], lrail.stations)  # 18919
}

# theft - calculate distances to park art
park.art = read.shapefile("ParkArt", "points", "+init=epsg:3435", "+init=epsg:26971")@coords  # 290
park.art.min.distance <- c()
for (i in 2:9) {
  park.art.min.distance[[i]] = get.min.distances(training.data[[i]][,c("x","y")],park.art)
}

# add predictor columns (theft density and factor.min.distance) to training data 
for (i in 2:9) {
  training.data[[i]] = cbind(training.data[[i]], theft.density[[i]], lrail.stations.min.distance[[i]], park.art.min.distance[[i]])  # 18919
  names(training.data[[i]]) = c("response", "x","y","theft.density","lrail.stations.min.distance","park.art.min.distance")
}

# fit GLM -- don't use x- and y-coordinates as predictors.
glm.fit <- c()
for (i in 2:9) {
  glm.fit[[i]] = glm(response ~. -x -y, data = training.data[[i]], family=binomial)
}

# build dataframe to predict, based on february's data
prediction.points = get.grid.points(city.boundary, prediction.resolution.meters, TRUE)
theft.density.prediction  <- c()
for (i in 2:9) {
  theft.density.prediction[[i]] = run.spatial.kde(theft.month[[i]][,c("x","y")], prediction.points, max.sample.size=1000)
}
# theft - calculate distances to CTA l'rail stations
lrail.stations.min.distance.prediction = get.min.distances(prediction.points, lrail.stations)  
# theft - calculate distances to park art
park.art.min.distance.prediction = get.min.distances(prediction.points,park.art)
# prediction data
prediction.data <- c()
for (i in 2:9) {
  prediction.data[[i]] = data.frame(cbind(prediction.points, theft.density.prediction[[i]], lrail.stations.min.distance.prediction, park.art.min.distance.prediction))  # 18919
  names(prediction.data[[i]]) = c("x","y","theft.density","lrail.stations.min.distance","park.art.min.distance")
}

# run prediction
threats <- c()
for (i in 2:9) {
  threats[[i]] = predict(glm.fit[[i]], prediction.data[[i]], type="response")
}

# build prediction dataframe
theft.prediction <- c()
for (i in 2:9) {
  theft.prediction[[i]] = cbind(prediction.points,threats[[i]])
  names(theft.prediction[[i]]) = c("x","y","threat")
}

# evaluate prediction on post.month crime records
par(mfrow=c(2,5))
for (i in 2:9) {
  plot.surveillance.curve(theft.prediction[[i]], theft.post.month[[i]][,c("x","y")], prediction.resolution.meters)
}

