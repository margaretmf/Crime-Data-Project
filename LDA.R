#########################################

## LDA ##
## AUTHOR: MARGARET FURR ##

#########################################

library(ks); library(rgdal); library(rgeos)
library(RColorBrewer); library(ggplot2); library(ggmap)
library(zoo); library(maptools); library(lubridate)library(XML)
library(tm); library(filehash)
library(topicmodels); library(SnowballC)library(slam)


source("CrimeUtil.R") 

# read chicago boundary
city.boundary = read.shapefile("City_Boundary/City_Boundary", "poly", "+init=epsg:3435", "+init=epsg:26971")

# set prediction resolution
prediction.resolution.meters = 150
# read a sample of tweets
tweets = read.csv("tweets_large_sample.csv", header = TRUE, stringsAsFactors = FALSE)
n.tweets = nrow(tweets)  # 100,000

# edit text
tweets$text = gsub("(RT|via)((?:\\b\\W*@\\w+)+)","",tweets$text)
tweets$text = gsub("http[^[:blank:]]+", "", tweets$text)
tweets$text = gsub("@\\w+", "", tweets$text)
tweets$text = gsub("[ \t]{2,}", "", tweets$text)
tweets$text = gsub("^\\s+|\\s+$", "", tweets$text)
tweets$text= gsub('\\d+', '', tweets$text)
tweets$text = gsub("[[:punct:]]", " ", tweets$text)

# reproject from degrees to meters
tweets.lonlat = cbind(tweets$longitude, tweets$latitude)
meters.locations = project(tweets.lonlat, proj="+init=epsg:26971")
tweets$x = meters.locations[,1]
tweets$y = meters.locations[,2]
tweets = tweets[is.finite(tweets$x) & is.finite(tweets$y),]  # some tweets are posted from outside the bounds of the chicago projection and their reprojections are infinite. remove such tweets.

# filter all tweets to be in the city
points.in.boundary = data.frame(filter.points.to.boundary(tweets[,c("x","y")], "+init=epsg:26971", city.boundary))
tweets = tweets[points.in.boundary$x==tweets$x&points.in.boundary$y==tweets$y,]

# times of tweets - convert strings to dates
tweets$timestamp = ymd_hm(tweets$timestamp)
tweets$month = month(tweets$timestamp)
tweets$hour = hour(tweets$timestamp)
tweets$minute = minute(tweets$timestamp)
tweets$year = year(tweets$timestamp)
tweets$wday = wday(tweets$timestamp)

# plot
plot(city.boundary)
points(tweets[,c("x", "y")], pch = ".")

# points
tweet.points <- tweets[,c("x","y")]

# estimation points of KDE
kde.est.points = get.grid.points(city.boundary, kde.resolution.meters, FALSE) 

# run KDE
kde.est = run.spatial.kde(tweet.points, kde.est.points, 1000) 

# plot KDE
plot.spatial.kde(kde.est, kde.est.points)
# add city boundary to KDE plot
plot(city.boundary, axes=TRUE, border="black", asp=1, add=TRUE)
title(main = "KDE Plot for Sample of Tweets")

set.seed(12345)

# SMALL SAMPLE

# read documents - weekday tweets
# corpus - documents = tweets grouped by days
tweets.text <- tweets$text
tweets.corpus = VCorpus(VectorSource(tweets.text))
inspect(tweets.corpus)

# stopwords
english.stopwords = c(stopwords("english")) 

# clean and compute tfidf - tweets
tweet.corpus.clean = tm_map(tweet.corpus, stripWhitespace)  # remove extra whitespace
tweet.corpus.clean = tm_map(tweet.corpus.clean, removeNumbers)  # remove numbers
tweet.corpus.clean = tm_map(tweet.corpus.clean, removePunctuation)  # remove punctuation
tweet.corpus.clean = tm_map(tweet.corpus.clean, content_transformer(tolower))  # ignore case
tweet.corpus.clean = tm_map(tweet.corpus.clean, stemDocument)  # stem all words
tweet.corpus.clean = tm_map(tweet.corpus.clean, removeWords, english.stopwords) # remove stopwords
inspect(tweet.corpus.clean)

# document term matrix - tweets
tweet.corpus.clean.tf = DocumentTermMatrix(tweet.corpus.clean, control = list(weighting = weightTf))
tweet.corpus.clean.tf

# remove empty documents - tweets 
row.sums.tweets = apply(tweet.corpus.clean.tf, 1, sum)  
tweets.emptied = tweets.text[row.sums.tweets > 0]
tweet.corpus.clean.tf.emptied = tweet.corpus.clean.tf[row.sums.tweets > 0,]
tweet.corpus.clean.tf.emptied

# train topic model - emptied - weekday grouped by day - 5 topics
tweet.corpus.clean.tf.emptied.model = LDA(tweet.corpus.clean.tf.emptied, 5, method = "Gibbs", control=list(iter=100, seed=12345, best=TRUE),pass=100)
terms(tweet.corpus.clean.tf.emptied.model, 7)

