# ===================================================================
# Replicate experiment in 8.3 section of Kasahara and Shimotsu (2015)
# to test the number of regimes in cross-country growth regression
# from Mankiew, Romer, and Weil (1992) based on the Solow-Swan model 
# ===================================================================
#install.packages("ggmap")
library(maps)
library(ggmap)
library(normalregMix)
library(ggplot2)

## Plots the world map with dots on countries whose colours are determined
## based on components they belong to 
PlotWorldMap <- function(map.data, m)
{
  if ((m != 2) && (m != 3))
    return (NULL)
  if (m == 2)
  {
    countries.one <- map.data[map.data$modelm2 == 1, ]
    countries.two <- map.data[map.data$modelm2 == 2, ]
    countries.three <- NULL
  }
  if (m == 3)
  {
    countries.one <- map.data[map.data$modelm3 == 1, ]
    countries.two <- map.data[map.data$modelm3 == 2, ]
    countries.three <- map.data[map.data$modelm3 == 3, ]
  }
  map.countries.one <- geom_point(aes(x=countries.one$LON, y=countries.one$LAT, label = countries.one$NAME), 
                                  color="red", size=4) 
  map.countries.two <- geom_point(aes(x=countries.two$LON, y=countries.two$LAT, label = countries.two$NAME), 
                                  color="green", size=4) 
  map.countries.three <- geom_point(aes(x=countries.three$LON, y=countries.three$LAT, label = countries.three$NAME), 
                                    color="blue", size=4) 
  
  map.world <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
  map <- ggplot() + map.world
  map <- map + map.countries.three
  map <- map + map.countries.two
  map <- map + map.countries.one 
  map 
}

## Plots the names of countries to the map for testing
PlotTextMap <- function(map.data, m)
{
  if ((m != 2) && (m != 3))
    return (NULL)
  if (m == 2)
  {
    countries.one <- map.data[map.data$modelm2 == 1, ]
    countries.two <- map.data[map.data$modelm2 == 2, ]
    countries.three <- NULL
  }
  if (m == 3)
  {
    countries.one <- map.data[map.data$modelm3 == 1, ]
    countries.two <- map.data[map.data$modelm3 == 2, ]
    countries.three <- map.data[map.data$modelm3 == 3, ]
  }
  map.countries.one <- geom_point(aes(x=countries.one$LON, y=countries.one$LAT), 
                                  color="red", size=4) 
  map.countries.two <- geom_point(aes(x=countries.two$LON, y=countries.two$LAT), 
                                  color="green", size=4) 
  map.countries.three <- geom_point(aes(x=countries.three$LON, y=countries.three$LAT), 
                                    color="blue", size=4) 
  
  map.world <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
  map <- ggplot() + map.world
  map <- map + map.countries.three
  map <- map + map.countries.two
  map <- map + map.countries.one 
  
  ggplot(map.data, aes(x= LON, y= LAT, colour="green", label=NAME)) +
     map.world + geom_point() +geom_text(aes(label=NAME),hjust=0, vjust=0) 
}

## Returns the longitude and latitude given a vector of characters of country names
GetLocation <- function(names) {
  df <- geocode(as.vector(names))
  colnames(df) <- c("LON", "LAT")
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- read.csv("MRWDataNameLoc.csv") # read the data (assuming it's in the same dir)
data.use <- data[data$INTERMEDIATE == 1, ] # only non-oil
# data.use <- data[data$INTERMEDIATE == 1, ] # only intermediate
# data.use <- data.use[complete.cases(data.use),] # filter NULL 

# model specification
y <- log(data.use$GDP85)
x <- cbind(IOVERY = log(data.use$IOVERY),
           POPGRO = log(data.use$POPGRO + 0.05),
           SCHOOL = log(data.use$SCHOOL))
# x <- cbind(GDP60 = log(data.use$GDP60),
#            IOVERY = log(data.use$IOVERY),
#            POPGRO = log(data.use$POPGRO + 0.05),
#            LIT60 = log(data.use$LIT60),
#            SCHOOL = log(data.use$SCHOOL),
#            NONOIL = data.use$NONOIL)

# Test H_0: m = 1 vs H_1: m = 2
# model.m2 <- regmixMEMtest(y = y, m = 1, x = x, parallel = 1, crit.method = "boot")  # took 6 mins in my i5 laptop
# Test H_0: m = 2 vs H_1: m = 3 
# model.m3 <- regmixMEMtest(y = y, m = 2, x = x, parallel = 1, crit.method = "boot")

# print(model.m2)
# print(model.m3)

## Draw the Map
# PMLE assuming m = 2
parlist.m2 <- regmixPMLE(y = y, x = x, m = 2)
parlist.m3 <- regmixPMLE(y = y, x = x, m = 3)
map.data <- as.data.frame(cbind(y, x))
map.data$modelm2 <- parlist.m2$components
map.data$modelm3 <- parlist.m3$components
map.data <- cbind(map.data, 
                  NAME = data.use$NAME, LON = data.use$LON, LAT =  data.use$LAT)

PlotWorldMap(map.data, m = 2)

## Draw residual plots (for m = 2)
y1 <- map.data[map.data$modelm2 == 1, ]$y
y2 <- map.data[map.data$modelm2 == 2, ]$y
x1 <- map.data[map.data$modelm2 == 1, 2:4]
x2 <- map.data[map.data$modelm2 == 2, 2:4]
#http://www.nuffield.ox.ac.uk/teaching/economics/bond/mrw.pdf
# On component 1
plotDiag(parlist.m1$components, y1, x1, 1)
plotDiag(parlist.m1$components, y1, x1, 2)
plotDiag(parlist.m1$components, y1, x1, 3)
# On component 2
plotDiag(parlist.m1$components, y2, x2, 1)
plotDiag(parlist.m1$components, y2, x2, 2)
plotDiag(parlist.m1$components, y2, x2, 3)
