# ===================================================================
# Replicate experiment in 8.3 section of Kasahara and Shimotsu (2015)
# to test the number of regimes in cross-country growth regression
# from Mankiew, Romer, and Weil (1992) based on the Solow-Swan model 
# ===================================================================
#install.packages("ggmap")
library(ggmap)
library(normalregMix)
library(ggplot2)

## Plots for residuals where pivot is the index of pivotted column
PlotRes <- function(y, x, pivot) {
  pivot.name <- as.character(pivot)
  if (!is.null(colnames(x)))
    pivot.name <- colnames(x)[pivot]
  ivs <- as.matrix(x)
  ivs.pivot <- ivs[,pivot]
  ivs.others <- ivs[,-pivot]
  lm.y.other <- lm(y ~ ivs.others)
  lm.pivot.other <- lm(ivs.pivot ~ ivs.others)
  plot.df <- data.frame(y.on.others = lm.y.other$residuals, 
                        pivot.on.others = lm.pivot.other$residuals) 
  plot <- ggplot(plot.df, aes(x=pivot.on.others, y=y.on.others))
  plot <- plot + geom_point(shape=1) + geom_smooth(method=lm) +
          xlab(paste("pivot.on.others (pivot on ", pivot.name, ")", sep = ""))
  plot
}

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

setwd("C:\\Users\\chiyahn\\Dropbox\\Work\\June06\\normalregMix\\experiments\\ApplicationMRW")

data <- read.csv("MRWDataNameLoc.csv") # read the data (assuming it's in the same dir)
data.use <- data[data$INTER == 1, ] # only intermediate
#data.use <- data[complete.cases(data),] # filter NULL only 

# model specification
y <- log(data.use$GDP85 / data.use$GDP60)
x <- cbind(GDP60 = log(data.use$GDP60),
           IONY = log(data.use$IONY),
           POPGRO = log(data.use$POPGRO + 0.05),
           SCHOOL = log(data.use$SCHOOL),
           LIT60 = log(data.use$LIT60),
           NONOIL = data.use$NONOIL)

# Test H_0: m = 1 vs H_1: m = 2
model.m2 <- regmixMEMtest(y = y, m = 1, x = x, parallel = TRUE, crit.method = "asy") # took 6 mins in my i5 laptop
# Test H_0: m = 2 vs H_1: m = 3 
# (used bootstrapping instead for this case. Check for number of bootstraps for actual experiment)
# Test for H_0: m = 2 has been rejected (with p* = 0.05) when other than intermediate countries are also included. 
model.m3 <- regmixMEMtest(y = y, m = 2, x = x, parallel = TRUE, crit.method = "boot", nbtsp = 200) 
# Test for H_0: m = 3 has been rejected (p = 0.30) when other than intermediate countries are also included. 
model.m4 <- regmixMEMtest(y = y, m = 3, x = x, parallel = TRUE, crit.method = "boot", nbtsp = 200)

## Draw the Map
# PMLE assuming m = 2
parlist.m2 <- regmixPMLE(y = y, x = x, m = 2, vcov.method = "OPG")
parlist.m3 <- regmixPMLE(y = y, x = x, m = 3, vcov.method = "OPG")
map.data <- as.data.frame(cbind(y, x))
map.data$modelm2 <- parlist.m2$indices
map.data$modelm3 <- parlist.m3$indices
map.data <- cbind(map.data, 
                  NAME = data.use$NAME, LON = data.use$LON, LAT =  data.use$LAT)
print(map.data)

PlotWorldMap(map.data, m = 2)
PlotWorldMap(map.data, m = 3)

## Draw residual plots (for m = 2)
y1 <- map.data[map.data$modelm2 == 1, ]$y
y2 <- map.data[map.data$modelm2 == 2, ]$y
x1 <- map.data[map.data$modelm2 == 1, 2:4]
x2 <- map.data[map.data$modelm2 == 2, 2:4]

# On component 1
PlotRes(y1,x1,1)
PlotRes(y1,x1,2)
PlotRes(y1,x1,3)
# On component 2
PlotRes(y2,x2,1)
PlotRes(y2,x2,2)
PlotRes(y2,x2,3)
