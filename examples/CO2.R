data(CO2)
attach(CO2)
library(normalregMix)
library(ggplot2)

treatment.bin <- sapply(Treatment, switch, "chilled" = 1, "unchilled" = 0)
plant.identifier <- as.numeric(Plant)
plant.identifier.matrix <- matrix(0, 
                                  ncol = max(plant.identifier),
                                  nrow = length(plant.identifier))
for (i in 1:length(plant.identifier))
  plant.identifier.matrix[i, plant.identifier[i]] <- 1
model.m2 <- regmixPMLE(y = uptake, 
                       x = cbind(conc, plant.identifier.matrix), 
                       m = 2)
plot.m2.df <- data.frame(y = uptake, x = conc, 
                         components = model.m2$components)
plot.m2 <- ggplot(plot.m2.df, aes(x, y)) + 
  geom_point(aes(colour = factor(components))) +
  scale_fill_discrete(name="Component")
plot.m2
