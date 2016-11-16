data(faithful)
attach(faithful)
library(normalregMix)
library(ggplot2)
regmixMEMtestSeq(y = eruptions, x = waiting)
regmixPMLE(y = eruptions, x = waiting, m = 1)
model.m2 <- regmixPMLE(y = Sepal.Length, x = waiting, m = 3)
plot.m2.df <- data.frame(y = eruptions, x = waiting, 
                        components = model.m2$components)
plot.m2 <- ggplot(plot.m2.df, aes(x, y)) + 
            geom_point(aes(colour = factor(components))) +
            scale_fill_discrete(name="Component")
plot.m2

regmixPMLE(y = eruptions, x = waiting, m = 3)

normalmixPMLE(y = eruptions, m = 1)
normalmixPMLE(y = eruptions, m = 2)

out<-normalmixPMLE(y = eruptions, m = 2)
summary(out)