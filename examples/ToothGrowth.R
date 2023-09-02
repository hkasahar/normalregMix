data(faithful)
attach(faithful)
library(normalregMix)
set.seed(1234)

normalmixPMLE(y = waiting, m = 2)
normalmixMEMtestSeq(y = waiting, maxm = 4)

data(ToothGrowth)
attach(ToothGrowth)
library(normalregMix)
len.log <- log(len)
dose.log <- log(dose)

set.seed(1234)
normalmixMEMtest(y = len, x = dose.log, m = 1, crit.method = 'asy')
normalmixMEMtest(y = len, x = dose.log, m = 2, crit.method = 'asy')

set.seed(1234)
normalmixMEMtestSeq(y = len, x = log(dose), maxm = 2)
beepr::beep(3)
normalmixPMLE(y = len.log, x = dose.log, m = 1)

model.m2 <- regmixPMLE(y = len, x = dose.log, m = 2)
plot.m2.df <- data.frame(y = len.log, x = dose.log, 
                         components = model.m2$components)
plot.m2 <- ggplot(plot.m2.df, aes(x, y)) + 
  geom_point(aes(colour = factor(components))) +
  scale_fill_discrete(name="Component")
plot.m2
#https://rstudio-pubs-static.s3.amazonaws.com/30490_f5f75fae348843819993bf1e0957f4f1.html
#https://rpubs.com/daniambrosio/tooth_growth_exploratory_data_analysis

supp.estimated <- as.character(model.m2$components)


avg <- aggregate(len.log~.,data=tooth.growth.df,mean)
g <- ggplot(aes(x=dose, y = len.log), data = tooth.growth.df) + 
  geom_point(aes(color = supp)) 
g <- g + geom_line(data=avg,aes(group=supp,colour=supp))
print(g)



tooth.growth.df <- data.frame(supp, len.log = log(len.log))
tooth.growth.df$dose.log <- log(dose)
avg <- aggregate(len.log~.,data=tooth.growth.df,mean)
g <- ggplot(aes(x=dose.log, y = len.log), data = tooth.growth.df) + 
  geom_point(aes(color = supp)) 
g <- g + geom_line(data=avg,aes(group=supp,colour=supp))
print(g)

tooth.growth.df <- data.frame(dose, len, treatment.estimated = as.character(model.m2$components))
treatment.estimated = as.character(model.m2$components)
ggplot(aes(x = treatment.estimated, y = len), data = tooth.growth.df) +
  geom_boxplot(aes(fill = treatment.estimated)) + facet_wrap(~ dose)

ggplot(aes(x = supp, y = len), data = ToothGrowth) +
  geom_boxplot(aes(fill = supp)) + facet_wrap(~ dose)
