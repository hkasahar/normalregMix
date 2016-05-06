rm(list=ls())
setwd("~/Desktop")

install.packages("minpack.lm")	# required for simulating the critcal values in regression mixtures
install.packages("normalregmix.tar.gz", repos = NULL, type = "source")
library(minpack.lm)
library(normalregMix)

data(faithful)
attach(faithful)

result2 <- normalmixPMLE(eruptions, m=2)
summary(result2)
normalmixMEMtest(eruptions, m=2)

result3 <- normalmixPMLE(eruptions, m=3)
summary(result3)
normalmixMEMtest(eruptions, m=3)

result4 <- normalmixPMLE(eruptions, m=4)
summary(result4)
normalmixMEMtest(eruptions, m=4)

# The following example shows that the modified EM
# test has a problem when alpha is small.

resultreg1 <- regmixPMLE(y=waiting, x=eruptions, m=1)
summary(resultreg1)
resultreg2 <- regmixPMLE(y=waiting, x=eruptions, m=2)
summary(resultreg2)

regmixMEMtest(y=waiting, x=eruptions, m=1, crit.method = "asy")
regmixMEMtest(y=waiting, x=eruptions, m=1, crit.method = "boot")

