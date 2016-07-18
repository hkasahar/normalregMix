data(faithful)
attach(faithful)
normalmixMEMtestSeq(y = eruptions)

regmixPMLE(y = eruptions, x = waiting, m = 1)
regmixPMLE(y = eruptions, x = waiting, m = 2)
normalmixPMLE(y = eruptions, m = 1)
normalmixPMLE(y = eruptions, m = 2)

out<-normalmixPMLE(y = eruptions, m = 2)
summary(out)