## Test environments
* local OS X install, R 3.3.1
* local Windows 10 Pro install, R 3.3.1

## R CMD check results
There were no ERRORs or WARNINGs. 
There was 1 NOTE:

* checking R code for possible problems ... NOTE
Found the following possibly unsafe calls:
File 'normalregMix/R/methods.R':
  unlockBinding("normalregMix.test.on", getNamespace("normalregMix"))
  unlockBinding("normalregMix.test.seed", getNamespace("normalregMix"))

  testMode uses global variables "normalregMix.test.on" and "normalregMix.test.seed" to set random seed uniformly when random number generation occurs. This is intended for reproducibility of experiments done with this package. 