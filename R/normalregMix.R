normalregMixtest.env <- new.env(parent = emptyenv())
normalregMixtest.env$normalregMix.test.on <- FALSE
normalregMixtest.env$normalregMix.test.seed <- 8888577

#' @useDynLib normalregMix
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @import stats utils doParallel RcppArmadillo foreach minpack.lm iterators parallel ggplot2
NULL