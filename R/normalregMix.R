#' @description Environment that determines whether \code{normalregMix}
#' is in the test mode.
#' @export
#' @title normalregMixtest.env
#' @name normalregMixtest.env
#' @format An object of class environment with two items. The first element, 
#' \code{normalregMix.test.on}, is logical and sets whether the test mode is on.
#' The second element, \code{normalregMix.test.seed}, is an integer that is the
#' seed for random number generation.
normalregMixtest.env <- new.env(parent = emptyenv())
normalregMixtest.env$normalregMix.test.on <- FALSE
normalregMixtest.env$normalregMix.test.seed <- 8888577

#' @useDynLib normalregMix
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @import stats utils RcppArmadillo minpack.lm parallel
NULL