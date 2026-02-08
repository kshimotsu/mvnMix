#' @keywords internal
"_PACKAGE"

#' @useDynLib mvnMix, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom mixtools rmvnorm
#' @importFrom parallel detectCores makeCluster parLapply stopCluster clusterSetRNGStream
#' @importFrom stats pchisq qchisq rnorm runif var
NULL
