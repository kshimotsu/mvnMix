#' @description Generates multivariate mixed normal random variables
#' @export
#' @title rmvnmix
#' @name rmvnmix
#' @param n The number of observations
#' @param alpha m by 1 vector that represents proportions of components
#' @param mu d by m matrix that represents mu
#' @param sigma d by d*m matrix that represents variance of components
#' @return n by d vector
rmvnmix <- function(n, alpha, mu, sigma){
  m <- length(alpha)
  d <- nrow(mu)
  Ind <- sample((1:m), n, replace=TRUE, prob=alpha)
  y <- matrix(0, nrow=n, ncol=d)

  for (j in (1:m)){
    nj <- sum(Ind==j)
    muj <- mu[,j]
    sigmaj <- sigma[,(d*(j-1)+1):(d*j)]
    yj <- rmvnorm(nj, mu = muj, sigma = sigmaj)
    y[Ind==j,] <- yj
  }
y

}

#' @description Convert a vech (lower triangular, column-major) vector to a symmetric matrix
#' @export
#' @title sigmavec2mat
#' @name sigmavec2mat
#' @param sigma.vec A vector of length d(d+1)/2 containing the lower triangle of a symmetric matrix
#' @param d The dimension of the matrix
#' @return d by d symmetric matrix
sigmavec2mat <- function(sigma.vec, d){
# sigma.vec is a vector of length d(d+1)/2
  sigma <- diag(d)
  sigma[lower.tri(sigma, diag=TRUE)] <- sigma.vec
  sigma <- t(sigma) + sigma
  diag(sigma) <- diag(sigma)/2
  sigma
} # end function sigmavec2mat


