rm(list=ls())

setwd('~/Dropbox/mv_normal/R')

library(Rcpp)
library(RcppArmadillo)
library(mvnMix)
library(parallel)
# source('~/Dropbox/mv_normal/R/mvn_plrt.R')

d <- 2
m <- 2
n <- 200
alpha <- c(0.3,0.7)
mu <- matrix(c(-1,0,1,0),nrow=2,ncol=2)
sigma <- cbind(diag(2),diag(2))

nrep <- 1000

outfilename <- "pmle_test.RData"

ncores <- detectCores()
cl <- makePSOCKcluster(rep('localhost',ncores),master='localhost')

sink("pmle_test.out", append=T)

DGPset <- c(1:1)
for (DGP in DGPset)
{
    time_start <- proc.time()
    set.seed(123456)
    clusterEvalQ(cl,set.seed(123456))
    # y <- rmvnmix(n,alpha,mu,sigma)
    if (m==1){
      y <- array(rnorm(n*d*nrep),dim=c(n,d,nrep))
    } else{
      y <- replicate(nrep, rmvnmix(n,alpha,mu,sigma))
    }
     # clusterEvalQ(cl, library(normalregMix))
    clusterEvalQ(cl, library(parallel))
    clusterEvalQ(cl, library(mixtools))
    clusterEvalQ(cl, library(Rcpp))
    clusterEvalQ(cl, library(RcppArmadillo))
    clusterEvalQ(cl, library(mvnMix))
    # clusterEvalQ(cl, source('~/Dropbox/mv_normal/R/mvn_plrt.R'))
    clusterExport(cl,varlist=c("y","m"))
    pmleout <- parLapply(cl,1:nrep, function(j) mvnmixPMLE(y=y[,,j], m=m))
    coefsum <- t(sapply(pmleout,"[[","coefficients"))
    logliksum <- t(sapply(pmleout,"[[","loglik"))

    time_end <- proc.time()
    runtime  <- time_end - time_start

} # end of DGP loop

rm(y)

stopCluster(cl)

sigma1 <- sigma[,c(1,2)]
sigma1vec <- sigma1[lower.tri(sigma1, diag=TRUE)]
sigma2 <- sigma[,c(3,4)]
sigma2vec <- sigma2[lower.tri(sigma2, diag=TRUE)]

bias <- colMeans(coefsum) - c(alpha,mu,sigma1vec,sigma2vec)

stdev <- apply(coefsum,2,sd)

save.image(file = outfilename)

sink()
