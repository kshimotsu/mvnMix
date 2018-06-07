rm(list=ls())

library(Rcpp)
library(RcppArmadillo)

setwd('~/Dropbox/mv_normal/R')
source('~/Dropbox/mv_normal/R/mvn_funcs.R')
source('~/Dropbox/mv_normal/R/mvn_tests.R')
source('~/Dropbox/mv_normal/R/mvn_methods.R')
sourceCpp("cppMVNmixPMLE.cpp")

library(parallel)
# library(normalregMix)
library(mixtools)

d <- 2
m <- 2
n <- 200
alpha <- c(0.3,0.7)
mu <- matrix(c(-1,1,0,0),nrow=2,ncol=2)
sigma <- cbind(diag(2),diag(2))

m <- 1

nrep <- 1000

outfilename <- "mvn.RData"

ncores <- detectCores()
cl <- makePSOCKcluster(rep('localhost',ncores),master='localhost')

sink("mvn.out", append=T)

# rejfreq5all <- matrix(0, nrow=length(DGPset), ncol=3*nnobs)
# rejfreq1all <- matrix(0, nrow=length(DGPset), ncol=3*nnobs)

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
    clusterEvalQ(cl, source('~/Dropbox/normal_mixture/R/mvn_funcs.R'))
    clusterEvalQ(cl, source('~/Dropbox/normal_mixture/R/mvn_tests.R'))
    clusterEvalQ(cl, source('~/Dropbox/normal_mixture/R/mvn_other_funcs.R'))
    clusterEvalQ(cl, sourceCpp("cppMVNmixPMLE.cpp"))
    # parallel <- 0
    # clusterExport(cl,varlist=c("y","m","parallel"))
    nbtsp <- 199
    clusterExport(cl,varlist=c("y","m","nbtsp"))
    # emout1 <- parLapply(cl,1:nrep, function(j) mvnmixPMLE(y=y[,,j], m=m,))
    # coefsum1 <- t(sapply(emout1,"[[","coefficients"))
    # logliksum1 <- t(sapply(emout1,"[[","loglik"))
    # emout <- parLapply(cl,1:nrep, function(j) mvnmixMEMtest(y=y[,,j], m=m,
    #                                                         crit.method="none", parallel=0))
    emout <- parLapply(cl,1:nrep, function(j) mvnmixMEMtest(y=y[,,j], m=m,
                crit.method="boot", nbtsp=nbtsp, parallel=0))
    # emout.0 <- lapply(1:nrep, function(j) mvnmixMEMtest(y=y[,,j], m=m, crit.method="boot", nbtsp=9, parallel=0))
    pvalsum <- t(sapply(emout,"[[","pvals"))
    print(pvalsum)
    rejfreq5 <- 100*colMeans(pvalsum < 0.05)
    rejfreq1 <- 100*colMeans(pvalsum < 0.01)
    # rejfreq5all[DGP,(3*inobs-2):(3*inobs)] <- rejfreq5
    # rejfreq1all[DGP,(3*inobs-2):(3*inobs)] <- rejfreq1

    time_end <- proc.time()
    runtime  <- time_end - time_start

    # print("DGP, nobs, nrep")
    # print(c(DGP, nobs, nrep))
    # print(runtime)
    # print("rejfreq5, rejfreq1")
    # print(c(rejfreq5,rejfreq1))

    save.image(file = outfilename)
  # system("mail -s 1v2_report kshimotsu@gmail.com < 1v2_output.out");
} # end of DGP loop

rm(y)
save.image(file = outfilename)

stopCluster(cl)

# sigma1 <- sigma[,c(1,2)]
# sigma1vec <- sigma1[lower.tri(sigma1, diag=TRUE)]
# sigma2 <- sigma[,c(3,4)]
# sigma2vec <- sigma2[lower.tri(sigma2, diag=TRUE)]
#
# bias <- colMeans(coefsum1) - c(alpha,t(mu),sigma1vec,sigma2vec)

# return(list(outall=outall,runtime=runtime))

sink()
