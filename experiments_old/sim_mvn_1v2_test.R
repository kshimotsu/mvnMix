rm(list=ls())

setwd('~/Dropbox/mv_normal/R')

library(Rcpp)
library(RcppArmadillo)
library(mvnMix)
library(parallel)
source('~/Dropbox/mv_normal/R/mvn_plrt.R')
# source('~/Dropbox/mv_normal/R/mvn_funcs.R')
# source('~/Dropbox/mv_normal/R/mvn_tests.R')
# source('~/Dropbox/mv_normal/R/mvn_methods.R')
# sourceCpp("cppMVNmixPMLE.cpp")

# library(parallel)
# library(normalregMix)
# library(mixtools)

d <- 2
m <- 1
n <- 200
alpha <- c(0.3,0.7)
mu <- matrix(c(-1,1,0,0),nrow=2,ncol=2)
sigma <- cbind(diag(2),diag(2))

nrep <- 10
nbtsp <- 199

outfilename <- "mvn_test.RData"

ncores <- detectCores()
cl <- makePSOCKcluster(rep('localhost',ncores),master='localhost')

# sink("mvn.out", append=T)

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
    clusterEvalQ(cl, library(mvnMix))
    clusterEvalQ(cl, source('~/Dropbox/mv_normal/R/mvn_plrt.R'))
    # clusterEvalQ(cl, source('~/Dropbox/normal_mixture/R/mvn_funcs.R'))
    # clusterEvalQ(cl, source('~/Dropbox/normal_mixture/R/mvn_tests.R'))
    # clusterEvalQ(cl, source('~/Dropbox/normal_mixture/R/mvn_other_funcs.R'))
    # clusterEvalQ(cl, sourceCpp("cppMVNmixPMLE.cpp"))
    # parallel <- 0
    # clusterExport(cl,varlist=c("y","m","parallel"))
    clusterExport(cl,varlist=c("y","m","nbtsp"))
    # emout1 <- parLapply(cl,1:nrep, function(j) mvnmixPMLE(y=y[,,j], m=m,))
    # coefsum1 <- t(sapply(emout1,"[[","coefficients"))
    # logliksum1 <- t(sapply(emout1,"[[","loglik"))
    emout <- parLapply(cl,1:nrep, function(j) mvnmixMEMtest(y=y[,,j], m=m,
                                                            crit.method="none", parallel=0))
    # emout <- parLapply(cl,1:nrep, function(j) mvnmixMEMtest(y=y[,,j], m=m,
    #             crit.method="boot", nbtsp=nbtsp, parallel=0))
    plrtout <- parLapply(cl,1:nrep, function(j) mvnmixPLRT(y=y[,,j], m=m,
                                                          crit.method="none", nbtsp=nbtsp, parallel=0))
    # plrtout <- parLapply(cl,1:nrep, function(j) mvnmixPLRT(y=y[,,j], m=m,
    #                                                         crit.method="boot", nbtsp=nbtsp, parallel=0))
    # emout.0 <- lapply(1:nrep, function(j) mvnmixMEMtest(y=y[,,j], m=m, crit.method="boot", nbtsp=9, parallel=0))
    emstatsum <- t(sapply(emout,"[[","emstat"))
    plrtsum <- t(sapply(plrtout,"[[","plrtstat"))
    ll0emsum <- t(sapply(emout,"[[","ll0"))
    ll0plsum <- t(sapply(plrtout,"[[","ll0"))
    ll0compare <- cbind(t(ll0emsum),t(ll0plsum))

    ll1emsum <- t(sapply(emout,"[[","ll1"))
    ll1plsum <- t(sapply(plrtout,"[[","ll1"))
    ll1compare <- cbind(ll1emsum[,3],t(ll1plsum))
    # pvalsum <- t(sapply(emout,"[[","pvals"))
    plrtpvalsum <- t(sapply(plrtout,"[[","pvals"))
    # print(pvalsum)
    # rejfreq5 <- 100*colMeans(pvalsum < 0.05)
    # rejfreq1 <- 100*colMeans(pvalsum < 0.01)
    plrtrejfreq5 <- 100*colMeans(plrtpvalsum < 0.05)
    plrtrejfreq1 <- 100*colMeans(plrtpvalsum < 0.01)

    outc <- cbind(emstatsum[,1],t(plrtsum))
    
    time_end <- proc.time()
    runtime  <- time_end - time_start

    # print("DGP, nobs, nrep")
    # print(c(DGP, nobs, nrep))
    # print(runtime)
    # print("rejfreq5, rejfreq1")
    # print(c(rejfreq5,rejfreq1))

    # save.image(file = outfilename)
  # system("mail -s 1v2_report kshimotsu@gmail.com < 1v2_output.out");
} # end of DGP loop

save.image(file = outfilename)

stopCluster(cl)

# sink()
