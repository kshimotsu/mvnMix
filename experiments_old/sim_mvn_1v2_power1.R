rm(list=ls())

setwd('~/Dropbox/mv_normal/R')

# library(Rcpp)
# library(RcppArmadillo)
library(mvnMix)
library(parallel)
library(mixtools)
# library(normalregMix)

outfilename <- "mvn_power1.RData"
sink_name <- "sim_mvn_1v2_power1.out"

m <- 1
d <- 2

# alpha <- c(0.3,0.7)
# mu <- matrix(c(-1,0,1,0),nrow=2,ncol=2)
# sigma <- cbind(diag(2),diag(2))
alphaset <- array(c(0.3,0.7), dim=c(2,1,2))
muset <- array(0,dim=c(2,2,2))
muset[,,1] <- matrix(c(-1,0,1,0),nrow=2,ncol=2)
muset[,,2] <- matrix(c(-0.5,0,0.5,0),nrow=2,ncol=2)
sigmaset <- array(0,dim=c(2,4,2))
sigmaset[,,1] <- cbind(diag(2),diag(2))
sigmaset[,,2] <- cbind(matrix(c(1,0,0,5),nrow=2,ncol=2),diag(2))

ninits <- 10
nrep <- 1000
nbtsp <- 199
nobsset <- c(200,400)
nnobs <- length(nobsset)

ncores <- detectCores()
cl <- makePSOCKcluster(rep('localhost',ncores),master='localhost')

sink(sink_name, append=T)

DGPset <- c(1:2)
ndgp <- length(DGPset)

rejfreq5all <- matrix(0, nrow=ndgp, ncol=6)
rejfreq1all <- matrix(0, nrow=ndgp, ncol=6)
out.all <- vector('list',length=2*ndgp)

for (DGP in DGPset){

  alpha <- alphaset[,,DGP]
  mu <- muset[,,DGP]
  sigma <-sigmaset[,,DGP]

  for (inobs in 1:nnobs)  {
    time_start <- proc.time()
    n <- nobsset[inobs]
    set.seed(123456)
    y <- replicate(nrep, rmvnmix(n,alpha,mu,sigma))
    clusterSetRNGStream(cl,123456)
    # clusterEvalQ(cl, library(normalregMix))
    clusterEvalQ(cl, library(mvnMix))
    clusterEvalQ(cl, library(mixtools))
    clusterEvalQ(cl, library(parallel))
    clusterEvalQ(cl, library(Rcpp))
    # clusterEvalQ(cl, library(RcppArmadillo))
    clusterExport(cl,varlist=c("y","m","nbtsp","ninits"))

    out <- parLapply(cl,1:nrep, function(j) mvnmixMEMtest(y=y[,,j], m=m,
                                                          ninits=ninits, crit.method="boot", nbtsp=nbtsp, parallel=0))

    out.all[[2*(DGP-1)+inobs]] <- out

    pval <- t(sapply(out,"[[","pvals"))
    rejfreq <- 100*rbind(colMeans(pval < 0.10),colMeans(pval < 0.05),colMeans(pval < 0.01))
    rejfreq5 <- rejfreq[2,]
    rejfreq1 <- rejfreq[3,]

    rejfreq5all[DGP,(3*inobs-2):(3*inobs)] <- rejfreq5
    rejfreq1all[DGP,(3*inobs-2):(3*inobs)] <- rejfreq1

    time_end <- proc.time()
    runtime  <- time_end - time_start

    print("DGP, n, nrep")
    print(c(DGP, n, nrep))
    print(runtime)
    print("rejfreq5, rejfreq1")
    print(c(rejfreq5,rejfreq1))

  } # end of inobs loop

} # end of DGP loop

stopCluster(cl)

rm(y)
save.image(file = outfilename)

sink()
