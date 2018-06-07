rm(list=ls())

library(Rcpp)
library(RcppArmadillo)
library(mvnMix)
library(mixtools)
library(parallel)

outfilename <- "mvn_power3.RData"
sink_name <- "sim_mvn_power3.out"

DGPset <- c(2:2)

m <- 2
d <- 2

alphaset <- array(c(0.15,0.35,0.5), dim=c(3,1,2))
muset <- array(0,dim=c(2,3,2))
muset[,,1] <- matrix(c(-2,0,0,2,2,0),nrow=2,ncol=3)
muset[,,2] <- matrix(c(-2,0,-1,4,2,0),nrow=2,ncol=3)
sigmaset <- array(0,dim=c(2,6,2))
sigmaset[,,1] <- cbind(diag(2),diag(2),diag(2))
sigmaset[,,2] <- cbind(diag(2),diag(2),diag(2))
# sigmaset[,,2] <- cbind(diag(2), matrix(c(5,0,0,1),nrow=2,ncol=2), 
#                        matrix(c(1,0,0,5),nrow=2,ncol=2))

ninits <- 10
nrep <- 100
nbtsp <- 19
nobsset <- c(200,400)
# nobsset <- c(200)
nnobs <- length(nobsset)

nworker <- min(nrep,detectCores())
cl <- makePSOCKcluster(nworker)

sink(sink_name, append=T)

rejfreq10all <- matrix(0, nrow=16, ncol=6)
rejfreq5all <- matrix(0, nrow=16, ncol=6)
rejfreq1all <- matrix(0, nrow=16, ncol=6)
out.all <- vector('list',length=32)

for (DGP in DGPset)
{

alpha <- alphaset[,,DGP]
mu <- muset[,,DGP]
sigma <-sigmaset[,,DGP]
  

	for (inobs in 1:nnobs)
	{
		time_start <- proc.time()
		nobs <- nobsset[inobs]
		set.seed(123456)
		Y <- replicate(nrep, rmvnmix(nobs,alpha,mu,sigma))
		# Y <- array(0, dim=c(nobs,d,nrep))
		# for (j in 1:nrep){
		# 	Y[,,j] <- matrix(rnorm(nobs*d),nrow=nobs) %*% chol(sigma)
		# }
		clusterSetRNGStream(cl, 123456)
    clusterEvalQ(cl, library(mvnMix))
    clusterEvalQ(cl, library(mixtools))
    clusterEvalQ(cl, library(Rcpp))
    clusterExport(cl=cl, varlist=c("Y", "m", "nbtsp", "ninits"))

		out <- parLapplyLB(cl,1:nrep, function(j) mvnmixMEMtest(y=Y[,,j],
              m=m, ninits=ninits, crit.method="boot", nbtsp=nbtsp, parallel=0))

		out.all[[2*(DGP-1)+inobs]] <- out

		pval <- t(sapply(out,"[[","pvals"))
		rejfreq <- 100*rbind(colMeans(pval < 0.10),colMeans(pval < 0.05),colMeans(pval < 0.01))
		rejfreq10 <- rejfreq[1,]
		rejfreq5 <- rejfreq[2,]
		rejfreq1 <- rejfreq[3,]

		rejfreq10all[DGP,(3*inobs-2):(3*inobs)] <- rejfreq10
		rejfreq5all[DGP,(3*inobs-2):(3*inobs)] <- rejfreq5
		rejfreq1all[DGP,(3*inobs-2):(3*inobs)] <- rejfreq1

		time_end <- proc.time()
		runtime  <- time_end - time_start

		print("DGP, nobs, nrep")
		print(c(DGP, nobs, nrep))
		print(runtime)
		print("rejfreq10, rejfreq5, rejfreq1")
		print(c(rejfreq10,rejfreq5,rejfreq1))

	} # end of inobs loop

} # end of DGP loop

stopCluster(cl)

rm(Y)
save.image(file = outfilename)

sink()
