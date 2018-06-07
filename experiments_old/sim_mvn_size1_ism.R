rm(list=ls())

library(Rcpp)
library(RcppArmadillo)
library(mvnMix)
library(mixtools)
library(parallel)

outfilename <- "mvn_size1.RData"
sink_name <- "sim_mvn_dgp1.out"

DGPset <- c(1:2)

m <- 1
d <- 2

sigmas_null <- list(diag(d),matrix(c(1,0.5,0.5,1), nrow=2))
param_null <- list()

ninits <- 10
nrep <- 10
nbtsp <- 19
nobsset <- c(200,400)
nnobs <- length(nobsset)

nworker <- min(nrep,detectCores())
cl <- makePSOCKcluster(nworker)

sink(sink_name, append=T)

rejfreq5all <- matrix(0, nrow=16, ncol=6)
rejfreq1all <- matrix(0, nrow=16, ncol=6)
out.all <- vector('list',length=32)


for (DGP in DGPset)
{

sigma <- sigmas_null[[DGP]]

	for (inobs in 1:nnobs)
	{
		time_start <- proc.time()
		nobs <- nobsset[inobs]
		set.seed(123456)
		Y <- array(0, dim=c(nobs,d,nrep))
		for (j in 1:nrep){
			Y[,,j] <- matrix(rnorm(nobs*d),nrow=nobs) %*% chol(sigma)
		}
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
		rejfreq5 <- rejfreq[2,]
		rejfreq1 <- rejfreq[3,]

		rejfreq5all[DGP,(3*inobs-2):(3*inobs)] <- rejfreq5
		rejfreq1all[DGP,(3*inobs-2):(3*inobs)] <- rejfreq1

		time_end <- proc.time()
		runtime  <- time_end - time_start

		print("DGP, nobs, nrep")
		print(c(DGP, nobs, nrep))
		print(runtime)
		print("rejfreq5, rejfreq1")
		print(c(rejfreq5,rejfreq1))

	} # end of inobs loop

} # end of DGP loop

stopCluster(cl)

rm(Y)
save.image(file = outfilename)

sink()
