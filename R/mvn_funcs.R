#' @useDynLib mvnMix
#' @importFrom Rcpp sourceCpp
#' @description Generate initial values used by the PMLE of multivariate normal mixture
#' @export
#' @title mvnmixPMLEinit
#' @name mvnmixPMLEinit
#' @param y n by d matrix of data
#' @param ninits number of initial values to be generated
#' @param m The number of components in the mixture
#' @return A list with the following items:
#' \item{alpha}{m by ninits matrix for alpha}
#' \item{mu}{d*m by ninits matrix for mu, each column is (mu_1',...,mu_m')'}
#' \item{sigma}{d(d+1)/2*m by ninits matrix for sigma, each column is (vech(sigma_1)',...,vech(sigma_m)')'}
mvnmixPMLEinit <- function (y, ninits = 1, m = 2)
{
  # if (normalregMix.test.on) # initial values controlled by normalregMix.test.on
  #   set.seed(normalregMix.test.seed)

  n <- nrow(y)
  d <- ncol(y)
  dsig <- d*(d+1)/2

  alpha <- matrix(runif(m * ninits), nrow=m)
  alpha <- t(t(alpha) / colSums(alpha))
  mu    <- matrix(0, nrow=d, ncol=m*ninits)
  variance    <- matrix(0, nrow=d, ncol=m*ninits)
  sigma <- matrix(0, nrow=dsig, ncol = m*ninits)
  corrmax <- 0.4 # maximum of correlation coefficient in randomly drawn sigma matrix
  # generate initial values for each element of y
  for (i in 1:d){
    y0 <- y[,i]
    mu.i  <- runif(m*ninits, min=min(y0), max =max(y0))
    variance.i <- runif(m*ninits, min=0.1, max=2)*var(y0)
    mu[i,] <- mu.i
    variance[i,] <- variance.i # d by m*ninits matrix
  }

  mu <- matrix(mu, nrow=d*m) # mu was d by m*ninits matrix
  sigma[1,] <- variance[1,]
  sigma[c(2:d),] <- sqrt( t(t(variance[c(2:d),]) * variance[1,]) ) *
    runif((d-1)*m*ninits, min=-corrmax, max=corrmax)
  if (d >=3){
    for (i in 1:(d-2)){
      sigma[i*d-i*(i-1)/2+1,] <- variance[i+1,]
      sigma[c((i*d-i*(i-1)/2+2):((i+1)*d-i*(i+1)/2)),] <-
        sqrt( t(t(variance[c((i+2):d),]) * variance[(i+1),])  ) *
        runif((d-i-1)*m*ninits, min=-corrmax, max=corrmax)
    }
  }
  sigma[dsig,] <- variance[d,]

  sigma <- matrix(sigma, nrow=dsig*m) # sigma was dsig by m*ninits matrix

  list(alpha = alpha, mu = mu, sigma = sigma)

}  # end function mvnmixPMLEinit

#' @description Estimates parameters of a finite mixture of multivariate normals by
#' penalized maximum log-likelhood functions.
#' @export
#' @title mvnmixPMLE
#' @name mvnmixPMLE
#' @param y n by d matrix of data
#' @param m The number of components in the mixture
#' @param ninits The number of randomly drawn initial values.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit The maximum number of iterations.
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param binit The initial value of parameter vector that is included as a candidate parameter vector
#' @return  A list of class \code{mvnmix} with items:
#' \item{coefficients}{A vector of parameter estimates. Ordered as \eqn{\alpha_1,\ldots,\alpha_m,\mu_1,\ldots,\mu_m,\sigma_1,\ldots,\sigma_m}.}
#' \item{parlist}{The parameter estimates as a list containing alpha, mu, and sigma.}
#' \item{vcov}{The estimated variance-covariance matrix.}
#' \item{loglik}{The maximized value of the log-likelihood.}
#' \item{penloglik}{The maximized value of the penalized log-likelihood.}
#' \item{aic}{Akaike Information Criterion of the fitted model.}
#' \item{bic}{Bayesian Information Criterion of the fitted model.}
#' \item{postprobs}{n by m matrix of posterior probabilities for observations}
#' \item{components}{n by 1 vector of integers that indicates the indices of components
#' each observation belongs to based on computed posterior probabilities}
#' \item{call}{The matched call.}
#' \item{m}{The number of components in the mixture.}
#' @note \code{mvnmixPMLE} maximizes the penalized log-likelihood function
#' using the EM algorithm with combining short and long runs of EM steps as in Biernacki et al. (2003).
#' \code{mvnmixPMLE} first runs the EM algorithm from \code{ninits}\eqn{* 4m(1 + p)} initial values
#' with the convertence criterion \code{epsilon.short} and \code{maxit.short}.
#' Then, \code{mvnmixPMLE} uses \code{ninits} best initial values to run the EM algorithm
#' with the convertence criterion \code{epsilon} and \code{maxit}.
#' @references Alexandrovich, G. (2014)
#' A Note on the Article `Inference for Multivariate Normal Mixtures' by J. Chen and X. Tan
#' \emph{Journal of Multivariate Analysis}, \bold{129}, 245--248.
#'
#' Biernacki, C., Celeux, G. and Govaert, G. (2003)
#' Choosing Starting Values for the EM Algorithm for Getting the
#' Highest Likelihood in Multivariate Gaussian Mixture Models,
#' \emph{Computational Statistics and Data Analysis}, \bold{41}, 561--575.
#'
#' Boldea, O. and Magnus, J. R. (2009)
#' Maximum Likelihood Estimation of the Multivariate Normal Mixture Model,
#' \emph{Journal of the American Statistical Association},
#' \bold{104}, 1539--1549.
#'
#' Chen, J. and Tan, X. (2009)
#' Inference for Multivariate Normal Mixtures,
#' \emph{Journal of Multivariate Analysis}, \bold{100}, 1367--1383.
#'
#' McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley \& Sons, Inc.
mvnmixPMLE <- function (y, m = 2, #vcov.method = c("Hessian", "OPG", "none"),
                        ninits = 100, epsilon = 1e-08, maxit = 2000,
                        epsilon.short = 1e-02, maxit.short = 500, binit = NULL) {

  y <- as.matrix(y)
  d <- ncol(y)
  if (d == 1) { stop("y must have more than one columns.") }
  dsig <- d*(d+1)/2
  n <- nrow(y)
  ninits.short <- ninits*10*m*d
  # vcov.method <- match.arg(vcov.method)
  vcov.method <- "none"

  var0   <- var(y) * (n-1)/n

  if (m == 1) {
    mu     <- colMeans(y)
    sigma  <- var0[lower.tri(var0, diag=TRUE)]

    loglik   <- - (n/2) *(d + log(2*pi) + log(det(var0)))
    aic      <- -2*loglik + 2*(m-1 + m*d + m*dsig)
    bic      <- -2*loglik + log(n)*(m-1 + m*d + m*dsig)
    penloglik <- loglik

    parlist <- list(alpha = 1, mu = mu, sigma = sigma)
    coefficients <- c(alpha = 1, mu = mu, sigma = sigma)
    postprobs <- rep(1, n)

  } else {  # m >= 2

    # generate initial values
    tmp <- mvnmixPMLEinit(y = y, ninits = ninits.short, m = m)

    # the following values for (h, k, tau, an) are given by default
    # h       <- 0  # setting h=0 gives PMLE
    # k       <- 0  # k is set to 0 because this is PMLE
    # tau     <- 0.5  # tau is set to 0.5 because this is PMLE
    an      <- 1/sqrt(n)  # penalty term for variance
    var0vec <- var0[lower.tri(var0, diag=TRUE)]
    sigma0  <- rep(var0vec, m)
    mu0     <- double(m+1) # dummy

    # short EM
    b0 <- as.matrix(rbind( tmp$alpha, tmp$mu, tmp$sigma))
    if (!is.null(binit)) {
      b0[ , 1] <- binit
    }
    out.short <- cppMVNmixPMLE(b0, y, mu0, sigma0, m, an, maxit.short,
                               ninits.short, epsilon.short)

    # long EM
    components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
    b1 <- b0[ ,components, drop=FALSE] # b0 has been updated
    out <- cppMVNmixPMLE(b1, y, mu0, sigma0, m, an, maxit, ninits, epsilon)

    index     <- which.max(out$penloglikset)
    alpha <- b1[1:m,index] # b0 has been updated
    mu <- b1[(m+1):(m+m*d),index]
    sigma <- b1[(m+m*d+1):(m+m*d+m*dsig),index]
    penloglik <- out$penloglikset[index]
    loglik    <- out$loglikset[index]
    # postprobs <- matrix(out$post[,index], nrow=n)

    aic     <- -2*loglik + 2*(m-1 + m*d + m*dsig)
    bic     <- -2*loglik + log(n)*(m-1 + m*d + m*dsig)

    mu.matrix <- matrix(mu, nrow=d, ncol=m) # d by m matrix
    sigma.matrix <- matrix(sigma, nrow=dsig, ncol=m) # dsig by m matrix
    mu.order  <- order(mu.matrix[1,])
    alpha     <- alpha[mu.order]
    mu.2        <- mu.matrix[,mu.order]
    sigma.2     <- sigma.matrix[,mu.order]
    mu <- c(mu.2)
    sigma <- c(sigma.2)

    # postprobs <- postprobs[, mu.order]
    # colnames(postprobs) <- c(paste("comp", ".", 1:m, sep = ""))

    parlist <- list(alpha = alpha, mu = mu, sigma = sigma)
    coefficients <- unlist(parlist)

  } # end m >= 2

  if (vcov.method == "none") {
    vcov <- NULL
  } else {
    vcov <- normalmixVcov(y = y, coefficients = coefficients, vcov.method = vcov.method)
  }

  a <- list(coefficients = coefficients, parlist = parlist, loglik = loglik,
            penloglik = penloglik, aic = aic, bic = bic,
            call = match.call(), m = m)
  # a <- list(coefficients = coefficients, parlist = parlist, vcov = vcov, loglik = loglik,
  #           penloglik = penloglik, aic = aic, bic = bic, postprobs = postprobs,
  #           components = getComponentcomponents(postprobs),
  #           call = match.call(), m = m, label = "PMLE")

  # class(a) <- "normalregMix"

  a
}  # end function mvnmixPMLE

#' @description Compute ordinary & penalized log-likelihood ratio resulting from
#' MEM algorithm at k=1,2,3.
#' @title mvnmixMaxPhi
#' @name mvnmixMaxPhi
#' @param y n by d matrix of data
#' @param parlist The parameter estimates as a list containing alpha, mu, and sigma
#' in the form of (alpha = (alpha_1,...,alpha_m), mu = (mu_1',...,mu_m'),
#' sigma = (vech(sigma_1)',...,vech(sigma_m)')
#' @param an a term used for penalty function
#' @param tauset A set of initial tau value candidates
#' @param ninits The number of randomly drawn initial values.
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param maxit The maximum number of iterations.
#' @param verb Determines whether to print a message if an error occurs.
#' @param parallel Determines what percentage of available cores are used, represented by a double in [0,1]. 0.75 is default.
#' @param cl Cluster used for parallelization; if it is \code{NULL}, the system will automatically
#' create a new one for computation accordingly.
#' @return A list with items:
#' \item{loglik}{Log-likelihood resulting from MEM algorithm at k=1,2,3.}
#' \item{penloglik}{Penalized log-likelihood resulting from MEM algorithm at k=1,2,3.}
mvnmixMaxPhi <- function (y, parlist, an, tauset = c(0.1,0.3,0.5),
                          ninits = 10, epsilon.short = 1e-02, epsilon = 1e-08,
                          maxit.short = 500, maxit = 2000,
                          verb = FALSE,
                          parallel = 0.75,
                          cl = NULL) {
  # Given a parameter estimate of an m component model and tuning paramter an,
  # maximize the objective function for computing the modified EM test statistic
  # for testing H_0 of m components against H_1 of m+1 for a univariate normal finite mixture

  warn  <- options(warn=-1) # Turn off warnings

  m <- length(parlist$alpha)
  d <- ncol(y)
  dsig <- d*(d+1)/2

  ninits.short <- ninits*10*m

  loglik.all <- matrix(0,nrow=m*length(tauset),ncol=3)
  penloglik.all <- matrix(0,nrow=m*length(tauset),ncol=3)
  coefficient.all <- matrix(0,nrow=m*length(tauset),ncol=((1+d+dsig)*(m+1)))

  # num.cores <- max(1,floor(detectCores()*parallel))
  # if (num.cores > 1) {
  #   if (is.null(cl))
  #     cl <- makeCluster(detectCores())
  #   registerDoParallel(cl)
  #   results <- foreach (t = 1:length(tauset),
  #                       .export = 'mvnmixMaxPhiStep', .combine = c)  %:%
  #     foreach (h = 1:m) %dopar% {
  #       mvnmixMaxPhiStep (c(h, tauset[t]), y, parlist, an,
  #                         ninits, ninits.short,
  #                         epsilon.short, epsilon,
  #                         maxit.short, maxit,
  #                         verb) }
  #   on.exit(cl)
  #   loglik.all <- t(sapply(results, "[[", "loglik"))
  #   penloglik.all <- t(sapply(results, "[[", "penloglik"))
  #   coefficient.all <- t(sapply(results, "[[", "coefficient"))
  # }
  # else
  for (h in 1:m)
    for (t in 1:length(tauset)) {
      rowindex <- (t-1)*m + h
      tau <- tauset[t]
      result <- mvnmixMaxPhiStep(c(h, tau), y, parlist, an,
                                 ninits, ninits.short,
                                 epsilon.short, epsilon,
                                 maxit.short, maxit,
                                 verb)
      loglik.all[rowindex,] <- result$loglik
      penloglik.all[rowindex,] <- result$penloglik
      coefficient.all[rowindex,] <- result$coefficient
    }

  loglik <- apply(loglik.all, 2, max)  # 3 by 1 vector
  penloglik <- apply(penloglik.all, 2, max)  # 3 by 1 vector
  index <- which.max(loglik.all[ ,3]) # a par (h,m) that gives the highest likelihood at k=3
  coefficient <- as.vector(coefficient.all[index,])

  out <- list(coefficient = coefficient, loglik = loglik, penloglik = penloglik)

  out

}  # end mvnmixMaxPhi

#' @description Given a pair of h and tau and data, compute ordinary &
#' penalized log-likelihood ratio resulting from MEM algorithm at k=1,2,3,
#' tailored for parallelization.
#' @title mvnmixMaxPhiStep
#' @name mvnmixMaxPhiStep
#' @param htaupair A set of h and tau
#' @param y n by d matrix of data
#' @param parlist The parameter estimates as a list containing alpha, mu, and sigma
#' in the form of (alpha = (alpha_1,...,alpha_m), mu = (mu_1',...,mu_m'),
#' sigma = (vech(sigma_1)',...,vech(sigma_m)')
#' @param an a term used for penalty function
#' @param ninits The number of randomly drawn initial values.
#' @param ninits.short The number of candidates used to generate an initial phi, in short MEM
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param maxit The maximum number of iterations.
#' @param verb Determines whether to print a message if an error occurs.
#' @return A list of phi, log-likelihood, and penalized log-likelihood resulting from MEM algorithm.
mvnmixMaxPhiStep <- function (htaupair, y, parlist, an,
                              ninits, ninits.short,
                              epsilon.short, epsilon,
                              maxit.short, maxit,
                              verb)
{
  alpha0 <- parlist$alpha

  m      <- length(alpha0)
  m1     <- m+1
  k      <- 1
  n      <- nrow(y)
  d      <- ncol(y)
  dsig   <- d*(d+1)/2
  h      <- as.numeric(htaupair[1])
  tau    <- as.numeric(htaupair[2])

  mu0    <- parlist$mu  # d*m by 1
  mu0matrix <- matrix(mu0, nrow=d, ncol=m)
  mu01   <- mu0matrix[1,]
  mu0h   <- c(-1e+10,mu01,1e+10)        # m+2 by 1
  sigma0 <- parlist$sigma # dsig*m by 1
  sigma0h<- c(sigma0[1:(h*dsig)],sigma0[((h-1)*dsig+1):(m*dsig)]) # (m+1)*dsig by 1

  # generate initial values
  tmp <- mvnmixPhiInit(y = y, parlist = parlist, h=h, tau = tau, ninits = ninits.short)

  # short EM
  b0 <- as.matrix(rbind(tmp$alpha, tmp$mu, tmp$sigma))
  out.short <- cppMVNmixPMLE(b0, y, mu0h, sigma0h, m1, an, maxit.short, ninits.short,
                             epsilon.short, tau, h, k)
  components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
  if (verb && any(out.short$notcg)) {
    cat(sprintf("non-convergence rate at short-EM = %.3f\n",mean(out.short$notcg)))
  }
  # long EM
  b1 <- as.matrix(b0[ ,components, drop=FALSE])
  out <- cppMVNmixPMLE(b1, y, mu0h, sigma0h, m1, an, maxit, ninits, epsilon, tau, h, k)

  index     <- which.max(out$penloglikset)
  alpha <- b1[1:m1,index]
  mu <- b1[(m1+1):(m1+d*m1),index]
  sigma <- b1[(m1+d*m1+1):(m1+d*m1+dsig*m1),index]

  mu.matrix <- matrix(mu, nrow=d, ncol=m1) # d by m1 matrix
  sigma.matrix <- matrix(sigma, nrow=dsig, ncol=m1) # dsig by m1 matrix
  sigma0h.matrix <- matrix(sigma0h, nrow=dsig, ncol=m1) # dsig by m1 matrix
  mu.order  <- order(mu.matrix[1,])
  alpha     <- alpha[mu.order]
  mu.2        <- mu.matrix[,mu.order]
  sigma.2     <- sigma.matrix[,mu.order]
  sigma0h.2     <- sigma0h.matrix[,mu.order]
  mu <- c(mu.2)
  sigma <- c(sigma.2)
  sigma0h <- c(sigma0h.2)
  b <- as.matrix( c(alpha, mu, sigma) )

  # initilization
  loglik <-  vector("double", 3)
  penloglik <-  vector("double", 3)
  coefficient <- vector("double", length(b))

  penloglik[1] <- out$penloglikset[[index]]
  loglik[1]    <- out$loglikset[[index]]
  for (k in 2:3) {
    ninits <- 1
    maxit <- 2
    # Two EM steps
    out <- cppMVNmixPMLE(b, y, mu0h, sigma0h, m1, an, maxit, ninits, epsilon, tau, h, k)
    alpha <- b[1:m1,1] # b has been updated
    mu <- b[(m1+1):(m1+d*m1),1]
    sigma <- b[(m1+d*m1+1):(m1+d*m1+dsig*m1),1]
    loglik[k]    <- out$loglikset[[1]]
    penloglik[k]   <- out$penloglikset[[1]]

    # Check singularity: if singular, break from the loop
    # compute determinants
    detsigma <- double(m1)
    for (j in 1:m1){
      sigma.jmat <- diag(d)
      sigma.j <- sigma[((j-1)*dsig+1):(j*dsig)]
      sigma.jmat[lower.tri(sigma.jmat, diag=TRUE)] <- sigma.j
      sigma.jmat <- t(sigma.jmat) + sigma.jmat
      diag(sigma.jmat) <- diag(sigma.jmat)/2
      detsigma[j] <- det(sigma.jmat)
    }

    if ( any(detsigma < 1e-06) || any(alpha < 1e-06) || is.na(sum(alpha)) ) {
      loglik[k]    <- -Inf
      penloglik[k]   <- -Inf
      break
    }

    mu.matrix <- matrix(mu, nrow=d, ncol=m1) # d by m1 matrix
    sigma.matrix <- matrix(sigma, nrow=dsig, ncol=m1) # dsig by m1 matrix
    sigma0h.matrix <- matrix(sigma0h, nrow=dsig, ncol=m1) # dsig by m1 matrix
    mu.order  <- order(mu.matrix[1,])
    alpha     <- alpha[mu.order]
    mu.2        <- mu.matrix[,mu.order]
    sigma.2     <- sigma.matrix[,mu.order]
    sigma0h.2     <- sigma0h.matrix[,mu.order]
    mu <- c(mu.2)
    sigma <- c(sigma.2)
    sigma0h <- c(sigma0h.2)
  }
  coefficient <- as.matrix( c(alpha, mu, sigma) ) # at k=3

  return (list(coefficient = coefficient, loglik = loglik, penloglik = penloglik))
}

#' @description Generates lists of parameters for initial candidates used by
#' the modified EM test for mixture of multivariate normals.
#' @title mvnmixPhiInit
#' @name mvnmixPhiInit
#' @param y n by d matrix of data
#' @param parlist The parameter estimates as a list containing alpha, mu, and sigma
#' in the form of (alpha = (alpha_1,...,alpha_m),
#' mu = (mu_1',...,mu_m'), sigma = (vech(sigma_1)',...,vech(sigma_m)')
#' @param h h used as index for pivoting
#' @param tau Tau used to split the h-th component
#' @param ninits number of initial values to be generated
#' @return A list with the following items:
#' \item{alpha}{m+1 by ninits matrix for alpha}
#' \item{mu}{d*(m+1) by ninits matrix for mu}
#' \item{sigma}{d*(d+1)/2*(m+1) by ninits matrix for sigma}
mvnmixPhiInit <- function (y, parlist, h, tau, ninits = 1)
{
  # if (normalregMix.test.on) # initial values controlled by normalregMix.test.on
  #   set.seed(normalregMix.test.seed)

  n       <- nrow(y)
  d       <- ncol(y)
  dsig    <- d*(d+1)/2
  mu0     <- parlist$mu
  sigma0  <- parlist$sigma
  alpha0  <- parlist$alpha
  m       <- length(alpha0)

  mu0matrix <- matrix(mu0, nrow=d, ncol=m)
  sigma0matrix <- matrix(sigma0, nrow=dsig, ncol=m)
  y1 <- y[,1]
  mu01 <- mu0matrix[1,]

  if (m>=2){
    mid <- (mu01[1:(m-1)]+mu01[2:m])/2  # m-1 by 1
    lb0 <- c(min(y1),mid)               # m by 1
    lb  <- c(lb0[1:h],lb0[h:m])         # m+1 by 1
    ub0 <- c(mid,max(y1))               # m by 1
    ub  <- c(ub0[1:h],ub0[h:m])         # m+1 by 1
  } else {
    lb  <- c(min(y1),min(y1))
    ub  <- c(max(y1),max(y1))
  }

  mu    <- matrix(0, nrow=d, ncol=(m+1)*ninits)
  sigma <- matrix(0, nrow=dsig, ncol=(m+1)*ninits)
  variance <- matrix(0, nrow=d, ncol=(m+1)*ninits)

  ninits1 <- floor(ninits/2)
  ninits2 <- ninits - ninits1

  mu[1,] <- runif((m+1)*ninits, min=lb, max=ub)
  sigma01 <- sigma0matrix[1,] # m by 1
  sigma.1.hyp <- c(sigma01[1:h],sigma01[h:m])  # m+1 by 1
  sigma.1 <- runif((m+1)*ninits,min=sigma.1.hyp*0.25,max=sigma.1.hyp*2)
  variance[1,] <- sigma.1

  for (i in 2:d){
    y.i <- y[,i]
    mu0i <- mu0matrix[i,] # m by 1
    sigma0i <- sigma0matrix[(i-1)*d-(i-1)*(i-2)/2+1,] # m by 1
    mu.i.hyp <- c(mu0i[1:h],mu0i[h:m])  # m+1 by 1
    sigma.i.hyp <- c(sigma0i[1:h],sigma0i[h:m])  # m+1 by 1
    sd.i.hyp <- sqrt(sigma.i.hyp)
    mu.i1  <- runif((m+1)*ninits1, min=mu.i.hyp-sd.i.hyp, max=mu.i.hyp+sd.i.hyp)
    mu.i2  <- runif((m+1)*ninits2, min=min(y.i), max=max(y.i))
    mu[i,] <- c(mu.i1, mu.i2)
    sigma.i <- runif((m+1)*ninits, min=sigma.i.hyp*0.25,max=sigma.i.hyp*2)
    variance[i,] <- sigma.i
  }
  mu <- matrix(mu, nrow=d*(m+1)) # mu was d by (m+1)*ninits matrix

  corrmax <- 0.4 # maximum of correlation coefficient in randomly drawn sigma matrix

  sigma[1,] <- variance[1,]
  sigma[c(2:d),] <- sqrt( t(t(variance[c(2:d),]) * variance[1,]) ) *
    runif((d-1)*(m+1)*ninits, min=-corrmax, max=corrmax)
  if (d >=3){
    for (i in 1:(d-2)){
      sigma[i*d-i*(i-1)/2+1,] <- variance[i+1,]
      sigma[c((i*d-i*(i-1)/2+2):((i+1)*d-i*(i+1)/2)),] <-
        sqrt( t(t(variance[c((i+2):d),]) * variance[(i+1),])  ) *
        runif((d-i-1)*(m+1)*ninits, min=-corrmax, max=corrmax)
    }
  }
  sigma[d*(d+1)/2,] <- variance[d,]

  sigma <- matrix(sigma, nrow=dsig*(m+1)) # sigma was dsig by (m+1)*ninits matrix

  alpha.hyp <- c(alpha0[1:h],alpha0[h:m])  # m+1 by 1
  alpha.hyp[h:(h+1)] <- c(alpha.hyp[h]*tau,alpha.hyp[h+1]*(1-tau))
  alpha <- matrix(rep.int(alpha.hyp,ninits),nrow=m+1)

  list(alpha = alpha, mu = mu, sigma = sigma)

}  # end function mvnmixPhiInit
