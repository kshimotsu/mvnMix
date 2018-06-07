rMVNmixPMLE <- function(b0, y, mu0, sigma0, m, p, an, maxit, ninits, tol, tau=0.5, h=0, k=0)
{
  # double oldpenloglik, s0j, diff, minr, w_j, sum_l_j, ssr_j, alphah, tauhat;
  # double ll = 0; // force initilization
  # double penloglik = 0; // force initialization
  
  ub <- rep(0,m)
  lb <- rep(0,m)
  # /* Lower and upper bound for mu */
  if (k==1) {  # If k==1, compute upper and lower bounds
    mu0[0] = -Inf;
    mu0[m] = Inf;
    for (j in (1:h)) {
      lb[j] <- (mu0[j]+mu0[j+1])/2.0
      ub[j] <- (mu0[j+1]+mu0[j+2])/2.0
    }
    for (j in ((h+1):m)) {
      lb(j) <- (mu0[j-1]+mu0[j])/2.0
      ub(j) <- (mu0[j]+mu0[j+1])/2.0
    }
  }
  
  penloglikset <- double(ninits)
  loglikset <- double(ninits)
  
  n <- nrow(y)
  d <- ncol(y)
  dsig <- d*(d+1)/2
  notcg <- integer(ninits)
  itersum <- double(ninits)
  l0 <- matrix(0, nrow=n, ncol=d)
  l0.acc <- matrix(0, nrow=n, ncol=d)
  pen <- double(m)
  
  nparam <- m + m*d + m*dsig # number of parameters in alpha, mu, sigma
  
  # /* iteration over ninits initial values of b */
  for (jn in (1:ninits)) {
    
    # /* initialize EM iteration */
    b_jn <- b0[,jn]
    # param.old <- b_jn
    alpha <- b_jn[1:m] 
    mu <- b_jn[(m+1):(m+m*d)]
    sigma <- b_jn[(m+m*d+1):(m+m*d+m*dsig)]
    sigma.mat <- matrix(0, nrow=d, ncol=d*m)  # matrix of sigmas
    sigma0.mat <- matrix(0, nrow=d, ncol=d*m)  # matrix of sigma0s
    for (j in (1:m)) {
      sigma.mat[,((j-1)*d+1):(j*d)] <- sigma.vec.to.mat(sigma[((j-1)*dsig+1):(j*dsig)],d)
      sigma0.mat[,((j-1)*d+1):(j*d)] <- 
        sigma.vec.to.mat(sigma0[((j-1)*dsig+1):(j*dsig)],d)
    }
    oldpenloglik <- -Inf
    diff <- 1.0
    sing <- 0
    det.sigma <- double(m)
    
    # /* EM loop begins */
    for (iter in (1:maxit)) {
      # /* Compute the penalized loglik. Note that penalized loglik uses old (not updated) sigma */
      # See Chan and Tan p. 1370 for the definition of the penalty term
      for (j in (1:m)){
        muj <- mu[((j-1)*d+1):(j*d)]
        sigmaj.mat <- sigma.mat[,((j-1)*d+1):(j*d)]
        l0[,j] <- dmvnorm(y,muj,sigmaj.mat)
        sigma0j.mat <-sigma0.mat[,((j-1)*d+1):(j*d)]
        s0j <- solve(sigmaj.mat,sigma0j.mat)
        pen[j] <- an*(sum(diag(s0j)) - 2*log(det(s0j)) -d)
      }
      ll <- sum(log(l0 %*% alpha))
      penloglik <- ll + log(2.0) + min(log(tau),log(1-tau)) - sum(pen)
      
      diff <- penloglik - oldpenloglik
      oldpenloglik <- penloglik
      
      # /* Exit from the loop if diff is NaN or NA */
      if (is.nan(diff) || is.na(diff)) {
        penloglik <- -Inf
        ll <- -Inf
        notcg[jn] <- 1
        itersum[jn] <- iter
        break
      }
      
      # /* Normal exit */
      if (diff < tol ){
        itersum[jn] <- iter
        break
      }
      
      # # store previous parameter estimtes for Anderson acceleration
      # 
      # param.older <- param.old
      # alpha.old <- alpha
      # mu.old <- mu
      # sigma.old <- sigma
      # param.old <- c(alpha.old,mu.old,sigma.old)
      # f.acc.old <- f.acc
      #       
      # /* update alpha, mu, and sigma */
      w <- t(t(l0)*alpha)/c(l0 %*% alpha)
      alpha <- colMeans(w)
      for (j in 1:m){
        muj <- colSums( y* w[,j] ) / (n*alpha[j])
        ydot <- t(t(y) - muj)
        ssrj <- t(ydot* w[,j]) %*% ydot
        sigmaj.mat <- (2*an*sigma0j.mat + ssrj)/ (2*an + n*alpha[j])
        # /* If k ==1, impose lower and upper bound */
        if (k==1) {
          muj[1] <- min( max(muj[1],lb[j]), ub[j])
        }
        
        mu[((j-1)*d+1):(j*d)] <- muj
        det.sigma[j] <- det(sigmaj.mat)
        sigma.mat[,((j-1)*d+1):(j*d)] <- sigmaj.mat
        sigma[((j-1)*dsig+1):(j*dsig)] <- sigmaj.mat[lower.tri(sigmaj.mat, diag=TRUE)]
      }
      
      # /* for PMLE, we set k=0 (default value) */
      #   /* for EM test, we start from k=1       */
      #   /*   if k==1, we don't update tau       */
      # /*   if k>1, we update tau              */
      # if (k==1){
      # alphah = (alpha(h-1)+alpha(h));
      # alpha(h-1) = alphah*tau;
      # alpha(h) = alphah*(1-tau);
      # } else if (k>1) {
      # alphah = (alpha(h-1)+alpha(h));
      # tauhat = alpha(h-1)/(alpha(h-1)+alpha(h));
      # if(tauhat <= 0.5) {
      # tau = fmin((alpha(h-1)*n + 1.0)/(alpha(h-1)*n + alpha(h)*n + 1.0), 0.5);
      # } else {
      # tau = fmax(alpha(h-1)*n /(alpha(h-1)*n + alpha(h)*n + 1.0), 0.5);
      # }
      # alpha(h-1) = alphah*tau;
      # alpha(h) = alphah*(1-tau);
      # }
      
      # param.new <- c(alpha,mu,sigma)
      # f.acc <- param.new - param.old # f in Anderson acceleration
      # fmatrix.old <- fmatrix
      # fmatrix[,c(1:(m.acc-1))] <- fmatrix[,c(2:m.acc)]
      # fmatrix[,m.acc] <- f.acc # collects parameter updates
      # Fmatrix <- fmatrix - fmatrix.old
      # Xmatrix[,c(1:(m.acc-1))] <- Xmatrix[,c(2:m.acc)]
      # Xmatrix[,m.acc] <- param.old - param.older
      # if (iter >=4){
      #   Gamma <- lsfit(x=Fmatrix, y=f.acc, intercept=FALSE)
      #   param.new.acc <- param.new - (Xmatrix + Fmatrix) %*% Gamma$coefficients
      #   alpha.acc <-  param.new.acc[1:m]
      #   mu.acc <-  param.new.acc[(m+1):(m+m*d)]
      #   sigma.acc <- param.new.acc[(m+m*d+1):(m+m*d+m*dsig)]
      #   # check if the penloglik increases
      #   for (j in (1:m)){
      #     muj <- mu[((j-1)*d+1):(j*d)]
      #     sigmaj.mat <- sigma.vec.to.mat(sigma[((j-1)*dsig+1):(j*dsig)],d)
      #     l0[,j] <- dmvnorm(y,muj,sigmaj.mat)
      #     sigma0j.mat <-sigma0.mat[,((j-1)*d+1):(j*d)]
      #     s0j <- solve(sigmaj.mat,sigma0j.mat)
      #     pen[j] <- -an*(sum(diag(s0j)) - 2*log(det(s0j)) -d)
      #   }
      #   ll <- sum(log(l0 %*% alpha))
      #   penloglik.new <- ll + log(2.0) + min(log(tau),log(1-tau)) - sum(pen)
      #   
      #   # check if the penloglik increases
      #   for (j in (1:m)){
      #     muj.acc <- mu.acc[((j-1)*d+1):(j*d)]
      #     sigmaj.mat.acc <- sigma.vec.to.mat(sigma.acc[((j-1)*dsig+1):(j*dsig)],d)
      #     l0.acc[,j] <- dmvnorm(y,muj.acc,sigmaj.mat.acc)
      #     sigma0j.mat <-sigma0.mat[,((j-1)*d+1):(j*d)]
      #     s0j.acc <- solve(sigmaj.mat.acc,sigma0j.mat)
      #     pen[j] <- -an*(sum(diag(s0j.acc)) - 2*log(det(s0j.acc)) -d)
      #   }
      #   ll.acc <- sum(log(l0.acc %*% alpha.acc))
      #   penloglik.acc <- ll.acc + log(2.0) + min(log(tau),log(1-tau)) - sum(pen)
      #   
      #   if (!is.nan(penloglik.acc) && !is.na(penloglik.acc) && (penloglik.acc > penloglik.new)){
      #     alpha <- alpha.acc
      #     mu <- mu.acc
      #     sigma <- sigma.acc
      #   }
      # }
      
      
      
      # /* Check singularity */
      if (min(alpha) < 1e-8 || any(is.nan(alpha)) || min(det.sigma) < 1e-8){
        sing <- 1
      }
      
      # /* Exit from the loop if singular */
      if (sing) {
        notcg[jn] <- 1
        itersum[jn] <- iter
        break
      }
      
    } #/* EM loop ends */
    
    penloglikset[jn] <- penloglik
    loglikset[jn] <- ll
    # update b0
    b0[1:m,jn] <- alpha
    b0[(m+1):(m+m*d),jn] <- mu
    b0[(m+m*d+1):(m+m*d+m*dsig),jn] <- sigma
    
  } # /* end for (jn=0; jn<ninits; jn++) loop */
  
  list(penloglikset = penloglikset, loglikset = loglikset, b0=b0, itersum = itersum, notcg = notcg)
  
}  # end function rMVNmixPMLE

# OLD version using R version of cppMVNmixPMLE
#' @description Estimates parameters of a finite mixture of multivariate normals by 
#' penalized maximum log-likelhood functions.
#' @export
#' @title mvnmixPMLE_R
#' @name mvnmixPMLE_R
#' @param y n by d matrix of data
#' @param m The number of components in the mixture
#' @param ninits The number of randomly drawn initial values.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit The maximum number of iterations.
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param binit The initial value of parameter vector that is included as a candidate parameter vector
#' @return  A list of class \code{normalMix} with items:
#' \item{coefficients}{A vector of parameter estimates. Ordered as \eqn{\alpha_1,\ldots,\alpha_m,\mu_1,\ldots,\mu_m,\sigma_1,\ldots,\sigma_m,\gamma}.}
#' \item{parlist}{The parameter estimates as a list containing alpha, mu, and sigma (and gam if z is included in the model).}
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
mvnmixPMLE_R <- function (y, m = 2, #vcov.method = c("Hessian", "OPG", "none"),
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
    sigma  <- var0
    
    loglik   <- - (n/2) *(d + log(2*pi) + log(det(sigma)))
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
    out.short <- rMVNmixPMLE(b0, y, mu0, sigma0, m, p, an, maxit.short,
                             ninits.short, epsilon.short)
    b01 <- b0[ ,1, drop=FALSE] # First column of b0
    out.0 <- rMVNmixPMLE(b01, y, mu0, sigma0, m, p, an, maxit.short,1, epsilon.short)
    # out.short <- cppMVNmixPMLE(b0, y, mu0, sigma0, m, p, an, maxit.short,
    #                            ninits.short, epsilon.short)
    # out.short <- cppNormalmixPMLE(b0, y, ztilde, mu0, sigma0, m, p, an, maxit.short,
    #                               ninits.short, epsilon.short)
    # long EM
    components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
    b1 <- out.short$b0[ ,components, drop=FALSE] # b0 has been updated
    out <- rMVNmixPMLE(b1, y, mu0, sigma0, m, p, an, maxit, ninits, epsilon)
    # out <- cppMVNmixPMLE(b1, y, mu0, sigma0, m, p, an, maxit, ninits, epsilon)
    
    index     <- which.max(out$penloglikset)
    alpha <- b1[1:m,index] # b0 has been updated
    mu <- b1[(m+1):(m+m*d),index]
    sigma <- b1[(m+m*d+1):(m+m*d+m*dsig),index]
    penloglik <- out$penloglikset[index]
    loglik    <- out$loglikset[index]
    # postprobs <- matrix(out$post[,index], nrow=n)
    
    aic     <- -2*loglik + 2*(m-1 + m*d + m*dsig)
    bic     <- -2*loglik + log(n)*(m-1 + m*d + m*dsig)
    
    mu.matrix <- t(matrix(mu, nrow=d, ncol=m)) # m by d matrix
    sigma.matrix <- t(matrix(sigma, nrow=dsig, ncol=m)) # m by dsig matrix
    mu.order  <- order(mu.matrix[,1])
    alpha     <- alpha[mu.order]
    mu.2        <- mu.matrix[mu.order,]
    sigma.2     <- sigma.matrix[mu.order,]
    mu <- c(t(mu.2))
    sigma <- c(t(sigma.2))
    
    # postprobs <- postprobs[, mu.order]
    # colnames(postprobs) <- c(paste("comp", ".", 1:m, sep = ""))
    
    parlist <- list(alpha = alpha, mu = mu, sigma = sigma)
    coefficients <- unlist(parlist)
    
  } # end m >= 2
  
  if (vcov.method == "none") {
    vcov <- NULL
  } else {
    vcov <- normalmixVcov(y = y, coefficients = coefficients, z = z , vcov.method = vcov.method)
  }
  
  a <- list(coefficients = coefficients, parlist = parlist,loglik = loglik,
            penloglik = penloglik, aic = aic, bic = bic, 
            call = match.call(), m = m)
  # a <- list(coefficients = coefficients, parlist = parlist, vcov = vcov, loglik = loglik,
  #           penloglik = penloglik, aic = aic, bic = bic, postprobs = postprobs,
  #           components = getComponentcomponents(postprobs),
  #           call = match.call(), m = m, label = "PMLE")
  
  # class(a) <- "normalregMix"
  
  a
}  # end function nvnmixPMLE_R


#' @description Estimates parameters of a finite mixture of multivariate normals by 
#' penalized maximum log-likelhood functions.
#' @export
#' @title mvnmixPMLE.old
#' @name mvnmixPMLE.old
#' @param y n by d matrix of data
#' @param m The number of components in the mixture
#' @param ninits The number of randomly drawn initial values.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit The maximum number of iterations.
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param binit The initial value of parameter vector that is included as a candidate parameter vector
#' @return  A list of class \code{normalMix} with items:
#' \item{coefficients}{A vector of parameter estimates. Ordered as \eqn{\alpha_1,\ldots,\alpha_m,\mu_1,\ldots,\mu_m,\sigma_1,\ldots,\sigma_m,\gamma}.}
#' \item{parlist}{The parameter estimates as a list containing alpha, mu, and sigma (and gam if z is included in the model).}
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
mvnmixPMLE.old <- function (y, m = 2, #vcov.method = c("Hessian", "OPG", "none"),
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
    sigma  <- var0
    
    loglik   <- - (n/2) *(d + log(2*pi) + log(det(sigma)))
    aic      <- -2*loglik + 2*(m-1 + m*d + m*dsig)
    bic      <- -2*loglik + log(n)*(m-1 + m*d + m*dsig)
    penloglik <- loglik
    
    parlist <- list(alpha = 1, mu = mu, sigma = sigma)
    coefficients <- c(alpha = 1, mu = mu, sigma = sigma)
    postprobs <- rep(1, n)
    
  } else {  # m >= 2
    
    # generate initial values
    tmp <- mvnmixPMLEinit.old(y = y, ninits = ninits.short, m = m)
    
    # the following values for (h, k, tau, an) are given by default
    # h       <- 0  # setting h=0 gives PMLE
    # k       <- 0  # k is set to 0 because this is PMLE
    # tau     <- 0.5  # tau is set to 0.5 because this is PMLE
    an      <- 1/sqrt(n)  # penalty term for variance
    var0vec <- var0[lower.tri(var0, diag=TRUE)]
    sigma0  <- rep(var0vec, m)
    mu0     <- double(m+1) # dummy
    
    # if (is.null(z)) {
    #   ztilde <- matrix(0) # dummy
    # } else {
    #   ztilde <- z
    # }
    # short EM
    b0 <- as.matrix(rbind( tmp$alpha, tmp$mu, tmp$sigma))
    if (!is.null(binit)) {
      b0[ , 1] <- binit
    }
    out.short <- cppMVNmixPMLE(b0, y, mu0, sigma0, m, p, an, maxit.short,
                               ninits.short, epsilon.short)
    # out.short <- cppNormalmixPMLE(b0, y, ztilde, mu0, sigma0, m, p, an, maxit.short,
    #                               ninits.short, epsilon.short)
    # long EM
    components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
    b1 <- out.short$b0[ ,components] # b0 has been updated
    out <- cppMVNmixPMLE(b1, y, mu0, sigma0, m, p, an, maxit, ninits, epsilon)
    
    index     <- which.max(out$penloglikset)
    alpha <- b1[1:m,index] # b0 has been updated
    mu <- b1[(m+1):(m+m*d),index]
    sigma <- b1[(m+m*d+1):(m+m*d+m*dsig),index]
    # if (!is.null(z)) {
    #   gam     <- b1[(3*m+1):(3*m+p),index]
    # }
    penloglik <- out$penloglikset[index]
    loglik    <- out$loglikset[index]
    # postprobs <- matrix(out$post[,index], nrow=n)
    
    aic     <- -2*loglik + 2*(m-1 + m*d + m*dsig)
    bic     <- -2*loglik + log(n)*(m-1 + m*d + m*dsig)
    
    mu.matrix <- t(matrix(mu, nrow=d, ncol=m)) # m by d matrix
    sigma.matrix <- t(matrix(sigma, nrow=dsig, ncol=m)) # m by dsig matrix
    mu.order  <- order(mu.matrix[,1])
    alpha     <- alpha[mu.order]
    mu.2        <- mu.matrix[mu.order,]
    sigma.2     <- sigma.matrix[mu.order,]
    mu <- c(t(mu.2))
    sigma <- c(t(sigma.2))
    
    # postprobs <- postprobs[, mu.order]
    # colnames(postprobs) <- c(paste("comp", ".", 1:m, sep = ""))
    
    parlist <- list(alpha = alpha, mu = mu, sigma = sigma)
    coefficients <- unlist(parlist)
    
  } # end m >= 2
  
  if (vcov.method == "none") {
    vcov <- NULL
  } else {
    vcov <- normalmixVcov(y = y, coefficients = coefficients, z = z , vcov.method = vcov.method)
  }
  
  a <- list(coefficients = coefficients, parlist = parlist,loglik = loglik,
            penloglik = penloglik, aic = aic, bic = bic, 
            call = match.call(), m = m)
  # a <- list(coefficients = coefficients, parlist = parlist, vcov = vcov, loglik = loglik,
  #           penloglik = penloglik, aic = aic, bic = bic, postprobs = postprobs,
  #           components = getComponentcomponents(postprobs),
  #           call = match.call(), m = m, label = "PMLE")
  
  # class(a) <- "normalregMix"
  
  a
}  # end function nvnmixPMLE2


#' @description Generate initial values used by the PMLE of multivariate normal mixture
#' @export
#' @title mvnmixPMLEinit.old
#' @name mvnmixPMLEinit.old
#' @param y n by d matrix of data
#' @param ninits number of initial values to be generated
#' @param m The number of components in the mixture
#' @return A list with the following items:
#' \item{alpha}{m by ninits matrix for alpha}
#' \item{mu}{d*m by ninits matrix for mu, each column is (mu_1',\ldots,mu_m')'}
#' \item{sigma}{d(d+1)/2*m by ninits matrix for sigma, each column is (vech(sigma_1)',\ldots,vech(sigma_m)')'}
mvnmixPMLEinit.old <- function (y, ninits = 1, m = 2)
{
  # if (normalregMix.test.on) # initial values controlled by normalregMix.test.on
  #   set.seed(normalregMix.test.seed)
  
  n <- nrow(y)
  d <- ncol(y)
  gam <- NULL
  
  alpha <- matrix(runif(m * ninits), nrow=m)
  alpha <- t(t(alpha) / colSums(alpha))
  mu    <- matrix(0, nrow=d, ncol=m*ninits)
  variance    <- matrix(0, nrow=d, ncol=m*ninits)
  sigma <- matrix(0,nrow=d*(d+1)/2, ncol = m*ninits)
  corrmax <- 0.4 # maximum of correlation coefficient in randomly drawn sigma matrix
  # generate initial values for each element of y
  for (i in 1:d){
    y0 <- y[,i]
    out <- normalmixPMLE(y0, m = m, vcov.method = "none")  # fit a univariate mixture to each column of y
    param.matrix <- matrix(out$coefficients,m,3)
    mu0 <- param.matrix[,2]
    sigma0 <- param.matrix[,3]
    msigma0 <- max(sigma0)
    mu.i  <- matrix(runif(m*ninits, min=-1, max=1), nrow=m) * 2*msigma0 + mu0
    variance.i  <- matrix(runif(m*ninits, min=1/2, max=2), nrow=m) * sigma0^2
    ind <- matrix(c(1:m*ninits), nrow=m)  # m x ninits matrix
    ind.s <- apply(ind,2,sample) # shuffule each column
    mu.i <- mu.i[ind.s]         # m*ninits vector
    variance.i <- variance.i[ind.s]   # m*ninits vector
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
  sigma[d*(d+1)/2,] <- variance[d,]
  
  sigma <- matrix(sigma, nrow=d*(d+1)/2*m) # sigma was d*(d+1)/2 by m*ninits matrix
  
  list(alpha = alpha, mu = mu, sigma = sigma)
  
}  # end function mvnmixPMLEinit.old
