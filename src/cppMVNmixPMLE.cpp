#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// const double SINGULAR_EPS = 10e-10; // criteria for matrix singularity

// [[Rcpp::export]]
List cppMVNmixPMLE(NumericMatrix bs,
                   NumericMatrix ys,
                   NumericVector mu0s,
                   NumericVector sigma0s,
                   int m,
                   double an,
                   int maxit = 2000,
                   int ninits = 10,
                   double tol = 1e-8,
                   double tau = 0.5,
                   int h = 0,
                   int k = 0) {
  int n = ys.nrow();
  int d = ys.ncol();
  int dsig = d*(d+1)/2;
  arma::mat b(bs.begin(), bs.nrow(), bs.ncol(), false);
  arma::mat y(ys.begin(), ys.nrow(), ys.ncol(), false);
  arma::mat sigmamat_j = arma::zeros(d, d);
  arma::mat sigma_j_inv = arma::zeros(d, d);
  arma::mat sigma0mat_j = arma::zeros(d, d);
  arma::mat s0j(d, d), ssr_j(d, d);
  arma::mat sigmamat = arma::zeros(d, m*d);
  arma::mat sigma0mat = arma::zeros(d, m*d);
  arma::vec mu0(mu0s.begin(), mu0s.size(), false);
  arma::vec sigma0(sigma0s.begin(), sigma0s.size(), false);
  arma::vec b_jn(bs.nrow());
  arma::vec lb(m),ub(m);
  arma::vec alpha(m), mu(m*d), sigma(m*dsig), alp_sig(m);
  arma::vec detsigma(m), pen(m);
  arma::mat l_m(n,m);
  arma::mat r(n,m);
  arma::vec minr(n), sum_l_m(n);
  arma::mat w(n,m);
  arma::mat post(m*n,ninits);
  arma::vec notcg = arma::zeros(ninits);
  arma::vec penloglikset(ninits), loglikset(ninits);
  arma::mat ydot(n,d);
  arma::mat rtilde(n,d);
  arma::vec mu_j(d);
  arma::vec wtilde(n);
  arma::mat ytilde(d,n);
  int sing;
  double oldpenloglik, w_j, diff, alphah, tauhat;
  // double  ;
  double ll = 0; // force initilization
  double penloglik = 0; // force initialization

  /* Lower and upper bound for the first elemetn of mu */
  if (k==1) {  // If k==1, compute upper and lower bounds
    mu0(0) = R_NegInf;
    mu0(m) = R_PosInf;
    for (int j=0; j<h; j++) {
      lb(j) = (mu0(j)+mu0(j+1))/2.0;
      ub(j) = (mu0(j+1)+mu0(j+2))/2.0;
    }
    for (int j=h; j<m; j++) {
      lb(j) = (mu0(j-1)+mu0(j))/2.0;
      ub(j) = (mu0(j)+mu0(j+1))/2.0;
    }
  }

  /* iteration over ninits initial values of b */
  for (int jn=0; jn<ninits; jn++) {

    /* initialize EM iteration */
    b_jn = b.col(jn);
    alpha = b_jn.subvec(0,m-1);
    mu = b_jn.subvec(m,m+m*d-1);
    sigma = b_jn.subvec(m+m*d,m+m*d+m*dsig-1);
    int dum=0;
    for (int j=0; j < m; ++j){
      for (int ii=0; ii<d; ++ii){
        for (int jj=ii; jj<d; ++jj){
          sigmamat_j(jj,ii) = sigma(dum);
          sigma0mat_j(jj,ii) = sigma0(dum);
          dum++;
        }
      }
      sigmamat_j = symmatl(sigmamat_j);
      sigma0mat_j = symmatl(sigma0mat_j);
      sigmamat.cols(j*d,(j+1)*d-1) = sigmamat_j;
      sigma0mat.cols(j*d,(j+1)*d-1) = sigma0mat_j;
    }
    oldpenloglik = R_NegInf;
    diff = 1.0;
    sing = 0;

// sigmamat.print();
// mu.print();

    /* EM loop begins */
    for (int iter = 0; iter < maxit; iter++) {
      /* standardized squared residual */
      for (int j=0; j < m; j++) {
        mu_j = mu.subvec(j*d,(j+1)*d-1);
        ydot = y.each_row() - mu_j.t();
        detsigma(j) = det(sigmamat.cols(j*d,(j+1)*d-1));
        // if ( detsigma(j) < 1e-8 || isnan(detsigma(j)) ) {
        //     sigma_j_inv = arma::eye(d,d);
        // } else {
        //     sigma_j_inv = inv_sympd(sigmamat.cols(j*d,(j+1)*d-1));
        // }
        // rtilde = 0.5*(ydot * sigma_j_inv) % ydot;
        ytilde = solve(sigmamat.cols(j*d,(j+1)*d-1), ydot.t());
        rtilde = 0.5* ytilde.t() % ydot;
        r.col(j) = sum(rtilde,1);
        // s0j = sigma0mat.cols(j*d,(j+1)*d-1) * sigma_j_inv;
        s0j = solve(sigmamat.cols(j*d,(j+1)*d-1), sigma0mat.cols(j*d,(j+1)*d-1)).t();
        pen(j) = trace(s0j) - log(det(s0j)) -d;
      }
      if ( detsigma.has_nan() || sigmamat.has_nan() ) {
      // if ( any(detsigma < 1e-8) || detsigma.has_nan() || sigmamat.has_nan() ) {
          penloglik = R_NegInf;
        break;
      }
      alp_sig = alpha / sqrt(detsigma);
      minr = min(r,1);
      /* posterior for i = 1,...,n */
      /* normalizing with minr avoids the problem of dividing by zero */
      l_m = exp(-(r.each_col() - minr));
      l_m.each_row() %= alp_sig.t();
      sum_l_m = sum(l_m,1);
      w = l_m.each_col() / sum_l_m; /* w(j,i) = alp_j*l_j / sum_j (alp_j*l_j) */
      /* loglikelihood*/
      ll = sum(log(sum_l_m) - minr) - (double)n *d * M_LN_SQRT_2PI;
      /* subtract back minr and subtract n/2 times log(2pi) */;

      /* Compute the penalized loglik. Note that penalized loglik uses old (not updated) sigma */
      penloglik = ll + log(2.0) + fmin(log(tau),log(1-tau)) - an*sum(pen);

      diff = penloglik - oldpenloglik;
      oldpenloglik = penloglik;

      /* Normal exit */
      if (diff < tol ){
        break;
      }

      /* update alpha, mu, and sigma */
      dum=0;
      for (int j = 0; j < m; j++) {
        w_j = sum( w.col(j) ); /* w_j = sum_i w(i,j) */
        alpha(j) = w_j / n;
        mu_j = trans(sum((y.each_col() % w.col(j)), 0)) / w_j;
        ydot =  y.each_row() - mu_j.t();
        ssr_j = trans((ydot.each_col() % w.col(j))) * ydot;
        sigmamat_j = (2*an*sigma0mat.cols(j*d,(j+1)*d-1) + ssr_j)/ (2*an + w_j);
  //       sigma(j) = fmax(sigma(j),0.01*sigma0(j));
        /* If k ==1, impose lower and upper bound on the first element of mu_j */
        if (k==1) {
          mu_j(0) = fmin( fmax(mu_j(0),lb(j)), ub(j));
        }
        mu.subvec(j*d,(j+1)*d-1) = mu_j;
        detsigma(j) = det(sigmamat_j);
        sigmamat.cols(j*d,(j+1)*d-1) = sigmamat_j;
        for (int ii=0; ii<d; ++ii){
          for (int jj=ii; jj<d; ++jj){
            sigma(dum) = sigmamat_j(jj,ii);
            dum++;
          }
        }
      // alpha.print();
      // mu.print();
      // sigmamat.print();
        // sigmamat_j.print();
        // sigma.print();
      }

      /* for PMLE, we set k=0 (default value) */
      /* for EM test, we start from k=1       */
      /*   if k==1, we don't update tau       */
      /*   if k>1, we update tau              */
      if (k==1){
        alphah = (alpha(h-1)+alpha(h));
        alpha(h-1) = alphah*tau;
        alpha(h) = alphah*(1-tau);
      } else if (k>1) {
        alphah = (alpha(h-1)+alpha(h));
        tauhat = alpha(h-1)/(alpha(h-1)+alpha(h));
        if(tauhat <= 0.5) {
            tau = fmin((alpha(h-1)*n + 1.0)/(alpha(h-1)*n + alpha(h)*n + 1.0), 0.5);
        } else {
            tau = fmax(alpha(h-1)*n /(alpha(h-1)*n + alpha(h)*n + 1.0), 0.5);
        }
        alpha(h-1) = alphah*tau;
        alpha(h) = alphah*(1-tau);
      }

    /* Check singularity */
      if (any(alpha < 1e-8) || alpha.has_nan() ) {
      // if (any(alpha < 1e-8) || alpha.has_nan() || any(detsigma < 1e-8)) {
          sing = 1;
      }

      /* Exit from the loop if singular */
      if (sing) {
        notcg(jn) = 1;
        break;
      }

    } /* EM loop ends */

    penloglikset(jn) = penloglik;
    loglikset(jn) = ll;
    /* update b_jn */
    b_jn.subvec(0,m-1) = alpha;
    b_jn.subvec(m,m+m*d-1) = mu;
    b_jn.subvec(m+m*d,m+m*d+m*dsig-1) = sigma;
    /* update b */
    b.col(jn) = b_jn;

    post.col(jn) = vectorise(w);

  } /* end for (jn=0; jn<ninits; jn++) loop */

  return Rcpp::List::create(Named("penloglikset") = wrap(penloglikset),
                            Named("loglikset") = wrap(loglikset),
                            Named("notcg") = wrap(notcg),
                            Named("post") = wrap(post)
                            );
}
