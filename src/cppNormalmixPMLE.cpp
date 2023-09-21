#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @description Updates parameter estimates of a finite mixture of
//' univariate normals by the EM algorithm.
//' @export
//' @title cppNormalmixPMLE
//' @name cppNormalmixPMLE
//' @param bs 3m + p by ninits matrix of initial values of (alpha,mu,sigma,gamma).
//' @param ys n by 1 vector of data.
//' @param mu0s m-1 vector of the estimate of mu from an m-1 component model.
//' @param sigma0s m-1 vector of the estimate of sigma from an m-1 component model.
//' @param m number of components in the mixture.
//' @param an tuning parameter.
//' @param cn tuning parameter.
//' @param maxit maximum number of iterations.
//' @param ninits number of initial values.
//' @param tol Convergence is declared when the penalized log-likelihood increases by less than \code{tol}.
//' @param tau tau used to split the h-th component.
//' @param h h used as index for pivoting.
//' @param k number of EM steps taken in computing the modified EM statisic.
//' @return  A list with items:
//' \item{penloglikset}{vector of the maximized value of the penalized log-likelihood.}
//' \item{loglikset}{vector of the maximized value of the log-likelihood.}
//' \item{notcg}{vector that records whether EM steps converged or not for each initial value.}
//' \item{post}{n*m by ininits matrix of posterior probabilities for observations.}
//'
// [[Rcpp::export]]
List cppNormalmixPMLE(NumericMatrix bs,
                      NumericVector ys,
                      NumericVector mu0s,
                      NumericVector sigma0s,
                      int m,
                      double an,
                      double cn,
                      int maxit = 2000,
                      int ninits = 10,
                      double tol = 1e-8,
                      double tau = 0.5,
                      int h = 0,
                      int k = 0) {
  int n = ys.size();
  arma::mat b(bs.begin(), bs.nrow(), bs.ncol(), false);
  arma::vec y(ys.begin(), ys.size(), false);
  arma::vec mu0(mu0s.begin(), mu0s.size(), false);
  arma::vec sigma0(sigma0s.begin(), sigma0s.size(), false);
  arma::vec b_jn(bs.nrow());
  arma::vec lb(m),ub(m);
  arma::vec alpha(m), mu(m), sigma(m), alp_sig(m);
  arma::vec r(m), l_j(m);
  arma::mat w(m,n);
  arma::mat post(m*n,ninits);
  arma::vec notcg(ninits), penloglikset(ninits), loglikset(ninits);
  arma::vec wtilde(n);
  int sing;
  double oldpenloglik, s0j, diff, minr, w_j, sum_l_j, ssr_j, betah, tauhat;
  double ll = 0; // force initilization
  double penloglik = 0; // force initialization
  notcg.zeros();  // initialization

  /* Lower and upper bound for mu */
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
    mu = b_jn.subvec(m,2*m-1);
    sigma = b_jn.subvec(2*m,3*m-1);
    oldpenloglik = R_NegInf;
    diff = 1.0;
    sing = 0;

    /* EM loop begins */
    for (int iter = 0; iter < maxit; iter++) {
      ll = - (double)n * M_LN_SQRT_2PI; /* n/2 times log(2pi) */
      alp_sig = alpha/sigma;

      for (int i = 0; i < n; i++) {
        /* standardized squared residual */
        r = (1.0/sigma) % (y(i) - mu);
        r = 0.5 * (r % r); /* This is faster than r = pow( r, 2.0 ) */
        minr = min(r);
        /* posterior for i */
        /* normalizing with minr avoids the problem of dividing by zero */
        l_j =  alp_sig % exp( minr-r );
        sum_l_j = sum( l_j );
        w.col(i) = l_j/sum_l_j; /* w(j,i) = alp_j*l_j / sum_j (alp_j*l_j) */
        /* loglikelihood*/
        ll +=  log(sum_l_j) - minr; /* subtract back minr */
      } /* end for (i=0; i<n; i++) loop */

      /* Compute the penalized loglik. */
      /* Here, we compute penalized loglik with old (not updated) sigma. */
      penloglik = ll + cn*log(2.0) + cn*fmin(log(tau),log(1-tau));
      for (int j=0; j<m; j++) {
        s0j = sigma0(j)/sigma(j);
        penloglik += -an*(s0j*s0j - 2.0*log(s0j) -1.0);
      }
      diff = penloglik - oldpenloglik;
      oldpenloglik = penloglik;

      /* Normal exit */
      if (diff < tol ){
        break;
      }

      /* update alpha, mu, and sigma */
      for (int j = 0; j < m; j++) {
        w_j = sum( w.row(j) ); /* w_j(j) = sum_i w(i,j) */
        alpha(j) = w_j / n;
        mu(j) = sum( trans(w.row(j)) % y ) / w_j;
        ssr_j = sum( trans(w.row(j)) % pow( y - mu(j), 2 ) );
        sigma(j) = sqrt( (ssr_j + 2.0*an*sigma0(j)*sigma0(j))  / (w_j + 2.0*an) );
        sigma(j) = fmax(sigma(j),0.01*sigma0(j));
        /* If k ==1, impose lower and upper bound */
        if (k==1) {
          mu(j) = fmin( fmax(mu(j),lb(j)), ub(j));
        }
      }
 
      /* for PMLE, we set k=0 (default value) */
      /* for EM test, we start from k=1       */
      /*   if k==1, we don't update tau       */
      /*   if k>1, we update tau              */
      if (k==1){
        betah = alpha(h-1)+alpha(h);
        alpha(h-1) = betah*tau;
        alpha(h) = betah*(1-tau);
      } else if (k>1) {
        betah = alpha(h-1)+alpha(h);
        tauhat = alpha(h-1)/betah;
        if(tauhat <= 0.5) {
            tau = fmin((alpha(h-1)*n + cn)/(betah*n + cn), 0.5);
        } else {
            tau = fmax(alpha(h-1)*n /(betah*n + cn), 0.5);
        }
        alpha(h-1) = betah*tau;
        alpha(h) = betah*(1-tau);
      }

      /* Check singularity */
      for (int j=0; j<m; j++) {
        if (alpha(j) < 1e-8 || std::isnan(alpha(j)) || sigma(j) < 1e-8){
          sing = 1;
        }
      }

      /* Exit from the loop if singular */
      if (sing) {
        notcg(jn) = 1;
        break;
      }

    } /* EM loop ends */

    /* Compute loglik and penalized loglik under the updated parameter value */
    alp_sig = alpha/sigma;
    
    ll = - (double)n * M_LN_SQRT_2PI; /* n/2 times log(2pi) */
    for (int i = 0; i < n; i++) {
      /* standardized squared residual */
      r = (1.0/sigma) % (y(i) - mu);
      r = 0.5 * (r % r); /* This is faster than r = pow( r, 2.0 ) */
      minr = min(r);
      /* normalizing with minr avoids the problem of dividing by zero */
      l_j =  alp_sig % exp( minr-r );
      sum_l_j = sum( l_j );
      /* loglikelihood*/
      ll +=  log(sum_l_j) - minr; /* subtract back minr */
    } /* end for (i=0; i<n; i++) loop */
      
    penloglik = ll + cn*log(2.0) + cn*fmin(log(tau),log(1-tau));
    for (int j=0; j<m; j++) {
      s0j = sigma0(j)/sigma(j);
      penloglik += -an*(s0j*s0j - 2.0*log(s0j) -1.0);
    }

    if (sing) {
      ll = R_NegInf;
      penloglik = R_NegInf;
    }
    
    penloglikset(jn) = penloglik;
    loglikset(jn) = ll;
    b_jn.subvec(0,m-1) = alpha;
    b_jn.subvec(m,2*m-1) = mu;
    b_jn.subvec(2*m,3*m-1) = sigma;
    
    b.col(jn) = b_jn; /* b is updated */
    post.col(jn) = vectorise(trans(w));

  } /* end for (jn=0; jn<ninits; jn++) loop */

  return Rcpp::List::create(Named("penloglikset") = wrap(penloglikset),
                            Named("loglikset") = wrap(loglikset),
                            Named("notcg") = wrap(notcg),
                            Named("post") = wrap(post)
  );
}
