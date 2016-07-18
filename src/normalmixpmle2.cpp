#include <R.h>
#include <Rmath.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List normalmixpmle2(NumericMatrix bs,
           NumericVector ys,
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
  int n = ys.size();
  arma::mat b(bs.begin(), bs.nrow(), bs.ncol(), false);
  arma::colvec y(ys.begin(), ys.size(), false);
  arma::colvec mu0(mu0s.begin(), mu0s.size(), false);
  arma::colvec sigma0(sigma0s.begin(), sigma0s.size(), false);
  arma::colvec b_jn(bs.nrow());
  arma::colvec lb(m),ub(m);
  arma::colvec alpha(m), mu(m), sigma(m), alp_sig(m);
  arma::colvec z(m), l_j(m);
  arma::mat w(m,n);
  arma::mat post(m*n,ninits);
  arma::colvec notcg(ninits), penloglikset(ninits), loglikset(ninits);
  int emit, sing;
  double oldpenloglik, penloglik, s0j, diff, minz, w_j, sum_l_j, ll, ssr_j, alphah, tauhat;

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

  for (int jn=0; jn<ninits; jn++) {

    /* initialize EM iteration */
    b_jn = b.col(jn);
    for (int j=0; j < m; ++j){
      alpha(j) = b_jn(j);
      mu(j) = b_jn(m+j);
      sigma(j) = b_jn(2*m+j);
    }
    oldpenloglik = R_NegInf;
    emit = 0;
    diff = 1.0;
    sing = 0;

    /* EM loop begins */
    for (int iter = 0; iter < maxit; iter++) {
      ll = - (double)n * M_LN_SQRT_2PI; /* n/2 times log(2pi) */
      alp_sig = alpha/sigma;
      for (int i = 0; i < n; i++) {
        /* standardized squared residual */
        z = (1.0/sigma) % (y(i)-mu);
        z = 0.5 * (z % z); /* This is faster than z = pow( z, 2 ) */
        minz = min(z);
        /* posterior for i */
        /* normalizing with minr avoids the problem of dividing by zero */
        l_j =  alp_sig % exp( minz-z );
        sum_l_j = sum( l_j );
        w.col(i) = l_j/sum_l_j; /* w(i,j) = alp_j*like_j / sum_j (alp_j*like_j) */
        /* loglikelihood*/
        ll += log(sum_l_j) - minz; /* subtract back minz */
      } /* end for (i=0; i<n; i++) loop */

      /* Compute the penalized loglik. Note that penalized loglik uses old (not updated) sigma */
      penloglik = ll + log(2.0) + fmin(log(tau),log(1-tau));
      for (int j=0; j<m; j++) {
        s0j = sigma0(j)/sigma(j);
        penloglik += -an*(s0j*s0j - 2.0*log(s0j) -1.0);
      }
      diff = penloglik - oldpenloglik;
      oldpenloglik = penloglik;

      /* Normal exit */
      if (diff < tol || emit>=maxit){
        break;
      }
      emit++;

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

      /* if k!=0, update alpha and/or tau */
      if (k!=0){
        alphah = (alpha(h-1)+alpha(h));
        /* If k!=1, update tau. If k==1, no update of tau. */
        if (k!=1) {
          tauhat = alpha(h-1)/(alpha(h-1)+alpha(h));
          if (tauhat <= 0.5) {
            tau = fmin((alpha(h-1)*n + 1.0)/(alpha(h-1)*n + alpha(h)*n + 1.0), 0.5);
          } else {
            tau = fmax(alpha(h-1)*n /(alpha(h-1)*n + alpha(h)*n + 1.0), 0.5);
          }
        }
        /* Using tau, revise the h and h+1 th element of alpha */
        alpha(h-1) = alphah*tau;
        alpha(h) = alphah*(1-tau);
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

    }/* EM loop ends */

    penloglikset(jn) = penloglik;
    loglikset(jn) = ll;
    for (int j=0; j < m; j++){
      b_jn(j) = alpha(j);
      b_jn(m+j) = mu(j);
      b_jn(2*m+j) = sigma(j);
    }
    b.col(jn) = b_jn;
    post.col(jn) = vectorise(trans(w));

  } /* end for (jn=0; jn<ninits; jn++) loop */

  return Rcpp::List::create(Named("penloglikset") = wrap(penloglikset),
                            Named("loglikset") = wrap(loglikset),
                            Named("notcg") = wrap(notcg),
                              Named("post") = wrap(post)
  );
}

