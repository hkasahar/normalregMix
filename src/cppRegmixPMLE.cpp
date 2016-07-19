#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

const double SINGULAR_EPS = 10e-10; // criteria for matrix singularity

// [[Rcpp::export]]
List cppRegmixPMLE(NumericMatrix bs,
                      NumericVector ys,
                      NumericMatrix xs,
                      NumericMatrix zs,
                      NumericVector mu0s,
                      NumericVector sigma0s,
                      int m,
                      int p,
                      double an,
                      int maxit = 2000,
                      int ninits = 10,
                      double tol = 1e-8,
                      double tau = 0.5,
                      int h = 0,
                      int k = 0) {
  int n = ys.size();
  int q = xs.ncol();
  int q1 = q + 1;
  arma::mat b(bs.begin(), bs.nrow(), bs.ncol(), false);
  arma::vec y(ys.begin(), ys.size(), false);
  arma::mat x(xs.begin(), xs.nrow(), xs.ncol(), false);
  arma::mat z(zs.begin(), zs.nrow(), zs.ncol(), false);
  arma::vec mu0(mu0s.begin(), mu0s.size(), false);
  arma::vec sigma0(sigma0s.begin(), sigma0s.size(), false);
  arma::vec b_jn(bs.nrow());
  arma::vec lb(m), ub(m);
  arma::vec alpha(m), sigma(m), alp_sig(m);
  arma::mat mubeta(q1,m);
  arma::vec r(m), l_j(m);
  arma::mat w(m,n);
  arma::mat post(m*n,ninits);
  arma::vec notcg(ninits), penloglikset(ninits), loglikset(ninits);
  arma::vec gamma(p);
  arma::vec ytilde(n);
  arma::vec wtilde(n);
  arma::mat xtilde(n,q1);
  arma::mat ztilde(n,p);
  arma::mat zz(p,p);
  arma::mat ze(p,1);
  int emit, sing;
  double oldpenloglik, penloglik, s0j, diff, minr, w_j, sum_l_j, ll, ssr_j, alphah, tauhat;
  arma::mat x1(n,q1);

  x1.col(0) = arma::ones(n);
  for (int i=1; i < q1; ++i) {
    x1.col(i) = x.col(i-1);
  }

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
    for (int j=0; j < m; ++j){
      alpha(j) = b_jn(j);
      for (int i=0; i<q1; ++i){
        mubeta(i,j) = b_jn(m+q1*j+i);
      }
      sigma(j) = b_jn((q1+1)*m+j);
    }
    if (p>0) {
      for (int j=0; j < p; j++){
        gamma(j) = b_jn((q1+2)*m+j);
      }
    }
    oldpenloglik = R_NegInf;
    emit = 0;
    diff = 1.0;
    sing = 0;

    /* EM loop begins */
    for (int iter = 0; iter < maxit; iter++) {
      ll = - (double)n * M_LN_SQRT_2PI; /* n/2 times log(2pi) */
      alp_sig = alpha/sigma;

      if (p==0) {
        ytilde = y;
      } else {
        ytilde = y - z*gamma;
      }

      for (int i = 0; i < n; i++) {
        /* standardized squared residual */
        r = (1.0/sigma) % ( ytilde(i) - trans(x1.row(i)*mubeta) );
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
        wtilde = trans(w.row(j));
        for (int ii = 0; ii < q1; ii++) {
          xtilde.col(ii) = wtilde % x1.col(ii);
        }
        
        arma::mat design_matrix = trans(xtilde) * x1;
        if (rcond(design_matrix) < SINGULAR_EPS)
        {
          sing = 1;
          break;
        }
        
        mubeta.col(j) = solve( design_matrix , trans(xtilde) * ytilde );
        //mubeta.col(j) = solve( trans(xtilde) * x1 , trans(xtilde) * ytilde );
        ssr_j = sum( trans(w.row(j)) % pow(  ytilde - x1*mubeta.col(j) , 2 ) );
        sigma(j) = sqrt( (ssr_j + 2.0*an*sigma0(j)*sigma0(j))  / (w_j + 2.0*an) );
        sigma(j) = fmax(sigma(j),0.05*sigma0(j));
        /* If k ==1, impose lower and upper bound */
        if (k==1) {
          mubeta(0,j) = fmin( fmax(mubeta(0,j),lb(j)), ub(j));
        }
        alpha(j) = fmin( fmax(alpha(j),0.01), 0.99);
      }
      
      /* if k!=0, update alpha and/or tau */
      // if (k!=0){
      //   alphah = (alpha(h-1)+alpha(h));
      //   /* If k!=1, update tau. If k==1, no update of tau. */
      //   if (k!=1) {
      //     tauhat = alpha(h-1)/(alpha(h-1)+alpha(h));
      //     if (tauhat <= 0.5) {
      //       tau = fmin((alpha(h-1)*n + 1.0)/(alpha(h-1)*n + alpha(h)*n + 1.0), 0.5);
      //     } else {
      //       tau = fmax(alpha(h-1)*n /(alpha(h-1)*n + alpha(h)*n + 1.0), 0.5);
      //     }
      //   }
      //   /* Using tau, revise the h and h+1 th element of alpha */
      //   alpha(h-1) = alphah*tau;
      //   alpha(h) = alphah*(1-tau);
      // }
      if (k==1){
        alphah = (alpha(h-1)+alpha(h));
        /* Using tau, revise the h and h+1 th element of alpha */
        alpha(h-1) = alphah*tau;
        alpha(h) = alphah*(1-tau);
      }

      if (p>0) { /* update gamma */
        zz.zeros();
        ze.zeros();
        for (int j = 0; j < m; j++) {
          wtilde = trans(w.row(j)) ;
          for (int ii = 0; ii < p; ii++) {
            ztilde.col(ii) = wtilde % z.col(ii);
          }
          zz = zz + ( trans(ztilde) * z ) / (sigma(j)*sigma(j));
          ze = ze + ( trans(ztilde) * (y - x1*mubeta.col(j)) ) / (sigma(j)*sigma(j));
        }
        
        if (rcond(zz) < SINGULAR_EPS)
        {
          sing = 1;
          break;
        }

        gamma = solve(zz,ze);
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
      for (int i=0; i<q1; ++i){
        b_jn(m+q1*j+i) = mubeta(i,j);
      }
      b_jn((q1+1)*m+j) = sigma(j);
    }
    if (p>0) {
      for (int j=0; j < p; j++){
        b_jn((q1+2)*m+j) = gamma(j);
      }
    }
    b.col(jn) = b_jn; /* b is updated */
    post.col(jn) = vectorise(trans(w));

  } /* end for (jn=0; jn<ninits; jn++) loop */

  return Rcpp::List::create(Named("penloglikset") = wrap(penloglikset),
                            Named("loglikset") = wrap(loglikset),
                            Named("notcg") = wrap(notcg),
                              Named("post") = wrap(post)
  );
}

