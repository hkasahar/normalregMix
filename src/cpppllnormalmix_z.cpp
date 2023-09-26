#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double emnormalmix_z_cpp(arma::vec gamma,
                     arma::vec y,
                     arma::mat z,
                     int m,
                     int p,
                     arma::vec theta,
                     double an,
                     arma::vec sigma0) {
  
  int n = y.n_elem;
  arma::vec alpha(m), mu(m), sigma(m), gamma0(p);
  arma::vec yhat(n), f(n), ytilde(n), wj(n);
  arma::mat w(n,m), ztilde(n,p);
  double ll, Wj;
  
  alpha = theta.subvec(0,m-1);
  mu    = theta.subvec(m,2*m-1);
  sigma = theta.subvec(2*m,3*m-1);
  gamma0 = theta.subvec(3*m,3*m+p-1);
  yhat = y - z * gamma0;
  
  // Initialize w matrix
  
  for (int j = 0; j < m; ++j) {
    w.col(j) = alpha(j)*normpdf(yhat, mu(j), sigma(j));
  }
  
  f = sum(w, 1);
  w.each_col()/=f;
    
  ll = 0.0;
  
  for (int j = 0; j < m; ++j) {
    wj = w.col(j);
    Wj = sum(wj);
    ytilde = y - sum(wj % y) / Wj;
    ztilde = z.each_row() - sum(z.each_col()%wj, 0)/Wj;
    ll -= (Wj+2.0*an)*log(sum(wj%(ytilde-ztilde*gamma)%(ytilde-ztilde*gamma))
                                 + 2.0*an*sigma0(j)*sigma0(j));
  }
  
  return -ll;
}

// [[Rcpp::export]]
arma::vec emnormalmix_z_grad_cpp(arma::vec gamma,
                         arma::vec y,
                         arma::mat z,
                         int m,
                         int p,
                         arma::vec theta,
                         double an,
                         arma::vec sigma0) {
  
  int n = y.n_elem;
  arma::vec alpha(m), mu(m), sigma(m), gamma0(p);
  arma::vec yhat(n), f(n), ytilde(n), wj(n), dll(p);
  arma::mat w(n,m), ztilde(n,p);
  double Wj;
  
  alpha = theta.subvec(0,m-1);
  mu    = theta.subvec(m,2*m-1);
  sigma = theta.subvec(2*m,3*m-1);
  gamma0 = theta.subvec(3*m,3*m+p-1);
  yhat = y - z * gamma0;
  
  // Initialize w matrix
  
  for (int j = 0; j < m; ++j) {
    w.col(j) = alpha(j)*normpdf(yhat, mu(j), sigma(j));
  }
  
  f = sum(w, 1);
  w.each_col()/=f;
  
  dll.fill(0);
  double denom;
  
  for (int j = 0; j < m; ++j) {
    wj = w.col(j);
    Wj = sum(wj);
    ytilde = y - sum(wj % y) / Wj;
    ztilde = z.each_row() - sum(z.each_col()%wj, 0)/Wj;
    denom = sum(wj%(ytilde-ztilde*gamma)%(ytilde-ztilde*gamma))
              + 2.0*an*sigma0(j)*sigma0(j);
    dll += (Wj+2.0*an)*2
            *sum((ztilde.each_col()%(wj%(ytilde-ztilde*gamma))),0)/denom;
  }
  
  return -dll;
}

// [[Rcpp::export]]
double pllnormalmix_z_cpp(arma::vec theta,
                          arma::vec y,
                          arma::mat z,
                          int m,
                          int p,
                          double an,
                          double cn,
                          arma::vec sigma0,
                          double tau = 0.5,
                          int h = 0) {
  int n = y.n_elem;
  arma::vec ytilde(n);
  arma::vec f(n);
  arma::vec alpha(m), mu(m), sigma(m), gamma(p), ss(m);
  double ll, pll;
  
  alpha = theta.subvec(0,m-1);
  mu = theta.subvec(m,2*m-1);
  sigma = theta.subvec(2*m,3*m-1);
  gamma = theta.subvec(3*m,3*m+p-1);
  ytilde = y - z * gamma;
  
  f.zeros();
  for (int j=0; j<m; j++) {
    f += alpha(j)*normpdf(ytilde, mu(j), sigma(j));
  }
  
  ll = sum(log(f));
  ss = sigma0/sigma;
  ss = pow(ss,2);
  pll = -(ll - an*sum(ss + log(1/ss) -1)  + cn*log(2.0) + cn*fmin(log(tau),log(1-tau)));
  
  return pll;
  
}

// [[Rcpp::export]]
arma::vec pllnormalmix_z_grad_cpp(arma::vec theta,
                          arma::vec y,
                          arma::mat z,
                          int m,
                          int p,
                          double an,
                          double cn,
                          arma::vec sigma0,
                          double tau = 0.5,
                          int h = 0) {
  int n = y.n_elem;
  arma::vec ytilde(n);
  arma::vec f(n);
  arma::vec alpha(m), mu(m), sigma(m), gamma(p), ss(m);
  arma::mat fall(n,m);
  arma::vec dll(3*m+p);
  
  alpha = theta.subvec(0,m-1);
  mu = theta.subvec(m,2*m-1);
  sigma = theta.subvec(2*m,3*m-1);
  gamma = theta.subvec(3*m,3*m+p-1);
  ytilde = y - z * gamma;
 
  for (int j=0; j<m; j++) {
     fall.col(j) = normpdf(ytilde, mu(j), sigma(j));
   }

  f = fall * alpha;
  
  arma::vec Dalphall(m), sss(m), Dmull(m), Dsigmall(m), Dgammall(p);
  arma::vec fj(n), ysj(n);
  arma::mat Dmuf(n,m), Dsigmaf(n,m), Dgammaf(n,p);
      
  for (int j=0; j<m; j++) {
        fj = fall.col(j);
        Dalphall(j) = sum(fj/f);
        Dmuf.col(j) = alpha(j)*(ytilde-mu(j))%fj/pow(sigma(j),2);
        ysj = (ytilde-mu(j))/sigma(j);
        Dsigmaf.col(j) = alpha(j)*(pow(ysj,2) - 1)%fj/sigma(j);
  }
    
    arma::colvec DumfSum = sum(Dmuf,1);
    Dgammaf  = z.each_col() % DumfSum;
      
    sss = pow(sigma0,2)/pow(sigma,3);
    Dmull = sum(Dmuf.each_col()/f, 0);
    Dsigmall = sum(Dsigmaf.each_col()/f, 0) - an*(-2*sss + 2/sigma);
        
    Dgammall = sum(Dgammaf.each_col()/f, 0);
        
    dll.subvec(0,m-1) = - Dalphall;
    dll.subvec(m,2*m-1) = - Dmull;
    dll.subvec(2*m,3*m-1) = - Dsigmall;
    dll.subvec(3*m,3*m+p-1) = - Dgammall;
        
    return dll;

}


// [[Rcpp::export]]
arma::vec normalmix_heq_cpp(arma::vec theta,
                         arma::vec y,
                         arma::mat z,
                         int m,
                         int p,
                         double an,
                         double cn,
                         arma::vec sigma0,
                         double tau = 0.5,
                         int h = 0) {
  arma::vec b(2);
  if (h==0) {
    b.subvec(0,0) = sum(theta.subvec(0,m-1)) - 1.0;
    b.shed_row(1);
  } else {
    b.subvec(0,0) = sum(theta.subvec(0,m-1)) - 1.0;
    b.subvec(1,1) = theta.subvec(h-1,h) 
                      - tau*sum(theta.subvec(h-1,h));
  }
  return b;
}

// [[Rcpp::export]]
arma::mat normalmix_heq_grad_cpp(arma::vec theta,
                               arma::vec y,
                               arma::mat z,
                               int m,
                               int p,
                               double an,
                               double cn,
                               arma::vec sigma0,
                               double tau = 0.5,
                               int h = 0) {

  arma::mat b(2,3*m+p);
  b.zeros();
  if (h==0) {
    b.submat(0,0,0,m-1).fill(1.0);
    b.shed_row(1);
  } else {
    b.submat(0,0,0,m-1).fill(1.0);
    b.submat(1,h-1,1,h-1) = 1-tau;
    b.submat(1,h,1,h) = -tau;
  }
  return b;
}
