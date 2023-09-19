#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double pllnormalmix_z_cpp(arma::vec theta,
                          arma::vec y,
                          arma::mat z,
                          int m,
                          int p,
                          double an,
                          arma::vec sigma0) {
  int n = y.size();
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
  pll = -(ll - an*sum(ss + log(1/ss) -1));
  
  return pll;
  
}

// [[Rcpp::export]]
arma::vec pllnormalmix_z_grad_cpp(arma::vec theta,
                          arma::vec y,
                          arma::mat z,
                          int m,
                          int p,
                          double an,
                          arma::vec sigma0) {
  int n = y.size();
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
double normalmix_heq_cpp(arma::vec theta,
                         arma::vec y,
                         arma::mat z,
                         int m,
                         int p,
                         double an,
                         arma::vec sigma0) {
  double h;
  h = sum(theta.subvec(0,m-1)) - 1.0;
  return h;
}

// [[Rcpp::export]]
arma::vec normalmix_heq_gr_cpp(arma::vec theta,
                               arma::vec y,
                               arma::mat z,
                               int m,
                               int p,
                               double an,
                               arma::vec sigma0) {

  arma::vec h(3*m+p);
  h.zeros();
  h.subvec(0, m-1).fill(1.0);
  
  return h;
}
