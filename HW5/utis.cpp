#include <RcppArmadillo.h>
#include <omp.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double log_likelihood(vec beta, vec Y, mat X){
  double ll = 0;
  double XB = 0;
  for(int i=0;i<X.n_rows;i++){
    XB = 0;
    for(int j=0;j<4;j++){
      XB = XB + X(i,j)*beta(j);
    }
    ll = ll+Y(i)*(XB-log(exp(XB)+1))-(1-Y(i))*log(exp(XB)+1);
  }
  return(ll);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double log_prior(vec beta){
  double lp = 0;
  for(int i=0;i<4;i++){
    lp = lp-0.5*pow(beta[i],2)/10;
  }
  lp = lp - 4*(0.5*log(2*M_PI)+0.5*log(10));
  return(lp);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
mat Rcpp_MCMC(int niter, vec beta_init, vec Y, mat X, double step){
  mat beta_chain(niter, 4);
  vec beta_prop(4);
  double a, b;
  for(int i=0;i<niter;i++){
    for(int j=0;j<4;j++){
      beta_prop(j) = beta_init(j)+randn()*step;
      beta_chain(i,j) = beta_init(j);
    }
    a = log_likelihood(beta_prop, Y, X)+log_prior(beta_prop);
    b = log_likelihood(beta_init, Y, X)+log_prior(beta_init);
    if(log(randu())<a-b){
      beta_init = beta_prop;
    }
  }
  return(beta_chain);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
cube Rcpp_MCMC_parallel(int niter, int ncore, vec y, mat X, double step){
  cube beta_chain(ncore, niter, 4);
  mat X_shard;
  vec y_shard;
  int len_shard = X.n_rows/ncore;
  int k;
#pragma omp parallel shared(beta_chain) private(k, X_shard, y_shard)
{
  
}
  for(int k=0;k<ncore;k++){
    X_shard = X.submat(len_shard*k, 0, len_shard*(k+1)-1, 3);
    y_shard = y.subvec(len_shard*k, len_shard*(k+1)-1);
    beta_chain.row(k) = Rcpp_MCMC(niter, randn(4)*0.1, y_shard, X_shard, step);
  }
  return(beta_chain);
}