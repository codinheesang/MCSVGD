// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-


// we only include RcppArmadillo.h which pulls Rcpp.h in for us
// [[Rcpp::depends("RcppArmadillo")]]

#include <RcppArmadillo.h>
#include <limits>
#include <omp.h>
#include <Rcpp.h>
#include <cmath>
//#define minimum(x,y) (((x) < (y)) ? (x) : (y))

using namespace std;
using namespace Rcpp;
using namespace arma;



// =============================================================================
//  Rejection sampling for COMP distribution
// =============================================================================
// Evaluate unnormalised density of the COM-poisson distribution 
// if fudge is set to lgammafn(mode+1) then the unnormalised density is one at the mode
// If mode and fudge are set to 0 then the usual unnormalised density is computed
// [[Rcpp::export]]
double unnorm_ldcpois(double x, double mu, double nu, double mode, double fudge) {
  return nu*((x-mode)*log(mu)-lgamma(x+1)+fudge);
}


// Sample from a geometric distribution truncated to {0,1,...,n}
// u is a U[0,1] realisation
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double truncated_geo_sample(double u, double logq, double n) {
  double C;
  if(logq > -DBL_EPSILON){
    return 0;
  } else {
    C = -expm1(logq*(n+1));
    return floor(log(1-C*u)/logq); 
  }
}


// Sample from a geometric distribution with range {0,1,2,...}
// u is a U[0,1] realisation
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double untruncated_geo_sample(double u, double logq) {
  if (logq > -DBL_EPSILON){
    return 0;
  } else {
    return floor(log(u)/logq); 
  }
}


// Compute the finite geometric series (1+q+...+q^deltax)
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double truncated_lweights(double deltax, double logq) {
  if (logq > -DBL_EPSILON)
    return log(deltax+1)+logq;
  return log1p(-exp((deltax+1)*logq)) - log1p(-exp(logq));
}


// Compute the geometric series (1+q+...)
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
vec Summary(mat X, vec y){
  int nrow = X.n_rows;
  int ncol = X.n_cols;
  vec result = zeros(ncol+1);
  
  //calculate summary statistic
  for(int j=0; j<ncol; j++){
    for(int i=0; i<nrow; i++){
      result(j) = result(j) + X(i,j)*y(i);
    }
  }
  for(int ii=0; ii<nrow; ii++){
    result(ncol) = result(ncol) + lgamma(y(ii)+1);
  }
  return result;
}







// Compute the log gradient of h(x|theta)
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
vec loggrah(mat X, vec y, double nu, vec theta){
  int nrow = X.n_rows;
  int ncol = X.n_cols;
  vec result = zeros(ncol+1);
  vec value = zeros(ncol+1);
  
  //calculate summary statistic
  for(int j=0; j<ncol; j++){
    for(int i=0; i<nrow; i++){
      result(j) = result(j) + X(i,j)*y(i);
    }
  }
  for(int ii=0; ii<nrow; ii++){
    result(ncol) = result(ncol) + lgamma(y(ii)+1);
  }
  
  for(int jj=0; jj<ncol; jj++){
    value(jj) = nu*result(jj);
  }
  
  for(int iii=0; iii<ncol; iii++){
    value(ncol) = value (ncol)+ theta(iii)*result(iii);
  }
  
  value(ncol) = value(ncol)-result(ncol);
  return value;
}


// Compute the log gradient of h(x|theta)
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
vec loggrah3(mat X, vec summary, double nu, vec theta){
  int nrow = X.n_rows;
  int ncol = X.n_cols;
  vec value = zeros(ncol+1);
  double result = 0;
  
  for(int i=0; i<(ncol); i++){
    value(i)=nu*summary(i);
  }
  
  for(int j=0; j<(ncol); j++){
    result=result+theta(j)*summary(j);
  }
  
  value(ncol)= result-summary(ncol);
  return value;
}
// Compute the log gradient of h(x|theta)
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
vec loggrah2(mat X, vec y, vec theta){
  int nrow = X.n_rows;
  int ncol = X.n_cols;
  vec result = zeros(ncol+1);
  vec value = zeros(ncol+1);
  
  //calculate summary statistic
  for(int j=0; j<ncol; j++){
    for(int i=0; i<nrow; i++){
      result(j) = result(j) + X(i,j)*y(i);
    }
  }
  for(int ii=0; ii<nrow; ii++){
    result(ncol) = result(ncol) + lgamma(y(ii)+1);
  }
  
  for(int jj=0; jj<ncol; jj++){
    value(jj) = theta(ncol)*result(jj);
  }
  
  for(int iii=0; iii<ncol; iii++){
    value(ncol) = value (ncol)+ theta(iii)*result(iii);
  }
  
  value(ncol) = value(ncol)-result(ncol);
  return value;
}



// Compute h(x|theta)
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double hfunc(vec suf, double nu, vec theta){
  int ncol = theta.size();
  double value = 0;
  
  for(int jj=0; jj<ncol-1; jj++){
    value = value + nu*theta(jj)*suf(jj);
  }
  
  value = value - nu*suf(ncol-1);
  
  value = value;
  return value;
}




// Compute h(x|theta)
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double hfunc2(vec suf, vec theta){
  int ncol = theta.size();
  double value = 0;
  
  for(int jj=0; jj<ncol-1; jj++){
    value = value + theta(ncol)*theta(jj)*suf(jj);
  }
  
  value = value - theta(ncol)*suf(ncol-1);
  
  value = exp(value);
  return value;
}


// Compute the geometric series (1+q+...)
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
vec Summary_1dim(mat X, vec y){
  int nrow = X.n_rows;
  int ncol = X.n_cols;
  vec result = zeros(ncol+1);
  
  //calculate summary statistic
  for(int j=0; j<ncol; j++){
    for(int i=0; i<nrow; i++){
      result(j) = result(j) + X(i,j)*y(i);
    }
  }
  for(int ii=0; ii<nrow; ii++){
    result(ncol) = result(ncol) + lgamma(y(ii)+1);
  }
  return result;
}


// Calculate sufficient statistic
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double untruncated_lweights(double logq) {
  return -log1p(-exp(logq));
}



// Sample from an COMP distribution
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
vec rCOMP(int n, double mu, double nu){ // mu = mode
  double negInf = -std::numeric_limits<float>::infinity();;
  
  double logmu, lmode, rmode, fudge, sd, lsd, rsd, maxlweight, logprob, x, u;
  int i, attempts;
  vec ldens(4), lweights(4), logq(4), sweights(4), result(n);
  logmu = log(mu);
  
  // Figure out mode and standard deviation
  lmode = ceil(mu)-1;
  fudge = lgamma(lmode+1);
  rmode = lmode+1;
  fudge = lgamma(lmode+1);
  sd = ceil(sqrt(mu)/sqrt(nu));
  if (sd < 5) {
    sd = 5;
  }
  
  // Set up two points at mode +/- sd
  lsd = round(lmode-sd);
  if (lsd < 0){
    lsd = -1;
  }
  rsd = round(rmode+sd);
  
  // Left most tail
  if (lsd == -1) {
    lweights[0] = negInf;
    logq[0] = 0;
    ldens[0] = negInf;
  } else {
    ldens[0] = unnorm_ldcpois(lsd, mu, nu, lmode, fudge);
    if (lsd == 0) {
      lweights[0] = ldens[0];
      logq[0] = 0;
    } else {
      logq[0] = nu * (-logmu + log(lsd));
      lweights[0] = ldens[0] + truncated_lweights(lsd, logq[0]);
    }
  }
  
  // within 1sd to the left of the mode
  ldens[1] = 0;
  if (lmode == 0) {
    lweights[1] = 0;
    logq[1] = 1;
  } else {
    logq[1] = nu * (-logmu + log(lmode));
    lweights[1] = truncated_lweights(lmode-lsd-1, logq[1]);
  }
  
  // within 1sd to the right of the mode
  logq[2] = nu * (logmu - log(rmode+1));
  ldens[2] = nu * (logmu - log(rmode));
  lweights[2] = ldens[2] + truncated_lweights(rsd-rmode-1, logq[2]);
  
  // right tail
  logq[3] = nu * (logmu - log(rsd+1));
  ldens[3] = unnorm_ldcpois(rsd, mu, nu, lmode, fudge);
  lweights[3] = ldens[3] + untruncated_lweights(logq[3]);
  
  // Find maximum log-weight
  maxlweight = lweights[0];
  for (i = 1; i < 4; i++){
    if (lweights[i] > maxlweight) { maxlweight = lweights[i]; }
  }
  // Compute the cumulative sum of the weights
  for (i = 0; i < 4; i++) {
    lweights[i] = lweights[i]-maxlweight;
    sweights[0] = exp(lweights[0]);
  }
  for (i = 1; i < 4; i++) {
    sweights[i]=sweights[i-1]+exp(lweights[i]);
  }
  
  // Draw the sample by rejection sampling
  attempts = 0;
  for (i = 0; i < n; i++) {
    while (TRUE) {
      attempts = attempts + 1;
      u = randu() * sweights[3];
      if (u < sweights[0]) {
        u = u / sweights[0];
        x = truncated_geo_sample(u, logq[0], lsd);
        logprob = ldens[0]+x*logq[0];
        x = lsd-x;
      } else {
        if (u < sweights[1]) {
          u = (u-sweights[0])/(sweights[1]-sweights[0]);
          x = truncated_geo_sample(u, logq[1], lmode-lsd-1);
          logprob = ldens[1]+x*logq[1];
          x = lmode - x;
        } else {
          if (u<sweights[2]) {
            u = (u-sweights[1])/(sweights[2]-sweights[1]);
            x = truncated_geo_sample(u, logq[2], rsd-rmode-1);
            logprob = ldens[2]+x*logq[2];
            x = rmode + x;
          } else {
            u = (u-sweights[2])/(sweights[3]-sweights[2]);
            x = untruncated_geo_sample(u, logq[3]);
            logprob = ldens[3]+x*logq[3];
            x = rsd + x;
          }
        }
      }
      if (log(randu()) < unnorm_ldcpois(x, mu, nu, lmode, fudge) - logprob) {
        result[i] = x;
        break;
      }
    }
  }
  return result;
}



// Auxiliary variable generation in parallel
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
mat rCOMP_parallel(int n, vec mode, double nu, int num){
  int N = mode.size();
  mat aux = zeros(n, N);
#pragma omp parallel num_threads(num)
{
#pragma omp for
  for(int i = 0; i < N; i++){
    aux.col(i) = rCOMP(n, mode[i], nu);
  }
}
return aux;
}


// Auxiliary variable generation in parallel
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
mat rCOMP_one(int n, vec mode, double nu){
  int N = mode.size();
  mat aux = zeros(n, N);
  for(int i = 0; i < N; i++){
    aux.col(i) = rCOMP(n, mode[i], nu);
  }
  return aux;
}

// Save one auxiliary dataset
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
mat rCOMP_one2(vec mode, double nu){
  int n = mode.size();
  mat aux = zeros(1, n);
  
  for(int i=0; i<n; i++){
    aux.col(i)=rCOMP(1, mode[i], nu)[0];
  }
  return aux;
}

// =============================================================================
// Distribution functions
// =============================================================================

// Unnormalized log likelihood of COMP
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double modeCOMP_logh(vec y, vec mode, double nu){
  double result = sum( nu * ( y % log(mode) -  lgamma(y+1) ) );
  return result;
}


// Unnormalized log multivariate normal
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double MVN_logh(vec y, vec mu, mat invSigma){
  vec result = - 0.5 * trans(y - mu) * invSigma * (y - mu);
  return result[0];
} 



// =============================================================================
// COMP Regression Models
// =============================================================================

// -----------------------------------------------------------------------------
// Exchange algorithm for inference
// -----------------------------------------------------------------------------
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List exactCOMPreg(int outer, vec y, mat X, vec beta, double alpha,
                  List block, vec sigma2, List COV, bool updateCOV,
                  int adaptInterval, double adaptFactorExponent, int adapIter,
                  int thin){
  
  int p = beta.size(), N = y.size(), iter = 0, m1 = block.size(), mj;
  mat posterior(outer, p+1), posterior_thined, accprob = zeros(outer, m1), prej, postj;
  vec betaprop(p), mode(N), modeprop(N), par(p+1), parprop(p+1), aux(N), dummy;
  vec rhat = zeros(m1), gamma1 = zeros(m1), gamma2 = zeros(m1), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  double logprob, u, alphaprop, nu, nuprop;
  uvec index_beta(p);
  for(int i = 0; i < p; i++){ index_beta[i] = i; }
  int index_alpha = p;
  par.elem(index_beta) = beta;
  par[index_alpha] = alpha;
  
  List cholCOV(m1);
  for(int i = 0; i < m1; i++){
    uvec index = block[i];
    mj = index.size();
    cholCOV[i] = trans( chol( sigma2[i] * ( as<arma::mat>(COV[i]) + 0.001 * diagmat(ones(mj)) ) ) );
  }
  
  
  // initial parameter values
  mode = exp(X * beta);
  nu = exp(alpha);
  
  // Start of MCMC Chain
  for(int k = 0; k < outer; k++){
    
    for(int j = 0; j < m1; j ++){
      
      uvec index = block[j];
      mj = index.size();
      
      // update proposal distribution
      if( updateCOV ){
        
        if( (k+1 >= adaptInterval) && (k+1 - (adaptInterval * trunc((k+1) / adaptInterval)) == 0) ){
          dummyaccprob = accprob.col(j);
          rhat[j] = sum( dummyaccprob.rows(k+1-adaptInterval, k-1) ) / (adaptInterval-1);
          gamma1[j] = 1 / pow(adapIter, c1);
          gamma2[j] = c0 * gamma1[j];
          sigma2[j] = exp( log(sigma2[j]) + gamma2[j] * (rhat[j] - ropt) );
          
          postj = posterior.cols(index);
          COV[j] = as<arma::mat>(COV[j]) + gamma1[j] * ( cov( postj.rows(k+1-adaptInterval, k-1) ) - as<arma::mat>(COV[j]) );
          
          cholCOV[j] = trans( chol( sigma2[j] * ( as<arma::mat>(COV[j]) + 0.001 * diagmat(ones(mj)) ) ) );
          
          if( j == m1-1 ){ adapIter = adapIter + 1; }
        }
      } 
      
      
      dummy = as<arma::mat>(cholCOV[j]) * randn(mj);
      parprop = par;
      for(int l = 0; l < mj; l ++){
        parprop[index[l]] = par[index[l]] + dummy[l];
      }
      
      betaprop = parprop.elem(index_beta);
      alphaprop = parprop[index_alpha];
      
      modeprop = exp( X * betaprop );
      nuprop = exp(alphaprop);
      
      // aux = trans( rCOMP_parallel(1, modeprop, nuprop, num) );
      for(int i = 0; i < N; i ++){
        aux[i] = ( rCOMP(1, modeprop[i], nuprop) )[0];
      }
      
      logprob = modeCOMP_logh(y, modeprop, nuprop) - modeCOMP_logh(y, mode, nu) +
        modeCOMP_logh(aux, mode, nu) - modeCOMP_logh(aux, modeprop, nuprop) +
        MVN_logh(parprop.elem(index), zeros(mj), diagmat(ones(mj))/100) -
        MVN_logh(par.elem(index), zeros(mj), diagmat(ones(mj))/100);
      
      u = log( randu() );
      if( u < logprob ){
        par = parprop;
        beta = betaprop;
        alpha = alphaprop;
        mode = modeprop;
        nu = nuprop;
        
        accprob(k, j) = 1;
      }
      
      for(int l = 0; l < mj; l ++){
        posterior(k,index[l]) = par[index[l]];
      }
    }
    
    // thining
    if(k - (thin * trunc(k / thin)) == 0){
      iter = iter + 1;
      posterior_thined.insert_rows(iter-1, posterior.row(k));
    }
    
  }
  
  return Rcpp::List::create(Rcpp::Named("Sample") = posterior_thined,
                            Rcpp::Named("Accprob") = accprob,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("COV") = COV);
}












// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List exactCOMPreg_fixnu(int outer, vec y, mat X, vec beta, double alpha,
                        List block, vec sigma2, List COV, bool updateCOV,
                        int adaptInterval, double adaptFactorExponent, int adapIter,
                        int thin){
  
  int p = beta.size(), N = y.size(), iter = 0, m1 = block.size(), mj;
  mat posterior(outer, p+1), posterior_thined, accprob = zeros(outer, m1), prej, postj;
  vec betaprop(p), mode(N), modeprop(N), par(p+1), parprop(p+1), aux(N), dummy;
  vec rhat = zeros(m1), gamma1 = zeros(m1), gamma2 = zeros(m1), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  double logprob, u, alphaprop, nu, nuprop;
  uvec index_beta(p);
  for(int i = 0; i < p; i++){ index_beta[i] = i; }
  int index_alpha = p;
  par.elem(index_beta) = beta;
  par[index_alpha] = alpha;
  
  List cholCOV(m1);
  for(int i = 0; i < m1; i++){
    uvec index = block[i];
    mj = index.size();
    cholCOV[i] = trans( chol( sigma2[i] * ( as<arma::mat>(COV[i]) + 0.001 * diagmat(ones(mj)) ) ) );
  }
  
  
  // initial parameter values
  mode = exp(X * beta);
  nu = 1.648721;
  
  // Start of MCMC Chain
  for(int k = 0; k < outer; k++){
    
    for(int j = 0; j < m1; j ++){
      
      uvec index = block[j];
      mj = index.size();
      
      // update proposal distribution
      if( updateCOV ){
        
        if( (k+1 >= adaptInterval) && (k+1 - (adaptInterval * trunc((k+1) / adaptInterval)) == 0) ){
          dummyaccprob = accprob.col(j);
          rhat[j] = sum( dummyaccprob.rows(k+1-adaptInterval, k-1) ) / (adaptInterval-1);
          gamma1[j] = 1 / pow(adapIter, c1);
          gamma2[j] = c0 * gamma1[j];
          sigma2[j] = exp( log(sigma2[j]) + gamma2[j] * (rhat[j] - ropt) );
          
          postj = posterior.cols(index);
          COV[j] = as<arma::mat>(COV[j]) + gamma1[j] * ( cov( postj.rows(k+1-adaptInterval, k-1) ) - as<arma::mat>(COV[j]) );
          
          cholCOV[j] = trans( chol( sigma2[j] * ( as<arma::mat>(COV[j]) + 0.001 * diagmat(ones(mj)) ) ) );
          
          if( j == m1-1 ){ adapIter = adapIter + 1; }
        }
      } 
      
      
      dummy = as<arma::mat>(cholCOV[j]) * randn(mj);
      parprop = par;
      for(int l = 0; l < mj; l ++){
        parprop[index[l]] = par[index[l]] + dummy[l];
      }
      
      betaprop = parprop.elem(index_beta);
      alphaprop = parprop[index_alpha];
      
      modeprop = exp( X * betaprop );
      nuprop = 1.648721;
      
      // aux = trans( rCOMP_parallel(1, modeprop, nuprop, num) );
      for(int i = 0; i < N; i ++){
        aux[i] = ( rCOMP(1, modeprop[i], nuprop) )[0];
      }
      
      logprob = modeCOMP_logh(y, modeprop, nuprop) - modeCOMP_logh(y, mode, nu) +
        modeCOMP_logh(aux, mode, nu) - modeCOMP_logh(aux, modeprop, nuprop) +
        MVN_logh(parprop.elem(index), zeros(mj), diagmat(ones(mj))/100) -
        MVN_logh(par.elem(index), zeros(mj), diagmat(ones(mj))/100);
      
      u = log( randu() );
      if( u < logprob ){
        par = parprop;
        beta = betaprop;
        alpha = alphaprop;
        mode = modeprop;
        nu = nuprop;
        
        accprob(k, j) = 1;
      }
      
      for(int l = 0; l < mj; l ++){
        posterior(k,index[l]) = par[index[l]];
      }
    }
    
    // thining
    if(k - (thin * trunc(k / thin)) == 0){
      iter = iter + 1;
      posterior_thined.insert_rows(iter-1, posterior.row(k));
    }
    
  }
  
  return Rcpp::List::create(Rcpp::Named("Sample") = posterior_thined,
                            Rcpp::Named("Accprob") = accprob,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("COV") = COV);
}







// =============================================================================
// Auxiliary variable simulation
// =============================================================================

// Save sufficient statistics only
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
mat rCOMP2_parallel(mat X, vec theta, int N, int num){
  int n = X.n_rows, p = theta.size();
  vec mode, beta = theta.rows(0, p-2);
  mat aux(N, n), sumstat(N, p);
  double alpha = theta[p-1], nu;
  mode = exp(X * beta);
  //nu = exp(alpha);
  nu = alpha;
  aux = rCOMP_parallel(N, mode, nu, num);
  sumstat.cols(0, p-2) = aux * X;
  sumstat.col(p-1) = sum(lgamma(aux+1), 1); 
  
  return sumstat;
}




// Save sufficient statistics only
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
mat rCOMP3_parallel(mat X, vec theta, int N, int num){
  int n = X.n_rows, p = theta.size();
  vec mode, beta = theta.rows(0, p-2);
  mat aux(N,n);
  double alpha = theta[p-1], nu;
  mode = exp(X * beta);
  //nu = exp(alpha);
  nu = alpha;
  aux = rCOMP_parallel(N, mode, nu, num);
  return aux;
}


// Save sufficient statistics only
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
mat rCOMP_one_suf(mat X, vec theta){
  int n = X.n_rows, p = theta.size();
  vec mode, beta = theta.rows(0, p-2);
  double alpha = theta[p-1], nu;
  mode = exp(X * beta);
  //nu = exp(alpha);
  nu = alpha;
  mat aux = rCOMP_one2(mode, nu);
  mat sumstat = zeros(1,p);
  sumstat.cols(0, p-2) = aux * X;
  sumstat.col(p-1) = sum(lgamma(aux+1), 1); 
  return sumstat;
}


// Save sufficient statistics only
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
mat rCOMP2_one(mat X, vec theta, int N){
  int n = X.n_rows, p = theta.size();
  vec mode, beta = theta.rows(0, p-2);
  mat aux(N,n);
  double alpha = theta[p-1], nu;
  mode = exp(X * beta);
  //nu = exp(alpha);
  nu = alpha;
  aux = rCOMP_one(N, mode, nu);
  return aux;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
mat rCOMP2(mat X, vec theta, int N){
  int n = X.n_rows, p = theta.size();
  vec mode, beta = theta.rows(0, p-2);
  mat aux(N, n), sumstat(N, p);
  double alpha = theta[p-1], nu;
  
  mode = exp(X * beta);
  //nu = exp(alpha);
  nu = alpha;
  
  for(int i = 0; i < n; i++){
    aux.col(i) = rCOMP(N, mode[i], nu);
  }
  sumstat.cols(0, p-2) = aux * X;
  sumstat.col(p-1) = sum(lgamma(aux+1), 1); 
  
  return sumstat;
}
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
vec rCOMP2_parallel_all(mat X, vec theta, int N, int num){
  int n = X.n_rows, p = theta.size();
  vec mode, beta = theta.rows(0, p-2);
  mat aux(N, n), sumstat(N, p);
  vec sumstatall = zeros(p);
  double alpha = theta[p-1], nu;
  
  mode = exp(X * beta);
  nu = exp(alpha);
  
  aux = rCOMP_parallel(N, mode, nu, num);
  sumstat.cols(0, p-2) = aux * X;
  sumstat.col(p-1) = sum(lgamma(aux+1), 1); 
  
  sumstatall = sum(sumstat, 0);
  sumstatall = sumstatall/N;
  return sumstatall;
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
vec rCOMP_aux(mat X, vec theta, int N){
  int n = X.n_rows, p = theta.size();
  vec mode, beta = theta.rows(0, p-2);
  vec aux(n);
  double alpha = theta[p-1], nu;
  
  mode = exp(X * beta);
  nu = exp(alpha);
  
  for(int i = 0; i < N; i ++){
    aux[i] = ( rCOMP(1, mode[i], nu) )[0];
  }
  vec sums = Summary_1dim(X,aux);
  
  return sums;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
vec rCOMP_aux_all(mat X, vec theta, int N){
  int n = X.n_rows, p = theta.size();
  vec mode, beta = theta.rows(0, p-2);
  mat aux(N, n), sumstat(N, p);
  vec sumstatall = zeros(p);
  double alpha = theta[p-1], nu;
  
  mode = exp(X * beta);
  nu = exp(alpha);
  
  for(int i = 0; i < N; i ++){
    aux[i] = ( rCOMP(1, mode[i], nu) )[0];
  }
  sumstat.cols(0, p-2) = aux * X;
  sumstat.col(p-1) = sum(lgamma(aux+1), 1); 
  
  sumstatall = sum(sumstat, 0);
  sumstatall = sumstatall/N;
  return sumstatall;
}



mat crossdist(mat m1, mat m2) {
  int nrow1 = m1.n_rows;
  int nrow2 = m2.n_rows;
  int ncol = m1.n_cols;
  
  if (ncol != m2.n_cols) {
    throw std::runtime_error("Incompatible number of dimensions");
  }
  
  mat out(nrow1, nrow2);
  
  for (int r1 = 0; r1 < nrow1; r1++) {
    for (int r2 = 0; r2 < nrow2; r2++) {
      double total = 0;
      for (int c12 = 0; c12 < ncol; c12++) {
        total += pow(m1(r1, c12) - m2(r2, c12), 2);
      }
      out(r1, r2) = sqrt(total);
    }
  }
  
  return out;
}

mat mmult( mat m , mat v, bool  byrow=true )
{
  if( ! m.n_rows == v.n_rows ) stop("Non-conformable arrays") ;
  if( ! m.n_cols == v.n_cols ) stop("Non-conformable arrays") ;
  
  mat out(m) ;
  
  for (int i = 0; i < (m.n_rows); i++) 
  {
    for (int j = 0; j < (m.n_cols); j++) 
    {
      out(i,j)=m(i,j) * v(i,j) ;
    }
  }
  return out ;
}



mat mpow( mat m , double v )
{
  mat out(m) ;
  
  for (int i = 0; i < m.n_rows; i++) 
  {
    for (int j = 0; j < m.n_cols; j++) 
    {
      out(i,j)=pow(m(i,j),v);
    }
  }
  return out ;
}

mat mexp( mat m )
{
  mat out(m) ;
  
  for (int i = 0; i < m.n_rows; i++) 
  {
    for (int j = 0; j < m.n_cols; j++) 
    {
      out(i,j)=exp(m(i,j));
    }
  }
  return out ;
}



mat minv( mat m )
{
  mat out(m) ;
  
  for (int i = 0; i < m.n_rows; i++) 
  {
    for (int j = 0; j < m.n_cols; j++) 
    {
      out(i,j)=(1/m(i,j));
    }
  }
  return out ;
}

mat msum( mat m, double t )
{
  mat out(m) ;
  
  for (int i = 0; i < m.n_rows; i++) 
  {
    for (int j = 0; j < m.n_cols; j++) 
    {
      out(i,j)=m(i,j)+t;
    }
  }
  return out ;
}


mat repsum( mat m, mat t )
{
  mat out(m) ;
  
  for (int i = 0; i < m.n_rows; i++) 
  {
    for (int j = 0; j < m.n_cols; j++) 
    {
      out(i,j)=m(i,j)+t(i,0);
    }
  }
  return out ;
}

mat msqrt( mat m )
{
  mat out(m) ;
  
  for (int i = 0; i < m.n_rows; i++) 
  {
    for (int j = 0; j < m.n_cols; j++) 
    {
      out(i,j)=sqrt(m(i,j));
    }
  }
  return out ;
}
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
vec colSum(mat x) {
  int ncol = x.n_cols;
  vec out(ncol);
  for (int i = 0; i < ncol; i++) {
    out[i] = sum(x.col(i));
  }
  return out;
}



// [[Rcpp::export]]
vec find_lowest_indices(vec vecs, int n = 10) {
  int vecs_size = vecs.size();
  vec indices(vecs_size);
  std::iota(indices.begin(), indices.end(), 0);
  
  std::partial_sort(
    indices.begin(), 
    indices.begin() + n, 
    indices.end(),
    [&](int i, int j) { return vecs[i] < vecs[j]; }
  );
  
  return indices.head(n);
}


arma::mat f_arma_mat_empty( arma::vec dim_mat ) {
  
  // Create an empty matrix with dimensions dim_mat
  return arma::mat(dim_mat[0], dim_mat[1]);
}
// [[Rcpp::export]]
double allsum(mat x) {
  int ncol = x.n_cols;
  vec out(ncol);
  for (int i = 0; i < ncol; i++) {
    out[i] = sum(x.col(i));
  }
  
  double allsum = sum(out);
  return allsum;
}

// [[Rcpp::export]]
mat rbind_mat(mat x, mat y) {
  int xcol = x.n_cols;
  int xrow = x.n_rows;
  int yrow = y.n_rows;
  
  mat out(xrow+yrow, xcol);
  for (int j = 0; j < xcol; j++) {
    for (int i = 0; i < xrow; i++){
      out(i,j) = x(i,j);
    }
    for (int ii = 0; ii < yrow; ii++){
      out(xrow+ii, j) = y(ii,j);
    }
  }
  return out;
}


// [[Rcpp::export]]
vec rbind_vec(vec x, double y) {
  int xcol = x.size();
  
  vec out = zeros(xcol+1);
  for (int i = 0; i < xcol; i++){
    out(i) = x(i);
  }
  out(xcol)=y;
  return out;
}

// [[Rcpp::export]]
mat int_to_mat(int x) {
  mat out(1,1);
  out(0,0)=x;
  return out;
}


// [[Rcpp::export]]
int vecminInd(vec x) {
  // Rcpp supports STL-style iterators
  vec::iterator it = std::min_element(x.begin(), x.end());
  // we want the value so dereference 
  return it - x.begin();
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
vec crossdist2(mat m1, mat m2) {
  int nrow1 = m1.n_rows;
  int nrow2 = m2.n_rows;
  int ncol = m1.n_cols;
  
  if (ncol != m2.n_cols) {
    throw std::runtime_error("Incompatible number of dimensions");
  }
  
  vec out(nrow2);
  
  for (int r1 = 0; r1 < nrow1; r1++) {
    for (int r2 = 0; r2 < nrow2; r2++) {
      double total = 0;
      for (int c12 = 0; c12 < ncol; c12++) {
        total += pow(m1(r1, c12) - m2(r2, c12), 2);
      }
      out[r2] = sqrt(total);
    }
  }
  
  return out;
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List rbf(mat theta, double h=-1){
  mat sqdist = crossdist(theta, theta);
  mat pairdist = mpow(sqdist,2);
  vec paird = vectorise(pairdist);
  if (h<0) {
    h=median(paird);
    h=sqrt(0.5*h/log(theta.n_rows+1));
  }
  
  pairdist = pairdist*(-1);
  pairdist = pairdist / (h*h);
  pairdist = pairdist / 2;
  mat kxy = mexp(pairdist);
  mat dxkxy = (kxy*theta)*(-1);
  vec sumkxy = colSum(kxy);
  for(int i=0 ; i < theta.n_cols; i++){
    dxkxy.col(i)=dxkxy.col(i)+mmult(theta.col(i),sumkxy);
  }
  dxkxy = dxkxy / (h*h);
  return List::create(Named("kxy")=kxy, Named("dxkxy")=dxkxy);
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
mat comp_fixnu_map(mat x0, mat X, vec y, int niter, int imporsize, double stepsize, double fixvalue, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  vec fixnu = rep(fixvalue, theta.n_rows);
  theta.col(theta.n_cols-1)=fixnu;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat zerotheta = zeros(theta.n_cols);
  mat alltheta(theta.n_rows*niter, theta.n_cols);
  int numthe = theta.n_rows;
  
  for(int iter=0; iter<niter; iter++){
    mat loggr(theta.n_cols, theta.n_rows);
    for(int tt=0; tt<theta.n_rows; tt++){
      vec vectheta = vectorise(theta.row(tt));
      mat auxs=rCOMP3_parallel(X, vectheta,imporsize,1);
      for(int ii=0; ii<imporsize; ii++){
        loggr.col(tt)=loggr.col(tt)+loggrah(X,y,fixvalue,vectheta)-loggrah(X,vectorise(auxs.row(ii)), fixvalue, vectheta);
      }
    }
    loggr = loggr/imporsize;
    
    mat lnpgrad = loggr;
    mat gradtheta = lnpgrad.t();
    gradtheta = gradtheta;
    gradtheta = gradtheta / theta.n_rows;
    
    mat pgradtheta = mpow(gradtheta,2);
    
    mat histo = zeros(pgradtheta.n_rows, pgradtheta.n_cols);
    mat fudge = zeros(pgradtheta.n_rows, pgradtheta.n_cols);
    fudge= msum(fudge, fudgefactor);
    if(iter==1){
      histo = histo + pgradtheta;
    } else {
      histo = alpha * histo + (1-alpha) * pgradtheta;
    }
    
    mat adjgrad = gradtheta/(fudge+msqrt(histo));
    
    theta = theta + stepsize * adjgrad;
    theta.col(theta.n_cols-1)=fixnu;
    for(int jjj=0; jjj<numthe; jjj++){
      alltheta.row(jjj+(iter*numthe)) = theta.row(jjj);
    }
  }
  return(alltheta);
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List comp_svi_fixnu_new(mat x0, mat X, vec y, mat map, int niter, int imporsize, int imporsizes, double stepsize, double thres, double fixvalue, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  vec fixnu = rep(fixvalue, theta.n_rows);
  theta.col(theta.n_cols-1)=fixnu;
  map.row(theta.n_cols-1)=fixvalue;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat zerotheta = zeros(theta.n_cols);
  mat alltheta(theta.n_rows*niter, theta.n_cols);
  int numthe = theta.n_rows;
  
  mat keeploggra(imporsize, theta.n_cols);
  vec keeph = zeros(imporsize);
  mat suf(imporsize, theta.n_cols);
  mat thetas = map.t();
  
  mat auxs=rCOMP3_parallel(X, map,imporsize,1);
  for(int line=0; line<imporsize; line++){
    keeploggra.row(line)=loggrah(X,vectorise(auxs.row(line)), fixvalue, map).t();
    keeph(line)=hfunc(Summary(X,vectorise(auxs.row(line))), fixvalue, map);
    suf.row(line)=Summary(X,vectorise(auxs.row(line))).t();
  }
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums = zeros(theta.n_rows);
  for(int iter=0; iter<niter; iter++){
    Rprintf("start %d, \n", iter+1);
    for(int tt=0; tt<theta.n_rows; tt++){
      vec vectheta = vectorise(theta.row(tt));
      double distmap = crossdist2(theta.row(tt), map.t())[0];
      double distmin = min( crossdist2(theta.row(tt), thetas) ); 
      if(distmap > distmin) {
        int minwhich = vecminInd(crossdist2(theta.row(tt), thetas));
        mat selip = zeros(imporsizes);
        mat selip2(imporsizes, theta.n_cols);
        mat norwei = zeros(imporsizes);
        vec rowthetas=vectorise(thetas.row(minwhich));
        for(int jj=0; jj<imporsizes; jj++){
          vec sufs = vectorise(suf.row(imporsize+(minwhich-1)*(imporsizes)+jj));
          selip.row(jj)=hfunc(sufs, fixvalue, vectheta)-keeph(imporsize+(minwhich-1)*(imporsizes)+jj); 
        }
        double mx = selip.max();
        selip= selip-mx;
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsizes; jj++){
          selip2.row(jj)=(expselip.row(jj)/selsum.t())*keeploggra.row(imporsize+(minwhich-1)*(imporsizes)+jj);
          mat nor = (expselip.row(jj)/selsum.t());
          norwei.row(jj)= nor;
        }
        vec selsum2 = colSum(selip2);
        norwei = mpow(norwei,2); 
        vec norweisum = colSum(norwei);
        norweisums(tt)=norweisum(0);
        parsum.row(tt)=loggrah(X,y,fixvalue,vectheta).t()-selsum2.t();
        ess = norweisums;
        mat midess=minv(ess);
        midess = midess.t();
        
        if(min(midess(tt))<(imporsizes/thres)){
          thetas = rbind_mat(thetas, theta.row(tt));
          parsum.row(tt)=zerotheta.t();
          mat gibs=rCOMP3_parallel(X, vectheta,imporsizes,1);
          mat savelog = zeros(theta.n_cols);
          for(int impors=0; impors<imporsizes; impors++){
            suf = rbind_mat(suf, Summary(X, vectorise(gibs.row(impors))).t());
            keeploggra = rbind_mat(keeploggra, loggrah(X, vectorise(gibs.row(impors)), fixvalue, vectheta).t());
            keeph = rbind_vec(keeph, hfunc(Summary(X, vectorise(gibs.row(impors))), fixvalue, vectheta));
            savelog = savelog + loggrah(X, vectorise(gibs.row(impors)), fixvalue, vectheta);
          }
          savelog = savelog/imporsizes;
          parsum.row(tt)=parsum.row(tt)+savelog.t();
          parsum.row(tt)=loggrah(X,y,fixvalue,vectheta).t()-parsum.row(tt);
        }
      }
      
      
      else{
        mat selip = zeros(imporsize);
        mat selip2(imporsize, theta.n_cols);
        mat norwei = zeros(imporsize);
        for(int jj=0; jj<imporsize; jj++){
          selip.row(jj)=hfunc(vectorise(suf.row(jj)), fixvalue, vectheta)-keeph(jj);
        }
        double mx = selip.max();
        selip= selip-mx;
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsize; jj++){
          selip2.row(jj)=(expselip.row(jj)/selsum.t())*keeploggra.row(jj);
          mat nor = (expselip.row(jj)/selsum.t());
          norwei.row(jj)= nor;
        }
        vec selsum2 = colSum(selip2);
        norwei = mpow(norwei,2); 
        vec norweisum = colSum(norwei);
        norweisums(tt)=norweisum(0);
        parsum.row(tt)=loggrah(X,y,fixvalue,vectheta).t()-selsum2.t();
        ess = norweisums;
        mat midess=minv(ess);
        midess = midess.t();
        if(min(midess(tt))<(imporsize/thres)){
          thetas = rbind_mat(thetas, theta.row(tt));
          parsum.row(tt)=zerotheta.t();
          mat gibs=rCOMP3_parallel(X, vectheta,imporsizes,1);
          mat savelog = zeros(theta.n_cols);
          for(int impors=0; impors<imporsizes; impors++){
            suf = rbind_mat(suf, Summary(X, vectorise(gibs.row(impors))).t());
            keeploggra = rbind_mat(keeploggra, loggrah(X, vectorise(gibs.row(impors)), fixvalue, vectheta).t());
            keeph = rbind_vec(keeph, hfunc(Summary(X, vectorise(gibs.row(impors))), fixvalue, vectheta));
            savelog = savelog + loggrah(X, vectorise(gibs.row(impors)), fixvalue, vectheta);
          }
          savelog = savelog/imporsizes;
          parsum.row(tt)=parsum.row(tt)+savelog.t();
          parsum.row(tt)=loggrah(X,y,fixvalue,vectheta).t()-parsum.row(tt);
        }
      }
    }  
    
    mat lnpgrad = parsum.t();
    List ker=rbf(theta, -1);
    mat kxy = ker["kxy"];
    mat dxkxy = ker["dxkxy"];
    mat gradtheta = kxy*lnpgrad.t();
    gradtheta = gradtheta + dxkxy;
    gradtheta = gradtheta / theta.n_rows;
    
    mat pgradtheta = mpow(gradtheta,2);
    
    mat histo = zeros(pgradtheta.n_rows, pgradtheta.n_cols);
    mat fudge = zeros(pgradtheta.n_rows, pgradtheta.n_cols);
    fudge= msum(fudge, fudgefactor);
    if(iter==1){
      histo = histo + pgradtheta;
    } else {
      histo = alpha * histo + (1-alpha) * pgradtheta;
    }
    
    mat adjgrad = gradtheta/(fudge+msqrt(histo));
    
    theta = theta + stepsize * adjgrad;
    theta.col(theta.n_cols-1)=fixnu;
    for(int jjj=0; jjj<numthe; jjj++){
      alltheta.row(jjj+(iter*numthe)) = theta.row(jjj);
    }
  }
  return List::create(Named("theta")=theta, Named("alltheta")=alltheta);
}








// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List comp_svi_fixnu_parallel(mat x0, mat X, vec y, mat map, int niter, int imporsize, int imporsizes, double stepsize, double thres, double fixvalue, int core, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  vec fixnu = rep(fixvalue, theta.n_rows);
  theta.col(theta.n_cols-1)=fixnu;
  map.row(theta.n_cols-1)=fixvalue;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat zerotheta = zeros(theta.n_cols);
  mat alltheta(theta.n_rows*niter, theta.n_cols);
  int numthe = theta.n_rows;
  
  mat keeploggra(imporsize, theta.n_cols);
  vec keeph = zeros(imporsize);
  mat suf(imporsize, theta.n_cols);
  mat thetas = map.t();
  vec summary = Summary(X,y);
  for(int line=0; line<imporsize; line++){
    mat summa=rCOMP_one_suf(X, map);
    keeploggra.row(line)=loggrah3(X, vectorise(summa), fixvalue, map).t();
    keeph(line)=hfunc(vectorise(summa), fixvalue, map);
    suf.row(line)=summa;
  }
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums = zeros(theta.n_rows);
  for(int iter=0; iter<niter; iter++){
    Rprintf("start %d, \n", iter+1);
    for(int tt=0; tt<theta.n_rows; tt++){
      vec vectheta = vectorise(theta.row(tt));
      double distmap = crossdist2(theta.row(tt), map.t())[0];
      double distmin = min( crossdist2(theta.row(tt), thetas) ); 
      if(distmap > distmin) {
        int minwhich = vecminInd(crossdist2(theta.row(tt), thetas));
        mat selip = zeros(imporsizes);
        mat selip2(imporsizes, theta.n_cols);
        mat norwei = zeros(imporsizes);
        vec rowthetas=vectorise(thetas.row(minwhich));
        for(int jj=0; jj<imporsizes; jj++){
          vec sufs = vectorise(suf.row(imporsize+(minwhich-1)*(imporsizes)+jj));
          selip.row(jj)=hfunc(sufs, fixvalue, vectheta)-keeph(imporsize+(minwhich-1)*(imporsizes)+jj); 
        }
        double mx = selip.max();
        selip= selip-mx;
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsizes; jj++){
          selip2.row(jj)=(expselip.row(jj)/selsum.t())*keeploggra.row(imporsize+(minwhich-1)*(imporsizes)+jj);
          mat nor = (expselip.row(jj)/selsum.t());
          norwei.row(jj)= nor;
        }
        vec selsum2 = colSum(selip2);
        norwei = mpow(norwei,2); 
        vec norweisum = colSum(norwei);
        norweisums(tt)=norweisum(0);
        parsum.row(tt)=loggrah3(X,summary,fixvalue,vectheta).t()-selsum2.t();
        ess = norweisums;
        mat midess=minv(ess);
        midess = midess.t();
        
        if(min(midess(tt))<(imporsizes/thres)){
          thetas = rbind_mat(thetas, theta.row(tt));
          parsum.row(tt)=zerotheta.t();
          mat savelog = zeros(theta.n_cols);
          for(int impors=0; impors<imporsizes; impors++){
            mat gibs=rCOMP_one_suf(X, vectheta);
            suf = rbind_mat(suf, gibs);
            keeploggra = rbind_mat(keeploggra, loggrah3(X, vectorise(gibs), fixvalue, vectheta).t());
            keeph = rbind_vec(keeph, hfunc(vectorise(gibs), fixvalue, vectheta));
            savelog = savelog + loggrah3(X, vectorise(gibs), fixvalue, vectheta);
          }
          savelog = savelog/imporsizes;
          parsum.row(tt)=parsum.row(tt)+savelog.t();
          parsum.row(tt)=loggrah3(X,summary,fixvalue,vectheta).t()-parsum.row(tt);
        }
      }
      
      
      else{
        mat selip = zeros(imporsize);
        mat selip2(imporsize, theta.n_cols);
        mat norwei = zeros(imporsize);
        for(int jj=0; jj<imporsize; jj++){
          selip.row(jj)=hfunc(vectorise(suf.row(jj)), fixvalue, vectheta)-keeph(jj);
        }
        double mx = selip.max();
        selip= selip-mx;
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsize; jj++){
          selip2.row(jj)=(expselip.row(jj)/selsum.t())*keeploggra.row(jj);
          mat nor = (expselip.row(jj)/selsum.t());
          norwei.row(jj)= nor;
        }
        vec selsum2 = colSum(selip2);
        norwei = mpow(norwei,2); 
        vec norweisum = colSum(norwei);
        norweisums(tt)=norweisum(0);
        parsum.row(tt)=loggrah3(X,summary,fixvalue,vectheta).t()-selsum2.t();
        ess = norweisums;
        mat midess=minv(ess);
        midess = midess.t();
        if(min(midess(tt))<(imporsize/thres)){
          thetas = rbind_mat(thetas, theta.row(tt));
          parsum.row(tt)=zerotheta.t();
          mat savelog = zeros(theta.n_cols);
          for(int impors=0; impors<imporsizes; impors++){
            mat gibs=rCOMP_one_suf(X, vectheta);
            suf = rbind_mat(suf, gibs);
            keeploggra = rbind_mat(keeploggra, loggrah3(X, vectorise(gibs), fixvalue, vectheta).t());
            keeph = rbind_vec(keeph, hfunc(vectorise(gibs), fixvalue, vectheta));
            savelog = savelog + loggrah3(X, vectorise(gibs), fixvalue, vectheta);
          }
          savelog = savelog/imporsizes;
          parsum.row(tt)=parsum.row(tt)+savelog.t();
          parsum.row(tt)=loggrah3(X,summary,fixvalue,vectheta).t()-parsum.row(tt);
        }
      }
    } 
    
    
    mat lnpgrad = parsum.t();
    List ker=rbf(theta, -1);
    mat kxy = ker["kxy"];
    mat dxkxy = ker["dxkxy"];
    mat gradtheta = kxy*lnpgrad.t();
    gradtheta = gradtheta + dxkxy;
    gradtheta = gradtheta / theta.n_rows;
    
    mat pgradtheta = mpow(gradtheta,2);
    
    mat histo = zeros(pgradtheta.n_rows, pgradtheta.n_cols);
    mat fudge = zeros(pgradtheta.n_rows, pgradtheta.n_cols);
    fudge= msum(fudge, fudgefactor);
    if(iter==1){
      histo = histo + pgradtheta;
    } else {
      histo = alpha * histo + (1-alpha) * pgradtheta;
    }
    
    mat adjgrad = gradtheta/(fudge+msqrt(histo));
    
    theta = theta + stepsize * adjgrad;
    theta.col(theta.n_cols-1)=fixnu;
    for(int jjj=0; jjj<numthe; jjj++){
      alltheta.row(jjj+(iter*numthe)) = theta.row(jjj);
    }
  }
  return List::create(Named("theta")=theta, Named("alltheta")=alltheta);
}









// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List comp_svi_fixnu_par(mat x0, mat X, vec y, mat map, int niter, int imporsize, int imporsizes, double stepsize, double thres, double fixvalue, int core, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  vec fixnu = rep(fixvalue, theta.n_rows);
  theta.col(theta.n_cols-1)=fixnu;
  map.row(theta.n_cols-1)=fixvalue;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat zerotheta = zeros(theta.n_cols);
  int numthe = theta.n_rows;
  vec mc = zeros(numthe);
  int percore = theta.n_rows/core;
  vec nummc = zeros(niter);
  
  mat keeploggra(imporsize, theta.n_cols);
  vec keeph = zeros(imporsize);
  mat suf(imporsize, theta.n_cols);
  mat thetas = map.t();
  vec summary = Summary(X,y);
  
  for(int line=0; line<imporsize; line++){
    mat summa=rCOMP_one_suf(X, map);
    keeploggra.row(line)=loggrah3(X, vectorise(summa), fixvalue, map).t();
    keeph(line)=hfunc(vectorise(summa), fixvalue, map);
    suf.row(line)=summa;
  }
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums = zeros(theta.n_rows);
  for(int iter=0; iter<niter; iter++){
    mc = zeros(theta.n_rows);
    for(int tt=0; tt<theta.n_rows; tt++){
      vec vectheta = vectorise(theta.row(tt));
      double distmap = crossdist2(theta.row(tt), map.t())[0];
      double distmin = min( crossdist2(theta.row(tt), thetas) ); 
      if(distmap > distmin) {
        int minwhich = vecminInd(crossdist2(theta.row(tt), thetas));
        mat selip = zeros(imporsizes);
        mat selip2(imporsizes, theta.n_cols);
        mat norwei = zeros(imporsizes);
        vec rowthetas=vectorise(thetas.row(minwhich));
        for(int jj=0; jj<imporsizes; jj++){
          vec sufs = vectorise(suf.row(imporsize+(minwhich-1)*(imporsizes)+jj));
          selip.row(jj)=hfunc(sufs, fixvalue, vectheta)-keeph(imporsize+(minwhich-1)*(imporsizes)+jj); 
        }
        double mx = selip.max();
        selip= selip-mx;
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsizes; jj++){
          selip2.row(jj)=(expselip.row(jj)/selsum.t())*keeploggra.row(imporsize+(minwhich-1)*(imporsizes)+jj);
          mat nor = (expselip.row(jj)/selsum.t());
          norwei.row(jj)= nor;
        }
        vec selsum2 = colSum(selip2);
        norwei = mpow(norwei,2); 
        vec norweisum = colSum(norwei);
        norweisums(tt)=norweisum(0);
        parsum.row(tt)=loggrah3(X,summary,fixvalue,vectheta).t()-selsum2.t();
        ess = norweisums;
        mat midess=minv(ess);
        midess = midess.t();
        
        if(min(midess(tt))<(imporsizes/thres)){
          mc[tt]=1;
        }
      }
      
      
      else{
        mat selip = zeros(imporsize);
        mat selip2(imporsize, theta.n_cols);
        mat norwei = zeros(imporsize);
        for(int jj=0; jj<imporsize; jj++){
          selip.row(jj)=hfunc(vectorise(suf.row(jj)), fixvalue, vectheta)-keeph(jj);
        }
        double mx = selip.max();
        selip= selip-mx;
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsize; jj++){
          selip2.row(jj)=(expselip.row(jj)/selsum.t())*keeploggra.row(jj);
          mat nor = (expselip.row(jj)/selsum.t());
          norwei.row(jj)= nor;
        }
        vec selsum2 = colSum(selip2);
        norwei = mpow(norwei,2); 
        vec norweisum = colSum(norwei);
        norweisums(tt)=norweisum(0);
        parsum.row(tt)=loggrah3(X,summary,fixvalue,vectheta).t()-selsum2.t();
        ess = norweisums;
        mat midess=minv(ess);
        midess = midess.t();
        
        if(min(midess(tt))<(imporsize/thres)){
          mc[tt]=1;
        }
        
      }
    }  
    
    mat suf2(theta.n_rows*imporsizes, theta.n_cols);
  
    
    int cor;
#pragma omp parallel shared(theta) private(cor)
{
#pragma omp for schedule(static)
  for(cor=0; cor<core; cor++){
    for(int tt2=0; tt2<percore; tt2++){
      if(mc[percore*cor+tt2]==1){
        vec vectheta = vectorise(theta.row(percore*cor+tt2));
        parsum.row(percore*cor+tt2)=zerotheta.t();
        mat savelog = zeros(theta.n_cols);
        for(int impors=0; impors<imporsizes; impors++){
          mat gibs=rCOMP_one_suf(X, vectheta);
          suf2.row(((percore*cor+tt2)*imporsizes)+impors)=gibs;
          savelog = savelog + loggrah3(X, vectorise(gibs), fixvalue, vectheta);
        }
        savelog = savelog/imporsizes;
        parsum.row(percore*cor+tt2)=parsum.row(percore*cor+tt2)+savelog.t();
        parsum.row(percore*cor+tt2)=loggrah3(X,summary,fixvalue,vectheta).t()-parsum.row(percore*cor+tt2);
        
      }
    }
  }
} 

for(int ttt=0; ttt<theta.n_rows; ttt++){    
  if(mc[ttt]==1){
    thetas = rbind_mat(thetas, theta.row(ttt));
    vec vectheta = vectorise(theta.row(ttt));
    for(int iii=0; iii<imporsizes; iii++){
      suf = rbind_mat(suf, suf2.row(ttt*imporsizes+iii));
      keeploggra = rbind_mat(keeploggra, loggrah3(X, vectorise(suf2.row(ttt*imporsizes+iii)), fixvalue, vectheta).t());
      keeph = rbind_vec(keeph, hfunc(vectorise(suf2.row(ttt*imporsizes+iii)), fixvalue, vectheta));
    }
  }
}

  
nummc[iter] = sum(mc);
mat lnpgrad = parsum.t();
List ker=rbf(theta, -1);
mat kxy = ker["kxy"];
mat dxkxy = ker["dxkxy"];
mat gradtheta = kxy*lnpgrad.t();
gradtheta = gradtheta + dxkxy;
gradtheta = gradtheta / theta.n_rows;

mat pgradtheta = mpow(gradtheta,2);

mat histo = zeros(pgradtheta.n_rows, pgradtheta.n_cols);
mat fudge = zeros(pgradtheta.n_rows, pgradtheta.n_cols);
fudge= msum(fudge, fudgefactor);
if(iter==1){
  histo = histo + pgradtheta;
} else {
  histo = alpha * histo + (1-alpha) * pgradtheta;
}

mat adjgrad = gradtheta/(fudge+msqrt(histo));

theta = theta + stepsize * adjgrad;
theta.col(theta.n_cols-1)=fixnu;

  }
  return List::create(Named("theta")=theta, Named("nummc")=nummc);
}




















// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List comp_svi_fixnu_par2(mat x0, mat X, vec y, mat map, int niter, int imporsize, int imporsizes, double stepsize, double thres, double fixvalue, int core, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  vec fixnu = rep(fixvalue, theta.n_rows);
  theta.col(theta.n_cols-1)=fixnu;
  map.row(theta.n_cols-1)=fixvalue;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat zerotheta = zeros(theta.n_cols);
  int numthe = theta.n_rows;
  vec mc = zeros(numthe);
  int percore = theta.n_rows/core;
  vec nummc = zeros(niter);
  
  mat keeploggra(imporsize, theta.n_cols);
  vec keeph = zeros(imporsize);
  mat suf(imporsize, theta.n_cols);
  mat thetas = map.t();
  vec summary = Summary(X,y);
  
  for(int line=0; line<imporsize; line++){
    mat summa=rCOMP_one_suf(X, map);
    keeploggra.row(line)=loggrah3(X, vectorise(summa), fixvalue, map).t();
    keeph(line)=hfunc(vectorise(summa), fixvalue, map);
    suf.row(line)=summa;
  }
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums = zeros(theta.n_rows);
  for(int iter=0; iter<niter; iter++){
    mc = zeros(theta.n_rows);
    for(int tt=0; tt<theta.n_rows; tt++){
      vec vectheta = vectorise(theta.row(tt));
      double distmap = crossdist2(theta.row(tt), map.t())[0];
      double distmin = min( crossdist2(theta.row(tt), thetas) ); 
      if(distmap > distmin) {
        int minwhich = vecminInd(crossdist2(theta.row(tt), thetas));
        mat selip = zeros(imporsizes);
        mat selip2(imporsizes, theta.n_cols);
        mat norwei = zeros(imporsizes);
        vec rowthetas=vectorise(thetas.row(minwhich));
        for(int jj=0; jj<imporsizes; jj++){
          vec sufs = vectorise(suf.row(imporsize+(minwhich-1)*(imporsizes)+jj));
          selip.row(jj)=hfunc(sufs, fixvalue, vectheta)-keeph(imporsize+(minwhich-1)*(imporsizes)+jj); 
        }
        double mx = selip.max();
        selip= selip-mx;
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsizes; jj++){
          selip2.row(jj)=(expselip.row(jj)/selsum.t())*keeploggra.row(imporsize+(minwhich-1)*(imporsizes)+jj);
          mat nor = (expselip.row(jj)/selsum.t());
          norwei.row(jj)= nor;
        }
        vec selsum2 = colSum(selip2);
        norwei = mpow(norwei,2); 
        vec norweisum = colSum(norwei);
        norweisums(tt)=norweisum(0);
        parsum.row(tt)=loggrah3(X,summary,fixvalue,vectheta).t()-selsum2.t();
        ess = norweisums;
        mat midess=minv(ess);
        midess = midess.t();
        
        if(min(midess(tt))<(imporsizes/thres)){
          mc[tt]=1;
        }
      }
      
      
      else{
        mat selip = zeros(imporsize);
        mat selip2(imporsize, theta.n_cols);
        mat norwei = zeros(imporsize);
        for(int jj=0; jj<imporsize; jj++){
          selip.row(jj)=hfunc(vectorise(suf.row(jj)), fixvalue, vectheta)-keeph(jj);
        }
        double mx = selip.max();
        selip= selip-mx;
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsize; jj++){
          selip2.row(jj)=(expselip.row(jj)/selsum.t())*keeploggra.row(jj);
          mat nor = (expselip.row(jj)/selsum.t());
          norwei.row(jj)= nor;
        }
        vec selsum2 = colSum(selip2);
        norwei = mpow(norwei,2); 
        vec norweisum = colSum(norwei);
        norweisums(tt)=norweisum(0);
        parsum.row(tt)=loggrah3(X,summary,fixvalue,vectheta).t()-selsum2.t();
        ess = norweisums;
        mat midess=minv(ess);
        midess = midess.t();
        
        if(min(midess(tt))<(imporsize/thres)){
          mc[tt]=1;
        }
        
      }
    }  
    
    mat suf2(theta.n_rows*imporsizes, theta.n_cols);
    
    int cor;
    
    std::vector<int> mc_indices;
    for (int i = 0; i < mc.size(); ++i) {
      if (mc[i] == 1) {
        mc_indices.push_back(i);  
      }
    }
    
    int total_tasks = mc_indices.size();
    int percore = (total_tasks + core - 1) / core;
    
#pragma omp parallel shared(theta, mc_indices) private(cor)
{
#pragma omp for schedule(static)
  for (cor = 0; cor < core; ++cor) {
    for (int tt2 = 0; tt2 < percore; ++tt2) {
      int task_idx = cor * percore + tt2;
      if (task_idx < total_tasks) {
        int idx = mc_indices[task_idx];  
        

        vec vectheta = vectorise(theta.row(idx));
        parsum.row(idx) = zerotheta.t();  
        mat savelog = zeros(theta.n_cols);  
        

        for (int impors = 0; impors < imporsizes; ++impors) {
          mat gibs = rCOMP_one_suf(X, vectheta);
          suf2.row(idx * imporsizes + impors) = gibs;
          savelog += loggrah3(X, vectorise(gibs), fixvalue, vectheta);  
        }
        
       
        savelog /= imporsizes;
        parsum.row(idx) += savelog.t();
        parsum.row(idx) = loggrah3(X, summary, fixvalue, vectheta).t() - parsum.row(idx); 
      }
    }
  }
}
  
 

for(int ttt=0; ttt<theta.n_rows; ttt++){    
  if(mc[ttt]==1){
    thetas = rbind_mat(thetas, theta.row(ttt));
    vec vectheta = vectorise(theta.row(ttt));
    for(int iii=0; iii<imporsizes; iii++){
      suf = rbind_mat(suf, suf2.row(ttt*imporsizes+iii));
      keeploggra = rbind_mat(keeploggra, loggrah3(X, vectorise(suf2.row(ttt*imporsizes+iii)), fixvalue, vectheta).t());
      keeph = rbind_vec(keeph, hfunc(vectorise(suf2.row(ttt*imporsizes+iii)), fixvalue, vectheta));
    }
  }
}


nummc[iter] = sum(mc);
mat lnpgrad = parsum.t();
List ker=rbf(theta, -1);
mat kxy = ker["kxy"];
mat dxkxy = ker["dxkxy"];
mat gradtheta = kxy*lnpgrad.t();
gradtheta = gradtheta + dxkxy;
gradtheta = gradtheta / theta.n_rows;

mat pgradtheta = mpow(gradtheta,2);

mat histo = zeros(pgradtheta.n_rows, pgradtheta.n_cols);
mat fudge = zeros(pgradtheta.n_rows, pgradtheta.n_cols);
fudge= msum(fudge, fudgefactor);
if(iter==1){
  histo = histo + pgradtheta;
} else {
  histo = alpha * histo + (1-alpha) * pgradtheta;
}

mat adjgrad = gradtheta/(fudge+msqrt(histo));

theta = theta + stepsize * adjgrad;
theta.col(theta.n_cols-1)=fixnu;

  }
  return List::create(Named("theta")=theta, Named("nummc")=nummc);
}





















// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List comp_weight_mcsvgd(mat x0, mat X, vec y, mat map, int niter, int imporsize, double close, double stepsize, double thres, double thres2,double fixvalue, int core, int prerun, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  vec fixnu = rep(fixvalue, theta.n_rows);
  theta.col(theta.n_cols-1)=fixnu;
  map.row(theta.n_cols-1)=fixvalue;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat zerotheta = zeros(theta.n_cols);
  int numthe = theta.n_rows;
  vec mc = zeros(numthe);
  int percore = theta.n_rows/core;
  vec nummc = zeros(niter);
  
  mat keeploggra(imporsize, theta.n_cols);
  vec keeph = zeros(imporsize);
  mat suf(imporsize, theta.n_cols);
  mat thetas = map.t();
  vec summary = Summary(X,y);
  
  for(int line=0; line<imporsize; line++){
    mat summa=rCOMP_one_suf(X, map);
    keeploggra.row(line)=loggrah3(X, vectorise(summa), fixvalue, map).t();
    keeph(line)=hfunc(vectorise(summa), fixvalue, map);
    suf.row(line)=summa;
  }
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums = zeros(theta.n_rows);
  for(int iter=0; iter<niter; iter++){
    mc = zeros(theta.n_rows);
    if(thetas.n_rows>prerun){
      for(int tt=0; tt<theta.n_rows; tt++){
        vec vectheta = vectorise(theta.row(tt));
        vec distval = crossdist2(theta.row(tt), thetas);
        vec minwhich = find_lowest_indices(distval, close);
        mat selip = zeros(imporsize*close);
        mat selip2(imporsize*close, theta.n_cols);
        mat norwei = zeros(imporsize*close);
        for(int near=0; near<close; near++){
          for(int jj=0; jj<imporsize; jj++){
            vec sufs = vectorise(suf.row((minwhich[near])*(imporsize)+jj));
            selip.row(imporsize*near+jj)=hfunc(sufs, fixvalue, vectheta)-keeph((minwhich[near])*(imporsize)+jj); 
          }
        }
        double mx = selip.max();
        selip= selip-mx;
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int near=0; near<close; near++){
          for(int jj=0; jj<imporsize; jj++){
            selip2.row(imporsize*near+jj)=(expselip.row(imporsize*near+jj)/selsum.t())*keeploggra.row((minwhich[near])*(imporsize)+jj);
            mat nor = (expselip.row(imporsize*near+jj)/selsum.t());
            norwei.row(imporsize*near+jj)= nor;
          }
        }
        vec selsum2 = colSum(selip2);
        norwei = mpow(norwei,2); 
        vec norweisum = colSum(norwei);
        norweisums(tt)=norweisum(0);
        parsum.row(tt)=loggrah3(X,summary,fixvalue,vectheta).t()-selsum2.t();
        ess = norweisums;
        mat midess=minv(ess);
        midess = midess.t();
        
        if(min(midess(tt))<(thres)){
          mc[tt]=1;
        }
      }
    }else{
      for(int tt=0; tt<theta.n_rows; tt++){
      vec vectheta = vectorise(theta.row(tt));
      int minwhich = vecminInd(crossdist2(theta.row(tt), thetas));
        mat selip = zeros(imporsize);
        mat selip2(imporsize, theta.n_cols);
        mat norwei = zeros(imporsize);
        for(int jj=0; jj<imporsize; jj++){
          selip.row(jj)=hfunc(vectorise(suf.row(jj)), fixvalue, vectheta)-keeph(jj);
        }
        double mx = selip.max();
        selip= selip-mx;
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsize; jj++){
          selip2.row(jj)=(expselip.row(jj)/selsum.t())*keeploggra.row(jj);
          mat nor = (expselip.row(jj)/selsum.t());
          norwei.row(jj)= nor;
        }
        vec selsum2 = colSum(selip2);
        norwei = mpow(norwei,2); 
        vec norweisum = colSum(norwei);
        norweisums(tt)=norweisum(0);
        parsum.row(tt)=loggrah3(X,summary,fixvalue,vectheta).t()-selsum2.t();
        ess = norweisums;
        mat midess=minv(ess);
        midess = midess.t();
        
        if(min(midess(tt))<(thres2)){
          mc[tt]=1;
        }
      } 
      }
      
    
    mat suf2(theta.n_rows*imporsize, theta.n_cols);
    
    
    int cor;
#pragma omp parallel shared(theta) private(cor)
{
#pragma omp for schedule(static)
  for(cor=0; cor<core; cor++){
    for(int tt2=0; tt2<percore; tt2++){
      if(mc[percore*cor+tt2]==1){
        vec vectheta = vectorise(theta.row(percore*cor+tt2));
        parsum.row(percore*cor+tt2)=zerotheta.t();
        mat savelog = zeros(theta.n_cols);
        for(int impors=0; impors<imporsize; impors++){
          mat gibs=rCOMP_one_suf(X, vectheta);
          suf2.row(((percore*cor+tt2)*imporsize)+impors)=gibs;
          savelog = savelog + loggrah3(X, vectorise(gibs), fixvalue, vectheta);
        }
        savelog = savelog/imporsize;
        parsum.row(percore*cor+tt2)=parsum.row(percore*cor+tt2)+savelog.t();
        parsum.row(percore*cor+tt2)=loggrah3(X,summary,fixvalue,vectheta).t()-parsum.row(percore*cor+tt2);
        
      }
    }
  }
} 

for(int ttt=0; ttt<theta.n_rows; ttt++){    
  if(mc[ttt]==1){
    thetas = rbind_mat(thetas, theta.row(ttt));
    vec vectheta = vectorise(theta.row(ttt));
    for(int iii=0; iii<imporsize; iii++){
      suf = rbind_mat(suf, suf2.row(ttt*imporsize+iii));
      keeploggra = rbind_mat(keeploggra, loggrah3(X, vectorise(suf2.row(ttt*imporsize+iii)), fixvalue, vectheta).t());
      keeph = rbind_vec(keeph, hfunc(vectorise(suf2.row(ttt*imporsize+iii)), fixvalue, vectheta));
    }
  }
}

nummc[iter] = sum(mc);
mat lnpgrad = parsum.t();
List ker=rbf(theta, -1);
mat kxy = ker["kxy"];
mat dxkxy = ker["dxkxy"];
mat gradtheta = kxy*lnpgrad.t();
gradtheta = gradtheta + dxkxy;
gradtheta = gradtheta / theta.n_rows;

mat pgradtheta = mpow(gradtheta,2);

mat histo = zeros(pgradtheta.n_rows, pgradtheta.n_cols);
mat fudge = zeros(pgradtheta.n_rows, pgradtheta.n_cols);
fudge= msum(fudge, fudgefactor);
if(iter==1){
  histo = histo + pgradtheta;
} else {
  histo = alpha * histo + (1-alpha) * pgradtheta;
}

mat adjgrad = gradtheta/(fudge+msqrt(histo));

theta = theta + stepsize * adjgrad;
theta.col(theta.n_cols-1)=fixnu;

  }
  return List::create(Named("theta")=theta, Named("nummc")=nummc);
}



