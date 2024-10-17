// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-


// we only include RcppArmadillo.h which pulls Rcpp.h in for us
// [[Rcpp::depends("RcppArmadillo")]]


#include <RcppArmadillo.h>
#include <limits>
#include <omp.h>
#include <Rcpp.h>
#include <cmath>
// #define min(x,y) (((x) < (y)) ? (x) : (y))

using namespace std;
using namespace Rcpp;
using namespace arma;




// // [[Rcpp::export]]
// mat potts_stat(RawVector foo, vec theta, int cycle){
//   
//   // Obtaining namespace of potts package
//   Environment pkg = Environment::namespace_env("potts");
//   
//   // Picking up calc_t_full() function from potts package
//   Function f = pkg["potts"];
//   
//   List result = f(foo, theta, cycle);
//   mat stat = as<arma::mat>( result[9] );
//   return stat.row(cycle - 1);
// }
// 
// 
// 
// // [[Rcpp::export]]
// mat potts_stat2(RawVector foo, vec theta, int cycle){
//   
//   // Obtaining namespace of potts package
//   Environment pkg = Environment::namespace_env("potts");
//   
//   // Picking up calc_t_full() function from potts package
//   Function f = pkg["potts"];
//   
//   List result = f(foo, theta, cycle);
//   return as<arma::mat>( result[9] );
// }
// 
// 
// 
// // [[Rcpp::export]]
// List pottsDMH(RawVector foo, mat stat, mat COV, vec theta, int outer, int cycle, 
//               bool updateCOV, double sigma2, int adaptInterval, double adaptFactorExponent, int adapIter){
// 
//   double logprob, u;
//   int p = theta.size();
//   vec thetaprev = zeros(p), thetaprop = zeros(p);
//   mat postSamples = zeros(outer, p), statprop = zeros(1, p);
//   double rhat, gamma1, gamma2, c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
//   vec accprob = zeros(outer);
//   mat cholCOV = trans( chol( sigma2 * ( COV + 0.00001 * diagmat(ones(p)) ) ) );
//   
// 
//   for(int l = 0; l< outer; l++){
//     
//     if(updateCOV){
//       if( (l+1 >= adaptInterval) && (l+1 - (adaptInterval * trunc((l+1) / adaptInterval)) == 0) ){
//         rhat = sum( accprob.rows(l+1-adaptInterval, l-1) ) / (adaptInterval-1);
//         gamma1 = 1 / pow(adapIter, c1);
//         gamma2 = c0 * gamma1;
//         sigma2 = exp( log(sigma2) + gamma2 * (rhat - ropt) );
//         
//         COV = COV + gamma1 * ( cov( postSamples.rows(l+1-adaptInterval, l-1) ) - COV );
//         cholCOV = trans( chol( sigma2 * ( COV + 0.00001 * diagmat(ones(p)) ) ) );
//         adapIter = adapIter + 1;
//       }
//     }
//     thetaprev = theta;
//     thetaprop = thetaprev + trans( cholCOV ) * randn(p);
//     
//     if(thetaprop[p-1] > 0){
//       statprop = potts_stat(foo, thetaprop, cycle);
//       
//       logprob = sum( -0.05 * trans(thetaprop) * thetaprop + 0.05 * trans(thetaprev) * thetaprev + 
//         (statprop - stat) * (thetaprev - thetaprop) );
//       
//       u = log( randu() );
//       if( u < logprob ){
//         theta = thetaprop;
//         accprob[l] = 1;
//       } 
//     }
//     postSamples.row(l) = trans(theta);
//     
//     if ( (l+1) % 100 == 0 ) {
//       Rprintf("Generated %d samples...\n", l+1);
//     } 
//   }
// 
//   return Rcpp::List::create(Rcpp::Named("postSamples") = postSamples,
//                             Rcpp::Named("accprob") = accprob,
//                             Rcpp::Named("adapIter") = adapIter,
//                             Rcpp::Named("sigma2") = sigma2,
//                             Rcpp::Named("COV") = COV);
// }




// [[Rcpp::export]]
double potts_stat(RawVector foo, int ncolor, double beta, int cycle){
  
  // Obtaining namespace of potts package
  Environment pkg = Environment::namespace_env("potts");
  
  // Picking up calc_t_full() function from potts package
  Function f = pkg["potts"];
  
  vec theta = zeros(ncolor+1);
  theta[ncolor] = beta;
  List result = f(foo, theta, cycle);
  NumericMatrix stat = result[9];
  return stat(cycle - 1, ncolor);
}



// [[Rcpp::export]]
vec potts_stat2(RawVector foo, int ncolor, double beta, int cycle){
  
  // Obtaining namespace of potts package
  Environment pkg = Environment::namespace_env("potts");
  
  // Picking up calc_t_full() function from potts package
  Function f = pkg["potts"];
  
  vec theta = zeros(ncolor+1);
  theta[ncolor] = beta;
  List result = f(foo, theta, cycle);
  mat stat = as<arma::mat>( result[9] );
  return stat.col(ncolor);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Summary statistics for potts model
int Summary(mat X){
  int nrow = X.n_rows, ncol = X.n_cols, s1 = 0, s2 = 0;
  int result;
  
  for(int j = 0; j< ncol; j++){
    for(int i = 0; i< nrow-1; i++){
      int diff = X(i,j)-X(i+1,j);
      if(diff == 0){
        s1 = s1 +1;
      }
      
    }
  }
  
  for(int i = 0; i< nrow; i++){
    for(int j = 0; j< ncol-1; j++){
      int diff = X(i,j)-X(i,j+1);
      if(diff == 0){
        s2 = s2 +1;
      }
      
    }
  }
  
  result = s1 + s2;
  
  return(result);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Summary statistics for potts model
int Summary2(mat X){
  int nrow = X.n_rows, ncol = X.n_cols, s1 = 0, s2 = 0, s3=0, s4=0;
  int result;
  
  for(int j = 0; j< ncol; j++){
    for(int i = 0; i< nrow-1; i++){
      int diff = X(i,j)-X(i+1,j);
      if(diff == 0){
        s1 = s1 +1;
      }
      
    }
  }
  
  for(int i = 0; i< nrow; i++){
    for(int j = 0; j< ncol-1; j++){
      int diff = X(i,j)-X(i,j+1);
      if(diff == 0){
        s2 = s2 +1;
      }
      
    }
  }
  
  for(int i = 0; i<nrow; i++){
    int diff4 = X(i, ncol-1)-X(i,0);
    if(diff4 ==0){
      s3 = s3 +1;
    }
  }
  
  for(int i = 0; i<ncol; i++){
    int diff5 = X(nrow-1, i)-X(0,i);
    if(diff5 ==0){
      s4 = s4 +1;
    }
  }
  result = s1 + s2 +s3+s4;
  
  return(result);
}




// [[Rcpp::export]]
List pottsDMH(RawVector foo, int ncolor, double stat, double COV, double beta, int outer, int cycle,
              bool updateCOV, double sigma2, int adaptInterval, double adaptFactorExponent, int adapIter){
  
  double logprob, u, betaprev, betaprop, statprop;
  vec postSamples = zeros(outer), accprob = zeros(outer);
  double rhat, gamma1, gamma2, c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  double cholCOV = sqrt( sigma2 * COV );;
  
  
  for(int l = 0; l< outer; l++){
    
    if(updateCOV){
      if( (l+1 >= adaptInterval) && (l+1 - (adaptInterval * trunc((l+1) / adaptInterval)) == 0) ){
        rhat = sum( accprob.rows(l+1-adaptInterval, l-1) ) / (adaptInterval-1);
        gamma1 = 1 / pow(adapIter, c1);
        gamma2 = c0 * gamma1;
        sigma2 = exp( log(sigma2) + gamma2 * (rhat - ropt) );
        
        COV = COV + gamma1 * ( var( postSamples.rows(l+1-adaptInterval, l-1) ) - COV );
        cholCOV = sqrt( sigma2 * ( COV + 0.00001 ) );
        adapIter = adapIter + 1;
      }
    }
    betaprev = beta;
    betaprop = betaprev + cholCOV * randn();
    
    if(betaprop > 0){
      statprop = potts_stat(foo, ncolor, betaprop, cycle);
      
      logprob = -0.05 * betaprop * betaprop + 0.05 * betaprev * betaprev +
        (statprop - stat) * (betaprev - betaprop);
      
      u = log( randu() );
      if( u < logprob ){
        beta = betaprop;
        accprob[l] = 1;
      }
    }
    postSamples[l] =  beta;
    
    if ( (l+1) % 1000 == 0 ) {
      Rprintf("Generated %d samples...\n", l+1);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("postSamples") = postSamples,
                            Rcpp::Named("accprob") = accprob,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("COV") = COV);
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Generate summary statistics of auxiliary variables for given particles
mat pAuxSamp(RawVector foo, int ncolor, int cycle, vec Designmat, int m, int num){
  
  int thnrow = Designmat.n_elem;             // number of design points   
  mat H(thnrow,m);                       // cube summary statistics will be stored. (thnrow by m)
  omp_set_num_threads(num);
  
  
  //////    START OF BIGGEST CHAIN (M)     ////// 
  for(int M = 0; M < m; M++){               // m is number of importance sampling estimate 
    
    int i;
#pragma omp parallel shared(Designmat) private(i)
{	
#pragma omp for schedule(static)  
  for(i = 0; i < thnrow; i++){
    double thetaprop = Designmat(i);        
    double sumstat = potts_stat(foo, ncolor, thetaprop, cycle);
    H(i,M) = sumstat;	 	
  }
}
  }
  return(H);        	
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Generate summary statistics of auxiliary variables for given particles to construct IS estimate

vec pResponse(RawVector foo, int ncolor,  int cycle, double hatparameter, int m, int num){
  
  vec H(m);                         // matrix where summary statistics will be stored. (m by 2)
  omp_set_num_threads(num);
  
  //////    START OF BIGGEST CHAIN (M)     ////// 
  int M;
#pragma omp parallel shared(H) private(M)
{	
#pragma omp for schedule(static)  
  for(M = 0; M < m; M++){                            // m is number of importance sampling estimate 
    // summary statistics
    H(M) = potts_stat(foo, ncolor, hatparameter, cycle);
  }
}

return(H);        	
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// LikEm for Potts
vec pottsLikEm(int Niter, vec theta, double COV, double lhXZ, vec betahat, vec phihat, vec Designmat, vec y, double stat){
  int thnrow = Designmat.n_elem;                                        // number of design points
  double thetaprev;                                                   // befor propose in MCMC, previous parameters
  double lhXZp,logprob,u;                                               // used in MCMC step
  double negativeInf = -std::numeric_limits<float>::infinity();;	               
  double phi1hat = phihat[0], sigmasqhat = phihat[1]; 
  
  mat h1(thnrow,thnrow);   // used to construct Sigma in Gaussian process
  vec h1dcross(thnrow);     // used to construct cross Sigma in Gaussian process
  
  for(int i = 0; i < thnrow; i++){
    for(int j = 0; j <= i; j++){
      h1(i,j) = h1(j,i) = fabs(Designmat(i)-Designmat(j));
    }
  }
  mat Sigma = sigmasqhat*(1+sqrt(3)*h1/phi1hat)%exp(-sqrt(3)*h1/phi1hat);
  mat InvSigma = inv(Sigma);	    // Inverse of Sigma in Gaussian process
  mat Xth = ones(thnrow,1);
  Xth.insert_cols(1,Designmat);	// make design matrix for linear model 
  
  
  // Start of MCMC Chain 
  for(int k = 0; k< Niter-1; k++){
    
    double Znormal = randn(); // multivariate proposal by using Cholesky factorization
    thetaprev = theta(k);
    
    // proposed parameter and corresponding rho coefficients
    double thetaprop = thetaprev + Znormal*COV;
    
    // constranits on prior space
    if( thetaprop < 0 ){
      logprob = negativeInf;	
      
    }else{			
      for(int i = 0; i< thnrow; i++){  // Caculating cross covaraince matrix
        h1dcross[i] =  fabs(thetaprop-Designmat(i));	
      }
      mat Sigmacross = sigmasqhat*(1+sqrt(3)*h1dcross/phi1hat)%exp(-sqrt(3)*h1dcross/phi1hat);
      vec xpoint = ones(1);
      xpoint.insert_rows(1,thetaprop);
      lhXZp = (trans(xpoint)*betahat + trans(Sigmacross)* InvSigma*(y-Xth*betahat))[0]; //Gaussian kriging for intractable term
      
      logprob = -0.05 * thetaprop * thetaprop + 0.05 * thetaprev * thetaprev + lhXZp - lhXZ;              // log probability ratio to determine acceptance of MCMC 
    } 
    
    u = log( randu() );
    if( u< logprob ){
      theta(k+1) = thetaprop;
      lhXZ = lhXZp;		
    }else{
      theta(k+1) = thetaprev;
    }
    
  }
  
  return theta;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// LikEm for Potts
vec pottsNormEm(int Niter, vec theta, double COV, double logconst, vec betahat, vec phihat, vec Designmat, vec y, double stat){
  int thnrow = Designmat.n_elem;                                        // number of design points
  double thetaprev, thetaprop;                                                   // befor propose in MCMC, previous parameters
  double logconstp,logprob,u;                                               // used in MCMC step
  double negativeInf = -std::numeric_limits<float>::infinity();;	               
  double phi1hat = phihat[0], sigmasqhat = phihat[1]; 
  
  mat h1(thnrow,thnrow);   // used to construct Sigma in Gaussian process
  vec h1dcross(thnrow);     // used to construct cross Sigma in Gaussian process
  
  for(int i = 0; i < thnrow; i++){
    for(int j = 0; j <= i; j++){
      h1(i,j) = h1(j,i) = fabs(Designmat(i)-Designmat(j));
    }
  }
  mat Sigma = sigmasqhat*(1+sqrt(3)*h1/phi1hat)%exp(-sqrt(3)*h1/phi1hat);
  mat InvSigma = inv(Sigma);	    // Inverse of Sigma in Gaussian process
  mat Xth = ones(thnrow,1);
  Xth.insert_cols(1,Designmat);	// make design matrix for linear model 
  
  
  // Start of MCMC Chain 
  for(int k = 0; k< Niter-1; k++){
    
    double Znormal = randn(); // multivariate proposal by using Cholesky factorization
    thetaprev = theta(k);
    
    // proposed parameter and corresponding rho coefficients
    thetaprop = thetaprev + Znormal*COV;
    
    // constranits on prior space
    if( thetaprop < 0 ){
      logprob = negativeInf;	
      
    }else{			
      for(int i = 0; i< thnrow; i++){  // Caculating cross covaraince matrix
        h1dcross[i] =  fabs(thetaprop-Designmat(i));	
      }
      mat Sigmacross = sigmasqhat*(1+sqrt(3)*h1dcross/phi1hat)%exp(-sqrt(3)*h1dcross/phi1hat);
      logconstp = (betahat(0) + betahat(1)*thetaprop + trans(Sigmacross)* InvSigma*(y-Xth*betahat))[0]; //Gaussian kriging for intractable term
      
      logprob = -0.05 * thetaprop * thetaprop + 0.05 * thetaprev * thetaprev + stat * (thetaprop - thetaprev) + logconst - logconstp;    // log probability ratio to determine acceptance of MCMC 
    } 
    
    u = log( randu() );
    if( u< logprob ){
      theta(k+1) = thetaprop;
      logconst = logconstp;			
    }else{
      theta(k+1) = thetaprev;
    }
    
  }
  
  return theta;
}










//rbf kernel
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
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


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
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
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
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


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
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
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
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

// [[Rcpp::export]]
mat vec2mat( vec m, int t )
{
  int n = m.size();
  int ns = n/t;
  mat out(ns, t) ;
  
  for (int i = 0; i < ns; i++) 
  {
    out(i,0)=m[i];
    for (int j = 1; j < t; j++) 
    {
      out(i,j)=m[i+(j*ns)];
    }
  }
  return out ;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
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

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
mat repint( mat m, double t )
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
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
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
double allsum(mat x) {
  int ncol = x.n_cols;
  vec out(ncol);
  for (int i = 0; i < ncol; i++) {
    out[i] = sum(x.col(i));
  }
  
  double allsum = sum(out);
  return allsum;
}
arma::mat f_arma_mat_empty( arma::vec dim_mat ) {
  
  // Create an empty matrix with dimensions dim_mat
  return arma::mat(dim_mat[0], dim_mat[1]);
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



// [[Rcpp::export]]
vec AIKSi(mat th, mat score, vec weight, double c, double beta, int i, int k){
  int p = score.n_cols;
  double temp, temp2, bthi, bthip, k0, k0thi, k0thip, k0thithip;
  mat knot = zeros(k,p);
  
  for(int ip = i; ip < k; ip++){
    temp = pow(c, 2);
    for(int j = 0; j < p; j++){
      temp += pow(th(i,j) - th(ip,j), 2);
    }
    
    for(int j = 0; j < p; j++){
      bthi = score(i,j);
      bthip = score(ip,j);
      temp2 = th(i,j) - th(ip,j);
      
      k0 = pow(temp, beta);
      k0thi = 2*beta*pow(temp, beta-1)*temp2;
      k0thip = -k0thi;
      k0thithip = -2*beta*pow(temp,beta-2)*(2*pow(temp2,2)*(beta-1)+temp);
      
      knot(ip,j) = bthi*bthip*k0 + bthi*k0thip + bthip*k0thi + k0thithip;
    }
  }
  
  vec H = zeros(p);
  
  for(int j = 0; j < p; j++){
    double wsq = pow(weight[i],2)*knot(i,j);
    if(i < k-1){
      for(int m = (i+1); m < k; m++){
        wsq += weight[i]*knot(m,j)*weight[m]*2;
      }
    }
    H[j] = wsq;
  }
  
  return(H);
}



// [[Rcpp::export]]
// Compute w^2 for each row
mat AIKS(mat th, mat score, vec weight, double c, double beta, int k){
  int p = score.n_cols;
  mat H(k,p);
  
  for(int i = 0; i < k; i++){
    vec res = AIKSi(th, score, weight, c, beta, i, k);
    
    for(int j = 0; j < p; j++){
      H(i,j) = res[j];
    }
  }
  return(H);
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
List rbf2(mat theta, double h=-1){
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
  mat dykxy = dxkxy*(-1);
  mat dxykxy = dxkxy*theta.t()*(-2);
  return List::create(Named("kxy")=kxy, Named("dxkxy")=dxkxy, Named("dykxy")=dykxy, Named("dxykxy")=dxykxy);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
mat potts_is(RawVector foo, int ncolor, mat x0, mat X, int niter, int inner, int imporsize, double stepsize, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  int summary = Summary2(X);
  
  for(int iter=0; iter<niter; iter++){
    vec suf(theta.n_rows);
    for(int ii=0; ii<theta.n_rows; ii++){  
      for(int insize=0;  insize<imporsize; insize++){
        suf[ii]=suf[ii]+potts_stat(foo, ncolor, theta[ii],inner);
      }
    }    
    
    suf = suf/imporsize;
    
    mat lnpgrad = summary - suf;
    List ker=rbf(theta, -1);
    mat kxy = ker["kxy"];
    mat dxkxy = ker["dxkxy"];
    mat gradtheta = kxy*lnpgrad;
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
  }
  return(theta);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
mat potts_map(RawVector foo, int ncolor, mat x0, mat X, int niter, int inner, int imporsize, double stepsize, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  
  for(int iter=0; iter<niter; iter++){
    vec suf(theta.n_rows);
    for(int ii=0; ii<theta.n_rows; ii++){  
      for(int insize=0;  insize<imporsize; insize++){
        suf[ii]=suf[ii]+potts_stat(foo, ncolor, theta[ii],inner);
      }
    }    
    
    suf = suf/imporsize;
    
    mat lnpgrad = Summary2(X) - suf;
    List ker=rbf(theta, -1);
    mat kxy = ker["kxy"];
    mat dxkxy = ker["dxkxy"];
    mat gradtheta = lnpgrad;
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
  }
  return(theta);
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List potts_snis(RawVector foo, int ncolor, mat x0, mat X, double map, int niter, int inner, int imporsize, double stepsize, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repint(summary,Summary2(X));
  
  mat suf(imporsize, theta.n_cols);
  mat esss(niter, theta.n_rows);
  
  
  for(int insize=0;  insize<imporsize; insize++){
    suf.row(insize)=potts_stat(foo, ncolor, map, inner);
  }
  
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums(theta.n_rows,theta.n_cols);
  
  for(int iter=0; iter<niter; iter++){
    for(int tt=0; tt<theta.n_rows; tt++){
      mat selip(imporsize, theta.n_cols);
      mat selip2(imporsize, theta.n_cols);
      mat norwei(imporsize, theta.n_cols);
      for(int jj=0; jj<imporsize; jj++){
        selip.row(jj)=mmult((theta.row(tt)-map), suf.row(jj));
      }

      
      mat expselip=mexp(selip);
      for(int seliprow=0; seliprow<selip.n_rows; seliprow++){
        for(int selipcol=0; selipcol<selip.n_cols; selipcol++){
          if(expselip(seliprow, selipcol)==arma::datum::inf){
            expselip(seliprow,selipcol)=1e300;
          }
          if(expselip(seliprow, selipcol)==0){
            expselip(seliprow,selipcol)=-1e300;
          }
        }
      }
      vec selsum = colSum(expselip);
      for(int jj=0; jj<imporsize; jj++){
        selip2.row(jj)=mmult((expselip.row(jj)/selsum.t()), suf.row(jj));
        mat nor = (expselip.row(jj)/selsum.t());
        norwei.row(jj)= nor;
      }
      vec selsum2 = colSum(selip2);
      norwei = mpow(norwei,2); 
      vec norweisum = colSum(norwei);
      norweisums.row(tt)=norweisum.t();
      parsum.row(tt)=selsum2.t();
    }      
    
    ess = norweisums;
    mat midess=minv(ess);
    midess = midess.t();
    esss.row(iter)=midess.row(0);
    
    
    mat lnpgrad = summary - parsum.t() - 0.1*theta.t();
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
  }
  return List::create(Named("theta")=theta, Named("esss")=esss);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List potts_svi(RawVector foo, int ncolor, mat x0, mat X, double map, int niter, int inner, int imporsize, int imporsizes, double stepsize, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repint(summary,Summary2(X));
  
  mat suf(imporsize, theta.n_cols);
  mat esss(niter, theta.n_rows);
  
  
  for(int insize=0;  insize<imporsize; insize++){
    suf.row(insize)=potts_stat(foo, ncolor, map,inner);
  }
  
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums(theta.n_rows,theta.n_cols);
  
  for(int iter=0; iter<niter; iter++){
    for(int tt=0; tt<theta.n_rows; tt++){
      mat selip(imporsize, theta.n_cols);
      mat selip2(imporsize, theta.n_cols);
      mat norwei(imporsize, theta.n_cols);
      for(int jj=0; jj<imporsize; jj++){
        selip.row(jj)=mmult((theta.row(tt)-map), suf.row(jj));
      }
      mat expselip=mexp(selip);
      vec selsum = colSum(expselip);
      for(int jj=0; jj<imporsize; jj++){
        selip2.row(jj)=mmult((expselip.row(jj)/selsum.t()), suf.row(jj));
        mat nor = (expselip.row(jj)/selsum.t());
        norwei.row(jj)= nor;
      }
      vec selsum2 = colSum(selip2);
      norwei = mpow(norwei,2); 
      vec norweisum = colSum(norwei);
      norweisums.row(tt)=norweisum.t();
      parsum.row(tt)=selsum2.t();
    }      
    
    ess = norweisums;
    mat midess=minv(ess);
    midess = midess.t();
    esss.row(iter)=midess.row(0);
   
    for(int ss=0; ss<theta.n_rows; ss++){
      if(esss(iter,ss)<(imporsize/3)){ 
        parsum.row(ss)=0;
        for(int insizes=0; insizes<imporsizes; insizes++){
        parsum.row(ss)=parsum.row(ss)+potts_stat(foo, ncolor, theta(ss,0),inner);
      }
      }    
    }
   
    mat lnpgrad = summary - parsum.t() - 0.1*theta.t();
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
  }
  return List::create(Named("theta")=theta, Named("esss")=esss);
}





// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List potts_svi_new(RawVector foo, int ncolor, mat x0, mat X, double map, int niter, int inner, int imporsize, int imporsizes, double stepsize, double thres, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repint(summary,Summary2(X));
  
  
  mat suf(imporsize, theta.n_cols);
  mat esss(niter, theta.n_rows);
  mat thetas = int_to_mat(map);
  
  for(int insize=0;  insize<imporsize; insize++){
    suf.row(insize)=potts_stat(foo, ncolor, map,inner);
  }
  
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums(theta.n_rows,theta.n_cols);
  mat suftheta = suf;
  
  for(int iter=0; iter<niter; iter++){
    for(int tt=0; tt<theta.n_rows; tt++){
      if(crossdist2(theta.row(tt), int_to_mat(map))[0]> min(crossdist2(theta.row(tt), thetas))){
      
      int minwhich = vecminInd(crossdist2(theta.row(tt), thetas));
      mat selip(imporsizes, theta.n_cols);
      mat selip2(imporsizes, theta.n_cols);
      mat norwei(imporsizes, theta.n_cols);
      for(int jj=0; jj<imporsizes; jj++){
        selip.row(jj)=mmult((theta.row(tt)-thetas.row(minwhich)), suftheta.row(imporsize+minwhich*imporsizes+jj));
      }
      mat expselip=mexp(selip);
      vec selsum = colSum(expselip);
      for(int jj=0; jj<imporsizes; jj++){
        selip2.row(jj)=mmult((expselip.row(jj)/selsum.t()), suftheta.row(imporsize+minwhich*imporsizes+jj));
        mat nor = (expselip.row(jj)/selsum.t());
        norwei.row(jj)= nor;
      }
      vec selsum2 = colSum(selip2);
      norwei = mpow(norwei,2); 
      vec norweisum = colSum(norwei);
      norweisums.row(tt)=norweisum.t();
      parsum.row(tt)=selsum2.t();
      
      
      ess = norweisums;
      mat midess=minv(ess);
      midess = midess.t();
      esss.row(iter)=midess.row(0);
      
      
        if(esss(iter,tt)<(imporsizes/thres)) {
          thetas = rbind_mat(thetas, int_to_mat(theta(tt,0)));
          parsum.row(tt)=0;
          for(int insizes=0; insizes<imporsizes; insizes++){
          mat suff = int_to_mat(potts_stat(foo, ncolor, theta(tt,0),inner));
          parsum.row(tt)=parsum.row(tt)+suff;
          suf = rbind_mat(suf, suff);
         }
        }
      
      }
      
      else{ 
      mat selip(imporsize, theta.n_cols);
      mat selip2(imporsize, theta.n_cols);
      mat norwei(imporsize, theta.n_cols);
      for(int jj=0; jj<imporsize; jj++){
        selip.row(jj)=mmult((theta.row(tt)-map), suf.row(jj));
      }
      mat expselip=mexp(selip);
      vec selsum = colSum(expselip);
      for(int jj=0; jj<imporsize; jj++){
        selip2.row(jj)=mmult((expselip.row(jj)/selsum.t()), suf.row(jj));
        mat nor = (expselip.row(jj)/selsum.t());
        norwei.row(jj)= nor;
      }
      vec selsum2 = colSum(selip2);
      norwei = mpow(norwei,2); 
      vec norweisum = colSum(norwei);
      norweisums.row(tt)=norweisum.t();
      parsum.row(tt)=selsum2.t();
      
      ess = norweisums;
      mat midess=minv(ess);
      midess = midess.t();
      esss.row(iter)=midess.row(0);
      
      
      if(esss(iter,tt)<(imporsize/thres)) {
        thetas = rbind_mat(thetas, int_to_mat(theta(tt,0)));
        for(int insizes=0; insizes<imporsizes; insizes++){
          parsum.row(tt)=0;
          mat suff = int_to_mat(potts_stat(foo, ncolor, theta(tt,0),inner));
          parsum.row(tt)=parsum.row(tt)+suff;
          suf = rbind_mat(suf, suff);
        }
      }
      }
    }      
    
    
    mat lnpgrad = summary - parsum.t() - 0.1*theta.t();
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
  }
  return List::create(Named("theta")=theta, Named("esss")=esss, Named("thetas")=thetas);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List potts_svi_aiks(RawVector foo, int ncolor, mat x0, mat X, double map, int niter, int inner, int imporsize, int imporsizes, double stepsize, double thres, vec weight, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repint(summary,Summary2(X));
  vec aiksall = zeros(niter);
  mat alltheta(theta.n_rows*niter, theta.n_cols);
  int numthe = theta.n_rows;
  
  mat suf(imporsize, theta.n_cols);
  mat esss(niter, theta.n_rows);
  mat thetas = int_to_mat(map);
  
  for(int insize=0;  insize<imporsize; insize++){
    suf.row(insize)=potts_stat(foo, ncolor, map,inner);
  }
  
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums(theta.n_rows,theta.n_cols);
  mat suftheta = suf;
  
  for(int iter=0; iter<niter; iter++){
    for(int tt=0; tt<theta.n_rows; tt++){
      if(crossdist2(theta.row(tt), int_to_mat(map))[0]> min(crossdist2(theta.row(tt), thetas))){
        
        int minwhich = vecminInd(crossdist2(theta.row(tt), thetas));
        mat selip(imporsizes, theta.n_cols);
        mat selip2(imporsizes, theta.n_cols);
        mat norwei(imporsizes, theta.n_cols);
        for(int jj=0; jj<imporsizes; jj++){
          selip.row(jj)=mmult((theta.row(tt)-thetas.row(minwhich)), suftheta.row(imporsize+minwhich*imporsizes+jj));
        }
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsizes; jj++){
          selip2.row(jj)=mmult((expselip.row(jj)/selsum.t()), suftheta.row(imporsize+minwhich*imporsizes+jj));
          mat nor = (expselip.row(jj)/selsum.t());
          norwei.row(jj)= nor;
        }
        vec selsum2 = colSum(selip2);
        norwei = mpow(norwei,2); 
        vec norweisum = colSum(norwei);
        norweisums.row(tt)=norweisum.t();
        parsum.row(tt)=selsum2.t();
        
        
        ess = norweisums;
        mat midess=minv(ess);
        midess = midess.t();
        esss.row(iter)=midess.row(0);
        
        
        if(esss(iter,tt)<(imporsizes/thres)) {
          thetas = rbind_mat(thetas, int_to_mat(theta(tt,0)));
          parsum.row(tt)=0;
          for(int insizes=0; insizes<imporsizes; insizes++){
            mat suff = int_to_mat(potts_stat(foo, ncolor, theta(tt,0),inner));
            parsum.row(tt)=parsum.row(tt)+suff;
            suf = rbind_mat(suf, suff);
          }
        }
        
      }
      
      else{ 
        mat selip(imporsize, theta.n_cols);
        mat selip2(imporsize, theta.n_cols);
        mat norwei(imporsize, theta.n_cols);
        for(int jj=0; jj<imporsize; jj++){
          selip.row(jj)=mmult((theta.row(tt)-map), suf.row(jj));
        }
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsize; jj++){
          selip2.row(jj)=mmult((expselip.row(jj)/selsum.t()), suf.row(jj));
          mat nor = (expselip.row(jj)/selsum.t());
          norwei.row(jj)= nor;
        }
        vec selsum2 = colSum(selip2);
        norwei = mpow(norwei,2); 
        vec norweisum = colSum(norwei);
        norweisums.row(tt)=norweisum.t();
        parsum.row(tt)=selsum2.t();
        
        ess = norweisums;
        mat midess=minv(ess);
        midess = midess.t();
        esss.row(iter)=midess.row(0);
        
        
        if(esss(iter,tt)<(imporsize/thres)) {
          thetas = rbind_mat(thetas, int_to_mat(theta(tt,0)));
          for(int insizes=0; insizes<imporsizes; insizes++){
            parsum.row(tt)=0;
            mat suff = int_to_mat(potts_stat(foo, ncolor, theta(tt,0),inner));
            parsum.row(tt)=parsum.row(tt)+suff;
            suf = rbind_mat(suf, suff);
          }
        }
      }
    }      
    
    
    mat lnpgrad = summary - parsum.t() - 0.1*theta.t();
    List ker=rbf(theta, -1);
    mat kxy = ker["kxy"];
    mat dxkxy = ker["dxkxy"];
    mat aiks = AIKS(theta, lnpgrad.t(), weight, 1, -1/2, numthe);
    double aikss = allsum(aiks);
    aiksall[iter]=aikss;
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
    for(int jjj=0; jjj<numthe; jjj++){
      alltheta.row(jjj+(iter*numthe)) = theta.row(jjj);
    }
  }
  return List::create(Named("theta")=theta, Named("alltheta")=alltheta, Named("aiksall")=aiksall);
}







// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List potts_svi_aiks_max(RawVector foo, int ncolor, mat x0, mat X, double map, int niter, int inner, int imporsize, int imporsizes, double stepsize, double thres, vec weight, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repint(summary,Summary2(X));
  vec aiksall = zeros(niter);
  mat alltheta(theta.n_rows*niter, theta.n_cols);
  int numthe = theta.n_rows;
  
  mat suf(imporsize, theta.n_cols);
  mat esss(niter, theta.n_rows);
  mat thetas = int_to_mat(map);
  
  for(int insize=0;  insize<imporsize; insize++){
    suf.row(insize)=potts_stat(foo, ncolor, map,inner);
  }
  
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums(theta.n_rows,theta.n_cols);
  mat suftheta = suf;
  
  for(int iter=0; iter<niter; iter++){
    for(int tt=0; tt<theta.n_rows; tt++){
      if(crossdist2(theta.row(tt), int_to_mat(map))[0]> min(crossdist2(theta.row(tt), thetas))){
        
        int minwhich = vecminInd(crossdist2(theta.row(tt), thetas));
        mat selip(imporsizes, theta.n_cols);
        mat selip2(imporsizes, theta.n_cols);
        mat norwei(imporsizes, theta.n_cols);
        for(int jj=0; jj<imporsizes; jj++){
          selip.row(jj)=mmult((theta.row(tt)-thetas.row(minwhich)), suftheta.row(imporsize+minwhich*imporsizes+jj));
        }
        double mx = selip.max();
        selip= selip-mx;
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsizes; jj++){
          selip2.row(jj)=mmult((expselip.row(jj)/selsum.t()), suftheta.row(imporsize+minwhich*imporsizes+jj));
          mat nor = (expselip.row(jj)/selsum.t());
          norwei.row(jj)= nor;
        }
        vec selsum2 = colSum(selip2);
        norwei = mpow(norwei,2); 
        vec norweisum = colSum(norwei);
        norweisums.row(tt)=norweisum.t();
        parsum.row(tt)=selsum2.t();
        
        
        ess = norweisums;
        mat midess=minv(ess);
        midess = midess.t();
        esss.row(iter)=midess.row(0);
        
        
        if(esss(iter,tt)<(imporsizes/thres)) {
          thetas = rbind_mat(thetas, int_to_mat(theta(tt,0)));
          parsum.row(tt)=0;
          for(int insizes=0; insizes<imporsizes; insizes++){
            mat suff = int_to_mat(potts_stat(foo, ncolor, theta(tt,0),inner));
            parsum.row(tt)=parsum.row(tt)+suff;
            suf = rbind_mat(suf, suff);
          }
        }
        
      }
      
      else{ 
        mat selip(imporsize, theta.n_cols);
        mat selip2(imporsize, theta.n_cols);
        mat norwei(imporsize, theta.n_cols);
        for(int jj=0; jj<imporsize; jj++){
          selip.row(jj)=mmult((theta.row(tt)-map), suf.row(jj));
        }
        double mx = selip.max();
        selip= selip-mx;
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsize; jj++){
          selip2.row(jj)=mmult((expselip.row(jj)/selsum.t()), suf.row(jj));
          mat nor = (expselip.row(jj)/selsum.t());
          norwei.row(jj)= nor;
        }
        vec selsum2 = colSum(selip2);
        norwei = mpow(norwei,2); 
        vec norweisum = colSum(norwei);
        norweisums.row(tt)=norweisum.t();
        parsum.row(tt)=selsum2.t();
        
        ess = norweisums;
        mat midess=minv(ess);
        midess = midess.t();
        esss.row(iter)=midess.row(0);
        
        
        if(esss(iter,tt)<(imporsize/thres)) {
          thetas = rbind_mat(thetas, int_to_mat(theta(tt,0)));
          for(int insizes=0; insizes<imporsizes; insizes++){
            parsum.row(tt)=0;
            mat suff = int_to_mat(potts_stat(foo, ncolor, theta(tt,0),inner));
            parsum.row(tt)=parsum.row(tt)+suff;
            suf = rbind_mat(suf, suff);
          }
        }
      }
    }      
    
    
    mat lnpgrad = summary - parsum.t();
    List ker=rbf(theta, -1);
    mat kxy = ker["kxy"];
    mat dxkxy = ker["dxkxy"];
    mat aiks = AIKS(theta, lnpgrad.t(), weight, 1, -1/2, numthe);
    double aikss = allsum(aiks);
    aiksall[iter]=aikss;
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
    for(int jjj=0; jjj<numthe; jjj++){
      alltheta.row(jjj+(iter*numthe)) = theta.row(jjj);
    }
  }
  return List::create(Named("theta")=theta, Named("alltheta")=alltheta, Named("aiksall")=aiksall);
}





// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
mat potts_ispar(RawVector foo, int ncolor, vec x0, mat X, int niter, int inner, int imporsize, double stepsize, int core, double bandwidth=-1, double alpha=0.9){
  vec theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  int summary = Summary(X);
  
  for(int iter=0; iter<niter; iter++){
    mat thetamat = vec2mat(theta,core);
    mat suf(thetamat.n_rows, core);
    
    int jj;
#pragma omp parallel shared(thetamat) private(jj)
{
#pragma omp for schedule(static)
  for(jj=0; jj<core; jj++){
    for(int ii=0; ii<thetamat.n_rows; ii++){  
      for(int insize=0;  insize<imporsize; insize++){
        suf(ii,jj)=suf(ii,jj)+potts_stat(foo, ncolor ,thetamat(ii,jj),inner);
      }
    }    
  }
}
suf = suf.as_col();
suf = suf/imporsize;
theta = thetamat.as_col();

mat lnpgrad = summary - suf ;
List ker=rbf(theta, -1);
mat kxy = ker["kxy"];
mat dxkxy = ker["dxkxy"];
mat gradtheta = kxy*lnpgrad;
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
  }
  return(theta);
}







// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List potts_svi_par(RawVector foo, int ncolor, mat x0, mat X, double map, int niter, int inner, int imporsize, int imporsizes, double stepsize, int core, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repint(summary,Summary2(X));
  int ttt = theta.n_rows/core;
  
  mat suf(imporsize, theta.n_cols);
  mat esss(niter, theta.n_rows);
  
  
  for(int insize=0;  insize<imporsize; insize++){
    suf.row(insize)=potts_stat(foo, ncolor, map,inner);
  }
  
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums(theta.n_rows,theta.n_cols);
  
  for(int iter=0; iter<niter; iter++){
    for(int tt=0; tt<theta.n_rows; tt++){
      mat selip(imporsize, theta.n_cols);
      mat selip2(imporsize, theta.n_cols);
      mat norwei(imporsize, theta.n_cols);
      for(int jj=0; jj<imporsize; jj++){
        selip.row(jj)=mmult((theta.row(tt)-map), suf.row(jj));
      }
      mat expselip=mexp(selip);
      for(int seliprow=0; seliprow<selip.n_rows; seliprow++){
        for(int selipcol=0; selipcol<selip.n_cols; selipcol++){
          if(expselip(seliprow, selipcol)==arma::datum::inf){
            expselip(seliprow,selipcol)=1e300;
          }
          if(expselip(seliprow, selipcol)==0){
            expselip(seliprow,selipcol)=-1e300;
          }
        }
      }
      vec selsum = colSum(expselip);
      for(int jj=0; jj<imporsize; jj++){
        selip2.row(jj)=mmult((expselip.row(jj)/selsum.t()), suf.row(jj));
        mat nor = (expselip.row(jj)/selsum.t());
        norwei.row(jj)= nor;
      }
      vec selsum2 = colSum(selip2);
      norwei = mpow(norwei,2); 
      vec norweisum = colSum(norwei);
      norweisums.row(tt)=norweisum.t();
      parsum.row(tt)=selsum2.t();
    }      
    
    ess = norweisums;
    mat midess=minv(ess);
    midess = midess.t();
    esss.row(iter)=midess.row(0);
    
    int jjj;
#pragma omp parallel shared(theta) private(jjj)
{
#pragma omp for schedule(static)
  for(jjj=0; jjj<core; jjj++){
    for(int ii=0; ii<ttt; ii++){
      if(esss(iter,(ttt*jjj)+ii)<(imporsize/3)) for(int insizes=0; insizes<imporsizes; insizes++){
        parsum.row((ttt*jjj)+ii)=0;
        parsum.row((ttt*jjj)+ii)=parsum.row((ttt*jjj)+ii)+potts_stat(foo, ncolor, theta((ttt*jjj)+ii,0),inner);
      }
    }
  }
}
    
    
    mat lnpgrad = summary - parsum.t();
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
  }
  return List::create(Named("theta")=theta, Named("esss")=esss);
}




