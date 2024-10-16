/////////////////////////MCSVGD_ERGM (24.05.29 ver)/////////////////////////////////

#include <RcppArmadillo.h>
#include <limits>
#include <omp.h>
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
//#define min(x,y) (((x) < (y)) ? (x) : (y))

using namespace std;
using namespace Rcpp;
using namespace arma;



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Choose function
double Choose(int x, int y){
  double result1 = 1, result2 = 1, result;
  int iter = y;
  
  if( x< y ){ result = 0;	}else{	
    for(int i = 0; i<iter; i++){
      result1 = result1*x;
      result2 = result2*y;
      y = y-1;
      x = x-1;   	
    }	
    result = result1/result2;
  }
  return(result);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Count degree of each node 
vec countDegree(vec rowsum){
  int nrow = rowsum.n_elem, ind;
  vec degree = zeros( nrow );
  for(int i = 0; i < nrow; i++ ){
    ind = rowsum(i); 
    if(ind>0){ degree(ind-1) = degree(ind-1) + 1; }
  }
  return(degree);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Count edgewise shared partners 
vec countShared(vec rowsum, mat X){
  
  int numedges = sum(rowsum)/2, nrow = rowsum.n_elem ,ind;
  vec hist = zeros(numedges);
  int num = 0;
  for(int k = 0; k< nrow; k++){
    for(int j = 0; j< k+1; j++){
      if( X(k,j) == 1){  for(int i = 0; i<nrow; i++){ hist(num) =  hist(num) + X(i,k)*X(k,j)*X(j,i); }
      num = num + 1; }
    }
  }
  vec shared = zeros(nrow);
  for(int i = 0; i < numedges; i++ ){
    ind = hist(i);	
    if(ind>0){ shared(ind-1) = shared(ind-1) + 1; }
  }
  return(shared);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Summary statistics edges gwd gwesp
vec Summary(mat X, vec grade, vec sex){
  int nrow = X.n_rows;
  double decay = 0.25;
  vec rowsum = sum(X,1), result = zeros(10);
  
  // count degree of each node
  vec degree = countDegree(rowsum);
  
  // count edgewise shared partners 
  vec shared = countShared(rowsum,X);
  
  // Calculate summary statistics
  for(int i = 0; i< nrow; i++){
    result(0) = result(0) + Choose(rowsum(i),1);
    result(8) = result(8) + (   1- pow( (1-exp(-decay)),i+1)  )*degree(i);
    result(9) = result(9) + (   1- pow( (1-exp(-decay)),i+1)  )*shared(i);
    
    for(int j = 0; j< i; j++){
      result(1) = result(1) + X(i,j)*( (grade[i]==7)*(grade[j]==7)   );
      result(2) = result(2) + X(i,j)*( (grade[i]==8)*(grade[j]==8)   );
      result(3) = result(3) + X(i,j)*( (grade[i]==9)*(grade[j]==9)   );	
      result(4) = result(4) + X(i,j)*( (grade[i]==10)*(grade[j]==10)   );
      result(5) = result(5) + X(i,j)*( (grade[i]==11)*(grade[j]==11)   );		
      result(6) = result(6) + X(i,j)*( (grade[i]==12)*(grade[j]==12)   );		   
      result(7) = result(7) + X(i,j)*( sex[i]==sex[j] );		   		   	    		   		   	    	
    }
  }
  result(0)=result(0)/2;
  result(8)=exp(decay)*result(8);
  result(9)=exp(decay)*result(9);   
  
  return result;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// random scan gibbs update
vec Gibbs(mat X, vec grade, vec sex, vec coef, int cycle){
  int nrow = X.n_rows,indi,indj,prei,prej,res;
  double decay = 0.25;
  vec star = sum(X,1), changestat = zeros(10), Sumstat = Summary(X,grade,sex) ;
  rowvec ivec = trans( zeros(nrow) ), jvec = trans( zeros(nrow) );
  vec geoweight = zeros(nrow);
  for(int i =0; i<nrow; i++){ geoweight(i) = 1- pow( (1-exp(-decay)),i+1); }
  
  
  for(int l = 0; l< cycle; l++){
    for(int i = 1; i< nrow; i++){
      for(int j = 0; j< i; j++){
        ivec(i) = 1; jvec(j) = 1;
        res = 0; 
        
        // When Xij is 0
        if(X(i,j)==0){  
          indi = star(i) + 1, indj = star(j) + 1; 
          // change statistics of edge
          changestat(0) = ( Choose(indi,1) + Choose(indj,1) - Choose(star(i),1) - Choose(star(j),1) )/2;
          
          // change statistics of nodematch
          changestat(1) = ( (grade[i]==7)*(grade[j]==7)   );
          changestat(2) = ( (grade[i]==8)*(grade[j]==8)   );
          changestat(3) = ( (grade[i]==9)*(grade[j]==9)   );	       
          changestat(4) = ( (grade[i]==10)*(grade[j]==10)   );
          changestat(5) = ( (grade[i]==11)*(grade[j]==11)   );            
          changestat(6) = ( (grade[i]==12)*(grade[j]==12)   );       		               
          changestat(7) = (  sex[i]==sex[j]  );  		               
          
          // change statistics of gwd
          changestat(8) = exp(decay)*( geoweight(indi-1) + geoweight(indj-1) );
          if(indi-1>0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indi-1-1); }
          if(indj-1>0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indj-1-1); }
          
          // change statistics of gwes
          changestat(9) = 0;     
          for(int k = 0; k< nrow; k++){
            if(   ( X(i,k) == 1 ) & ( X(j,k) == 1  )   ){
              prei = sum( X.row(i)%X.row(k) );
              prej = sum( X.row(j)%X.row(k) );
              
              changestat(9) = changestat(9) + exp(decay)*geoweight(prei+1-1) ; 
              if(prei >0){ changestat(9) = changestat(9) - exp(decay)*geoweight(prei-1);} 
              changestat(9) = changestat(9) + exp(decay)*geoweight(prej+1-1) ; 
              if(prej >0){ changestat(9) = changestat(9) - exp(decay)*geoweight(prej-1);} 
              res = res + 1; // X(k,i) multiply X(i,j) +1 multiply X(j,k)
            }
          }
          if(res > 0){ changestat(9) = changestat(9) + exp(decay)*geoweight(res-1); }
          
          // accept reject step  	       
          // probability of Xij = 1 for given others are fixed 
          vec r =  exp(trans(coef)*changestat) ;         
          double p = r(0)/(1+r(0));   		   
          if( randu() < p  ){
            X(i,j) = X(j,i) = 1; 
            star(i) = indi; star(j) = indj;		 
            Sumstat = Sumstat + changestat;
          }
          
          
          // When Xij is 1  
        }else{
          indi = star(i) - 1, indj = star(j) - 1;	
          // change statistics of edge
          changestat(0) = ( Choose(star(i),1) + Choose(star(j),1) - Choose(indi,1) - Choose(indj,1) )/2;
          
          // change statistics of nodematch
          changestat(1) = ( (grade[i]==7)*(grade[j]==7)   );
          changestat(2) = ( (grade[i]==8)*(grade[j]==8)   );
          changestat(3) = ( (grade[i]==9)*(grade[j]==9)   );	       
          changestat(4) = ( (grade[i]==10)*(grade[j]==10)   );
          changestat(5) = ( (grade[i]==11)*(grade[j]==11)   );           
          changestat(6) = ( (grade[i]==12)*(grade[j]==12)   );     		              
          changestat(7) = (  sex[i]==sex[j]  );     		              
          
          // change statistics of gwd
          changestat(8) = exp(decay)*( geoweight(indi-1+1) + geoweight(indj-1+1)  );
          if(indi-1>=0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indi-1); }       
          if(indj-1>=0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indj-1); }        		   
          
          // change statistics of gwesp 
          changestat(9) = 0;
          for(int k = 0; k< nrow; k++){
            if(   ( X(i,k) == 1 ) & ( X(j,k) == 1  )   ){
              prei = sum( X.row(i)%X.row(k) );
              prej = sum( X.row(j)%X.row(k) );
              if(prei-1 >0){changestat(9) = changestat(9) - exp(decay)*geoweight(prei-1-1); } 
              changestat(9) = changestat(9) + exp(decay)*geoweight(prei-1); 
              if(prej-1 >0){changestat(9) = changestat(9) - exp(decay)*geoweight(prej-1-1); } 
              changestat(9) = changestat(9) + exp(decay)*geoweight(prej-1); 
              res = res + 1; // X(k,i) multiply X(i,j)  multiply X(j,k)
            }
          }
          if(res > 0){ changestat(9) = changestat(9) + exp(decay)*geoweight(res-1); }
          
          // accept reject step  	       
          // probability of Xij = 1 for given others are fixed 
          vec r =  exp(trans(coef)*changestat) ;         
          double p = r(0)/(1+r(0));   		   
          if( randu() > p  ){
            X(i,j) = X(j,i) = 0; 
            star(i) = indi; star(j) = indj;  
            Sumstat = Sumstat - changestat;
          }             
        }
        // End of a single update  
        ivec(i) = 0; jvec(j) = 0;
      }
    }   
  }
  
  return(Sumstat); 
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// random scan gibbs update
mat Gibbs2(mat X, vec grade, vec sex, vec coef, int cycle){
  mat allinner(10, cycle);
  int nrow = X.n_rows,indi,indj,prei,prej,res;
  double decay = 0.25;
  vec star = sum(X,1), changestat = zeros(10), Sumstat = Summary(X,grade,sex) ;
  rowvec ivec = trans( zeros(nrow) ), jvec = trans( zeros(nrow) );
  vec geoweight = zeros(nrow);
  for(int i =0; i<nrow; i++){ geoweight(i) = 1- pow( (1-exp(-decay)),i+1); }
  
  for(int l = 0; l< cycle; l++){
    for(int i = 1; i< nrow; i++){
      for(int j = 0; j< i; j++){
        ivec(i) = 1; jvec(j) = 1;
        res = 0; 
        
        // When Xij is 0
        if(X(i,j)==0){  
          indi = star(i) + 1, indj = star(j) + 1; 
          // change statistics of edge
          changestat(0) = ( Choose(indi,1) + Choose(indj,1) - Choose(star(i),1) - Choose(star(j),1) )/2;
          
          // change statistics of nodematch
          changestat(1) = ( (grade[i]==7)*(grade[j]==7)   );
          changestat(2) = ( (grade[i]==8)*(grade[j]==8)   );
          changestat(3) = ( (grade[i]==9)*(grade[j]==9)   );	       
          changestat(4) = ( (grade[i]==10)*(grade[j]==10)   );
          changestat(5) = ( (grade[i]==11)*(grade[j]==11)   );            
          changestat(6) = ( (grade[i]==12)*(grade[j]==12)   );       		               
          changestat(7) = (  sex[i]==sex[j]  );  		               
          
          // change statistics of gwd
          changestat(8) = exp(decay)*( geoweight(indi-1) + geoweight(indj-1) );
          if(indi-1>0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indi-1-1); }
          if(indj-1>0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indj-1-1); }
          
          // change statistics of gwes
          changestat(9) = 0;     
          for(int k = 0; k< nrow; k++){
            if(   ( X(i,k) == 1 ) & ( X(j,k) == 1  )   ){
              prei = sum( X.row(i)%X.row(k) );
              prej = sum( X.row(j)%X.row(k) );
              
              changestat(9) = changestat(9) + exp(decay)*geoweight(prei+1-1) ; 
              if(prei >0){ changestat(9) = changestat(9) - exp(decay)*geoweight(prei-1);} 
              changestat(9) = changestat(9) + exp(decay)*geoweight(prej+1-1) ; 
              if(prej >0){ changestat(9) = changestat(9) - exp(decay)*geoweight(prej-1);} 
              res = res + 1; // X(k,i) multiply X(i,j) +1 multiply X(j,k)
            }
          }
          if(res > 0){ changestat(9) = changestat(9) + exp(decay)*geoweight(res-1); }
          
          // accept reject step  	       
          // probability of Xij = 1 for given others are fixed 
          vec r =  exp(trans(coef)*changestat) ;         
          double p = r(0)/(1+r(0));   		   
          if( randu() < p  ){
            X(i,j) = X(j,i) = 1; 
            star(i) = indi; star(j) = indj;		 
            Sumstat = Sumstat + changestat;
          }
          
          
          // When Xij is 1  
        }else{
          indi = star(i) - 1, indj = star(j) - 1;	
          // change statistics of edge
          changestat(0) = ( Choose(star(i),1) + Choose(star(j),1) - Choose(indi,1) - Choose(indj,1) )/2;
          
          // change statistics of nodematch
          changestat(1) = ( (grade[i]==7)*(grade[j]==7)   );
          changestat(2) = ( (grade[i]==8)*(grade[j]==8)   );
          changestat(3) = ( (grade[i]==9)*(grade[j]==9)   );	       
          changestat(4) = ( (grade[i]==10)*(grade[j]==10)   );
          changestat(5) = ( (grade[i]==11)*(grade[j]==11)   );           
          changestat(6) = ( (grade[i]==12)*(grade[j]==12)   );     		              
          changestat(7) = (  sex[i]==sex[j]  );     		              
          
          // change statistics of gwd
          changestat(8) = exp(decay)*( geoweight(indi-1+1) + geoweight(indj-1+1)  );
          if(indi-1>=0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indi-1); }       
          if(indj-1>=0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indj-1); }        		   
          
          // change statistics of gwesp 
          changestat(9) = 0;
          for(int k = 0; k< nrow; k++){
            if(   ( X(i,k) == 1 ) & ( X(j,k) == 1  )   ){
              prei = sum( X.row(i)%X.row(k) );
              prej = sum( X.row(j)%X.row(k) );
              if(prei-1 >0){changestat(9) = changestat(9) - exp(decay)*geoweight(prei-1-1); } 
              changestat(9) = changestat(9) + exp(decay)*geoweight(prei-1); 
              if(prej-1 >0){changestat(9) = changestat(9) - exp(decay)*geoweight(prej-1-1); } 
              changestat(9) = changestat(9) + exp(decay)*geoweight(prej-1); 
              res = res + 1; // X(k,i) multiply X(i,j)  multiply X(j,k)
            }
          }
          if(res > 0){ changestat(9) = changestat(9) + exp(decay)*geoweight(res-1); }
          
          // accept reject step  	       
          // probability of Xij = 1 for given others are fixed 
          vec r =  exp(trans(coef)*changestat) ;         
          double p = r(0)/(1+r(0));   		   
          if( randu() > p  ){
            X(i,j) = X(j,i) = 0; 
            star(i) = indi; star(j) = indj;  
            Sumstat = Sumstat - changestat;
          }
          
        }
        // End of a single update  
        ivec(i) = 0; jvec(j) = 0;
      }
    }
    allinner.col(l)=Sumstat;
  }
  
  return(allinner); 
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// random scan gibbs update
mat Gibbs3(mat X, vec grade, vec sex, vec coef, int cycle, int impor){
  mat allinner(10, impor);
  int nrow = X.n_rows,indi,indj,prei,prej,res;
  double decay = 0.25;
  vec star = sum(X,1), changestat = zeros(10), Sumstat = Summary(X,grade,sex) ;
  rowvec ivec = trans( zeros(nrow) ), jvec = trans( zeros(nrow) );
  vec geoweight = zeros(nrow);
  for(int i =0; i<nrow; i++){ geoweight(i) = 1- pow( (1-exp(-decay)),i+1); }
  int numm = 0;
  
  for(int l = 0; l< cycle; l++){
    for(int i = 1; i< nrow; i++){
      numm=numm+1;
      for(int j = 0; j< i; j++){
        ivec(i) = 1; jvec(j) = 1;
        res = 0; 
        
        // When Xij is 0
        if(X(i,j)==0){  
          indi = star(i) + 1, indj = star(j) + 1; 
          // change statistics of edge
          changestat(0) = ( Choose(indi,1) + Choose(indj,1) - Choose(star(i),1) - Choose(star(j),1) )/2;
          
          // change statistics of nodematch
          changestat(1) = ( (grade[i]==7)*(grade[j]==7)   );
          changestat(2) = ( (grade[i]==8)*(grade[j]==8)   );
          changestat(3) = ( (grade[i]==9)*(grade[j]==9)   );	       
          changestat(4) = ( (grade[i]==10)*(grade[j]==10)   );
          changestat(5) = ( (grade[i]==11)*(grade[j]==11)   );            
          changestat(6) = ( (grade[i]==12)*(grade[j]==12)   );       		               
          changestat(7) = (  sex[i]==sex[j]  );  		               
          
          // change statistics of gwd
          changestat(8) = exp(decay)*( geoweight(indi-1) + geoweight(indj-1) );
          if(indi-1>0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indi-1-1); }
          if(indj-1>0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indj-1-1); }
          
          // change statistics of gwes
          changestat(9) = 0;     
          for(int k = 0; k< nrow; k++){
            if(   ( X(i,k) == 1 ) & ( X(j,k) == 1  )   ){
              prei = sum( X.row(i)%X.row(k) );
              prej = sum( X.row(j)%X.row(k) );
              
              changestat(9) = changestat(9) + exp(decay)*geoweight(prei+1-1) ; 
              if(prei >0){ changestat(9) = changestat(9) - exp(decay)*geoweight(prei-1);} 
              changestat(9) = changestat(9) + exp(decay)*geoweight(prej+1-1) ; 
              if(prej >0){ changestat(9) = changestat(9) - exp(decay)*geoweight(prej-1);} 
              res = res + 1; // X(k,i) multiply X(i,j) +1 multiply X(j,k)
            }
          }
          if(res > 0){ changestat(9) = changestat(9) + exp(decay)*geoweight(res-1); }
          
          // accept reject step  	       
          // probability of Xij = 1 for given others are fixed 
          vec r =  exp(trans(coef)*changestat) ;         
          double p = r(0)/(1+r(0));   		   
          if( randu() < p  ){
            X(i,j) = X(j,i) = 1; 
            star(i) = indi; star(j) = indj;		 
            Sumstat = Sumstat + changestat;
          }
          
          
          // When Xij is 1  
        }else{
          indi = star(i) - 1, indj = star(j) - 1;	
          // change statistics of edge
          changestat(0) = ( Choose(star(i),1) + Choose(star(j),1) - Choose(indi,1) - Choose(indj,1) )/2;
          
          // change statistics of nodematch
          changestat(1) = ( (grade[i]==7)*(grade[j]==7)   );
          changestat(2) = ( (grade[i]==8)*(grade[j]==8)   );
          changestat(3) = ( (grade[i]==9)*(grade[j]==9)   );	       
          changestat(4) = ( (grade[i]==10)*(grade[j]==10)   );
          changestat(5) = ( (grade[i]==11)*(grade[j]==11)   );           
          changestat(6) = ( (grade[i]==12)*(grade[j]==12)   );     		              
          changestat(7) = (  sex[i]==sex[j]  );     		              
          
          // change statistics of gwd
          changestat(8) = exp(decay)*( geoweight(indi-1+1) + geoweight(indj-1+1)  );
          if(indi-1>=0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indi-1); }       
          if(indj-1>=0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indj-1); }        		   
          
          // change statistics of gwesp 
          changestat(9) = 0;
          for(int k = 0; k< nrow; k++){
            if(   ( X(i,k) == 1 ) & ( X(j,k) == 1  )   ){
              prei = sum( X.row(i)%X.row(k) );
              prej = sum( X.row(j)%X.row(k) );
              if(prei-1 >0){changestat(9) = changestat(9) - exp(decay)*geoweight(prei-1-1); } 
              changestat(9) = changestat(9) + exp(decay)*geoweight(prei-1); 
              if(prej-1 >0){changestat(9) = changestat(9) - exp(decay)*geoweight(prej-1-1); } 
              changestat(9) = changestat(9) + exp(decay)*geoweight(prej-1); 
              res = res + 1; // X(k,i) multiply X(i,j)  multiply X(j,k)
            }
          }
          if(res > 0){ changestat(9) = changestat(9) + exp(decay)*geoweight(res-1); }
          
          // accept reject step  	       
          // probability of Xij = 1 for given others are fixed 
          vec r =  exp(trans(coef)*changestat) ;         
          double p = r(0)/(1+r(0));   		   
          if( randu() > p  ){
            X(i,j) = X(j,i) = 0; 
            star(i) = indi; star(j) = indj;  
            Sumstat = Sumstat - changestat;
          }
          
        }
        // End of a single update  
        ivec(i) = 0; jvec(j) = 0;
      }
      int seat = numm%impor;
      allinner.col(seat)=Sumstat;
    }
  }
  
  return(allinner); 
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// random scan gibbs update
mat Gibbs4(mat X, vec grade, vec sex, vec coef, int cycle, int impor){
  mat allinner(10, impor);
  int nrow = X.n_rows,indi,indj,prei,prej,res;
  double decay = 0.25;
  vec star = sum(X,1), changestat = zeros(10), Sumstat = Summary(X,grade,sex) ;
  rowvec ivec = trans( zeros(nrow) ), jvec = trans( zeros(nrow) );
  vec geoweight = zeros(nrow);
  for(int i =0; i<nrow; i++){ geoweight(i) = 1- pow( (1-exp(-decay)),i+1); }
  int numm = 0;
  int numm2 = 0;
  
  for(int l = 0; l< cycle; l++){
    for(int i = 1; i< nrow; i++){
      for(int j = 0; j< i; j++){
        numm=numm+1;
        if((numm&1000)==0){numm2=numm2+1;
        }
        ivec(i) = 1; jvec(j) = 1;
        res = 0; 
        
        // When Xij is 0
        if(X(i,j)==0){  
          indi = star(i) + 1, indj = star(j) + 1; 
          // change statistics of edge
          changestat(0) = ( Choose(indi,1) + Choose(indj,1) - Choose(star(i),1) - Choose(star(j),1) )/2;
          
          // change statistics of nodematch
          changestat(1) = ( (grade[i]==7)*(grade[j]==7)   );
          changestat(2) = ( (grade[i]==8)*(grade[j]==8)   );
          changestat(3) = ( (grade[i]==9)*(grade[j]==9)   );	       
          changestat(4) = ( (grade[i]==10)*(grade[j]==10)   );
          changestat(5) = ( (grade[i]==11)*(grade[j]==11)   );            
          changestat(6) = ( (grade[i]==12)*(grade[j]==12)   );       		               
          changestat(7) = (  sex[i]==sex[j]  );  		               
          
          // change statistics of gwd
          changestat(8) = exp(decay)*( geoweight(indi-1) + geoweight(indj-1) );
          if(indi-1>0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indi-1-1); }
          if(indj-1>0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indj-1-1); }
          
          // change statistics of gwes
          changestat(9) = 0;     
          for(int k = 0; k< nrow; k++){
            if(   ( X(i,k) == 1 ) & ( X(j,k) == 1  )   ){
              prei = sum( X.row(i)%X.row(k) );
              prej = sum( X.row(j)%X.row(k) );
              
              changestat(9) = changestat(9) + exp(decay)*geoweight(prei+1-1) ; 
              if(prei >0){ changestat(9) = changestat(9) - exp(decay)*geoweight(prei-1);} 
              changestat(9) = changestat(9) + exp(decay)*geoweight(prej+1-1) ; 
              if(prej >0){ changestat(9) = changestat(9) - exp(decay)*geoweight(prej-1);} 
              res = res + 1; // X(k,i) multiply X(i,j) +1 multiply X(j,k)
            }
          }
          if(res > 0){ changestat(9) = changestat(9) + exp(decay)*geoweight(res-1); }
          
          // accept reject step  	       
          // probability of Xij = 1 for given others are fixed 
          vec r =  exp(trans(coef)*changestat) ;         
          double p = r(0)/(1+r(0));   		   
          if( randu() < p  ){
            X(i,j) = X(j,i) = 1; 
            star(i) = indi; star(j) = indj;		 
            Sumstat = Sumstat + changestat;
          }
          
          
          // When Xij is 1  
        }else{
          indi = star(i) - 1, indj = star(j) - 1;	
          // change statistics of edge
          changestat(0) = ( Choose(star(i),1) + Choose(star(j),1) - Choose(indi,1) - Choose(indj,1) )/2;
          
          // change statistics of nodematch
          changestat(1) = ( (grade[i]==7)*(grade[j]==7)   );
          changestat(2) = ( (grade[i]==8)*(grade[j]==8)   );
          changestat(3) = ( (grade[i]==9)*(grade[j]==9)   );	       
          changestat(4) = ( (grade[i]==10)*(grade[j]==10)   );
          changestat(5) = ( (grade[i]==11)*(grade[j]==11)   );           
          changestat(6) = ( (grade[i]==12)*(grade[j]==12)   );     		              
          changestat(7) = (  sex[i]==sex[j]  );     		              
          
          // change statistics of gwd
          changestat(8) = exp(decay)*( geoweight(indi-1+1) + geoweight(indj-1+1)  );
          if(indi-1>=0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indi-1); }       
          if(indj-1>=0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indj-1); }        		   
          
          // change statistics of gwesp 
          changestat(9) = 0;
          for(int k = 0; k< nrow; k++){
            if(   ( X(i,k) == 1 ) & ( X(j,k) == 1  )   ){
              prei = sum( X.row(i)%X.row(k) );
              prej = sum( X.row(j)%X.row(k) );
              if(prei-1 >0){changestat(9) = changestat(9) - exp(decay)*geoweight(prei-1-1); } 
              changestat(9) = changestat(9) + exp(decay)*geoweight(prei-1); 
              if(prej-1 >0){changestat(9) = changestat(9) - exp(decay)*geoweight(prej-1-1); } 
              changestat(9) = changestat(9) + exp(decay)*geoweight(prej-1); 
              res = res + 1; // X(k,i) multiply X(i,j)  multiply X(j,k)
            }
          }
          if(res > 0){ changestat(9) = changestat(9) + exp(decay)*geoweight(res-1); }
          
          // accept reject step  	       
          // probability of Xij = 1 for given others are fixed 
          vec r =  exp(trans(coef)*changestat) ;         
          double p = r(0)/(1+r(0));   		   
          if( randu() > p  ){
            X(i,j) = X(j,i) = 0; 
            star(i) = indi; star(j) = indj;  
            Sumstat = Sumstat - changestat;
          }
          
        }
        // End of a single update  
        ivec(i) = 0; jvec(j) = 0;
      }
      if((numm%1000)==0){
        int seat = numm2%impor;
        allinner.col(seat)=Sumstat;
      }
    }
  }
  
  return(allinner); 
}







// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// random scan gibbs update
mat Gibbs5(mat X, vec grade, vec sex, vec coef, int cycle, int impor, int cyc, int thin){
  mat allinner(10, impor);
  int nrow = X.n_rows,indi,indj,prei,prej,res;
  double decay = 0.25;
  vec star = sum(X,1), changestat = zeros(10), Sumstat = Summary(X,grade,sex) ;
  rowvec ivec = trans( zeros(nrow) ), jvec = trans( zeros(nrow) );
  vec geoweight = zeros(nrow);
  for(int i =0; i<nrow; i++){ geoweight(i) = 1- pow( (1-exp(-decay)),i+1); }
  int numm = 0;
  int numm2 = 0;
  
  for(int l = 0; l< (cycle+cyc); l++){
    for(int i = 1; i< nrow; i++){
      for(int j = 0; j< i; j++){
        ivec(i) = 1; jvec(j) = 1;
        res = 0; 
        
        // When Xij is 0
        if(X(i,j)==0){  
          indi = star(i) + 1, indj = star(j) + 1; 
          // change statistics of edge
          changestat(0) = ( Choose(indi,1) + Choose(indj,1) - Choose(star(i),1) - Choose(star(j),1) )/2;
          
          // change statistics of nodematch
          changestat(1) = ( (grade[i]==7)*(grade[j]==7)   );
          changestat(2) = ( (grade[i]==8)*(grade[j]==8)   );
          changestat(3) = ( (grade[i]==9)*(grade[j]==9)   );	       
          changestat(4) = ( (grade[i]==10)*(grade[j]==10)   );
          changestat(5) = ( (grade[i]==11)*(grade[j]==11)   );            
          changestat(6) = ( (grade[i]==12)*(grade[j]==12)   );       		               
          changestat(7) = (  sex[i]==sex[j]  );  		               
          
          // change statistics of gwd
          changestat(8) = exp(decay)*( geoweight(indi-1) + geoweight(indj-1) );
          if(indi-1>0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indi-1-1); }
          if(indj-1>0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indj-1-1); }
          
          // change statistics of gwes
          changestat(9) = 0;     
          for(int k = 0; k< nrow; k++){
            if(   ( X(i,k) == 1 ) & ( X(j,k) == 1  )   ){
              prei = sum( X.row(i)%X.row(k) );
              prej = sum( X.row(j)%X.row(k) );
              
              changestat(9) = changestat(9) + exp(decay)*geoweight(prei+1-1) ; 
              if(prei >0){ changestat(9) = changestat(9) - exp(decay)*geoweight(prei-1);} 
              changestat(9) = changestat(9) + exp(decay)*geoweight(prej+1-1) ; 
              if(prej >0){ changestat(9) = changestat(9) - exp(decay)*geoweight(prej-1);} 
              res = res + 1; // X(k,i) multiply X(i,j) +1 multiply X(j,k)
            }
          }
          if(res > 0){ changestat(9) = changestat(9) + exp(decay)*geoweight(res-1); }
          
          // accept reject step  	       
          // probability of Xij = 1 for given others are fixed 
          vec r =  exp(trans(coef)*changestat) ;         
          double p = r(0)/(1+r(0));   		   
          if( randu() < p  ){
            X(i,j) = X(j,i) = 1; 
            star(i) = indi; star(j) = indj;		 
            Sumstat = Sumstat + changestat;
          }
          
          
          // When Xij is 1  
        }else{
          indi = star(i) - 1, indj = star(j) - 1;	
          // change statistics of edge
          changestat(0) = ( Choose(star(i),1) + Choose(star(j),1) - Choose(indi,1) - Choose(indj,1) )/2;
          
          // change statistics of nodematch
          changestat(1) = ( (grade[i]==7)*(grade[j]==7)   );
          changestat(2) = ( (grade[i]==8)*(grade[j]==8)   );
          changestat(3) = ( (grade[i]==9)*(grade[j]==9)   );	       
          changestat(4) = ( (grade[i]==10)*(grade[j]==10)   );
          changestat(5) = ( (grade[i]==11)*(grade[j]==11)   );           
          changestat(6) = ( (grade[i]==12)*(grade[j]==12)   );     		              
          changestat(7) = (  sex[i]==sex[j]  );     		              
          
          // change statistics of gwd
          changestat(8) = exp(decay)*( geoweight(indi-1+1) + geoweight(indj-1+1)  );
          if(indi-1>=0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indi-1); }       
          if(indj-1>=0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indj-1); }        		   
          
          // change statistics of gwesp 
          changestat(9) = 0;
          for(int k = 0; k< nrow; k++){
            if(   ( X(i,k) == 1 ) & ( X(j,k) == 1  )   ){
              prei = sum( X.row(i)%X.row(k) );
              prej = sum( X.row(j)%X.row(k) );
              if(prei-1 >0){changestat(9) = changestat(9) - exp(decay)*geoweight(prei-1-1); } 
              changestat(9) = changestat(9) + exp(decay)*geoweight(prei-1); 
              if(prej-1 >0){changestat(9) = changestat(9) - exp(decay)*geoweight(prej-1-1); } 
              changestat(9) = changestat(9) + exp(decay)*geoweight(prej-1); 
              res = res + 1; // X(k,i) multiply X(i,j)  multiply X(j,k)
            }
          }
          if(res > 0){ changestat(9) = changestat(9) + exp(decay)*geoweight(res-1); }
          
          // accept reject step  	       
          // probability of Xij = 1 for given others are fixed 
          vec r =  exp(trans(coef)*changestat) ;         
          double p = r(0)/(1+r(0));   		   
          if( randu() > p  ){
            X(i,j) = X(j,i) = 0; 
            star(i) = indi; star(j) = indj;  
            Sumstat = Sumstat - changestat;
          }
          
        }
        // End of a single update  
        ivec(i) = 0; jvec(j) = 0;
        if((l+1)>cycle){
          numm=numm+1;
          if((numm%thin)==0){
            numm2=numm2+1;
            int seat = (numm2%impor);
            allinner.col(seat)=Sumstat;
          }
        }
      }
    }
  }
  
  return(allinner); 
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// random scan gibbs update
vec Gibbsall(mat X, vec grade, vec sex, vec coef, int cycle){
  vec allinner = zeros(10);
  int nrow = X.n_rows,indi,indj,prei,prej,res;
  double decay = 0.25;
  vec star = sum(X,1), changestat = zeros(10), Sumstat = Summary(X,grade,sex) ;
  rowvec ivec = trans( zeros(nrow) ), jvec = trans( zeros(nrow) );
  vec geoweight = zeros(nrow);
  for(int i =0; i<nrow; i++){ geoweight(i) = 1- pow( (1-exp(-decay)),i+1); }
  
  
  for(int l = 0; l< cycle; l++){
    for(int i = 1; i< nrow; i++){
      for(int j = 0; j< i; j++){
        ivec(i) = 1; jvec(j) = 1;
        res = 0; 
        
        // When Xij is 0
        if(X(i,j)==0){  
          indi = star(i) + 1, indj = star(j) + 1; 
          // change statistics of edge
          changestat(0) = ( Choose(indi,1) + Choose(indj,1) - Choose(star(i),1) - Choose(star(j),1) )/2;
          
          // change statistics of nodematch
          changestat(1) = ( (grade[i]==7)*(grade[j]==7)   );
          changestat(2) = ( (grade[i]==8)*(grade[j]==8)   );
          changestat(3) = ( (grade[i]==9)*(grade[j]==9)   );	       
          changestat(4) = ( (grade[i]==10)*(grade[j]==10)   );
          changestat(5) = ( (grade[i]==11)*(grade[j]==11)   );            
          changestat(6) = ( (grade[i]==12)*(grade[j]==12)   );       		               
          changestat(7) = (  sex[i]==sex[j]  );  		               
          
          // change statistics of gwd
          changestat(8) = exp(decay)*( geoweight(indi-1) + geoweight(indj-1) );
          if(indi-1>0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indi-1-1); }
          if(indj-1>0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indj-1-1); }
          
          // change statistics of gwes
          changestat(9) = 0;     
          for(int k = 0; k< nrow; k++){
            if(   ( X(i,k) == 1 ) & ( X(j,k) == 1  )   ){
              prei = sum( X.row(i)%X.row(k) );
              prej = sum( X.row(j)%X.row(k) );
              
              changestat(9) = changestat(9) + exp(decay)*geoweight(prei+1-1) ; 
              if(prei >0){ changestat(9) = changestat(9) - exp(decay)*geoweight(prei-1);} 
              changestat(9) = changestat(9) + exp(decay)*geoweight(prej+1-1) ; 
              if(prej >0){ changestat(9) = changestat(9) - exp(decay)*geoweight(prej-1);} 
              res = res + 1; // X(k,i) multiply X(i,j) +1 multiply X(j,k)
            }
          }
          if(res > 0){ changestat(9) = changestat(9) + exp(decay)*geoweight(res-1); }
          
          // accept reject step  	       
          // probability of Xij = 1 for given others are fixed 
          vec r =  exp(trans(coef)*changestat) ;         
          double p = r(0)/(1+r(0));   		   
          if( randu() < p  ){
            X(i,j) = X(j,i) = 1; 
            star(i) = indi; star(j) = indj;		 
            Sumstat = Sumstat + changestat;
          }
          
          
          // When Xij is 1  
        }else{
          indi = star(i) - 1, indj = star(j) - 1;	
          // change statistics of edge
          changestat(0) = ( Choose(star(i),1) + Choose(star(j),1) - Choose(indi,1) - Choose(indj,1) )/2;
          
          // change statistics of nodematch
          changestat(1) = ( (grade[i]==7)*(grade[j]==7)   );
          changestat(2) = ( (grade[i]==8)*(grade[j]==8)   );
          changestat(3) = ( (grade[i]==9)*(grade[j]==9)   );	       
          changestat(4) = ( (grade[i]==10)*(grade[j]==10)   );
          changestat(5) = ( (grade[i]==11)*(grade[j]==11)   );           
          changestat(6) = ( (grade[i]==12)*(grade[j]==12)   );     		              
          changestat(7) = (  sex[i]==sex[j]  );     		              
          
          // change statistics of gwd
          changestat(8) = exp(decay)*( geoweight(indi-1+1) + geoweight(indj-1+1)  );
          if(indi-1>=0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indi-1); }       
          if(indj-1>=0){ changestat(8) = changestat(8) - exp(decay)*geoweight(indj-1); }        		   
          
          // change statistics of gwesp 
          changestat(9) = 0;
          for(int k = 0; k< nrow; k++){
            if(   ( X(i,k) == 1 ) & ( X(j,k) == 1  )   ){
              prei = sum( X.row(i)%X.row(k) );
              prej = sum( X.row(j)%X.row(k) );
              if(prei-1 >0){changestat(9) = changestat(9) - exp(decay)*geoweight(prei-1-1); } 
              changestat(9) = changestat(9) + exp(decay)*geoweight(prei-1); 
              if(prej-1 >0){changestat(9) = changestat(9) - exp(decay)*geoweight(prej-1-1); } 
              changestat(9) = changestat(9) + exp(decay)*geoweight(prej-1); 
              res = res + 1; // X(k,i) multiply X(i,j)  multiply X(j,k)
            }
          }
          if(res > 0){ changestat(9) = changestat(9) + exp(decay)*geoweight(res-1); }
          
          // accept reject step  	       
          // probability of Xij = 1 for given others are fixed 
          vec r =  exp(trans(coef)*changestat) ;         
          double p = r(0)/(1+r(0));   		   
          if( randu() > p  ){
            X(i,j) = X(j,i) = 0; 
            star(i) = indi; star(j) = indj;  
            Sumstat = Sumstat - changestat;
          }             
        }
        // End of a single update  
        ivec(i) = 0; jvec(j) = 0;
      }
    } 
    allinner=allinner+Sumstat;
  }
  allinner=allinner/cycle;
  return(allinner); 
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Double metropolis hastings algorithm
mat ergmDMH(mat X, vec grade, vec sex, mat COV, mat theta, int outer, int cycle){
  
  // Initializing part
  double logprob,u;                          // used in Outer MCMC
  int nCOVcols = COV.n_cols;                 // number of parameters	
  vec thetaprev(nCOVcols);                   // before propose in Outer MCMC, previous parameters
  vec stat = Summary(X, grade,sex), statprop(nCOVcols); // sufficient statistics
  
  //// Start of OUTER MCMC Chain 
  for(int l = 0; l< outer; l++){
    
    //	if( (l > 1000) && (l <= 10000) ){ // adaptively update COV until 10000 iterations 
    //	COV = cov(theta);
    //    }	
    
    for(int i = 0; i< nCOVcols; i++){
      thetaprev[i] = theta(l,i);
    }
    
    vec Znormal = randn(nCOVcols);                                           // multivariate proposal by using Cholesky factorization
    vec thetaprop = trans(  trans(thetaprev) + trans(Znormal)*chol(COV)  );  // proposed parameter
    
    // proposed auxiliary variable		
    statprop = Gibbs(X, grade, sex, thetaprop, cycle);
    
    // log probability ratio to determine acceptance of Outer MCMC 	   
    vec dummy = ( -0.05*trans(thetaprop)*thetaprop + 0.05*trans(thetaprev)*thetaprev + trans(thetaprev - thetaprop)*(statprop - stat) );
    logprob = dummy[0];
    u = log( randu() );
    if( u< logprob ){
      theta.insert_rows(l+1,trans(thetaprop));
    }else{
      theta.insert_rows(l+1,trans(thetaprev));
    }
    
  }
  
  return(theta);	
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Double metropolis hastings algorithm
mat ergmDMH_update(mat X, vec grade, vec sex, mat COV, mat theta, int outer, int cycle, int adaptInterval, double adaptFactorExponent, int adapIter,
                   int thin){
  
  // Initializing part
  double logprob,u;                         // used in Outer MCMC
  double rhat = 0, gamma1 = 0, gamma2 = 0;
  vec accprob = zeros(outer);
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  int nCOVcols = COV.n_cols;                 // number of parameters	
  vec thetaprev(nCOVcols);                   // before propose in Outer MCMC, previous parameters
  vec stat = Summary(X, grade,sex), statprop(nCOVcols); // sufficient statistics
  double sigma2 = (5.76)/nCOVcols;
  
  mat cholCOV=COV;
  int mj = nCOVcols;
  cholCOV = trans( chol( sigma2 * (COV + 0.001 * diagmat(ones(mj)) ) )) ;
  
  
  //// Start of OUTER MCMC Chain 
  for(int l = 0; l< outer; l++){
    
    //	if( (l > 1000) && (l <= 10000) ){ // adaptively update COV until 10000 iterations 
    //	COV = cov(theta);
    //    }	
    if( (l+1 >= adaptInterval) && (l+1 - (adaptInterval * trunc((l+1) / adaptInterval)) == 0) ){
      double dummyacc=0;
      for(int acc=l+1-adaptInterval; acc<l; acc++){
        dummyacc = dummyacc + accprob(acc);
      }
      rhat =  dummyacc / (adaptInterval-1);
      gamma1 = 1 / pow(adapIter, c1);
      gamma2 = c0 * gamma1;
      sigma2 = exp( log(sigma2) + gamma2 * (rhat - ropt) );
      
      COV = COV + gamma1 * ( cov( theta.rows(l+1-adaptInterval, l-1) ) - COV );
      cholCOV = trans( chol( sigma2 * (COV + 0.001 * diagmat(ones(mj)) ) )) ;
      
      adapIter = adapIter + 1;
    }
    
    
    for(int i = 0; i< nCOVcols; i++){
      thetaprev[i] = theta(l,i);
    }
    
    vec Znormal = randn(nCOVcols);                                           // multivariate proposal by using Cholesky factorization
    vec thetaprop = trans(  trans(thetaprev) + trans(Znormal)*cholCOV  );  // proposed parameter
    
    // proposed auxiliary variable		
    statprop = Gibbs(X, grade, sex, thetaprop, cycle);
    
    // log probability ratio to determine acceptance of Outer MCMC 	   
    vec dummy = ( -0.05*trans(thetaprop)*thetaprop + 0.05*trans(thetaprev)*thetaprev + trans(thetaprev - thetaprop)*(statprop - stat) );
    logprob = dummy[0];
    u = log( randu() );
    if( u< logprob ){
      theta.insert_rows(l+1,trans(thetaprop));
      accprob(l)=1;
    }else{
      theta.insert_rows(l+1,trans(thetaprev));
    }
    
  }
  
  return(theta);	
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Simulate pcf on grid and measure distance from data points 
List pAuxgen(mat X, vec grade, vec sex, int cycle, int num, int numcore, mat Designmat){
  
  omp_set_num_threads(numcore);
  int numpoint = Designmat.n_rows;       // number of designpoints
  int ndim = Designmat.n_cols;           // parameter dimension
  
  mat Mu(numpoint, ndim);
  mat Cov(numpoint*ndim, ndim);
  
  int i;
#pragma omp parallel shared(Designmat) private(i)
{	
#pragma omp for schedule(static)  
  for(i = 0; i < numpoint; i++){
    mat pseudo(num,ndim);
    for(int j = 0; j < num; j++){ pseudo.row(j) =  trans( Gibbs(X, grade, sex, trans( Designmat.row( i ) ), cycle) ); }
    
    Mu.row(i) =  mean(pseudo,0) ;    
    Cov.submat( span(ndim*i, ndim*i+ndim-1), span(0, ndim-1) ) = cov(pseudo);
    
    
  }
}

return List::create(Named("Cov") = Cov, Named("Mu") = Mu);	
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// generate cross distance matrix
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
/// product with each element in matrix
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
// power on each element of matrix
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
// taking exp on each element of matrix
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
// taking inverse on each element of matrix
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



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// RBF kernel
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
mat ergmvgdc2(mat x0, mat X, vec grade, vec sex, int niter, int inner, int imporsize, double stepsize, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repsum(summary,Summary(X, grade, sex));
  
  for(int iter=0; iter<niter; iter++){
    mat suf(theta.n_cols, theta.n_rows);
    for(int ii=0; ii<theta.n_rows; ii++){
      vec vectheta = vectorise(theta.row(ii));
      for(int insize=0;  insize<imporsize; insize++){
        suf.col(ii)=suf.col(ii)+(Gibbs(X,grade, sex, vectheta,inner));
      }
    }    
    
    suf = suf/imporsize;
    mat lnpgrad = summary - suf - 0.1*theta.t();
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
  return(theta);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
mat ergmispar2(mat x0, mat X, vec grade, vec sex, int niter, int inner, int imporsize, double stepsize, int core, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repsum(summary,Summary(X, grade, sex));
  int tt = theta.n_rows/core;
  
  for(int iter=0; iter<niter; iter++){
    mat suf(theta.n_cols, theta.n_rows);
    
    int jj;
#pragma omp parallel shared(theta) private(jj)
{
#pragma omp for schedule(static)
  for(jj=0; jj<core; jj++){
    for(int ii=0; ii<tt; ii++){
      vec vectheta = vectorise(theta.row((tt*jj)+ii));
      for(int insize=0;  insize<imporsize; insize++){
        suf.col((tt*jj)+ii)=suf.col((tt*jj)+ii)+(Gibbs(X,grade, sex, vectheta,inner));
      }
    }
  }
}

suf = suf/imporsize;
mat lnpgrad = summary - suf - 0.1*theta.t();
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
  return(theta);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
mat ergmispar3(mat x0, mat X, vec grade, vec sex, int niter, int inner, int imporsize, double stepsize, int core, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  mat alltheta(theta.n_rows*niter, theta.n_cols);
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repsum(summary,Summary(X, grade, sex));
  int tt = theta.n_rows/core;
  int numthe = theta.n_rows;
  
  for(int iter=0; iter<niter; iter++){
    mat suf(theta.n_cols, theta.n_rows);
    
    int jj;
#pragma omp parallel shared(theta) private(jj)
{
#pragma omp for schedule(static)
  for(jj=0; jj<core; jj++){
    for(int ii=0; ii<tt; ii++){
      vec vectheta = vectorise(theta.row((tt*jj)+ii));
      for(int insize=0;  insize<imporsize; insize++){
        suf.col((tt*jj)+ii)=suf.col((tt*jj)+ii)+(Gibbs(X,grade, sex, vectheta,inner));
      }
    }
  }
}

suf = suf/imporsize;
mat lnpgrad = summary - suf - 0.1*theta.t();
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

for(int jjj=0; jjj<numthe; jjj++){
  alltheta.row(jjj+(iter*numthe)) = theta.row(jjj);
}
  }
  return(alltheta);
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
mat ergmispar4(mat x0, mat X, vec grade, vec sex, int niter, int inner, int imporsize, double stepsize, int core, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repsum(summary,Summary(X, grade, sex));
  int tt = theta.n_rows/core;
  
  for(int iter=0; iter<niter; iter++){
    mat suf(theta.n_cols, theta.n_rows);
    
    int jj;
#pragma omp parallel shared(theta) private(jj)
{
#pragma omp for schedule(static)
  for(jj=0; jj<core; jj++){
    for(int ii=0; ii<tt; ii++){
      vec vectheta = vectorise(theta.row((tt*jj)+ii));
      for(int insize=0;  insize<imporsize; insize++){
        suf.col((tt*jj)+ii)=suf.col((tt*jj)+ii)+(Gibbsall(X,grade, sex, vectheta,inner));
      }
    }
  }
}

suf = suf/imporsize;
mat lnpgrad = summary - suf - 0.1*theta.t();
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
  return(theta);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// get the MAP
mat ergm_map(mat x0, mat X, vec grade, vec sex, int niter, int inner, int imporsize, double stepsize, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repsum(summary,Summary(X, grade, sex));
  mat alltheta(theta.n_rows*niter, theta.n_cols);
  int numthe = theta.n_rows;
  
  for(int iter=0; iter<niter; iter++){
    mat suf(theta.n_cols, theta.n_rows);
    for(int ii=0; ii<theta.n_rows; ii++){
      vec vectheta = vectorise(theta.row(ii));
      for(int insize=0;  insize<imporsize; insize++){
        suf.col(ii)=suf.col(ii)+(Gibbs(X,grade,sex,vectheta,inner));
      }
    }    
    
    suf = suf/imporsize;
    mat lnpgrad = summary - suf - 0.1*theta.t();
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
    for(int jjj=0; jjj<numthe; jjj++){
      alltheta.row(jjj+(iter*numthe)) = theta.row(jjj);
    }
  }
  return(alltheta);
}






// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List ergm_mesa_svi_par_all(mat x0, mat X, vec grade, vec sex, mat map, int niter, int inner, int imporsize, int imporsizes, double stepsize, double thres, int core, int burn, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repsum(summary,Summary(X, grade, sex));
  mat zerotheta = zeros(theta.n_cols);
  int percore = theta.n_rows/core;
  vec mc = zeros(theta.n_rows);
  mat alltheta(theta.n_rows*niter, theta.n_cols);
  mat allscore(theta.n_rows*niter, theta.n_cols);
  int numthe = theta.n_rows;
  vec nummc = zeros(niter);
  
  mat suf(imporsize*inner, theta.n_cols);
  mat thetas = map.t();
  
  for(int insize=0;  insize<imporsize; insize++){
    mat gib = Gibbs2(X,grade, sex, map,inner+burn);
    gib = gib.t();
    for(int innersize=0; innersize<inner; innersize++){
      suf.row(insize*inner+innersize)=gib.row(innersize+burn);
    }
  }
  
  
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums(theta.n_rows,theta.n_cols);
  mat suftheta = suf;
  
  for(int iter=0; iter<niter; iter++){
    mc = zeros(theta.n_rows);
    for(int tt=0; tt<theta.n_rows; tt++){
      double distmap = crossdist2(theta.row(tt), map.t())[0];
      double distmin = min( crossdist2(theta.row(tt), thetas) ); 
      if(distmap > distmin) {
        int minwhich = vecminInd(crossdist2(theta.row(tt), thetas));
        mat selip(imporsizes*inner, theta.n_cols);
        mat selip2(imporsizes*inner, theta.n_cols);
        mat norwei(imporsizes*inner, theta.n_cols);
        for(int jj=0; jj<imporsizes*inner; jj++){
          selip.row(jj)=mmult((theta.row(tt)-thetas.row(minwhich)), suf.row(imporsize*inner+(minwhich-1)*(imporsizes*inner)+jj));
        }
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsizes*inner; jj++){
          selip2.row(jj)=mmult((expselip.row(jj)/selsum.t()), suf.row(imporsize*inner+(minwhich-1)*(imporsizes*inner)+jj));
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
        
        
        if(min(midess.col(tt))<(imporsizes*inner/thres)){
          mc[tt]=1;
        }
      }else{
        mat selip(imporsize*inner, theta.n_cols);
        mat selip2(imporsize*inner, theta.n_cols);
        mat norwei(imporsize*inner, theta.n_cols);
        for(int jj=0; jj<imporsize*inner; jj++){
          selip.row(jj)=mmult((theta.row(tt)-map.t()), suf.row(jj));
        }
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsize*inner; jj++){
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
        
        if(min(midess.col(tt))<(imporsize*inner/thres)){
          mc[tt]=1;
        }
      }
    }  
    
    mat suf2(theta.n_rows*inner*imporsizes, theta.n_cols);
    
    
    
    int cor;
#pragma omp parallel shared(theta) private(cor)
{
#pragma omp for schedule(static)
  for(cor=0; cor<core; cor++){
    for(int tt2=0; tt2<percore; tt2++){
      if(mc[percore*cor+tt2]==1){
        vec vectheta = vectorise(theta.row(percore*cor+tt2));
        parsum.row(percore*cor+tt2)=zerotheta.t();
        for(int insizes=0; insizes<imporsizes; insizes++){
          mat gibs=Gibbs2(X,grade, sex, vectheta,inner+burn);
          gibs = gibs.t();
          for(int innersize=0; innersize<inner; innersize++){
            parsum.row(percore*cor+tt2)=parsum.row(percore*cor+tt2)+gibs.row(innersize+burn);
            suf2.row(((percore*cor+tt2)*inner*imporsizes)+(insizes*inner+innersize))=gibs.row(innersize+burn);
          }
          parsum.row(percore*cor+tt2) = parsum.row(percore*cor+tt2)/inner;
        }
        parsum.row(percore*cor+tt2) = parsum.row(percore*cor+tt2)/imporsizes; 
      }
    }
  }
}

for(int ttt=0; ttt<theta.n_rows; ttt++){    
  if(mc[ttt]==1){
    thetas = rbind_mat(thetas, theta.row(ttt));
    for(int iii=0; iii<inner*imporsizes; iii++){
      suf = rbind_mat(suf, suf2.row(ttt*inner*imporsizes+iii));
    }
  }
}
nummc[iter] = sum(mc);
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
for(int jjj=0; jjj<numthe; jjj++){
  alltheta.row(jjj+(iter*numthe)) = theta.row(jjj);
}

for(int jjjj=0; jjjj<numthe; jjjj++){
  mat lnp = lnpgrad.t();
  allscore.row(jjjj+(iter*numthe)) = lnp.row(jjjj);
}

  }
  return List::create(Named("theta")=theta, Named("alltheta")=alltheta, Named("allscore")=allscore, Named("thetas")=thetas, Named("suf")=suf, Named("nummc")=nummc);
}






// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List ergm_mesa_svi_par_all2(mat x0, mat X, vec grade, vec sex, mat map, int niter, int inner, int imporsize, int imporsizes, double stepsize, double thres, int core, int burn, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repsum(summary,Summary(X, grade, sex));
  mat zerotheta = zeros(theta.n_cols);
  int percore = theta.n_rows/core;
  vec mc = zeros(theta.n_rows);
  mat alltheta(theta.n_rows*niter, theta.n_cols);
  int numthe = theta.n_rows;
  vec nummc = zeros(niter);
  
  mat suf(imporsize, theta.n_cols);
  mat thetas = map.t();
  
  mat gib = Gibbs3(X,grade, sex, map,inner+burn, imporsize);
  gib = gib.t();
  for(int innersize=0; innersize<imporsize; innersize++){
    suf.row(innersize)=gib.row(innersize);
  }
  
  
  
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums(theta.n_rows,theta.n_cols);
  mat suftheta = suf;
  
  for(int iter=0; iter<niter; iter++){
    mc = zeros(theta.n_rows);
    for(int tt=0; tt<theta.n_rows; tt++){
      double distmap = crossdist2(theta.row(tt), map.t())[0];
      double distmin = min( crossdist2(theta.row(tt), thetas) ); 
      if(distmap > distmin) {
        int minwhich = vecminInd(crossdist2(theta.row(tt), thetas));
        mat selip(imporsizes, theta.n_cols);
        mat selip2(imporsizes, theta.n_cols);
        mat norwei(imporsizes, theta.n_cols);
        for(int jj=0; jj<imporsizes; jj++){
          selip.row(jj)=mmult((theta.row(tt)-thetas.row(minwhich)), suf.row(imporsize+(minwhich-1)*(imporsizes)+jj));
        }
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsizes; jj++){
          selip2.row(jj)=mmult((expselip.row(jj)/selsum.t()), suf.row(imporsize+(minwhich-1)*(imporsizes)+jj));
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
        
        
        if(min(midess.col(tt))<(imporsizes/thres)){
          mc[tt]=1;
        }
      }else{
        mat selip(imporsize, theta.n_cols);
        mat selip2(imporsize, theta.n_cols);
        mat norwei(imporsize, theta.n_cols);
        for(int jj=0; jj<imporsize; jj++){
          selip.row(jj)=mmult((theta.row(tt)-map.t()), suf.row(jj));
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
        
        if(min(midess.col(tt))<(imporsize/thres)){
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
        mat gibs=Gibbs3(X,grade, sex, vectheta,inner+burn, imporsize);
        gibs = gibs.t();
        for(int innersize=0; innersize<imporsize; innersize++){
          parsum.row(percore*cor+tt2)=parsum.row(percore*cor+tt2)+gibs.row(innersize);
          suf2.row(((percore*cor+tt2)*imporsizes)+(innersize))=gibs.row(innersize);
        }
        parsum.row(percore*cor+tt2) = parsum.row(percore*cor+tt2)/imporsizes; 
      }
    }
  }
}

for(int ttt=0; ttt<theta.n_rows; ttt++){    
  if(mc[ttt]==1){
    thetas = rbind_mat(thetas, theta.row(ttt));
    for(int iii=0; iii<imporsizes; iii++){
      suf = rbind_mat(suf, suf2.row(ttt*imporsizes+iii));
    }
  }
}
nummc[iter] = sum(mc);
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
for(int jjj=0; jjj<numthe; jjj++){
  alltheta.row(jjj+(iter*numthe)) = theta.row(jjj);
}

  }
  return List::create(Named("theta")=theta, Named("alltheta")=alltheta, Named("thetas")=thetas, Named("nummc")=nummc);
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List ergm_mesa_svi_par_all3(mat x0, mat X, vec grade, vec sex, mat map, int niter, int inner, int imporsize, int imporsizes, double stepsize, double thres, int core, int burn, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repsum(summary,Summary(X, grade, sex));
  mat zerotheta = zeros(theta.n_cols);
  int percore = theta.n_rows/core;
  vec mc = zeros(theta.n_rows);
  mat alltheta(theta.n_rows*niter, theta.n_cols);
  int numthe = theta.n_rows;
  
  mat suf(imporsize, theta.n_cols);
  mat thetas = map.t();
  
  mat gib = Gibbs4(X,grade, sex, map,inner+burn, imporsize);
  gib = gib.t();
  for(int innersize=0; innersize<imporsize; innersize++){
    suf.row(innersize)=gib.row(innersize);
  }
  
  
  
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums(theta.n_rows,theta.n_cols);
  mat suftheta = suf;
  
  for(int iter=0; iter<niter; iter++){
    mc = zeros(theta.n_rows);
    for(int tt=0; tt<theta.n_rows; tt++){
      double distmap = crossdist2(theta.row(tt), map.t())[0];
      double distmin = min( crossdist2(theta.row(tt), thetas) ); 
      if(distmap > distmin) {
        int minwhich = vecminInd(crossdist2(theta.row(tt), thetas));
        mat selip(imporsizes, theta.n_cols);
        mat selip2(imporsizes, theta.n_cols);
        mat norwei(imporsizes, theta.n_cols);
        for(int jj=0; jj<imporsizes; jj++){
          selip.row(jj)=mmult((theta.row(tt)-thetas.row(minwhich)), suf.row(imporsize+(minwhich-1)*(imporsizes)+jj));
        }
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsizes; jj++){
          selip2.row(jj)=mmult((expselip.row(jj)/selsum.t()), suf.row(imporsize+(minwhich-1)*(imporsizes)+jj));
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
        
        
        if(min(midess.col(tt))<(imporsizes/thres)){
          mc[tt]=1;
        }
      }else{
        mat selip(imporsize, theta.n_cols);
        mat selip2(imporsize, theta.n_cols);
        mat norwei(imporsize, theta.n_cols);
        for(int jj=0; jj<imporsize; jj++){
          selip.row(jj)=mmult((theta.row(tt)-map.t()), suf.row(jj));
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
        
        if(min(midess.col(tt))<(imporsize/thres)){
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
        mat gibs=Gibbs4(X,grade, sex, vectheta,inner+burn, imporsize);
        gibs = gibs.t();
        for(int innersize=0; innersize<imporsize; innersize++){
          parsum.row(percore*cor+tt2)=parsum.row(percore*cor+tt2)+gibs.row(innersize);
          suf2.row(((percore*cor+tt2)*imporsizes)+(innersize))=gibs.row(innersize);
        }
        parsum.row(percore*cor+tt2) = parsum.row(percore*cor+tt2)/imporsizes; 
      }
    }
  }
}

for(int ttt=0; ttt<theta.n_rows; ttt++){    
  if(mc[ttt]==1){
    thetas = rbind_mat(thetas, theta.row(ttt));
    for(int iii=0; iii<imporsizes; iii++){
      suf = rbind_mat(suf, suf2.row(ttt*imporsizes+iii));
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
for(int jjj=0; jjj<numthe; jjj++){
  alltheta.row(jjj+(iter*numthe)) = theta.row(jjj);
}

  }
  return List::create(Named("theta")=theta, Named("alltheta")=alltheta, Named("thetas")=thetas, Named("suf")=suf);
}











// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List ergm_mesa_svi_par_all4(mat x0, mat X, vec grade, vec sex, mat map, int niter, int inner, int imporsize, int imporsizes, double stepsize, double thres, int core, int burn, int cyc, int thin, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repsum(summary,Summary(X, grade, sex));
  mat zerotheta = zeros(theta.n_cols);
  int percore = theta.n_rows/core;
  vec mc = zeros(theta.n_rows);
  int numthe = theta.n_rows;
  vec nummc = zeros(niter);
  
  mat suf(imporsize, theta.n_cols);
  mat thetas = map.t();
  
  mat gib = Gibbs5(X,grade, sex, map,inner+burn, imporsize, cyc, thin);
  gib = gib.t();
  for(int innersize=0; innersize<imporsize; innersize++){
    suf.row(innersize)=gib.row(innersize);
  }
  
  
  
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums(theta.n_rows,theta.n_cols);
  mat suftheta = suf;
  
  for(int iter=0; iter<niter; iter++){
    mc = zeros(theta.n_rows);
    for(int tt=0; tt<theta.n_rows; tt++){
      double distmap = crossdist2(theta.row(tt), map.t())[0];
      double distmin = min( crossdist2(theta.row(tt), thetas) ); 
      if(distmap > distmin) {
        int minwhich = vecminInd(crossdist2(theta.row(tt), thetas));
        mat selip(imporsizes, theta.n_cols);
        mat selip2(imporsizes, theta.n_cols);
        mat norwei(imporsizes, theta.n_cols);
        for(int jj=0; jj<imporsizes; jj++){
          selip.row(jj)=mmult((theta.row(tt)-thetas.row(minwhich)), suf.row(imporsize+(minwhich-1)*(imporsizes)+jj));
        }
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsizes; jj++){
          selip2.row(jj)=mmult((expselip.row(jj)/selsum.t()), suf.row(imporsize+(minwhich-1)*(imporsizes)+jj));
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
        
        
        if(min(midess.col(tt))<(imporsizes/thres)){
          mc[tt]=1;
        }
      }else{
        mat selip(imporsize, theta.n_cols);
        mat selip2(imporsize, theta.n_cols);
        mat norwei(imporsize, theta.n_cols);
        for(int jj=0; jj<imporsize; jj++){
          selip.row(jj)=mmult((theta.row(tt)-map.t()), suf.row(jj));
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
        
        if(min(midess.col(tt))<(imporsize/thres)){
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
        mat gibs=Gibbs5(X,grade, sex, vectheta,inner+burn, imporsize, cyc, thin);
        gibs = gibs.t();
        for(int innersize=0; innersize<imporsize; innersize++){
          parsum.row(percore*cor+tt2)=parsum.row(percore*cor+tt2)+gibs.row(innersize);
          suf2.row(((percore*cor+tt2)*imporsizes)+(innersize))=gibs.row(innersize);
        }
        parsum.row(percore*cor+tt2) = parsum.row(percore*cor+tt2)/imporsizes; 
      }
    }
  }
}

for(int ttt=0; ttt<theta.n_rows; ttt++){    
  if(mc[ttt]==1){
    thetas = rbind_mat(thetas, theta.row(ttt));
    for(int iii=0; iii<imporsizes; iii++){
      suf = rbind_mat(suf, suf2.row(ttt*imporsizes+iii));
    }
  }
}
nummc[iter] = sum(mc);
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
  return List::create(Named("theta")=theta, Named("thetas")=thetas, Named("nummc")=nummc);
}




















// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List ergm_mesa_svi_par_alls(mat x0, mat X, vec grade, vec sex, mat map, int niter, int inner, int imporsize, int imporsizes, double stepsize, double thres, int core, int burn, int cyc, int thin, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repsum(summary,Summary(X, grade, sex));
  mat zerotheta = zeros(theta.n_cols);
  int percore = theta.n_rows/core;
  vec mc = zeros(theta.n_rows);
  int numthe = theta.n_rows;
  vec nummc = zeros(niter);
  
  mat suf(imporsize, theta.n_cols);
  mat thetas = map.t();
  
  mat gib = Gibbs5(X,grade, sex, map,inner+burn, imporsize, cyc, thin);
  gib = gib.t();
  for(int innersize=0; innersize<imporsize; innersize++){
    suf.row(innersize)=gib.row(innersize);
  }
  
  
  
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums(theta.n_rows,theta.n_cols);
  mat suftheta = suf;
  
  for(int iter=0; iter<niter; iter++){
    mc = zeros(theta.n_rows);
    for(int tt=0; tt<theta.n_rows; tt++){
      double distmap = crossdist2(theta.row(tt), map.t())[0];
      double distmin = min( crossdist2(theta.row(tt), thetas) ); 
      if(distmap > distmin) {
        int minwhich = vecminInd(crossdist2(theta.row(tt), thetas));
        mat selip(imporsizes, theta.n_cols);
        mat selip2(imporsizes, theta.n_cols);
        mat norwei(imporsizes, theta.n_cols);
        for(int jj=0; jj<imporsizes; jj++){
          selip.row(jj)=mmult((theta.row(tt)-thetas.row(minwhich)), suf.row(imporsize+(minwhich-1)*(imporsizes)+jj));
        }
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsizes; jj++){
          selip2.row(jj)=mmult((expselip.row(jj)/selsum.t()), suf.row(imporsize+(minwhich-1)*(imporsizes)+jj));
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
        
        
        if(min(midess.col(tt))<(imporsizes/thres)){
          mc[tt]=1;
        }
      }else{
        mat selip(imporsize, theta.n_cols);
        mat selip2(imporsize, theta.n_cols);
        mat norwei(imporsize, theta.n_cols);
        for(int jj=0; jj<imporsize; jj++){
          selip.row(jj)=mmult((theta.row(tt)-map.t()), suf.row(jj));
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
        
        if(min(midess.col(tt))<(imporsize/thres)){
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
        
        mat gibs = Gibbs5(X, grade, sex, vectheta, inner + burn, imporsize, cyc, thin);
        gibs = gibs.t();  

        for (int innersize = 0; innersize < imporsize; ++innersize) {
          parsum.row(idx) += gibs.row(innersize);
          suf2.row((idx * imporsizes) + innersize) = gibs.row(innersize);  
        }

        parsum.row(idx) /= imporsize;
      }
    }
  }
}

for(int ttt=0; ttt<theta.n_rows; ttt++){    
  if(mc[ttt]==1){
    thetas = rbind_mat(thetas, theta.row(ttt));
    for(int iii=0; iii<imporsizes; iii++){
      suf = rbind_mat(suf, suf2.row(ttt*imporsizes+iii));
    }
  }
}
nummc[iter] = sum(mc);
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
  return List::create(Named("theta")=theta, Named("thetas")=thetas, Named("nummc")=nummc);
}
















// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List ergm_weight_mcsvgd(mat x0, mat X, vec grade, vec sex, mat map, int niter, int inner, int imporsize, double close, double stepsize, double thres, int core, int burn, int prerun, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repsum(summary,Summary(X, grade, sex));
  mat zerotheta = zeros(theta.n_cols);
  int percore = theta.n_rows/core;
  vec mc = zeros(theta.n_rows);
  mat alltheta(theta.n_rows*niter, theta.n_cols);
  mat allscore(theta.n_rows*niter, theta.n_cols);
  int numthe = theta.n_rows;
  vec nummc = zeros(niter);
  
  mat suf(imporsize*inner, theta.n_cols);
  mat thetas = map.t();
  
  for(int insize=0;  insize<imporsize; insize++){
    mat gib = Gibbs2(X,grade, sex, map,inner+burn);
    gib = gib.t();
    for(int innersize=0; innersize<inner; innersize++){
      suf.row(insize*inner+innersize)=gib.row(innersize+burn);
    }
  }
  
  
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums(theta.n_rows,theta.n_cols);
  mat suftheta = suf;
  
  for(int iter=0; iter<niter; iter++){
    mc = zeros(theta.n_rows);
    if(thetas.n_rows>prerun){
      for(int tt=0; tt<theta.n_rows; tt++){
        vec distval = crossdist2(theta.row(tt), thetas);
        vec minwhich = find_lowest_indices(distval, close);
        vec weights(close);
        weights.fill(1.0/close);
        mat selip2(imporsize*inner*close, theta.n_cols);
        mat norwei(imporsize*inner*close, theta.n_cols);
        mat selip(imporsize*inner*close, theta.n_cols);
        for(int near=0; near<close; near++){
          for(int jj=0; jj<imporsize*inner; jj++){
            selip.row((imporsize*inner)*near+jj)=mmult((theta.row(tt)-thetas.row(minwhich[near])), suf.row((minwhich[near])*(imporsize*inner)+jj));
          }
        }
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int near=0; near<close; near++){
          for(int jj=0; jj<imporsize*inner; jj++){
            selip2.row((imporsize*inner)*near+jj)=mmult((expselip.row((imporsize*inner)*near+jj)/selsum.t()), suf.row((minwhich[near])*(imporsize*inner)+jj));
            mat nor = (expselip.row((imporsize*inner)*near+jj)/selsum.t());
            norwei.row((imporsize*inner)*near+jj)= nor;
          }
        }
        vec selsum2 = colSum(selip2);
        norwei = mpow(norwei,2); 
        vec norweisum = colSum(norwei);
        norweisums.row(tt)=norweisum.t();
        parsum.row(tt)=selsum2.t();
        
        ess = norweisums;
        mat midess=minv(ess);
        midess = midess.t();
        
        
        if(min(midess.col(tt))<(thres)){
          mc[tt]=1;
        }
      }  
    }else{
      for(int tt=0; tt<theta.n_rows; tt++){
        int minwhich = vecminInd(crossdist2(theta.row(tt), thetas));
        mat selip(imporsize*inner, theta.n_cols);
        mat selip2(imporsize*inner, theta.n_cols);
        mat norwei(imporsize*inner, theta.n_cols);
        for(int jj=0; jj<imporsize*inner; jj++){
          selip.row(jj)=mmult((theta.row(tt)-thetas.row(minwhich)), suf.row((minwhich)*(imporsize*inner)+jj));
        }
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsize*inner; jj++){
          selip2.row(jj)=mmult((expselip.row(jj)/selsum.t()), suf.row((minwhich)*(imporsize*inner)+jj));
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
        
        
        if(min(midess.col(tt))<(thres)){
          mc[tt]=1;
        }
      }  
    }
    mat suf2(theta.n_rows*inner*imporsize, theta.n_cols);
    
    
    
    int cor;
#pragma omp parallel shared(theta) private(cor)
{
#pragma omp for schedule(static)
  for(cor=0; cor<core; cor++){
    for(int tt2=0; tt2<percore; tt2++){
      if(mc[percore*cor+tt2]==1){
        vec vectheta = vectorise(theta.row(percore*cor+tt2));
        parsum.row(percore*cor+tt2)=zerotheta.t();
        for(int insizes=0; insizes<imporsize; insizes++){
          mat gibs=Gibbs2(X,grade, sex, vectheta,inner+burn);
          gibs = gibs.t();
          for(int innersize=0; innersize<inner; innersize++){
            parsum.row(percore*cor+tt2)=parsum.row(percore*cor+tt2)+gibs.row(innersize+burn);
            suf2.row(((percore*cor+tt2)*inner*imporsize)+(insizes*inner+innersize))=gibs.row(innersize+burn);
          }
          parsum.row(percore*cor+tt2) = parsum.row(percore*cor+tt2)/inner;
        }
        parsum.row(percore*cor+tt2) = parsum.row(percore*cor+tt2)/imporsize; 
      }
    }
  }
}

for(int ttt=0; ttt<theta.n_rows; ttt++){    
  if(mc[ttt]==1){
    thetas = rbind_mat(thetas, theta.row(ttt));
    for(int iii=0; iii<inner*imporsize; iii++){
      suf = rbind_mat(suf, suf2.row(ttt*inner*imporsize+iii));
    }
  }
}
nummc[iter] = sum(mc);
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
for(int jjj=0; jjj<numthe; jjj++){
  alltheta.row(jjj+(iter*numthe)) = theta.row(jjj);
}

for(int jjjj=0; jjjj<numthe; jjjj++){
  mat lnp = lnpgrad.t();
  allscore.row(jjjj+(iter*numthe)) = lnp.row(jjjj);
}

  }
  return List::create(Named("theta")=theta, Named("alltheta")=alltheta, Named("allscore")=allscore, Named("thetas")=thetas, Named("suf")=suf, Named("nummc")=nummc);
}









// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List ergm_uniweight_mcsvgd(mat x0, mat X, vec grade, vec sex, mat map, int niter, int inner, int imporsize, double close, double stepsize, double thres, int core, int burn, int prerun, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repsum(summary,Summary(X, grade, sex));
  mat zerotheta = zeros(theta.n_cols);
  int percore = theta.n_rows/core;
  vec mc = zeros(theta.n_rows);
  mat alltheta(theta.n_rows*niter, theta.n_cols);
  mat allscore(theta.n_rows*niter, theta.n_cols);
  int numthe = theta.n_rows;
  vec nummc = zeros(niter);
  
  mat suf(imporsize*inner, theta.n_cols);
  mat thetas = map.t();
  
  for(int insize=0;  insize<imporsize; insize++){
    mat gib = Gibbs2(X,grade, sex, map,inner+burn);
    gib = gib.t();
    for(int innersize=0; innersize<inner; innersize++){
      suf.row(insize*inner+innersize)=gib.row(innersize+burn);
    }
  }
  
  
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums(theta.n_rows,theta.n_cols);
  mat suftheta = suf;
  
  for(int iter=0; iter<niter; iter++){
    mc = zeros(theta.n_rows);
    if(thetas.n_rows>prerun){
      for(int tt=0; tt<theta.n_rows; tt++){
        vec distval = crossdist2(theta.row(tt), thetas);
        vec minwhich = find_lowest_indices(distval, close);
        vec weights(close);
        weights.fill(1.0/close);
        mat selip2(imporsize*inner, theta.n_cols);
        mat norwei(imporsize*inner, theta.n_cols);
        for(int near=0; near<close; near++){
          mat selip(imporsize*inner, theta.n_cols);
          for(int jj=0; jj<imporsize*inner; jj++){
            selip.row(jj)=mmult((theta.row(tt)-thetas.row(minwhich[near])), suf.row((minwhich[near])*(imporsize*inner)+jj));
          }
          mat expselip=mexp(selip);
          vec selsum = colSum(expselip);
          for(int jj=0; jj<imporsize*inner; jj++){
            selip2.row(jj)=selip2.row(jj)+weights[near]*mmult((expselip.row(jj)/selsum.t()), suf.row((minwhich[near])*(imporsize*inner)+jj));
            mat nor = (expselip.row(jj)/selsum.t());
            norwei.row(jj)= norwei.row(jj)+ weights[near]*nor;
          }
        }
        vec selsum2 = colSum(selip2);
        norwei = mpow(norwei,2); 
        vec norweisum = colSum(norwei);
        norweisums.row(tt)=norweisum.t();
        parsum.row(tt)=selsum2.t();
        
        ess = norweisums;
        mat midess=minv(ess);
        midess = midess.t();
        
        
        if(min(midess.col(tt))<(thres)){
          mc[tt]=1;
        }
      }  
    }else{
      for(int tt=0; tt<theta.n_rows; tt++){
        int minwhich = vecminInd(crossdist2(theta.row(tt), thetas));
        mat selip(imporsize*inner, theta.n_cols);
        mat selip2(imporsize*inner, theta.n_cols);
        mat norwei(imporsize*inner, theta.n_cols);
        for(int jj=0; jj<imporsize*inner; jj++){
          selip.row(jj)=mmult((theta.row(tt)-thetas.row(minwhich)), suf.row((minwhich)*(imporsize*inner)+jj));
        }
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsize*inner; jj++){
          selip2.row(jj)=mmult((expselip.row(jj)/selsum.t()), suf.row((minwhich)*(imporsize*inner)+jj));
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
        
        
        if(min(midess.col(tt))<(thres)){
          mc[tt]=1;
        }
      }  
    }
    mat suf2(theta.n_rows*inner*imporsize, theta.n_cols);
    
    
    
    int cor;
#pragma omp parallel shared(theta) private(cor)
{
#pragma omp for schedule(static)
  for(cor=0; cor<core; cor++){
    for(int tt2=0; tt2<percore; tt2++){
      if(mc[percore*cor+tt2]==1){
        vec vectheta = vectorise(theta.row(percore*cor+tt2));
        parsum.row(percore*cor+tt2)=zerotheta.t();
        for(int insizes=0; insizes<imporsize; insizes++){
          mat gibs=Gibbs2(X,grade, sex, vectheta,inner+burn);
          gibs = gibs.t();
          for(int innersize=0; innersize<inner; innersize++){
            parsum.row(percore*cor+tt2)=parsum.row(percore*cor+tt2)+gibs.row(innersize+burn);
            suf2.row(((percore*cor+tt2)*inner*imporsize)+(insizes*inner+innersize))=gibs.row(innersize+burn);
          }
          parsum.row(percore*cor+tt2) = parsum.row(percore*cor+tt2)/inner;
        }
        parsum.row(percore*cor+tt2) = parsum.row(percore*cor+tt2)/imporsize; 
      }
    }
  }
}

for(int ttt=0; ttt<theta.n_rows; ttt++){    
  if(mc[ttt]==1){
    thetas = rbind_mat(thetas, theta.row(ttt));
    for(int iii=0; iii<inner*imporsize; iii++){
      suf = rbind_mat(suf, suf2.row(ttt*inner*imporsize+iii));
    }
  }
}
nummc[iter] = sum(mc);
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
for(int jjj=0; jjj<numthe; jjj++){
  alltheta.row(jjj+(iter*numthe)) = theta.row(jjj);
}

for(int jjjj=0; jjjj<numthe; jjjj++){
  mat lnp = lnpgrad.t();
  allscore.row(jjjj+(iter*numthe)) = lnp.row(jjjj);
}

  }
  return List::create(Named("theta")=theta, Named("alltheta")=alltheta, Named("allscore")=allscore, Named("thetas")=thetas, Named("suf")=suf, Named("nummc")=nummc);
}








// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List ergm_distweight_mcsvgd(mat x0, mat X, vec grade, vec sex, mat map, int niter, int inner, int imporsize, double close, double stepsize, double thres, int core, int burn, int prerun, double bandwidth=-1, double alpha=0.9){
  mat theta = x0;
  double fudgefactor = 1e-6;
  double historicalgrad = 0;
  mat summary(theta.n_cols, theta.n_rows);
  summary = repsum(summary,Summary(X, grade, sex));
  mat zerotheta = zeros(theta.n_cols);
  int percore = theta.n_rows/core;
  vec mc = zeros(theta.n_rows);
  mat alltheta(theta.n_rows*niter, theta.n_cols);
  mat allscore(theta.n_rows*niter, theta.n_cols);
  int numthe = theta.n_rows;
  vec nummc = zeros(niter);
  
  mat suf(imporsize*inner, theta.n_cols);
  mat thetas = map.t();
  
  for(int insize=0;  insize<imporsize; insize++){
    mat gib = Gibbs2(X,grade, sex, map,inner+burn);
    gib = gib.t();
    for(int innersize=0; innersize<inner; innersize++){
      suf.row(insize*inner+innersize)=gib.row(innersize+burn);
    }
  }
  
  
  mat parsum(theta.n_rows,theta.n_cols);
  mat ess(theta.n_rows,theta.n_cols);
  mat norweisums(theta.n_rows,theta.n_cols);
  mat suftheta = suf;
  
  for(int iter=0; iter<niter; iter++){
    mc = zeros(theta.n_rows);
    if(thetas.n_rows>prerun){
      for(int tt=0; tt<theta.n_rows; tt++){
        vec distval = crossdist2(theta.row(tt), thetas);
        vec minwhich = find_lowest_indices(distval, close);
        vec weights(close);
        for(int aa=0; aa<close; aa++){
          weights[aa]=distval[minwhich[aa]];
        }
        weights=weights/sum(weights);
        mat selip2(imporsize*inner, theta.n_cols);
        mat norwei(imporsize*inner, theta.n_cols);
        for(int near=0; near<close; near++){
          mat selip(imporsize*inner, theta.n_cols);
          for(int jj=0; jj<imporsize*inner; jj++){
            selip.row(jj)=mmult((theta.row(tt)-thetas.row(minwhich[near])), suf.row((minwhich[near])*(imporsize*inner)+jj));
          }
          mat expselip=mexp(selip);
          vec selsum = colSum(expselip);
          for(int jj=0; jj<imporsize*inner; jj++){
            selip2.row(jj)=selip2.row(jj)+weights[near]*mmult((expselip.row(jj)/selsum.t()), suf.row((minwhich[near])*(imporsize*inner)+jj));
            mat nor = (expselip.row(jj)/selsum.t());
            norwei.row(jj)= norwei.row(jj)+ weights[near]*nor;
          }
        }
        vec selsum2 = colSum(selip2);
        norwei = mpow(norwei,2); 
        vec norweisum = colSum(norwei);
        norweisums.row(tt)=norweisum.t();
        parsum.row(tt)=selsum2.t();
        
        ess = norweisums;
        mat midess=minv(ess);
        midess = midess.t();
        
        
        if(min(midess.col(tt))<(thres)){
          mc[tt]=1;
        }
      }  
    }else{
      for(int tt=0; tt<theta.n_rows; tt++){
        int minwhich = vecminInd(crossdist2(theta.row(tt), thetas));
        mat selip(imporsize*inner, theta.n_cols);
        mat selip2(imporsize*inner, theta.n_cols);
        mat norwei(imporsize*inner, theta.n_cols);
        for(int jj=0; jj<imporsize*inner; jj++){
          selip.row(jj)=mmult((theta.row(tt)-thetas.row(minwhich)), suf.row((minwhich)*(imporsize*inner)+jj));
        }
        mat expselip=mexp(selip);
        vec selsum = colSum(expselip);
        for(int jj=0; jj<imporsize*inner; jj++){
          selip2.row(jj)=mmult((expselip.row(jj)/selsum.t()), suf.row((minwhich)*(imporsize*inner)+jj));
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
        
        
        if(min(midess.col(tt))<(thres)){
          mc[tt]=1;
        }
      }  
    }
    mat suf2(theta.n_rows*inner*imporsize, theta.n_cols);
    
    
    
    int cor;
#pragma omp parallel shared(theta) private(cor)
{
#pragma omp for schedule(static)
  for(cor=0; cor<core; cor++){
    for(int tt2=0; tt2<percore; tt2++){
      if(mc[percore*cor+tt2]==1){
        vec vectheta = vectorise(theta.row(percore*cor+tt2));
        parsum.row(percore*cor+tt2)=zerotheta.t();
        for(int insizes=0; insizes<imporsize; insizes++){
          mat gibs=Gibbs2(X,grade, sex, vectheta,inner+burn);
          gibs = gibs.t();
          for(int innersize=0; innersize<inner; innersize++){
            parsum.row(percore*cor+tt2)=parsum.row(percore*cor+tt2)+gibs.row(innersize+burn);
            suf2.row(((percore*cor+tt2)*inner*imporsize)+(insizes*inner+innersize))=gibs.row(innersize+burn);
          }
          parsum.row(percore*cor+tt2) = parsum.row(percore*cor+tt2)/inner;
        }
        parsum.row(percore*cor+tt2) = parsum.row(percore*cor+tt2)/imporsize; 
      }
    }
  }
}

for(int ttt=0; ttt<theta.n_rows; ttt++){    
  if(mc[ttt]==1){
    thetas = rbind_mat(thetas, theta.row(ttt));
    for(int iii=0; iii<inner*imporsize; iii++){
      suf = rbind_mat(suf, suf2.row(ttt*inner*imporsize+iii));
    }
  }
}
nummc[iter] = sum(mc);
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
for(int jjj=0; jjj<numthe; jjj++){
  alltheta.row(jjj+(iter*numthe)) = theta.row(jjj);
}

for(int jjjj=0; jjjj<numthe; jjjj++){
  mat lnp = lnpgrad.t();
  allscore.row(jjjj+(iter*numthe)) = lnp.row(jjjj);
}

  }
  return List::create(Named("theta")=theta, Named("alltheta")=alltheta, Named("allscore")=allscore, Named("thetas")=thetas, Named("suf")=suf, Named("nummc")=nummc);
}






