rm(list=ls())
library(potts); library(Rcpp); library(RcppArmadillo)
library(sp)
library(gstat)
library(fields)
library(classInt)
library(spatstat)
library(MASS)
library(DiceKriging)
library(DiceDesign)
library(rdist)
library(doParallel)
library(doRNG)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
#======================================================================
# call data and functions
#======================================================================
load("ice_potts_data.RData")

sourceCpp("potts_RcppFtns.cpp")

set.seed(2024)
load("potts_neighbor.RData")
po<-glm(as.vector(x)~as.vector(same), family="poisson")
summary(po)



cycle = 30


COV = 1
sigma2 = 1

burn = 1000
outer = burn + 10000
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
adapIter = 1
outers = unique(c(0, seq(10000, outer, by = 10000), outer))

# initialize theta
beta = rgamma(1, 1, 1)

start = 1; postSamples = c(); Accprob = 0; rtime = 0


pottsDMHr2 <- function(foo, ncolor, stat, COV, beta, outer, cycle, updateCOV, sigma2, adaptInterval, adaptFactorExponent, adapIter){
  postSamples=rep(0, outer)
  accprob=rep(0, outer)
  c0=1
  c1=adaptFactorExponent
  ropt=0.234
  cholCOV=sqrt(sigma2*COV)
  
  for(l in 1:outer){
    if(updateCOV){
      if((l>=adaptInterval) && (l - (adaptInterval * trunc(((l)/adaptInterval))==0))){
        rhat = sum(accprob[(l-adaptInterval):(l-2)])/(adaptInterval-1)
        gamma1 = 1/(adapIter^c1)
        gamma2 = c0*gamma1
        sigma2 = exp(log(sigma2) + gamma2 * (rhat - ropt))
        
        COV = COV + gamma1 * ( var(postSamples[(l-adaptInterval):(l-2)]) - COV)
        cholCOV = sqrt(sigma2 * (COV + 0.00001))
        adapIter = adapIter + 1
      }
    }
    betaprev = beta
    betaprop = betaprev + cholCOV * rnorm(1)
    
    if(betaprop > 0){
      out <- potts(foo, c(rep(0, ncolor), betaprop), nbatch=cycle)
      xx <- unpackPotts(out$final)
      statprop <- calc_t(xx, ncolor)[ncolor+1]
      
      
      logprob = -0.05 * betaprop * betaprop + 0.05 * betaprev * betaprev +
        (statprop - stat) * (betaprev - betaprop)
      
      u = log(runif(1))
      if( u < logprob ){
        beta = betaprop
        accprob[l] = 1
      }
    }
    postSamples[l] = beta
  }
  list(postSamples=postSamples, accprob=accprob, adapIter=adapIter, sigma2=sigma2, COV=COV)
}



for(i in (start+1):length(outers) ){
  outeri = outers[i]-outers[i-1]
  
  ptm = proc.time()[3]
  dummy = pottsDMHr2(foo, ncolor, stat, COV, beta, outeri, cycle, updateCOV, 
                     sigma2, adaptInterval, adaptFactorExponent, adapIter)
  rtime = rtime + proc.time()[3] - ptm
  
  postSamples = c(postSamples, dummy$postSamples)
  Accprob = ( Accprob * outers[i-1] + sum(dummy$accprob) ) / outers[i]
  
  nSamples = length(dummy$postSamples)
  beta = dummy$postSamples[nSamples]
  COV = dummy$COV
  adapIter = dummy$adapIter
  sigma2 = dummy$sigma2
  
  save(postSamples, Accprob, rtime, beta, COV, adapIter, sigma2, file = "potts_ice_dmh_cycle30_r2_final.RData")
  
}




load("potts_ice_dmh_cycle30_r2_final.RData")




## RBF kernel in R
rbf_kernel <- function(theta, h=-1){
  sq_dist <- pdist(theta)
  pair_dist <- sq_dist^2
  if(h<0){
    h <- median(pair_dist)
    h <- sqrt(0.5*h / log(dim(theta)[1]+1))
  }
  
  # compute kernel
  
  kxy <- exp(-pair_dist / (h^2) /2)
  kxy <- as.matrix(kxy)
  theta <- as.matrix(theta)
  dxkxy <- -(kxy%*%theta)
  sumkxy <- colSums(kxy)
  for(i in 1:dim(theta)[2]){
    dxkxy[,i] <- dxkxy[,i] + theta[,i]*sumkxy
  }
  dxkxy <- dxkxy/(h^2)
  list(kxy=kxy, dxkxy=dxkxy)
}


potts_svi_r = function(foo, ncolor, x0, X, map, niter, inner, imporsize, imporsizes, stepsize, thres, bandwidth=-1, alpha=0.9){
  theta = x0
  fudgefactor = 1e-6
  historicalgrad = 0
  summary = matrix(rep(Summary2(X), nrow(theta)))
  suf = matrix(rep(0, imporsize*ncol(theta)), ncol=ncol(theta))
  thetas = map
  alltheta = matrix(rep(0, nrow(theta)*niter*ncol(theta)), ncol=ncol(theta))
  
  for(insize in 1:imporsize){
    suf[insize,] = potts_stat(foo, ncolor, map, inner)
  }
  
  parsum = matrix(rep(0, nrow(theta)*ncol(theta)), nrow=nrow(theta))
  ess = matrix(rep(0, nrow(theta)*ncol(theta)), nrow=nrow(theta))
  norweisums = rep(0, nrow(theta))
  for(iter in 1:niter){
    for(tt in 1:nrow(theta)){
      distmap = cdist(t(as.matrix(theta[tt,])), map)[1]
      distmin = min( cdist(t(as.matrix(theta[tt,])), thetas) )
      if(distmap > distmin) {
        minwhich = which(cdist(t(as.matrix(theta[tt,])), thetas)==min(cdist(t(as.matrix(theta[tt,])), thetas)))
        selip = rep(0, imporsizes)
        selip2= matrix(rep(0, imporsizes*ncol(theta)), nrow=imporsizes)
        norwei= rep(0, imporsizes)
        for(jj in 1:imporsizes){
          selip[jj]=(theta[tt,]-thetas[minwhich,])*suf[imporsize+(minwhich-2)*imporsizes+jj,]
        }
        mx = max(selip)
        selip=selip-mx
        expselip = exp(selip)
        selsum = sum(expselip)
        for(jj in 1:imporsizes){
          selip2[jj,]=(expselip[jj]/selsum)*suf[imporsize+(minwhich-2)*(imporsizes)+jj,]
          nor = (expselip[jj]/t(selsum))
          norwei[jj]= nor
        }
        selsum2 = colSums(selip2)
        norwei = norwei^2 
        norweisum = sum(norwei)
        norweisums[tt]=norweisum
        parsum[tt,]=t(selsum2)
        
        ess = norweisums
        midess=(ess^(-1))
        
        if( (min(midess[tt])) <(imporsizes/thres) ){
          thetas = rbind(thetas, theta[tt,])
          parsum[tt,]=0
          
          for(insizes in 1:imporsizes){
            suff = potts_stat(foo, ncolor, theta[tt,], inner)
            parsum[tt,]=parsum[tt,]+suff
            suf = rbind(suf, suff)
          }
          parsum[tt,] = parsum[tt,]/imporsizes
        }
      }else{
        selip = rep(0, imporsize)
        selip2= matrix(rep(0, imporsize*ncol(theta)), nrow=imporsize)
        norwei= rep(0, imporsize)
        for(jj in 1:imporsize){
          selip[jj] = (theta[tt,]-map)*suf[jj,]
        }
        mx = max(selip)
        selip = selip-mx
        expselip = exp(selip)
        selsum = sum(expselip)
        for(jj in 1:imporsize){
          selip2[jj,]=(expselip[jj]/selsum)*suf[jj,]
          nor = (expselip[jj]/selsum)
          norwei[jj]= nor
        }
        
        selsum2 = colSums(selip2)
        norwei = norwei^2
        norweisum = sum(norwei)
        norweisums[tt]=norweisum
        parsum[tt,]=t(selsum2)
        
        ess = norweisums
        midess=(ess^(-1))
        
        if( (min(midess[tt]))<(imporsize/thres)){
          thetas = rbind(thetas, theta[tt,])
          parsum[tt,]=0
          for(insizes in 1:imporsizes){
            suff = potts_stat(foo, ncolor, theta[tt,], inner)
            parsum[tt,]=parsum[tt,]+suff
            suf = rbind(suf, suff)
          }
          parsum[tt,] = parsum[tt,]/imporsizes
        }
      }
    }  
    lnpgrad = t(summary) - t(parsum) 
    ker=rbf_kernel(theta, -1)
    kxy = ker$kxy
    dxkxy = ker$dxkxy
    gradtheta = kxy%*%t(lnpgrad)
    gradtheta = gradtheta + dxkxy
    gradtheta = gradtheta / nrow(theta)
    
    pgradtheta = gradtheta^2
    
    histo = matrix(rep(0, nrow(pgradtheta)*ncol(pgradtheta)), nrow=nrow(pgradtheta))
    fudge = matrix(rep(0, nrow(pgradtheta)*ncol(pgradtheta)), nrow=nrow(pgradtheta))
    fudge= fudge+ fudgefactor
    if(iter==1){
      histo = histo + pgradtheta
    } else {
      histo = alpha * histo + (1-alpha) * pgradtheta
    }
    adjgrad = gradtheta/(fudge+sqrt(histo))
    theta = theta + stepsize*adjgrad
    
    for(jjj in 1:nrow(theta)){
      alltheta[jjj+((iter-1)*nrow(theta)),] = theta[jjj,]
    }
  }
  list(theta=theta, thetas=thetas, alltheta = alltheta)
}

map=mean(postSamples)
M=3
x0<-rnorm(M,mean=beta, sd=0.01)
x0<-as.matrix(x0)
X<-as.matrix(x)
ptm <-proc.time()
potts_svi_result<-potts_svi_r(foo, ncolor, x0, X, map, 300, 30, 10, 3, 0.0001, 0)
MCMCtime=proc.time()-ptm
MCMCtime





potts_svi_r_par = function(foo, ncolor, x0, X, map, niter, inner, imporsize, imporsizes, stepsize, thres, bandwidth=-1, alpha=0.9){
  theta = x0
  fudgefactor = 1e-6
  historicalgrad = 0
  summary = matrix(rep(Summary2(X), nrow(theta)))
  suf = matrix(rep(0, imporsize*ncol(theta)), ncol=ncol(theta))
  thetas = map
  alltheta = matrix(rep(0, nrow(theta)*niter*ncol(theta)), ncol=ncol(theta))
  mc = rep(0, nrow(theta))
  nummc = rep(0, niter)
  
  #for(insize in 1:imporsize){
  #  suf[insize,] = potts_stat(foo, ncolor, map, inner)
  #}
  out<-potts(foo, c(rep(0,ncolor), map), nbatch=inner+imporsize)
  suff<-out$batch[,ncolor+1]
  for(insize in 1:imporsize){
    suf[insize,] = suff[inner+insize]
  }
  
  parsum = matrix(rep(0, nrow(theta)*ncol(theta)), nrow=nrow(theta))
  ess = matrix(rep(0, nrow(theta)*ncol(theta)), nrow=nrow(theta))
  norweisums = rep(0, nrow(theta))
  for(iter in 1:niter){
    mc = rep(0, nrow(theta))
    for(tt in 1:nrow(theta)){
      distmap = cdist(t(as.matrix(theta[tt,])), map)[1]
      distmin = min( cdist(t(as.matrix(theta[tt,])), thetas) )
      if(distmap > distmin) {
        minwhich = which(cdist(t(as.matrix(theta[tt,])), thetas)==min(cdist(t(as.matrix(theta[tt,])), thetas)))
        selip = rep(0, imporsizes)
        selip2= matrix(rep(0, imporsizes*ncol(theta)), nrow=imporsizes)
        norwei= rep(0, imporsizes)
        for(jj in 1:imporsizes){
          selip[jj]=(theta[tt,]-thetas[minwhich,])*suf[imporsize+(minwhich-2)*imporsizes+jj,]
        }
        mx = max(selip)
        selip=selip-mx
        expselip = exp(selip)
        selsum = sum(expselip)
        for(jj in 1:imporsizes){
          selip2[jj,]=(expselip[jj]/selsum)*suf[imporsize+(minwhich-2)*(imporsizes)+jj,]
          nor = (expselip[jj]/t(selsum))
          norwei[jj]= nor
        }
        selsum2 = colSums(selip2)
        norwei = norwei^2 
        norweisum = sum(norwei)
        norweisums[tt]=norweisum
        parsum[tt,]=t(selsum2)
        
        ess = norweisums
        midess=(ess^(-1))
        
        if( (min(midess[tt])) <(imporsizes/thres) ){
          mc[tt]=1
        }
      }else{
        selip = rep(0, imporsize)
        selip2= matrix(rep(0, imporsize*ncol(theta)), nrow=imporsize)
        norwei= rep(0, imporsize)
        for(jj in 1:imporsize){
          selip[jj] = (theta[tt,]-map)*suf[jj,]
        }
        mx = max(selip)
        selip = selip-mx
        expselip = exp(selip)
        selsum = sum(expselip)
        for(jj in 1:imporsize){
          selip2[jj,]=(expselip[jj]/selsum)*suf[jj,]
          nor = (expselip[jj]/selsum)
          norwei[jj]= nor
        }
        
        selsum2 = colSums(selip2)
        norwei = norwei^2
        norweisum = sum(norwei)
        norweisums[tt]=norweisum
        parsum[tt,]=t(selsum2)
        
        ess = norweisums
        midess=(ess^(-1))
        
        
        if( (min(midess[tt])) <(imporsize/thres) ){
          mc[tt]=1
        }
        
      }
    }  
    
    
    # suf2 = matrix(rep(0, nrow(theta)*imporsizes*ncol(theta)), ncol=ncol(theta))
    
    
    rr<-foreach(i=1:nrow(theta)) %dopar% {
      library(potts)
      if(mc[i]==1){
        parsum[i,]=0
        out<-potts(foo, c(rep(0,ncolor), theta[i,]), nbatch=inner+imporsizes)
        suff<-out$batch[,ncolor+1]
        suff
      }else{parsum[i,]}
    }
    for(i in 1:nrow(theta)){
      if(mc[i]==1){
        thetas = rbind(thetas, theta[i,])
        
        suf=rbind(suf,as.matrix(rr[[i]][(inner+1):(inner+imporsizes)]))
        parsum[i,]<-mean(rr[[i]][(inner+1):(inner+imporsizes)])
      }
    }
    nummc[iter]=sum(mc)
    lnpgrad = t(summary) - t(parsum) 
    ker=rbf_kernel(theta, -1)
    kxy = ker$kxy
    dxkxy = ker$dxkxy
    gradtheta = kxy%*%t(lnpgrad)
    gradtheta = gradtheta + dxkxy
    gradtheta = gradtheta / nrow(theta)
    
    pgradtheta = gradtheta^2
    
    histo = matrix(rep(0, nrow(pgradtheta)*ncol(pgradtheta)), nrow=nrow(pgradtheta))
    fudge = matrix(rep(0, nrow(pgradtheta)*ncol(pgradtheta)), nrow=nrow(pgradtheta))
    fudge= fudge+ fudgefactor
    if(iter==1){
      histo = histo + pgradtheta
    } else {
      histo = alpha * histo + (1-alpha) * pgradtheta
    }
    adjgrad = gradtheta/(fudge+sqrt(histo))
    theta = theta + stepsize*adjgrad
    
    for(jjj in 1:nrow(theta)){
      alltheta[jjj+((iter-1)*nrow(theta)),] = theta[jjj,]
    }
  }
  list(theta=theta, thetas=thetas, alltheta = alltheta, nummc=nummc)
}


registerDoParallel(32)
getDoParWorkers()
map=mean(potts_svi_result$theta)
M=64
x0<-rnorm(M,mean=map, sd=0.01)
x0<-as.matrix(x0)
X<-as.matrix(x)
ptm <-proc.time()
potts_svi_result<-potts_svi_r_par(foo, ncolor, x0, X, map, 500, 30, 50, 50, 0.0001, 3)
MCMCtime=proc.time()-ptm
MCMCtime


library(HDInterval)
print(hdi(potts_svi_result$theta))

plot(density(postSamples[1001:11000]), col="red")
lines(density(potts_svi_result$theta), col="black")
save(potts_svi_result, MCMCtime, file="MCSVGD_POTTS_R_par64_iter500_inner30_impor50_50_step0.0001_thres3_new.RData")


