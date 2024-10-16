rm(list=ls())
library(COMPoissonReg)
library(Rcpp);library(RcppArmadillo);library(MASS)
library(rdist);library(fields);library(lattice)
library(Bergm);library(coda);library(xtable)
library(DiceKriging);library(DiceDesign)
library(dplyr);library(sp);library(gstat)
library(classInt);library(spatstat)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

sourceCpp('MCSVGD_COMP.cpp')


set.seed(2024)
###########################   n=10   ############################################################

load("simData2dim_n10_inter.RData")
betaInd = 1:p
beta_initial = rnorm(p)
alpha_initial = log(1) # log(nu) = alpha

b = c(rep(1, p), 2)
block = sapply(1:length(unique(b)), function(i) which(b == i) - 1, simplify = F)
k = length(block)

sigma2 = rep(1^2, k)
COV = sapply(1:length(block), function(i) diag( length(block[[i]]) ), simplify = F)

### run real
thin = 1
burn = 1000
niters = 50000 * thin + burn
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
adapIter = 1

### run exchange algorithm
ptm = proc.time()[3]
Murray = exactCOMPreg_fixnu(niters, y, X, beta_initial, alpha_initial, block, sigma2, COV, updateCOV, 
                            adaptInterval, adaptFactorExponent, adapIter, thin)
rtime = proc.time()[3] - ptm

postSamp = Murray$Sample
accprob = colMeans(Murray$Accprob)

save(Murray, rtime, postSamp, accprob, file = paste0('simuexchange_50000_dim2_n10_inter_final2.RData'))

#load exchange algorithm result
load("simuexchange_50000_dim2_n10_inter_final2.RData")

M=1
x0<-rmvnorm(M,mean=c(beta_initial, alpha_initial), sigma=diag(4))
x0<-as.matrix(x0)
pt <- proc.time()
stein_comp_svi<-comp_fixnu_map(x0, X, y, 300, 10, 0.005, exp(0.5))
ptFinal_glm <- proc.time()-pt
ptFinal_glm
map<-as.matrix(colMeans(stein_comp_svi[291:300,]))


colnames(X)<-c("x1","x2","x3")
datas<-cbind(y,X)
datas<-data.frame(datas)
glm_re<-glm.cmp(y~x1+x2, data=datas)
cov1<-vcov(glm_re)[2:3, 2:3]
cov2<-rbind(cov1,vcov(glm_re)[1, 2:3])
cov3<-cbind(cov2,vcov(glm_re)[c(2,3,1),1])
cov4<-rbind(cov3,vcov(glm_re)[4,c(2,3,1)])
cov5<-cbind(cov4, vcov(glm_re)[c(2,3,1,4),4])

M=96
x0<-rmvnorm(M,mean=map, sigma=cov5)
x0<-as.matrix(x0)

pt <- proc.time()
stein_comp_svi<-comp_svi_fixnu_par2(x0, X, y, map, 500, 50, 50, 0.0005, 3, exp(0.5),32)
ptFinal_glm <- proc.time()-pt
ptFinal_glm

pos_in<-matrix(rep(0,9), nrow=3)
pos_in[,1]<-colMeans(stein_comp_svi$theta)[1:3]
library(HDInterval)
for(i in 1:3){
  pos_in[i,2:3]<-hdi(stein_comp_svi$theta[,i])
}
save(stein_comp_svi,ptFinal_glm, pos_in,file="comp_svi_dim2_n10_inter_par96_iter500_impor50_50_thres3_step0.0005_cov_glmcmp_final_new.RData")

par(mfrow=c(1,3))
for(jj in 1:3){
  plot(density(postSamp[1001:50000,jj]), col="red")
  lines(density(stein_comp_svi$theta[,jj]),col="black")
  lines(density(x0[,jj]), col="blue")
}


rm(list=ls())
library(COMPoissonReg)
library(Rcpp);library(RcppArmadillo);library(MASS)
library(rdist);library(fields);library(lattice)
library(Bergm);library(coda);library(xtable)
library(DiceKriging);library(DiceDesign)
library(dplyr);library(sp);library(gstat)
library(classInt);library(spatstat)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

sourceCpp('MCSVGD_COMP.cpp')

set.seed(2024)
###################################  n=15 ####################################################

load("simData2dim_n15_inter.RData")
betaInd = 1:p
beta_initial = rnorm(p)
alpha_initial = log(1) # log(nu) = alpha

b = c(rep(1, p), 2)
block = sapply(1:length(unique(b)), function(i) which(b == i) - 1, simplify = F)
k = length(block)

sigma2 = rep(1^2, k)
COV = sapply(1:length(block), function(i) diag( length(block[[i]]) ), simplify = F)

### run real
thin = 1
burn = 1000
niters = 50000 * thin + burn
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
adapIter = 1


ptm = proc.time()[3]
Murray = exactCOMPreg_fixnu(niters, y, X, beta_initial, alpha_initial, block, sigma2, COV, updateCOV, 
                            adaptInterval, adaptFactorExponent, adapIter, thin)
rtime = proc.time()[3] - ptm

postSamp = Murray$Sample
accprob = colMeans(Murray$Accprob)

save(Murray, rtime, postSamp, accprob, file = paste0('simuexchange_50000_dim2_n15_inter_final2.RData'))

load("simuexchange_50000_dim2_n15_inter_final2.RData")

colnames(X)<-c("x1","x2","x3")
datas<-cbind(y,X)
datas<-data.frame(datas)
glm_re<-glm.cmp(y~x1+x2, data=datas)
cov1<-vcov(glm_re)[2:3, 2:3]
cov2<-rbind(cov1,vcov(glm_re)[1, 2:3])
cov3<-cbind(cov2,vcov(glm_re)[c(2,3,1),1])
cov4<-rbind(cov3,vcov(glm_re)[4,c(2,3,1)])
cov5<-cbind(cov4, vcov(glm_re)[c(2,3,1,4),4])

M=1
x0<-rmvnorm(M,mean=c(beta_initial, alpha_initial), sigma=diag(4))
x0<-as.matrix(x0)
pt <- proc.time()
stein_comp_svi<-comp_fixnu_map(x0, X, y, 300, 10, 0.005, exp(0.5))
ptFinal_glm <- proc.time()-pt
ptFinal_glm
map<-as.matrix(colMeans(stein_comp_svi[291:300,]))


M=96
x0<-rmvnorm(M,mean=map, sigma=cov5)
x0<-as.matrix(x0)

pt <- proc.time()
stein_comp_svi<-comp_svi_fixnu_par2(x0, X, y, map, 500, 50, 50, 0.0005, 3, exp(0.5),32)
ptFinal_glm <- proc.time()-pt
ptFinal_glm

par(mfrow=c(1,3))
for(jj in 1:3){
  plot(density(postSamp[1001:50000,jj]), col="red")
  lines(density(stein_comp_svi$theta[,jj]),col="black")
  #lines(density(x0[,jj]), col="blue")
}

pos_in<-matrix(rep(0,9), nrow=3)
pos_in[,1]<-colMeans(stein_comp_svi$theta)[1:3]
library(HDInterval)
for(i in 1:3){
  pos_in[i,2:3]<-hdi(stein_comp_svi$theta[,i])
}
save(stein_comp_svi,ptFinal_glm, pos_in, file="comp_svi_dim2_n15_inter_par96_iter500_impor50_thres3_step0.0005_hoho.RData")

par(mfrow=c(1,3))
for(jj in 1:3){
  plot(density(postSamp[1001:50000,jj]), col="red")
  lines(density(stein_comp_svi$theta[,jj]),col="black")
  #lines(density(x0[,jj]), col="blue")
}


rm(list=ls())
library(Rcpp);library(RcppArmadillo);library(MASS)
library(rdist);library(fields);library(lattice)
library(Bergm);library(coda);library(xtable)
library(DiceKriging);library(DiceDesign)
library(dplyr);library(sp);library(gstat)
library(classInt);library(spatstat)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

sourceCpp('MCSVGD_COMP.cpp')





set.seed(2024)
###########################   n=20   ############################################################



load("simData2dim_n20_inter.RData")
betaInd = 1:p
beta_initial = rnorm(p)
alpha_initial = log(1) # log(nu) = alpha

b = c(rep(1, p), 2)
block = sapply(1:length(unique(b)), function(i) which(b == i) - 1, simplify = F)
k = length(block)

sigma2 = rep(1^2, k)
COV = sapply(1:length(block), function(i) diag( length(block[[i]]) ), simplify = F)


### run real
thin = 1
burn = 1000
niters = 50000 * thin + burn
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
adapIter = 1


ptm = proc.time()[3]
Murray = exactCOMPreg_fixnu(niters, y, X, beta_initial, alpha_initial, block, sigma2, COV, updateCOV, 
                            adaptInterval, adaptFactorExponent, adapIter, thin)
rtime = proc.time()[3] - ptm

postSamp = Murray$Sample
accprob = colMeans(Murray$Accprob)

save(Murray, rtime, postSamp, accprob, file = paste0('simuexchange_50000_dim2_n20_inter_final2.RData'))

load("simuexchange_50000_dim2_n20_inter_final2.RData")


M=1
x0<-rmvnorm(M,mean=c(beta_initial, alpha_initial), sigma=diag(4))
x0<-as.matrix(x0)
pt <- proc.time()
stein_comp_svi<-comp_fixnu_map(x0, X, y, 300, 10, 0.005, exp(0.5))
ptFinal_glm <- proc.time()-pt
ptFinal_glm
map<-as.matrix(colMeans(stein_comp_svi[291:300,]))

colnames(X)<-c("x1","x2","x3")
datas<-cbind(y,X)
datas<-data.frame(datas)
glm_re<-glm.cmp(y~x1+x2, data=datas)
cov1<-vcov(glm_re)[2:3, 2:3]
cov2<-rbind(cov1,vcov(glm_re)[1, 2:3])
cov3<-cbind(cov2,vcov(glm_re)[c(2,3,1),1])
cov4<-rbind(cov3,vcov(glm_re)[4,c(2,3,1)])
cov5<-cbind(cov4, vcov(glm_re)[c(2,3,1,4),4])

M=96
x0<-rmvnorm(M,mean=map, sigma=cov5)
x0<-as.matrix(x0)


pt <- proc.time()
stein_comp_svi<-comp_svi_fixnu_par2(x0, X, y, map, 500, 50, 50, 0.0005, 3, exp(0.5),32)
ptFinal_glm <- proc.time()-pt
ptFinal_glm

pos_in<-matrix(rep(0,9), nrow=3)
pos_in[,1]<-colMeans(stein_comp_svi$theta)[1:3]
library(HDInterval)
for(i in 1:3){
  pos_in[i,2:3]<-hdi(stein_comp_svi$theta[,i])
}

save(stein_comp_svi,ptFinal_glm,  pos_in, file="comp_svi_dim2_n20_inter_par96_iter500_impor50_50_thres3_step0.0005_cov_glmcmp_final_new.RData")
par(mfrow=c(1,3))
for(jj in 1:3){
  plot(density(postSamp[1001:50000,jj]), col="red")
  lines(density(stein_comp_svi$theta[,jj]),col="black")
  lines(density(x0[,jj]), col="blue")
}





rm(list=ls())
library(Rcpp);library(RcppArmadillo);library(MASS)
library(rdist);library(fields);library(lattice)
library(Bergm);library(coda);library(xtable)
library(DiceKriging);library(DiceDesign)
library(dplyr);library(sp);library(gstat)
library(classInt);library(spatstat)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

sourceCpp('MCSVGD_COMP.cpp')




set.seed(2024)
###########################   n=25   ############################################################



load("simData2dim_n25_inter.RData")
betaInd = 1:p
beta_initial = rnorm(p)
alpha_initial = log(1) # log(nu) = alpha

b = c(rep(1, p), 2)
block = sapply(1:length(unique(b)), function(i) which(b == i) - 1, simplify = F)
k = length(block)

sigma2 = rep(1^2, k)
COV = sapply(1:length(block), function(i) diag( length(block[[i]]) ), simplify = F)
### run real
thin = 1
burn = 1000
niters = 50000 * thin + burn
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
adapIter = 1


ptm = proc.time()[3]
Murray = exactCOMPreg_fixnu(niters, y, X, beta_initial, alpha_initial, block, sigma2, COV, updateCOV, 
                            adaptInterval, adaptFactorExponent, adapIter, thin)
rtime = proc.time()[3] - ptm

postSamp = Murray$Sample
accprob = colMeans(Murray$Accprob)

save(Murray, rtime, postSamp, accprob, file = paste0('simuexchange_50000_dim2_n25_inter_final2.RData'))

load("simuexchange_50000_dim2_n25_inter_final2.RData")

M=1
x0<-rmvnorm(M,mean=c(beta_initial, alpha_initial), sigma=diag(4))
x0<-as.matrix(x0)
pt <- proc.time()
stein_comp_svi<-comp_fixnu_map(x0, X, y, 300, 10, 0.005, exp(0.5))
ptFinal_glm <- proc.time()-pt
ptFinal_glm
map<-as.matrix(colMeans(stein_comp_svi[291:300,]))

colnames(X)<-c("x1","x2","x3")
datas<-cbind(y,X)
datas<-data.frame(datas)
glm_re<-glm.cmp(y~x1+x2, data=datas)
cov1<-vcov(glm_re)[2:3, 2:3]
cov2<-rbind(cov1,vcov(glm_re)[1, 2:3])
cov3<-cbind(cov2,vcov(glm_re)[c(2,3,1),1])
cov4<-rbind(cov3,vcov(glm_re)[4,c(2,3,1)])
cov5<-cbind(cov4, vcov(glm_re)[c(2,3,1,4),4])

M=96
x0<-rmvnorm(M,mean=map, sigma=cov5)
x0<-as.matrix(x0)


pt <- proc.time()
stein_comp_svi<-comp_svi_fixnu_par2(x0, X, y, map, 500, 50, 50, 0.0001, 3, exp(0.5),32)
ptFinal_glm <- proc.time()-pt
ptFinal_glm

pos_in<-matrix(rep(0,9), nrow=3)
pos_in[,1]<-colMeans(stein_comp_svi$theta)[1:3]
library(HDInterval)
for(i in 1:3){
  pos_in[i,2:3]<-hdi(stein_comp_svi$theta[,i])
}
save(stein_comp_svi,ptFinal_glm, pos_in,  file="comp_svi_dim2_n25_inter_par96_iter500_impor50_50_thres3_step0.0001_cov_glmcmp_final_new.RData")
par(mfrow=c(1,3))
for(jj in 1:3){
  plot(density(postSamp[1001:50000,jj]), col="red")
  lines(density(stein_comp_svi$theta[,jj]),col="black")
  lines(density(x0[,jj]), col="blue")
}




rm(list=ls())
library(Rcpp);library(RcppArmadillo);library(MASS)
library(rdist);library(fields);library(lattice)
library(Bergm);library(coda);library(xtable)
library(DiceKriging);library(DiceDesign)
library(dplyr);library(sp);library(gstat)
library(classInt);library(spatstat)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

sourceCpp('MCSVGD_COMP.cpp')





set.seed(2024)
###########################   n=30   ############################################################



load("simData2dim_n30_inter.RData")
betaInd = 1:p
beta_initial = rnorm(p)
alpha_initial = log(1) # log(nu) = alpha

b = c(rep(1, p), 2)
block = sapply(1:length(unique(b)), function(i) which(b == i) - 1, simplify = F)
k = length(block)

sigma2 = rep(1^2, k)
COV = sapply(1:length(block), function(i) diag( length(block[[i]]) ), simplify = F)
### run real
thin = 1
burn = 1000
niters = 50000 * thin + burn
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
adapIter = 1


ptm = proc.time()[3]
Murray = exactCOMPreg_fixnu(niters, y, X, beta_initial, alpha_initial, block, sigma2, COV, updateCOV, 
                            adaptInterval, adaptFactorExponent, adapIter, thin)
rtime = proc.time()[3] - ptm

postSamp = Murray$Sample
accprob = colMeans(Murray$Accprob)

save(Murray, rtime, postSamp, accprob, file = paste0('simuexchange_50000_dim2_n30_inter_final2.RData'))

load("simuexchange_50000_dim2_n30_inter_final2.RData")



M=1
x0<-rmvnorm(M,mean=c(beta_initial, alpha_initial), sigma=diag(4))
x0<-as.matrix(x0)
pt <- proc.time()
stein_comp_svi<-comp_fixnu_map(x0, X, y, 300, 10, 0.005, exp(0.5))
ptFinal_glm <- proc.time()-pt
ptFinal_glm
map<-as.matrix(colMeans(stein_comp_svi[291:300,]))

colnames(X)<-c("x1","x2","x3")
datas<-cbind(y,X)
datas<-data.frame(datas)
glm_re<-glm.cmp(y~x1+x2, data=datas)
cov1<-vcov(glm_re)[2:3, 2:3]
cov2<-rbind(cov1,vcov(glm_re)[1, 2:3])
cov3<-cbind(cov2,vcov(glm_re)[c(2,3,1),1])
cov4<-rbind(cov3,vcov(glm_re)[4,c(2,3,1)])
cov5<-cbind(cov4, vcov(glm_re)[c(2,3,1,4),4])

M=96
x0<-rmvnorm(M,mean=map, sigma=cov5)
x0<-as.matrix(x0)


pt <- proc.time()
stein_comp_svi<-comp_svi_fixnu_par2(x0, X, y, map, 500, 50, 50, 0.0001, 3, exp(0.5),32)
ptFinal_glm <- proc.time()-pt
ptFinal_glm

pos_in<-matrix(rep(0,9), nrow=3)
pos_in[,1]<-colMeans(stein_comp_svi$theta)[1:3]
library(HDInterval)
for(i in 1:3){
  pos_in[i,2:3]<-hdi(stein_comp_svi$theta[,i])
}
save(stein_comp_svi,ptFinal_glm,  pos_in, file="comp_svi_dim2_n30_inter_par96_iter500_impor50_50_thres3_step0.0001_cov_glmcmp_final_new.RData")
par(mfrow=c(1,3))
for(jj in 1:3){
  plot(density(postSamp[1001:50000,jj]), col="red")
  lines(density(stein_comp_svi$theta[,jj]),col="black")
  lines(density(x0[,jj]), col="blue")
}


rm(list=ls())
library(Rcpp);library(RcppArmadillo);library(MASS)
library(rdist);library(fields);library(lattice)
library(Bergm);library(coda);library(xtable)
library(DiceKriging);library(DiceDesign)
library(dplyr);library(sp);library(gstat)
library(classInt);library(spatstat)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

sourceCpp('MCSVGD_COMP.cpp')





set.seed(2024)
###########################   n=35   ############################################################



load("simData2dim_n35_inter.RData")
betaInd = 1:p
beta_initial = rnorm(p)
alpha_initial = log(1) # log(nu) = alpha

b = c(rep(1, p), 2)
block = sapply(1:length(unique(b)), function(i) which(b == i) - 1, simplify = F)
k = length(block)

sigma2 = rep(1^2, k)
COV = sapply(1:length(block), function(i) diag( length(block[[i]]) ), simplify = F)
### run real
thin = 1
burn = 1000
niters = 50000 * thin + burn
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
adapIter = 1


ptm = proc.time()[3]
Murray = exactCOMPreg_fixnu(niters, y, X, beta_initial, alpha_initial, block, sigma2, COV, updateCOV, 
                            adaptInterval, adaptFactorExponent, adapIter, thin)
rtime = proc.time()[3] - ptm

postSamp = Murray$Sample
accprob = colMeans(Murray$Accprob)

save(Murray, rtime, postSamp, accprob, file = paste0('simuexchange_50000_dim2_n35_inter_final2.RData'))

load("simuexchange_50000_dim2_n35_inter_final2.RData")


M=1
x0<-rmvnorm(M,mean=c(beta_initial, alpha_initial), sigma=diag(4))
x0<-as.matrix(x0)
pt <- proc.time()
stein_comp_svi<-comp_fixnu_map(x0, X, y, 300, 10, 0.005, exp(0.5))
ptFinal_glm <- proc.time()-pt
ptFinal_glm
map<-as.matrix(colMeans(stein_comp_svi[291:300,]))

colnames(X)<-c("x1","x2","x3")
datas<-cbind(y,X)
datas<-data.frame(datas)
glm_re<-glm.cmp(y~x1+x2, data=datas)
cov1<-vcov(glm_re)[2:3, 2:3]
cov2<-rbind(cov1,vcov(glm_re)[1, 2:3])
cov3<-cbind(cov2,vcov(glm_re)[c(2,3,1),1])
cov4<-rbind(cov3,vcov(glm_re)[4,c(2,3,1)])
cov5<-cbind(cov4, vcov(glm_re)[c(2,3,1,4),4])

M=96
x0<-rmvnorm(M,mean=map, sigma=cov5)
x0<-as.matrix(x0)


pt <- proc.time()
stein_comp_svi<-comp_svi_fixnu_par2(x0, X, y, map, 500, 50, 50, 0.0001, 3, exp(0.5),32)
ptFinal_glm <- proc.time()-pt
ptFinal_glm

pos_in<-matrix(rep(0,9), nrow=3)
pos_in[,1]<-colMeans(stein_comp_svi$theta)[1:3]
library(HDInterval)
for(i in 1:3){
  pos_in[i,2:3]<-hdi(stein_comp_svi$theta[,i])
}
save(stein_comp_svi,ptFinal_glm,  pos_in, file="comp_svi_dim2_n35_inter_par96_iter500_impor50_50_thres3_step0.0001_cov_glmcmp_final_new.RData")
par(mfrow=c(1,3))
for(jj in 1:3){
  plot(density(postSamp[1001:50000,jj]), col="red")
  lines(density(stein_comp_svi$theta[,jj]),col="black")
  lines(density(x0[,jj]), col="blue")
}




rm(list=ls())
library(Rcpp);library(RcppArmadillo);library(MASS)
library(rdist);library(fields);library(lattice)
library(Bergm);library(coda);library(xtable)
library(DiceKriging);library(DiceDesign)
library(dplyr);library(sp);library(gstat)
library(classInt);library(spatstat)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

sourceCpp('MCSVGD_COMP.cpp')





set.seed(2032)
###########################   n=40   ############################################################



load("simData2dim_n40_inter.RData")
betaInd = 1:p
beta_initial = rnorm(p)
alpha_initial = log(1) # log(nu) = alpha

b = c(rep(1, p), 2)
block = sapply(1:length(unique(b)), function(i) which(b == i) - 1, simplify = F)
k = length(block)

sigma2 = rep(1^2, k)
COV = sapply(1:length(block), function(i) diag( length(block[[i]]) ), simplify = F)
### run real
thin = 1
burn = 1000
niters = 50000 * thin + burn
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
adapIter = 1


ptm = proc.time()[3]
Murray = exactCOMPreg_fixnu(niters, y, X, beta_initial, alpha_initial, block, sigma2, COV, updateCOV, 
                            adaptInterval, adaptFactorExponent, adapIter, thin)
rtime = proc.time()[3] - ptm

postSamp = Murray$Sample
accprob = colMeans(Murray$Accprob)

save(Murray, rtime, postSamp, accprob, file = paste0('simuexchange_50000_dim2_n40_inter_final2.RData'))

load("simuexchange_50000_dim2_n40_inter_final2.RData")

M=1
x0<-rmvnorm(M,mean=c(beta_initial, alpha_initial), sigma=diag(4))
x0<-as.matrix(x0)
pt <- proc.time()
stein_comp_svi<-comp_fixnu_map(x0, X, y, 300, 10, 0.005, exp(0.5))
ptFinal_glm <- proc.time()-pt
ptFinal_glm
map<-as.matrix(colMeans(stein_comp_svi[291:300,]))

colnames(X)<-c("x1","x2","x3")
datas<-cbind(y,X)
datas<-data.frame(datas)
glm_re<-glm.cmp(y~x1+x2, data=datas)
cov1<-vcov(glm_re)[2:3, 2:3]
cov2<-rbind(cov1,vcov(glm_re)[1, 2:3])
cov3<-cbind(cov2,vcov(glm_re)[c(2,3,1),1])
cov4<-rbind(cov3,vcov(glm_re)[4,c(2,3,1)])
cov5<-cbind(cov4, vcov(glm_re)[c(2,3,1,4),4])

M=96
x0<-rmvnorm(M,mean=map, sigma=cov5)
x0<-as.matrix(x0)


pt <- proc.time()
stein_comp_svi<-comp_svi_fixnu_par2(x0, X, y, map, 500, 50, 50, 0.0001, 3, exp(0.5),32)
ptFinal_glm <- proc.time()-pt
ptFinal_glm

pos_in<-matrix(rep(0,9), nrow=3)
pos_in[,1]<-colMeans(stein_comp_svi$theta)[1:3]
library(HDInterval)
for(i in 1:3){
  pos_in[i,2:3]<-hdi(stein_comp_svi$theta[,i])
}
save(stein_comp_svi,ptFinal_glm,  pos_in, file="comp_svi_dim2_n40_inter_par96_iter500_impor50_50_thres3_step0.0001_cov_glmcmp_final_new.RData")
par(mfrow=c(1,3))
for(jj in 1:3){
  plot(density(postSamp[1001:50000,jj]), col="red")
  lines(density(stein_comp_svi$theta[,jj]),col="black")
  lines(density(x0[,jj]), col="blue")
}





rm(list=ls())
library(Rcpp);library(RcppArmadillo);library(MASS)
library(rdist);library(fields);library(lattice)
library(Bergm);library(coda);library(xtable)
library(DiceKriging);library(DiceDesign)
library(dplyr);library(sp);library(gstat)
library(classInt);library(spatstat)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

sourceCpp('MCSVGD_COMP.cpp')





set.seed(2100)
###########################   n=45   ############################################################



load("simData2dim_n45_inter.RData")
betaInd = 1:p
beta_initial = rnorm(p)
alpha_initial = log(1) # log(nu) = alpha

b = c(rep(1, p), 2)
block = sapply(1:length(unique(b)), function(i) which(b == i) - 1, simplify = F)
k = length(block)

sigma2 = rep(1^2, k)
COV = sapply(1:length(block), function(i) diag( length(block[[i]]) ), simplify = F)
### run real
thin = 1
burn = 1000
niters = 50000 * thin + burn
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
adapIter = 1


ptm = proc.time()[3]
Murray = exactCOMPreg_fixnu(niters, y, X, beta_initial, alpha_initial, block, sigma2, COV, updateCOV, 
                            adaptInterval, adaptFactorExponent, adapIter, thin)
rtime = proc.time()[3] - ptm

postSamp = Murray$Sample
accprob = colMeans(Murray$Accprob)

save(Murray, rtime, postSamp, accprob, file = paste0('simuexchange_50000_dim2_n45_inter_final2.RData'))

load("simuexchange_50000_dim2_n45_inter_final2.RData")


M=1
x0<-rmvnorm(M,mean=c(beta_initial, alpha_initial), sigma=diag(4))
x0<-as.matrix(x0)
pt <- proc.time()
stein_comp_svi<-comp_fixnu_map(x0, X, y, 300, 10, 0.005, exp(0.5))
ptFinal_glm <- proc.time()-pt
ptFinal_glm
map<-as.matrix(colMeans(stein_comp_svi[291:300,]))

colnames(X)<-c("x1","x2","x3")
datas<-cbind(y,X)
datas<-data.frame(datas)
glm_re<-glm.cmp(y~x1+x2, data=datas)
cov1<-vcov(glm_re)[2:3, 2:3]
cov2<-rbind(cov1,vcov(glm_re)[1, 2:3])
cov3<-cbind(cov2,vcov(glm_re)[c(2,3,1),1])
cov4<-rbind(cov3,vcov(glm_re)[4,c(2,3,1)])
cov5<-cbind(cov4, vcov(glm_re)[c(2,3,1,4),4])

M=96
x0<-rmvnorm(M,mean=map, sigma=cov5)
x0<-as.matrix(x0)


pt <- proc.time()
stein_comp_svi<-comp_svi_fixnu_par2(x0, X, y, map, 500, 50, 50, 0.0001, 3, exp(0.5),32)
ptFinal_glm <- proc.time()-pt
ptFinal_glm

pos_in<-matrix(rep(0,9), nrow=3)
pos_in[,1]<-colMeans(stein_comp_svi$theta)[1:3]
library(HDInterval)
for(i in 1:3){
  pos_in[i,2:3]<-hdi(stein_comp_svi$theta[,i])
}
save(stein_comp_svi,ptFinal_glm, pos_in,  file="comp_svi_dim2_n45_inter_par96_iter500_impor50_50_thres3_step0.0001_cov_glmcmp_final_new.RData")
par(mfrow=c(1,3))
for(jj in 1:3){
  plot(density(postSamp[1001:50000,jj]), col="red")
  lines(density(stein_comp_svi$theta[,jj]),col="black")
  lines(density(x0[,jj]), col="blue")
}








rm(list=ls())
library(Rcpp);library(RcppArmadillo);library(MASS)
library(rdist);library(fields);library(lattice)
library(Bergm);library(coda);library(xtable)
library(DiceKriging);library(DiceDesign)
library(dplyr);library(sp);library(gstat)
library(classInt);library(spatstat)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

sourceCpp('MCSVGD_COMP.cpp')





set.seed(2025)
###########################   n=50   ############################################################



load("simData2dim_n50_inter.RData")
betaInd = 1:p
beta_initial = rnorm(p)
alpha_initial = log(1) # log(nu) = alpha

b = c(rep(1, p), 2)
block = sapply(1:length(unique(b)), function(i) which(b == i) - 1, simplify = F)
k = length(block)

sigma2 = rep(1^2, k)
COV = sapply(1:length(block), function(i) diag( length(block[[i]]) ), simplify = F)
### run real
thin = 1
burn = 1000
niters = 50000 * thin + burn
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
adapIter = 1


ptm = proc.time()[3]
Murray = exactCOMPreg_fixnu(niters, y, X, beta_initial, alpha_initial, block, sigma2, COV, updateCOV, 
                            adaptInterval, adaptFactorExponent, adapIter, thin)
rtime = proc.time()[3] - ptm

postSamp = Murray$Sample
accprob = colMeans(Murray$Accprob)

save(Murray, rtime, postSamp, accprob, file = paste0('simuexchange_50000_dim2_n50_inter_final2.RData'))

load("simuexchange_50000_dim2_n50_inter_final2.RData")


M=1
x0<-rmvnorm(M,mean=c(beta_initial, alpha_initial), sigma=diag(4))
x0<-as.matrix(x0)
pt <- proc.time()
stein_comp_svi<-comp_fixnu_map(x0, X, y, 300, 10, 0.005, exp(0.5))
ptFinal_glm <- proc.time()-pt
ptFinal_glm
map<-as.matrix(colMeans(stein_comp_svi[291:300,]))

colnames(X)<-c("x1","x2","x3")
datas<-cbind(y,X)
datas<-data.frame(datas)
glm_re<-glm.cmp(y~x1+x2, data=datas)
cov1<-vcov(glm_re)[2:3, 2:3]
cov2<-rbind(cov1,vcov(glm_re)[1, 2:3])
cov3<-cbind(cov2,vcov(glm_re)[c(2,3,1),1])
cov4<-rbind(cov3,vcov(glm_re)[4,c(2,3,1)])
cov5<-cbind(cov4, vcov(glm_re)[c(2,3,1,4),4])
cov5

M=96
x0<-rmvnorm(M,mean=map, sigma=cov5)
x0<-as.matrix(x0)


pt <- proc.time()
stein_comp_svi<-comp_svi_fixnu_par2(x0, X, y, map, 500, 50, 50, 0.0001, 3, exp(0.5),32)
ptFinal_glm <- proc.time()-pt
ptFinal_glm


pos_in<-matrix(rep(0,9), nrow=3)
pos_in[,1]<-colMeans(stein_comp_svi$theta)[1:3]
library(HDInterval)
for(i in 1:3){
  pos_in[i,2:3]<-hdi(stein_comp_svi$theta[,i])
}
save(stein_comp_svi,ptFinal_glm, pos_in,   file="comp_svi_dim2_n50_inter_par96_iter500_impor50_50_thres3_step0.0001_cov_glmcmp_final_new.RData")
par(mfrow=c(1,3))
for(jj in 1:3){
  plot(density(postSamp[1001:50000,jj]), col="red")
  lines(density(stein_comp_svi$theta[,jj]),col="black")
  lines(density(x0[,jj]), col="blue")
}
