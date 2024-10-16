########## ERGM exmpales  ##########
## change statistics of gwdegree: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2031865/#FD24 ##
## ergm basics: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2743438/

rm(list=ls())
library(coda);library(Bergm);library(Rcpp)
library(RcppArmadillo);library(DiceKriging);library(DiceDesign)
library(mvtnorm);library(MASS)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

########## Call functions ##########
sourceCpp("MCSVGD_ERGM.cpp")

set.seed(2024)
########## Run step ##########

###  faus.mesa network  ###
data(faux.mesa.high)
faux.mesa.high
plot(faux.mesa.high)


### data processing

grade <- faux.mesa.high %v% "Grade" 
sex <- faux.mesa.high %v% "Sex" 
sex[sex=="M"] = 1
sex[sex=="F"] = 0
sex <- as.numeric(sex)

# summary statistics for data X
X = as.matrix(faux.mesa.high)

# Parameter estimation via MLE or MPLE #
formula <- faux.mesa.high ~ edges + nodematch("Grade",diff=T) + nodematch("Sex",diff=F) + gwdegree(0.25,fixed=TRUE) + gwesp(0.25,fixed=TRUE)
summary(formula)
Xsummary = Summary(X, grade, sex)
Xsummary

m <-ergm(formula,estimate="MPLE")
m$coef
theta <- matrix(m$coef,1)
COV <- solve(-m$hessian)

### 1. Make Design points  via MPLE ###

num.point <- 400
p <- length(Xsummary)
th <- rmvt(num.point, sigma = 3*solve(-m$hessian), df = 8, delta = m$coef )

par(mfrow=c(2,5))
for(i in 1:p){
  plot(density(th[,i]))
}


###### DMH ######
thin = 1
outer = 80000 * thin
adaptInterval = 200
adaptFactorExponent = 0.8
adapIter = 1

pt <- proc.time()
ergmdmh<-ergmDMH_update(X, grade, sex, COV, theta, outer, 10,  adaptInterval, adaptFactorExponent, adapIter, thin)
ptFinal_glm <- proc.time()-pt
ptFinal_glm
save(ergmdmh, file="mesa_ergm_dmh_80000_update_inner10_com.RData")



load("mesa_ergm_dmh_80000_update_inner10.RData")


###map
M=1 #number of particle
x0<-rmvnorm(M, m$coef, COV)
x0<-as.matrix(x0)
pt <- proc.time()
ergm_mesa_map<-ergm_map(x0,X, grade, sex, 500, 10,1, 0.001)
ptFinal_glm <- proc.time()-pt
ptFinal_glm

save(ergm_mesa_map, file="ergm_map_500iter_final.RData")


M=320 #number of particle
x0<-rmvnorm(M, m$coef, COV)
x0<-as.matrix(x0)
#map<-as.matrix(colMeans(ergmdmh))
map<-as.matrix(colMeans(ergm_mesa_map[491:500,]))


pt <- proc.time()
ergm_svi<-ergm_mesa_svi_par_alls(x0, X, grade, sex, map, 500, 5, 50,50, 0.0005,1.1,32,5, 2, 500)
ptFinal_glm <- proc.time()-pt
ptFinal_glm

colMeans(ergm_svi$theta)
library(HDInterval)
for(i in 1:10){
  print(hdi(ergm_svi$theta[,i]))
}

save(ergm_svi,ptFinal_glm, file="ergm_svi_par320_iter1000_impo50_step0.0005_thres1.5_12cycle_10burn_500thin_final.RData")
par(mfrow=c(2,5))
for(jj in 1:10){
  plot(density(ergm_svi$theta[,jj]))
  lines(density(ergmdmh[,jj]), col="red")
}
