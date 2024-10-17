rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(MASS)

library(rdist)
library(fields)
library(lattice)
library(Bergm)
library(coda)
library(xtable)
library(DiceKriging)
library(DiceDesign)
library(dplyr)
library(sp)
library(gstat)
library(classInt)
library(spatstat)

sourceCpp('MCSVGD_COMP.cpp')


rm(list=ls())
library(ngspatial); library(MASS); library(foreach)
library(Rcpp); library(RcppArmadillo)

sourceCpp('RcppFtns.cpp')
# ==============================================================================
# data simulation
# ==============================================================================
set.seed(1012)

### adjacency matrix and design matrix
n = 10
A = adjacency.matrix(n)
x = y = seq(0, 1, length.out = n)
coord = cbind(rep(x, times = n), rep(y, each = n))
t = 1
X = foreach(j = 1:t, .combine = 'rbind') %do% { coord }
p = ncol(X)+1; N = nrow(X)
n = rep(nrow(A), t); nt = c(0, cumsum(n))

inter<-rep(1,N)
X<-cbind(X,inter)
#load("automul_mpleresult.RData")
### true parameters
beta=c(1,1,0.1)
#beta = c(0.1,-0.1,0.3,0.5,0.05,0.5,0.3,0.2,0.01,-0.02)
alpha = 0.5
trpar = list(beta = beta, alpha = alpha)


### detection probability, mode, and dispersion parameters
mode = exp( as.vector( X %*% beta ) )
nu = exp(alpha)


### generate outcome
y = sapply(1:N, function(i)  rCOMP(1, mode[i], nu))


save(coord, t, n, nt, y, X, N, p, trpar, mode, nu, file = paste0("simData2dim_n10_inter.RData"))


