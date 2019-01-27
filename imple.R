library(copula)
library(ReIns)
library(tailDepFun)
library(cubature)
setwd('E:/Files/Toronto Courses/RA')
source('MEstimatorSummary.R')
source('factor.R')
#########################################################
n = 1500
d = 2
theta0 = 0.5
k=40
##used to generate data
gumbel = gumbelCopula(param=1/theta0, dim=d)
md= mvdc(gumbel,margins = c('frechet','frechet'),
         paramMargins = list(list(loc=0,scale=1,shape=1),list(loc=0,scale=1,shape=1)))


##################single parameter estimation#########################
kv = seq(40,320,40)
estimation = matrix(NaN,nrow = 200, ncol = length(kv))
estimation_x1 = matrix(NaN,nrow = 200, ncol = length(kv))
estimation_xsum = matrix(NaN,nrow = 200, ncol = length(kv))
estimation_xsumSquare = matrix(NaN,nrow = 200, ncol = length(kv))
for (i in 1:length(kv)) {
  k = kv[i];
  for (j in 1:200) {
    set.seed(j)
    da = rMvdc(n,md)
    estimation[j,i] = Mestimator_Gumbel(da, k)
    estimation_x1[j,i] = Mestimator_Gumbel_x1(da, k)
    estimation_xsum[j,i] = Mestimator_Gumbel_xsum(da,k)
    estimation_xsumSquare[j,i] = Mestimator_Gumbel_xsquare_sum(da,k)
  }
}
### for Censured MLE ###
kv1 = round(seq(8,100,length.out=8))
estimation_MLE = matrix(NaN,nrow = 200, ncol = length(kv1))
#estimation_OneSmallk = matrix(NaN,nrow = 200, ncol = length(kv1))
for (i in 1:length(kv)) {
  k = kv1[i];
  for (j in 1:200) { # 56 
    set.seed(j)
    if(j>=56)
      set.seed(j+1)
    da = rMvdc(n,md)
    model = POT::fitbvgpd(da, rep(k,2),model = 'log',cscale = T,cshape = T)
    estimation_MLE[j,i] = model$param['alpha']
    #estimation_OneSmallk[j,i] = Mestimator_Gumbel(da, k)
  }
}

####################five dimension###############################
d =5
gumbel = gumbelCopula(param=1/theta0, dim=d)
md= mvdc(gumbel,margins = c('frechet','frechet','frechet','frechet','frechet'),
         paramMargins = list(list(loc=0,scale=1,shape=1),list(loc=0,scale=1,shape=1),list(loc=0,scale=1,shape=1),
                             list(loc=0,scale=1,shape=1),list(loc=0,scale=1,shape=1)))
estimation_5d = matrix(NaN,nrow = 200, ncol = length(kv))

for (i in 1:length(kv)) {
  k = kv[i];
  for (j in 1:200) {
    set.seed(j)
    da = rMvdc(n,md)
    estimation_5d[j,i] = Mestimator_Gumbel(da, k)
  }
}


######################factor Model###################
source('factor.R')
###############generate data#########################
n= 5000
theta = c(0.2,0.5,0.7,0.9)
threshold = 40
Z1 = rfrechet(n,shape = 1)
Z2 = rfrechet(n,shape = 1)
A1 = matrix(theta,nrow = 4,ncol = 1)
A2 = 1 - A1
A = cbind(A1,A2)

tem1 = A1 %*% t(Z1)
tem2 = A2 %*% t(Z2)
X = t(ifelse(tem1>=tem2, tem1,tem2)) ##rows n, dimension d

kv = c(80,round(seq(40,280,40) /1500 * 5000))
CoefS = list()
for (j in 1:length(kv)) {
  threshold = kv[j]
  coefs = c()
  for (i in 1:200) {
    set.seed(i*3)
    Z1 = rfrechet(n,shape = 1)
    set.seed(i*7)
    Z2 = rfrechet(n,shape = 1)
    A1 = matrix(theta,nrow = 4,ncol = 1)
    A2 = 1 - A1
    A = cbind(A1,A2)
    tem1 = A1 %*% t(Z1)
    tem2 = A2 %*% t(Z2)
    X = t(ifelse(tem1>=tem2, tem1,tem2))
    int1 = c()
    int1[1] = InteStdf_gx(X,k = threshold,gx= 'x',order = 1)
    int1[2] = InteStdf_gx(X,k = threshold,gx= 'x',order = 2)
    int1[3] = InteStdf_gx(X,k = threshold,gx= 'x',order = 3)
    int1[4] = InteStdf_gx(X,k = threshold,gx= 'x',order = 4)
    int1[5] = InteStdf_gx(X,k = threshold,gx= '1',order = 1)
    Matrix_Mestimator = function(V, X,threshold,r){
      d = ncol(X) 
      sum_square = (factorTailDepIntegral(V,1,0,d,r) - int1[5])^2
      for (j in 1:d) {
        sum_square = sum_square + 
          (factorTailDepIntegral(V,j,1,d,r)-int1[j])^2
      }
      return( sum_square )
    }
    
    result = stats::optim(c(0.2,0.5,0.7,0.9),Matrix_Mestimator,X= X, r= 2,threshold = threshold,lower = 0.05,upper = 0.99,
                          method = 'L-BFGS-B')$par
    print(paste(j,i,sep='---'))
    print(result)
    coefs = rbind( coefs, result )
  }
  CoefS[[j]] = coefs
}


