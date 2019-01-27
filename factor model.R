library(copula)
library(cubature)
library(ReIns)
library(tailDepFun)

setwd('E:/Files/Toronto Courses/RA')
source('MEstimatorSummary.R')
#####max factor model####
fx = function(x, vector,j,k,s,d){
  product = ( min( vector[j]/vector[k] * x,1) )^s
  for (l in 1:d) {
    product = product * ( min( vector[j]/vector[l] * x,1) )
  }
  return(product)
}
factorTailDepIntegral = function(V,k,s,d,r){
  ##k,s determines g(x) as x[k]^s, V is a vector with length d*(r-1), 
  ##coefficients are in V, d is row, r is column
  #V1 = c()
  #for (i in 1:d) {
  #    V1 = c(V1,1 - sum(V[numbers::mod(1:length(V),d)==i]))
  #}
  #V1[d] = 1 - sum(V[numbers::mod(1:length(V),d)==0])
  V1 = 1-V
  B = c(V,V1)
  Sum = 0
  for (i in 1:r) {
    Ci = B[((i-1)*d+1) : (i * d)]
    for (j in 1:d) {
      integral = adaptIntegrate(fx,lower = 0,upper = 1,vector=Ci,j=j,k=k,s=s,d=d,maxEval = 500)$integral
      if(j == k){
        Sum = Sum + Ci[j] * integral
      }
      else
        Sum = Sum + Ci[j]/(1+s) * integral
    }
  }
  return(Sum)
}

##very slow by putting empirical integral function in it## I change it by separating the integral out in the simulation
Matrix_Mestimator = function(V, X,threshold,r){
  d = ncol(X) 
  sum_square = 0
  for (j in 1:d) {
    sum_square = sum_square + 
      (factorTailDepIntegral(V,j,1,d,r)-InteStdf_gx(X,k = threshold,gx= 'x',order = j))^2
  }
  sum_square  = sum_square + 
    ( factorTailDepIntegral(V,1,0,d,r)-InteStdf_gx(X,k = threshold,gx= '1') )^2

  return( sum_square )
}
##
haha = function(X, threshold,ini,r,accuracy){
  ini = ini
  temp = stats::optim(ini,Matrix_Mestimator,X= X, r= r,threshold = threshold,lower = 0,upper = 1,
                      method = 'L-BFGS-B')
  return(temp$par)
}


  