library(cubature)
library(copula)
library(ReIns)

##Empirical tail dependence function
Emp_stdf = function(da, k, cst){
  rank = apply(da,2,rank, ties.method = 'first')
  d = ncol(rank); n = nrow(rank)
  if(length(cst)!=d){
    warning('The length of cst')
  }
  bool_matrix = t(apply(rank, 1, '>', n + 0.5 - k * cst))
  temp = apply(bool_matrix, 1, any)
  return(sum(temp)/k)
}

#####compute empirical tail function integral with different g(x)#################
InteStdf_gx= function(da,k,gx,order=1){
  ## gx = '1' or 'x' or 'sum' or 'sum_square'
  ranks = apply(da,2,rank)
  d = ncol(da);n=nrow(da)
  if(gx == '1')
    integrand1 = function(x) return(1)
  else if(gx == 'x')
    integrand1 = function(x,th=order) return(x[th])
  else if(gx == 'sum')
    integrand1 = function(x) return(sum(x))
  else if(gx == 'sum_square')
    integrand1 = function(x) return(sum(x^2))
  Sum = 0
  for (i in 1:n) {
    upp = (n + 0.5 - ranks[i,])/k
    if(any(upp<=1)){
      upper = ifelse(rep(1,d)>upp,upp,1 )
      integral1 = adaptIntegrate(integrand1, lowerLimit = rep(0,d),upperLimit = rep(1,d))$integral
      Sum = Sum + integral1 - adaptIntegrate(integrand1, lowerLimit = rep(0,d),upperLimit = upper)$integral
    }
  }
  return(Sum/k)
}

############## M Estimator for g(x) = 1####################
LogisTailFunction=function(x,theta){
  return ((sum(x^(1/theta)))^theta) ## x is an vector
}

Gumbel_integral = function(da,params, k){
  d = ncol(da)
  result = cubature::adaptIntegrate(LogisTailFunction, lowerLimit = rep(0,d),upperLimit = rep(1,d),
                                    theta = params,absError = 10^-4)$integral
  result_emp = InteStdf_gx(da,k,'1')
  result = (result- result_emp)^2
  return(result)
}

Mestimator_Gumbel = function(da, k){
  #x is the data p*d matrix
  
  theta = stats::optimize(Gumbel_integral,interval = c(0.05,0.99), da=da, k=k, tol= 10^-4)
  
  return(theta$minimum)
}

################## M Estimator for g(x) = x####################
Gumbel_integral_x1 = function(params, da, k){
  d = ncol(da)
  GumbelIntegrand = function(x,theta) return(x[1] *((sum(x^(1/theta)))^theta))
  result_x1 = cubature::adaptIntegrate(GumbelIntegrand, lowerLimit = rep(0,d),upperLimit = rep(1,d),
                                       theta = params,absError = 10^-4)$integral
  emp_result = InteStdf_gx(da,k,'x')
  result = (result_x1 - emp_result )^2
  return(result)
}

Mestimator_Gumbel_x1 = function(da, k){
  #x is the data p*d matrix
  theta = stats::optimize(Gumbel_integral_x1,interval = c(0.05,0.99), da = da, k = k,tol = 10^-4)
  
  return(theta$minimum)
}

#################### M Estimator for g(x) = sum(x)####################
Gumbel_integral_xsum = function(params, da, k){
  d = ncol(da)
  GumbelIntegrand_xsum = function(x,theta) return(sum(x) *((sum(x^(1/theta)))^theta))
  result_x1 = cubature::adaptIntegrate(GumbelIntegrand_xsum, lowerLimit = rep(0,d),upperLimit = rep(1,d),
                                       theta = params,absError = 10^-4)$integral
  emp_result = InteStdf_gx(da,k,'sum')
  result = (result_x1 - emp_result )^2
  return(result)
}

Mestimator_Gumbel_xsum = function(da, k){
  #x is the data p*d matrix
  theta = stats::optimize(Gumbel_integral_xsum,interval = c(0.05,0.99), da = da, k = k,tol = 10^-4)
  
  return(theta$minimum)
}
#################### M Estimator for g(x) = sum(x^2)####################
Gumbel_integral_xsquare_sum = function(params, da, k){
  d = ncol(da)
  GumbelIntegrand_xsquare_sum = function(x,theta) return(sum(x^2) * ((sum(x^(1/theta)))^theta))
  result_x1 = cubature::adaptIntegrate(GumbelIntegrand_xsquare_sum, lowerLimit = rep(0,d),upperLimit = rep(1,d),
                                       theta = params,absError = 10^-4)$integral
  emp_result =  InteStdf_gx(da,k,'sum_square')
  result = (result_x1 - emp_result )^2
  return(result)
}

Mestimator_Gumbel_xsquare_sum = function(da, k){
  #x is the data p*d matrix
  theta = stats::optimize(Gumbel_integral_xsquare_sum,interval = c(0.05,0.99), da = da, k = k,tol=10^-4)
  
  return(theta$minimum)
}
