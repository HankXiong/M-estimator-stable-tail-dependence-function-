library(graphics)
########################some plot 2d########################################
RMSE = function(series, realTheta){
  return( sqrt(mean((series-realTheta)^2))  )
}
bias = function(series,realTheta){
  return( mean(series-realTheta) )
}
par(mar=c(4,4,1,0.5))

rmse_estimation = apply(estimation,2,RMSE, realTheta = 0.5)
rmse_estimation_x1 = apply(estimation_x1,2,RMSE,realTheta = 0.5)
rmse_estimation_xsum = apply(estimation_xsum,2,RMSE,realTheta = 0.5)
rmse_estimation_xsquare_sum = apply(estimation_xsumSquare,2,RMSE,realTheta = 0.5)
bias_estimation = apply(estimation,2,bias,realTheta = 0.5)
bias_estimation_x1 = apply(estimation_x1,2,bias,realTheta = 0.5)
bias_estimation_xsum = apply(estimation_xsum,2,bias,realTheta = 0.5)
bias_estimation_xsquare_sum = apply(estimation_xsumSquare,2,bias,realTheta = 0.5)

plot(kv,rmse_estimation,type='b',xlab = 'Threshold',ylab='RMSE',lty = 1,pch=1,ylim = c(0.035,0.1))
lines(kv,rmse_estimation_x1,type = 'b',lty=2,pch=2)
lines(kv,rmse_estimation_xsum,type = 'b',lty=3,pch=3)
lines(kv,rmse_estimation_xsquare_sum,type = 'b',lty=4,pch=4)
lines(kv[1:4],rmse_mle,type = 'b',col = 'green',lty=5,pch=5)
legend('topright',legend=c('g(x)=1','g(x)= x[1]','sum x','sum x^2'),lty=1:4,pch=1:4)

plot(kv,bias_estimation,type='b',xlab = 'Threshold',ylab='Bias',lty = 1,pch=1,ylim=c(-0.05,0.02))
lines(kv,bias_estimation_x1,type = 'b',lty=2,pch=2)
lines(kv,bias_estimation_xsum,type = 'b',lty=3,pch=3)
lines(kv,bias_estimation_xsquare_sum,type = 'b',lty=4,pch=4)
lines(kv[1:4],bias_mle,type = 'b',col = 'green', lty=5,pch =5)
legend('topright',legend=c('g(x)=1','g(x)= x[1]','sum x','sum x^2'),lty=1:4,pch=1:4)

rmse_mle = apply(estimation_MLE[,1:4],2,RMSE,realTheta = 0.5)
bias_mle = apply(estimation_MLE[,1:4],2,bias,realTheta = 0.5)

##smaller k##
rmse_mle = apply(estimation_MLE,2,RMSE,realTheta=0.5)
bias_mle = apply(estimation_MLE,2,bias,realTheta=0.5)
rmse_smaller = apply(estimation_OneSmallk,2,RMSE,realTheta=0.5)
bias_smaller = apply(estimation_OneSmallk,2,bias,realTheta=0.5)

plot(kv1,bias_mle,type='b',xlab = 'Threshold',ylab='Bias',lty = 1,pch=1,ylim=c(-0.01,0.01))
lines(kv1,bias_smaller,type = 'b',lty=2,pch=2)
legend('bottomleft',legend=c('Censured MLE','g(x)=1'),lty=1:2,pch=1:2)

plot(kv1,rmse_mle,type='b',xlab = 'Threshold',ylab='RMSE',lty = 1,pch=1)
lines(kv1,rmse_smaller,type = 'b',lty=2,pch=2)
legend('bottom',legend=c('Censured MLE','g(x)=1'),lty=1:2,pch=1:2)


########################plot 5d###########################################
par(mfrow=c(2,2),mar=c(2,4,2,0.5))
rmse_5d = apply(estimation_5d,2,RMSE,realTheta = 0.5)
bias_5d = apply(estimation_5d,2,bias,realTheta = 0.5)

emp_tail_FiveOne = matrix(NaN,nrow =200, ncol = 8)
for (i in 1:length(kv)) {
  k = kv[i];
  for (j in 1:200) {
    set.seed(j)
    da = rMvdc(n,md)
    emp_tail_FiveOne[j,i] = Emp_stdf(da,k,c(1,1,1,1,1))
  }
}
emp_tail_FiveOne_bias = apply(emp_tail_FiveOne, 2, bias,realTheta = sqrt(5))
emp_tail_FiveOne_RMSE = apply(emp_tail_FiveOne, 2, RMSE,realTheta = sqrt(5))

d5_tail_value = 5^estimation_5d
d5_tail_FiveOne_bias = apply(d5_tail_value, 2, bias,realTheta = sqrt(5))
d5_tail_FiveOne_RMSE = apply(d5_tail_value, 2, RMSE,realTheta = sqrt(5))
##
plot(kv,bias_5d,type= 'b',ylab = 'Bias',xlab = 'k',main = 'Mestimator 5d')
plot(kv,rmse_5d, type= 'b',ylab = 'RMSE',xlab ='k',main = 'Mestimator 5d')

plot(kv,emp_tail_FiveOne_bias,type='b',lty = 2,ylim=c(-0.3,0),ylab ='Bias',
     main = expression(paste(plain(theta),'=0.5, ',plain(l),'(1,1,1,1,1;',plain(theta),') = ',sqrt(5))),xlab = 'k')
lines(kv,d5_tail_FiveOne_bias,type='b',lty=1)
legend('topright', legend= c(expression(paste(plain(l),'(1,1,1,1,1;',plain(hat(theta)),')')),expression(paste(plain(hat(l)),'(1,1,1,1,1)'))),lty=1:2)

plot(kv,emp_tail_FiveOne_RMSE,type='b',lty =2,ylim=c(0.12,0.3),ylab ='RMSE',
     main = expression(paste(plain(theta),'=0.5, ',plain(l),'(1,1,1,1,1;',plain(theta),') = ',sqrt(5))),xlab = 'k')
lines(kv,d5_tail_FiveOne_RMSE,type='b',lty=1)
legend('topleft', legend= c(expression(paste(plain(l),'(1,1,1,1,1;',plain(hat(theta)),')')),expression(paste(plain(hat(l)),'(1,1,1,1,1)'))),lty=1:2)

###############Factor Model plot #################
fac1 = matrix(NA,nrow = 200,ncol = 8)
fac2 = matrix(NA,nrow = 200,ncol = 8)
fac3 = matrix(NA,nrow = 200,ncol = 8)
fac4 = matrix(NA,nrow = 200,ncol = 8)
for (i in 1:8) {
  fac1[,i] = CoefS[[i]][,1]
  fac2[,i] = CoefS[[i]][,2]
  fac3[,i] = CoefS[[i]][,3]
  fac4[,i] = CoefS[[i]][,4]
}
bias1 = apply(fac1,2,bias,realTheta = 0.2)
bias2 = apply(fac2,2,bias,realTheta = 0.5)
bias3 = apply(fac3,2,bias,realTheta = 0.7)
bias4 = apply(fac4,2,bias,realTheta = 0.9)
RMSE1 = apply(fac1,2,RMSE,realTheta = 0.2)
RMSE2 = apply(fac2,2,RMSE,realTheta = 0.5)
RMSE3 = apply(fac3,2,RMSE,realTheta = 0.7)
RMSE4 = apply(fac4,2,RMSE,realTheta = 0.9)
kv_factor = c(80,round(seq(40,280,40) /1500 * 5000))

par(mfrow=c(2,2),mar=c(4,4,1,0.5))
plot(kv_factor,bias1,type='b',ylab = 'Bias',ylim=c(-0.01,0.038),xlab = 'k',
     main = expression( paste(plain(theta[1]),'=0.2, ',plain(theta[2]),'=0.5') ))
lines(kv_factor,bias2,type='b',lty=2,pch=2)
legend('topleft',legend=c(expression(plain(theta[1])),expression(plain(theta[2]))),lty=1:2,pch=1:2)

plot(kv_factor,bias3,type='b',ylab = 'Bias',ylim=c(-0.04,0.003),xlab='k',
     main = expression( paste(plain(theta[3]),'=0.7, ',plain(theta[4]),'=0.9') ) )
lines(kv_factor,bias4,type='b',lty=2,pch=2)
legend('topright',legend=c(expression(plain(theta[3])),expression(plain(theta[4]))),lty=1:2,pch=1:2)

plot(kv_factor,RMSE1,type='b',ylab = 'RMSE',ylim=c(0.00,0.05),xlab = 'k',
     main = expression( paste(plain(theta[1]),'=0.2, ',plain(theta[2]),'=0.5') ))
lines(kv_factor,RMSE2,type='b',lty=2,pch=2)
legend('topleft',legend=c(expression(plain(theta[1])),expression(plain(theta[2]))),lty=1:2,pch=1:2)

plot(kv_factor,RMSE3,type='b',ylab = 'RMSE',ylim=c(0.00,0.05),xlab='k',
     main = expression( paste(plain(theta[3]),'=0.7, ',plain(theta[4]),'=0.9') ) )
lines(kv_factor,RMSE4,type='b',lty=2,pch=2)
legend('topright',legend=c(expression(plain(theta[3])),expression(plain(theta[4]))),lty=1:2,pch=1:2)

