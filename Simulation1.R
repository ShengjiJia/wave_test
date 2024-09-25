library(lars)

###Simulation1
######true null distribution
set.seed(9)
True=-log(-log(runif(100000, min=0, max=1)))
q=quantile(True,prob=0.95)
######data generating
n=500
x=1:n
signal=1*(x>=0.3*n)-2*(x>=0.4*n)+2*(x>=0.8*n)-1*(x>=0.9*n)
#matrix X
x1=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1))
  x1[i,(i+1):n]=0
TT0=1:500
TT1=1:500
TT2=1:500
prob=1:500          #probability cotaining 450 for s=5
power0=matrix(0,nrow=500,ncol=5)
power1=matrix(0,nrow=500,ncol=5)
power2=matrix(0,nrow=500,ncol=5)
beta=c(0.025,0.05,0.075,0.10,0.125)
for(i in 1:500){        
  e=rnorm(n,mean=0,sd=0.5)     
  fi=runif(1, min = 0, max = 2*pi)
  thi=runif(1, min = 0, max = 2*pi)
  wave=sin(2*pi*x/96+fi)+2*sin(2*pi*x/240+thi) 
  y0=signal+e              #under Null
  #oracle test
  r0=lm(y0~x1[,c(0.3*n,0.4*n,0.8*n,0.9*n)])$residuals
  #Fourier transform
  fft0=NULL
  for(j in 1:(n/2)){
    a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r0)
    b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r0)
    fft0=c(fft0, a, b)
  }
  v0=var(fft0[floor(n/4):n])        #sigma^2
  T0=NULL
  for(m in 1:n){
    t=(2*m*v0^2)^(-0.5)*sum((fft0[1:m])^2-v0)
    T0=c(T0,t)
  }
  TT0[i]=sqrt(2*log(log(n)))*max(T0)-2*log(log(n))-0.5*log(log(log(n)))+0.5*log(4*pi)            #test statistic
  for(k in 1:5){
    y1=y0+beta[k]*wave        #Under alternative
    #oracle
    r=lm(y1~x1[,c(0.3*n,0.4*n,0.8*n,0.9*n)])$residuals
    #Fourier transform
    fft=NULL
    for(j in 1:(n/2)){
      a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r)
      b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r)
      fft=c(fft, a, b)
    }
    v=var(fft[floor(n/4):n])        #sigma^2
    T=NULL
    for(m in 1:n){
      t=(2*m*v^2)^(-0.5)*sum((fft[1:m])^2-v)
      T=c(T,t)
    }
    power0[i,k]=(sqrt(2*log(log(n)))*max(T)-2*log(log(n))-0.5*log(log(log(n)))+0.5*log(4*pi) >q)         #test statistic
  }
  
  #lasso with s=10
  model1=lars(x1[1:n,2:n], y =y0, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=10)
  sub=which((which(model1$entry>0)[2:10]-which(model1$entry>0)[1:9])==1)
  if (length(sub)==0)
    r1=lm(y0~x1[,(which(model1$entry>0)+1)])$residuals
  if (length(sub)>0)
    r1=lm(y0~x1[,(which(model1$entry>0)[-sub]+1)])$residuals
  #Fourier transform
  fft1=NULL
  for(j in 1:(n/2)){
    a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r1)
    b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r1)
    fft1=c(fft1, a, b)
  }
  v1=var(fft1[floor(n/4):n])        #sigma^2
  T1=NULL
  for(m in 1:n){
    t=(2*m*v1^2)^(-0.5)*sum((fft1[1:m])^2-v1)
    T1=c(T1,t)
  }
  TT1[i]=sqrt(2*log(log(n)))*max(T1)-2*log(log(n))-0.5*log(log(log(n)))+0.5*log(4*pi)           #test statistic
  for(k in 1:5){
    y1=y0+beta[k]*wave        #Under alternative
    #lasso with s=10
    model=lars(x1[1:n,2:n], y =y1, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=10)
    sub=which((which(model$entry>0)[2:10]-which(model$entry>0)[1:9])==1)
    if (length(sub)==0)
      r=lm(y1~x1[,(which(model$entry>0)+1)])$residuals
    if (length(sub)>0)
      r=lm(y1~x1[,(which(model$entry>0)[-sub]+1)])$residuals
    #Fourier transform
    fft=NULL
    for(j in 1:(n/2)){
      a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r)
      b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r)
      fft=c(fft, a, b)
    }
    v=var(fft[floor(n/4):n])        #sigma^2
    T=NULL
    for(m in 1:n){
      t=(2*m*v^2)^(-0.5)*sum((fft[1:m])^2-v)
      T=c(T,t)
    }
    power1[i,k]=(sqrt(2*log(log(n)))*max(T)-2*log(log(n))-0.5*log(log(log(n)))+0.5*log(4*pi)>q)         #test statistic
  }
  
  #lasso with s=5
  model2=lars(x1[1:n,2:n], y =y0, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=5)
  sub=which((which(model2$entry>0)[2:5]-which(model2$entry>0)[1:4])==1)
  if (length(sub)==0)
    r2=lm(y0~x1[,(which(model2$entry>0)+1)])$residuals
  if (length(sub)>0)
    r2=lm(y0~x1[,(which(model2$entry>0)[-sub]+1)])$residuals
  prob[i]=sum((which(model2$entry>0)<453)*(which(model2$entry>0)>447))
  #Fourier transform
  fft2=NULL
  for(j in 1:(n/2)){
    a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r2)
    b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r2)
    fft2=c(fft2, a, b)
  }
  v2=var(fft2[floor(n/4):n])        #sigma^2
  T2=NULL
  for(m in 1:n){
    t=(2*m*v2^2)^(-0.5)*sum((fft2[1:m])^2-v2)
    T2=c(T2,t)
  }
  TT2[i]=sqrt(2*log(log(n)))*max(T2)-2*log(log(n))-0.5*log(log(log(n)))+0.5*log(4*pi)            #test statistic
  for(k in 1:5){
    y1=y0+beta[k]*wave        #Under alternative
    #lasso with s=5
    model=lars(x1[1:n,2:n], y =y1, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=5)
    sub=which((which(model$entry>0)[2:5]-which(model$entry>0)[1:4])==1)
    if (length(sub)==0)
      r=lm(y1~x1[,(which(model$entry>0)+1)])$residuals
    if (length(sub)>0)
      r=lm(y1~x1[,(which(model$entry>0)[-sub]+1)])$residuals
    #Fourier transform
    fft=NULL
    for(j in 1:(n/2)){
      a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r)
      b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r)
      fft=c(fft, a, b)
    }
    v=var(fft[floor(n/4):n])        #sigma^2
    T=NULL
    for(m in 1:n){
      t=(2*m*v^2)^(-0.5)*sum((fft[1:m])^2-v)
      T=c(T,t)
    }
    power2[i,k]=(sqrt(2*log(log(n)))*max(T)-2*log(log(n))-0.5*log(log(log(n)))+0.5*log(4*pi) >q)         #test statistic
  }
}

par(mfrow=c(2,3))
plot(quantile(True,probs=(1:199)/200),quantile(TT0,probs=(1:199)/200),ylim=c(-5,5),xlab="Theoritical Quantiles", ylab="Sample Quantiles",main="Oracle")
lines(quantile(True,probs=(1:199)/200),quantile(True,probs=(1:199)/200),lwd=2)
plot(quantile(True,probs=(1:199)/200),quantile(TT1,probs=(1:199)/200),ylim=c(-5,5),xlab="Theoritical Quantiles", ylab="Sample Quantiles",main="Overestimate")
lines(quantile(True,probs=(1:199)/200),quantile(True,probs=(1:199)/200),lwd=2)
plot(quantile(True,probs=(1:199)/200),quantile(TT2,probs=(1:199)/200),xlab="Theoritical Quantiles", ylab="Sample Quantiles",main="Underestimate")
lines(quantile(True,probs=(1:199)/200),quantile(True,probs=(1:199)/200),lwd=2)
plot(x=c(0,beta),y=c(sum(TT0>q)/500, mean(power0[,1]), mean(power0[,2]), mean(power0[,3]), mean(power0[,4]), mean(power0[,5])), 
     xlab="theta",ylab="power",ylim=c(0,1),main="Oracle", type="b")
plot(x=c(0,beta),y=c(sum(TT1>q)/500, mean(power1[,1]), mean(power1[,2]), mean(power1[,3]), mean(power1[,4]), mean(power1[,5])), 
     xlab="theta",ylab="power",ylim=c(0,1),main="Overestimate", type="b")
plot(x=c(0,beta),y=c(sum(TT2>q)/500, mean(power2[,1]), mean(power2[,2]), mean(power2[,3]), mean(power2[,4]), mean(power2[,5])), 
     xlab="theta",ylab="power",ylim=c(0,1),main="Underestimate", type="b")

rbind(c(sum(TT0>q)/500, mean(power0[,1]), mean(power0[,2]), mean(power0[,3]), mean(power0[,4]), mean(power0[,5])),
      c(sum(TT1>q)/500, mean(power1[,1]), mean(power1[,2]), mean(power1[,3]), mean(power1[,4]), mean(power1[,5])),
      c(sum(TT2>q)/500, mean(power2[,1]), mean(power2[,2]), mean(power2[,3]), mean(power2[,4]), mean(power2[,5])))

sum(prob>0)/500
q
