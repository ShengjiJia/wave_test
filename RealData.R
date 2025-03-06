library(lars)
library(grpreg)
library(imputeTS)
library(cumSeg)
library(DNAcopy)    #download from https://bioconductor.org/packages/release/bioc/html/DNAcopy.html

####real data 
SNPdata <- read.delim("C:/Users/PC/Desktop/我的文档/Research/Projects/change points(wave test)/SNPdata.txt") 

#################################define some functions
estimateSigma<-function (Y, h = 10) {  
  n = length(Y)
  YBar = rep(0, n)
  for (i in 1:n) {
    a = min(n, i + h)
    b = max(1, i - h)
    YBar[i] = mean(Y[b:a])
  }
  return(sqrt(var(Y - YBar) * (2 * h + 1)/(2 * h)))
}

localDiagnostic<-function (y, h) { 
  yy = c(rep(0, h - 1), y, rep(0, h))
  n = length(y)
  z = rep(0, n)
  for(i in 1:n){
    z[i]=sum(yy[i:(h+i-1)])/h-sum(yy[(h+i):(2*h-1+i)])/h
  }
  return(z)
}

localMax<-function (y, span = 5) {  
  if (length(y) < span * 2 + 1) 
    return(NULL)
  n = length(y)
  index = NULL
  for (i in (span + 1):(n - span)) {
    if (y[i] == max(y[(i - span):(i + span)])) 
      index = c(index, i)
  }
  return(index)
}

MultiScan<-function(y, h=20){
  m=dim(y)[1]
  n=dim(y)[2]
  s=matrix(0,nr=m,nc=n)
  sigma=1:m
  for (i in 1:m){
    sigma[i]=estimateSigma(y[i,], h = max(3, 2*floor(log(n))))
    s[i,]=(localDiagnostic(y[i,],h=h))^2*h/(2*sigma[i])
  }
  S=apply(s,2,sum)
  index=localMax(S,span=2*h)
  return(data.frame(S[index],index))
}

threshold<-function(y, alpha=0.05, h=20) {
  m=dim(y)[1]
  n=dim(y)[2]
  empirical=NULL
  for (k in 1:100) {
    Y=matrix(rnorm(n*m),nr=m)
    s=matrix(0,nr=m,nc=n)
    for (i in 1:m){ s[i,]=(localDiagnostic(Y[i,],h=h))^2*h/2 }
    S=apply(s,2,sum)
    index=localMax(S,span=2*h)
    empirical=c(empirical, S[index])
  }
  return(quantile(empirical, probs=1-alpha))
}

screening<-function(x, y, h=10){       #local linear smoothing (rightlimit-leftlimit)
  n=length(x)
  xx=1:(n+2*h)
  yy=c(rep(y[1],h),y,rep(y[n],h))       #create data outside the boundaries
  right=rep(0,n)     #rightlimit for xx[1:n]
  left=rep(0,n)       #leftlimit for xx[(1+2*h):(n+2*h)]
  for (i in 1:n){
    model=locpoly(xx[i:(i+2*h)],yy[i:(i+2*h)],kernel="epanech",bandwidth=h,gridsize=1+2*h)
    right[i]=model$y[1]
    left[i]=model$y[1+2*h]
  }
  L=c(rep(0,h),right[(2*h+2):n]-left[1:(n-2*h-1)],rep(0,h+1))
  return(L)
}

####Chro 21
y=rbind(SNPdata$X99HI0697A.Log.R.Ratio[which(SNPdata$Chr==21)],
        SNPdata$X99HI0698C.Log.R.Ratio[which(SNPdata$Chr==21)],
        SNPdata$X99HI0700A.Log.R.Ratio[which(SNPdata$Chr==21)]) 
n=ncol(y)
d=nrow(y)
x=1:n
y[which(abs(y)>1.5)]=NA     #delete outliers
for (i in 1:d) {y[i,]=na_ma(y[i,],k=5,weighting="linear")}        #imputation
h=10
num1=1:4        
yy=rbind(SNPdata$X99HI0697A.Log.R.Ratio[which(SNPdata$Chr==22)],
         SNPdata$X99HI0698C.Log.R.Ratio[which(SNPdata$Chr==22)],
         SNPdata$X99HI0700A.Log.R.Ratio[which(SNPdata$Chr==22)])     
rm(SNPdata)
#####CBS
set.seed(2)
CBS=DNAcopy::segment(CNA(t(y), chrom=rep(1,n), maploc=1:n))
estimate1=sort(unique(CBS$output[,4]))
estimate1=estimate1[-length(estimate1)]                
num1[1]=length(estimate1)
######SaRa
sara=MultiScan(y,h=10)
thres=threshold(y,alpha=0.1,h=10)
estimate2=sara$index[which(sara$S.index.>thres)]
num1[2]=length(estimate2)
######group Lasso
y0=as.vector(y)
x1=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1))
  x1[i,(i+1):n]=0            
X=kronecker(x1,diag(1, d))                          #contain intercept
group=NULL
for (i in 1:n)
  group<-c(group,rep(i-1,d))                        #contain intercept    
model1=grpreg(X=X[,2:length(y0)], y=y0, group=group[2:length(y0)], penalty="grLasso",family="gaussian", gmax=10)  
loss=NULL    #residual sum of squares
for(i in 1:length(model1$lambda)){
  loss=c(loss, sum(residuals(model1, lambda=model1$lambda[i])^2))
}
opt=which.min(log(loss/length(y0)) + log(length(y0))*model1$df*0.4*log(log(length(y0)))/length(y0))    
estimate3=as.vector((which(model1$beta[,opt]!=0)[d*(1:(length(which(model1$beta[,opt]!=0))/d))]/d-1)[-1])     
if(length(which(diff(estimate3)<5))>0){estimate33=estimate3[-which(diff(estimate3)<5)]} else{estimate33=estimate3}        
num1[3]=length(estimate33)
rm(model1)
#######partial penalized with m=5
z=matrix(1,nrow=n,ncol=11)                     
for (i in 1:5){
  z[,(2*i)]=cos(2*pi*i*x/n)
  z[,(2*i+1)]=sin(2*pi*i*x/n)             
}
Z=kronecker(z[,2:11],diag(1, d))
X=cbind(X,Z)
group=c(group, rep(0,(10*d)))
model1<-grpreg(X=X[,-1], y=y0, group=group[-1], penalty="grLasso",family="gaussian", gmax=10)
loss=NULL    #residual sum of squares
for(i in 1:length(model1$lambda)){
  loss=c(loss, sum(residuals(model1, lambda=model1$lambda[i])^2))
}
opt=which.min(log(loss/length(y0)) + log(length(y0))*model1$df*0.4*log(log(length(y0)))/length(y0))    
estimate4=as.vector((which(model1$beta[1:(d*n),opt]!=0)[d*(1:(length(which(model1$beta[1:(d*n),opt]!=0))/d))]/d-1)[-1])     
if(length(which(diff(estimate4)<5))>0){estimate44=estimate4[-which(diff(estimate4)<5)]} else{estimate44=estimate4}        
num1[4]=length(estimate44)
rm(model1)
######adaptive Neyman test 
r1=lm(y[1,]~x1[,(estimate3+1)])$residuals
r2=lm(y[2,]~x1[,(estimate3+1)])$residuals
r3=lm(y[3,]~x1[,(estimate3+1)])$residuals
#Fourier transform
fft1=NULL
for(j in 1:floor(n/2)){
  a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r1)
  b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r1)
  fft1=c(fft1, a, b)
}
v1=var(fft1[floor(n/4):length(fft1)])        #sigma^2
fft2=NULL
for(j in 1:floor(n/2)){
  a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r2)
  b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r2)
  fft2=c(fft2, a, b)
}
v2=var(fft2[floor(n/4):length(fft2)])        #sigma^2
fft3=NULL
for(j in 1:floor(n/2)){
  a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r3)
  b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r3)
  fft3=c(fft3, a, b)
}
v3=var(fft3[floor(n/4):length(fft3)])        #sigma^2
T=NULL
for(m in 1:length(fft1)){
  t=sqrt((sum((fft1[1:m])^2-v1)/sqrt(2*m*v1^2))^2+(sum((fft2[1:m])^2-v2)/sqrt(2*m*v2^2))^2+(sum((fft3[1:m])^2-v3)/sqrt(2*m*v3^2))^2)
  T=c(T,t)
}
T1n=sqrt(2*log(log(n)))*max(T)-2*log(log(n))-2*log(log(log(n)))+log(gamma(2)) 

#####Chro 22               
n=ncol(yy)
d=nrow(yy)
x=1:n
y[which(abs(yy)>1.5)]=NA     #delete outliers
for (i in 1:d) {yy[i,]=na_ma(yy[i,],k=5,weighting="linear")}        #imputation
h=10
num2=1:4
#####CBS
set.seed(2)
CBS=DNAcopy::segment(CNA(t(yy), chrom=rep(1,n), maploc=1:n))
Estimate1=sort(unique(CBS$output[,4]))
Estimate1=Estimate1[-length(Estimate1)]                
num2[1]=length(Estimate1)
######SaRa
sara=MultiScan(yy,h=10)
thres=threshold(yy,alpha=0.1,h=10)
Estimate2=sara$index[which(sara$S.index.>thres)]
num2[2]=length(Estimate2)
######group Lasso
y0=as.vector(yy)
x1=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1))
  x1[i,(i+1):n]=0            
X=kronecker(x1,diag(1, d))                          #contain intercept
group=NULL
for (i in 1:n)
  group<-c(group,rep(i-1,d))                        #contain intercept    
model1=grpreg(X=X[,2:length(y0)], y=y0, group=group[2:length(y0)], penalty="grLasso",family="gaussian", gmax=10)  
loss=NULL    #residual sum of squares
for(i in 1:length(model1$lambda)){
  loss=c(loss, sum(residuals(model1, lambda=model1$lambda[i])^2))
}
opt=which.min(log(loss/length(y0)) + log(length(y0))*model1$df*0.4*log(log(length(y0)))/length(y0))    
Estimate3=as.vector((which(model1$beta[,opt]!=0)[d*(1:(length(which(model1$beta[,opt]!=0))/d))]/d-1)[-1])     
if(length(which(diff(Estimate3)<5))>0){Estimate33=Estimate3[-which(diff(Estimate3)<5)]} else{Estimate33=Estimate3}     
num2[3]=length(Estimate33)
rm(model1)
#######partial penalized with m=5
z=matrix(1,nrow=n,ncol=11)                     
for (i in 1:5){
  z[,(2*i)]=cos(2*pi*i*x/n)
  z[,(2*i+1)]=sin(2*pi*i*x/n)             
}
Z=kronecker(z[,2:11],diag(1, d))
X=cbind(X,Z)
group=c(group, rep(0,(10*d)))
model1<-grpreg(X=X[,-1], y=y0, group=group[-1], penalty="grLasso",family="gaussian", gmax=10)
loss=NULL    #residual sum of squares
for(i in 1:length(model1$lambda)){
  loss=c(loss, sum(residuals(model1, lambda=model1$lambda[i])^2))
}
opt=which.min(log(loss/length(y0)) + log(length(y0))*model1$df*0.4*log(log(length(y0)))/length(y0))    
Estimate4=as.vector((which(model1$beta[1:(d*n),opt]!=0)[d*(1:(length(which(model1$beta[1:(d*n),opt]!=0))/d))]/d-1)[-1])     
if(length(which(diff(Estimate4)<5))>0){Estimate44=Estimate4[-which(diff(Estimate4)<5)]} else{Estimate44=Estimate4}        
num2[4]=length(Estimate44)
rm(model1)
######adaptive Neyman test 
r1=lm(yy[1,]~x1[,(Estimate3+1)])$residuals
r2=lm(yy[2,]~x1[,(Estimate3+1)])$residuals
r3=lm(yy[3,]~x1[,(Estimate3+1)])$residuals
#Fourier transform
fft1=NULL
for(j in 1:floor(n/2)){
  a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r1)
  b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r1)
  fft1=c(fft1, a, b)
}
v1=var(fft1[floor(n/4):length(fft1)])        #sigma^2
fft2=NULL
for(j in 1:floor(n/2)){
  a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r2)
  b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r2)
  fft2=c(fft2, a, b)
}
v2=var(fft2[floor(n/4):length(fft2)])        #sigma^2
fft3=NULL
for(j in 1:floor(n/2)){
  a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r3)
  b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r3)
  fft3=c(fft3, a, b)
}
v3=var(fft3[floor(n/4):length(fft3)])        #sigma^2
T=NULL
for(m in 1:length(fft1)){
  t=sqrt((sum((fft1[1:m])^2-v1)/sqrt(2*m*v1^2))^2+(sum((fft2[1:m])^2-v2)/sqrt(2*m*v2^2))^2+(sum((fft3[1:m])^2-v3)/sqrt(2*m*v3^2))^2)
  T=c(T,t)
}
T2n=sqrt(2*log(log(n)))*max(T)-2*log(log(n))-2*log(log(log(n)))+log(gamma(2)) 

#####show results
rbind(c(num1,T1n),c(num2,T2n))
