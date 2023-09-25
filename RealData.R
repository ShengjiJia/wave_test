library(lars)
library(penalized)
library(grpreg)
library(imputeTS)
library(cumSeg)
library(DNAcopy)

####real data 
CGHdata <- read.csv("C:/Users/Acer/Desktop/My Document/Documents/Research/Projects/change points(wave test)/CGHdataset.csv",header=T)
index=c(11, 19, 45, 56)
data=CGHdata[2:2301,(1+3*index)]               
n=nrow(data)
d=ncol(data)
x=1:n
for (i in 1:d){
  data[,i]=na_ma(data[,i],k=5,weighting="linear")        #imputation
} 

####matrix X_1 and Z
x1=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1))
  x1[i,(i+1):n]=0                
zz=matrix(1,nrow=n,ncol=21)                     
for (i in 1:10){
  zz[,(2*i)]=cos(2*pi*i*x/n)
  zz[,(2*i+1)]=sin(2*pi*i*x/n)             
}
X=kronecker(x1,diag(1, d))                          #with intercept term
group=NULL
for (i in 1:n)
  group<-c(group,rep(i-1,d))                        #with intercept term

#####1st sequence  
y1=data[,1]
num1=1:4
#####CBS
set.seed(2)
CBS=DNAcopy::segment(CNA(y1, rep(1,n), 1:n))                
num1[1]=length(CBS$output[,4])-1
######cumSeg
num1[2]=jumpoints(y1,k=66,output="2")$n.psi
######LB
model1<-lars(x1[,2:n], y =y1, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=66)   
n3=which.min(log(model1$RSS/n) + log(n)*model1$df*2*log(log(n))/n)                 
model11<-lars(x1[,2:n], y =y1, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=n3)
num1[3]=n3-sum(diff(sort(abs(unlist(model11$actions))))<=2)              #delete adjacent estimators
######adaptive Neyman test 
sub=which(diff(sort(abs(unlist(model1$actions))))<2)
r1=lm(y1~x1[,which(model1$entry>0)[-sub]+1])$residuals 
#Adaptive Neyman test
fft1=NULL
for(j in 1:n/2){
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
T1n=sqrt(2*log(log(n)))*max(T1)-2*log(log(n))-0.5*log(log(log(n)))+0.5*log(4*pi) 
#######select m1
gBIC=1:10
for (j in 1:10){
  y11=lm(y1~0+zz[,1:(1+2*j)])$residuals
  RX1=matrix(1,nrow=n,ncol=n-1)
  for(i in 1:(n-1)){
    RX1[,i]=lm(x1[,(i+1)]~0+zz[,1:(1+2*j)])$residuals
  }
  model3=lars(RX1,y=y11,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=66)
  gBIC[j]=min(log(model3$RSS/n) + log(n)*(model3$df+2*j)*2*log(log(n))/n)
}
m1=which.min(gBIC)
#######partial penalized with selected m1
yy1=lm(y1~0+zz[,1:(1+2*m1)])$residuals
RRX1=matrix(1,nrow=n,ncol=n-1)
for(i in 1:(n-1)){
  RRX1[,i]=lm(x1[,(i+1)]~0+zz[,1:(1+2*m1)])$residuals
}
model2<-lars(RRX1,y=yy1,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=66)
n4=which.min(log(model2$RSS/n) + log(n)*model2$df*2*log(log(n))/n)                
model22=lars(RRX1,y=yy1,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=n4)
sub=which(diff(sort(abs(unlist(model22$actions))))<=2)
loc1=sort(abs(unlist(model22$actions)))[-sub]
num1[4]=length(loc1)

#####2nd sequence  
y2=data[,2]
num2=1:4
#####CBS
CBS=DNAcopy::segment(CNA(y2, rep(1,n), 1:n))                
num2[1]=length(CBS$output[,4])-1
######cumSeg
num2[2]=jumpoints(y2,k=66,output="2")$n.psi
######LB
model1<-lars(x1[,2:n], y =y2, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=66)   
n3=which.min(log(model1$RSS/n) + log(n)*model1$df*2*log(log(n))/n)
model11<-lars(x1[,2:n], y =y2, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=n3)
num2[3]=n3-sum(diff(sort(abs(unlist(model11$actions))))<=2)
######adaptive Neyman test 
sub=which(diff(sort(abs(unlist(model1$actions))))<2)
r1=lm(y2~x1[,which(model1$entry>0)[-sub]+1])$residuals 
#Adaptive Neyman test
fft1=NULL
for(j in 1:n/2){
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
T2n=sqrt(2*log(log(n)))*max(T1)-2*log(log(n))-0.5*log(log(log(n)))+0.5*log(4*pi) 
#######select m2
gBIC=1:10
for (j in 1:10){
  y22=lm(y2~0+zz[,1:(1+2*j)])$residuals
  RX1=matrix(1,nrow=n,ncol=n-1)
  for(i in 1:(n-1)){
    RX1[,i]=lm(x1[,(i+1)]~0+zz[,1:(1+2*j)])$residuals
  }
  model3=lars(RX1,y=y22,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=66)
  gBIC[j]=min(log(model3$RSS/n) + log(n)*(model3$df+2*j)*2*log(log(n))/n)
}
m2=which.min(gBIC)
#######partial penalized with selected m2
yy2=lm(y2~0+zz[,1:(1+2*m2)])$residuals
RRX1=matrix(1,nrow=n,ncol=n-1)
for(i in 1:(n-1)){
  RRX1[,i]=lm(x1[,(i+1)]~0+zz[,1:(1+2*m2)])$residuals
}
model2<-lars(RRX1,y=yy2,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=66)
n4=which.min(log(model2$RSS/n) + log(n)*model2$df*2*log(log(n))/n)
model22<-lars(RRX1,y=yy2,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=n4)
sub=which(diff(sort(abs(unlist(model22$actions))))<=2)
loc2=sort(abs(unlist(model22$actions)))[-sub]
num2[4]=length(loc2)

#####3nd sequence  
y3=data[,3]
num3=1:4
#####CBS
CBS=DNAcopy::segment(CNA(y3, rep(1,n), 1:n))                
num3[1]=length(CBS$output[,4])-1
######cumSeg
num3[2]=jumpoints(y3,k=66,output="2")$n.psi
######LB
model1<-lars(x1[,2:n], y =y3, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=66)   
n3=which.min(log(model1$RSS/n) + log(n)*model1$df*2*log(log(n))/n)
model11<-lars(x1[,2:n], y =y3, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=n3)   
num3[3]=n3-sum(diff(sort(abs(unlist(model11$actions))))<=2)
######adaptive Neyman test 
sub=which(diff(sort(abs(unlist(model1$actions))))<2)
r1=lm(y3~x1[,which(model1$entry>0)[-sub]+1])$residuals 
#Adaptive Neyman test
fft1=NULL
for(j in 1:n/2){
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
T3n=sqrt(2*log(log(n)))*max(T1)-2*log(log(n))-0.5*log(log(log(n)))+0.5*log(4*pi) 
#######select m3
gBIC=1:10
for (j in 1:10){
  y33=lm(y3~0+zz[,1:(1+2*j)])$residuals
  RX1=matrix(1,nrow=n,ncol=n-1)
  for(i in 1:(n-1)){
    RX1[,i]=lm(x1[,(i+1)]~0+zz[,1:(1+2*j)])$residuals
  }
  model3=lars(RX1,y=y33,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=66)
  gBIC[j]=min(log(model3$RSS/n) + log(n)*(model3$df+2*j)*2*log(log(n))/n)
}
m3=which.min(gBIC)
#######partial penalized with selected m3
yy3=lm(y3~0+zz[,1:(1+2*m3)])$residuals
RRX1=matrix(1,nrow=n,ncol=n-1)
for(i in 1:(n-1)){
  RRX1[,i]=lm(x1[,(i+1)]~0+zz[,1:(1+2*m3)])$residuals
}
model2<-lars(RRX1,y=yy3,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=66)
n4=which.min(log(model2$RSS/n) + log(n)*model2$df*2*log(log(n))/n)
model22<-lars(RRX1,y=yy3,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=n4)
sub=which(diff(sort(abs(unlist(model22$actions))))<=2)
loc3=sort(abs(unlist(model22$actions)))[-sub]
num3[4]=length(loc3)

#####4nd sequence  
y4=data[,4]
num4=1:4
#####CBS
CBS=DNAcopy::segment(CNA(y4, rep(1,n), 1:n))                
num4[1]=length(CBS$output[,4])-1
######cumSeg
num4[2]=jumpoints(y4,k=66,output="2")$n.psi
######LB
model1<-lars(x1[,2:n], y =y4, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=66)   
n3=which.min(log(model1$RSS/n) + log(n)*model1$df*2*log(log(n))/n)
model11<-lars(x1[,2:n], y =y4, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=n3)   
num4[3]=n3-sum(diff(sort(abs(unlist(model11$actions))))<=2)
######adaptive Neyman test 
sub=which(diff(sort(abs(unlist(model1$actions))))<2)
r1=lm(y4~x1[,which(model1$entry>0)[-sub]+1])$residuals 
#Adaptive Neyman test
fft1=NULL
for(j in 1:n/2){
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
T4n=sqrt(2*log(log(n)))*max(T1)-2*log(log(n))-0.5*log(log(log(n)))+0.5*log(4*pi) 
#######select m4
gBIC=1:10
for (j in 1:10){
  y44=lm(y4~0+zz[,1:(1+2*j)])$residuals
  RX1=matrix(1,nrow=n,ncol=n-1)
  for(i in 1:(n-1)){
    RX1[,i]=lm(x1[,(i+1)]~0+zz[,1:(1+2*j)])$residuals
  }
  model3=lars(RX1,y=y44,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=66)
  gBIC[j]=min(log(model3$RSS/n) + log(n)*(model3$df+2*j)*2*log(log(n))/n)
}
m4=which.min(gBIC)
#######partial penalized with selected m4
yy4=lm(y4~0+zz[,1:(1+2*m4)])$residuals
RRX1=matrix(1,nrow=n,ncol=n-1)
for(i in 1:(n-1)){
  RRX1[,i]=lm(x1[,(i+1)]~0+zz[,1:(1+2*m4)])$residuals
}
model2<-lars(RRX1,y=yy4,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=66)
n4=which.min(log(model2$RSS/n) + log(n)*model2$df*2*log(log(n))/n)
model22<-lars(RRX1,y=yy4,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=n4)
sub=which(diff(sort(abs(unlist(model22$actions))))<=2)
loc4=sort(abs(unlist(model22$actions)))[-sub]
num4[4]=length(loc4)


#####show Table 
rbind(c(num1,T1n),c(num2,T2n),c(num3,T3n),c(num4,T4n))      


#####show Figure
par(mfrow=c(1,2))
plot(y3[1:500], x=1:500, xlab="locations", ylab="Log 2 ratio", main="X506",pch=20,col=8)
lines(lm(y3~0+zz[,1:(1+2*m3)]+x1[,(loc3+1)])$fitted.values[1:500],x=1:500,col=2,lwd=2)
plot(y4[1:500], x=1:500, xlab="locations", ylab="Log 2 ratio", main="X2259-1",pch=20,col=8)
lines(lm(y4~0+zz[,1:(1+2*m4)]+x1[,(loc4+1)])$fitted.values[1:500],x=1:500,col=2,lwd=2)
m=max(c(m1,m2,m3,m4))
ZZ=kronecker(zz[,2:(1+2*m)],diag(1, d))
XX=cbind(X,ZZ)
Group=c(group, rep(0,(2*m*d)))
y0=as.vector(t(data))
model4=grpreg(X=XX[,-1], y=y0, group=Group[-1], penalty="grLasso",family="gaussian", gmax=66)     #group lasso
nonzero=as.vector(which(model4$beta[1:(n*d),length(model4$lambda)]!=0))
sub=which(diff(nonzero[4*(2:(length(nonzero)/4))]/4)<=15)
loc=(nonzero[4*(2:(length(nonzero)/4))]/4)[-sub]
par(mfrow=c(4,1))
plot(y1, x=1:n, xlab="locations", ylab="Log 2 ratio", main="X1410",pch=20,col=8)
lines(lm(y1~0+zz[,1:(1+2*m1)]+x1[,(loc1+1)])$fitted.values,x=1:n,col=2)
abline(v=loc, lty=2)
plot(y2, x=1:n, xlab="locations", ylab="Log 2 ratio", ylim=c(-1,3), main="X1533-1",pch=20,col=8)
lines(lm(y2~0+zz[,1:(1+2*m2)]+x1[,(loc2+1)])$fitted.values,x=1:n,col=2)
abline(v=loc, lty=2)
plot(y3, x=1:n, xlab="locations", ylab="Log 2 ratio", main="X506",pch=20,col=8)
lines(lm(y3~0+zz[,1:(1+2*m3)]+x1[,(loc3+1)])$fitted.values,x=1:n,col=2)
abline(v=loc, lty=2)
plot(y4, x=1:n, xlab="locations", ylab="Log 2 ratio", main="X2259-1",pch=20,col=8)
lines(lm(y4~0+zz[,1:(1+2*m4)]+x1[,(loc4+1)])$fitted.values,x=1:n,col=2)
abline(v=loc, lty=2)

