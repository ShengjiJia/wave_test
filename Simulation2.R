library(lars)
library(DNAcopy)     #download from https://bioconductor.org/packages/release/bioc/html/DNAcopy.html

set.seed(23)
######Simulation2
######Hausdorff distance
hausdorff<-function(x, y){
  if (length(y)*length(x)==0){  
    return(NA)
  }   
  else{
    A=rep(0, length(y))
    for(i in 1:length(y)){
      A[i]=min(abs(y[i]-x))  
    }
    B=rep(0, length(x))
    for(i in 1:length(x)){
      B[i]=min(abs(x[i]-y))
    }
    return(max(max(A),max(B)))
  }
}

######matrix Z (m=5,10)
n=500
x=1:n
ZZ=matrix(1,nrow=n,ncol=21)
for (i in 1:10){
  ZZ[,(2*i)]=cos(2*pi*i*x/n)
  ZZ[,(2*i+1)]=sin(2*pi*i*x/n)
}
Z=ZZ[,1:11]
######regress X on Z (m=5,10)
X1=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1))
  X1[i,(i+1):n]=0
X1=X1[,2:n]
RX1=matrix(1,nrow=n,ncol=n-1)
for(i in 1:(n-1)){
  RX1[,i]=lm(X1[,i]~0+Z)$residuals
}
RRX1=matrix(1,nrow=n,ncol=n-1)
for(i in 1:(n-1)){
  RRX1[,i]=lm(X1[,i]~0+ZZ)$residuals
}
######data generating
####wave1
num1=rep(0,10)         #m=0
num2=rep(0,10)         #m=5
num3=rep(0,10)         #m=10
Num1=1:500              #CBS correct
haus1=matrix(0,nrow=500,ncol=10)         #m=0
haus2=matrix(0,nrow=500,ncol=10)         #m=5
haus3=matrix(0,nrow=500,ncol=10)         #m=10
Haus1=1:500
signal=0.5*(x>=0.3*n)-1*(x>=0.4*n)+1*(x>=0.8*n)-0.5*(x>=0.9*n)   
true=c(0.3,0.4,0.8,0.9)*n
wave1=0.1*(cos(2*pi*x/n)+sin(2*pi*x/n)+cos(2*pi*2*x/n)+sin(2*pi*2*x/n)) 
for(i in 1:500){
  e=rnorm(n,mean=0,sd=0.15)  
  y1=signal+wave1+e  
  ###m=0
  model1<-lars(X1, y =y1, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=10)       
  s1=unlist(model1$actions) 
  for(j in 1:10){
    num1[j]=num1[j]+(sum((s1[1:j]>147)*(s1[1:j]<153))>0)+(sum((s1[1:j]>197)*(s1[1:j]<203))>0)+(sum((s1[1:j]>397)*(s1[1:j]<403))>0)+(sum((s1[1:j]>447)*(s1[1:j]<453))>0)
    haus1[i,j]=hausdorff(s1[1:j], true)
  }
  ###m=5
  ry1=lm(y1~0+Z)$residuals
  model2<-lars(RX1,y=ry1,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=10)
  s2=unlist(model2$actions)
  for(j in 1:10){
    num2[j]=num2[j]+(sum((s2[1:j]>147)*(s2[1:j]<153))>0)+(sum((s2[1:j]>197)*(s2[1:j]<203))>0)+(sum((s2[1:j]>397)*(s2[1:j]<403))>0)+(sum((s2[1:j]>447)*(s2[1:j]<453))>0)
    haus2[i,j]=hausdorff(s2[1:j], true)
  }
  ###m=10
  rry1=lm(y1~0+ZZ)$residuals
  model3<-lars(RRX1,y=rry1,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=10)
  s3=unlist(model3$actions)
  for(j in 1:10){
    num3[j]=num3[j]+(sum((s3[1:j]>147)*(s3[1:j]<153))>0)+(sum((s3[1:j]>197)*(s3[1:j]<203))>0)+(sum((s3[1:j]>397)*(s3[1:j]<403))>0)+(sum((s3[1:j]>447)*(s3[1:j]<453))>0)
    haus3[i,j]=hausdorff(s3[1:j], true)
  }
  ####CBS
  CBS=DNAcopy::segment(CNA(y1, rep(1,n), 1:n))
  estimate1=CBS$output[-length(CBS$output[,4]),4]
  Num1[i]=(sum((estimate1>147)*(estimate1<153))>0)+(sum((estimate1>197)*(estimate1<203))>0)+(sum((estimate1>397)*(estimate1<403))>0)+(sum((estimate1>447)*(estimate1<453))>0)
  Haus1[i]=hausdorff(estimate1, true)
}

####wave2
num4=rep(0,10)
num5=rep(0,10)
num6=rep(0,10)
Num4=1:500
haus4=matrix(0,nrow=500,ncol=10)
haus5=matrix(0,nrow=500,ncol=10)
haus6=matrix(0,nrow=500,ncol=10)
Haus4=1:500
wave2=wave1+0.1*(cos(2*pi*6*x/n)+sin(2*pi*6*x/n)+cos(2*pi*7*x/n)+sin(2*pi*7*x/n)) 
for(i in 1:500){
  e=rnorm(n,mean=0,sd=0.15)  
  y2=signal+wave2+e  
  model4<-lars(X1, y =y2, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=10)       
  s4=unlist(model4$actions)
  for(j in 1:10){
    num4[j]=num4[j]+(sum((s4[1:j]>147)*(s4[1:j]<153))>0)+(sum((s4[1:j]>197)*(s4[1:j]<203))>0)+(sum((s4[1:j]>397)*(s4[1:j]<403))>0)+(sum((s4[1:j]>447)*(s4[1:j]<453))>0)
    haus4[i,j]=hausdorff(s4[1:j], true)
  }
  ry2=lm(y2~0+Z)$residuals
  model5<-lars(RX1,y=ry2,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=10)
  s5=unlist(model5$actions)
  for(j in 1:10){
    num5[j]=num5[j]+(sum((s5[1:j]>147)*(s5[1:j]<153))>0)+(sum((s5[1:j]>197)*(s5[1:j]<203))>0)+(sum((s5[1:j]>397)*(s5[1:j]<403))>0)+(sum((s5[1:j]>447)*(s5[1:j]<453))>0)
    haus5[i,j]=hausdorff(s5[1:j], true)
  }
  rry2=lm(y2~0+ZZ)$residuals
  model6<-lars(RRX1,y=rry2,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=10)
  s6=unlist(model6$actions)
  for(j in 1:10){
    num6[j]=num6[j]+(sum((s6[1:j]>147)*(s6[1:j]<153))>0)+(sum((s6[1:j]>197)*(s6[1:j]<203))>0)+(sum((s6[1:j]>397)*(s6[1:j]<403))>0)+(sum((s6[1:j]>447)*(s6[1:j]<453))>0)
    haus6[i,j]=hausdorff(s6[1:j], true)
  }
  ####CBS
  CBS=DNAcopy::segment(CNA(y2, rep(1,n), 1:n))
  estimate1=CBS$output[-length(CBS$output[,4]),4]
  Num4[i]=(sum((estimate1>147)*(estimate1<153))>0)+(sum((estimate1>197)*(estimate1<203))>0)+(sum((estimate1>397)*(estimate1<403))>0)+(sum((estimate1>447)*(estimate1<453))>0)
  Haus4[i]=hausdorff(estimate1, true)
}

####wave3
num7=rep(0,10)
num8=rep(0,10)
num9=rep(0,10)
Num7=1:500
haus7=matrix(0,nrow=500,ncol=10)
haus8=matrix(0,nrow=500,ncol=10)
haus9=matrix(0,nrow=500,ncol=10)
Haus7=1:500
for(i in 1:500){
  e=rnorm(n,mean=0,sd=0.15)  
  fi=runif(1, min = 0, max = 2*pi)
  thi=runif(1, min = 0, max = 2*pi)
  wave3=0.1*(sin(2*pi*x/96+fi)+2*sin(2*pi*x/240+thi)) 
  y3=signal+wave3+e  
  model7<-lars(X1, y =y3, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=10)       
  s7=unlist(model7$actions) 
  for(j in 1:10){
    num7[j]=num7[j]+(sum((s7[1:j]>147)*(s7[1:j]<153))>0)+(sum((s7[1:j]>197)*(s7[1:j]<203))>0)+(sum((s7[1:j]>397)*(s7[1:j]<403))>0)+(sum((s7[1:j]>447)*(s7[1:j]<453))>0)
    haus7[i,j]=hausdorff(s7[1:j], true)
  }
  ry3=lm(y3~0+Z)$residuals
  model8<-lars(RX1,y=ry3,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=10)
  s8=unlist(model8$actions) 
  for(j in 1:10){
    num8[j]=num8[j]+(sum((s8[1:j]>147)*(s8[1:j]<153))>0)+(sum((s8[1:j]>197)*(s8[1:j]<203))>0)+(sum((s8[1:j]>397)*(s8[1:j]<403))>0)+(sum((s8[1:j]>447)*(s8[1:j]<453))>0)
    haus8[i,j]=hausdorff(s8[1:j], true)
  }
  rry3=lm(y3~0+ZZ)$residuals
  model9<-lars(RRX1,y=rry3,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=10)
  s9=unlist(model9$actions)
  for(j in 1:10){
    num9[j]=num9[j]+(sum((s9[1:j]>147)*(s9[1:j]<153))>0)+(sum((s9[1:j]>197)*(s9[1:j]<203))>0)+(sum((s9[1:j]>397)*(s9[1:j]<403))>0)+(sum((s9[1:j]>447)*(s9[1:j]<453))>0)
    haus9[i,j]=hausdorff(s9[1:j], true)
  }
  ####CBS
  CBS=DNAcopy::segment(CNA(y3, rep(1,n), 1:n))
  estimate1=CBS$output[-length(CBS$output[,4]),4]
  Num7[i]=(sum((estimate1>147)*(estimate1<153))>0)+(sum((estimate1>197)*(estimate1<203))>0)+(sum((estimate1>397)*(estimate1<403))>0)+(sum((estimate1>447)*(estimate1<453))>0)
  Haus7[i]=hausdorff(estimate1, true)
}

par(mfrow=c(2,3))
plot(x=1:10,y=num1/500,pch=2,col=3, type="b", ylim=c(0,4), xlab = "Number of selected change points",ylab ="Correct selection",main="Scenario I")
lines(x=1:10,y=num2/500, col=2, type="b")
lines(x=1:10,y=num3/500,pch=4, col=4, type="b")
legend("bottomright", legend=c("LASSO(m=0)","m=5","m=10"),pch=c(2,1,4), col=c(3,2,4),bty="n")
plot(x=1:10,y=num4/500,pch=2,col=3, type="b", ylim=c(0,4), xlab = "Number of selected change points",ylab ="Correct selection",main="Scenario II")
lines(x=1:10,y=num5/500, col=2, type="b")
lines(x=1:10,y=num6/500,pch=4, col=4, type="b")
legend("bottomright", legend=c("LASSO(m=0)","m=5","m=10"),pch=c(2,1,4), col=c(3,2,4),bty="n")
plot(x=1:10,y=num7/500,pch=2,col=3, type="b", ylim=c(0,4), xlab = "Number of selected change points",ylab ="Correct selection",main="Scenario III")
lines(x=1:10,y=num8/500,col=2, type="b")
lines(x=1:10,y=num9/500,pch=4, col=4, type="b")
legend("bottomright", legend=c("LASSO(m=0)","m=5","m=10"),pch=c(2,1,4), col=c(3,2,4),bty="n")

plot(x=1:10,y=apply(haus1,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,400), xlab = "Number of selected change points",ylab ="Hausdorff distance",main="Scenario I")
lines(x=1:10,y=apply(haus2,2,mean,na.rm=T), col=2, type="b")
lines(x=1:10,y=apply(haus3,2,mean,na.rm=T),pch=4, col=4, type="b")
lines(x=1:10, y=rep(mean(Haus1),10), pch=3, type="b")
legend("topright", legend=c("LASSO(m=0)","m=5","m=10","CBS"), pch=c(2,1,4,3), col=c(3,2,4,1), bty="n")
plot(x=1:10,y=apply(haus4,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,400), xlab = "Number of selected change points",ylab ="Hausdorff distance",main="Scenario II")
lines(x=1:10,y=apply(haus5,2,mean,na.rm=T), col=2, type="b")
lines(x=1:10,y=apply(haus6,2,mean,na.rm=T),pch=4, col=4, type="b")
lines(x=1:10, y=rep(mean(Haus4),10), pch=3, type="b")
legend("topright", legend=c("LASSO(m=0)","m=5","m=10", "CBS"), pch=c(2,1,4,3), col=c(3,2,4,1), bty="n")
plot(x=1:10,y=apply(haus7,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,400), xlab = "Number of selected change points",ylab ="Hausdorff distance",main="Scenario III")
lines(x=1:10,y=apply(haus8,2,mean,na.rm=T), col=2, type="b")
lines(x=1:10,y=apply(haus9,2,mean,na.rm=T),pch=4, col=4, type="b")
lines(x=1:10, y=rep(mean(Haus7),10), pch=3, type="b")
legend("topright", legend=c("LASSO(m=0)","m=5","m=10", "CBS"), pch=c(2,1,4,3), col=c(3,2,4,1), bty="n")

