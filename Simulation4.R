library(lars)
library(grpreg)

set.seed(23)
#################################define some functions
TruePositive<-function(estimate, true, error=2){         
  #calculate the number of true positives
  Num=0
  if(length(estimate)>0){
    for(l in 1:length(true)){
      if(min(abs(estimate-true[l]))<=error)
        Num=Num+1
    }
  }
  return(Num)
}

######additional simulation in supplementary material
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
num2=matrix(NA,nr=100,nc=10)         #m=5 Lasso
Num2=matrix(NA,nr=100,nc=10)         #m=5 SCAD
num3=matrix(NA,nr=100,nc=10)         #m=10 Lasso
Num3=matrix(NA,nr=100,nc=10)        #m=10 SCAD
signal=0.5*(x>=0.3*n)-1*(x>=0.4*n)+1*(x>=0.8*n)-0.5*(x>=0.9*n)   
true=c(0.3,0.4,0.8,0.9)*n
wave1=0.1*(cos(2*pi*x/n)+sin(2*pi*x/n)+cos(2*pi*2*x/n)+sin(2*pi*2*x/n)) 
for(i in 1:100){
  e=rnorm(n,mean=0,sd=0.15)  
  y1=signal+wave1+e  
  ###m=5
  ry1=lm(y1~0+Z)$residuals
  model2<-grpreg(X=RX1, y=ry1, group=1:(n-1), penalty="grLasso",family="gaussian", gmax=10)
  for(j in 1:dim(model2$beta)[2]){
    estimate=as.vector((which(model2$beta[,j]!=0))[-1])
    if(length(estimate)<=10)  num2[i,length(estimate)]=TruePositive(estimate, true)
  }
  Model2<- grpreg(X=RX1, y=ry1, group=1:(n-1), penalty="grSCAD",family="gaussian", gmax=10)
  for(j in 1:dim(Model2$beta)[2]){
    estimate=as.vector((which(Model2$beta[,j]!=0))[-1])
    if(length(estimate)<=10)  Num2[i,length(estimate)]=TruePositive(estimate, true)
  }
  ###m=10
  rry1=lm(y1~0+ZZ)$residuals
  model3<-grpreg(X=RRX1, y=rry1, group=1:(n-1), penalty="grLasso",family="gaussian", gmax=10)
  for(j in 1:dim(model3$beta)[2]){
    estimate=as.vector((which(model3$beta[,j]!=0))[-1])
    if(length(estimate)<=10)  num3[i,length(estimate)]=TruePositive(estimate, true)
  }
  Model3<- grpreg(X=RRX1, y=rry1, group=1:(n-1), penalty="grSCAD",family="gaussian", gmax=10)
  for(j in 1:dim(Model3$beta)[2]){
    estimate=as.vector((which(Model3$beta[,j]!=0))[-1])
    if(length(estimate)<=10)  Num3[i,length(estimate)]=TruePositive(estimate, true)
  }
}

####wave2
num5=matrix(NA,nr=100,nc=10)      #m=5 Lasso
Num5=matrix(NA,nr=100,nc=10)      #m=5 SCAD
num6=matrix(NA,nr=100,nc=10)       #m=10 Lasso
Num6=matrix(NA,nr=100,nc=10)      #m=10 SCAD
wave2=wave1+0.1*(cos(2*pi*6*x/n)+sin(2*pi*6*x/n)+cos(2*pi*7*x/n)+sin(2*pi*7*x/n)) 
for(i in 1:100){
  e=rnorm(n,mean=0,sd=0.15)  
  y2=signal+wave2+e  
  ###m=5
  ry2=lm(y2~0+Z)$residuals
  model5<-grpreg(X=RX1, y=ry2, group=1:(n-1), penalty="grLasso",family="gaussian", gmax=10)
  for(j in 1:dim(model5$beta)[2]){
    estimate=as.vector((which(model5$beta[,j]!=0))[-1])
    if(length(estimate)<=10)  num5[i,length(estimate)]=TruePositive(estimate, true)
  }
  Model5<-grpreg(X=RX1, y=ry2, group=1:(n-1), penalty="grSCAD",family="gaussian", gmax=10)
  for(j in 1:dim(Model5$beta)[2]){
    estimate=as.vector((which(Model5$beta[,j]!=0))[-1])
    if(length(estimate)<=10)  Num5[i,length(estimate)]=TruePositive(estimate, true)
  }
  ###m=10
  rry2=lm(y2~0+ZZ)$residuals
  model6<-grpreg(X=RRX1, y=rry2, group=1:(n-1), penalty="grLasso",family="gaussian", gmax=10)
  for(j in 1:dim(model6$beta)[2]){
    estimate=as.vector((which(model6$beta[,j]!=0))[-1])
    if(length(estimate)<=10)  num6[i,length(estimate)]=TruePositive(estimate, true)
  }
  Model6<-grpreg(X=RRX1, y=rry2, group=1:(n-1), penalty="grSCAD",family="gaussian", gmax=10)
  for(j in 1:dim(Model6$beta)[2]){
    estimate=as.vector((which(Model6$beta[,j]!=0))[-1])
    if(length(estimate)<=10)  Num6[i,length(estimate)]=TruePositive(estimate, true)
  }
}

####wave3
num8=matrix(NA,nr=100,nc=10)        #m=5 Lasso
Num8=matrix(NA,nr=100,nc=10)        #m=5 SCAD
num9=matrix(NA,nr=100,nc=10)        #m=10 Lasso
Num9=matrix(NA,nr=100,nc=10)        #m=10 SCAD
for(i in 1:100){
  e=rnorm(n,mean=0,sd=0.15)  
  fi=runif(1, min = 0, max = 2*pi)
  thi=runif(1, min = 0, max = 2*pi)
  wave3=0.1*(sin(2*pi*x/96+fi)+2*sin(2*pi*x/240+thi)) 
  y3=signal+wave3+e  
  ###m=5
  ry3=lm(y3~0+Z)$residuals
  model8<-grpreg(X=RX1, y=ry3, group=1:(n-1), penalty="grLasso",family="gaussian", gmax=10)
  for(j in 1:dim(model8$beta)[2]){
    estimate=as.vector((which(model8$beta[,j]!=0))[-1])
    if(length(estimate)<=10)  num8[i,length(estimate)]=TruePositive(estimate, true)
  }
  Model8<-grpreg(X=RX1, y=ry3, group=1:(n-1), penalty="grSCAD",family="gaussian", gmax=10)
  for(j in 1:dim(Model8$beta)[2]){
    estimate=as.vector((which(Model8$beta[,j]!=0))[-1])
    if(length(estimate)<=10)  Num8[i,length(estimate)]=TruePositive(estimate, true)
  }
  ###m=10
  rry3=lm(y3~0+ZZ)$residuals
  model9<-grpreg(X=RRX1, y=ry3, group=1:(n-1), penalty="grLasso",family="gaussian", gmax=10)
  for(j in 1:dim(model9$beta)[2]){
    estimate=as.vector((which(model9$beta[,j]!=0))[-1])
    if(length(estimate)<=10)  num9[i,length(estimate)]=TruePositive(estimate, true)
  }
  Model9<-grpreg(X=RRX1, y=ry3, group=1:(n-1), penalty="grSCAD",family="gaussian", gmax=10)
  for(j in 1:dim(Model9$beta)[2]){
    estimate=as.vector((which(Model9$beta[,j]!=0))[-1])
    if(length(estimate)<=10)  Num9[i,length(estimate)]=TruePositive(estimate, true)
  }
}

par(mfrow=c(2,3))
plot(x=1:10,y=apply(num2,2,mean,na.rm=T)[1:10],xlab="Number of detected change-points",ylab="Correct selection",ylim=c(0,4),main="Scenario I, m=5", type="b")
lines(x=1:10,y=apply(Num2,2,mean,na.rm=T)[1:10],type="b",pch=4)
legend("bottomright", legend=c("SCAD","Lasso"),pch=c(4,1),bty="n")
plot(x=1:10,y=apply(num5,2,mean,na.rm=T)[1:10],xlab="Number of detected change-points",ylab="Correct selection",ylim=c(0,4),main="Scenario II, m=5", type="b")
lines(x=1:10,y=apply(Num5,2,mean,na.rm=T)[1:10],type="b",pch=4)
legend("bottomright", legend=c("SCAD","Lasso"),pch=c(4,1),bty="n")
plot(x=1:10,y=apply(num8,2,mean,na.rm=T)[1:10],xlab="Number of detected change-points",ylab="Correct selection",ylim=c(0,4),main="Scenario III, m=5", type="b")
lines(x=1:10,y=apply(Num8,2,mean,na.rm=T)[1:10],type="b",pch=4)
legend("bottomright", legend=c("SCAD","Lasso"),pch=c(4,1),bty="n")
plot(x=1:10,y=apply(num3,2,mean,na.rm=T)[1:10],xlab="Number of detected change-points",ylab="Correct selection",ylim=c(0,4),main="Scenario I, m=10", type="b")
lines(x=1:10,y=apply(Num3,2,mean,na.rm=T)[1:10],type="b",pch=4)
legend("bottomright", legend=c("SCAD","Lasso"),pch=c(4,1),bty="n")
plot(x=1:10,y=apply(num6,2,mean,na.rm=T)[1:10],xlab="Number of detected change-points",ylab="Correct selection",ylim=c(0,4),main="Scenario II, m=10", type="b")
lines(x=1:10,y=apply(Num6,2,mean,na.rm=T)[1:10],type="b",pch=4)
legend("bottomright", legend=c("SCAD","Lasso"),pch=c(4,1),bty="n")
plot(x=1:10,y=apply(num9,2,mean,na.rm=T)[1:10],xlab="Number of detected change-points",ylab="Correct selection",ylim=c(0,4),main="Scenario III, m=10", type="b")
lines(x=1:10,y=apply(Num9,2,mean,na.rm=T)[1:10],type="b",pch=4)
legend("bottomright", legend=c("SCAD","Lasso"),pch=c(4,1),bty="n")
