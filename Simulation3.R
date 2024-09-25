library(lars)
library(grpreg)

###Simulation3
set.seed(1005)
#######matrix X and Z
n=500
d=4
x=1:n
x1=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1))
  x1[i,(i+1):n]=0
X=kronecker(x1,diag(1, d))                          #contain intercept
group=NULL
for (i in 1:n)
  group<-c(group,rep(i-1,d))                        #contain intercept
zz=matrix(1,nrow=n,ncol=21)     
for (i in 1:10){
  zz[,(2*i)]=cos(2*pi*i*x/n)
  zz[,(2*i+1)]=sin(2*pi*i*x/n)
}
z=zz[,1:11]
Z=kronecker(z[,2:11],diag(1, d))
ZZ=kronecker(zz[,2:21],diag(1, d))
XX=cbind(X,Z)
XXX=cbind(X,ZZ)
Group=c(group, rep(0,(10*4)))
GGroup=c(group, rep(0,(20*4)))
RX1=matrix(1,nrow=n,ncol=n-1)
for(i in 1:(n-1)){
  RX1[,i]=lm(x1[,(i+1)]~0+z)$residuals
}
RRX1=matrix(1,nrow=n,ncol=n-1)
for(i in 1:(n-1)){
  RRX1[,i]=lm(x1[,(i+1)]~0+zz)$residuals
}

####change points estimation
loc1=matrix(0,nrow=200,ncol=4)           #m=0, d=1
loc2=matrix(0,nrow=200,ncol=4)           #m=0, d=4
loc3=matrix(0,nrow=200,ncol=4)           #m=5, d=1
loc4=matrix(0,nrow=200,ncol=4)           #m=5, d=4
signal=1*(x>=0.2*n)-2*(x>=0.4*n)+2*(x>=0.7*n)-1*(x>=0.9*n)
wave=cos(2*pi*x/n)+sin(2*pi*x/n)+cos(2*pi*2*x/n)+sin(2*pi*2*x/n)
for (i in 1:200){
  e1=rnorm(n,mean=0,sd=0.3)     
  e2=rnorm(n,mean=0,sd=0.3) 
  e3=rnorm(n,mean=0,sd=0.3)     
  e4=rnorm(n,mean=0,sd=0.3)   
  y1=signal+e1+sample(c(-0.3,0.3),size=1,replace=TRUE)*wave 
  y2=signal+e2+sample(c(-0.3,0.3),size=1,replace=TRUE)*wave
  y3=signal+e3+sample(c(-0.3,0.3),size=1,replace=TRUE)*wave  
  y4=signal+e4+sample(c(-0.3,0.3),size=1,replace=TRUE)*wave  
  y0=NULL
  for (j in 1:n)
    y0<-c(y0,c(y1[j],y2[j],y3[j],y4[j]))
  #m=0
  model1=lars(x1[,2:n], y =y1, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=5) 
  s1=unlist(model1$actions)  
  loc1[i,1]=s1[which.min(abs(s1-0.2*n))]
  loc1[i,2]=s1[which.min(abs(s1-0.4*n))]
  loc1[i,3]=s1[which.min(abs(s1-0.7*n))]
  loc1[i,4]=s1[which.min(abs(s1-0.9*n))]
  model2=grpreg(X=X[,2:length(y0)], y=y0, group=group[2:length(y0)], penalty="grLasso",family="gaussian", gmax=4)      #group lasso
  a=which(model2$beta[,length(model2$lambda)]!=0)
  s2=(a[4*(1:(length(a)/4))]/4-1)[-1]
  loc2[i,1]=s2[which.min(abs(s2-0.2*n))]
  loc2[i,2]=s2[which.min(abs(s2-0.4*n))]
  loc2[i,3]=s2[which.min(abs(s2-0.7*n))]
  loc2[i,4]=s2[which.min(abs(s2-0.9*n))]
  #m=5
  ry1=lm(y1~0+z)$residuals
  model3=lars(RX1,y=ry1,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=5)
  s3=unlist(model3$actions) 
  loc3[i,1]=s3[which.min(abs(s3-0.2*n))]
  loc3[i,2]=s3[which.min(abs(s3-0.4*n))]
  loc3[i,3]=s3[which.min(abs(s3-0.7*n))]
  loc3[i,4]=s3[which.min(abs(s3-0.9*n))]
  model4=grpreg(X=XX[,-1], y=y0, group=Group[-1], penalty="grLasso",family="gaussian", gmax=4)      #group lasso
  a=which(model4$beta[1:(4*n),length(model4$lambda)]!=0)
  s4=(a[4*(1:(length(a)/4))]/4-1)[-1]
  loc4[i,1]=s4[which.min(abs(s4-0.2*n))]
  loc4[i,2]=s4[which.min(abs(s4-0.4*n))]
  loc4[i,3]=s4[which.min(abs(s4-0.7*n))]
  loc4[i,4]=s4[which.min(abs(s4-0.9*n))]
}

#show results
CoverPr1=matrix(0,nrow=4,ncol=4)
CoverPr1[1,]=c(sum((loc1[,1]>0.2*n-3)*(loc1[,1]<0.2*n+3))/200,sum((loc1[,2]>0.4*n-3)*(loc1[,2]<0.4*n+3))/200,sum((loc1[,3]>0.7*n-3)*(loc1[,3]<0.7*n+3))/200,sum((loc1[,4]>0.9*n-3)*(loc1[,4]<0.9*n+3))/200)
CoverPr1[2,]=c(sum((loc2[,1]>0.2*n-3)*(loc2[,1]<0.2*n+3))/200,sum((loc2[,2]>0.4*n-3)*(loc2[,2]<0.4*n+3))/200,sum((loc2[,3]>0.7*n-3)*(loc2[,3]<0.7*n+3))/200,sum((loc2[,4]>0.9*n-3)*(loc2[,4]<0.9*n+3))/200)
CoverPr1[3,]=c(sum((loc3[,1]>0.2*n-3)*(loc3[,1]<0.2*n+3))/200,sum((loc3[,2]>0.4*n-3)*(loc3[,2]<0.4*n+3))/200,sum((loc3[,3]>0.7*n-3)*(loc3[,3]<0.7*n+3))/200,sum((loc3[,4]>0.9*n-3)*(loc3[,4]<0.9*n+3))/200)
CoverPr1[4,]=c(sum((loc4[,1]>0.2*n-3)*(loc4[,1]<0.2*n+3))/200,sum((loc4[,2]>0.4*n-3)*(loc4[,2]<0.4*n+3))/200,sum((loc4[,3]>0.7*n-3)*(loc4[,3]<0.7*n+3))/200,sum((loc4[,4]>0.9*n-3)*(loc4[,4]<0.9*n+3))/200)

####change points estimation
loc1=matrix(0,nrow=200,ncol=4)           #m=0, d=1
loc2=matrix(0,nrow=200,ncol=4)           #m=0, d=4
loc3=matrix(0,nrow=200,ncol=4)           #m=5, d=1
loc4=matrix(0,nrow=200,ncol=4)           #m=5, d=4
signal=1*(x>=0.2*n)-2*(x>=0.4*n)+2*(x>=0.7*n)-1*(x>=0.9*n)
for (i in 1:200){
  e1=rnorm(n,mean=0,sd=0.3)     
  e2=rnorm(n,mean=0,sd=0.3) 
  e3=rnorm(n,mean=0,sd=0.3)     
  e4=rnorm(n,mean=0,sd=0.3)   
  wave1=sin(2*pi*x/96+runif(1, min = 0, max = 2*pi))+2*sin(2*pi*x/240+runif(1, min = 0, max = 2*pi)) 
  wave2=sin(2*pi*x/96+runif(1, min = 0, max = 2*pi))+2*sin(2*pi*x/240+runif(1, min = 0, max = 2*pi)) 
  wave3=sin(2*pi*x/96+runif(1, min = 0, max = 2*pi))+2*sin(2*pi*x/240+runif(1, min = 0, max = 2*pi)) 
  wave4=sin(2*pi*x/96+runif(1, min = 0, max = 2*pi))+2*sin(2*pi*x/240+runif(1, min = 0, max = 2*pi)) 
  y1=signal+e1+0.2*wave1
  y2=signal+e2+0.2*wave1
  y3=signal+e3+0.2*wave1
  y4=signal+e4+0.2*wave1
  y0=NULL
  for (j in 1:n)
    y0<-c(y0,c(y1[j],y2[j],y3[j],y4[j]))
  #m=0
  model1=lars(x1[,2:n], y =y1, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=5) 
  s1=unlist(model1$actions)  
  loc1[i,1]=s1[which.min(abs(s1-0.2*n))]
  loc1[i,2]=s1[which.min(abs(s1-0.4*n))]
  loc1[i,3]=s1[which.min(abs(s1-0.7*n))]
  loc1[i,4]=s1[which.min(abs(s1-0.9*n))]
  model2=grpreg(X=X[,2:length(y0)], y=y0, group=group[2:length(y0)], penalty="grLasso",family="gaussian", gmax=4)      #group lasso
  a=which(model2$beta[,length(model2$lambda)]!=0)
  s2=(a[4*(1:(length(a)/4))]/4-1)[-1]
  loc2[i,1]=s2[which.min(abs(s2-0.2*n))]
  loc2[i,2]=s2[which.min(abs(s2-0.4*n))]
  loc2[i,3]=s2[which.min(abs(s2-0.7*n))]
  loc2[i,4]=s2[which.min(abs(s2-0.9*n))]
  #m=5
  ry1=lm(y1~0+z)$residuals
  model3=lars(RX1,y=ry1,type="lasso", normalize = FALSE, intercept = FALSE, trace = FALSE, max.steps=5)
  s3=unlist(model3$actions) 
  loc3[i,1]=s3[which.min(abs(s3-0.2*n))]
  loc3[i,2]=s3[which.min(abs(s3-0.4*n))]
  loc3[i,3]=s3[which.min(abs(s3-0.7*n))]
  loc3[i,4]=s3[which.min(abs(s3-0.9*n))]
  model4=grpreg(X=XX[,-1], y=y0, group=Group[-1], penalty="grLasso",family="gaussian", gmax=4)      #group lasso
  a=which(model4$beta[1:(4*n),length(model4$lambda)]!=0)
  s4=(a[4*(1:(length(a)/4))]/4-1)[-1]
  loc4[i,1]=s4[which.min(abs(s4-0.2*n))]
  loc4[i,2]=s4[which.min(abs(s4-0.4*n))]
  loc4[i,3]=s4[which.min(abs(s4-0.7*n))]
  loc4[i,4]=s4[which.min(abs(s4-0.9*n))]
}

#show results
CoverPr2=matrix(0,nrow=4,ncol=4)
CoverPr2[1,]=c(sum((loc1[,1]>0.2*n-3)*(loc1[,1]<0.2*n+3))/200,sum((loc1[,2]>0.4*n-3)*(loc1[,2]<0.4*n+3))/200,sum((loc1[,3]>0.7*n-3)*(loc1[,3]<0.7*n+3))/200,sum((loc1[,4]>0.9*n-3)*(loc1[,4]<0.9*n+3))/200)
CoverPr2[2,]=c(sum((loc2[,1]>0.2*n-3)*(loc2[,1]<0.2*n+3))/200,sum((loc2[,2]>0.4*n-3)*(loc2[,2]<0.4*n+3))/200,sum((loc2[,3]>0.7*n-3)*(loc2[,3]<0.7*n+3))/200,sum((loc2[,4]>0.9*n-3)*(loc2[,4]<0.9*n+3))/200)
CoverPr2[3,]=c(sum((loc3[,1]>0.2*n-3)*(loc3[,1]<0.2*n+3))/200,sum((loc3[,2]>0.4*n-3)*(loc3[,2]<0.4*n+3))/200,sum((loc3[,3]>0.7*n-3)*(loc3[,3]<0.7*n+3))/200,sum((loc3[,4]>0.9*n-3)*(loc3[,4]<0.9*n+3))/200)
CoverPr2[4,]=c(sum((loc4[,1]>0.2*n-3)*(loc4[,1]<0.2*n+3))/200,sum((loc4[,2]>0.4*n-3)*(loc4[,2]<0.4*n+3))/200,sum((loc4[,3]>0.7*n-3)*(loc4[,3]<0.7*n+3))/200,sum((loc4[,4]>0.9*n-3)*(loc4[,4]<0.9*n+3))/200)

######null distribution n=500，d=1,4
n=500   
True=1:100000
sTrue=1:100000
for(i in 1:100000){
  e1=rnorm(n,mean=0,sd=1)     
  e2=rnorm(n,mean=0,sd=1)  
  e3=rnorm(n,mean=0,sd=1)     
  e4=rnorm(n,mean=0,sd=1)      
  T=NULL
  for(k in 1:n){
    t=sqrt((sum((e1[1:k])^2-1)/sqrt(2*k))^2+(sum((e2[1:k])^2-1)/sqrt(2*k))^2+(sum((e3[1:k])^2-1)/sqrt(2*k))^2+(sum((e4[1:k])^2-1)/sqrt(2*k))^2)
    T=c(T,t)
  }
  True[i]=sqrt(2*log(log(n)))*max(T)-2*log(log(n))-2*log(log(log(n)))+log(gamma(2))         #test statistic for d=4
}
for(i in 1:100000){
  e=rnorm(n,mean=0,sd=1)     
  sT=NULL
  for(m in 1:n){
    st=(2*m)^(-0.5)*sum((e[1:m])^2-1)
    sT=c(sT,st)
  }
  sTrue[i]=sqrt(2*log(log(n)))*max(sT)-2*log(log(n))-0.5*log(log(log(n)))+0.5*log(4*pi)         #test statistic for d=1
}
q=quantile(True,prob=0.95)                              #d=4                  
sq=quantile(sTrue,prob=0.95)                            #d=1
#q=4.581187
#sq=3.869114

######data generating for power
power0=matrix(0,nrow=200,ncol=6)       #for d=4 
power1=matrix(0,nrow=200,ncol=6)
power2=matrix(0,nrow=200,ncol=6)
spower0=matrix(0,nrow=200,ncol=6)        #for d=1
spower1=matrix(0,nrow=200,ncol=6)
spower2=matrix(0,nrow=200,ncol=6)
prob=matrix(0,nrow=200,ncol=6)
sprob=matrix(0,nrow=200,ncol=6)
beta=c(0,0.02,0.04,0.06,0.08,0.10)
signal=1*(x>=0.2*n)-2*(x>=0.4*n)+2*(x>=0.7*n)-1*(x>=0.9*n)
for(i in 1:200){        
  e1=rnorm(n,mean=0,sd=0.3)     
  e2=rnorm(n,mean=0,sd=0.3) 
  e3=rnorm(n,mean=0,sd=0.3)     
  e4=rnorm(n,mean=0,sd=0.3)   
  wave1=sin(2*pi*x/96+runif(1, min = 0, max = 2*pi))+2*sin(2*pi*x/240+runif(1, min = 0, max = 2*pi)) 
  wave2=sin(2*pi*x/96+runif(1, min = 0, max = 2*pi))+2*sin(2*pi*x/240+runif(1, min = 0, max = 2*pi)) 
  wave3=sin(2*pi*x/96+runif(1, min = 0, max = 2*pi))+2*sin(2*pi*x/240+runif(1, min = 0, max = 2*pi)) 
  wave4=sin(2*pi*x/96+runif(1, min = 0, max = 2*pi))+2*sin(2*pi*x/240+runif(1, min = 0, max = 2*pi)) 
  y01=signal+e1              #under Null
  y02=signal+e2
  y03=signal+e3
  y04=signal+e4
  ########power 
  for(k in 1:6){
    y11=y01+beta[k]*wave1      
    y12=y02+beta[k]*wave2 
    y13=y03+beta[k]*wave3 
    y14=y04+beta[k]*wave4 
    y1=NULL
    for (j in 1:n)
      y1<-c(y1,c(y11[j],y12[j],y13[j],y14[j]))
    r1=lm(y11~x1[,c(0.2*n,0.4*n,0.7*n,0.9*n)])$residuals           #oracle, d=4
    r2=lm(y12~x1[,c(0.2*n,0.4*n,0.7*n,0.9*n)])$residuals
    r3=lm(y13~x1[,c(0.2*n,0.4*n,0.7*n,0.9*n)])$residuals
    r4=lm(y14~x1[,c(0.2*n,0.4*n,0.7*n,0.9*n)])$residuals
    #Fourier transform
    fft1=NULL
    for(j in 1:(n/2)){
      a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r1)
      b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r1)
      fft1=c(fft1, a, b)
    }
    v1=var(fft1[floor(n/4):n])        #sigma^2
    fft2=NULL
    for(j in 1:(n/2)){
      a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r2)
      b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r2)
      fft2=c(fft2, a, b)
    }
    v2=var(fft2[floor(n/4):n])        #sigma^2
    fft3=NULL
    for(j in 1:(n/2)){
      a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r3)
      b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r3)
      fft3=c(fft3, a, b)
    }
    v3=var(fft3[floor(n/4):n])        #sigma^2
    fft4=NULL
    for(j in 1:(n/2)){
      a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r4)
      b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r4)
      fft4=c(fft4, a, b)
    }
    v4=var(fft4[floor(n/4):n])        #sigma^2
    T=NULL
    for(m in 1:n){
      t=sqrt((sum((fft1[1:m])^2-v1)/sqrt(2*m*v1^2))^2+(sum((fft2[1:m])^2-v2)/sqrt(2*m*v2^2))^2+(sum((fft3[1:m])^2-v3)/sqrt(2*m*v3^2))^2+(sum((fft4[1:m])^2-v4)/sqrt(2*m*v4^2))^2)
      T=c(T,t)
    }
    power0[i,k]=(sqrt(2*log(log(n)))*max(T)-2*log(log(n))-2*log(log(log(n)))+log(gamma(2))>q)   
    sT=NULL                #oracle, d=1
    for(m in 1:n){
      t=(2*m*v1^2)^(-0.5)*sum((fft1[1:m])^2-v1)
      sT=c(sT,t)
    }     
    spower0[i,k]=(sqrt(2*log(log(n)))*max(sT)-2*log(log(n))-0.5*log(log(log(n)))+0.5*log(4*pi)>sq) 
  }
  #lasso with s=10, d=4
  for(k in 1:6){
    y11=y01+beta[k]*wave1        #Under alternative
    y12=y02+beta[k]*wave2 
    y13=y03+beta[k]*wave3 
    y14=y04+beta[k]*wave4 
    y1=NULL
    for (j in 1:n)
      y1<-c(y1,c(y11[j],y12[j],y13[j],y14[j]))
    model1=grpreg(X=X[1:length(y1),2:length(y1)], y=y1, group=group[2:length(y1)], penalty="grLasso",family="gaussian", gmax=9)      #group lasso
    a=which(model1$beta[,length(model1$lambda)]!=0)
    b=(a[4*(1:(length(a)/4))]/4-1)[-1]
    sub=which((b[2:length(b)]-b[1:(length(b)-1)])==1)
    if (length(sub)==0){
      r1=lm(y11~x1[,(b+1)])$residuals
      r2=lm(y12~x1[,(b+1)])$residuals
      r3=lm(y13~x1[,(b+1)])$residuals
      r4=lm(y14~x1[,(b+1)])$residuals
    }
    if (length(sub)>0){
      r1=lm(y11~x1[,(b[-sub]+1)])$residuals
      r2=lm(y12~x1[,(b[-sub]+1)])$residuals
      r3=lm(y13~x1[,(b[-sub]+1)])$residuals
      r4=lm(y14~x1[,(b[-sub]+1)])$residuals
    }
    #Fourier transform
    fft1=NULL
    for(j in 1:(n/2)){
      a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r1)
      b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r1)
      fft1=c(fft1, a, b)
    }
    v1=var(fft1[floor(n/4):n])        #sigma^2
    fft2=NULL
    for(j in 1:(n/2)){
      a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r2)
      b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r2)
      fft2=c(fft2, a, b)
    }
    v2=var(fft2[floor(n/4):n])        #sigma^2
    fft3=NULL
    for(j in 1:(n/2)){
      a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r3)
      b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r3)
      fft3=c(fft3, a, b)
    }
    v3=var(fft3[floor(n/4):n])        #sigma^2
    fft4=NULL
    for(j in 1:(n/2)){
      a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r4)
      b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r4)
      fft4=c(fft4, a, b)
    }
    v4=var(fft4[floor(n/4):n])        #sigma^2
    T=NULL
    for(m in 1:n){
      t=sqrt((sum((fft1[1:m])^2-v1)/sqrt(2*m*v1^2))^2+(sum((fft2[1:m])^2-v2)/sqrt(2*m*v2^2))^2+(sum((fft3[1:m])^2-v3)/sqrt(2*m*v3^2))^2+(sum((fft4[1:m])^2-v4)/sqrt(2*m*v4^2))^2)
      T=c(T,t)
    }
    power1[i,k]=(sqrt(2*log(log(n)))*max(T)-2*log(log(n))-2*log(log(log(n)))+log(gamma(2))>q)   
    smodel1=lars(x1[1:n,2:n], y =y11, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=10)      #lasso with s=10,d=1
    sub=which((which(smodel1$entry>0)[2:10]-which(smodel1$entry>0)[1:9])==1)
    if (length(sub)==0)
      r1=lm(y11~x1[,(which(smodel1$entry>0)+1)])$residuals
    if (length(sub)>0)
      r1=lm(y11~x1[,(which(smodel1$entry>0)[-sub]+1)])$residuals
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
    spower1[i,k]=(sqrt(2*log(log(n)))*max(T1)-2*log(log(n))-0.5*log(log(log(n)))+0.5*log(4*pi)>sq)
  }
  #lasso with s=5, d=4
  for(k in 1:6){
    y11=y01+beta[k]*wave1        #Under alternative
    y12=y02+beta[k]*wave2 
    y13=y03+beta[k]*wave3 
    y14=y04+beta[k]*wave4 
    y1=NULL
    for (j in 1:n)
      y1<-c(y1,c(y11[j],y12[j],y13[j],y14[j]))
    model2=grpreg(X=X[1:length(y1),2:length(y1)], y=y1, group=group[2:length(y1)], penalty="grLasso",family="gaussian", gmax=4)      #group lasso
    a=which(model2$beta[,length(model2$lambda)]!=0)
    b=(a[4*(1:(length(a)/4))]/4-1)[-1]
    sub=which((b[2:length(b)]-b[1:(length(b)-1)])==1)
    if (length(sub)==0){
      r1=lm(y11~x1[,(b+1)])$residuals
      r2=lm(y12~x1[,(b+1)])$residuals
      r3=lm(y13~x1[,(b+1)])$residuals
      r4=lm(y14~x1[,(b+1)])$residuals
    }
    if (length(sub)>0){
      r1=lm(y11~x1[,(b[-sub]+1)])$residuals
      r2=lm(y12~x1[,(b[-sub]+1)])$residuals
      r3=lm(y13~x1[,(b[-sub]+1)])$residuals
      r4=lm(y14~x1[,(b[-sub]+1)])$residuals
    }
    prob[i,k]=sum((b<453)*(b>447))
    #Fourier transform
    fft1=NULL
    for(j in 1:(n/2)){
      a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r1)
      b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r1)
      fft1=c(fft1, a, b)
    }
    v1=var(fft1[floor(n/4):n])        #sigma^2
    fft2=NULL
    for(j in 1:(n/2)){
      a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r2)
      b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r2)
      fft2=c(fft2, a, b)
    }
    v2=var(fft2[floor(n/4):n])        #sigma^2
    fft3=NULL
    for(j in 1:(n/2)){
      a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r3)
      b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r3)
      fft3=c(fft3, a, b)
    }
    v3=var(fft3[floor(n/4):n])        #sigma^2
    fft4=NULL
    for(j in 1:(n/2)){
      a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r4)
      b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r4)
      fft4=c(fft4, a, b)
    }
    v4=var(fft4[floor(n/4):n])        #sigma^2
    T=NULL
    for(m in 1:n){
      t=sqrt((sum((fft1[1:m])^2-v1)/sqrt(2*m*v1^2))^2+(sum((fft2[1:m])^2-v2)/sqrt(2*m*v2^2))^2+(sum((fft3[1:m])^2-v3)/sqrt(2*m*v3^2))^2+(sum((fft4[1:m])^2-v4)/sqrt(2*m*v4^2))^2)
      T=c(T,t)
    }
    power2[i,k]=(sqrt(2*log(log(n)))*max(T)-2*log(log(n))-2*log(log(log(n)))+log(gamma(2))>q)   
    smodel2=lars(x1[1:n,2:n], y =y11, type ="lasso", normalize = FALSE, intercept = TRUE, trace = FALSE, max.steps=5)      #lasso with s=5, d=1
    sub=which((which(smodel2$entry>0)[2:5]-which(smodel2$entry>0)[1:4])==1)
    if (length(sub)==0)
      r1=lm(y11~x1[,(which(smodel2$entry>0)+1)])$residuals
    if (length(sub)>0)
      r1=lm(y11~x1[,(which(smodel2$entry>0)[-sub]+1)])$residuals
    sprob[i,k]=sum((which(smodel2$entry>0)<453)*(which(smodel2$entry>0)>447))
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
    spower2[i,k]=(sqrt(2*log(log(n)))*max(T1)-2*log(log(n))-0.5*log(log(log(n)))+0.5*log(4*pi)>sq)
  }
}
sum(prob[,1]>0)/200
sum(sprob[,1]>0)/200

####show the plots for powers
par(mfrow=c(1,3))
plot(x=beta,y=c(mean(power0[,1]), mean(power0[,2]), mean(power0[,3]), mean(power0[,4]), mean(power0[,5]), mean(power0[,6])), xlab="theta",ylab="power",ylim=c(0,1),main="Case I", type="b")
lines(x=beta,y=c(mean(spower0[,1]), mean(spower0[,2]), mean(spower0[,3]), mean(spower0[,4]), mean(spower0[,5]), mean(spower0[,6])), xlab="theta",ylab="power",main="Case I", type="b",pch=4)
legend("bottomright", legend=c("d=1","d=4"),pch=c(4,1),bty="n")
plot(x=beta,y=c(mean(power1[,1]), mean(power1[,2]), mean(power1[,3]), mean(power1[,4]), mean(power1[,5]), mean(power1[,6])), xlab="theta",ylab="power",ylim=c(0,1),main="Case II", type="b")
lines(x=beta,y=c(mean(spower1[,1]), mean(spower1[,2]), mean(spower1[,3]), mean(spower1[,4]), mean(spower1[,5]), mean(spower1[,6])), xlab="theta",ylab="power",main="Case II", type="b",pch=4)
legend("bottomright", legend=c("d=1","d=4"),pch=c(4,1),bty="n")
plot(x=beta,y=c(mean(power2[,1]), mean(power2[,2]), mean(power2[,3]), mean(power2[,4]), mean(power2[,5]), mean(power2[,6])), xlab="theta",ylab="power",ylim=c(0,1),main="Case III", type="b")
lines(x=beta,y=c(mean(spower2[,1]), mean(spower2[,2]), mean(spower2[,3]), mean(spower2[,4]), mean(spower2[,5]), mean(spower2[,6])), xlab="theta",ylab="power",main="Case III", type="b",pch=4)
legend("bottomright", legend=c("d=1","d=4"),pch=c(4,1),bty="n")

####show the table
cbind(CoverPr1,CoverPr2)
