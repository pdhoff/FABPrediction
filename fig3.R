## functions for z-regions 
source("FABmvnPredict.R") 

## global variables
k<-1 ; alpha<-.1 ;lambdas<-c(.1,1,10,100) ; nmc<-1000

## p=1
try(load("mvnSim.rdata"),silent=TRUE)
if(!exists("V1")){ 
p<-1 ; sigma<-1 ; mu<-0 

V1<-NULL
for(lambda in lambdas){ 
for(theta in seq(0,4,length=25)){ 
  aw<-0  
  x<-rnorm(nmc,theta,sigma) 
  for(s in 1:nmc){ 
    aw<-aw+diff(predZR1(x[s],sigma,k,mu,lambda,alpha))
  } 
  V1<-rbind(V1,c(lambda,theta,aw/nmc)) 
cat(lambda,theta,"\n") 
}}
}


V1M<-t(tapply(V1[,3],list(lambda=V1[,1],theta=V1[,2]),mean) ) 
rtheta1<-as.numeric(rownames(V1M))
V1S<-apply(V1M,2,function(y){ lm(y~poly(rtheta1,3))$fitted } ) 


## p=2
if(!exists("V2")){
p<-2 ; iSh<-diag(p) ; mu<-rep(0,p) 

V2<-NULL
for(lambda in lambdas){ 
for(r in seq(0,4,length=25)){  
  theta<-r*c(1,1)/sqrt(2)
  aw<-0  
  for(s in 1:nmc){ 
    x<-theta+rnorm(2)
    aw<-aw+predZR2(x,iSh,k,mu,lambda,alpha)$vol
  } 
  V2<-rbind(V2,c(lambda,theta,aw/nmc)) 
cat(lambda,theta,"\n") 
}}
}

r<-V2[,2]*sqrt(2) 
V2M<-t(tapply(V2[,4],list(lambda=V2[,1],ntheta=r),mean) )
rtheta2<-as.numeric(rownames(V2M)) 
V2S<-apply(V2M,2,function(y){ lm(y~poly(rtheta2,3))$fitted } ) 


pdf("fig3.pdf",height=4,width=8,family="Times") 

par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))

cols<-gray( (1:(1+length(lambdas))-1)/(1+length(lambdas))) 


v1inf<-sqrt(k+1)*2*sqrt(qchisq(1-alpha,1))
plot(range(rtheta1),range(V1S),type="n",
     xlab=expression("|"*italic(x)*"|"),ylab="expected interval width") 
for(j in 1:ncol(V1S)){ lines(rtheta1,V1S[,j],lwd=3,col=cols[j] ) } 
abline(h=v1inf ,col="gray",lwd=3) 
legend(-.15,7.6,lwd=c(3,3),col=c(cols),bty="n",
legend=c(
  expression(lambda==1/10),
  expression(lambda==1),
  expression(lambda==10),
  expression(lambda==100),
  expression(lambda==infinity)))

v2inf<-pi*qchisq(1-alpha,2)*(k+1)
plot(range(rtheta2),range(V2S),type="n",
     xlab=expression("||"*italic(x)*"||"),ylab="expected region area") 
for(j in 1:ncol(V2S)){ lines(rtheta2,V2S[,j],lwd=3,col=cols[j] ) } 
abline(h=v2inf ,col="gray",lwd=3) 
legend(-.15,7.6,lwd=c(3,3),col=c(cols),bty="n",
legend=c(
  expression(lambda==1/10),
  expression(lambda==1),
  expression(lambda==10),
  expression(lambda==100),
  expression(lambda==infinity)))

dev.off()

save(V1,V2,file="mvnSim.rdata") 

