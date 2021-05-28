## functions for z-regions 
source("FABmvnPredict.R") 


## test
test<-FALSE
if(test){
par(mfrow=c(1,3),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
for(p in 1:3){ 
alpha<-.1
ACC<-NULL
for(i in 1:100){ 
  set.seed(i)
 
  ## conditions
  Sh<-crossprod(matrix(rnorm(2*p^2),2*p,p))/sqrt(2*p)
  iSh<-solve(Sh) 
  theta<-Sh%*%rnorm(p) 
  mu<-Sh%*%rnorm(p) 
  k<-rexp(1,1) 
  lambda<-rexp(1,2) 

  ## mc coverage approximation
  acc<-0
  for(j in 1:1000){
    y<-theta + Sh%*%rnorm(p)
    x<-theta + sqrt(k)*Sh%*%rnorm(p)
    acc<-acc+1*(pcheck(y,x,iSh,k,mu,lambda,alpha)<0 )   
  }
  cat(i,acc/j,"\n") 
  ACC<-c(ACC,acc/j) 
} 

mci<-1-alpha+c(-1,0,1)*qnorm(1-.5/i)*sqrt(alpha*(1-alpha)/j) 
hist(ACC,xlim=range(c(ACC,mci)),main="") 
abline(v=mci) 
}
}

## global variables
k<-1 ; alpha<-.1 ;lambdas<-c(.1,1,10,100) 
cols<-gray( (1:(1+length(lambdas))-1)/(1+length(lambdas))) 


## p=1 widths as a function of x
p<-1 ; sigma<-1 ; mu<-0 
lambda<-1 
PYXF<-PYX0<-NULL
xs<-seq(-4,4,length=100)
for(x in xs){
  PYXF<-rbind(PYXF,predZR1(x,sigma,k,mu,lambda,alpha))
  PYX0<-rbind(PYX0,x+sqrt(k+1)*c(-1,1)*qnorm(alpha/2))
 } 

V1<-NULL
for(lambda in lambdas){ 
for(r in seq(0,4,length=100)){ 
  x<-r
  pyx<-predZR1(x,sigma,k,mu,lambda,alpha) 
  V1<-rbind(V1,c(lambda,r, diff(pyx))) 
}}

pdf("prZp1.pdf",height=4,width=8,family="Times")
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(range(xs),range(c(PYXF,PYX0)),type="n",xlab=expression(italic(x)),
     ylab="interval endpoints") 
abline(h=0,lty=2,col="gray") 
abline(v=0,lty=2,col="gray") 
legend(-4.25,6.75,lwd=c(3,3),col=cols[c(2,5)],bty="n",
legend=c(expression(lambda==1),expression(lambda==infinity)))
matplot(xs,PYXF,type="l",lwd=2,lty=1,col=cols[2],add=TRUE) 
matplot(xs,PYX0,type="l",lwd=2,lty=1,col=cols[5],add=TRUE) 


v1inf<-sqrt(k+1)*2*sqrt(qchisq(1-alpha,1))
plot(range(V1[,2]),range(V1[,3]),type="n",
     xlab=expression("|"*italic(x)*"|"),ylab="interval width") 
for(lambda in lambdas){ lines(V1[V1[,1]==lambda,2:3], lwd=3,
    col=cols[match(lambda,lambdas)]) }
abline(h=v1inf ,col=cols[5],lwd=3) 
legend(-.15,7.6,lwd=c(3,3),col=c(cols),bty="n",
legend=c(
  expression(lambda==1/10),
  expression(lambda==1),
  expression(lambda==10),
  expression(lambda==100),
  expression(lambda==infinity)))

dev.off()

# p=2 areas as a function of x
p<-2 ; iSh<-diag(p) ; mu<-rep(0,p) 
lambda<-1
PXYF0<-list()
for(r in 0:4){ 
  x<-r*c(1,1)/sqrt(2) 
  PXYF0[[length(PXYF0)+1]]<-predZR2(x,iSh,k,mu,lambda,alpha)  
} 

V2<-NULL
for(lambda in lambdas){ 
for(r in seq(0,4,length=100)){ 
  x<-r*c(1/sqrt(2),1/sqrt(2))  
  pyx<-predZR2(x,iSh,k,mu,lambda,alpha) 
  V2<-rbind(V2,c(lambda,r,pyx$vol)) 
}}

pdf("prZp2.pdf",height=4,width=8,family="Times")
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(c(-4,7),c(-3.75,7),
     type="n",xlab=expression(italic(x)[1]),ylab=expression(italic(x)[2]))     
abline(v=0,lty=2,col="gray")    
abline(h=0,lty=2,col="gray")   
for(i in 1:length(PXYF0)){ lines(PXYF0[[i]]$vertInf,col=cols[5]) }
for(i in 1:length(PXYF0)){ lines(PXYF0[[i]]$vert,lwd=2,col=cols[2]) }
legend(-4.5,7.3,lwd=c(3,3),col=cols[c(2,5)],bty="n",
legend=c(expression(lambda==1),expression(lambda==infinity)))

v2inf<-pi*prod(1/eigen(iSh)$val)*qchisq(1-alpha,2)*(k+1)
plot(range(V2[,2]),range(V2[,3]),type="n",
     xlab=expression("||"*italic(x)*"||"),ylab="region area") 
for(lambda in lambdas){ lines(V2[V2[,1]==lambda,2:3], lwd=3,
    col=cols[match(lambda,lambdas)]) }
abline(h=v2inf ,col=cols[5],lwd=3) 
legend(-.15,50,lwd=c(3,3),col=cols,bty="n",
legend=c(
  expression(lambda==1/10),
  expression(lambda==1),
  expression(lambda==10),
  expression(lambda==100),
  expression(lambda==infinity)))
dev.off()


