## critical points for optimal statistic
qX<-function(x,y,iSh,k,mu,lambda,alpha=.05){ 
  # x observed value
  # y candidate predicted value 
  # iSh inverse square root of covariance matrix 
  # k scale for variance of x 
  # mu prior mean of theta
  # lambda scale for prior variance of theta
  # alpha error rate
  z<-(x+k*y)/(1+k) 
  vxz<-k^2/(k+1) 
  b<-iSh%*%(mu-z)*sqrt(k+1)/(lambda*(1+1/k)+1)
  qchisq(1-alpha,length(x),ncp=apply(b^2,2,sum))

} 


## check if y is in Ax
pcheck<-function(y,x,iSh,k,mu,lambda,alpha=.05){ 
  vinf<-k+1 
  vlam<-1+1/(1/lambda+1/k) 
  theta<-(x/k +mu/lambda)/(1/k+1/lambda)  
  sum( (iSh%*%(y-theta))^2/vlam )*vinf/vlam - 
  qX(x,y,iSh,k,mu,lambda,alpha) 
}


## p=1 prediction region 
predZR1<-function(x,sigma,k,mu,lambda,alpha=.05){  
  # x observed value 
  # sigma standard deviation of Y
  # k k*sigma^2 is variance of x 
  # mu prior mean of theta 
  # lambda lambda*sigma^2 is prior variance of theta  
  # alpha error rate 

  ## setup
  vinf<-k+1 
  vlam<-1+1/(1/lambda+1/k) 
  theta<-(x/k+mu/lambda)/(1/k+1/lambda)  

  ## endpoint solution
  f<-function(y){ 
     (vinf/vlam)*(y-theta)^2/(sigma^2*vlam) -
      qX(x,y,1/sigma,k,mu,lambda,alpha) }

  ## usual interval 
  pinf<-x+c(-1,1)*sigma*sqrt(1+k)*qnorm(1-alpha/2) 

  ## bounds for uniroot
  lb<-pinf[1]-sigma*sqrt(1+k)*qnorm(1-alpha/2)
  while(f(lb)*f(x)>0){ lb<-lb-.1*sigma } 
    
  ub<-pinf[1]+sigma*sqrt(1+k)*qnorm(1-alpha/2)
  while(f(ub)*f(x)>0){ ub<-ub+.1*sigma } 

  ## fab interval
  c(uniroot(f,c(lb,x))$root,uniroot(f,c(x,ub))$root) 

}


## p=2 prediction region
predZR2<-function(x,iSh,k,mu,lambda,alpha=.05,nch=100){  

  # x observed value 
  # iSh inverse square root of covariance of Y
  # k k*Sigma is variance of x 
  # mu prior mean of theta 
  # lambda lambda*Sigma is prior variance of theta  
  # alpha error rate 

 
  ## setup
  vinf<-k+1 
  vlam<-1+1/(1/lambda+1/k) 
  theta<-(x/k+mu/lambda)/(1/k+1/lambda)  

  ## make usual region
  w<-seq(0,2*pi,length=nch) 
  cs<-cbind(cos(w),sin(w))
  r<-sqrt( qchisq(1-alpha,2)*(k+1)/apply((cs%*%iSh)^2,1,sum) )
  Y0<-sweep(cs*r,2,x,"+") 
  vol0<-pi*prod(1/eigen(iSh)$val)*qchisq(1-alpha,2)*(k+1)

  ## deform usual region
  y0<-(x+theta)/2
  YL<-Y0
  for(i in 1:nrow(Y0)){
    g<-function(v){ pcheck(y0+v*(Y0[i,]-y0),x,iSh,k,mu,lambda,alpha)  }  
    # bounds for uniroot 
    ub<-1 
    while( g(0)*g(ub)>0 ){ ub<-ub+.1 }
    YL[i,]<- y0+(Y0[i,]-y0)*uniroot(g,c(0,ub))$root 
  }  

  ## compute convex hull and area
  chull<-geometry::convhulln(YL,output.options=TRUE)  

  ## results 
  list(vert=YL,vol=chull$vol,
       vertInf=Y0,volInf=vol0,
       volInfcheck=geometry::convhulln(Y0,output.options=TRUE)$vol) 
}

