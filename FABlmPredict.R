sigmasplit<-function(y,X,dfsplit=3/4){

  # split RSS to obtain two independent estimates of sigma^2
  G<-MASS::Null(X)
  e<-t(G)%*%y   
  df0<-round((nrow(X)-ncol(X))*dfsplit)
  e0<-e[1:df0] 
  e1<-e[-(1:df0)] 
  sigma0<-sqrt(mean(e0^2))
  sigma1<-sqrt(mean(e1^2))
  list(sigma0=sigma0,df0=df0,sigma1=sigma1,df1=length(e1)) 
}



FABlmPredict<-function(y0,X0,X1,Psi=diag(1/sigmadelta^2,nrow=ncol(X0)),
   sigma=sqrt( sum(lm(y0~-1+X0)$res^2)/(nrow(X0)-ncol(X0))), 
   df=nrow(X0)-ncol(X0),
   sigmadelta=1,
   sigmasplit=NA,
   alpha=.05)
{ 
  ## y0 vector of outcomes in original sample 
  ## X0 matrix of covariates in original sample
  ## X1 matrix of covariates at which to make prediction intervals
  ## Psi penalty parameter - inverse precision scaled by s2
  
  ## sigma estimate of error std dev
  ## df df for sigma 
  ## sigmadelta independent estimate of error std dev
  ## alpha error rate/non-coverage level
  ## sigmasplit split ss for sigma estimates  

  sigmaHat<-sqrt( sum(lm(y0~-1+X0)$res^2)/(nrow(X0)-ncol(X0)) )
  
  if(!is.na(sigmasplit)){
    if(sigmasplit==TRUE){ dfsplit<-3/4 } 
    if(sigmasplit!=TRUE){ dfsplit<-sigmasplit }
    sigmas<-sigmasplit(y0,X0,dfsplit)  
    sigma<-sigmas$sigma0
    df<-sigmas$df0
    sigmadelta<-sigmas$sigma1 
    if(sum(abs(Psi))==0){ sigma<-sigmaHat ; df<-nrow(X0)-ncol(X0) } 
  }

  qTdelta<-function(delta,alpha,df){
    fq<-function(q) { pt(q-delta,df) - pt(-q-delta,df) - (1-alpha) } 
    lb<- abs(delta) + qnorm(1-alpha) - 1e-5
    ub<-lb 
    while(fq(ub)<0){ ub<-ub*1.1 } 
    uniroot( fq,  interval= c(lb,ub) )$root 
  }

  S0<-solve(crossprod(X0))
  Sg<-solve(crossprod(X0)+Psi)
  tX0y0<-t(X0)%*%y0
  beta<-S0%*%tX0y0
  X1beta<-X1%*%beta

  FAB<-PIV<-NULL
  for(i in 1:nrow(X1)){
    x1<-X1[i,]
    w0<-1+c(t(x1)%*%S0%*%x1)
    wg<-1+c(t(x1)%*%Sg%*%x1)
    deltaprez<-t(x1)%*%(S0/w0 - Sg/wg)*sqrt(w0)/sigmadelta
    x1beta<-X1beta[i,] 

    piv<-x1beta+sigmaHat*sqrt(w0)*qt(1-alpha/2,nrow(X0)-ncol(X0))*c(-1,1)
    PIV<-rbind(PIV,piv)

    ## upper endpoint 
    fy<-function(y){
      z<-tX0y0+x1*y 
      deltaz<-c(deltaprez%*%z)
      qz<-qTdelta(deltaz,alpha,df)
      (y - x1beta)/(sigma*sqrt(w0)) + deltaz - qz
    }

    lu<- sort(piv[2]*c(.9,1.1))
    while( fy(lu[1])*fy(lu[2])>0 ){ lu<-piv[2]+(lu-piv[2])*1.1 }
    yu<-uniroot(fy,interval=lu)$root

    ## lower endpoint 
    fy<-function(y){
      z<-tX0y0+x1*y 
      deltaz<-c(deltaprez%*%z)
      qz<-qTdelta(deltaz,alpha,df)
      (y - x1beta)/(sigma*sqrt(w0)) + deltaz + qz
    }
    
    lu<-sort(piv[1]*c(.9,1.1))
    while( fy(lu[1])*fy(lu[2])>0 ){ lu<-piv[1]+(lu-piv[1])*1.1 }
    yl<-uniroot(fy,interval=lu)$root 
    FAB<-rbind(FAB,c(yl,yu) ) 
  }

  list(FAB=FAB,PIV=PIV,yfit=c(X1beta))
}


