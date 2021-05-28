source("FABlmPredict.R") 

set.seed(1)

n<-100
p<-75
sigma2<-1
alpha<-.1

X0<-matrix(rnorm(n*p),n,p)%*%matrix(rnorm(p*p),p,p)
X0<-sweep(X0,2,apply(X0,2,mean),"-") 
X0<-sweep(X0,2,apply(X0,2,sd),"/") 

taus<-seq(.1,sqrt(10),length=15) 

FRES1<-FRES2<-PRES<-NULL

for(tau in taus){

    fres1<-fres2<-pres<-NULL
    for(s in 1:2500){
      beta<-rnorm(p,0,tau) 
      mu<-X0%*%beta
      e<-rnorm(n)*sqrt(sigma2) 

      PR1<-FABlmPredict(1*mu+e,X0,X0,Psi=diag(p)*sigma2/tau^2,
                         sigmasplit=TRUE,alpha=alpha) 
      PR2<-FABlmPredict(2*mu+e,X0,X0,Psi=diag(p)*sigma2/tau^2,
                         sigmasplit=TRUE,alpha=alpha)  
      fres1<-c(fres1,mean(PR1$FAB[,2]-PR1$FAB[,1]))
      fres2<-c(fres2,mean(PR2$FAB[,2]-PR2$FAB[,1]))
      pres<-c(pres,mean(PR1$PIV[,2]-PR1$PIV[,1]))
      cat(s,tau," ",mean(fres1),mean(fres2)," ",mean(pres),"\n") 
    }
  FRES1<-cbind(FRES1,fres1)   
  FRES2<-cbind(FRES2,fres2) 
  PRES<-cbind(PRES,pres) 
} 


### calculate EPW without simulation
w0<-NULL 
iXX0<-solve(t(X0)%*%X0) 
for(i in 1:n){ w0<-c(w0,1+t(X0[i,])%*%iXX0%*%c(X0[i,])) }
EPR<-2*mean(sqrt(w0))*sqrt(sigma2)*qt(1-alpha/2,n-p)*
    (sqrt(2)*gamma((n-p+1)/2)/gamma((n-p)/2))/sqrt(n-p)

save(taus,FRES1,FRES2,PRES,file="lmSim.rdata") 




### plots 
frisk1mc<-apply(FRES1,2,mean) 
frisk2mc<-apply(FRES2,2,mean) 

frisk1<-smooth.spline(taus,isoreg(taus,frisk1mc)$yf,nknots=5)$y
frisk2<-smooth.spline(taus,isoreg(taus,frisk2mc)$yf,nknots=5)$y

pdf("fig4.pdf",height=3,width=6,family="Times") 
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(taus,frisk1,type="l",lwd=3,ylim=range(c(frisk1,frisk2)),
     xlab=expression(tau),ylab="average Bayes risk",col=gray(.5))
lines(taus,frisk2,lwd=3,col="black")  
#lines(taus,frisk1mc) 
#lines(taus,frisk2mc)
abline(h=EPR,lwd=3,col=gray(.9)) 
legend(1.5,4.3,lwd=3,col=c(gray(.9),gray(.5),"black"),
   legend=c(expression(tau[pi]==infinity),
            expression(tau[pi]==tau),
            expression(tau[pi]==tau/2 )),bty="n") 
dev.off()





