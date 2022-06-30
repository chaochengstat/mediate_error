####################################################################################
#                                                                                  #
#                                R Code                                            #
#                                                                                  #
#       Mediation analysis in the presence of exposure measurement error under     #
#                   main study/validation study designs                            #
#                            from Multiple Studies                                 #
#                                                                                  #
#                   For any question please contact Chao Cheng                     #   
#                         Email: c.cheng@yale.edu                                  #
#                                                                                  #
#                                                                                  #
####################################################################################


################################################################################################
#
# There are two sections, the 1st section is the functions, the 2nd is an illustrative example
#
#################################################################################################

#############################################################################################
#   
# SECTION 1: FUNCTIONS
#
#############################################################################################

rm(list=ls())
library("spatstat")
library("nleqslv")


# data generation function
# n total sample size
# n2 validation study sample size
# rho correlation between X* and X
# cormx correlation between M and X
gen_data=function(n=200,n2=50,rho=0.6,cormx=0.4,TE=log(2),MP=0.5,PY=0.5,outcome="binary") {
  # n=5000
  # n2=200
  # rho=0.8
  # cormx=0.6
  # TE = log(2)
  # MP=0.5
  # PY=0.1
  # outcome="binary"
  
  n1=n-n2
  NIE = TE*MP
  
  XMW.mat=MASS::mvrnorm(n=n,mu=c(0,0,0),Sigma=matrix(c(1,cormx,0.2,cormx,1,0.2,0.2,0.2,1),ncol=3,nrow=3,byrow=T))
  X=XMW.mat[,1]
  M=XMW.mat[,2]
  W=XMW.mat[,3]
  Xs=X+rnorm(n,mean=0,sd=sqrt(1/rho^2 - 1))
  gamma1=c(solve(matrix(c(1,0.2,0.2,1),ncol=2,nrow=2,byrow=T)) %*% c(cormx,0.2))[1]
  beta1=TE*(1-MP)
  beta2= (TE*MP)/gamma1
  
  if (outcome=="binary") {
    beta0.f=function(beta0) {
      m.f=beta0+beta1*0+beta2*0+log(1.2)*0;
      v.f=t(c(beta1,beta2,log(1.2))) %*% matrix(c(1,cormx,0.2,cormx,1,0.2,0.2,0.2,1),ncol=3,nrow=3,byrow=T) %*% c(beta1,beta2,log(1.2))
      z=gauss.hermite(f=function(x) {exp(x)}, mu = m.f, sd = sqrt(v.f[1,1]),  order = 40)
      PY-z
    }
    beta0.res=nleqslv::nleqslv(c(-1),beta0.f)
    beta0=beta0.res$x
    PYMXW=exp(beta0+beta1*X+beta2*M+log(1.2)*W)#/(1+exp(beta0+beta1*X+beta2*M))
    Y=as.numeric(runif(n)<PYMXW)
  } else {
    beta0=2
    Y=beta0+beta1*X+beta2*M+log(1.2)*W+rnorm(n)
  }
  I = rep(0,n); I[(n1+1):n]=1
  X[which(I==0)]=NA
  data.frame(X=X,M=M,Y=Y,Xs=Xs,W=W,I=I)
}

my.ghq = function(f, mu = 0, sd = 1, ..., order = 5) {
  Hn <- HermiteCoefs(order)
  Hn1 <- HermiteCoefs(order - 1)
  x <- sort(Re(polyroot(Hn)))
  Hn1x <- matrix(Hn1, nrow = 1) %*% t(outer(x, 0:(order - 1),"^"))
  w <- 2^(order - 1) * factorial(order) * sqrt(pi)/(order * Hn1x)^2
  ww <- w/sqrt(pi)
  xx <- mu + sqrt(2)*sd %*% t(x)
  fall=f(xx[,1],...)
  for (i in (2:order)) {
    fall=cbind(fall,f(xx[,i],...))
  }
  c(fall %*% t(ww))
}

# Functions for the full EEE
get_eee1=function(data,Xname="X",Xsname="Xs",Yname="Y",Mname="M",Wname="W",Zname="W",Iname="I") {
  #Xname="X";Xsname="Xs";Yname="Y";Mname="M";Wname="W";Zname="W";Iname="I"
  mydata=data
  
  gradU=function(name="Ua",a,sa,b,sb,c,sc,d,sd,model=1) {
    mygrad=function (f, x0, heps = 1e-6, ...) {
      fun <- match.fun(f)
      f <- function(x) fun(x, ...)
      p =length(f(x0))
      n <- length(x0)
      hh <- rep(0, n)
      gr <- matrix(0,nrow=n,ncol=p)
      for (i in 1:n) {
        hh[i] <- heps
        gr[i,] <- (f(x0 + hh) - f(x0 - hh))/(2 * heps)
        hh[i] <- 0
      }
      return(gr)
    }
    if (model==1) {
      x=c(a,c,sc);loc1=1:length(a);
      loc2=(length(a)+1):(length(a)+length(c));loc3=length(x)
      myf=function(x,data.main,data.cali) {
        get(name)(a=x[loc1],c=x[loc2],sc=x[loc3],data.main=data.main,data.cali=data.cali)
      }
    } else {
      x=c(b,d,sd);loc1=1:length(b);
      loc2=(length(b)+1):(length(b)+length(d));loc3=length(x)
      myf=function(x,data.main,data.cali) {
        get(name)(b=x[loc1],d=x[loc2],sd=x[loc3],data.main=data.main,data.cali=data.cali)
      }
    }
    mygrad(myf,x0=x,data.main=subset(mydata,I==0),data.cali=subset(mydata,I==1))
  }
  
  # number of covariates in outcome model and measurement error model
  nw=length(Wname);nz=length(Zname)
  # regression models for the outcome
  fcon=as.formula(paste(Yname,"~",Xsname,"+",Mname,ifelse(nw==0,"",paste("+",paste(Wname,collapse="+")))))
  fmar=as.formula(paste(Yname,"~",Xsname,ifelse(nw==0,"",paste("+",paste(Wname,collapse="+")))))
  # regression models for the measurement error
  fecon=as.formula(paste(Xname,"~",Xsname,"+",Mname,ifelse(nz==0,"",paste("+",paste(Zname,collapse="+")))))
  femar=as.formula(paste(Xname,"~",Xsname,ifelse(nz==0,"",paste("+",paste(Zname,collapse="+")))))
  
  # Marginal Outcome Model
  m1=glm(fmar,data=mydata,family=poisson(link = "log"))
  ahat0 = m1$coefficients
  # Conditional Outcome Model
  m2=glm(fcon,data=mydata,family=poisson(link = "log"))
  bhat0 = m2$coefficients
  # Mariginal measurement error model
  m3=lm(femar,data=subset(mydata,mydata[,Iname]==1))
  chat = m3$coefficients
  sigma2c = var(m3$residuals)
  # Conditional measurement error model
  m4=lm(fecon,data=subset(mydata,I==1))
  dhat = m4$coefficients
  sigma2d = var(m4$residuals)
  
  NIE.naive = ahat0[2]-bhat0[2]
  MP.naive = (ahat0[2]-bhat0[2])/ahat0[2]
  
  # corrected a 
  ahat=m1$coefficients
  ahat[2] = ahat0[2]/chat[2]
  ahat[1] = ahat0[1] - ahat[2]*chat[1] - 0.5*(ahat[2]^2)*sigma2c
  # corrected b
  bhat=m2$coefficients
  bhat[2] = bhat0[2]/dhat[2]
  bhat[3] = bhat0[3] - bhat[2]*dhat[3]
  bhat[1] = bhat0[1] - bhat[2]*dhat[1] - 0.5*(bhat[2]^2)*sigma2d
  
  NIE.3s=  ahat[2]-bhat[2]
  MP.3s = (ahat[2]-bhat[2])/ahat[2]
  
  ############ Ua
  Ua=function(a,c,sc,data.main,data.cali,type=1) {
    data.main=subset(mydata,mydata[,Iname]==0)
    data.cali = subset(mydata,mydata[,Iname]==1)
    n1=dim(data.main)[1];n2=dim(data.cali)[1]
    Wterm=`if`(nw==0,0, c(as.matrix(data.main[,Wname,drop=F]) %*% a[-c(1:2)]) )
    Zterm=`if`(nz==0,0, c(as.matrix(data.main[,Zname,drop=F]) %*% c[-c(1:2)]) )
    Wterm_cali=`if`(nw==0,0, c(as.matrix(data.cali[,Wname,drop=F]) %*% a[-c(1:2)]) )
    Zterm_cali=`if`(nz==0,0, c(as.matrix(data.cali[,Zname,drop=F]) %*% c[-c(1:2)]) )
    f1=function(x,V) {
      P=exp(a[1]+a[2]*x+Wterm)
      Pa=ifelse(P<1,P,0.9)
      if (V[1]=="x") {
        x*(data.main[,Yname]-P)*(Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
      } else {
        V*(data.main[,Yname]-P)*(Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
      }
    }
    f2=function(x) {
      P=exp(a[1]+a[2]*x+Wterm)
      Pa=ifelse(P<1,P,0.9)
      (Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    }
    # main dataset
    out = array(0,dim=c(2+nw,n1+n2))
    out[1,1:n1]=my.ghq(f1,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),V=1,order=15)/my.ghq(f2,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),order=15)
    out[2,1:n1]=my.ghq(f1,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),V="x",order=15)/my.ghq(f2,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),order=15)
    if (nw>0) {
      for (m in (1:nw)) {
        out[2+m,1:n1]= my.ghq(f1,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),V=data.main[,Wname[m]],order=15)/my.ghq(f2,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),order=15)
      }
    }
    P=exp(a[1]+a[2]*data.cali[,Xname]+Wterm_cali)
    out[1,(n1+1):(n1+n2)]=1*(data.cali[,Yname]-P)
    out[2,(n1+1):(n1+n2)]=data.cali[,Xname]*(data.cali$Y-P)
    if (nw>0) {
      for (m in (1:nw)) {
        out[2+m,(n1+1):(n1+n2)]=data.cali[,Wname[m]]*(data.cali$Y-P)
      }
    }
    if (type==2) {return(out)}
    apply(out,1,sum)
  }
  
  ############ Uc
  Uc=function(c,a,sc,data.main,data.cali,type=1) {
    data.main=subset(mydata,mydata[,Iname]==0)
    data.cali = subset(mydata,mydata[,Iname]==1)
    n1=dim(data.main)[1];n2=dim(data.cali)[1]
    Wterm=`if`(nw==0,0, c(as.matrix(data.main[,Wname,drop=F]) %*% a[-c(1:2)]) )
    Zterm=`if`(nz==0,0, c(as.matrix(data.main[,Zname,drop=F]) %*% c[-c(1:2)]) )
    Wterm_cali=`if`(nw==0,0, c(as.matrix(data.cali[,Wname,drop=F]) %*% a[-c(1:2)]) )
    Zterm_cali=`if`(nz==0,0, c(as.matrix(data.cali[,Zname,drop=F]) %*% c[-c(1:2)]) )
    f1=function(x,V) {
      P=exp(a[1]+a[2]*x+Wterm)
      Pa=ifelse(P<1,P,0.9)
      V*(x-c[1]-c[2]*data.main[,Xsname]-Zterm)*(Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    }
    f2=function(x) {
      P=exp(a[1]+a[2]*x+Wterm)
      Pa=ifelse(P<1,P,0.9)
      (Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    }
    # main dataset
    out = array(0,dim=c(2+nz,n1+n2))
    out[1,1:n1]=my.ghq(f1,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),V=1,order=15)/my.ghq(f2,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),order=15)
    out[2,1:n1]=my.ghq(f1,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),V=data.main[,Xsname],order=15)/my.ghq(f2,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),order=15)
    if (nz>0) {
      for (m in (1:nz)) {
        out[2+m,1:n1]=my.ghq(f1,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),V=data.main[,Zname[m]],order=15)/my.ghq(f2,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),order=15)
      }
    }
    out[1,(n1+1):(n1+n2)]=1*(data.cali[,Xname]-c[1]-c[2]*data.cali[,Xsname]-Zterm_cali)
    out[2,(n1+1):(n1+n2)]=data.cali[,Xsname]*(data.cali[,Xname]-c[1]-c[2]*data.cali[,Xsname]-Zterm_cali)
    if (nz>0) {
      for (m in (1:nz)) {
        out[2+m,(n1+1):(n1+n2)]=data.cali[,Wname[m]]*(data.cali[,Xname]-c[1]-c[2]*data.cali[,Xsname]-Zterm_cali)
      }
    }
    if (type==2) {return(out)}
    apply(out,1,sum)
  }
  
  ############ Usc
  Usc=function(sc,a,c,sc.o,data.main,data.cali,type=1) {
    data.main=subset(mydata,mydata[,Iname]==0)
    data.cali = subset(mydata,mydata[,Iname]==1)
    n1=dim(data.main)[1];n2=dim(data.cali)[1]
    Wterm=`if`(nw==0,0, c(as.matrix(data.main[,Wname,drop=F]) %*% a[-c(1:2)]) )
    Zterm=`if`(nz==0,0, c(as.matrix(data.main[,Zname,drop=F]) %*% c[-c(1:2)]) )
    Wterm_cali=`if`(nw==0,0, c(as.matrix(data.cali[,Wname,drop=F]) %*% a[-c(1:2)]) )
    Zterm_cali=`if`(nz==0,0, c(as.matrix(data.cali[,Zname,drop=F]) %*% c[-c(1:2)]) )
    f1=function(x) {
      P=exp(a[1]+a[2]*x+Wterm)
      Pa=ifelse(P<1,P,0.9)
      ((x-c[1]-c[2]*data.main[,Xsname]-Zterm)^2 - sc) * (Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    }
    f2=function(x) {
      P=exp(a[1]+a[2]*x+Wterm)
      Pa=ifelse(P<1,P,0.9)
      (Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    }
    # main dataset
    out = array(0,dim=c(1,n1+n2))
    out[1,1:n1]=my.ghq(f1,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),order=15)/my.ghq(f2,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),order=15)
    out[1,(n1+1):(n1+n2)]=(data.cali[,Xname]-c[1]-c[2]*data.cali[,Xsname]-Zterm_cali)^2 - sc
    if (type==2) {return(out)}
    apply(out,1,sum)
  }
  
  
  ############ Ub
  Ub=function(b,d,sd,data.main,data.cali,type=1) {
    data.main=subset(mydata,mydata[,Iname]==0)
    data.cali = subset(mydata,mydata[,Iname]==1)
    n1=dim(data.main)[1];n2=dim(data.cali)[1]
    Wterm=`if`(nw==0,0, c(as.matrix(data.main[,Wname,drop=F]) %*% b[-c(1:3)]) )
    Zterm=`if`(nz==0,0, c(as.matrix(data.main[,Zname,drop=F]) %*% d[-c(1:3)]) )
    Wterm_cali=`if`(nw==0,0, c(as.matrix(data.cali[,Wname,drop=F]) %*% b[-c(1:3)]) )
    Zterm_cali=`if`(nz==0,0, c(as.matrix(data.cali[,Zname,drop=F]) %*% d[-c(1:3)]) )
    f1=function(x,V) {
      P=exp(b[1]+b[2]*x+b[3]*data.main[,Mname]+Wterm)
      Pa=ifelse(P<1,P,0.9)
      if (V[1]=="x") {
        x*(data.main[,Yname]-P)*(Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
      } else {
        V*(data.main[,Yname]-P)*(Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
      }
    }
    f2=function(x) {
      P=exp(b[1]+b[2]*x+b[3]*data.main[,Mname]+Wterm)
      Pa=ifelse(P<1,P,0.9)
      (Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    }
    # main dataset
    out = array(0,dim=c(3+nw,n1+n2))
    mux = d[1]+d[2]*data.main[,Xsname] + d[3]*data.main[,Mname] + Zterm
    sdx = rep(sqrt(sd),n1)
    out[1,1:n1]=my.ghq(f1,mu=mux,sd=sdx,V=1,order=15)/my.ghq(f2,mu=mux,sd=sdx,order=15)
    out[2,1:n1]=my.ghq(f1,mu=mux,sd=sdx,V="x",order=15)/my.ghq(f2,mu=mux,sd=sdx,order=15)
    out[3,1:n1]=my.ghq(f1,mu=mux,sd=sdx,V=data.main[,Mname],order=15)/my.ghq(f2,mu=mux,sd=sdx,order=15)
    
    P=exp(b[1]+b[2]*data.cali[,Xname]+b[3]*data.cali[,Mname]+Wterm_cali)
    out[1,(n1+1):(n1+n2)]=1*(data.cali[,Yname]-P)
    out[2,(n1+1):(n1+n2)]=data.cali[,Xname]*(data.cali[,Yname]-P)
    out[3,(n1+1):(n1+n2)]=data.cali[,Mname]*(data.cali[,Yname]-P)
    if (nw>0) {
      for (m in (1:nw)) {
        out[3+m,1:n1]=my.ghq(f1,mu=mux,sd=sdx,V=data.main[,Wname[m]],order=15)/my.ghq(f2,mu=mux,sd=sdx,order=15)
        out[3+m,(n1+1):(n1+n2)]=data.cali[,Wname[m]]*(data.cali[,Yname]-P)
      }
    }
    if (type==2) {return(out)}
    apply(out,1,sum)
  }
  
  ############ Ud
  Ud=function(d,b,sd,data.main,data.cali,type=1) {
    data.main=subset(mydata,mydata[,Iname]==0)
    data.cali = subset(mydata,mydata[,Iname]==1)
    n1=dim(data.main)[1];n2=dim(data.cali)[1]
    Wterm=`if`(nw==0,0, c(as.matrix(data.main[,Wname,drop=F]) %*% b[-c(1:3)]) )
    Zterm=`if`(nz==0,0, c(as.matrix(data.main[,Zname,drop=F]) %*% d[-c(1:3)]) )
    Wterm_cali=`if`(nw==0,0, c(as.matrix(data.cali[,Wname,drop=F]) %*% b[-c(1:3)]) )
    Zterm_cali=`if`(nz==0,0, c(as.matrix(data.cali[,Zname,drop=F]) %*% d[-c(1:3)]) )
    f1=function(x,V) {
      P=exp(b[1]+b[2]*x+b[3]*data.main[,Mname]+Wterm)
      Pa=ifelse(P<1,P,0.9)
      V*(x-d[1]-d[2]*data.main[,Xsname]-d[3]*data.main[,Mname]-Zterm)*(Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    }
    f2=function(x) {
      P=exp(b[1]+b[2]*x+b[3]*data.main[,Mname]+Wterm)
      Pa=ifelse(P<1,P,0.9)
      (Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    }
    # main dataset
    out = array(0,dim=c(3+nz,n1+n2))
    mux = d[1]+d[2]*data.main[,Xsname] + d[3]*data.main[,Mname] + Zterm
    sdx = rep(sqrt(sd),n1)
    out[1,1:n1]=my.ghq(f1,mu=mux,sd=sdx,V=1,order=15)/my.ghq(f2,mu=mux,sd=sdx,order=15)
    out[2,1:n1]=my.ghq(f1,mu=mux,sd=sdx,V=data.main[,Xsname],order=15)/my.ghq(f2,mu=mux,sd=sdx,order=15)
    out[3,1:n1]=my.ghq(f1,mu=mux,sd=sdx,V=data.main[,Mname],order=15)/my.ghq(f2,mu=mux,sd=sdx,order=15)
    out[1,(n1+1):(n1+n2)]=1*(data.cali[,Xname]-d[1]-d[2]*data.cali[,Xsname]-d[3]*data.cali[,Mname]-Zterm_cali)
    out[2,(n1+1):(n1+n2)]=data.cali[,Xsname] * (data.cali[,Xname]-d[1]-d[2]*data.cali[,Xsname]-d[3]*data.cali[,Mname]-Zterm_cali)
    out[3,(n1+1):(n1+n2)]=data.cali[,Mname]  * (data.cali[,Xname]-d[1]-d[2]*data.cali[,Xsname]-d[3]*data.cali[,Mname]-Zterm_cali)
    if (nz>0) {
      for (m in (1:nz)) {
        out[3+m,1:n1]=my.ghq(f1,mu=mux,sd=sdx,V=data.main[,Zname[m]],order=15)/my.ghq(f2,mu=mux,sd=sdx,order=15)
        out[3+m,(n1+1):(n1+n2)]=data.cali[,Zname[m]]  * (data.cali[,Xname]-d[1]-d[2]*data.cali[,Xsname]-d[3]*data.cali[,Mname]-Zterm_cali)
      }
    }
    if (type==2) {return(out)}
    apply(out,1,sum)
  }
  
  ############ Usd
  Usd=function(sd,b,d,data.main,data.cali,type=1) {
    data.main=subset(mydata,mydata[,Iname]==0)
    data.cali = subset(mydata,mydata[,Iname]==1)
    n1=dim(data.main)[1];n2=dim(data.cali)[1]
    Wterm=`if`(nw==0,0, c(as.matrix(data.main[,Wname,drop=F]) %*% b[-c(1:3)]) )
    Zterm=`if`(nz==0,0, c(as.matrix(data.main[,Zname,drop=F]) %*% d[-c(1:3)]) )
    Wterm_cali=`if`(nw==0,0, c(as.matrix(data.cali[,Wname,drop=F]) %*% b[-c(1:3)]) )
    Zterm_cali=`if`(nz==0,0, c(as.matrix(data.cali[,Zname,drop=F]) %*% d[-c(1:3)]) )
    f1=function(x) {
      P=exp(b[1]+b[2]*x+b[3]*data.main[,Mname]+Wterm)
      Pa=ifelse(P<1,P,0.9)
      ((x-d[1]-d[2]*data.main[,Xsname] - d[3]*data.main[,Mname]-Zterm)^2 - sd) * (Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    }
    f2=function(x) {
      P=exp(b[1]+b[2]*x+b[3]*data.main[,Mname]+Wterm)
      Pa=ifelse(P<1,P,0.9)
      (Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    }
    # main dataset
    out = array(0,dim=c(1,n1+n2))
    out[1,1:n1]=my.ghq(f1,mu=d[1]+d[2]*data.main[,Xsname]+d[3]*data.main[,Mname]+Zterm,sd=rep(sqrt(sd),n1),order=15)/my.ghq(f2,mu=d[1]+d[2]*data.main[,Xsname]+d[3]*data.main[,Mname]+Zterm,sd=rep(sqrt(sd),n1),order=15)
    out[1,(n1+1):(n1+n2)]=(data.cali[,Xname]-d[1]-d[2]*data.cali[,Xsname] - d[3]*data.cali[,Mname]-Zterm_cali)^2 - sd
    if (type==2) {return(out)}
    apply(out,1,sum)
  }
  
  
  ############ estimation of parameters in marginal model
  par.vec=c(ahat,chat,sigma2c)
  out = matrix(NA,ncol=length(par.vec),nrow=10)
  out[1,]=par.vec
  for (i in (2:10)) {
    ahat=nleqslv(x=ahat,fn=Ua,c=chat,sc=sigma2c,
                 data.main=subset(mydata,mydata[,Iname]==0),data.cali=subset(mydata,mydata[,Iname]==1),
                 control=list(maxit=10))$x
    chat=nleqslv(x=chat,fn=Uc,a=ahat,sc=sigma2c,
                 data.main=subset(mydata,mydata[,Iname]==0),
                 data.cali=subset(mydata,mydata[,Iname]==1),control=list(maxit=10))$x
    sigma2c=nleqslv(x=sigma2c,fn=Usc,a=ahat,c=chat,sc.o=sigma2c,
                    data.main=subset(mydata,mydata[,Iname]==0),
                    data.cali=subset(mydata,mydata[,Iname]==1),control=list(maxit=10))$x
    par.vec=c(ahat,chat,sigma2c)
    out[i,]=par.vec
    k=i
    if (max(abs(out[i,]-out[i-1,]))<0.001) {break}
  }
  par_mar=out[k,]
  aW=`if`(nw==0,NULL,paste("aW",1:nw,sep=""))
  cZ=`if`(nz==0,NULL,paste("cZ",1:nz,sep=""))
  names(par_mar)=c("a0","a1",aW,"c0","c1",cZ,"sc")
  
  ############ estimation of parameters in conditional model
  par.vec=c(bhat,dhat,sigma2d)
  out = matrix(NA,ncol=length(par.vec),nrow=10)
  out[1,]=par.vec
  for (i in (2:10)) {
    bhat=nleqslv(x=bhat,fn=Ub,d=dhat,sd=sigma2d,
                 data.main=subset(mydata,mydata[,Iname]==0),data.cali=subset(mydata,mydata[,Iname]==1),
                 control=list(maxit=10))$x
    dhat=nleqslv(x=dhat,fn=Ud,b=bhat,sd=sigma2d,
                 data.main=subset(mydata,mydata[,Iname]==0),
                 data.cali=subset(mydata,mydata[,Iname]==1),
                 control=list(maxit=10))$x
    sigma2d=nleqslv(x=sigma2d,fn=Usd,b=bhat,d=dhat,
                    data.main=subset(mydata,mydata[,Iname]==0),
                    data.cali=subset(mydata,mydata[,Iname]==1),
                    control=list(maxit=10))$x
    par.vec=c(bhat,dhat,sigma2d)
    out[i,]=par.vec
    k=i
    if (max(abs(out[i,]-out[i-1,]))<0.001) {break}
  }
  par_con=out[k,]
  bW=`if`(nw==0,NULL,paste("bW",1:nw,sep=""))
  dZ=`if`(nz==0,NULL,paste("dZ",1:nz,sep=""))
  names(par_con)=c("b0","b1","b2",bW,"d0","d1","d2",dZ,"sd")
  
  NIE.eee = par_mar[2] - par_con[2]
  MP.eee  =  (par_mar[2] - par_con[2])/par_mar[2]
  
  
  
  
  ### sandwich variance A^-1 B A^-1
  # (1) calculate B
  Ba=Ua(a=ahat,c=chat,sc=sigma2c,data.main=subset(mydata,mydata[,Iname]==0),data.cali=subset(mydata,mydata[,Iname]==1),type=2)
  Bc=Uc(c=chat,a=ahat,sc=sigma2c,data.main=subset(mydata,mydata[,Iname]==0),data.cali=subset(mydata,mydata[,Iname]==1),type=2)
  Bsc=Usc(sc=sigma2c,a=ahat,c=chat,data.main=subset(mydata,mydata[,Iname]==0),data.cali=subset(mydata,mydata[,Iname]==1),type=2)
  B1=rbind(Ba,Bc,Bsc)
  Bb =Ub(b=bhat,d=dhat,sd=sigma2d,data.main=subset(mydata,mydata[,Iname]==0),data.cali=subset(mydata,mydata[,Iname]==1),type=2)
  Bd =Ud(d=dhat,b=bhat,sd=sigma2d,data.main=subset(mydata,mydata[,Iname]==0),data.cali=subset(mydata,mydata[,Iname]==1),type=2)
  Bsd=Usd(sd=sigma2d,b=bhat,d=dhat,data.main=subset(mydata,mydata[,Iname]==0),data.cali=subset(mydata,mydata[,Iname]==1),type=2)
  B2=rbind(Bb,Bd,Bsd)
  Bs=rbind(B1,B2)
  B=Bs %*% t(Bs)/dim(mydata)[1]
  # (2) calculate A
  Aa=gradU(name="Ua",a=ahat,b=bhat,c=chat,sc=sigma2c,d=dhat,sd=sigma2d,model=1)
  Ac=gradU(name="Uc",a=ahat,b=bhat,c=chat,sc=sigma2c,d=dhat,sd=sigma2d,model=1)
  Asc=gradU(name="Usc",a=ahat,b=bhat,c=chat,sc=sigma2c,d=dhat,sd=sigma2d,model=1)
  Ab=gradU(name="Ub",a=ahat,b=bhat,c=chat,sc=sigma2c,d=dhat,sd=sigma2d,model=2)
  Ad=gradU(name="Ud",a=ahat,b=bhat,c=chat,sc=sigma2c,d=dhat,sd=sigma2d,model=2)
  Asd=gradU(name="Usd",a=ahat,b=bhat,c=chat,sc=sigma2c,d=dhat,sd=sigma2d,model=2)
  A1=cbind(Aa,Ac,Asc)
  A2=cbind(Ab,Ad,Asd)
  A=rbind(cbind(t(A1), matrix(0, nrow=nrow(A1), ncol=ncol(A2))),
          cbind(matrix(0, nrow=nrow(A2), ncol=ncol(A1)), t(A2)))/dim(mydata)[1]
  Ai=solve(A)
  covP=(Ai %*% B %*% t(Ai))
  colnames(covP)=rownames(covP)=c("a0","a1",aW,"c0","c1",cZ,"sc",
                                  "b0","b1","b2",bW,"d0","d1","d2",dZ,"sd")
  covb<-covP[c("a1","b1"),c("a1","b1")]/dim(mydata)[1]
  
  
  # SE of MP and NIE
  se.MPeee=sqrt(as.numeric(t(c(bhat[2]/(ahat[2]^2),-1/ahat[2])) %*% covb %*% c(bhat[2]/(ahat[2]^2),-1/ahat[2])))
  se.NIEeee=sqrt(as.numeric(t(c(1,-1)) %*% covb %*% c(1,-1)))
  output=c(NIE.eee,MP.eee,se.NIEeee,se.MPeee)
  names(output)=c("NIE_eee","MP_eee","seNIE_eee","seMP_eee")
  output
}

# Functions for the reduced EEE
get_eee2=function(data,Xname="X",Xsname="Xs",Yname="Y",Mname="M",Wname="W",Zname="W",Iname="I") {
  #Xname="X";Xsname="Xs";Yname="Y";Mname="M";Wname="W";Zname="W";Iname="I"
  mydata=data
  
  gradU=function(name="Ua",a,sa,b,sb,c,sc,d,sd,model=1) {
    mygrad=function (f, x0, heps = 1e-6, ...) {
      fun <- match.fun(f)
      f <- function(x) fun(x, ...)
      p =length(f(x0))
      n <- length(x0)
      hh <- rep(0, n)
      gr <- matrix(0,nrow=n,ncol=p)
      for (i in 1:n) {
        hh[i] <- heps
        gr[i,] <- (f(x0 + hh) - f(x0 - hh))/(2 * heps)
        hh[i] <- 0
      }
      return(gr)
    }
    if (model==1) {
      x=c(a,c,sc);loc1=1:length(a);
      loc2=(length(a)+1):(length(a)+length(c));loc3=length(x)
      myf=function(x,data.main,data.cali) {
        get(name)(a=x[loc1],c=x[loc2],sc=x[loc3],data.main=data.main,data.cali=data.cali)
      }
    } else {
      x=c(b,d,sd);loc1=1:length(b);
      loc2=(length(b)+1):(length(b)+length(d));loc3=length(x)
      myf=function(x,data.main,data.cali) {
        get(name)(b=x[loc1],d=x[loc2],sd=x[loc3],data.main=data.main,data.cali=data.cali)
      }
    }
    mygrad(myf,x0=x,data.main=subset(mydata,I==0),data.cali=subset(mydata,I==1))
  }
  
  # number of covariates in outcome model and measurement error model
  nw=length(Wname);nz=length(Zname)
  # regression models for the outcome
  fcon=as.formula(paste(Yname,"~",Xsname,"+",Mname,ifelse(nw==0,"",paste("+",paste(Wname,collapse="+")))))
  fmar=as.formula(paste(Yname,"~",Xsname,ifelse(nw==0,"",paste("+",paste(Wname,collapse="+")))))
  # regression models for the measurement error
  fecon=as.formula(paste(Xname,"~",Xsname,"+",Mname,ifelse(nz==0,"",paste("+",paste(Zname,collapse="+")))))
  femar=as.formula(paste(Xname,"~",Xsname,ifelse(nz==0,"",paste("+",paste(Zname,collapse="+")))))
  
  # Marginal Outcome Model
  m1=glm(fmar,data=mydata,family=poisson(link = "log"))
  ahat0 = m1$coefficients
  # Conditional Outcome Model
  m2=glm(fcon,data=mydata,family=poisson(link = "log"))
  bhat0 = m2$coefficients
  # Mariginal measurement error model
  m3=lm(femar,data=subset(mydata,mydata[,Iname]==1))
  chat = m3$coefficients
  sigma2c = var(m3$residuals)
  # Conditional measurement error model
  m4=lm(fecon,data=subset(mydata,I==1))
  dhat = m4$coefficients
  sigma2d = var(m4$residuals)
  
  NIE.naive = ahat0[2]-bhat0[2]
  MP.naive = (ahat0[2]-bhat0[2])/ahat0[2]
  
  # corrected a 
  ahat=m1$coefficients
  ahat[2] = ahat0[2]/chat[2]
  ahat[1] = ahat0[1] - ahat[2]*chat[1] - 0.5*(ahat[2]^2)*sigma2c
  # corrected b
  bhat=m2$coefficients
  bhat[2] = bhat0[2]/dhat[2]
  bhat[3] = bhat0[3] - bhat[2]*dhat[3]
  bhat[1] = bhat0[1] - bhat[2]*dhat[1] - 0.5*(bhat[2]^2)*sigma2d
  
  NIE.3s=  ahat[2]-bhat[2]
  MP.3s = (ahat[2]-bhat[2])/ahat[2]
  
  ############ Ua
  Ua=function(a,c,sc,data.main,data.cali,type=1) {
    data.main=subset(mydata,mydata[,Iname]==0)
    data.cali = subset(mydata,mydata[,Iname]==1)
    n1=dim(data.main)[1];n2=dim(data.cali)[1]
    Wterm=`if`(nw==0,0, c(as.matrix(data.main[,Wname,drop=F]) %*% a[-c(1:2)]) )
    Zterm=`if`(nz==0,0, c(as.matrix(data.main[,Zname,drop=F]) %*% c[-c(1:2)]) )
    Wterm_cali=`if`(nw==0,0, c(as.matrix(data.cali[,Wname,drop=F]) %*% a[-c(1:2)]) )
    Zterm_cali=`if`(nz==0,0, c(as.matrix(data.cali[,Zname,drop=F]) %*% c[-c(1:2)]) )
    f1=function(x,V) {
      P=exp(a[1]+a[2]*x+Wterm)
      Pa=ifelse(P<1,P,0.9)
      if (V[1]=="x") {
        x*(data.main[,Yname]-P)*(Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
      } else {
        V*(data.main[,Yname]-P)*(Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
      }
    }
    f2=function(x) {
      P=exp(a[1]+a[2]*x+Wterm)
      Pa=ifelse(P<1,P,0.9)
      (Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    }
    # main dataset
    out = array(0,dim=c(2+nw,n1+n2))
    out[1,1:n1]=my.ghq(f1,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),V=1,order=15)/my.ghq(f2,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),order=15)
    out[2,1:n1]=my.ghq(f1,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),V="x",order=15)/my.ghq(f2,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),order=15)
    if (nw>0) {
      for (m in (1:nw)) {
        out[2+m,1:n1]= my.ghq(f1,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),V=data.main[,Wname[m]],order=15)/my.ghq(f2,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),order=15)
      }
    }
    P=exp(a[1]+a[2]*data.cali[,Xname]+Wterm_cali)
    out[1,(n1+1):(n1+n2)]=1*(data.cali[,Yname]-P)
    out[2,(n1+1):(n1+n2)]=data.cali[,Xname]*(data.cali$Y-P)
    if (nw>0) {
      for (m in (1:nw)) {
        out[2+m,(n1+1):(n1+n2)]=data.cali[,Wname[m]]*(data.cali$Y-P)
      }
    }
    if (type==2) {return(out)}
    apply(out,1,sum)
  }
  
  ############ Uc
  Uc=function(c,a,sc,data.main,data.cali,type=1) {
    data.main=subset(mydata,mydata[,Iname]==0)
    data.cali = subset(mydata,mydata[,Iname]==1)
    n1=dim(data.main)[1];n2=dim(data.cali)[1]
    #Wterm=`if`(nw==0,0, c(as.matrix(data.main[,Wname,drop=F]) %*% a[-c(1:2)]) )
    #Zterm=`if`(nz==0,0, c(as.matrix(data.main[,Zname,drop=F]) %*% c[-c(1:2)]) )
    Wterm_cali=`if`(nw==0,0, c(as.matrix(data.cali[,Wname,drop=F]) %*% a[-c(1:2)]) )
    Zterm_cali=`if`(nz==0,0, c(as.matrix(data.cali[,Zname,drop=F]) %*% c[-c(1:2)]) )
    # f1=function(x,V) {
    #   P=exp(a[1]+a[2]*x+Wterm)
    #   Pa=ifelse(P<1,P,0.9)
    #   V*(x-c[1]-c[2]*data.main[,Xsname]-Zterm)*(Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    # }
    # f2=function(x) {
    #   P=exp(a[1]+a[2]*x+Wterm)
    #   Pa=ifelse(P<1,P,0.9)
    #   (Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    # }
    # main dataset
    out = array(0,dim=c(2+nz,n1+n2))
    #out[1,1:n1]=my.ghq(f1,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),V=1,order=15)/my.ghq(f2,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),order=15)
    #out[2,1:n1]=my.ghq(f1,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),V=data.main[,Xsname],order=15)/my.ghq(f2,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),order=15)
    # if (nz>0) {
    #   for (m in (1:nz)) {
    #     out[2+m,1:n1]=my.ghq(f1,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),V=data.main[,Zname[m]],order=15)/my.ghq(f2,mu=c[1]+c[2]*data.main[,Zname[m]]+Zterm,sd=rep(sqrt(sc),n1),order=15)
    #   }
    # }
    out[1,(n1+1):(n1+n2)]=1*(data.cali[,Xname]-c[1]-c[2]*data.cali[,Xsname]-Zterm_cali)
    out[2,(n1+1):(n1+n2)]=data.cali[,Xsname]*(data.cali[,Xname]-c[1]-c[2]*data.cali[,Xsname]-Zterm_cali)
    if (nz>0) {
      for (m in (1:nz)) {
        out[2+m,(n1+1):(n1+n2)]=data.cali[,Wname[m]]*(data.cali[,Xname]-c[1]-c[2]*data.cali[,Xsname]-Zterm_cali)
      }
    }
    if (type==2) {return(out)}
    apply(out,1,sum)
  }
  
  ############ Usc
  Usc=function(sc,a,c,sc.o,data.main,data.cali,type=1) {
    data.main=subset(mydata,mydata[,Iname]==0)
    data.cali = subset(mydata,mydata[,Iname]==1)
    n1=dim(data.main)[1];n2=dim(data.cali)[1]
    #Wterm=`if`(nw==0,0, c(as.matrix(data.main[,Wname,drop=F]) %*% a[-c(1:2)]) )
    #Zterm=`if`(nz==0,0, c(as.matrix(data.main[,Zname,drop=F]) %*% c[-c(1:2)]) )
    Wterm_cali=`if`(nw==0,0, c(as.matrix(data.cali[,Wname,drop=F]) %*% a[-c(1:2)]) )
    Zterm_cali=`if`(nz==0,0, c(as.matrix(data.cali[,Zname,drop=F]) %*% c[-c(1:2)]) )
    # f1=function(x) {
    #   P=exp(a[1]+a[2]*x+Wterm)
    #   Pa=ifelse(P<1,P,0.9)
    #   ((x-c[1]-c[2]*data.main[,Xsname]-Zterm)^2 - sc) * (Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    # }
    # f2=function(x) {
    #   P=exp(a[1]+a[2]*x+Wterm)
    #   Pa=ifelse(P<1,P,0.9)
    #   (Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    # }
    # main dataset
    out = array(0,dim=c(1,n1+n2))
    #out[1,1:n1]=my.ghq(f1,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),order=15)/my.ghq(f2,mu=c[1]+c[2]*data.main[,Xsname]+Zterm,sd=rep(sqrt(sc),n1),order=15)
    out[1,(n1+1):(n1+n2)]=(data.cali[,Xname]-c[1]-c[2]*data.cali[,Xsname]-Zterm_cali)^2 - sc
    if (type==2) {return(out)}
    apply(out,1,sum)
  }
  
  
  ############ Ub
  Ub=function(b,d,sd,data.main,data.cali,type=1) {
    data.main=subset(mydata,mydata[,Iname]==0)
    data.cali = subset(mydata,mydata[,Iname]==1)
    n1=dim(data.main)[1];n2=dim(data.cali)[1]
    Wterm=`if`(nw==0,0, c(as.matrix(data.main[,Wname,drop=F]) %*% b[-c(1:3)]) )
    Zterm=`if`(nz==0,0, c(as.matrix(data.main[,Zname,drop=F]) %*% d[-c(1:3)]) )
    Wterm_cali=`if`(nw==0,0, c(as.matrix(data.cali[,Wname,drop=F]) %*% b[-c(1:3)]) )
    Zterm_cali=`if`(nz==0,0, c(as.matrix(data.cali[,Zname,drop=F]) %*% d[-c(1:3)]) )
    f1=function(x,V) {
      P=exp(b[1]+b[2]*x+b[3]*data.main[,Mname]+Wterm)
      Pa=ifelse(P<1,P,0.9)
      if (V[1]=="x") {
        x*(data.main[,Yname]-P)*(Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
      } else {
        V*(data.main[,Yname]-P)*(Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
      }
    }
    f2=function(x) {
      P=exp(b[1]+b[2]*x+b[3]*data.main[,Mname]+Wterm)
      Pa=ifelse(P<1,P,0.9)
      (Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    }
    # main dataset
    out = array(0,dim=c(3+nw,n1+n2))
    mux = d[1]+d[2]*data.main[,Xsname] + d[3]*data.main[,Mname] + Zterm
    sdx = rep(sqrt(sd),n1)
    out[1,1:n1]=my.ghq(f1,mu=mux,sd=sdx,V=1,order=15)/my.ghq(f2,mu=mux,sd=sdx,order=15)
    out[2,1:n1]=my.ghq(f1,mu=mux,sd=sdx,V="x",order=15)/my.ghq(f2,mu=mux,sd=sdx,order=15)
    out[3,1:n1]=my.ghq(f1,mu=mux,sd=sdx,V=data.main[,Mname],order=15)/my.ghq(f2,mu=mux,sd=sdx,order=15)
    
    P=exp(b[1]+b[2]*data.cali[,Xname]+b[3]*data.cali[,Mname]+Wterm_cali)
    out[1,(n1+1):(n1+n2)]=1*(data.cali[,Yname]-P)
    out[2,(n1+1):(n1+n2)]=data.cali[,Xname]*(data.cali[,Yname]-P)
    out[3,(n1+1):(n1+n2)]=data.cali[,Mname]*(data.cali[,Yname]-P)
    if (nw>0) {
      for (m in (1:nw)) {
        out[3+m,1:n1]=my.ghq(f1,mu=mux,sd=sdx,V=data.main[,Wname[m]],order=15)/my.ghq(f2,mu=mux,sd=sdx,order=15)
        out[3+m,(n1+1):(n1+n2)]=data.cali[,Wname[m]]*(data.cali[,Yname]-P)
      }
    }
    if (type==2) {return(out)}
    apply(out,1,sum)
  }
  
  ############ Ud
  Ud=function(d,b,sd,data.main,data.cali,type=1) {
    data.main=subset(mydata,mydata[,Iname]==0)
    data.cali = subset(mydata,mydata[,Iname]==1)
    n1=dim(data.main)[1];n2=dim(data.cali)[1]
    #Wterm=`if`(nw==0,0, c(as.matrix(data.main[,Wname,drop=F]) %*% b[-c(1:3)]) )
    #Zterm=`if`(nz==0,0, c(as.matrix(data.main[,Zname,drop=F]) %*% d[-c(1:3)]) )
    Wterm_cali=`if`(nw==0,0, c(as.matrix(data.cali[,Wname,drop=F]) %*% b[-c(1:3)]) )
    Zterm_cali=`if`(nz==0,0, c(as.matrix(data.cali[,Zname,drop=F]) %*% d[-c(1:3)]) )
    # f1=function(x,V) {
    #   P=exp(b[1]+b[2]*x+b[3]*data.main[,Mname]+Wterm)
    #   Pa=ifelse(P<1,P,0.9)
    #   V*(x-d[1]-d[2]*data.main[,Xsname]-d[3]*data.main[,Mname]-Zterm)*(Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    # }
    # f2=function(x) {
    #   P=exp(b[1]+b[2]*x+b[3]*data.main[,Mname]+Wterm)
    #   Pa=ifelse(P<1,P,0.9)
    #   (Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    # }
    # main dataset
    out = array(0,dim=c(3+nz,n1+n2))
    # mux = d[1]+d[2]*data.main[,Xsname] + d[3]*data.main[,Mname] + Zterm
    # sdx = rep(sqrt(sd),n1)
    # out[1,1:n1]=my.ghq(f1,mu=mux,sd=sdx,V=1,order=15)/my.ghq(f2,mu=mux,sd=sdx,order=15)
    # out[2,1:n1]=my.ghq(f1,mu=mux,sd=sdx,V=data.main[,Xsname],order=15)/my.ghq(f2,mu=mux,sd=sdx,order=15)
    # out[3,1:n1]=my.ghq(f1,mu=mux,sd=sdx,V=data.main[,Mname],order=15)/my.ghq(f2,mu=mux,sd=sdx,order=15)
    out[1,(n1+1):(n1+n2)]=1*(data.cali[,Xname]-d[1]-d[2]*data.cali[,Xsname]-d[3]*data.cali[,Mname]-Zterm_cali)
    out[2,(n1+1):(n1+n2)]=data.cali[,Xsname] * (data.cali[,Xname]-d[1]-d[2]*data.cali[,Xsname]-d[3]*data.cali[,Mname]-Zterm_cali)
    out[3,(n1+1):(n1+n2)]=data.cali[,Mname]  * (data.cali[,Xname]-d[1]-d[2]*data.cali[,Xsname]-d[3]*data.cali[,Mname]-Zterm_cali)
    if (nz>0) {
      for (m in (1:nz)) {
        #out[3+m,1:n1]=my.ghq(f1,mu=mux,sd=sdx,V=data.main[,Zname[m]],order=15)/my.ghq(f2,mu=mux,sd=sdx,order=15)
        out[3+m,(n1+1):(n1+n2)]=data.cali[,Zname[m]]  * (data.cali[,Xname]-d[1]-d[2]*data.cali[,Xsname]-d[3]*data.cali[,Mname]-Zterm_cali)
      }
    }
    if (type==2) {return(out)}
    apply(out,1,sum)
  }
  
  ############ Usd
  Usd=function(sd,b,d,data.main,data.cali,type=1) {
    data.main=subset(mydata,mydata[,Iname]==0)
    data.cali = subset(mydata,mydata[,Iname]==1)
    n1=dim(data.main)[1];n2=dim(data.cali)[1]
    #Wterm=`if`(nw==0,0, c(as.matrix(data.main[,Wname,drop=F]) %*% b[-c(1:3)]) )
    #Zterm=`if`(nz==0,0, c(as.matrix(data.main[,Zname,drop=F]) %*% d[-c(1:3)]) )
    Wterm_cali=`if`(nw==0,0, c(as.matrix(data.cali[,Wname,drop=F]) %*% b[-c(1:3)]) )
    Zterm_cali=`if`(nz==0,0, c(as.matrix(data.cali[,Zname,drop=F]) %*% d[-c(1:3)]) )
    # f1=function(x) {
    #   P=exp(b[1]+b[2]*x+b[3]*data.main[,Mname]+Wterm)
    #   Pa=ifelse(P<1,P,0.9)
    #   ((x-d[1]-d[2]*data.main[,Xsname] - d[3]*data.main[,Mname]-Zterm)^2 - sd) * (Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    # }
    # f2=function(x) {
    #   P=exp(b[1]+b[2]*x+b[3]*data.main[,Mname]+Wterm)
    #   Pa=ifelse(P<1,P,0.9)
    #   (Pa^data.main[,Yname])*((1-Pa)^(1-data.main[,Yname]))
    # }
    # main dataset
    out = array(0,dim=c(1,n1+n2))
    #out[1,1:n1]=my.ghq(f1,mu=d[1]+d[2]*data.main[,Xsname]+d[3]*data.main[,Mname]+Zterm,sd=rep(sqrt(sd),n1),order=15)/my.ghq(f2,mu=d[1]+d[2]*data.main[,Xsname]+d[3]*data.main[,Mname]+Zterm,sd=rep(sqrt(sd),n1),order=15)
    out[1,(n1+1):(n1+n2)]=(data.cali[,Xname]-d[1]-d[2]*data.cali[,Xsname] - d[3]*data.cali[,Mname]-Zterm_cali)^2 - sd
    if (type==2) {return(out)}
    apply(out,1,sum)
  }
  
  
  ############ estimation of parameters in marginal model
  par.vec=c(ahat,chat,sigma2c)
  out = matrix(NA,ncol=length(par.vec),nrow=10)
  out[1,]=par.vec
  for (i in (2:10)) {
    chat=nleqslv(x=chat,fn=Uc,a=ahat,sc=sigma2c,
                 data.main=subset(mydata,mydata[,Iname]==0),
                 data.cali=subset(mydata,mydata[,Iname]==1),control=list(maxit=50))$x
    sigma2c=nleqslv(x=sigma2c,fn=Usc,a=ahat,c=chat,sc.o=sigma2c,
                    data.main=subset(mydata,mydata[,Iname]==0),
                    data.cali=subset(mydata,mydata[,Iname]==1),control=list(maxit=50))$x
    ahat=nleqslv(x=ahat,fn=Ua,c=chat,sc=sigma2c,
                 data.main=subset(mydata,mydata[,Iname]==0),data.cali=subset(mydata,mydata[,Iname]==1),
                 control=list(maxit=50))$x
    par.vec=c(ahat,chat,sigma2c)
    out[i,]=par.vec
    k=i
    if (max(abs(out[i,]-out[i-1,]))<0.001) {break}
  }
  #print(out)
  par_mar=out[k,]
  aW=`if`(nw==0,NULL,paste("aW",1:nw,sep=""))
  cZ=`if`(nz==0,NULL,paste("cZ",1:nz,sep=""))
  names(par_mar)=c("a0","a1",aW,"c0","c1",cZ,"sc")
  
  ############ estimation of parameters in conditional model
  par.vec=c(bhat,dhat,sigma2d)
  out = matrix(NA,ncol=length(par.vec),nrow=10)
  out[1,]=par.vec
  for (i in (2:10)) {
    sigma2d=nleqslv(x=sigma2d,fn=Usd,b=bhat,d=dhat,
                    data.main=subset(mydata,mydata[,Iname]==0),
                    data.cali=subset(mydata,mydata[,Iname]==1),
                    control=list(maxit=50))$x
    bhat=nleqslv(x=bhat,fn=Ub,d=dhat,sd=sigma2d,
                 data.main=subset(mydata,mydata[,Iname]==0),data.cali=subset(mydata,mydata[,Iname]==1),
                 control=list(maxit=50))$x
    dhat=nleqslv(x=dhat,fn=Ud,b=bhat,sd=sigma2d,
                 data.main=subset(mydata,mydata[,Iname]==0),
                 data.cali=subset(mydata,mydata[,Iname]==1),
                 control=list(maxit=50))$x
    par.vec=c(bhat,dhat,sigma2d)
    out[i,]=par.vec
    k=i
    if (max(abs(out[i,]-out[i-1,]))<0.001) {break}
  }
  #print(out)
  par_con=out[k,]
  bW=`if`(nw==0,NULL,paste("bW",1:nw,sep=""))
  dZ=`if`(nz==0,NULL,paste("dZ",1:nz,sep=""))
  names(par_con)=c("b0","b1","b2",bW,"d0","d1","d2",dZ,"sd")
  
  NIE.eee = par_mar[2] - par_con[2]
  MP.eee  =  (par_mar[2] - par_con[2])/par_mar[2]
  
  
  
  
  ### sandwich variance A^-1 B A^-1
  # (1) calculate B
  Ba=Ua(a=ahat,c=chat,sc=sigma2c,data.main=subset(mydata,mydata[,Iname]==0),data.cali=subset(mydata,mydata[,Iname]==1),type=2)
  Bc=Uc(c=chat,a=ahat,sc=sigma2c,data.main=subset(mydata,mydata[,Iname]==0),data.cali=subset(mydata,mydata[,Iname]==1),type=2)
  Bsc=Usc(sc=sigma2c,a=ahat,c=chat,data.main=subset(mydata,mydata[,Iname]==0),data.cali=subset(mydata,mydata[,Iname]==1),type=2)
  B1=rbind(Ba,Bc,Bsc)
  Bb =Ub(b=bhat,d=dhat,sd=sigma2d,data.main=subset(mydata,mydata[,Iname]==0),data.cali=subset(mydata,mydata[,Iname]==1),type=2)
  Bd =Ud(d=dhat,b=bhat,sd=sigma2d,data.main=subset(mydata,mydata[,Iname]==0),data.cali=subset(mydata,mydata[,Iname]==1),type=2)
  Bsd=Usd(sd=sigma2d,b=bhat,d=dhat,data.main=subset(mydata,mydata[,Iname]==0),data.cali=subset(mydata,mydata[,Iname]==1),type=2)
  B2=rbind(Bb,Bd,Bsd)
  Bs=rbind(B1,B2)
  B=Bs %*% t(Bs)/dim(mydata)[1]
  # (2) calculate A
  Aa=gradU(name="Ua",a=ahat,b=bhat,c=chat,sc=sigma2c,d=dhat,sd=sigma2d,model=1)
  Ac=gradU(name="Uc",a=ahat,b=bhat,c=chat,sc=sigma2c,d=dhat,sd=sigma2d,model=1)
  Asc=gradU(name="Usc",a=ahat,b=bhat,c=chat,sc=sigma2c,d=dhat,sd=sigma2d,model=1)
  Ab=gradU(name="Ub",a=ahat,b=bhat,c=chat,sc=sigma2c,d=dhat,sd=sigma2d,model=2)
  Ad=gradU(name="Ud",a=ahat,b=bhat,c=chat,sc=sigma2c,d=dhat,sd=sigma2d,model=2)
  Asd=gradU(name="Usd",a=ahat,b=bhat,c=chat,sc=sigma2c,d=dhat,sd=sigma2d,model=2)
  A1=cbind(Aa,Ac,Asc)
  A2=cbind(Ab,Ad,Asd)
  A=rbind(cbind(t(A1), matrix(0, nrow=nrow(A1), ncol=ncol(A2))),
          cbind(matrix(0, nrow=nrow(A2), ncol=ncol(A1)), t(A2)))/dim(mydata)[1]
  Ai=solve(A)
  covP=(Ai %*% B %*% t(Ai))
  colnames(covP)=rownames(covP)=c("a0","a1",aW,"c0","c1",cZ,"sc",
                                  "b0","b1","b2",bW,"d0","d1","d2",dZ,"sd")
  covb<-covP[c("a1","b1"),c("a1","b1")]/dim(mydata)[1]
  
  
  # SE of MP and NIE
  se.MPeee=sqrt(as.numeric(t(c(bhat[2]/(ahat[2]^2),-1/ahat[2])) %*% covb %*% c(bhat[2]/(ahat[2]^2),-1/ahat[2])))
  se.NIEeee=sqrt(as.numeric(t(c(1,-1)) %*% covb %*% c(1,-1)))
  output=c(NIE.eee,MP.eee,se.NIEeee,se.MPeee)
  names(output)=c("NIE_eee","MP_eee","seNIE_eee","seMP_eee")
  output
}

# Functions for the regression calibration
get_rc=function(data=mydata,Xname="X",Xsname="Xs",Yname="Y",Mname="M",Wname="W",Zname="W",Iname="I",outcome="binary") {
  #data=mydata;outcome="binary"
  mydata=data
  # number of covariates in outcome model and measurement error model
  nw=length(Wname);nz=length(Zname)
  # regression models for the outcome
  fcon=as.formula(paste(Yname,"~",Xsname,"+",Mname,ifelse(nw==0,"",paste("+",paste(Wname,collapse="+")))))
  fmar=as.formula(paste(Yname,"~",Xsname,ifelse(nw==0,"",paste("+",paste(Wname,collapse="+")))))
  
  fconreal=as.formula(paste(Yname,"~",Xname,"+",Mname,ifelse(nw==0,"",paste("+",paste(Wname,collapse="+")))))
  fmarreal=as.formula(paste(Yname,"~",Xname,ifelse(nw==0,"",paste("+",paste(Wname,collapse="+")))))
  # regression models for the measurement error
  fecon=as.formula(paste(Xname,"~",Xsname,"+",Mname,ifelse(nz==0,"",paste("+",paste(Wname,collapse="+")))))
  femar=as.formula(paste(Xname,"~",Xsname,ifelse(nz==0,"",paste("+",paste(Wname,collapse="+")))))
  
  # measurement error parameters
  calib.m1=lm(femar,data=subset(data,data[,Iname]==1))
  calib.m2=lm(fecon,data=subset(data,data[,Iname]==1))
  cd = c(calib.m1$coefficients[2],calib.m2$coefficients[2])
  
  # main study parameters
  if (outcome=="binary") {
    y.m1 = glm(fmar,family=poisson(link="log"),data=subset(data,data[,Iname]==0))
    y.m2 = glm(fcon,family=poisson(link="log"),data=subset(data,data[,Iname]==0))
  } else {
    y.m1 = lm(fmar,data=subset(data,data[,Iname]==0))
    y.m2 = lm(fcon,data=subset(data,data[,Iname]==0))
  }
  ab.m = c(y.m1$coefficients[2],y.m2$coefficients[2])
  
  # internal validation study parameters
  if (outcome=="binary") {
    if (sum(subset(data,data[,Iname]==1)[,Yname])>0) {
      y.i1 = glm(fmarreal,family=poisson(link="log"),data=subset(data,data[,Iname]==1))
      y.i2 = glm(fconreal,family=poisson(link="log"),data=subset(data,data[,Iname]==1))
      ab.i = c(y.i1$coefficients[2],y.i2$coefficients[2])
    }
  } else {
    y.i1 = lm(fmarreal,data=subset(data,data[,Iname]==1))
    y.i2 = lm(fconreal,data=subset(data,data[,Iname]==1))
    ab.i = c(y.i1$coefficients[2],y.i2$coefficients[2])
  }
  # variance of measurement error paramaters
  Uc_f=function(c=c(1,2),data=mydata,outcome="binary") {
    data = subset(data,data[,Iname]==1)
    n=dim(data)[1]
    out = array(0,dim=c(2+nz,n))
    Zterm=`if`(nz==0,0,c(as.matrix(data[,Zname,drop=F]) %*% c[-c(1:2)])  )
    eps = data[,Xname] - as.vector(cbind(1,data[,Xsname]) %*% c[1:2]) - Zterm 
    out[1,]=1*eps
    out[2,]=1*data[,Xsname]*eps
    if (nz>0) {
      for (m in (1:nz)) {
        out[2+m,]=data[,Zname[m]]*eps
      }
    }
    out
  }
  
  dUc_f = function(c,data=mydata,outcome="binary") {
    Uc_f2 = function(c=c(1,2),k=1) {
      apply(Uc_f(c=c,data=data,outcome=outcome),1,sum)[k]
    }
    out=matrix(0,ncol=2+nz,nrow=2+nz)
    out[,1]=numDeriv::grad(Uc_f2,x=c,k=1)
    out[,2]=numDeriv::grad(Uc_f2,x=c,k=2)
    if (nz>0) {
      for (m in (1:nz)) {
        out[,2+m]=numDeriv::grad(Uc_f2,x=c,k=2+m)
      }
    }
    rbind(out)
  }
  
  Ud_f=function(d=c(1,2,3),data=mydata) {
    data = subset(data,data[,Iname]==1)
    n=dim(data)[1]
    out = array(0,dim=c(3+nz,n))
    Zterm=`if`(nz==0,0,c(as.matrix(data[,Zname,drop=F]) %*% d[-c(1:3)])  )
    eps = data[,Xname] - as.vector(cbind(1,data[,Xsname],data[,Mname]) %*% d[1:3]) -Zterm
    out[1,]=1*eps
    out[2,]=data[,Xsname]*eps
    out[3,]=data[,Mname]*eps
    if (nz>0) {
      for (m in (1:nz)) {
        out[3+m,]=data[,Zname[m]]*eps
      }
    }
    out
  }
  
  dUd_f = function(d,data=mydata) {
    Ud_f2 = function(d=d(1,2,3),k=1) {
      apply(Ud_f(d=d,data=data),1,sum)[k]
    }
    out=matrix(0,ncol=3+nz,nrow=3+nz)
    out[,1]=numDeriv::grad(Ud_f2,x=d,k=1)
    out[,2]=numDeriv::grad(Ud_f2,x=d,k=2)
    out[,3]=numDeriv::grad(Ud_f2,x=d,k=3)
    if (nz>0) {
      for (m in (1:nz)) {
        out[,3+m]=numDeriv::grad(Ud_f2,x=d,k=3+m)
      }
    }
    out
  }
  Uc = Uc_f(c=calib.m1$coefficients,data=data)
  Qc = matrix(0,ncol=2+nz+3+nz,nrow=2+nz)
  Qc[1:(2+nz),1:(2+nz)] = dUc_f(c=calib.m1$coefficients,data=data) 
  
  Ud = Ud_f(d=calib.m2$coefficients,data=data)
  Qd = matrix(0,ncol=2+nz+3+nz,nrow=3+nz)
  Qd[1:(3+nz),(2+nz+1):(2+nz+3+nz)] = dUd_f(d=calib.m2$coefficients,data=data) 
  Q=rbind(Qc,Qd)
  U = rbind(Uc,Ud)
  UUT = U %*% t(U)
  V.cd= (solve(Q) %*% UUT %*% solve(Q))[c(2,2+nz+2),c(2,2+nz+2)]
  # variance of outcome model parameters in the main study
  Ua_f = function(a,data=mydata,outcome="binary") {
    n=dim(data)[1]
    out = array(0,dim=c(2+nw,n))
    Wterm=`if`(nw==0,0,c(as.matrix(data[,Wname,drop=F]) %*% a[-c(1:2)])  )
    if (outcome=="binary") {
      eps = data[,Yname] - exp(as.vector(cbind(1,data[,Xsname]) %*% a[1:2])+Wterm )
    } else {
      eps = data[,Yname] - (as.vector(cbind(1,data[,Xsname]) %*% a[1:2])+Wterm)
    }
    out[1,]= eps
    out[2,]= data[,Xsname]*eps
    if (nw>0) {
      for (m in (1:nw)) {
        out[2+m,]=data[,Wname[m]]*eps
      }
    }
    out
  }
  
  dUa_f = function(a,data=mydata,outcome="binary") {
    Ua_f2 = function(x=c(1,2),k=1) {
      apply(Ua_f(a=x,data=data,outcome=outcome),1,sum)[k]
    }
    out=matrix(0,ncol=2+nw,nrow=2+nw)
    out[,1]=numDeriv::grad(Ua_f2,x=c(a),k=1)
    out[,2]=numDeriv::grad(Ua_f2,x=c(a),k=2)
    if (nw>0) {
      for (m in (1:nw)) {
        out[,2+m]=numDeriv::grad(Ua_f2,x=a,k=2+m)
      }
    }
    out
  }
  
  Ub_f = function(b,data=mydata,outcome="binary") {
    n=dim(data)[1]
    out = array(0,dim=c(3+nw,n))
    Wterm=`if`(nw==0,0,c(as.matrix(data[,Wname,drop=F]) %*% b[-c(1:3)])  )
    if (outcome=="binary") {
      eps = data[,Yname] - exp(as.vector(cbind(1,data[,Xsname],data[,Mname]) %*% b[1:3])+Wterm)
    } else {
      eps = data[,Yname] - (as.vector(cbind(1,data[,Xsname],data[,Mname]) %*% b[1:3])+Wterm)
    }
    out[1,]= eps
    out[2,]= data[,Xsname]*eps
    out[3,]= data[,Mname]*eps
    if (nw>0) {
      for (m in (1:nw)) {
        out[3+m,]=data[,Wname[m]]*eps
      }
    }
    out
  }
  
  dUb_f = function(b,data=mydata,outcome="binary") {
    Ub_f2 = function(x=c(1,2,3),k=1) {
      apply(Ub_f(b=x,data=data,outcome=outcome),1,sum)[k]
    }
    out=matrix(0,ncol=3+nw,nrow=3+nw)
    out[,1]=numDeriv::grad(Ub_f2,x=c(b),k=1)
    out[,2]=numDeriv::grad(Ub_f2,x=c(b),k=2)
    out[,3]=numDeriv::grad(Ub_f2,x=c(b),k=3)
    if (nw>0) {
      for (m in (1:nw)) {
        out[,3+m]=numDeriv::grad(Ub_f2,x=b,k=3+m)
      }
    }
    out
  }
  Ua = Ua_f(a=y.m1$coefficients,data=subset(data,data[,Iname]==0),outcome=outcome)
  Qa = matrix(0,ncol=2+nw+3+nw,nrow=2+nw)
  Qa.tmp = dUa_f(a=y.m1$coefficients,data=subset(data,data[,Iname]==0),outcome=outcome)
  Qa[1:(2+nw),1:(2+nw)] = Qa.tmp[1:(2+nw),1:(2+nw)]
  
  Ub = Ub_f(b=y.m2$coefficients,data=subset(data,data[,Iname]==0),outcome=outcome)
  Qb = matrix(0,ncol=2+nw+3+nw,nrow=3+nw)
  Qb.tmp = dUb_f(b=y.m2$coefficients,data=subset(data,data[,Iname]==0),outcome=outcome)
  Qb[1:(3+nw),(2+nw+1):(2+nw+3+nw)] = Qb.tmp[1:(3+nw),1:(3+nw)]
  Q=rbind(Qa,Qb)
  U = rbind(Ua,Ub)
  UUT = U %*% t(U)
  V.abm= (solve(Q) %*% UUT %*% solve(Q))[c(2,2+nw+2),c(2,2+nw+2)]
  # delta method
  XX = matrix(0,ncol=4,nrow=4); XX[1:2,1:2]=V.abm; XX[3:4,3:4] = V.cd
  V11 = 1/(cd[1]^2) * V.abm[1,1] + ab.m[1]^2/(cd[1]^4) * V.cd[1,1]
  V22 = 1/(cd[2]^2) * V.abm[2,2] + ab.m[2]^2/(cd[2]^4) * V.cd[2,2]
  V12 = as.vector(t(c(1/cd[1],0,-ab.m[1]/(cd[1]^2),0)) %*% XX %*% c(0,1/cd[2],0,-ab.m[2]/(cd[2]^2)))
  V.abm=matrix(c(V11,V12,V12,V22),ncol=2,byrow=T)
  ab.m = ab.m/cd
  # variance of outcome model parameters in the validation study
  if (outcome=="binary") {
    if (sum(subset(data,data[,Iname]==1)[,Yname])>0) {
      subdat = subset(data,data[,Iname]==1); subdat[,Xsname]=subdat[,Xname]
      Ua = Ua_f(a=y.i1$coefficients,data=subdat,outcome=outcome)
      Qa = matrix(0,ncol=2+nw+3+nw,nrow=2+nw)
      Qa.tmp = dUa_f(a=y.i1$coefficients,data=subdat,outcome=outcome)
      Qa[1:(2+nw),1:(2+nw)] = Qa.tmp[1:(2+nw),1:(2+nw)]
      Ub = Ub_f(b=y.i2$coefficients,data=subdat,outcome=outcome)
      Qb = matrix(0,ncol=2+nw+3+nw,nrow=3+nw)
      Qb.tmp = dUb_f(b=y.i2$coefficients,data=subdat,outcome=outcome)
      Qb[1:(3+nw),(2+nw+1):(2+nw+3+nw)] = Qb.tmp[1:(3+nw),1:(3+nw)]
      Q=rbind(Qa,Qb)
      U = rbind(Ua,Ub)
      UUT = U %*% t(U)
      V.abi= (solve(Q) %*% UUT %*% solve(Q))[c(2,2+nw+2),c(2,2+nw+2)]
      ab.com = as.vector(solve(solve(V.abm)+solve(V.abi))  %*%  (solve(V.abm) %*% ab.m + solve(V.abi) %*% ab.i))
      V.abcom = solve(solve(V.abm)+solve(V.abi))
    } else {
      ab.com = ab.m
      V.abcom = V.abm
    }
  } else {
    subdat = subset(data,data[,Iname]==1); subdat[,Xsname]=subdat[,Xname]
    Ua = Ua_f(a=y.i1$coefficients,data=subdat,outcome=outcome)
    Qa = matrix(0,ncol=2+nw+3+nw,nrow=2+nw)
    Qa.tmp = dUa_f(a=y.i1$coefficients,data=subdat,outcome=outcome)
    Qa[1:(2+nw),1:(2+nw)] = Qa.tmp[1:(2+nw),1:(2+nw)]
    Ub = Ub_f(b=y.i2$coefficients,data=subdat,outcome=outcome)
    Qb = matrix(0,ncol=2+nw+3+nw,nrow=3+nw)
    Qb.tmp = dUb_f(b=y.i2$coefficients,data=subdat,outcome=outcome)
    Qb[1:(3+nw),(2+nw+1):(2+nw+3+nw)] = Qb.tmp[1:(3+nw),1:(3+nw)]
    Q=rbind(Qa,Qb)
    U = rbind(Ua,Ub)
    UUT = U %*% t(U)
    V.abi= (solve(Q) %*% UUT %*% solve(Q))[c(2,2+nw+2),c(2,2+nw+2)]
    ab.com = as.vector(solve(solve(V.abm)+solve(V.abi))  %*%  (solve(V.abm) %*% ab.m + solve(V.abi) %*% ab.i))
    V.abcom = solve(solve(V.abm)+solve(V.abi))
  }
  
  TE = ab.com[1]
  NIE= ab.com[1]-ab.com[2]
  MP = 1-ab.com[2]/ab.com[1]
  # SE of TE
  se.TE = sqrt(V.abcom[1,1])
  # SE of MP
  a1=ab.com[1]
  b1=ab.com[2]
  se.MP=sqrt(as.numeric(t(c(b1/(a1^2),-1/a1)) %*% V.abcom %*% c(b1/(a1^2),-1/a1)))
  # SE of NIE
  se.NIE=sqrt(as.numeric(t(c(1,-1)) %*% V.abcom %*% c(1,-1)))
  
  c(TE=TE,NIE=NIE,MP=MP,se.TE=se.TE,se.MP=se.MP,se.NIE=se.NIE)
}

# Functions for the naive estimator
get_naive=function(data=mydata,Xname="X",Xsname="Xs",Yname="Y",Mname="M",Wname="W",Zname="W",Iname="I",type="true",outcome="binary") {
  #data=mydata
  #Xname="X";Xsname="Xs";Yname="Y";Mname="M";Wname="W";Zname="W";Iname="I"
  #outcome="continuous"
  #type="true"
  nw=length(Wname);nz=length(Zname)
  
  if (type=="true") {
    data[,Xsname]=data[,Xname]
  }
  
  fcon=as.formula(paste(Yname,"~",Xsname,"+",Mname,ifelse(nw==0,"",paste("+",paste(Wname,collapse="+")))))
  fmar=as.formula(paste(Yname,"~",Xsname,ifelse(nw==0,"",paste("+",paste(Wname,collapse="+")))))
  
  #data=mydata
  if (outcome=="binary") {
    y.m1 = glm(fmar,family=poisson(link="log"),data=data)
    y.m2 = glm(fcon,family=poisson(link="log"),data=data)
  } else {
    y.m1 = lm(fmar,data=data)
    y.m2 = lm(fcon,data=data)
  }
  
  TE = y.m1$coefficients[2]
  NIE = y.m1$coefficients[2] - y.m2$coefficients[2]
  MP = NIE/TE
  
  Ua_f = function(a,data=mydata,outcome="binary") {
    n=dim(data)[1]
    out = array(0,dim=c(2+nw,n))
    Wterm=`if`(nw==0,0,c(as.matrix(data[,Wname,drop=F]) %*% a[-c(1:2)])  )
    if (outcome=="binary") {
      eps = data[,Yname] - exp(as.vector(cbind(1,data[,Xsname]) %*% a[1:2])+Wterm )
    } else {
      eps = data[,Yname] - (as.vector(cbind(1,data[,Xsname]) %*% a[1:2])+Wterm)
    }
    out[1,]= eps
    out[2,]= data[,Xsname]*eps
    if (nw>0) {
      for (m in (1:nw)) {
        out[2+m,]=data[,Wname[m]]*eps
      }
    }
    out
  }
  
  dUa_f = function(a,data=mydata,outcome="binary") {
    Ua_f2 = function(x=c(1,2),k=1) {
      apply(Ua_f(a=x,data=data,outcome=outcome),1,sum)[k]
    }
    out=matrix(0,ncol=2+nw,nrow=2+nw)
    out[,1]=numDeriv::grad(Ua_f2,x=c(a),k=1)
    out[,2]=numDeriv::grad(Ua_f2,x=c(a),k=2)
    if (nw>0) {
      for (m in (1:nw)) {
        out[,2+m]=numDeriv::grad(Ua_f2,x=a,k=2+m)
      }
    }
    out
  }
  
  Ub_f = function(b,data=mydata,outcome="binary") {
    n=dim(data)[1]
    out = array(0,dim=c(3+nw,n))
    Wterm=`if`(nw==0,0,c(as.matrix(data[,Wname,drop=F]) %*% b[-c(1:3)])  )
    if (outcome=="binary") {
      eps = data[,Yname] - exp(as.vector(cbind(1,data[,Xsname],data[,Mname]) %*% b[1:3])+Wterm)
    } else {
      eps = data[,Yname] - (as.vector(cbind(1,data[,Xsname],data[,Mname]) %*% b[1:3])+Wterm)
    }
    out[1,]= eps
    out[2,]= data[,Xsname]*eps
    out[3,]= data[,Mname]*eps
    if (nw>0) {
      for (m in (1:nw)) {
        out[3+m,]=data[,Wname[m]]*eps
      }
    }
    out
  }
  
  dUb_f = function(b,data=mydata,outcome="binary") {
    Ub_f2 = function(x=c(1,2,3),k=1) {
      apply(Ub_f(b=x,data=data,outcome=outcome),1,sum)[k]
    }
    out=matrix(0,ncol=3+nw,nrow=3+nw)
    out[,1]=numDeriv::grad(Ub_f2,x=c(b),k=1)
    out[,2]=numDeriv::grad(Ub_f2,x=c(b),k=2)
    out[,3]=numDeriv::grad(Ub_f2,x=c(b),k=3)
    if (nw>0) {
      for (m in (1:nw)) {
        out[,3+m]=numDeriv::grad(Ub_f2,x=b,k=3+m)
      }
    }
    out
  }
  
  # parameter a
  Ua = Ua_f(a=y.m1$coefficients,data=data,outcome=outcome)
  Qa = matrix(0,ncol=2+nw+3+nw,nrow=2+nw)
  Qa.tmp = dUa_f(a=y.m1$coefficients,data=data,outcome=outcome)
  Qa[1:(2+nw),1:(2+nw)] = Qa.tmp[1:(2+nw),1:(2+nw)]
  # parameter b
  Ub = Ub_f(b=y.m2$coefficients,data=data,outcome=outcome)
  Qb = matrix(0,ncol=2+nw+3+nw,nrow=3+nw)
  Qb.tmp = dUb_f(b=y.m2$coefficients,data=data,outcome=outcome)
  Qb[1:(3+nw),(2+nw+1):(2+nw+3+nw)] = Qb.tmp[1:(3+nw),1:(3+nw)]
  # Q
  Q=rbind(Qa,Qb)
  # UUT
  U = rbind(Ua,Ub)
  UUT = U %*% t(U)
  # vcov
  V.par= solve(Q) %*% UUT %*% solve(Q)
  aW=`if`(nw==0,NULL,paste("aW",1:nw,sep=""))
  bW=`if`(nw==0,NULL,paste("bW",1:nw,sep=""))
  colnames(V.par) = rownames(V.par) =c("a0","a1",aW,"b0","b1","b2",bW)
  # SE of TE
  se.TE = sqrt(V.par["a1","a1"])
  # SE of MP
  V.sub = V.par[c("a1","b1"),c("a1","b1")]
  a1=y.m1$coefficients[2]
  b1=y.m2$coefficients[2]
  se.MP=sqrt(as.numeric(t(c(b1/(a1^2),-1/a1)) %*% V.sub %*% c(b1/(a1^2),-1/a1)))
  # SE of NIE
  se.NIE=sqrt(as.numeric(t(c(1,-1)) %*% V.sub %*% c(1,-1)))
  
  c(TE=TE,NIE=NIE,MP=MP,se.TE=se.TE,se.MP=se.MP,se.NIE=se.NIE)
}

#############################################################################################
#   
# SECTION 2: ILLUSTRATIVE EXAMPLE
#
#############################################################################################

#######################
# Step 1: generate data
#######################
set.seed(2022)
mydata=gen_data(n=1000,n2=200,rho=0.7,cormx=0.4,TE=log(1.5),MP=0.5,PY=0.5,outcome="binary")
# an overview of the dataset. Here, X is the exposure, M is mediator, Y is outcome,
# Xs is the surrogate exposure, W is the confounder, I is the indicator for validation study
# (I=0 for main study; I=0 for validation study)
#   X           M  Y         Xs          W I
# 1 NA -0.07861102 0 -0.2002656 -1.4038913 0
# 2 NA  0.22702352 0  0.2734961  1.3763708 0
# 3 NA  0.41977096 1  0.4318740  0.8339203 0
# 4 NA  0.27585670 0  2.4026005  0.8805160 0
# 5 NA -0.39751675 1 -0.2560545  0.6151019 0
# 6 NA  1.85265051 1  3.3037772  2.4865786 0

#######################
# Step 2: Full EEE
#######################

get_eee1(data=mydata,Xname="X",Xsname="Xs",Yname="Y",Mname="M",Wname="W",Zname="W",Iname="I")


#######################
# Step 3: Reduced EEE
#######################

get_eee2(data=mydata,Xname="X",Xsname="Xs",Yname="Y",Mname="M",Wname="W",Zname="W",Iname="I")


################################
# Step 4: Regression Calibration
################################

get_rc(data=mydata,Xname="X",Xsname="Xs",Yname="Y",Mname="M",Wname="W",Zname="W",Iname="I")

################################
# Step 5: Naive Estimation
################################
get_naive(data=mydata,Xname="X",Xsname="Xs",Yname="Y",Mname="M",Wname="W",Zname="W",Iname="I",
          type="naive",outcome="binary")
