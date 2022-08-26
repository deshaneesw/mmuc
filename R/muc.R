#' @title Mediation analysis with multiple mediators and covariates under unmeasured mediator-outcome confounding
#'
#' @description Mediation analysis with multiple mediators and covariates under unmeasured mediator-outcome confounding
#'
#' @param y a character string indicating the name of the outcome variable.
#' @param a a character string indicating the name of the exposure, a binary or continuous variable.
#' @param m a vector of character strings indicating the names of the mediators.
#' @param cv a vector of character strings indicating the names of the covariates.
#' @param Data a data frame containing all the variables. If data frame name is missing, the variables are taken from the enclosing environment.
#' @param conf.level level of confidence for confidence intervals.
#' @param a.link a suitable link function which transforms the expectation of the exposure to the linear predictor. This accepts the links "logit", "identity", "log" and "inverse".
#'
#' @details This function is for estimating mediation effects when there are multiple mediators and covariates with unmeasured mediator-outcome confounding.
#' It reports the estimates of Natural Direct Effect (NDE) and Natural Indirect Effect (NIE) together with the corresponding standard errors and confidence intervals.
#' It can handle both binary and continuous type exposure variables.
#'
#' @return `muc` returns a list containing the following components
#' \item{nde}{point estimate of the natural direct effect of the exposure on the outcome}
#' \item{nie}{point estimate of the natural indirect effect of the exposure on the outcome via both mediators}
#' \item{se.nde}{Standard error of nde}
#' \item{se.nie}{Standard error of nie}
#' \item{conf.level}{the level of confidence used}
#' \item{y}{the outcome variable}
#' \item{cv}{the covariate}
#' \item{a}{the exposure variable}
#' \item{m}{the mediators}
#' \item{n0}{the sample size before removing observations with missing values}
#' \item{n}{the sample size after removing observations with missing values}
#' \item{nmed}{the number of mediators}
#' \item{nc}{the number of covariates}
#' \item{conv}{State of convergence}
#'
#' @references
#' Wickramarachchi, D.S., Lim, L.H.M. and Sun, B. (2022). \href{https://doi.org/10.48550/arXiv.2205.15206}{Mediation Analysis with Multiple Mediators under Unmeasured Mediator-Outcome Confounding}. arXiv e-prints.
#'
#' @seealso \code{\link{summary.muc}}, \code{\link{print.muc}}
#'
#' @examples
#' #A simulation with one binary exposure, one continuous covariate,
#' #one continuous outcome and three continuous mediators.
#' library("mmuc")
#'
#' #The following package is needed to simulate data and not a requirement to use the `mmuc` package
#' library("MASS")
#'
#' #Parameter specification for the simulation run
#' n=400
#' iter=1000
#' eta=0.5
#' delta=5
#' t.beta_1<-c(1.2,1.5,1.1,eta)
#' t.beta_2<-c(0.5,1.2,1.8,eta)
#' t.beta_3<-c(1.3,1,0.5,eta)
#' t.theta<-c(1.3,2.5,1.5,eta)
#' t.theta_2<-c(1.2,0.8,1)
#'
#'
#' #Actual NDE and NIE values
#' t.theta[2]  #2.5
#' c(t.beta_1[2],t.beta_2[2],t.beta_3[2])%*%t.theta_2  #3.76
#'
#' #Data simulation
#' set.seed(234)
#' X<-rnorm(n,0,1)
#' expit<- function(x) {1/(1+exp(-x))}
#' A<-rbinom(n,1,expit(as.vector(cbind(1,X)%*%c(0.8,1.2))))
#' U<-rnorm(n,1+0.5*X,1)
#' M<-c()
#' for (j in 1:n) {
#'   a<-A[j]
#'   R<-(1+delta*a)*(matrix(c(1,2^(-1),2^(-2),2^(-1),1,2^(-1),2^(-2),2^(-1),1),nrow = 3, ncol = 3))
#'   mu<-c(cbind(1,a,X[j],U[j])%*%t.beta_1,cbind(1,a,X[j],U[j])%*%t.beta_2,cbind(1,a,X[j],U[j])%*%t.beta_3)
#'   m<-mvrnorm(n=1,mu,Sigma = R)
#'   M<-rbind(M,m)
#' }
#' colnames(M)<-c("M1","M2","M3")
#' Y<-rnorm(n,cbind(1,A,X,U)%*%t.theta+M%*%t.theta_2,1)
#' simdata<-cbind(A,X,Y,M)
#'
#'
#' #Application of the 'muc' function on the simulated data
#' fit<-muc(y="Y",a="A",m=c("M1","M2","M3"),cv="X",simdata)
#' print(fit)
#' summary(fit)
#'
#' @export
#'
#'
muc <- function(y,a,m,cv,Data,conf.level=0.95,a.link="logit") {
  if(missing(Data)) {
    if(missing(cv)){
      Data<-data.frame(get(y),get(a),mget(m,inherits = T))
      names(Data)<-c(y,a,m)
      } else {
    Data<-data.frame(get(y),get(a),mget(m,inherits = T),mget(cv,inherits = T))
    names(Data)<-c(y,a,m,cv)
      }
    }
  if(missing(cv)){
    var<-c(y,a,m)
    nc<-0
  } else {
    var<-c(y,a,m,cv)
    nc<-length(cv)
  }
  Dt <- stats::na.omit(Data[,var])
  n0<-nrow(Data)
  n<-nrow(Dt)
  nmed<-length(m)
  if(missing(cv)) {
    C<-as.matrix(rep(1,n))
  } else C<-cbind(1,as.matrix(Dt[,cv]))

  ## estimating functions
  inv.link<-function(x,link){
    switch(link, logit=1/(1+exp(-x)), identity=x, log=exp(x),
            inverse=x^(-1), neginverse=-(x^(-1)), sqroot=x^2, invsq=x^(-1/2),
           probit=, loglog=exp(-exp(-x)), cloglog=1-exp(-exp(x)) )
  }

  ng<-1+nc+(2+nc)*nmed+1+nmed
  est<-function(g){
    h<-rep(0,ng);
    pi<-inv.link(as.vector(C%*%g[1:(1+nc)]),link=a.link)
    beta<-g[(2+nc):(1+nc+(2+nc)*nmed)]
    theta_1<-g[1+nc+(2+nc)*nmed+1]
    theta_2<-g[(nc+(2+nc)*nmed+3):(nc+(2+nc)*nmed+2+nmed)]
    h[1:(1+nc)]<-t(C)%*%(Dt[,a]-pi)
    h[(2+nc):(1+nc+(2+nc)*nmed)]<-as.vector(t(cbind(Dt[,a],C))%*%as.matrix(Dt[,m]-cbind(Dt[,a],C)%*%matrix(beta,ncol=nmed)))
    h[1+nc+(2+nc)*nmed+1]<-t(Dt[,a]-pi)%*%(Dt[,y]-Dt[,a]*theta_1-as.matrix(Dt[,m])%*%theta_2)
    h[(2+nc+(2+nc)*nmed+1):(2+nc+(2+nc)*nmed+nmed)]<- t(Dt[,a]-pi)%*%as.matrix((Dt[,m]-cbind(Dt[,a],C)%*%matrix(beta,ncol=nmed))*as.vector(Dt[,y]-as.matrix(Dt[,m])%*%theta_2))
    h
  }

  est.fun<-function(g){
    pi<-inv.link(as.vector(C%*%g[1:(1+nc)]),link=a.link)
    beta<-g[(2+nc):(1+nc+(2+nc)*nmed)]
    theta_1<-g[1+nc+(2+nc)*nmed+1]
    theta_2<-g[(nc+(2+nc)*nmed+3):(nc+(2+nc)*nmed+2+nmed)]
    f<-C*(Dt[,a]-pi)
    for (i in 1:nmed) {
      f<-cbind(f,cbind(Dt[,a],C)*as.vector(Dt[,m[i]]-cbind(Dt[,a],C)%*%(matrix(beta,ncol=nmed)[,i])))
    }
    f<-cbind(f,(Dt[,a]-pi)*(Dt[,y]-Dt[,a]*theta_1-as.matrix(Dt[,m])%*%theta_2))
    for (i in 1:nmed) {
      f<-cbind(f,(Dt[,a]-pi)*((Dt[,m[i]]-cbind(Dt[,a],C)%*%(matrix(beta,ncol=nmed)[,i]))*(Dt[,y]-as.matrix(Dt[,m])%*%theta_2)))
    }
    f
  }

  #Proposed method
  run<-BB::BBsolve(par = rep(0,ng), fn = est,quiet=T,control=list(maxit=100000,noimp=50000))
  run$residual
  conv<-run$message
  run$fn.reduction
  sol<-run[1]
  nde<-sol$par[1+nc+(2+nc)*nmed+1]
  nie.m<-sol$par[seq(2+nc,(1+nc+(2+nc)*nmed),2+nc)]*sol$par[seq(3+nc+(2+nc)*nmed,1+nc+(2+nc)*nmed+1+nmed,1)]
  nie<-sum(nie.m)

  G<-numDeriv::jacobian(func=est,x=sol$par)/n
  om<-est.fun(sol$par)
  gmm.var<-solve(G)%*%(t(om)%*%om)%*%t(solve(G))/n
  var_nde<-diag(solve(G)%*%(t(om)%*%om/n)%*%t(solve(G))/n)[1+nc+(2+nc)*nmed+1]
  se.nde<-sqrt(var_nde)

  #variance of nie
  h<-rep(0,1+nc)
  for (i in seq(3+nc+(2+nc)*nmed,1+nc+(2+nc)*nmed+1+nmed,1)) {
    h<-append(h,c(sol$par[i],rep(0,1+nc)))
  }
  h<-append(h,0)
  for (i in seq(2+nc,(1+nc+(2+nc)*nmed),2+nc)) {
    h<-append(h,sol$par[i])
  }
  var_nie<-t(h)%*%(solve(G)%*%(t(om)%*%om/n)%*%t(solve(G))/n)%*%t(t(h))
  se.nie<-sqrt(var_nie)
  if(missing(cv)){
    out<-list(nde=nde,nie=nie,se.nde=se.nde,se.nie=se.nie,
              conf.level=conf.level,y=y,cv="None",a=a,m=m,n0=n0,n=n,nmed=nmed,nc=nc,conv=conv)

  } else {out<-list(nde=nde,nie=nie,se.nde=se.nde,se.nie=se.nie,
                   conf.level=conf.level,y=y,cv=cv,a=a,m=m,n0=n0,n=n,nmed=nmed,nc=nc,conv=conv)}

  class(out)<-"muc"
  out
}


## Summary

#' Summarizing output from \code{\link{muc}}
#'
#' This function is a method for class \code{\link{muc}} objects.
#'
#' @aliases summary.muc print.muc
#'
#' @param object an object of class `\code{muc}`.
#' @param x an object of class \code{summary.muc}.
#' @param ...  additional arguments affecting the summary produced.
#'
#' @seealso \code{\link{muc}}
#'
#' @examples
#' #For examples see example(muc)
#'
#' @export
#'
summary.muc<-function(object,...){
  x<-object
  cat("Mediation analysis under unmeasured confounding","\n","\n")
  cat("Outcome:",x$y,"\n","Exposure:",x$a,"\n","Mediator(s):",x$m,"\n","Covariate(s):",x$cv,"\n","\n",sep = " ")
  if(x$n0-x$n>0) cat("(",x$n0-x$n," observations deleted due to missingness)",sep = "")
  cip <- 100 * x$conf.level
  q<-stats::qnorm((1-x$conf.level)/2,lower.tail = F)
  est<-c(x$nde,x$nie)
  se<-c(x$se.nde,x$se.nie)
  ci_l<-c(x$nde-q*x$se.nde,x$nie-q*x$se.nie)
  ci_u<-c(x$nde+q*x$se.nde,x$nie+q*x$se.nie)
  est<-cbind(est,se,ci_l,ci_u)
  colnames(est)<-c("Estimate","S.E.",paste(cip, "% CI Lower", sep=""),
                   paste(cip, "% CI Upper", sep=""))
  rownames(est)<-c("NDE","NIE")
  stats::printCoefmat(est,digits = 5)
}
#' @rdname summary.muc
#' @export
print.muc<-function(x,...){
  cat("Mediation analysis under unmeasured confounding","\n","\n")
  cip <- 100 * x$conf.level
  q<-stats::qnorm((1-x$conf.level)/2,lower.tail = F)
  est<-c(x$nde,x$nie)
  se<-c(x$se.nde,x$se.nie)
  est<-cbind(est,se)
  colnames(est)<-c("Estimate","S.E.")
  rownames(est)<-c("NDE","NIE")
  stats::printCoefmat(est,digits = 5)
}
