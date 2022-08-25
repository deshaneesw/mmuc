# mmuc
An R toolkit for mediation analysis with multiple mediators under unmeasured mediator-outcome confounding.

# How to install
This package can be installed using devtools package:
```
devtools::install_github("deshaneesw/mmuc")
```
and then loaded into the R enviorenment using:
```
library("mmuc")
```

# Overview
The main function in this package `muc` performs mediation analysis with multiple mediators under unmeasured mediator-outcome confounding. It reports the estimates of Natural Direct Effect (NDE) and Natural Indirect Effect (NIE) together with the corresponding standard errors and confidence intervals. It can handle both binary and continuous type exposure variables. 
```
#y          :a character string indicating the name of the outcome variable, a continuous variable.
#a          :a character string indicating the name of the exposure, a binary or continuous variable.
#m          :a vector of character strings indicating the names of the mediators.
#cv         :a vector of character strings indicating the names of the covariates.
#Data       :a data frame containing all the variables. If data frame name is missing, the variables are taken from the enclosing environment.
#conf.level :level of confidence for confidence intervals.
#a.link     :a suitable link function which transforms the expectation of the exposure to the linear predictor. This accepts the links "logit", "identity", "log" and "inverse".

muc(y,a,m,cv,Data,conf.level=0.95,a.link="logit")
```

Two different layouts of the results from the `muc` function can be obtained by applying the `print` or `summary` function on the objects belonging to the class `muc`.
```
fit<-muc(y,a,m,cv,Data,conf.level=0.95,a.link="logit")
summary(fit)
print(fit)
```

# Example
```
#A simulation with one binary exposure, one continuous covariate, one continuous outcome and three continuous mediators.
library("mmuc")

#The following package is required to simulate data, but not a requirement to use the `mmuc` package
library("MASS")

#Parameter specification for the simulation run
n=400
iter=1000
eta=0.5
delta=5
t.beta_1<-c(1.2,1.5,1.1,eta)
t.beta_2<-c(0.5,1.2,1.8,eta)
t.beta_3<-c(1.3,1,0.5,eta)
t.theta<-c(1.3,2.5,1.5,eta)
t.theta_2<-c(1.2,0.8,1)


#Actual NDE and NIE values
t.theta[2]  #2.5
c(t.beta_1[2],t.beta_2[2],t.beta_3[2])%*%t.theta_2  #3.76

#Data simulation
set.seed(234)
X<-rnorm(n,0,1)
expit<- function(x) {1/(1+exp(-x))}
A<-rbinom(n,1,expit(as.vector(cbind(1,X)%*%c(0.8,1.2))))
U<-rnorm(n,1+0.5*X,1)
M<-c()
for (j in 1:n) {
  a<-A[j]
  R<-(1+delta*a)*(matrix(c(1,2^(-1),2^(-2),2^(-1),1,2^(-1),2^(-2),2^(-1),1),nrow = 3, ncol = 3))
  mu<-c(cbind(1,a,X[j],U[j])%*%t.beta_1,cbind(1,a,X[j],U[j])%*%t.beta_2,cbind(1,a,X[j],U[j])%*%t.beta_3)
  m<-mvrnorm(n=1,mu,Sigma = R)
  M<-rbind(M,m)
}
colnames(M)<-c("M1","M2","M3")
Y<-rnorm(n,cbind(1,A,X,U)%*%t.theta+M%*%t.theta_2,1)
simdata<-cbind(A,X,Y,M)


#Application of the 'muc' function on the simulated data
fit<-muc(y="Y",a="A",m=c("M1","M2","M3"),cv="X",simdata)
print(fit)
summary(fit)

```

# References
Wickramarachchi, D.S., Lim, L.H.M. and Sun, B. (2022). [Mediation Analysis with Multiple Mediators under Unmeasured Mediator-Outcome Confounding](
https://doi.org/10.48550/arXiv.2205.15206). arXiv e-prints. 
