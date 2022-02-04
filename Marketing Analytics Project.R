install.packages("tidyverse")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("lubridate")
install.packages("moments")
install.packages("ggthemes")
install.packages("bayesm")

library(tidyverse)
library(ggplot2)
library(dplyr)
library(lubridate)
library(moments)
library(ggthemes)
library(bayesm)

data <- list(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385)

result <- array(c(monitor_data_X_list_R, monitor_data_y_list_R), dim = c(385,18,3))
print(result)

N <- nlevels(monitor_data_X_list_R$Respondent)
data <- vector(mode = "list", length =N)
for(i in 1:4620){
  data[[i]]$y <- monitor_data_y_list_R[monitor_data_y_list_R$Respondent==levels(monitor_data_y_list_R$Respondent)[i], "Choice"]
  data[[i]]$x <- as.matrix(monitor_data_X_list_R[monitor_data_X_list_R$Respondent==levels(monitor_data_X_list_R$Respondent)[i], 4:18])
}

##monitor_data <- read_excel("C:\MSM McCombs\Fall Term\Marketing Analytics - I\Class Project")
m1_data <- pull(monitor_data_X_list_R,Alternative)
m2_data <- pull(monitor_data_y_list_R,Choice)
m3_data <- m1_data + m2_data
m3_data

lgtdata[[i]]$y

Data = list(monitor_data_y_list_R[[385]]$Choice,monitor_data_X_list_R[[18480]]$) 

Prior = list(0, 0.01*385, list(1, min(50,0.1*48), 0.8), list(c(0.01, 2),  c(0.001, 3), c(0.1, 4))) 

Mcmc = list((385*12), 1, 100, (2.93/sqrt(4)), 0.1, 200, 20)

rhierMnlDP(Data, Prior, Mcmc)



?nvar


if(nchar(Sys.getenv("monitor_data_y_list_R")) != 0) {R=20000} else {R=10}
set.seed(66)

p = 4                                # num of choice alterns
ncoef = 48  
nlgt = 300                           # num of cross sectional units
nz = 2
Z = matrix(runif(nz*nlgt),ncol=nz)
Z = t(t(Z)-apply(Z,2,mean))          # demean Z
ncomp = 3                            # no of mixture components
Delta=matrix(c(1,0,1,0,1,2),ncol=2)

comps = NULL
comps[[1]] = list(mu=c(0,-1,-2),   rooti=diag(rep(2,3)))
comps[[2]] = list(mu=c(0,-1,-2)*2, rooti=diag(rep(2,3)))
comps[[3]] = list(mu=c(0,-1,-2)*4, rooti=diag(rep(2,3)))
pvec=c(0.4, 0.2, 0.4)

##  simulate from MNL model conditional on X matrix
simmnlwX = function(n,X,beta) {
  k = length(beta)
  Xbeta = X%*%beta
  j = nrow(Xbeta) / n
  Xbeta = matrix(Xbeta, byrow=TRUE, ncol=j)
  Prob = exp(Xbeta)
  iota = c(rep(1,j))
  denom = Prob%*%iota
  Prob = Prob / as.vector(denom)
  y = vector("double", n)
  ind = 1:j
  for (i in 1:n) {
    yvec = rmultinom(1, 1, Prob[i,])
    y[i] = ind%*%yvec}
  return(list(y=y, X=X, beta=beta, prob=Prob))
}

## simulate data with a mixture of 3 normals
simlgtdata = NULL
ni = rep(50,300)
for (i in 1:nlgt) {
  betai = Delta%*%Z[i,] + as.vector(rmixture(1,pvec,comps)$x)
  Xa = matrix(runif(ni[i]*p,min=-1.5,max=0), ncol=p)
  X = createX(p, na=1, nd=NULL, Xa=Xa, Xd=NULL, base=1)
  outa = simmnlwX(ni[i], X, betai)
  simlgtdata[[i]] = list(y=outa$y, X=X, beta=betai)
}

## plot betas
if(0){
  bmat = matrix(0, nlgt, ncoef)
  for(i in 1:nlgt) { bmat[i,] = simlgtdata[[i]]$beta }
  par(mfrow = c(ncoef,1))
  for(i in 1:ncoef) { hist(bmat[,i], breaks=30, col="magenta") }
}

## set Data and Mcmc lists
keep = 5
Mcmc1 = list(R=R, keep=keep)
Data1 = list(p=p, lgtdata=simlgtdata, Z=Z)

out = rhierMnlDP(Data=Data1, Mcmc=Mcmc1)

cat("Summary of Delta draws", fill=TRUE)
summary(out$Deltadraw, tvalues=as.vector(Delta))

## plotting examples
if(0) {
  plot(out$betadraw)
  plot(out$nmix)
}

