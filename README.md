# metaFlexB
This is a README file of the R package _metaFlexB_. In our paper, we develop new Bayesian procedures for estimating and testing the overall treatment effect and inter-study heterogeneity in random-effects meta-analysis of rare binary events.

## Installation of the package

To install our package, you may execute the following codes:

```{r, eval = FALSE}
install.packages("devtools")
devtools::install_github("chriszhangm/metaFlexB")
library(metaFlexB)
```
For Mac Users who cannot compile the code, please refer [this answer](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/).

## A Example of Using _metaFlexB_

We show a toy example to apply the function `main_draw`, which produces random posterior draws of global parameters of interest, and the function `metatest`, which will report a model selection result using one of criteria {AIC, BIC, DIC, BST}.

### Generate Data

We will generate a dataset containing 50 independent trials with rare background incident rates.
```{r, eval = FALSE}
----------------------------------------------------------
#nct(k): generate k studies(sample size 50~1000).
#Pct(mu,k,theta,tauS,w): probabilities developing rare events given parameters.
##mu: the background incidence rate.
##k: the number of studies.
##theta: the overall treatment effect.
##tauS: the heterogeneity parameter.
##w: the direction paramter.

#Xct(nct,Pct): generate total events for k studies.
----------------------------------------------------------

nct<-function(k)
{nt<-rep(NA,k)
  nc<-rep(NA,k)
  nc<-runif(k,min=50,max=1000)
  nt<- runif(k,min=50,max=1000)         
  nct<-rbind(round(nc),round(nt))
  return(nct)
}
#####################################################################

Pct<-function(mu,k,theta,tauS,w)
{
  epsi1<-rnorm(k,0,sqrt(0.5))  #sigmaS=0.5
  epsi2<-rnorm(k,0,sqrt(tauS))
  pt<-exp(mu+epsi1+(1-w)*(theta+epsi2))/(1+exp(mu+epsi1+(1-w)*(theta+epsi2)))
  pc<-exp(mu+epsi1-w*(theta+epsi2))/(1+exp(mu+epsi1-w*(theta+epsi2)))
  p<-rbind(pc,pt)
  return(p)
}
#####################################################################
Xct = function(nct,Pct){
  k<-length(nct[1,])
  Xt<-rep(NA,k)
  Xc<-rep(NA,k)
  for(i in 1:k)
  {
    Xc[i]<-rbinom(1,nct[1,i],Pct[1,i])
    Xt[i]<-rbinom(1,nct[2,i],Pct[2,i])
  }
  Xct<-rbind(Xc,Xt)
  return(Xct)
}

----------------------------------------------------------
set.seeds(1234)
#generate large size data (k=50)
nnn = nct(k = 50)
ppp = Pct(mu = -5,k = 50,theta = 1,tauS = 0.8,w = 0)
xxx = Xct(nnn,ppp)
```
### FlexB Estimator
```{r, eval = FALSE}
rrr = main_draw(10000,xxx,nnn)
#treatment effect
mean(rrr$theta0[5001:10000]) #1.037
#inter-study heterogeneity
median(rrr$tauS[5001:10000]) #0.925
#background incident rate
mean(rrr$mu0[5001:10000]) #-4.93
#variance of background incident rate
median(rrr$sigmaS[5001:10000]) #0.53
#ratio variablity
table(rrr$w[5001:10000]) #prob(w==0) = 0.99

#computation time (based on 64-bit operating system; i7-7700HQ CPU @2.8GHz 2.8GHz; 8-cores)
system.time(expr=main_draw(10000,xxx,nnn)) #2.13s
```
### Markov Chain Monte Carlo Convergence Diagnostics
### Posterior Densities of All Global Parameters
### Model Selection


## Notes
Alternative way to access the function:

Step1: Download src folder entirely

Step2: Install Rcpp,RcppArmadillo and RcppDist packages.

Step3: Let your R locate to src folder.

Step4: Use the following R code

```{r, eval = FALSE}
library(Rcpp)
SourceCpp('META.cpp')
```
Then, you should see the function named main_draw() in the section of "Environment".
