\name{main_draw}
\alias{main_draw}

\title{Parameter estimation of the flexible random-effects model in Meta-Analysis}
\description{
  Fit a flexible random-effects model, and generate all posterior estimates (FlexB), including the overall treatment effect, the heterogeneity parameter, the background incidence rate, the residual part, and the direction parameter. 
  
}
\usage{
main_draw(M,xct,nct)
}
\arguments{
  \item{M}{M iterations for the Markov chain Monte Carlo algorithm. 4000-9000 iterations are recommanded for the number of studies < 50, and more than 10,000 iterations are suggested for more studies to achieve the convergence.}
  \item{xct}{a 2 by k matrix of total (rare) events, where k is the number of studies. The first row are total events of the control group, and the second row are cases for the treatment group.}
  \item{nct}{a 2 by k matrix of total subjects, where k is the number of studies. The first row are total subjects of the control group, and the second row are subjects for the treatment group.}
}
\value{
  \code{main_draw(M,xct,nct)} return a full list of posterior samples.\cr
  \code{main_draw(M,xct,nct)$theta0} return M posterior samples of the overall treatment effect.
  \code{main_draw(M,xct,nct)$tauS} return M posterior samples of the heterogeneity parameter.
  \code{main_draw(M,xct,nct)$mu0} return M posterior samples of the background incidence rate.
  \code{main_draw(M,xct,nct)$sigmaS} return M posterior samples of the residual part.
  \code{main_draw(M,xct,nct)$w} return M posterior samples of the direction parameter.

}
\details{
Users are recommanded to use M/2 as a burn-in period. For instance, one may use M = 10000 and estimate the overall treatment effect by calculating the posterior mean of \code{main_draw(M,xct,nct)$theta0 [5001:10000]}.\cr\cr
Then, we use posterior mean to estimate the overall treatment effect and background incidence rate, posterior median to estimate the heterogeneity paramter and residual, and posterior mode to estimate the direction parameter.
}
\references{
  See the documentation for Armadillo, and RcppArmadillo, for more details.
}
\examples{
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
rrr = main_draw(10000,xxx,nnn)

#FlexB estimates
mean(rrr$theta0[5001:10000]) #1.037
median(rrr$tauS[5001:10000]) #0.925
mean(rrr$mu0[5001:10000]) #-4.93
median(rrr$sigmaS[5001:10000]) #0.53
table(rrr$w[5001:10000]) #prob(w==0) = 0.99

#computation time (based on 64-bit operating system; i7-7700HQ CPU @2.8GHz 2.8GHz; 8-cores)
system.time(expr=main_draw(10000,xxx,nnn)) #2.13s

}
\author{Ming Zhang <mingz@smu.edu>}

