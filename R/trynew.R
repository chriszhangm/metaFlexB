metatest=function(xxx,nnn,criterion='AIC'){
  #data generating process nct, Pct, Xct
  nct<-function(k){
    nt<-rep(NA,k)
    nc<-rep(NA,k)
    nc<-runif(k,min=50,max=1000)
    #R = sample(c(1,1.5,2),1)
    nt<- runif(k,min=50,max=1000)         
    nct<-rbind(round(nc),round(nt))
    return(nct)
  }
  Pct<-function(mu,k,theta,tauS,w){ #####
    
    epsi1<-rnorm(k,0,sqrt(0.5))  #sigmaS=0.5
    epsi2<-rnorm(k,0,sqrt(tauS))
    pt<-exp(mu+epsi1+(1-w)*(theta+epsi2))/(1+exp(mu+epsi1+(1-w)*(theta+epsi2)))
    pc<-exp(mu+epsi1-w*(theta+epsi2))/(1+exp(mu+epsi1-w*(theta+epsi2)))
    p<-rbind(pc,pt)
    return(p)}
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
  if (criterion=='AIC') {
    k=dim(nnn)[2]
    study = rep(1:k,2)
    treat1 = rep(c(0,1),each=k)
    treat2 = rep(c(1,0),each=k)
    treat12 = rep(c(-0.5,0.5),each=k)
    n = c(nnn[2,],nnn[1,])
    event = c(xxx[2,],xxx[1,])
    dff = data.frame(study=factor(study),treat1=treat1,treat2 = treat2,treat12 = treat12, n=n,event=event)
    #m1 #both =0
    ss1 = glmer(cbind(event,n-event)~(1|study),data=dff,family=binomial(link="logit"),
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5),check.conv.grad= .makeCC("warning", tol = 2e-2, relTol = NULL)
                                     ,check.conv.singular = .makeCC(action = "warning", tol = formals(isSingular)$tol)),nAGQ = 1)
    #m2-m4 #theta=0
    #w=0
    ss3 = glmer(cbind(event,n-event)~(1|study)+(treat1-1|study),data=dff,family=binomial(link="logit"),
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5),check.conv.grad= .makeCC("warning", tol = 2e-2, relTol = NULL),check.conv.singular = .makeCC(action = "warning", tol = formals(isSingular)$tol)),nAGQ = 1)
    #w=1
    ss32 = glmer(cbind(event,n-event)~(1|study)+(treat2-1|study),data=dff,family=binomial(link="logit"),
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5),check.conv.grad= .makeCC("warning", tol = 2e-2, relTol = NULL),check.conv.singular = .makeCC(action = "warning", tol = formals(isSingular)$tol)),nAGQ = 1)
    #w=0.5
    ss31 = glmer(cbind(event,n-event)~(1|study)+(treat12-1|study),data=dff,family=binomial(link="logit"),
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5),check.conv.grad= .makeCC("warning", tol = 2e-2, relTol = NULL),check.conv.singular = .makeCC(action = "warning", tol = formals(isSingular)$tol)),nAGQ = 1)
    
    #m5 #tau=0
    ss2 = glmer(cbind(event,n-event)~(1|study)+factor(treat2),data=dff,family=binomial(link="logit"),
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5),check.conv.grad= .makeCC("warning", tol = 2e-2, relTol = NULL),check.conv.singular = .makeCC(action = "warning", tol = formals(isSingular)$tol)),nAGQ = 1)
    
    #m6-m8 #both !=0
    ss4 = glmer(cbind(event,n-event)~(1|study)+(treat1-1|study)+factor(treat1),data=dff,family=binomial(link="logit"),
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5),check.conv.grad= .makeCC("warning", tol = 2e-2, relTol = NULL),check.conv.singular = .makeCC(action = "warning", tol = formals(isSingular)$tol)),nAGQ = 1)
    ss42 = glmer(cbind(event,n-event)~(1|study)+(treat2-1|study)+factor(treat2),data=dff,family=binomial(link="logit"),
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5),check.conv.grad= .makeCC("warning", tol = 2e-2, relTol = NULL),check.conv.singular = .makeCC(action = "warning", tol = formals(isSingular)$tol)),nAGQ = 1)
    ss41 = glmer(cbind(event,n-event)~(1|study)+(treat12-1|study)+factor(treat12),data=dff,family=binomial(link="logit"),
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5),check.conv.grad= .makeCC("warning", tol = 2e-2, relTol = NULL),check.conv.singular = .makeCC(action = "warning", tol = formals(isSingular)$tol)),nAGQ = 1)
    #AIC
    aic_all = c(AIC(ss1),AIC(ss2),AIC(ss3),AIC(ss31),AIC(ss32),AIC(ss4),AIC(ss41),AIC(ss42))
    if (which.min(aic_all)==3|which.min(aic_all)==4|which.min(aic_all)==5) {
      return(cat('Treatment effects: No.\nInter-study heterogeneity of treatment effects: Yes.'))
    }
    else if (which.min(aic_all)==6|which.min(aic_all)==7|which.min(aic_all)==8) {
      return(cat('Treatment effects: Yes.\nInter-study heterogeneity of treatment effects: Yes.'))
    }
    else if(which.min(aic_all)==1){return(cat('Treatment effects: No.\nInter-study heterogeneity of treatment effects: No.'))}
    else{return(cat('Treatment effects: Yes.\nInter-study heterogeneity of treatment effects: No.'))}
  }
  else if (criterion=='BIC') {
    k=dim(nnn)[2]
    study = rep(1:k,2)
    treat1 = rep(c(0,1),each=k)
    treat2 = rep(c(1,0),each=k)
    treat12 = rep(c(-0.5,0.5),each=k)
    n = c(nnn[2,],nnn[1,])
    event = c(xxx[2,],xxx[1,])
    dff = data.frame(study=factor(study),treat1=treat1,treat2 = treat2,treat12 = treat12, n=n,event=event)
    #m1 #both =0
    ss1 = glmer(cbind(event,n-event)~(1|study),data=dff,family=binomial(link="logit"),
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5),check.conv.grad= .makeCC("warning", tol = 2e-2, relTol = NULL)
                                     ,check.conv.singular = .makeCC(action = "warning", tol = formals(isSingular)$tol)),nAGQ = 1)
    #m2-m4 #theta=0
    #w=0
    ss3 = glmer(cbind(event,n-event)~(1|study)+(treat1-1|study),data=dff,family=binomial(link="logit"),
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5),check.conv.grad= .makeCC("warning", tol = 2e-2, relTol = NULL),check.conv.singular = .makeCC(action = "warning", tol = formals(isSingular)$tol)),nAGQ = 1)
    #w=1
    ss32 = glmer(cbind(event,n-event)~(1|study)+(treat2-1|study),data=dff,family=binomial(link="logit"),
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5),check.conv.grad= .makeCC("warning", tol = 2e-2, relTol = NULL),check.conv.singular = .makeCC(action = "warning", tol = formals(isSingular)$tol)),nAGQ = 1)
    #w=0.5
    ss31 = glmer(cbind(event,n-event)~(1|study)+(treat12-1|study),data=dff,family=binomial(link="logit"),
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5),check.conv.grad= .makeCC("warning", tol = 2e-2, relTol = NULL),check.conv.singular = .makeCC(action = "warning", tol = formals(isSingular)$tol)),nAGQ = 1)
    
    #m5 #tau=0
    ss2 = glmer(cbind(event,n-event)~(1|study)+factor(treat2),data=dff,family=binomial(link="logit"),
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5),check.conv.grad= .makeCC("warning", tol = 2e-2, relTol = NULL),check.conv.singular = .makeCC(action = "warning", tol = formals(isSingular)$tol)),nAGQ = 1)
    
    #m6-m8 #both !=0
    ss4 = glmer(cbind(event,n-event)~(1|study)+(treat1-1|study)+factor(treat1),data=dff,family=binomial(link="logit"),
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5),check.conv.grad= .makeCC("warning", tol = 2e-2, relTol = NULL),check.conv.singular = .makeCC(action = "warning", tol = formals(isSingular)$tol)),nAGQ = 1)
    ss42 = glmer(cbind(event,n-event)~(1|study)+(treat2-1|study)+factor(treat2),data=dff,family=binomial(link="logit"),
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5),check.conv.grad= .makeCC("warning", tol = 2e-2, relTol = NULL),check.conv.singular = .makeCC(action = "warning", tol = formals(isSingular)$tol)),nAGQ = 1)
    ss41 = glmer(cbind(event,n-event)~(1|study)+(treat12-1|study)+factor(treat12),data=dff,family=binomial(link="logit"),
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5),check.conv.grad= .makeCC("warning", tol = 2e-2, relTol = NULL),check.conv.singular = .makeCC(action = "warning", tol = formals(isSingular)$tol)),nAGQ = 1)
    #BIC
    bic_all = c(BIC(ss1),BIC(ss2),BIC(ss3),BIC(ss31),BIC(ss32),BIC(ss4),BIC(ss41),BIC(ss42))
    if (which.min(bic_all)==3|which.min(bic_all)==4|which.min(bic_all)==5) {
      return(cat('Treatment effects: No.\nInter-study heterogeneity of treatment effects: Yes.'))
    }
    else if (which.min(bic_all)==6|which.min(bic_all)==7|which.min(bic_all)==8) {
      return(cat('Treatment effects: Yes.\nInter-study heterogeneity of treatment effects: Yes.'))
    }
    else if(which.min(bic_all)==1){return(cat('Treatment effects: No.\nInter-study heterogeneity of treatment effects: No.'))}
    else{return(cat('Treatment effects: Yes.\nInter-study heterogeneity of treatment effects: No.'))}
  }
  else if (criterion=='DIC'){
    model.lm1 = "model{

  # Likelihood
  for (i in 1:n){
    xt[i] ~ dbin(pt[i],nt[i])
    xc[i] ~ dbin(pc[i],nc[i])
    logit(pt[i]) = mu[i] 
    logit(pc[i]) = mu[i]
    mu[i] ~ dnorm(mu0,sigmaSinv)
  }
  
  mu0 ~ dunif(-10,0)
  #Prior on inv.var
  sigmaSinv ~ dgamma(0.01,0.01)
  sigmaS=1/sigmaSinv
}"
    #m2:tau!=0,theta=0 w=0
    model.lm2 = "model{

  # Likelihood
  for (i in 1:n){
    xt[i] ~ dbin(pt[i],nt[i])
    xc[i] ~ dbin(pc[i],nc[i])
    logit(pt[i]) = mu[i] + theta[i]
    logit(pc[i]) = mu[i] 
    mu[i] ~ dnorm(mu0,sigmaSinv)
    theta[i] ~dnorm(0,tauSinv)
  }
  
  #Priors on  beta
  mu0 ~ dunif(-10,0)
  #Prior on inv.var
  sigmaSinv ~ dgamma(0.01,0.01)
  tauSinv ~ dgamma(0.01,0.01)
  sigmaS=1/sigmaSinv
  tauS = 1/tauSinv
}"
    #m3:tau!=0,theta=0 w=0.5
    model.lm3 = "model{

  # Likelihood
  for (i in 1:n){
    xt[i] ~ dbin(pt[i],nt[i])
    xc[i] ~ dbin(pc[i],nc[i])
    logit(pt[i]) = mu[i] + 0.5*theta[i]
    logit(pc[i]) = mu[i] -0.5*theta[i]
    mu[i] ~ dnorm(mu0,sigmaSinv)
    theta[i] ~dnorm(0,tauSinv)
  }
  
  #Priors on  beta
  mu0 ~ dunif(-10,0)
  #Prior on inv.var
  sigmaSinv ~ dgamma(0.01,0.01)
  tauSinv ~ dgamma(0.01,0.01)
  sigmaS=1/sigmaSinv
  tauS = 1/tauSinv
}"
    #m4:tau!=0,theta=0 w=1
    model.lm4 = "model{

  # Likelihood
  for (i in 1:n){
    xt[i] ~ dbin(pt[i],nt[i])
    xc[i] ~ dbin(pc[i],nc[i])
    logit(pt[i]) = mu[i] 
    logit(pc[i]) = mu[i] -theta[i]
    mu[i] ~ dnorm(mu0,sigmaSinv)
    theta[i] ~dnorm(0,tauSinv)
  }
  
  #Priors on  beta
  mu0 ~ dunif(-10,0)
  #Prior on inv.var
  sigmaSinv ~ dgamma(0.01,0.01)
  tauSinv ~ dgamma(0.01,0.01)
  sigmaS=1/sigmaSinv
  tauS = 1/tauSinv
}"
    
    #m5:tau=0,theta!=0 w=0
    model.lm5 = "model{

  # Likelihood
  for (i in 1:n){
    xt[i] ~ dbin(pt[i],nt[i])
    xc[i] ~ dbin(pc[i],nc[i])
    logit(pt[i]) = mu[i] + theta0
    logit(pc[i]) = mu[i] 
    mu[i] ~ dnorm(mu0,sigmaSinv)
  }
  
  #Priors on  beta
  mu0 ~ dunif(-10,0)
  theta0 ~ dunif(-5,5)
  #Prior on inv.var
  sigmaSinv ~ dgamma(0.01,0.01)
  sigmaS=1/sigmaSinv
}"
    #m6:tau=!0,theta!=0 w=0
    model.lm6 = "model{

  # Likelihood
  for (i in 1:n){
    xt[i] ~ dbin(pt[i],nt[i])
    xc[i] ~ dbin(pc[i],nc[i])
    logit(pt[i]) = mu[i] + theta[i]
    logit(pc[i]) = mu[i] 
    mu[i] ~ dnorm(mu0,sigmaSinv)
    theta[i] ~ dnorm(theta0,tauSinv)
  }
  
  #Priors on  beta
  mu0 ~ dunif(-10,0)
  theta0 ~ dunif(-5,5)
  #Prior on inv.var
  sigmaSinv ~ dgamma(0.01,0.01)
  tauSinv ~ dgamma(0.01,0.01)
  sigmaS=1/sigmaSinv
  tauS = 1/tauSinv
}"
    #m7:tau=!0,theta!=0 w=0.5
    model.lm7 = "model{

  # Likelihood
  for (i in 1:n){
    xt[i] ~ dbin(pt[i],nt[i])
    xc[i] ~ dbin(pc[i],nc[i])
    logit(pt[i]) = mu[i] + 0.5*theta[i]
    logit(pc[i]) = mu[i] -0.5*theta[i]
    mu[i] ~ dnorm(mu0,sigmaSinv)
    theta[i] ~ dnorm(theta0,tauSinv)
  }
  
  #Priors on  beta
  mu0 ~ dunif(-10,0)
  theta0 ~ dunif(-5,5)
  #Prior on inv.var
  sigmaSinv ~ dgamma(0.01,0.01)
  tauSinv ~ dgamma(0.01,0.01)
  sigmaS=1/sigmaSinv
  tauS = 1/tauSinv
}"
    #m8:tau=!0,theta!=0 w=1
    model.lm8 = "model{

  # Likelihood
  for (i in 1:n){
    xt[i] ~ dbin(pt[i],nt[i])
    xc[i] ~ dbin(pc[i],nc[i])
    logit(pt[i]) = mu[i]
    logit(pc[i]) = mu[i] -theta[i]
    mu[i] ~ dnorm(mu0,sigmaSinv)
    theta[i] ~ dnorm(theta0,tauSinv)
  }
  
  #Priors on  beta
  mu0 ~ dunif(-10,0)
  theta0 ~ dunif(-5,5)
  #Prior on inv.var
  sigmaSinv ~ dgamma(0.01,0.01)
  tauSinv ~ dgamma(0.01,0.01)
  sigmaS=1/sigmaSinv
  tauS = 1/tauSinv
}"
    xt <- xxx[1,]
    xc <- xxx[2,]
    nt <- nnn[1,]
    nc <- nnn[2,]
    n = dim(nnn)[2]
    
    model1 = jags.model(textConnection(model.lm1),data = list(xt=xt,xc=xc,nc=nc,nt=nt,n=n),n.chains = 2,quiet=TRUE)
    model2 = jags.model(textConnection(model.lm2),data = list(xt=xt,xc=xc,nc=nc,nt=nt,n=n),n.chains = 2,quiet=TRUE)
    model3 = jags.model(textConnection(model.lm3),data = list(xt=xt,xc=xc,nc=nc,nt=nt,n=n),n.chains = 2,quiet=TRUE)
    model4 = jags.model(textConnection(model.lm4),data = list(xt=xt,xc=xc,nc=nc,nt=nt,n=n),n.chains = 2,quiet=TRUE)
    model5 = jags.model(textConnection(model.lm5),data = list(xt=xt,xc=xc,nc=nc,nt=nt,n=n),n.chains = 2,quiet=TRUE)
    model6 = jags.model(textConnection(model.lm6),data = list(xt=xt,xc=xc,nc=nc,nt=nt,n=n),n.chains = 2,quiet=TRUE)
    model7 = jags.model(textConnection(model.lm7),data = list(xt=xt,xc=xc,nc=nc,nt=nt,n=n),n.chains = 2,quiet=TRUE)
    model8 = jags.model(textConnection(model.lm8),data = list(xt=xt,xc=xc,nc=nc,nt=nt,n=n),n.chains = 2,quiet=TRUE)
    
    dic1 = dic.samples(model1,variable.names=c('mu0','sigmaS'),n.iter = 5000,progress.bar=NULL)
    dic2 = dic.samples(model2,variable.names=c('mu0','sigmaS','tauS'),n.iter = 5000,progress.bar=NULL)
    dic3 = dic.samples(model3,variable.names=c('mu0','sigmaS','tauS'),n.iter = 5000,progress.bar=NULL)
    dic4 = dic.samples(model4,variable.names=c('mu0','sigmaS','tauS'),n.iter = 5000,progress.bar=NULL)
    dic5 = dic.samples(model5,variable.names=c('mu0','theta0','sigmaS'),n.iter = 5000,progress.bar=NULL)
    dic6 = dic.samples(model6,variable.names=c('mu0','theta0','sigmaS'),n.iter = 5000,progress.bar=NULL)
    dic7 = dic.samples(model7,variable.names=c('mu0','theta0','sigmaS'),n.iter = 5000,progress.bar=NULL)
    dic8 = dic.samples(model8,variable.names=c('mu0','theta0','sigmaS','tauS'),n.iter = 5000,progress.bar=NULL)
    
    
    dic =c(sum(dic1$deviance) + sum(dic1$penalty),
            sum(dic2$deviance) + sum(dic2$penalty),
            sum(dic3$deviance) + sum(dic3$penalty),
            sum(dic4$deviance) + sum(dic4$penalty),
            sum(dic5$deviance) + sum(dic5$penalty),
            sum(dic6$deviance) + sum(dic6$penalty),
            sum(dic7$deviance) + sum(dic7$penalty),
            sum(dic8$deviance) + sum(dic8$penalty))
    if (which.min(dic)==3|which.min(dic)==4|which.min(dic)==5) {
      return(cat('Treatment effects: No.\nInter-study heterogeneity of treatment effects: Yes.'))
    }
    else if (which.min(dic)==6|which.min(dic)==7|which.min(dic)==8) {
      return(cat('Treatment effects: Yes.\nInter-study heterogeneity of treatment effects: Yes.'))
    }
    else if(which.min(dic)==1){return(cat('Treatment effects: No.\nInter-study heterogeneity of treatment effects: No.'))}
    else{return(cat('Treatment effects: Yes.\nInter-study heterogeneity of treatment effects: No.'))}
  }
  else if (criterion=='BST') {
    #Simple Average Method
    k=dim(nnn)[2]
    mu0<-function(Xc,Xt,nc,nt){
      return(mean(log((Xc+1/2)/(nc-Xc+1/2))))
    }
    mu00<-function(Xc,Xt,nc,nt)
    {
      return(log((Xc+1/2)/(nc-Xc+1/2)))
    }
    SA<-function(Xc,Xt,nc,nt){
      thetaa<-log((Xt+1/2)/(nt-Xt+1/2))-log((Xc+1/2)/(nc-Xc+1/2))
      return(mean(thetaa))
    }
    #IPM Method
    tIPM<-function(Xc,Xt,nc,nt){
      pt<-(Xt+1/2)/(nt+2*1/2)
      pc<-(Xc+1/2)/(nc+2*1/2)
      muhat<-mean(log((Xc+1/2)/(nc-Xc+1/2)))
      thetahat<-log((Xt+1/2)/(nt-Xt+1/2))-log((Xc+1/2)/(nc-Xc+1/2))
      thetasa<-mean(thetahat)
      tauS<-numeric()
      tauS[1]<-0
      stSstar<-1/(nt+1)*(exp(-muhat-thetasa+tauS[1]/2)+2+exp(muhat+thetasa+tauS[1]/2))+1/(nc+1)*(exp(-muhat)+2+exp(muhat))
      W<-1/(stSstar+tauS[1])
      thetawta<-sum(W*thetahat)/sum(W) 
      FtauS<-sum(W*(thetahat-thetawta)^2)-(length(Xc)-1)
      if(FtauS<=0)
      {
        tauSPM=0
      }
      else
      {             
        i<-1
        while(FtauS!=0)
        {
          W<-1/(1/(nt+1)*(exp(-muhat-thetasa+tauS[i]/2)+2+exp(muhat+thetasa+tauS[i]/2))+1/(nc+1)*(exp(-muhat)+2+exp(muhat))
                +tauS[i])
          thetawta<-sum(W*thetahat)/sum(W) 
          FtauS<-sum(W*(thetahat-thetawta)^2)-(length(Xc)-1)
          #print(FtauS)
          #flush.console()
          
          if (abs(FtauS)<=(10^(-13)*2))
          {
            tauSPM<-tauS[i]
            break
          }
          else
          {
            deltatauS<-(sum(W*(thetahat-thetawta)^2)-(length(Xc)-1))/sum(W^2*(thetahat-thetawta)^2)
            tauS[i+1]<-tauS[i]+deltatauS
            i<-i+1
          }
        }
      }
      return(tauSPM)
    } 
    #a function to test the heterogeneity
    fre_test_statistc = function(Xc,Xt,nc,nt){
      thetaa<-log((Xt+1/2)/(nt-Xt+1/2))-log((Xc+1/2)/(nc-Xc+1/2))
      k = length(Xc)
      pt = numeric(k)
      pc = numeric(k)
      #y = numeric(k)
      sa = numeric(k)
      sigma0 = numeric(k)
      for (i in 1:k) {
        pt[i]<-(Xt[i]+1/2)/(nt[i]+1)
        pc[i]<-(Xc[i]+1/2)/(nc[i]+1)
        sa[i] =  log((Xt[i]+1/2)/(nt[i]-Xt[i]+1/2))-log((Xc[i]+1/2)/(nc[i]-Xc[i]+1/2))
        sigma0[i] = 1/(nt[i]*pt[i]*(1-pt[i])) + 1/(nc[i]*pc[i]*(1-pc[i]))
      }
      y = (sa-mean(sa))^2
      B = ((k-2)/k)*sigma0 + sum(sigma0)/k^2
      return(sum(y-B)/sqrt(sum(2*B)))
    }
    #test theta0
    r1=NA
    r2=NA
    pichat = (xxx[1,]+1/2)/(nnn[1,]+1)
    pithat = (xxx[2,]+1/2)/(nnn[2,]+1)
    SASE= 1/(nnn[1,]*pichat*(1-pichat))+1/(nnn[2,]*pithat*(1-pithat))+tIPM(xxx[1,],xxx[2,],nnn[1,],nnn[2,])
    tsa<-k*SA(xxx[1,],xxx[2,],nnn[1,],nnn[2,])/sqrt(sum(SASE))
    if(1-pnorm(tsa)>0.05){r1=0}
    else{r1=1}
    #test tauS
    t4 = fre_test_statistc(xxx[1,],xxx[2,],nnn[1,],nnn[2,])
    #PB
    thetahat=SA(xxx[1,],xxx[2,],nnn[1,],nnn[2,])
    muhat=mu00(xxx[1,],xxx[2,],nnn[1,],nnn[2,])
    pcpb=exp(muhat)/(1+exp(muhat))
    ptpb=exp(muhat+thetahat)/(1+exp(muhat+thetahat))
    t4s = numeric(1000)
    for (j in 1:1000) {
      xctnew = Xct(nnn,Pct = rbind(pcpb,ptpb))
      t4s[j]=fre_test_statistc(xctnew[1,],xctnew[2,],nnn[1,],nnn[2,])
    }
    if (t4 > quantile(t4s,probs = 0.95)[[1]]) {
      r2 = 1
    }
    else{r2=0}
    #BST test
    if (r1==0&r2==1) {
      return(cat('Treatment effects: No.\nInter-study heterogeneity of treatment effects: Yes.'))
    }
    else if (r1==1&r2==1) {
      return(cat('Treatment effects: Yes.\nInter-study heterogeneity of treatment effects: Yes.'))
    }
    else if(r1==0&r2==0){return(cat('Treatment effects: No.\nInter-study heterogeneity of treatment effects: No.'))}
    else{return(cat('Treatment effects: Yes.\nInter-study heterogeneity of treatment effects: No.'))}
  }
  else{return(cat('Please choose one of critera: {AIC,BIC,DIC,BST}.'))}
}