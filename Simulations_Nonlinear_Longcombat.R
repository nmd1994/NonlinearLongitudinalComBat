library(invgamma)
library(lme4)
library('pbkrtest')
library('mgcv')
library(gamm4)
library('matrixStats')
nlinpvalues=matrix(0,nrow=50,ncol=20)
longcombatpvalues=matrix(0,nrow=50,ncol=20)
multlongcombat=matrix(0,nrow=50,ncol=20)
nlinmultlongcombat=matrix(0,nrow=50,ncol=20)

#simulat_data=list(0)
#NOTE: These simulations also use the NeuroCombat script as a dependency



for(q in 1:50){
  #################################
  # simulate data to run the functions
  #################################
  # 100 subjects with 5 time points each
  # 2 scanners / batches
  # 20 features
  #################################
  # set random seed
  # simulate the covariates
  simdata <- data.frame(
    subid=rep(1:200, each=5),
    age=c(rep(runif(100,10,90), each=5),rep(runif(100,10,90), each=5)),
    #age=c(rep(sample(c(10:35), 100, replace=TRUE), each=5),rep(sample(c(30:90), 100, replace=TRUE), each=5)),
    #age=c(rep(sample(c(1:40), 50, replace=TRUE), each=5),rep(sample(c(15:90), 50, replace=TRUE), each=5)),
    diagnosis=rep(c(0,1), each=500),
    time=rep(0:4, times=200)
  )
  simdata$batch=c(rep(1,500),rep(2,500))
  # simulate the brain features 
  features <- matrix(rnorm(200*5*20), nrow=1000)
  # simulate additive batch effects (normally distributed)
  gamma <- runif(n=2, min=-30, max=30)
  #tau <- runif(n=2, min=0.1, max=0.3)
  #gamma <- runif(n=2, min=-100, max=100)
  tau <- runif(n=2, min=0.1, max=0.3)
  batch.add <- matrix(c(
    rnorm(mean=gamma[1], sd=tau[1], n=20),
    rnorm(mean=gamma[2], sd=tau[2], n=20)),
    ncol=2) 
  # simulate multiplicative batch effects (inverse gamma distributed)
  lambda <- sample(c(2, 3), 2, replace=TRUE)
  theta <- sample(c(2, 5), 2, replace=TRUE)
  batch.mult <- matrix(c(
    rinvgamma(n=20, shape=lambda[1], scale=theta[1]),
    rinvgamma(n=20, shape=lambda[2], scale=theta[2])),
    ncol=2)
  # add / multiply batch effects to the features
  for(i in 1:1000){
    features[i,] <- features[i,]*batch.mult[,simdata$batch[i]] + batch.add[,simdata$batch[i]]
  }
  # add covariate effects to the features (CHOOSE FUNCTIONAL FORM HERE)
  #features <- features + 0.005*((simdata$age-50)^2)
  #features <- features + 0.005*((simdata$age-50)^2) - 0.1*simdata$time 
  features <- features - 0.0003*((simdata$age)^3) 
  #features <- features - 5*log(0.1*simdata$age-0.99) 
  #features <- features + 5*log(0.1*simdata$age-0.99) 
  #features <- features +  100*dlnorm(0.003*simdata$age,meanlog=0, sdlog=1.5)  
  #features <- features + (1/(sqrt(2*pi)*sigma))*exp(-(simdata$age-40)^(2)/(2*sigma^(2))) - 0.1*simdata$time 
  #features =  features - 1*simdata$age 
  # add subject random effect to the features (will be the same across features in this case)
  features <- features + rep(rnorm(n=200), each=5)
  # save feature names
  featurenames <- paste0('feature', 1:20)
  colnames(features) <- featurenames
  # combine into one data frame
  simdata <- data.frame(simdata, features)
  # remove some stuff no longer needed
  #rm('features', 'batch.patterns', 'batch.pattern.sample', 'i')
  #simulat_data[[q]]=simdata
  
  
  ########################################
  #Do Both Linear and Non-Linear ComBat
  
  #Then, refit both approaches with non-linear fit (separating by batch)
  
  #If intuition is correct - non-linear combat should look to be the same fit
  
  #Linear ComBat should show a huge batch difference (i.e. a shift)
  ########################################
  
  stand=c()
  subid=as.factor(simdata$subid)
  maxi=matrix(0,nrow=nrow(simdata),ncol=20)
  Fitted_site=matrix(0,nrow=nrow(simdata),ncol=20)
  Residual=matrix(0,nrow=nrow(simdata),ncol=20)
  
  for(i in 1:20){
    y=simdata[,i+5]
    Xe=cbind(y,simdata[,2:5])
    colnames(Xe)=c('y',"Age","Diagnosis",
                   "timeyrs","site")
    yur=gamm4(y~Age+as.factor(site)+1,
              random=~(1|subid),data=data.frame(Xe))
    maxi[,i]=as.matrix(yur$gam$fitted.values) 
    Fitted_site[,i]=model.matrix(yur$gam)[,3]%*%as.matrix(yur$gam$coefficients[3])
    stand[i]=as.data.frame(VarCorr(yur$mer))[2,5]
    Residual[,i]=(y-(maxi[,i]-Fitted_site[,i]))/(as.data.frame(VarCorr(yur$mer))[2,5])
    print(i)
  }
  
  data.harmonized=neuroCombat(t(Residual),as.factor(Xe[,5]))
  dat=t(data.harmonized$dat.combat)
  
  for(i in 1:20){
    dat[,i]=dat[,i]*stand[i]
    print(i)
  }
  harmonized_lin=dat+(maxi-Fitted_site)
  
  
  for(i in 1:20){
    y=simdata[,i+5]
    Xe=cbind(y,simdata[,2:5])
    colnames(Xe)=c('y',"Age","Diagnosis",
                   "timeyrs","site")
    yur=gamm4(y~s(Age,bs='cr', k= 45)+as.factor(site)+1,
              random=~(1|subid),data=data.frame(Xe))
  
    maxi[,i]=as.matrix(yur$gam$fitted.values) 
    Fitted_site[,i]=model.matrix(yur$gam)[,2]%*%as.matrix(yur$gam$coefficients[2])
    stand[i]=as.data.frame(VarCorr(yur$mer))[3,5]
    Residual[,i]=(y-(maxi[,i]-Fitted_site[,i]))/(as.data.frame(VarCorr(yur$mer))[3,5])
    print(i)
  }
  
  data.harmonized=neuroCombat(t(Residual),as.factor(Xe[,5]))
  dat=t(data.harmonized$dat.combat)
  
  for(i in 1:20){
    dat[,i]=dat[,i]*stand[i]
    print(i)
  }
  harmonized_nlin=dat+(maxi-Fitted_site)
  
  nlinlongcombatpval=c()
  mult2=c()
  resid=matrix(0,nrow = 1000, ncol = 20)
  for(i in 1:20){
    y=harmonized_nlin[,i]
    Xe=cbind(y,simdata[,2:5])
    colnames(Xe)=c('y',"Age","Diagnosis",
                   "timeyrs","site")
    fmLarge=gamm4(y~s(Age,bs='cr', k = 45)+as.factor(site)+1,
                  random=~(1|subid),data=data.frame(Xe))
    fmSmall1=gamm4(y~s(Age,bs='cr', k= 45)+1,
                  random=~(1|subid),data=data.frame(Xe))
    red=KRmodcomp(fmLarge$mer, fmSmall1$mer)$stats$p.value
    nlinlongcombatpval[i]=red
    
    #fmSmall1=gamm4(y~s(Age,bs='cr',k=45) + as.factor(site)+1,
    #               random=~(1|subid),data=data.frame(Xe))
    
    maxi[,i]=as.matrix(fmLarge$gam$fitted.values) 
    #Fitted_site[,i]=model.matrix(fmSmall1$gam)[,2]%*%as.matrix(fmSmall1$gam$coefficients[2])
    resid[,i]= y - fmLarge$gam$fitted.values
    
    mult2[i]=as.numeric(fligner.test((y - maxi[,i]) ~ Xe$site)[3])
    print(i)
  }
  
  longcombatpval=c()
  mult=c()
  resid2=matrix(0,nrow = 1000, ncol = 20)
  for(i in 1:20){
    y=harmonized_lin[,i]
    Xe=cbind(y,simdata[,2:5])
    colnames(Xe)=c('y',"Age","Diagnosis",
                   "timeyrs","site")
    fmLarge=gamm4(y~s(Age,bs='cr',k= 45)+as.factor(site)+1,
                  random=~(1|subid),data=data.frame(Xe))
    fmSmall1=gamm4(y~s(Age,bs='cr',k= 45)+1,
                  random=~(1|subid),data=data.frame(Xe))
    red=KRmodcomp(fmLarge$mer, fmSmall1$mer)$stats$p.value
    
    #fmSmall1=gamm4(y~s(Age,bs='cr',k=45) + as.factor(site)+1,
    #               random=~(1|subid),data=data.frame(Xe))
    
    maxi[,i]=as.matrix(fmLarge$gam$fitted.values) 
    #Fitted_site[,i]=model.matrix(fmLarge$gam)[,2]%*%as.matrix(fmLarge$gam$coefficients[2])
    
    mult[i]=as.numeric(fligner.test((y - maxi[,i]) ~ Xe$site)[3])
    resid2[,i]=y - fmLarge$gam$fitted.values
    longcombatpval[i]=red
    print(i)
  }
  
  
  nlinpvalues[q,]=nlinlongcombatpval
  longcombatpvalues[q,]=longcombatpval
  multlongcombat[q,]=mult
  nlinmultlongcombat[q,]=mult2
  print(q)
}

#---------------------------
#Key For Relevant Variables
#---------------------------
#nlinpvalues = pvalues for additive effects using nonlinear longcombat
#longcombatpvalues = pvalues for additive effects using linear longcombat
#multlongcombat = pvalues for multiplicative effects using linear longcombat
#nlinmultlongcombat = pvalues for multiplicative effects using nonlinear longcombat











  
  