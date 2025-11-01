#--------------------------
#REAL DATA ANALYSIS 
#--------------------------
library(invgamma)
library(lme4)
library('pbkrtest')
library('mgcv')
library(gamm4)
library('matrixStats')
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
#--------------------------
#STEP 1: LOAD ADNI DATASET (ON_GITHUB PAGE) AND SUBSET TO FEATURES CHOSEN FOR CONSIDERATION
#--------------------------
COMB=ADNIDataset

Feat= c(18,19,26,27,30,31,39,40,51,52,62,63,66,67,116,117,132,133,148,149)
COMB_USE=COMB[,c(1:15,Feat)]
#-------------------------
#STEP 2: CONDUCT TESTS FOR ADDITIVE AND MULTIPLICATIVE EFFECTS ON RAW DATA
#-------------------------
COMB=COMB_USE[-c(which(COMB_USE$Diagnosis_nearest_2.0=="")),]
subid=as.factor(COMB$participant_id)

longcombatpval=c()
for(i in 16:35){
  y=COMB[,i]
  Xe=cbind(y,COMB[,2:15])
  colnames(Xe)=c('y',colnames(COMB)[2:15])
  fmLarge=gamm4(y~s(Age,bs='cr',k=100)+DLICV_baseline+as.factor(Sex)+as.factor(Study)
                +Diagnosis_nearest_2.0+1,
                random=~(1|subid),data=data.frame(Xe))$mer
  fmSmall=gamm4(y~s(Age,bs='cr',k=100)+DLICV_baseline+as.factor(Sex)
                +Diagnosis_nearest_2.0+1,
                random=~(1|subid),data=data.frame(Xe))$mer
  red=KRmodcomp(fmLarge, fmSmall)$stats$p.value
  longcombatpval[i]=red
  print(i)
}
longcombatpval=longcombatpval[16:35]
sum(longcombatpval<0.05/20)
#10 Bonferonni corrected 

#Test Multiplicative Effects 
maxi=matrix(0,nrow=1445,ncol=20)
for(i in 16:35){
  y=COMB[,i]
  Xe=cbind(y,COMB[,2:15])
  colnames(Xe)=c('y',colnames(COMB)[2:15])
  fmSmall=gamm4(y~s(Age,bs='cr',k=100)+DLICV_baseline+as.factor(Sex)
                +as.factor(Study)+Diagnosis_nearest_2.0+1,
                random=~(1|subid),data=data.frame(Xe))$mer
  maxi[,i-15]=fitted(fmSmall)
}

Residual= COMB[,16:35]-maxi
mult=c()
for(i in 1:20){
  mult[i]=fligner.test(Residual[,i] ~ COMB$Study)[3]
}
length(which(unlist(mult)<(0.05/20)))
#16 out of 20 Bonferonni Corrected Remaining Batch effects 

#-------------------------
#STEP 3: GET HARMONIZED DATA
#-------------------------
stand=c()
#subid=as.factor(simdata$subid)
subid=as.factor(COMB$participant_id)
maxi=matrix(0,nrow=nrow(COMB),ncol=20)
Fitted_site=matrix(0,nrow=nrow(COMB),ncol=20)
Residual=matrix(0,nrow=nrow(COMB),ncol=20)

for(i in 16:35){
  y=COMB[,i]
  Xe=cbind(y,COMB[,2:15])
  colnames(Xe)=c('y',colnames(COMB)[2:15])
  yur=gamm4(y~s(Age,bs='cr', k = 100)+as.factor(Study)+DLICV_baseline+as.factor(Sex)
            +Diagnosis_nearest_2.0+1,
            random=~(1|subid),data=data.frame(Xe))
  maxi[,i-15]=as.matrix(yur$gam$fitted.values) 
  Fitted_site[,i-15]=model.matrix(yur$gam)[,2]%*%as.matrix(yur$gam$coefficients[2])
  stand[i-15]=sigma(yur$mer)
  Residual[,i-15]=(y-(maxi[,i-15]-Fitted_site[,i-15]))/sigma(yur$mer)
  print(i)
}
Xd$Sex

data.harmonized=neuroCombat(t(Residual),as.factor(Xe[,3]))
dat=t(data.harmonized$dat.combat)

for(i in 1:20){
  dat[,i]=dat[,i]*stand[i]
  print(i)
}

harmonized=dat+(maxi-Fitted_site)

#---------------------------------
#STEP 4: TEST FOR ADDITIVE AND MULTIPLICATIVE BATCH EFFECTS ON HARMONIZED DATA
#---------------------------------
#Additive Effects 
nlinlongcombatpval=c()
for(i in 1:20){
  y=harmonized[,i]
  Xe=cbind(y,COMB[,2:15])
  colnames(Xe)=c('y',colnames(COMB)[2:15])
  fmLarge=gamm4(y~s(Age,bs='cr', k = 100)+DLICV_baseline+as.factor(Sex)+as.factor(Study)
                +Diagnosis_nearest_2.0+1,
                random=~(1|subid),data=data.frame(Xe))$mer
  fmSmall=gamm4(y~s(Age,bs='cr', k = 100)+DLICV_baseline+as.factor(Sex)
                +Diagnosis_nearest_2.0+1,
                random=~(1|subid),data=data.frame(Xe))$mer
  red=KRmodcomp(fmLarge, fmSmall)$stats$p.value
  nlinlongcombatpval[i]=red
  print(i)
}

length(which(unlist(nlinlongcombatpval)<0.05/20))
#No remaining additive batch effects Bonferonni Corrected 


#Multiplicative Effects 
maxi=matrix(0,nrow=1445,ncol=20)
for(i in 16:35){
  y=harmonized[,i-15]
  Xe=cbind(y,COMB[,2:15])
  colnames(Xe)=c('y',colnames(COMB)[2:15])
  fmSmall=gamm4(y~s(Age,bs='cr',k=100)+DLICV_baseline+as.factor(Sex)
                +Diagnosis_nearest_2.0+as.factor(Study)+1,
                random=~(1|subid),data=data.frame(Xe))
  maxi[,i-15]=fitted(fmSmall$gam)
  print(i)
}

Residual= harmonized - maxi
colnames(Residual)=colnames(COMB)[16:35]


mult=c()
for(i in 1:20){
  mult[i]=fligner.test(dat[,i] ~ COMB$Study)[3]
}
length(which(unlist(mult)<0.05/20))
#No remaining batch effects Bonferonni Corrected 


#-------------------------------------
#STEP 5: Look At Quick at Residuals Post Harmonization
#-------------------------------------
yuhu=list(0)
colnames(Residual)=colnames(COMB)[16:35]

par(mfrow=c(4,5))
for(i in 1:20){plot(Residual[,i],col=c(rep("blue",629),rep('red',823)),
                    xlab = "", ylab = "",xaxt='n',yaxt='n')
  abline(h = 0, col = "black", lty = 2, lwd= 3) }


