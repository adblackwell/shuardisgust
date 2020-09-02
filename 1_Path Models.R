library(mice)
library(brms)
library(HDInterval)

#Data is available by request from the Shuar Health and Life History Project
#See https://www.shuarproject.org/data-sharing
#d<-read.delim("Disgust Paper Data.csv",stringsAsFactors = TRUE, header = TRUE)

#For the purposes of working through the code, data can be simulated with the code in "0_Simulate data as example.R". Note that the simulated data is randomized and won't reproduce the same associations or factor loadings.


#set undetectable IL6 at low value
d$IL6[d$IL6==0]<-min(d$IL6[d$IL6>0],na.rm=TRUE)

#log and then noramlize infection variables
d$NlnAscaris_EPG<-(log(d$Ascaris_EPG+1)-mean(log(d$Ascaris_EPG+1)))/sd(log(d$Ascaris_EPG+1),na.rm=TRUE)
d$NlnTrich_EPG<-(log(d$Trich_EPG+1)-mean(log(d$Trich_EPG+1)))/sd(log(d$Trich_EPG+1),na.rm=TRUE)
d$NlnCRP<-(log(d$CRP)-mean(log(d$CRP),na.rm=TRUE))/sd(log(d$CRP),na.rm=TRUE)
d$NlnIgE <-(log(d$IgE) -mean(log(d$IgE) ,na.rm=TRUE))/sd(log(d$IgE) ,na.rm=TRUE)
d$NlnIL6 <-(log(d$IL6) -mean(log(d$IL6) ,na.rm=TRUE))/sd(log(d$IL6) ,na.rm=TRUE)

#Normalize Market Variables
d$MSOL<-(d$MSOL-mean(d$MSOL,na.rm=TRUE))/sd(d$MSOL,na.rm=TRUE)
d$HOUSE<-(d$HOUSE-mean(d$HOUSE,na.rm=TRUE))/sd(d$HOUSE,na.rm=TRUE)
d$TSOL<-(d$TSOL-mean(d$TSOL,na.rm=TRUE))/sd(d$TSOL,na.rm=TRUE)


#produce disgust factors
fitd<-d[,grep("D[0-9]+_",names(d))] #Find the columns that are disgust items
library(nFactors)
ev <- eigen(cor(fitd)) # get eigenvalues
ap <- parallel(subject=nrow(fitd),var=ncol(fitd),rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)

library(psych)
fit0<-principal(fitd,nfactors=3,rotate="oblimin")
d$PDSCont<-fit0$scores[,1]
d$PDSFood<-fit0$scores[,2]
d$PDSPest<-fit0$scores[,3]
##In most analyses, excluding PDSPest, since it is not very interesting.

#Infection factors
#Subset to use for imputation
dsub<-d[,c("Age","Sex","Village","PDSCont","PDSFood","MSOL","HOUSE","TSOL","NlnAscaris_EPG","NlnTrich_EPG","NlnIgE","NlnCRP","NlnIL6")]
##Impute missing
set.seed(8274)
dimp<-mice(dsub,m=10,method="rf") 
dsets<-complete(dimp,"all")
#Add additional vaiables back onto dataset
for(i in 1:length(dsets)){
  dsets[[i]]<-cbind(dsets[[i]],d[,c("PID","Family","CRP","IgE","IL6","Ascaris_EPG","Trich_EPG","PDSPest")])
}
#Long format (all 10 datasets together) for principal components
dlong<-complete(dimp,"long")

#Functions to calculate values for other household members and community members outside the household
otherFam<-function(x){
  HH<-ave(x,dd$Family,FUN=function(x) mean(x,na.rm=TRUE))
  HH<-(HH*dd$inHH - x)/(dd$inHH-1)
  HH[dd$inHH==1]<-0
  HH
}
otherVillage<-function(x){
  HH<-ave(x,dd$Family,FUN=function(x) mean(x,na.rm=TRUE))
  Vi<-ave(x,dd$Village,FUN=function(x) mean(x,na.rm=TRUE))
  Vi<-(Vi*dd$inVi - HH*dd$inHH)/(dd$inVi-dd$inHH)
  Vi
}

#do the principal on all imputed combined so don't get different loadings across imputations.
#scree plot
dlongfit<-dlong[,c("NlnAscaris_EPG","NlnTrich_EPG","NlnIgE","NlnCRP","NlnIL6")]
ev <- eigen(cor(dlongfit)) # get eigenvalues
ap <- parallel(subject=nrow(dlongfit),var=ncol(dlongfit),rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)

fit2<-principal(dlongfit,nfactors=2,rotate="oblimin")
dlong$Parasites <- fit2$scores[,1]
dlong$Inflam <- fit2$scores[,2]

#Add component loadings and household and village values to the 10 imputed datasets.
##include orignal unimputed for some purposes
dsets[[length(dsets)+1]]<-d

for(i in 1:length(dsets)){
  dd<-dsets[[i]]
  pp<-predict(fit2,dsets[[i]][,c("NlnAscaris_EPG","NlnTrich_EPG","NlnIgE","NlnCRP","NlnIL6")],dlong[,c("NlnAscaris_EPG","NlnTrich_EPG","NlnIgE","NlnCRP","NlnIL6")])
  dd$Parasites<-pp[,1]
  dd$Inflam<-pp[,2]
  #Mean values for other family members and other village members
  dd$inHH<-ave(rep(1,nrow(dd)),dd$Family,FUN=sum)
  dd$inVi<-ave(rep(1,nrow(dd)),dd$Village,FUN=sum)
  
  dd$HHPDSCont<-otherFam(dd$PDSCont)
  dd$ViPDSCont<-otherVillage(dd$PDSCont)
  dd$HHPDSFood<-otherFam(dd$PDSFood)
  dd$ViPDSFood<-otherVillage(dd$PDSFood)
  dd$HHPDSPest<-otherFam(dd$PDSPest)
  dd$ViPDSPest<-otherVillage(dd$PDSPest)
  dd$HHMSOL<-otherFam(dd$MSOL)
  dd$ViMSOL<-otherVillage(dd$MSOL)
  dd$HHHOUSE<-otherFam(dd$HOUSE)
  dd$ViHOUSE<-otherVillage(dd$HOUSE)
  dd$HHTSOL<-otherFam(dd$TSOL)
  dd$ViTSOL<-otherVillage(dd$TSOL)
  #dd$HHMarket<-otherFam(dd$Market)
  #dd$ViMarket<-otherVillage(dd$Market)
  #dd$HHTraditional<-otherFam(dd$Traditional)
  #dd$ViTraditional<-otherVillage(dd$Traditional)
  dd$HHParasites<-otherFam(dd$Parasites)
  dd$ViParasites<-otherVillage(dd$Parasites)
  dd$HHInflam<-otherFam(dd$Inflam)
  dd$ViInflam<-otherVillage(dd$Inflam)
  dsets[[i]]<-dd
}
class(dsets)<-"list"

#Seperate unimputed datset into it's own variable
d2<-dsets[[length(dsets)]]
dsets<-dsets[1:(length(dsets)-1)]

#Path analysis models
brmhelp<-function(f1,...) brm_multiple(f1, data=dsets, family=gaussian(), chains=2,cores=4,iter=4000,control = list(adapt_delta = 0.99, max_treedepth=12),...)

#Inflam - PDSCont
#Step 1 Direct
mod1<-brmhelp(bf(Inflam ~ Age + Sex + PDSCont))
mod2<-brmhelp(bf(Inflam ~ Age + Sex + HHInflam))
mod3<-brmhelp(bf(Inflam ~ Age + Sex + HHPDSCont))
mod4<-brmhelp(bf(HHInflam ~ Age + Sex + HHPDSCont + (1|Family)))
mod1b<-brmhelp(bf(PDSCont ~ Age + Sex + Inflam))
mod2b<-brmhelp(bf(PDSCont ~ Age + Sex + HHInflam))
mod3b<-brmhelp(bf(PDSCont ~ Age + Sex + HHPDSCont))
mod4b<-brmhelp(bf(HHPDSCont ~ Age + Sex + HHInflam + (1|Family)))

#Step 2 Mediation on one side
mod6<-brmhelp(bf(Inflam ~ Age + Sex + PDSCont + HHPDSCont) + bf(PDSCont ~ Age + Sex + HHPDSCont) + set_rescor(FALSE))
mod7<-brmhelp(bf(Inflam ~ Age + Sex + HHInflam + HHPDSCont) + bf(HHInflam ~ Age + Sex + HHPDSCont +(1|Family)) + set_rescor(FALSE))
mod6b<-brmhelp(bf(PDSCont ~ Age + Sex + Inflam + HHInflam) + bf(Inflam ~ Age + Sex + HHInflam) + set_rescor(FALSE))
mod7b<-brmhelp(bf(PDSCont ~ Age + Sex + HHInflam + HHPDSCont) + bf(HHPDSCont ~ Age + Sex + HHInflam +(1|Family)) + set_rescor(FALSE))

#Step 3 Full Mediation one way
mod8<-brmhelp(bf(Inflam ~ Age + Sex + PDSCont + HHInflam + HHPDSCont) + bf(PDSCont ~ Age + Sex + HHPDSCont) + bf(HHInflam ~ Age + Sex + HHPDSCont +(1|Family)) + set_rescor(FALSE))
mod8b<-brmhelp(bf(PDSCont ~ Age + Sex + Inflam + HHInflam + HHPDSCont) + bf(Inflam ~ Age + Sex + HHInflam) + bf(HHPDSCont ~ Age + Sex + HHInflam +(1|Family)) + set_rescor(FALSE))

#Save models so I don't need to rerun them every time
save(mod1,mod2,mod3,mod4,mod6,mod7,mod8,file="modCI.Rdata")
save(mod1b,mod2b,mod3b,mod4b,mod6b,mod7b,mod8b,file="modCIb.Rdata")

#Inflam - PDSFood
#Step 1 Direct
modF1<-brmhelp(bf(Inflam ~ Age + Sex + PDSFood))
modF2<-brmhelp(bf(Inflam ~ Age + Sex + HHInflam))
modF3<-brmhelp(bf(Inflam ~ Age + Sex + HHPDSFood))
modF4<-brmhelp(bf(HHInflam ~ Age + Sex + HHPDSFood + (1|Family)))
modF1b<-brmhelp(bf(PDSFood ~ Age + Sex + Inflam))
modF2b<-brmhelp(bf(PDSFood ~ Age + Sex + HHInflam))
modF3b<-brmhelp(bf(PDSFood ~ Age + Sex + HHPDSFood))
modF4b<-brmhelp(bf(HHPDSFood ~ Age + Sex + HHInflam + (1|Family)))

#Step 2 Mediation on one side
modF6<-brmhelp(bf(Inflam ~ Age + Sex + PDSFood + HHPDSFood) + bf(PDSFood ~ Age + Sex + HHPDSFood) + set_rescor(FALSE))
modF7<-brmhelp(bf(Inflam ~ Age + Sex + HHInflam + HHPDSFood) + bf(HHInflam ~ Age + Sex + HHPDSFood +(1|Family)) + set_rescor(FALSE))
modF6b<-brmhelp(bf(PDSFood ~ Age + Sex + Inflam + HHInflam) + bf(Inflam ~ Age + Sex + HHInflam) + set_rescor(FALSE))
modF7b<-brmhelp(bf(PDSFood ~ Age + Sex + HHInflam + HHPDSFood) + bf(HHPDSFood ~ Age + Sex + HHInflam +(1|Family)) + set_rescor(FALSE))

#Step 3 Full Mediation one way
modF8<-brmhelp(bf(Inflam ~ Age + Sex + PDSFood + HHInflam + HHPDSFood) + bf(PDSFood ~ Age + Sex + HHPDSFood) + bf(HHInflam ~ Age + Sex + HHPDSFood +(1|Family)) + set_rescor(FALSE))
modF8b<-brmhelp(bf(PDSFood ~ Age + Sex + Inflam + HHInflam + HHPDSFood) + bf(Inflam ~ Age + Sex + HHInflam) + bf(HHPDSFood ~ Age + Sex + HHInflam +(1|Family)) + set_rescor(FALSE))

save(modF1,modF2,modF3,modF4,modF6,modF7,modF8,file="modF.Rdata")
save(modF1b,modF2b,modF3b,modF4b,modF6b,modF7b,modF8b,file="modFb.Rdata")

#Parasites - PDSCont
#Step 1 Direct
modP1<-brmhelp(bf(Parasites ~ Age + Sex + PDSCont))
modP2<-brmhelp(bf(Parasites ~ Age + Sex + HHParasites))
modP3<-brmhelp(bf(Parasites ~ Age + Sex + HHPDSCont))
modP4<-brmhelp(bf(HHParasites ~ Age + Sex + HHPDSCont +(1|Family)))
modP1b<-brmhelp(bf(PDSCont ~ Age + Sex + Parasites))
modP2b<-brmhelp(bf(PDSCont ~ Age + Sex + HHParasites))
modP3b<-brmhelp(bf(PDSCont ~ Age + Sex + HHPDSCont))
modP4b<-brmhelp(bf(HHPDSCont ~ Age + Sex + HHParasites +(1|Family)))

#Step 2 Mediation on one side
modP6<-brmhelp(bf(Parasites ~ Age + Sex + PDSCont + HHPDSCont) + bf(PDSCont ~ Age + Sex + HHPDSCont) + set_rescor(FALSE))
modP7<-brmhelp(bf(Parasites ~ Age + Sex + HHParasites + HHPDSCont) + bf(HHParasites ~ Age + Sex + HHPDSCont +(1|Family)) + set_rescor(FALSE))
modP6b<-brmhelp(bf(PDSCont ~ Age + Sex + Parasites + HHParasites) + bf(Parasites ~ Age + Sex + HHParasites) + set_rescor(FALSE))
modP7b<-brmhelp(bf(PDSCont ~ Age + Sex + HHParasites + HHPDSCont) + bf(HHPDSCont ~ Age + Sex + HHParasites +(1|Family)) + set_rescor(FALSE))

#Step 3 Full Mediation one way
modP8<-brmhelp(bf(Parasites ~ Age + Sex + PDSCont + HHParasites + HHPDSCont) + bf(PDSCont ~ Age + Sex + HHPDSCont) + bf(HHParasites ~ Age + Sex + HHPDSCont +(1|Family)) + set_rescor(FALSE))
modP8b<-brmhelp(bf(PDSCont ~ Age + Sex + Parasites + HHParasites + HHPDSCont) + bf(Parasites ~ Age + Sex + HHParasites) + bf(HHPDSCont ~ Age + Sex + HHParasites +(1|Family)) + set_rescor(FALSE))

save(modP1,modP2,modP3,modP4,modP6,modP7,modP8,file="modP.Rdata")
save(modP1b,modP2b,modP3b,modP4b,modP6b,modP7b,modP8b,file="modPb.Rdata")

#Parasites - PDSFood
#Step 1 Direct
modPF1<-brmhelp(bf(Parasites ~ Age + Sex + PDSFood))
modPF2<-brmhelp(bf(Parasites ~ Age + Sex + HHParasites))
modPF3<-brmhelp(bf(Parasites ~ Age + Sex + HHPDSFood))
modPF4<-brmhelp(bf(HHParasites ~ Age + Sex + HHPDSFood +(1|Family)))
modPF1b<-brmhelp(bf(PDSFood ~ Age + Sex + Parasites))
modPF2b<-brmhelp(bf(PDSFood ~ Age + Sex + HHParasites))
modPF3b<-brmhelp(bf(PDSFood ~ Age + Sex + HHPDSFood))
modPF4b<-brmhelp(bf(HHPDSFood ~ Age + Sex + HHParasites +(1|Family)))

#Step 2 Mediation on one side
modPF6<-brmhelp(bf(Parasites ~ Age + Sex + PDSFood + HHPDSFood) + bf(PDSFood ~ Age + Sex + HHPDSFood) + set_rescor(FALSE))
modPF7<-brmhelp(bf(Parasites ~ Age + Sex + HHParasites + HHPDSFood) + bf(HHParasites ~ Age + Sex + HHPDSFood +(1|Family)) + set_rescor(FALSE))
modPF6b<-brmhelp(bf(PDSFood ~ Age + Sex + Parasites + HHParasites) + bf(Parasites ~ Age + Sex + HHParasites) + set_rescor(FALSE))
modPF7b<-brmhelp(bf(PDSFood ~ Age + Sex + HHParasites + HHPDSFood) + bf(HHPDSFood ~ Age + Sex + HHParasites +(1|Family)) + set_rescor(FALSE))

#Step 3 Full Mediation one way
modPF8<-brmhelp(bf(Parasites ~ Age + Sex + PDSFood + HHParasites + HHPDSFood) + bf(PDSFood ~ Age + Sex + HHPDSFood) + bf(HHParasites ~ Age + Sex + HHPDSFood +(1|Family)) + set_rescor(FALSE))
modPF8b<-brmhelp(bf(PDSFood ~ Age + Sex + Parasites + HHParasites + HHPDSFood) + bf(Parasites ~ Age + Sex + HHParasites) + bf(HHPDSFood ~ Age + Sex + HHParasites +(1|Family)) + set_rescor(FALSE))

save(modPF1,modPF2,modPF3,modPF4,modPF6,modPF7,modPF8,file="modPF.Rdata")
save(modPF1b,modPF2b,modPF3b,modPF4b,modPF6b,modPF7b,modPF8b,file="modPFb.Rdata")

#Step 4 Single model
#This model has a regularizing prior for village level household effect, since in a few models when controlling for disgust these are estimated as negative effects (with) with wide intervals. Cross-validation suggests model with prior is a better fit.
pr<-c(prior(normal(0,0.15), class = b, coef = ViInflam, resp = HHInflam),
      prior(normal(0,0.15), class = b, coef = ViInflam, resp = Inflam))

mod11<-brmhelp(
  bf(Inflam ~ Age + Sex + PDSCont+HHPDSCont+HHInflam+ViInflam) +
    bf(PDSCont ~ Age + Sex + HHPDSCont+ViPDSCont) +
    bf(HHInflam ~ HHPDSCont + ViInflam + (1|p|Family)) + 
    bf(HHPDSCont ~ ViPDSCont + (1|p|Family)) +
    set_rescor(FALSE),prior=pr)

modF11<-brmhelp(
  bf(Inflam ~ Age + Sex + PDSFood+HHPDSFood+HHInflam+ViInflam) +
    bf(PDSFood ~ Age + Sex + HHPDSFood+ViPDSFood) +
    bf(HHInflam ~ HHPDSFood + ViInflam + (1|p|Family)) + 
    bf(HHPDSFood ~ ViPDSFood + (1|p|Family)) +
    set_rescor(FALSE),prior=pr)

modP11<-brmhelp(
  bf(Parasites ~ Age + Sex + PDSCont+HHPDSCont+HHParasites+ViParasites) +
    bf(PDSCont ~ Age + Sex + HHPDSCont+ViPDSCont) +
    bf(HHParasites ~ HHPDSCont + ViParasites + (1|p|Family)) + 
    bf(HHPDSCont ~ ViPDSCont + (1|p|Family)) +
    set_rescor(FALSE))

modPF11<-brmhelp(
  bf(Parasites ~ Age + Sex + PDSFood+HHPDSFood+HHParasites+ViParasites) +
    bf(PDSFood ~ Age + Sex + HHPDSFood+ViPDSFood) +
    bf(HHParasites ~ HHPDSFood +ViParasites +(1|p|Family)) + 
    bf(HHPDSFood ~ ViPDSFood+(1|p|Family)) +
    set_rescor(FALSE))

save(mod11,modP11,modF11,modPF11,file="mod11.Rdata")

#b is reverse causal direction
mod11b<-brmhelp(
  bf(Inflam ~ Age + Sex + HHInflam+ViInflam) +
    bf(PDSCont ~ Age + Sex + HHPDSCont+Inflam+HHInflam+ViPDSCont) +
    bf(HHInflam ~ ViInflam + (1|p|Family)) + 
    bf(HHPDSCont ~ HHInflam + ViPDSCont + (1|p|Family)) +
    set_rescor(FALSE))

modF11b<-brmhelp(
  bf(Inflam ~ Age + Sex + HHInflam+ViInflam) +
    bf(PDSFood ~ Age + Sex + HHPDSFood+Inflam+HHInflam+ViPDSFood) +
    bf(HHInflam ~ ViInflam + (1|p|Family)) + 
    bf(HHPDSFood ~ HHInflam + ViPDSFood + (1|p|Family)) +
    set_rescor(FALSE))

modP11b<-brmhelp(
  bf(Parasites ~ Age + Sex + HHParasites+ViParasites) +
    bf(PDSCont ~ Age + Sex + HHPDSCont+Parasites+HHParasites+ViPDSCont) +
    bf(HHParasites ~ ViParasites + (1|p|Family)) +
    bf(HHPDSCont ~ HHParasites + ViPDSCont + (1|p|Family)) +
    set_rescor(FALSE))

modPF11b<-brmhelp(
  bf(Parasites ~ Age + Sex + HHParasites+ViParasites) +
    bf(PDSFood ~ Age + Sex + HHPDSFood+Parasites+HHParasites+ViPDSFood) +
    bf(HHParasites ~ ViParasites +(1|p|Family)) +
    bf(HHPDSFood ~ HHParasites +ViPDSFood+(1|p|Family)) +
    set_rescor(FALSE))

save(mod11b,modP11b,modF11b,modPF11b,file="mod11b.Rdata")

#MI models
mod11MI<-brmhelp(
  bf(Inflam ~ Age + Sex + PDSCont+HHPDSCont+HHInflam+(1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(PDSCont ~ Age + Sex + HHPDSCont+(1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(HHInflam ~ HHPDSCont + (1|q|Village) + (1|p|Family) + MSOL + TSOL + HOUSE) + 
    bf(HHPDSCont ~ (1|q|Village) + (1|p|Family) + MSOL + TSOL + HOUSE) +
    set_rescor(FALSE))

modF11MI<-brmhelp(
  bf(Inflam ~ Age + Sex + PDSFood+HHPDSFood+HHInflam+(1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(PDSFood ~ Age + Sex + HHPDSFood+(1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(HHInflam ~ HHPDSFood + (1|q|Village) + (1|p|Family) + MSOL + TSOL + HOUSE) + 
    bf(HHPDSFood ~ (1|q|Village) + (1|p|Family) + MSOL + TSOL + HOUSE) +
    set_rescor(FALSE))

modP11MI<-brmhelp(
  bf(Parasites ~ Age + Sex + PDSCont+HHPDSCont+HHParasites++ (1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(PDSCont ~ Age + Sex + HHPDSCont+(1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(HHParasites ~ HHPDSCont + (1|q|Village) + (1|p|Family) + MSOL + TSOL + HOUSE) + 
    bf(HHPDSCont ~ (1|q|Village) + (1|p|Family) + MSOL + TSOL + HOUSE) +
    set_rescor(FALSE))

modPF11MI<-brmhelp(
  bf(Parasites ~ Age + Sex + PDSFood+HHPDSFood+HHParasites + (1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(PDSFood ~ Age + Sex + HHPDSFood + (1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(HHParasites ~ HHPDSFood + (1|q|Village) +(1|p|Family) + MSOL + TSOL + HOUSE) + 
    bf(HHPDSFood ~ (1|q|Village) + (1|p|Family) + MSOL + TSOL + HOUSE) +
    set_rescor(FALSE))

save(mod11MI,modP11MI,modF11MI,modPF11MI,file="mod11MI.Rdata")
#with family as RE

mod11MI2<-brmhelp(
  bf(Inflam ~ Age + Sex + PDSCont + (1|p|Family) +(1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(PDSCont ~ Age + Sex + Inflam + (1|p|Family) +(1|q|Village) + MSOL + TSOL + HOUSE) +
    set_rescor(FALSE))

modF11MI2<-brmhelp(
  bf(Inflam ~ Age + Sex + PDSFood + (1|p|Family)+(1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(PDSFood ~ Age + Sex + Inflam + (1|p|Family) +(1|q|Village) + MSOL + TSOL + HOUSE) +
    set_rescor(FALSE))

modP11MI2<-brmhelp(
  bf(Parasites ~ Age + Sex + PDSCont + (1|p|Family) + (1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(PDSCont ~ Age + Sex + Parasites + (1|p|Family) +(1|q|Village) + MSOL + TSOL + HOUSE) +
    set_rescor(FALSE))

modPF11MI2<-brmhelp(
  bf(Parasites ~ Age + Sex + PDSFood + (1|p|Family) + (1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(PDSFood ~ Age + Sex + Parasites + (1|p|Family) + (1|q|Village) + MSOL + TSOL + HOUSE) +
    set_rescor(FALSE))
save(mod11MI2,modP11MI2,modF11MI2,modPF11MI2,file="mod11MI2.Rdata")

# 
# 
# betaCI<-function(beta){
#   paste0(format(round(beta[1],2),nsmall=2,digits=2)," (",format(round(beta[2],2),nsmall=2,digits=2), " - ", format(round(beta[3],2),nsmall=2,digits=2),")")
# }
# 
# modeltable<-function(models,names){
#   names(models)<-names
#   df<-lapply(models,function(mod) data.frame(Variable=row.names(summary(mod)$fixed),Beta=apply(summary(mod)$fixed[,c(1,3,4)],1,betaCI)))
#   for(i in 1:length(df)) colnames(df[[i]])[2]<-names[i]
#   df2<-Reduce(function(m,y) merge(m,y,all=TRUE),df)
# }
# 
# betaCI<-function(beta){
#   paste0(format(round(beta[1],2),nsmall=2,digits=2)," (",format(round(beta[2],2),nsmall=2,digits=2), " - ", format(round(beta[3],2),nsmall=2,digits=2),")")
# }

#descriptive tables
output<-function(x) paste0(format(round(mean(x),2),nsmall=2,digits=2)," (",format(round(sd(x),2),nsmall=2,digits=2),")")
ms<-aggregate(cbind(Inflam,Parasites,PDSCont,PDSFood,MSOL,HOUSE,TSOL,Age,Sex)~Village,data=dlong,FUN=output)

#simple tests
modC1x<-brmhelp(bf(Inflam ~ Age + Sex + PDSCont + (1|Family) + (1|Village)))
modF1x<-brmhelp(bf(Inflam ~ Age + Sex + PDSFood + (1|Family) + (1|Village)))
modP1x<-brmhelp(bf(Parasites ~ Age + Sex + PDSCont + (1|Family) + (1|Village)))
modPF1x<-brmhelp(bf(Parasites ~ Age + Sex + PDSFood + (1|Family) + (1|Village)))
save(modC1x,modF1x,modP1x,modPF1x,file="simpleREmodels.Rdata")

modI2x<-brmhelp(bf(Inflam ~ Age + Sex + PDSCont + PDSFood + (1|Family) + (1|Village)))
modP2x<-brmhelp(bf(Parasites ~ Age + Sex + PDSCont + PDSFood + (1|Family) + (1|Village)))
modF2x<-brmhelp(bf(PDSFood ~ Age + Sex + Inflam + Parasites + (1|Family) + (1|Village)))
modC2x<-brmhelp(bf(PDSCont ~ Age + Sex + Inflam + Parasites + (1|Family) + (1|Village)))
save(modI2x,modP2x,modF2x,modC2x,file="REmodels2.Rdata")

modI3x<-brmhelp(bf(Inflam ~ Age + Sex + PDSCont + PDSFood + MSOL + HOUSE + TSOL + (1|Family) + (1|Village)))
modP3x<-brmhelp(bf(Parasites ~ Age + Sex + PDSCont + PDSFood + MSOL + HOUSE + TSOL + (1|Family) + (1|Village)))
modP3x2<-brmhelp(bf(Parasites ~ Age + Sex + PDSCont + PDSFood + HOUSE + (1|Family)))

#Pest models
mod1pestx<-brmhelp(bf(Inflam ~ Age + Sex + PDSPest + (1|Family) + (1|Village)))
modP1pestx<-brmhelp(bf(Parasites ~ Age + Sex + PDSPest + (1|Family) + (1|Village)))
mod2pestx<-brmhelp(bf(PDSPest ~ Age + Sex + Inflam + Parasites + (1|Family) + (1|Village)))
save(mod1pestx,modP1pestx,mod2pestx,file="Pest Models.Rdata")

#correlations
cors<-cor(d2[,c("PDSCont","PDSFood","Inflam","Parasites","MSOL","HOUSE","TSOL","HHPDSCont","HHPDSFood","HHInflam","HHParasites","ViPDSCont","ViPDSFood","ViInflam","ViParasites")],use="pairwise.complete.obs")
heatmap(cors,Rowv=NA,Colv=NA,)
