

##Rstan 2.21 is causing crashing when running multiple models. Use this to revert to 2.19 if needed
#require(devtools)
#install_version("rstan", version = "2.19.3", repos = "http://cran.us.r-project.org")

library(mice)
library(brms)
library(HDInterval)
library(psych)
library(nFactors)

#sessionInfo()

# R version 4.0.3 (2020-10-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
# [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] nFactors_2.4.1   lattice_0.20-41  psych_2.0.9      HDInterval_0.2.2 brms_2.14.0      Rcpp_1.0.5 
# [7] mice_3.11.0     
# 
# loaded via a namespace (and not attached):
# [1] nlme_3.1-149         matrixStats_0.57.0   xts_0.12.1           threejs_0.3.3       
# [5] rstan_2.19.3         tools_4.0.3          backports_1.1.10     R6_2.5.0            
# [9] DT_0.16              colorspace_1.4-1     mnormt_2.0.2         tidyselect_1.1.0    
# [13] gridExtra_2.3        prettyunits_1.1.1    processx_3.4.4       Brobdingnag_1.2-6   
# [17] emmeans_1.5.1        compiler_4.0.3       cli_2.1.0            shinyjs_2.0.0       
# [21] sandwich_3.0-0       colourpicker_1.1.0   scales_1.1.1         dygraphs_1.1.1.6    
# [25] mvtnorm_1.1-1        ggridges_0.5.2       callr_3.5.1          stringr_1.4.0       
# [29] digest_0.6.27        StanHeaders_2.21.0-6 base64enc_0.1-3      pkgconfig_2.0.3     
# [33] htmltools_0.5.0      fastmap_1.0.1        htmlwidgets_1.5.2    rlang_0.4.8         
# [37] rstudioapi_0.11      shiny_1.5.0          generics_0.0.2       zoo_1.8-8           
# [41] crosstalk_1.1.0.1    gtools_3.8.2         dplyr_1.0.2          inline_0.3.16       
# [45] magrittr_1.5         loo_2.3.1            bayesplot_1.7.2      Matrix_1.2-18       
# [49] munsell_0.5.0        fansi_0.4.1          abind_1.4-5          lifecycle_0.2.0     
# [53] stringi_1.5.3        multcomp_1.4-14      MASS_7.3-53          pkgbuild_1.1.0      
# [57] plyr_1.8.6           grid_4.0.3           parallel_4.0.3       promises_1.1.1      
# [61] crayon_1.3.4         miniUI_0.1.1.1       splines_4.0.3        tmvnsim_1.0-2       
# [65] ps_1.4.0             pillar_1.4.6         igraph_1.2.6         markdown_1.1        
# [69] estimability_1.3     shinystan_2.5.0      reshape2_1.4.4       codetools_0.2-16    
# [73] stats4_4.0.3         rstantools_2.1.1     glue_1.4.2           RcppParallel_5.0.2  
# [77] vctrs_0.3.4          httpuv_1.5.4         gtable_0.3.0         purrr_0.3.4         
# [81] tidyr_1.1.2          assertthat_0.2.1     ggplot2_3.3.2        mime_0.9            
# [85] xtable_1.8-4         broom_0.7.2          coda_0.19-4          later_1.1.0.1       
# [89] rsconnect_0.8.16     survival_3.2-7       tibble_3.0.4         shinythemes_1.1.2   
# [93] TH.data_1.0-10       ellipsis_0.3.1       bridgesampling_1.0-0


# stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


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

#put the disgust items in a new dataframe for manipulation
fitd<-d[,grep("D[0-9]+_",names(d))] #Find the columns that are disgust items

#First reduce to a single component
fit<-principal(fitd,nfactors=1)
#single factor scores are essentially the same as just summing the items and z-scoring
plot(fit$scores,rowSums(fitd))
d$PDSTotal<-fit$scores[,1]

#Now check for subcomponents and create loadings
ev <- eigen(cor(fitd)) # get eigenvalues
ap <- parallel(subject=nrow(fitd),var=ncol(fitd),rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)
#Parallel analysis and optimal coordinates suggest 3 factors
fit0<-principal(fitd,nfactors=3,rotate="oblimin")
d$PDSCont<-fit0$scores[,1]
d$PDSFood<-fit0$scores[,2]
d$PDSPest<-fit0$scores[,3]
##In most analyses, excluding PDSPest, since it is not very interesting.

#the three "raw meat" questions are relatively similar to each other so merge into one average, and check if loadings are similar
cor(fitd[,13:15]) #correlation between raw questions
fitd2<-fitd
fitd2$D13_15_Raw<-rowMeans(fitd2[,13:15])
fitd2<-fitd2[,-(13:15)]
ev <- eigen(cor(fitd2)) # get eigenvalues
ap <- parallel(subject=nrow(fitd2),var=ncol(fitd2),rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)
#both yield essentially the same loadings and three factors. Only Factor 2 loading shifts a bit, but still highly correlated
fit1<-principal(fitd2,nfactors=3,rotate="oblimin")
plot(fit0$scores[,1],fit1$scores[,1],main=cor(fit0$scores[,1],fit1$scores[,1]))
plot(fit0$scores[,2],fit1$scores[,2],main=cor(fit0$scores[,2],fit1$scores[,2]))
plot(fit0$scores[,3],fit1$scores[,3],main=cor(fit0$scores[,3],fit1$scores[,3]))

#Correlations between the single component and three components
cor(d[,c("PDSTotal","PDSCont","PDSFood","PDSPest")])


#Infection factors
#Subset to use for imputation
dsub<-d[,c("Age","Sex","Village","PDSTotal","PDSCont","PDSFood","PDSPest","MSOL","HOUSE","TSOL","NlnAscaris_EPG","NlnTrich_EPG","NlnIgE","NlnCRP","NlnIL6")]
##Impute missing
set.seed(8274)
dimp<-mice(dsub,m=10,method="rf") 
dsets<-complete(dimp,"all")
#Add additional variables back onto dataset
for(i in 1:length(dsets)){
  dsets[[i]]<-cbind(dsets[[i]],d[,c("PID","Family","CRP","IgE","IL6","Ascaris_EPG","Trich_EPG")])
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
#Two factors supported
fit2<-principal(dlongfit,nfactors=2,rotate="oblimin")
dlong$Parasites <- fit2$scores[,1]
dlong$Inflam <- fit2$scores[,2]

#Add component loadings and household and village values to the 10 imputed datasets.
##include original unimputed for some purposes
dsets[[length(dsets)+1]]<-d

for(i in 1:length(dsets)){
  dd<-dsets[[i]]
  pp<-predict(fit2,dsets[[i]][,c("NlnAscaris_EPG","NlnTrich_EPG","NlnIgE","NlnCRP","NlnIL6")],dlong[,c("NlnAscaris_EPG","NlnTrich_EPG","NlnIgE","NlnCRP","NlnIL6")])
  dd$Parasites<-pp[,1]
  dd$Inflam<-pp[,2]
  #Mean values for other family members and other village members
  dd$inHH<-ave(rep(1,nrow(dd)),dd$Family,FUN=sum)
  dd$inVi<-ave(rep(1,nrow(dd)),dd$Village,FUN=sum)
  dd$HHPDSTotal<-otherFam(dd$PDSTotal)
  dd$ViPDSTotal<-otherVillage(dd$PDSTotal)
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
#Due to running the models on multiple imputed datasets and then combining, the Rhats and ESS will look bad. Run models with the combine=FALSE option to individually inspect the models for imputed dataset and see that individually they are fine, which is what we should be concerned with.

#Inflam - PDSTotal
#Step 1 Direct
mod1T<-brmhelp(bf(Inflam ~ Age + Sex + PDSTotal))
mod2<-brmhelp(bf(Inflam ~ Age + Sex + HHInflam))
mod3T<-brmhelp(bf(Inflam ~ Age + Sex + HHPDSTotal))
mod4T<-brmhelp(bf(HHInflam ~ Age + Sex + HHPDSTotal + (1|Family)))

#Step 2 Mediation on one side
mod6T<-brmhelp(bf(Inflam ~ Age + Sex + PDSTotal + HHPDSTotal) + bf(PDSTotal ~ Age + Sex + HHPDSTotal) + set_rescor(FALSE))
mod7T<-brmhelp(bf(Inflam ~ Age + Sex + HHInflam + HHPDSTotal) + bf(HHInflam ~ Age + Sex + HHPDSTotal +(1|Family)) + set_rescor(FALSE))

#Step 3 Full Mediation one way
mod8T<-brmhelp(bf(Inflam ~ Age + Sex + PDSTotal + HHInflam + HHPDSTotal) + bf(PDSTotal ~ Age + Sex + HHPDSTotal) + bf(HHInflam ~ Age + Sex + HHPDSTotal +(1|Family)) + set_rescor(FALSE))

#Save models so I don't need to rerun them every time
save(mod1T,mod2,mod3T,mod4T,mod6T,mod7T,mod8T,file="modTI.Rdata")

#Parasites - PDSTotal
#Step 1 Direct
modP1T<-brmhelp(bf(Parasites ~ Age + Sex + PDSTotal))
modP2<-brmhelp(bf(Parasites ~ Age + Sex + HHParasites))
modP3T<-brmhelp(bf(Parasites ~ Age + Sex + HHPDSTotal))
modP4T<-brmhelp(bf(HHParasites ~ Age + Sex + HHPDSTotal + (1|Family)))

#Step 2 Mediation on one side
modP6T<-brmhelp(bf(Parasites ~ Age + Sex + PDSTotal + HHPDSTotal) + bf(PDSTotal ~ Age + Sex + HHPDSTotal) + set_rescor(FALSE))
modP7T<-brmhelp(bf(Parasites ~ Age + Sex + HHParasites + HHPDSTotal) + bf(HHParasites ~ Age + Sex + HHPDSTotal +(1|Family)) + set_rescor(FALSE))

#Step 3 Full Mediation one way
modP8T<-brmhelp(bf(Parasites ~ Age + Sex + PDSTotal + HHParasites + HHPDSTotal) + bf(PDSTotal ~ Age + Sex + HHPDSTotal) + bf(HHParasites ~ Age + Sex + HHPDSTotal +(1|Family)) + set_rescor(FALSE))

#Save models so I don't need to rerun them every time
save(modP1T,modP2,modP3T,modP4T,modP6T,modP7T,modP8T,file="modTP.Rdata")


#Step 4 Single model
#This model has a regularizing prior for village level household effect, since in a few models when controlling for disgust these are estimated as negative effects (with) with wide intervals. Cross-validation suggests model with prior is a better fit.
pr<-c(prior(normal(0,0.15), class = b, coef = ViInflam, resp = HHInflam),
      prior(normal(0,0.15), class = b, coef = ViInflam, resp = Inflam))

#with PDSTotal
modT11<-brmhelp(
  bf(Inflam ~ Age + Sex + PDSTotal+HHPDSTotal+HHInflam+ViInflam) +
    bf(PDSTotal ~ Age + Sex + HHPDSTotal+ViPDSTotal) +
    bf(HHInflam ~ HHPDSTotal + ViInflam + (1|p|Family)) + 
    bf(HHPDSTotal ~ ViPDSTotal + (1|p|Family)) +
    set_rescor(FALSE),prior=pr)

modTP11<-brmhelp(
  bf(Parasites ~ Age + Sex + PDSTotal+HHPDSTotal+HHParasites+ViParasites) +
    bf(PDSTotal ~ Age + Sex + HHPDSTotal+ViPDSTotal) +
    bf(HHParasites ~ HHPDSTotal + ViParasites + (1|p|Family)) + 
    bf(HHPDSTotal ~ ViPDSTotal + (1|p|Family)) +
    set_rescor(FALSE))

save(modT11,modTP11,file="modT11.Rdata")

#b is reverse causal direction
modT11b<-brmhelp(
  bf(Inflam ~ Age + Sex + HHInflam+ViInflam) +
    bf(PDSTotal ~ Age + Sex + HHPDSTotal+Inflam+HHInflam+ViPDSTotal) +
    bf(HHInflam ~ ViInflam + (1|p|Family)) + 
    bf(HHPDSTotal ~ HHInflam + ViPDSTotal + (1|p|Family)) +
    set_rescor(FALSE))

modTP11b<-brmhelp(
  bf(Parasites ~ Age + Sex + HHParasites+ViParasites) +
    bf(PDSTotal ~ Age + Sex + HHPDSTotal+Parasites+HHParasites+ViPDSTotal) +
    bf(HHParasites ~ ViParasites + (1|p|Family)) +
    bf(HHPDSTotal ~ HHParasites + ViPDSTotal + (1|p|Family)) +
    set_rescor(FALSE))

save(modT11b,modTP11b,file="modT11b.Rdata")

#MI models
#with PDSTotal
modT11MI<-brmhelp(
  bf(Inflam ~ Age + Sex + PDSTotal+HHPDSTotal+HHInflam+(1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(PDSTotal ~ Age + Sex + HHPDSTotal+(1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(HHInflam ~ HHPDSTotal + (1|q|Village) + (1|p|Family) + MSOL + TSOL + HOUSE) + 
    bf(HHPDSTotal ~ (1|q|Village) + (1|p|Family) + MSOL + TSOL + HOUSE) +
    set_rescor(FALSE))

modTP11MI<-brmhelp(
  bf(Parasites ~ Age + Sex + PDSTotal+HHPDSTotal+HHParasites++ (1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(PDSTotal ~ Age + Sex + HHPDSTotal+(1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(HHParasites ~ HHPDSTotal + (1|q|Village) + (1|p|Family) + MSOL + TSOL + HOUSE) + 
    bf(HHPDSTotal ~ (1|q|Village) + (1|p|Family) + MSOL + TSOL + HOUSE) +
    set_rescor(FALSE))

save(modT11MI,modTP11MI,file="modT11MI.Rdata")

#with family as RE

#PDSTotal
modT11MI2<-brmhelp(
  bf(Inflam ~ Age + Sex + PDSTotal + (1|p|Family) +(1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(PDSTotal ~ Age + Sex + Inflam + (1|p|Family) +(1|q|Village) + MSOL + TSOL + HOUSE) +  set_rescor(FALSE))

modTP11MI2<-brmhelp(
  bf(Parasites ~ Age + Sex + PDSTotal + (1|p|Family) + (1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(PDSTotal ~ Age + Sex + Parasites + (1|p|Family) +(1|q|Village) + MSOL + TSOL + HOUSE) +
    set_rescor(FALSE))

save(modT11MI2,modTP11MI2,file="modT11MI2.Rdata")

#simple tests
modTI1x<-brmhelp(bf(Inflam ~ Age + Sex + PDSTotal + (1|Family) + (1|Village)))
modTP1x<-brmhelp(bf(Parasites ~ Age + Sex + PDSTotal + (1|Family) + (1|Village)))
save(modTI1x,modTP1x, modC1x,modF1x,modP1x,modPF1x,file="simpleREmodels.Rdata")

#correlations
cors<-cor(d2[,c("PDSCont","PDSFood","Inflam","Parasites","MSOL","HOUSE","TSOL","HHPDSCont","HHPDSFood","HHInflam","HHParasites","ViPDSCont","ViPDSFood","ViInflam","ViParasites")],use="pairwise.complete.obs")
heatmap(cors,Rowv=NA,Colv=NA,)



#Models with three component scores
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

modV11<-brmhelp(
  bf(Inflam ~ Age + Sex + PDSPest+HHPDSPest+HHInflam+ViInflam) +
    bf(PDSPest ~ Age + Sex + HHPDSPest+ViPDSPest) +
    bf(HHInflam ~ HHPDSCont + ViInflam + (1|p|Family)) + 
    bf(HHPDSPest ~ ViPDSPest + (1|p|Family)) +
    set_rescor(FALSE),prior=pr)

modPV11<-brmhelp(
  bf(Parasites ~ Age + Sex + PDSPest+HHPDSPest+HHParasites+ViParasites) +
    bf(PDSPest ~ Age + Sex + HHPDSPest+ViPDSPest) +
    bf(HHParasites ~ HHPDSPest + ViParasites + (1|p|Family)) + 
    bf(HHPDSPest ~ ViPDSPest + (1|p|Family)) +
    set_rescor(FALSE))


save(mod11,modP11,modF11,modPF11,modV11,modPV11,file="mod11.Rdata")




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

modV11MI<-brmhelp(
  bf(Inflam ~ Age + Sex + PDSPest+HHPDSPest+HHInflam+(1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(PDSPest ~ Age + Sex + HHPDSPest+(1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(HHInflam ~ HHPDSPest + (1|q|Village) + (1|p|Family) + MSOL + TSOL + HOUSE) + 
    bf(HHPDSPest ~ (1|q|Village) + (1|p|Family) + MSOL + TSOL + HOUSE) +
    set_rescor(FALSE))

modPV11MI<-brmhelp(
  bf(Parasites ~ Age + Sex + PDSPest+HHPDSPest+HHParasites++ (1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(PDSPest ~ Age + Sex + HHPDSPest+(1|q|Village) + MSOL + TSOL + HOUSE) +
    bf(HHParasites ~ HHPDSPest + (1|q|Village) + (1|p|Family) + MSOL + TSOL + HOUSE) + 
    bf(HHPDSPest ~ (1|q|Village) + (1|p|Family) + MSOL + TSOL + HOUSE) +
    set_rescor(FALSE))


save(mod11MI,modP11MI,modF11MI,modPF11MI,modV11MI,modPV11MI, file="mod11MI.Rdata")
