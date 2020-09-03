#Simulate dataset to serve as an example
#Variable correlations
cors<-read.csv("pairwise correlations.csv")
varnames<-colnames(cors)
#variables are all normalized except age and sex. Mean and sd approximate from original data
means<-c(20,0.41,0,0,0,0,0,0,0,0, rep(3,19))
sds<-c(15,0.5,1,1,1,1,1,1,1,1, rep(1.5,19))

#Simulate variables
covMat <- sds %*% t(sds) * cors
set.seed(13256)
library(MASS)
dat1 <- mvrnorm(n = 75, mu = means, Sigma = covMat, empirical = TRUE)
dat1<-data.frame(dat1)
names(dat1)<-varnames
#Round Sex variable
dat1$Sex<-round(dat1$Sex)

#add household and village and add clustering
dat1$PID<-row.names(dat1)
dat1$Village<-sample(1:3,nrow(dat1),replace=TRUE)
#roughly 9 families per village
dat1$Family<-sample(1:9,nrow(dat1),replace=TRUE)+dat1$Village*10

#sds for random effects for each variable. Note, no attempt to recreate correlation patterns in household or village clustering, and same sd applied to all disgust variables.
famsd<-c(0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.3,0.3)
vilsd<-c(0,0,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2)
families<-unique(dat1$Family)
fameffects<-data.frame(sapply(famsd,function(x) rnorm(length(families),0,x)))
names(fameffects)<-varnames[1:10]
fameffects$Family<-unique(dat1$Family)

fameffects<-merge(dat1[,c("PID","Family")],fameffects)
dat1<-dat1[order(dat1$Family),]
dat1[varnames[1:10]]<-dat1[varnames[1:10]] + fameffects[varnames[1:10]]

vileffects<-data.frame(sapply(vilsd,function(x) rnorm(1:3,0,x)))
names(vileffects)<-varnames[1:10]
vileffects$Village<-1:3
vileffects<-merge(dat1[,c("PID","Village")],vileffects)
dat1<-dat1[order(dat1$Village),]
dat1[varnames[1:10]]<-dat1[varnames[1:10]] + vileffects[varnames[1:10]]

famdisgust<-data.frame(Family=unique(dat1$Family),famef=rnorm(length(families),0,0.1))
famdisgust<-merge(dat1[,c("PID","Family")],famdisgust)
dat1[varnames[11:19]]<-apply(dat1[varnames[11:19]],2,function(x) x + famdisgust$famef)

vildisgust<-data.frame(Family=unique(dat1$Village),ef=rnorm(3,0,0.2))
vildisgust<-merge(dat1[,c("PID","Village")],vildisgust)
dat1[varnames[11:19]]<-apply(dat1[varnames[11:19]],2,function(x) x + vildisgust$ef)

#make some data missing
dat1$NlnIgE[sample(nrow(dat1),15)]<-NA
dat1$NlnIL6[sample(nrow(dat1),19)]<-NA
dat1$NlnCRP[sample(nrow(dat1),11)]<-NA

#Exponentiate biomarkers back onto original scale
dat1$Ascaris_EPG<-exp(dat1$NlnAscaris_EPG)-1
dat1$Trich_EPG<-exp(dat1$NlnTrich_EPG)-1
dat1$CRP<-exp(dat1$NlnCRP)
dat1$IL6<-exp(dat1$NlnIL6)
dat1$IgE<-exp(dat1$NlnIgE)

#rename for use in code
d<-dat1

