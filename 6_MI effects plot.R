library(brms)
library(HDInterval)
source('Plotting Functions.R', echo = TRUE)

#posterior parameter plot
cols<-c("#ffbf47","#dd6031","#0262a2","#4daa57","#311e10")[c(3,1)]


#single factor
load("modT11MI.Rdata")
dhT<-getdens(modT11MI,modTP11MI)
dhTP<-getdens(modTP11MI,modTP11MI)
dhTI<-getdens(modT11MI,modT11MI)

dhC<-getdens(mod11MI,modP11MI)
dhF<-getdens(modF11MI,modPF11MI)
dhV<-getdens(modV11MI,modPV11MI)

tiff("Figure 4.tif",width=2300,height=1500,res = 300,pointsize=13,compression="lzw")
par(mfrow=c(1,2),mar=c(5,4,1,1))
names<-c("Disgust","Inflammation","Parasites")
plot(0,0,xlim=c(-1,1),ylim=c(0,10.2),type="n",main=NA,yaxt="n",ylab=NA,xlab=NA)
shifts<-c(9.2,8.1,7.0, 5.7,4.6,3.5, 2.2,1.1,0.0)

#MSOL - Total
denspoly(dhT$den[[11]],dhT$hdi[[11]],cols[2],shifts[1])
denspoly(dhT$den[[8]],dhT$hdi[[8]],cols[1],shifts[1])
#HSOL - Total
denspoly(dhT$den[[10]],dhT$hdi[[10]],cols[2],shifts[2])
denspoly(dhT$den[[7]],dhT$hdi[[7]],cols[1],shifts[2])
#TSOL - Total
denspoly(dhT$den[[12]],dhT$hdi[[12]],cols[2],shifts[3])
denspoly(dhT$den[[9]],dhT$hdi[[9]],cols[1],shifts[3])

#MSOL - Inflam
denspoly(dhTI$den[[5]],dhTI$hdi[[5]],cols[2],shifts[4])
denspoly(dhTI$den[[2]],dhTI$hdi[[2]],cols[1],shifts[4])

#HSOL - Inflam
denspoly(dhTI$den[[4]],dhTI$hdi[[4]],cols[2],shifts[5])
denspoly(dhTI$den[[1]],dhTI$hdi[[1]],cols[1],shifts[5])

#TSOL - Inflam
denspoly(dhTI$den[[6]],dhTI$hdi[[6]],cols[2],shifts[6])
denspoly(dhTI$den[[3]],dhTI$hdi[[3]],cols[1],shifts[6])

#MSOL - Parasites
denspoly(dhTP$den[[5]],dhTP$hdi[[5]],cols[2],shifts[7])
denspoly(dhTP$den[[2]],dhTP$hdi[[2]],cols[1],shifts[7])

#HSOL - Parasites
denspoly(dhTP$den[[4]],dhTP$hdi[[4]],cols[2],shifts[8])
denspoly(dhTP$den[[1]],dhTP$hdi[[1]],cols[1],shifts[8])

#TSOL - Parasites
denspoly(dhTP$den[[6]],dhTP$hdi[[6]],cols[2],shifts[9])
denspoly(dhTP$den[[3]],dhTP$hdi[[3]],cols[1],shifts[9])

abline(v=0,lty=2)
abline(h=shifts[3]-0.15,lty=3)
abline(h=shifts[6]-0.15,lty=3)
axis(2,at=shifts+0.5,labels=rep(c("M","H","T"),3),las=2,tick=FALSE,line=-2,padj=0)
axis(2,at=shifts[c(2,5,8)]+0.5,labels=names,tick=FALSE,line=0)

mtext("A",side=3,at=-1.4,cex=1.5,line=-0.5)


cols<-c("#ffbf47","#dd6031","#0262a2","#4daa57","#311e10")[c(3,1)]
#par(mfrow=c(1,1),mar=c(5,4,3,1))
names<-c("C1:Contagion","C2:Food","C3:Various")
plot(0,0,xlim=c(-1,1),ylim=c(0,10.2),type="n",main=NA,yaxt="n",ylab=NA,xlab=NA)
shifts<-c(9.2,8.1,7.0, 5.7,4.6,3.5, 2.2,1.1,0.0)
#MSOL - Cont
denspoly(dhC$den[[11]],dhC$hdi[[11]],cols[2],shifts[1])
denspoly(dhC$den[[8]],dhC$hdi[[8]],cols[1],shifts[1])

#HSOL - Cont
denspoly(dhC$den[[10]],dhC$hdi[[10]],cols[2],shifts[2])
denspoly(dhC$den[[7]],dhC$hdi[[7]],cols[1],shifts[2])

#TSOL - Cont
denspoly(dhC$den[[12]],dhC$hdi[[12]],cols[2],shifts[3])
denspoly(dhC$den[[9]],dhC$hdi[[9]],cols[1],shifts[3])

#MSOL - Food
denspoly(dhF$den[[11]],dhF$hdi[[11]],cols[2],shifts[4])
denspoly(dhF$den[[8]],dhF$hdi[[8]],cols[1],shifts[4])

#HSOL - Food
denspoly(dhF$den[[10]],dhF$hdi[[10]],cols[2],shifts[5])
denspoly(dhF$den[[7]],dhF$hdi[[7]],cols[1],shifts[5])

#TSOL - Food
denspoly(dhF$den[[12]],dhF$hdi[[12]],cols[2],shifts[6])
denspoly(dhF$den[[9]],dhF$hdi[[9]],cols[1],shifts[6])

#MSOL - C3
denspoly(dhV$den[[11]],dhV$hdi[[11]],cols[2],shifts[7])
denspoly(dhV$den[[8]],dhV$hdi[[8]],cols[1],shifts[7])

#HSOL - C3
denspoly(dhV$den[[10]],dhV$hdi[[10]],cols[2],shifts[8])
denspoly(dhV$den[[7]],dhV$hdi[[7]],cols[1],shifts[8])

#TSOL - C3
denspoly(dhV$den[[12]],dhV$hdi[[12]],cols[2],shifts[9])
denspoly(dhV$den[[9]],dhV$hdi[[9]],cols[1],shifts[9])

abline(v=0,lty=2)
abline(h=shifts[3]-0.15,lty=3)
abline(h=shifts[6]-0.15,lty=3)
#axis(2,at=shifts+0.5,labels=rep(names,3),crt=2,tick=FALSE)
#axis(2,at=shifts[c(2,5,8)]+0.5,labels=c("MSOL","HSOL","TSOL"),crt=2,tick=FALSE,line=1.5)
axis(2,at=shifts+0.5,labels=rep(c("M","H","T"),3),las=2,tick=FALSE,line=-2,padj=0)
axis(2,at=shifts[c(2,5,8)]+0.5,labels=names,tick=FALSE,line=0)
mtext("B",side=3,at=-1.4,cex=1.5,line=-0.5)
dev.off()


