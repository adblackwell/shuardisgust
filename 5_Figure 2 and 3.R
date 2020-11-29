library(brms)
#Figure 2 and 3
source('3_Plotting Functions.R', echo = TRUE)

#Marginal values for Village Plot
ml0<-brmhelp(Inflam~PDSTotal+Age+Sex)
ml1<-brmhelp(Inflam~PDSTotal+Age+Sex+(PDSTotal|Village))
ml2<-brmhelp(Parasites~PDSTotal+Age+Sex)
ml3<-brmhelp(Parasites~PDSTotal+Age+Sex+(PDSTotal|Village))


#data conditional on age and sex
p<-posterior_samples(ml0)
d2$Inflamcond<-d2$Inflam-d2$Age*mean(p[,"b_Age"])-d2$Sex*mean(p[,"b_Sex"])
Int<-p[,"b_Intercept"]
Slope<-p[,"b_PDSTotal"]
out<-data.frame(V=0,I=mean(Int),S=mean(Slope),IL=hdi(Int)[1],IH=hdi(Int)[2],SL=hdi(Slope)[1],SH=hdi(Slope)[2])
xseq<-seq(-2.2,1.2,0.1)
mu.ci <- list(sapply( xseq , function(x) hdi( Int + Slope*x ) ))
p<-posterior_samples(ml1)
for(i in 1:3){
  Int<-p[,"b_Intercept"]+p[,paste0("r_Village[",i,",Intercept]")]
  Slope<-p[,"b_PDSTotal"]+p[,paste0("r_Village[",i,",PDSTotal]")]
  out<-rbind(out,data.frame(V=i,I=mean(Int),S=mean(Slope),IL=hdi(Int)[1],IH=hdi(Int)[2],SL=hdi(Slope)[1],SH=hdi(Slope)[2]))
  mu.ci[[i+1]] <- sapply( xseq , function(x) hdi( Int + Slope*x ) )
}

p<-posterior_samples(ml2)
d2$Parasitescond<-d2$Parasites-d2$Age*mean(p[,"b_Age"])-d2$Sex*mean(p[,"b_Sex"])
Int<-p[,"b_Intercept"]
Slope<-p[,"b_PDSTotal"]
outP<-data.frame(V=0,I=mean(Int),S=mean(Slope),IL=hdi(Int)[1],IH=hdi(Int)[2],SL=hdi(Slope)[1],SH=hdi(Slope)[2])
mu.ciP <- list(sapply( xseq , function(x) hdi( Int + Slope*x ) ))
p<-posterior_samples(ml3)
for(i in 1:3){
  Int<-p[,"b_Intercept"]+p[,paste0("r_Village[",i,",Intercept]")]
  Slope<-p[,"b_PDSTotal"]+p[,paste0("r_Village[",i,",PDSTotal]")]
  outP<-rbind(outP,data.frame(V=i,I=mean(Int),S=mean(Slope),IL=hdi(Int)[1],IH=hdi(Int)[2],SL=hdi(Slope)[1],SH=hdi(Slope)[2]))
  mu.ciP[[i+1]] <- sapply( xseq , function(x) hdi( Int + Slope*x ) )
}

load("modT11.Rdata")
load("mod11MI.Rdata")
#Path plot
cols<-c("#ffbf47","#dd6031","#0262a2","#4daa57","#311e10")
threecols<-cols[c(3,5,1)]

boxes1<-data.frame(
  names=c("Inflam","HHInflam","ViInflam","PDSTotal","HHPDSTotal","ViPDSTotal"),
  xs=c(0.4,0.4,0.4,-0.4,-0.4,-0.4),
  ys=c(-0.8,0,0.8,-0.8,0,0.8),
  labpos=c(0.5,0.5,0.5,0.5,0.5,0.5),
  labels=c("Inflam","Household\nInflam","Community\nInflam","Disgust","Household\nDisgust","Community\nDisgust"),stringsAsFactors = FALSE)

boxes3<-data.frame(
  names=c("Parasites","HHParasites","ViParasites","PDSTotal","HHPDSTotal","ViPDSTotal"),
  xs=c(0.4,0.4,0.4,-0.4,-0.4,-0.4),
  ys=c(-0.8,0,0.8,-0.8,0,0.8),
  labpos=c(0.5,0.5,0.5,0.5,0.5,0.5),
  labels=c("Parasites","Household\nParasites","Community\nParasites","Disgust","Household\nDisgust","Community\nDisgust"),stringsAsFactors = FALSE)


#Complete model plots

tiff("Figure 2.tif",width=3600,height=1300,res=400,pointsize=11,compression="lzw")
layout(matrix(1:3,nrow=1),widths=c(1,1,0.7))
par(mar=c(0,0,0,0))
pathplot(modT11,boxes1,"A")
pathplot(modTP11,boxes3,"B")


#Parameters for component models
posts<-list(
  p1<-posterior_samples(mod11MI),
  p2<-posterior_samples(modF11MI),
  p3<-posterior_samples(modV11MI),
  p4<-posterior_samples(modP11MI),
  p5<-posterior_samples(modPF11MI),
  p6<-posterior_samples(modPV11MI))

par(mar=c(4,3,3,1))
cols<-c("#ffbf47","#dd6031","#0262a2","#4daa57","#311e10")[c(3,1,4)]
plot(0,0,xlim=c(-0.8,0.8),ylim=c(0.5,6.5),type="n",xlab="Posterior",yaxt="n",ylab=NA)
abline(v=0,lty=2)
off<-0.15
for(i in 1:6){
  pInd<-posts[[i]][,7]
  pIHH<-posts[[i]][,8]
  pHH<-posts[[i]][,19]
  lines(hdi(pInd,0.95),c(i+off,i+off),col=cols[1])
  lines(hdi(pInd,0.80),c(i+off,i+off),col=cols[1],lwd=3)
  lines(hdi(pInd,0.25),c(i+off,i+off),col=cols[1],lwd=7)
  lines(hdi(pIHH,0.95),c(i,i),col=cols[3])
  lines(hdi(pIHH,0.80),c(i,i),col=cols[3],lwd=3)
  lines(hdi(pIHH,0.25),c(i,i),col=cols[3],lwd=7)
  lines(hdi(pHH,0.95),c(i-off,i-off),col=cols[2])
  lines(hdi(pHH,0.80),c(i-off,i-off),col=cols[2],lwd=3)
  lines(hdi(pHH,0.25),c(i-off,i-off),col=cols[2],lwd=7)
}
abline(h=3.5,lty=3)
axis(2,at=1:6,labels=c("C1","C2","C3","C1","C2","C3"),las=2,tick=FALSE)
axis(2,at=c(2,5),labels=c("Inflam","Parasites"),line=-2,tick=FALSE,hadj=0.5)
box()
mtext("C",side=3,line=0.5,at=-1.2,cex=1.45)
dev.off()


#By Village plot
cols<-c("#ffbf47","#dd6031","#0262a2","#4daa57","#311e10")[c(1,3,4)]

pchs<-c(19,17,15)

tiff("Figure 3.tif",width=2500,height=1300,res=350,pointsize=11,compression="lzw")
par(mfrow=c(1,2),mar=c(5,5,1,1))
plot(Inflamcond~PDSTotal,data=d2,col=col.alpha(cols[Village],0.8),pch=pchs[Village],cex=2,xlab="Disgust",ylab="Inflammation")

library(rethinking)
shade(mu.ci[[1]],xseq,col=col.alpha("black",0.3))
shade(mu.ci[[2]],xseq,col=col.alpha(cols[1],0.3))
shade(mu.ci[[3]],xseq,col=col.alpha(cols[2],0.3))
shade(mu.ci[[4]],xseq,col=col.alpha(cols[3],0.3))

abline(out[1,2],out[1,3],col="black",lwd=2,lty=2)
abline(out[2,2],out[2,3],col=cols[1],lwd=2)
abline(out[3,2],out[3,3],col=cols[2],lwd=2)
abline(out[4,2],out[4,3],col=cols[3],lwd=2)

plot(Parasitescond~PDSCont,data=d2,col=col.alpha(cols[Village],0.8),pch=pchs[Village],cex=2,xlab="Disgust",ylab="Parasites")

shade(mu.ciP[[1]],xseq,col=col.alpha("black",0.3))
shade(mu.ciP[[2]],xseq,col=col.alpha(cols[1],0.3))
shade(mu.ciP[[3]],xseq,col=col.alpha(cols[2],0.3))
shade(mu.ciP[[4]],xseq,col=col.alpha(cols[3],0.3))

abline(outP[1,2],outP[1,3],col="black",lwd=2,lty=2)
abline(outP[2,2],outP[2,3],col=cols[1],lwd=2)
abline(outP[3,2],outP[3,3],col=cols[2],lwd=2)
abline(outP[4,2],outP[4,3],col=cols[3],lwd=2)
mtext(c("A","C","E"),side=2, line=-2,at=c(0.95,0.66,0.33),outer=TRUE,las=2,cex=2)
mtext(c("B","D","F"),side=2, line=-25,at=c(0.95,0.66,0.33),outer=TRUE,las=2,cex=2)

dev.off()