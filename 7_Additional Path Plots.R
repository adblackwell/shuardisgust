#Path plots

library(mice)
library(brms)
library(HDInterval)

#plotting of path models

cols<-c("#ffbf47","#dd6031","#0262a2","#4daa57","#311e10")
threecols<-cols[c(3,5,1)]


#Single component plots

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

#Complete model reversed plot

load("modT11b.Rdata")
tiff("complete Modelsb_1F.tif",width=3200,height=1300,res=400,pointsize=11,compression="lzw")
par(mfrow=c(1,2),mar=c(0,0,0,0))
pathplot(modT11b,boxes1,"A")
pathplot(modTP11b,boxes3,"B")
dev.off()

#large mediation models
load("modTI.Rdata")
load("modTP.Rdata")

tx <- -0.8
tiff("Mediation_1F.tif",width=1700,height=3000,res=360,pointsize=8,compression="lzw")
par(mfrow=c(1,2),mar=c(0,0,0,0))
#Inflam-Total
boxesT<-boxes1
boxesT$ys<-boxesT$ys+2
boxesT$xs<-boxesT$xs-0.15
pathplot(mod1T,boxesT,ylim=c(-4.8,2.4))

boxesT$ys<-boxesT$ys+0.2
pathplot(mod3T,boxesT,new=FALSE)
boxesT$xs<-boxesT$xs+0.35
pathplot(mod2,boxesT,new=FALSE)
boxesT$xs<-boxesT$xs+0.55

boxesT<-boxes1
pathplot(mod6T,boxesT,new=FALSE)
boxesT$ys<-boxesT$ys+0.4
pathplot(mod7T,boxesT,new=FALSE)

boxesT<-boxes1
boxesT$ys<-boxesT$ys-1.6
pathplot(mod8T,boxesT,new=FALSE)

boxesT<-boxes1
boxesT$ys<-boxesT$ys-4.0
pathplot(modT11,boxesT,new=FALSE)

l<-1
text(tx,2.4,LETTERS[l],cex=2)
text(tx,0.6,LETTERS[l+2],cex=2)
text(tx,-1.4,LETTERS[l+4],cex=2)
text(tx,-3.0,LETTERS[l+6],cex=2)
l<-l+1


#Parasites-Total
boxesT<-boxes3
boxesT$ys<-boxesT$ys+2
boxesT$xs<-boxesT$xs-0.15
pathplot(modP1T,boxesT,ylim=c(-4.8,2.4))

boxesT$ys<-boxesT$ys+0.2
pathplot(modP3T,boxesT,new=FALSE)
boxesT$xs<-boxesT$xs+0.35
pathplot(modP2,boxesT,new=FALSE)
boxesT$xs<-boxesT$xs+0.55

boxesT<-boxes3
pathplot(modP6T,boxesT,new=FALSE)
boxesT$ys<-boxesT$ys+0.4
pathplot(modP7T,boxesT,new=FALSE)

boxesT<-boxes3
boxesT$ys<-boxesT$ys-1.6
pathplot(modP8T,boxesT,new=FALSE)

boxesT<-boxes3
boxesT$ys<-boxesT$ys-4.0
pathplot(modTP11,boxesT,new=FALSE)

text(tx,2.4,LETTERS[l],cex=2)
text(tx,0.6,LETTERS[l+2],cex=2)
text(tx,-1.4,LETTERS[l+4],cex=2)
text(tx,-3.0,LETTERS[l+6],cex=2)
l<-l+1

dev.off()

#Market Integration Path Plots
boxes1MI<-data.frame(
  names=c("Inflam","HHInflam","ViInflam","PDSTotal","HHPDSTotal","ViPDSTotal","MSOL","HOUSE","TSOL","MSOL2","HOUSE2","TSOL2"),
  xs=c(0.5,0.5,0.5,-0.5,-0.5,-0.5,1.5,1.5,1.5,-1.5,-1.5,-1.5),
  ys=c(-0.8,0,0.8,-0.8,0,0.8,0.4,-0.4,-1.2,0.4,-0.4,-1.2),
  labpos=c(0.5,0.6,0.4,0.5,0.6,0.4,0.3,0.3,0.3,0.3,0.3,0.3),
  labels=c("Inflam","Household\nInflam","Community\nInflam","Disgust","Household\nDisgust","Community\nDisgust","MSOL","HSOL","TSOL","MSOL","HSOL","TSOL"),stringsAsFactors = FALSE)

boxes3MI<-data.frame(
  names=c("Parasites","HHParasites","ViParasites","PDSTotal","HHPDSTotal","ViPDSTotal","MSOL","HOUSE","TSOL","MSOL2","HOUSE2","TSOL2"),
  xs=c(0.5,0.5,0.5,-0.5,-0.5,-0.5,1.5,1.5,1.5,-1.5,-1.5,-1.5),
  ys=c(-0.8,0,0.8,-0.8,0,0.8,0.4,-0.4,-1.2,0.4,-0.4,-1.2),
  labpos=c(0.5,0.6,0.4,0.5,0.6,0.4,0.3,0.3,0.3,0.3,0.3,0.3),
  labels=c("Parasites","Household\nParasites","Community\nParasites","Disgust","Household\nDisgust","Community\nDisgust","MSOL","HSOL","TSOL","MSOL","HSOL","TSOL"),stringsAsFactors = FALSE)


load("modT11MI.Rdata")
tiff("MI Models_1F.tif",width=1850,height=3000,res=350,pointsize=11,compression="lzw")
par(mfrow=c(2,1),mar=c(0,0,0,0))
tx <- -1.5
ty <- 0.8
#Inflam-Total
pathplot(modT11MI,boxes1MI,xlim=c(-1.8,1.8),ylim=c(-1.5,1),split=TRUE)
l<-1
text(tx,ty,LETTERS[l],cex=3); l <- l +1
pathplot(modTP11MI,boxes3MI,xlim=c(-1.8,1.8),ylim=c(-1.5,1),split=TRUE)
text(tx,ty,LETTERS[l],cex=3); l <- l +1
dev.off()