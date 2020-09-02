#Models to examine variance components
#Needs dsets data from code #1
ParasitesV<-brm_multiple(Parasites ~ 1 + (1|Family) + (1|Village), data=dsets, family=gaussian(), chains=4,cores=4,control = list(adapt_delta = 0.99, max_treedepth=12))

InflamV<-brm_multiple(Inflam ~ 1 + (1|Family) + (1|Village), data=dsets, family=gaussian(), chains=4,cores=4,control = list(adapt_delta = 0.99, max_treedepth=12))

PDSContModV<-brm(PDSCont ~ 1 + (1|Family) + (1|Village), data=d2, family=gaussian(), chains=4,cores=4,control = list(adapt_delta = 0.99, max_treedepth=12))

PDSFoodModV<-brm(PDSFood ~ 1 + (1|Family) + (1|Village), data=d2, family=gaussian(), chains=4,cores=4,control = list(adapt_delta = 0.99, max_treedepth=12))

MSOLModV<-brm(MSOL ~ 1 + (1|Family) + (1|Village), data=d2, family=gaussian(), chains=4,cores=4,control = list(adapt_delta = 0.99, max_treedepth=12))

HOUSEModV<-brm(HOUSE ~ 1 + (1|Family) + (1|Village), data=d2, family=gaussian(), chains=4,cores=4,control = list(adapt_delta = 0.99, max_treedepth=12))

TSOLModV<-brm(TSOL ~ 1 + (1|Family) + (1|Village), data=d2, family=gaussian(), chains=4,cores=4,control = list(adapt_delta = 0.99, max_treedepth=12))

AgeModV<-brm(Age ~ 1 + (1|Family) + (1|Village), data=d2, family=gaussian(), chains=4,cores=4,control = list(adapt_delta = 0.99, max_treedepth=12))

SexModV<-brm(Age ~ 1 + (1|Family) + (1|Village), data=d2, family=gaussian(), chains=4,cores=4,control = list(adapt_delta = 0.99, max_treedepth=12))

save.image("Variancemodels.RData")

#Uncomment if starting with models already run
#load("Variancemodels.RData")
library(brms)
Sds<-rbind(
  unlist(VarCorr(PDSFoodModV))[c(5,1,9)],
  unlist(VarCorr(PDSContModV))[c(5,1,9)],
  unlist(VarCorr(ParasitesV))[c(5,1,9)],
  unlist(VarCorr(InflamV))[c(5,1,9)],
  unlist(VarCorr(MSOLModV))[c(5,1,9)],
  unlist(VarCorr(HOUSEModV))[c(5,1,9)],
  unlist(VarCorr(TSOLModV))[c(5,1,9)]

)
Vars<-apply(Sds,c(1,2),function(x) x^2)
rowSums(Vars)
#Table of variance components
Comps<-Vars/rowSums(Vars)
          
varlabs<-c("Contagion Disgust","Food Disgust","Parasites","Inflammation","MSOL","HSOL","TSOL")

#make table
library(flextable)
library(officer)
colnames(Comps)<-c("Community","Household","Individual")
rownames(Comps)<-varlabs
CompsT<-as.data.frame(Comps[c(4,1,2,3,7,5,6),])
CompsT$Variable<-row.names(CompsT)
ft<-flextable(CompsT,col_keys = c("Variable","Village","(officer)","Individual"))
ft<-colformat_num(ft,j=c("Village","Family","Individual"))
ft<-autofit(ft)
ft<-font(ft,fontname="Calibri",part="all")
ft<-fontsize(ft,size=11)
ft<-bold(ft,part="header")
ft
save_as_docx(ft, path = "Variance Table.docx")

#Plot of variance components and village means
tiff("VarianceParams.tif",width=3000,height=1200,res=300,pointsize=12,compression="lzw")
layout(matrix(c(4,1,5,2,3),nrow=1),widths=c(0.1,1,0.1,0.85,0.15))
par(mar=c(5,1,5,1))
#library(wesanderson)
#cols<-wes_palette("Rushmore1",5)[c(2,3,4)]
#<-c("#157A6E","#F6AE2D","#A63446")
cols<-c("#ffbf47","#dd6031","#0262a2","#4daa57","#311e10")[c(4,1,3)]
colnames(Comps)<-c("Community","Household","Individual")
rownames(Comps)<-varlabs
Comps2<-Comps[c(4,1,2,3,7,5,6),]
mosaicplot(Comps2,color=cols,dir=c("h","v"),las=1,cex=1,main=NA)
axis(1,at=c(0.23,0.61,0.99),labels=c("0%","50%","100%"),line=0.5)
axis(1,at=0.61,labels=c("Proportion of Variance"),line=2.5,tick=FALSE)

par(mar=c(5,7,6,1))
ms<-aggregate(cbind(PDSCont,PDSFood,Parasites,Inflam,MSOL,HOUSE,TSOL)~Village,data=dlong,FUN=mean)
ms<-as.matrix(ms[,-1])
colnames(ms)<-varlabs
ms<-ms[,rev(row.names(Comps2))]
cols2<-colorRampPalette(rev(cols[c(2,1,3)]))(9)
image(ms,zlim=c(-1.1,1.1),col=cols2,xaxt="n",yaxt="n",main="Community")
axis(2,at=seq(0,1,length=ncol(ms)),labels=colnames(ms),las=2,tick=FALSE)
axis(3,at=seq(0,1,length=3),labels=c(1,2,3),las=1,tick=FALSE,line=0)
grid(3,7,lty=1,col="black")
box()
par(mar=c(7,2,8,1))
image(y=seq(-1,1,length=10),x=c(0,1),z=matrix(seq(-1,1,length=10),nrow=1),zlim=c(-1.1,1.1),col=cols2,xaxt="n",yaxt="n",main="Z-score",xlab=NA)
axis(2,at=seq(-1,1,length=6),labels=seq(-1,1,length=6),las=2)
box()

par(mar=c(1,1,1,1))
plot(0,1,xlim=c(0,1),ylim=c(0,1),type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA,bty="n")
text(0.5,0.9,"A",cex=3)
plot(0,1,xlim=c(0,1),ylim=c(0,1),type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA,bty="n")
text(0.5,0.9,"B",cex=3)
dev.off()


#Similarity Plot
SimPlot<-function(v,cols,cols2,g1,g2,line1,line2,data=NA,...){
  if(!class(data)=="data.frame"){
    mats<-lapply(data,function(d) 
    {d<-d[order(d[,g1],d[,g2],d$Age),]
    outer(d[,v],d[,v],FUN="-") })
    vals<-lapply(data,function(d) 
    {d<-d[order(d[,g1],d[,g2],d$Age),]
    d[,v] })
    cc<-Reduce("+",mats)/length(mats)
    vals<-Reduce("+",vals)/length(vals)
    g1<-data[[1]][,g1]
    g2<-data[[1]][,g2]
  } else {
    data<-data[order(data[,g1],data[,g2],data$Age),]
    vals<-data[,v]
    cc<-outer(vals,vals,FUN="-")
    g1<-data[,g1]
    g2<-data[,g2]
  }
  xy<-seq(0,nrow(cc),by=1)
  par(mar=c(1,1,4,1))
  image(x=xy,y=xy,abs(cc), col = cols,xaxt="n",yaxt="n",xlab=NA,ylab=NA,...)
  lg<-c(0,cumsum(table(g1)))
  for(i in 1:(length(lg)-1)){
    polygon(c(lg[i],lg[i],lg[i+1],lg[i+1]),c(lg[i],lg[i+1],lg[i+1],lg[i]),lwd=2,border=line1)
  }
  lg<-c(0,cumsum(aggregate(rep(1,length(g2))~g2,FUN=sum)[,2]))
  for(i in 1:(length(lg)-1)){
    polygon(c(lg[i],lg[i],lg[i+1],lg[i+1]),c(lg[i],lg[i+1],lg[i+1],lg[i]),lwd=2,border=line2)
  }
  #ValPlot
  par(mar=c(1,1,0,1))
  vals<-as.matrix(vals,ncol=1)
  xy<-seq(0,length(vals),by=1)
  image(x=xy,y=c(1,2),vals, col = cols2,xaxt="n",yaxt="n",xlab=NA,ylab=NA,...)
}


#Similarity Plot
tiff("Similarity Matrix.tif",width=3000,height=2000,res=300,pointsize=13,compression="lzw")
laym<-matrix(1:8,byrow=FALSE,ncol=4)
laym<-rbind(laym,matrix(9:16,byrow=FALSE,ncol=4))
layout(laym,heights=rep(c(1,0.1),4))

#cols<-rev(brewer.pal(9,"Blues"))
cols3<-colorRampPalette(c(cols[3],"#FFFFFF"))(9)
cols2<-colorRampPalette(rev(cols[c(2,1,3)]))(9)
linecol<-cols[2]
SimPlot("Parasites",cols3,cols2,"Village","Family",linecol,linecol,main="Parasites",data=dsets)

SimPlot("Inflam",cols3,cols2,"Village","Family",linecol,linecol,main="Inflammation",data=dsets)

SimPlot("PDSCont",cols3,cols2,"Village","Family",linecol,linecol,main="Contagion Disgust",data=dsets)

SimPlot("PDSFood",cols3,cols2,"Village","Family",linecol,linecol,main="Food Disgust",data=dsets)

SimPlot("MSOL",cols3,cols2,"Village","Family",linecol,linecol,main="MSOL",data=dsets)

SimPlot("HOUSE",cols3,cols2,"Village","Family",linecol,linecol,main="HSOL",data=dsets)

SimPlot("TSOL",cols3,cols2,"Village","Family",linecol,linecol,main="TSOL",data=dsets)

dev.off()


#Compare villages
#simple Statistics for village comparisons
t(sapply(c("Inflam","Parasites","PDSCont","PDSFood","MSOL","HOUSE","TSOL","Age"), function(v) {
  m<-lm(formula(paste0(v,"~Village")),data=d2)
  a<-anova(m)
  c(a$`F value`[1],a$`Pr(>F)`[1])
}))

m<-glm(Sex~Village,data=d2,family=binomial)
a<-anova(m)

ms<-aggregate(~Village,data=dlong,FUN=mean)
chisq.test(d2$Sex,d2$Village)


#Multidimensional scaling of varianbles plotted with village
#Uses d2 dataset from code #1
ddist <- dist(d2[,c("Inflam","Parasites","PDSCont","PDSFood","MSOL","HOUSE","TSOL")]) # euclidean distances between the rows
fit <- cmdscale(ddist,eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS",col=c("red","blue","purple")[d2$Village])

cols<-c("#1B9E77","#D95F02","#7570B3")

col.alpha <- function (acol, alpha = 0.5){  # color funtion for plotting
  acol <- col2rgb(acol)
  acol.red <- acol["red",]/255
  acol.green <- acol["green",]/255
  acol.blue <- acol["blue",]/255
  acol <- mapply(function(red, green, blue, alphas) rgb(red, green, blue, alphas), acol.red, acol.green, acol.blue, alpha)
  as.character(acol)
}

library(MASS)
k<-kde2d(x[d2$Village==1],y[d2$Village==1],lims=c(-3,3,-3,3),n=50)
image(k,zlim=c(0.05,0.2),col=col.alpha(colorRampPalette(c("#FFFFFF", cols[1]))(12)[3:12]))
k<-kde2d(x[d2$Village==2],y[d2$Village==2],lims=c(-3,3,-3,3),n=50)
image(k,zlim=c(0.05,0.2),add=TRUE,col=col.alpha(colorRampPalette(c("#FFFFFF", cols[2]))(12)[3:12]))
k<-kde2d(x[d2$Village==3],y[d2$Village==3],lims=c(-3,3,-3,3),n=50)
image(k,zlim=c(0.05,0.2),add=TRUE,col=col.alpha(colorRampPalette(c("#FFFFFF", cols[3]))(12)[3:12]))
points(x,y,col=cols[d2$Village],pch=19)

xy<-cor(cbind(x,y,d2[,c("Inflam","Parasites","PDSCont","PDSFood","MSOL","HOUSE","TSOL")]),use="pairwise.complete.obs")[-c(1,2),c(1,2)]
yoff<-c(0,0,0.1,0,0,-0.1,0)
lab<-c("Inflam","Parasites","Cont Disgust","Food Disgust","MSOL","HSOL","TSOL")
for(i in 1:nrow(xy)){
  lines(c(0,xy[i,1]*3),c(0,xy[i,2]*3))
  text(xy[i,1]*3,xy[i,2]*3+yoff[i],lab[i])
}

