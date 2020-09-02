#Plots of effects by village

cols<-c("#ffbf47","#dd6031","#0262a2","#4daa57","#311e10")[c(1,3,4)]
col.alpha <- function (acol, alpha = 0.5){  # color funtion for plotting
  acol <- col2rgb(acol)
  acol.red <- acol["red",]/255
  acol.green <- acol["green",]/255
  acol.blue <- acol["blue",]/255
  acol <- mapply(function(red, green, blue, alphas) rgb(red, green, blue, alphas), acol.red, acol.green, acol.blue, alpha)
  as.character(acol)
}

partline<-function(m,...){
  nd<-data.frame(range(m$model[,2]))
  names(nd)<-names(m$model)[2]
  p<-predict(m,newdata=nd)
  lines(nd[,1],p,...)
  
}

#Use the unimputed dataset from code #1
d<-d2
#Replace Inflam and Parasites with the mean of the imputed datasets for each individual. Only has an effect on those with missing values.
avp<-rowMeans(data.frame(lapply(dsets,function(ds) ds$Parasites)))
avi<-rowMeans(data.frame(lapply(dsets,function(ds) ds$Inflam)))
d$Inflam<-avi
d$Parasites<-avp

pchs<-c(19,17,15)

tiff("Village plot.tif",width=2500,height=2500,res=350,pointsize=11,compression="lzw")
par(mfrow=c(2,2),mar=c(5,5,1,1))
plot(Inflam~PDSCont,data=d,col=col.alpha(cols[Village],0.8),pch=pchs[Village],cex=2,xlab="Contagion Disgust",ylab="Inflammation")
ml0<-lm(Inflam~PDSCont,data=d)
ml1<-lm(Inflam~PDSCont,data=d[d$Village==1,])
ml2<-lm(Inflam~PDSCont,data=d[d$Village==2,])
ml3<-lm(Inflam~PDSCont,data=d[d$Village==3,])

partline(ml0,col="black",lwd=2,lty=2)
partline(ml1,col=cols[1],lwd=2)
partline(ml2,col=cols[2],lwd=2)
partline(ml3,col=cols[3],lwd=2)

plot(Inflam~PDSFood,data=d,col=col.alpha(cols[Village],0.8),pch=pchs[Village],cex=2,xlab="Food Disgust",ylab="Inflammation")
ml0<-lm(Inflam~PDSFood,data=d)
ml1<-lm(Inflam~PDSFood,data=d[d$Village==1,])
ml2<-lm(Inflam~PDSFood,data=d[d$Village==2,])
ml3<-lm(Inflam~PDSFood,data=d[d$Village==3,])

partline(ml0,col="black",lwd=2,lty=2)
partline(ml1,col=cols[1],lwd=2)
partline(ml2,col=cols[2],lwd=2)
partline(ml3,col=cols[3],lwd=2)

plot(Parasites~PDSCont,data=d,col=col.alpha(cols[Village],0.8),pch=pchs[Village],cex=2,xlab="Contagion Disgust",ylab="Parasites")
ml0<-lm(Parasites~PDSCont,data=d)
ml1<-lm(Parasites~PDSCont,data=d[d$Village==1,])
ml2<-lm(Parasites~PDSCont,data=d[d$Village==2,])
ml3<-lm(Parasites~PDSCont,data=d[d$Village==3,])

partline(ml0,col="black",lwd=2,lty=2)
partline(ml1,col=cols[1],lwd=2)
partline(ml2,col=cols[2],lwd=2)
partline(ml3,col=cols[3],lwd=2)

plot(Parasites~PDSFood,data=d,col=col.alpha(cols[Village],0.8),pch=pchs[Village],cex=2,xlab="Food Disgust",ylab="Parasites")
ml0<-lm(Parasites~PDSFood,data=d)
ml1<-lm(Parasites~PDSFood,data=d[d$Village==1,])
ml2<-lm(Parasites~PDSFood,data=d[d$Village==2,])
ml3<-lm(Parasites~PDSFood,data=d[d$Village==3,])

partline(ml0,col="black",lwd=2,lty=2)
partline(ml1,col=cols[1],lwd=2)
partline(ml2,col=cols[2],lwd=2)
partline(ml3,col=cols[3],lwd=2)
dev.off()
