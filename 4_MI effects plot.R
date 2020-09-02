library(brms)
library(HDInterval)

col.alpha <- function (acol, alpha = 0.5){  # color funtion for plotting
  acol <- col2rgb(acol)
  acol.red <- acol["red",]/255
  acol.green <- acol["green",]/255
  acol.blue <- acol["blue",]/255
  acol <- mapply(function(red, green, blue, alphas) rgb(red, green, blue, alphas), acol.red, acol.green, acol.blue, alpha)
  as.character(acol)
}

#posterior parameter plot
cols<-c("#ffbf47","#dd6031","#0262a2","#4daa57","#311e10")[c(1,3)]

getdens<-function(model,model2){  
  p<-posterior_samples(model)
  p2<-posterior_samples(model2)
  n1<-strsplit(names(p[1]),"_")[[1]][2]
  n2<-strsplit(names(p[2]),"_")[[1]][2]
  n3<-strsplit(names(p2[1]),"_")[[1]][2]
  n4<-strsplit(names(p2[2]),"_")[[1]][2]
  
  d1HSOL<-density(c(p[,paste0("b_",n1,"_HOUSE")],p2[,paste0("b_",n3,"_HOUSE")]))
  d1MSOL<-density(c(p[,paste0("b_",n1,"_MSOL")],p2[,paste0("b_",n3,"_MSOL")]))
  d1TSOL<-density(c(p[,paste0("b_",n1,"_TSOL")],p2[,paste0("b_",n3,"_TSOL")]))
  d1HHHSOL<-density(c(p[,paste0("b_HH",n1,"_HOUSE")],p2[,paste0("b_HH",n3,"_HOUSE")]))
  d1HHMSOL<-density(c(p[,paste0("b_HH",n1,"_MSOL")],p2[,paste0("b_HH",n3,"_MSOL")]))
  d1HHTSOL<-density(c(p[,paste0("b_HH",n1,"_TSOL")],p2[,paste0("b_HH",n3,"_TSOL")]))
  d2HSOL<-density(c(p[,paste0("b_",n2,"_HOUSE")],p2[,paste0("b_",n4,"_HOUSE")]))
  d2MSOL<-density(c(p[,paste0("b_",n2,"_MSOL")],p2[,paste0("b_",n4,"_MSOL")]))
  d2TSOL<-density(c(p[,paste0("b_",n2,"_TSOL")],p2[,paste0("b_",n4,"_TSOL")]))
  d2HHHSOL<-density(c(p[,paste0("b_HH",n2,"_HOUSE")],p2[,paste0("b_HH",n4,"_HOUSE")]))
  d2HHMSOL<-density(c(p[,paste0("b_HH",n2,"_MSOL")],p2[,paste0("b_HH",n4,"_MSOL")]))
  d2HHTSOL<-density(c(p[,paste0("b_HH",n2,"_TSOL")],p2[,paste0("b_HH",n4,"_TSOL")]))
  
  d1HSOLhdi<-hdi(c(p[,paste0("b_",n1,"_HOUSE")],p2[,paste0("b_",n3,"_HOUSE")]))
  d1MSOLhdi<-hdi(c(p[,paste0("b_",n1,"_MSOL")],p2[,paste0("b_",n3,"_MSOL")]))
  d1TSOLhdi<-hdi(c(p[,paste0("b_",n1,"_TSOL")],p2[,paste0("b_",n3,"_TSOL")]))
  d1HHHSOLhdi<-hdi(c(p[,paste0("b_HH",n1,"_HOUSE")],p2[,paste0("b_HH",n3,"_HOUSE")]))
  d1HHMSOLhdi<-hdi(c(p[,paste0("b_HH",n1,"_MSOL")],p2[,paste0("b_HH",n3,"_MSOL")]))
  d1HHTSOLhdi<-hdi(c(p[,paste0("b_HH",n1,"_TSOL")],p2[,paste0("b_HH",n3,"_TSOL")]))
  d2HSOLhdi<-hdi(c(p[,paste0("b_",n2,"_HOUSE")],p2[,paste0("b_",n4,"_HOUSE")]))
  d2MSOLhdi<-hdi(c(p[,paste0("b_",n2,"_MSOL")],p2[,paste0("b_",n4,"_MSOL")]))
  d2TSOLhdi<-hdi(c(p[,paste0("b_",n2,"_TSOL")],p2[,paste0("b_",n4,"_TSOL")]))
  d2HHHSOLhdi<-hdi(c(p[,paste0("b_HH",n2,"_HOUSE")],p2[,paste0("b_HH",n4,"_HOUSE")]))
  d2HHMSOLhdi<-hdi(c(p[,paste0("b_HH",n2,"_MSOL")],p2[,paste0("b_HH",n4,"_MSOL")]))
  d2HHTSOLhdi<-hdi(c(p[,paste0("b_HH",n2,"_TSOL")],p2[,paste0("b_HH",n4,"_TSOL")]))
  return(list(den=list(d1HSOL,d1MSOL,d1TSOL,d1HHHSOL,d1HHMSOL,d1HHTSOL,d2HSOL,d2MSOL,d2TSOL,d2HHHSOL,d2HHMSOL,d2HHTSOL),hdi=list(d1HSOLhdi,d1MSOLhdi,d1TSOLhdi,d1HHHSOLhdi,d1HHMSOLhdi,d1HHTSOLhdi,d2HSOLhdi,d2MSOLhdi,d2TSOLhdi,d2HHHSOLhdi,d2HHMSOLhdi,d2HHTSOLhdi),names=c(n1,n2,n3,n4)))
}

load("mod11MI.Rdata")
dhC<-getdens(mod11MI,modP11MI)
dhF<-getdens(modF11MI,modPF11MI)
dhP<-getdens(modP11MI,modPF11MI)
dhI<-getdens(mod11MI,modF11MI)
  
denspoly<-function(den,hdi,col,shift){
  y<-den$y[den$x >= hdi[1] & den$x <= hdi[2]]
  y<-y/max(y) + shift
  x<-den$x[den$x >= hdi[1] & den$x <= hdi[2]]
  lines(x,y)
  polygon(c(x[1],x,x[length(x)]),c(shift,y,shift),col=col.alpha(col,0.7),border=NA)
}

tiff("MI Posterior plot.tif",width=3000,height=1400,res = 300,pointsize=13,compression="lzw")
par(mfrow=c(1,3),mar=c(5,4,3,1))
names<-c("Inflammation","Parasites","Contagion","Food")
plot(0,0,xlim=c(-1,1),ylim=c(0,4.5),type="n",main="MSOL",yaxt="n",ylab=NA,xlab=NA)
  shifts<-c(0,1.1,2.4,3.5)
  #MSOL - Inflam
  denspoly(dhI$den[[2]],dhI$hdi[[2]],cols[1],shifts[1])
  denspoly(dhI$den[[5]],dhI$hdi[[5]],cols[2],shifts[1])
  #MSOL - Parasites
  denspoly(dhP$den[[2]],dhP$hdi[[2]],cols[1],shifts[2])
  denspoly(dhP$den[[5]],dhP$hdi[[5]],cols[2],shifts[2])
  #MSOL - PDS
  denspoly(dhC$den[[8]],dhC$hdi[[8]],cols[1],shifts[3])
  denspoly(dhC$den[[11]],dhC$hdi[[11]],cols[2],shifts[3])
  #MSOL - Food
  denspoly(dhF$den[[8]],dhF$hdi[[8]],cols[1],shifts[4])
  denspoly(dhF$den[[11]],dhF$hdi[[11]],cols[2],shifts[4])
  abline(v=0,lty=2)

  axis(2,at=shifts+0.5,labels=names,crt=2,tick=FALSE)

plot(0,0,xlim=c(-1,1),ylim=c(0,4.5),type="n",main="HSOL",yaxt="n",ylab=NA,xlab="Posterior Parameter")
  #HSOL - Inflam
  denspoly(dhI$den[[1]],dhI$hdi[[1]],cols[1],shifts[1])
  denspoly(dhI$den[[4]],dhI$hdi[[4]],cols[2],shifts[1])
  #HSOL - Parasites
  denspoly(dhP$den[[1]],dhP$hdi[[1]],cols[1],shifts[2])
  denspoly(dhP$den[[4]],dhP$hdi[[4]],cols[2],shifts[2])
  #HSOL - PDS
  denspoly(dhC$den[[7]],dhC$hdi[[7]],cols[1],shifts[3])
  denspoly(dhC$den[[10]],dhC$hdi[[10]],cols[2],shifts[3])
  #HSOL - Food
  denspoly(dhF$den[[7]],dhF$hdi[[7]],cols[1],shifts[4])
  denspoly(dhF$den[[10]],dhF$hdi[[10]],cols[2],shifts[4])
  abline(v=0,lty=2)
  axis(2,at=shifts+0.5,labels=names,crt=2,tick=FALSE)
  
plot(0,0,xlim=c(-1,1),ylim=c(0,4.5),type="n",main="TSOL",yaxt="n",ylab=NA,xlab=NA)
  #HSOL - Inflam
  denspoly(dhI$den[[3]],dhI$hdi[[3]],cols[1],shifts[1])
  denspoly(dhI$den[[6]],dhI$hdi[[6]],cols[2],shifts[1])
  #HSOL - Parasites
  denspoly(dhP$den[[3]],dhP$hdi[[3]],cols[1],shifts[2])
  denspoly(dhP$den[[6]],dhP$hdi[[6]],cols[2],shifts[2])
  #HSOL - PDS
  denspoly(dhC$den[[9]],dhC$hdi[[9]],cols[1],shifts[3])
  denspoly(dhC$den[[12]],dhC$hdi[[12]],cols[2],shifts[3])
  #HSOL - Food
  denspoly(dhF$den[[9]],dhF$hdi[[9]],cols[1],shifts[4])
  denspoly(dhF$den[[12]],dhF$hdi[[12]],cols[2],shifts[4])
  abline(v=0,lty=2)
  axis(2,at=shifts+0.5,labels=names,crt=2,tick=FALSE)
  
dev.off()
