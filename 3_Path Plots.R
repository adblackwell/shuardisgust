#Path plots
col.alpha <- function (acol, alpha = 0.5){  # color funtion for plotting
  acol <- col2rgb(acol)
  acol.red <- acol["red",]/255
  acol.green <- acol["green",]/255
  acol.blue <- acol["blue",]/255
  acol <- mapply(function(red, green, blue, alphas) rgb(red, green, blue, alphas), acol.red, acol.green, acol.blue, alpha)
  as.character(acol)
}

beta<-function(beta){
  format(round(beta[1],2),nsmall=2,digits=2)
}

betaCI<-function(beta){
  paste0(format(round(beta[1],2),nsmall=2,digits=2)," (",format(round(beta[2],2),nsmall=2,digits=2), " - ", format(round(beta[3],2),nsmall=2,digits=2),")")
}

betaCI2<-function(beta){
  paste0(format(round(beta[1],2),nsmall=2,digits=2),"\n(",format(round(beta[2],2),nsmall=2,digits=2), " - ", format(round(beta[3],2),nsmall=2,digits=2),")")
}

library(mice)
library(brms)
library(HDInterval)


#plotting of path models

cols<-c("#ffbf47","#dd6031","#0262a2","#4daa57","#311e10")
threecols<-cols[c(3,5,1)]


## Note: curvedarrow source code is From the 'diagram' package
##==============================================================================
# curvedarrow: Plot curved arrow at certain distance between two points
##==============================================================================

curvedarrow <- function(from, to, lwd=2, lty=1, lcol="black", arr.col=lcol, 
                        arr.pos=0.5, lab=NA, lab.pos=0.8, curve=1, dr=0.01, endhead=FALSE, segment = c(0,1),tpos=NULL, ...)   {
  require(shape)
  dpos  <- to-from
  angle <- atan(dpos[2]/dpos[1])*180/pi         # angle between both
  if (is.nan(angle)) return
  mid   <- 0.5*(to+from)                        # midpoint of ellipsoid arrow
  dst   <- dist(rbind(to, from))                # distance from-to
  ry    <- curve*dst                            # small radius of ellepsoid
  aFrom<-0                                      # angle to and from
  aTo<-pi
  if ( from[1] <= to[1]) {
    aFrom <- pi
    aTo <- 2*pi
  }
  
  if (segment [1] != 0)
    From <- segment[1] * aTo + (1-segment[1]) * aFrom
  else
    From <- aFrom
  
  if (segment [2] != 1)
    To <- segment[2] * aTo + (1-segment[2]) * aFrom
  else
    To <- aTo
  
  meanpi <- arr.pos * aTo + (1-arr.pos) * aFrom
  if (endhead) To <- meanpi
  meanpilab <- lab.pos * aTo + (1-lab.pos) * aFrom
  
  plotellipse(rx=dst/2,  ry=ry, mid=mid, angle=angle, from = From, to = To,
              lwd=lwd, lty=lty, lcol=lcol)
  ell <- getellipse(rx=dst/2, ry=ry, mid=mid, angle=angle,
                    from=1.001*meanpi, to=0.999*meanpi, dr= 0.002)       #Changed from -0.002
  Arrows(ell[1,1], ell[1,2], ell[nrow(ell),1], ell[nrow(ell),2],
         code=1, lcol=lcol, arr.col=arr.col, ...)
  ell2 <- getellipse(rx=dst/2, ry=ry, mid=mid, angle=angle,
                     from=1.001*meanpilab, to=0.999*meanpilab, dr= 0.002)
  text(ell2[nrow(ell2),1], ell2[nrow(ell2),2],lab,pos=tpos)
  curvedarrow <- c(ell[nrow(ell),1], ell[nrow(ell),2])
}


pathplot<-function(model,boxes,letter=NA,xlim=c(-1,1),ylim=c(-1,1),new=TRUE,split=FALSE){
  if(class(model)[1] %in% c("brmsfit","brmsfit_multiple")) {
    segments<-data.frame(summary(model)$fixed,stringsAsFactors = FALSE)
    bps<-apply(posterior_samples(model)[,1:nrow(segments)],2,function(m) sum(m>0)/length(m))
    segments$bps<-sapply(bps,function(n) min(n,1-n))
    if(!all(grepl("_",row.names(segments)))) row.names(segments)<-paste0(all.vars(model$formula[[1]])[1],"_",row.names(segments))
  } else segments<-Reduce(rbind,lapply(model,function(x) {
    segments<-data.frame(summary(x)$fixed,stringsAsFactors = FALSE)
    bps<-apply(posterior_samples(x)[,1:nrow(segments)],2,function(m) sum(m>0)/length(m))
    segments$bps<-sapply(bps,function(n) min(n,1-n))
    if(!all(grepl("_",row.names(segments)))) row.names(segments)<-paste0(all.vars(x$formula[[1]])[1],"_",row.names(segments))
    segments}))
  
  Pair<-data.frame(Reduce(rbind,strsplit(row.names(segments),"_")),stringsAsFactors = FALSE)
  names(Pair)<-c("to","from")
  segments<-cbind(segments,Pair)
  if(split) {
    segments$from[(segments$from %in% c("MSOL","TSOL","HOUSE")) & (segments$to %in% c("PDSCont","HHPDSCont","PDSFood","HHPDSFood"))]<-paste0(segments$from[(segments$from %in% c("MSOL","TSOL","HOUSE")) & (segments$to %in% c("PDSCont","HHPDSCont","PDSFood","HHPDSFood"))],"2")
  }
  segs<-segments[(segments$to %in% boxes$names) & (segments$from %in% boxes$names),]
  boxes<-boxes[boxes$names %in% c(segments$to,segments$from),]
  segs$curve<-(grepl("Vi",segs$from) & !grepl("HH|Vi",segs$to))*0.2*sign(boxes$xs[sapply(segs$from,function(x) which(x==boxes$names))])*-1
  segs$lty<-c(2,2,1)[(segs$bps<=0.05) + (segs$bps<=0.10) + 1]
  segs$bidir<-paste0(segs$from,segs$to) %in% paste0(segs$to,segs$from)
  segs$labpos<-boxes$labpos[sapply(segs$from,function(x) which(boxes$name==x))]
  segs$lwd<-sapply((abs(segs$Estimate)*3)^2+2,function(x) min(x,4))
  labs<-apply(segs[,c(1,3,4)],1,beta)

  segs$Estimate[segs$bps>=0.20]<-0
  segs$lcol<-threecols[sign(segs$Estimate)+2]
  segs$lcol<-col.alpha(segs$lcol,sapply(segs$bps,function(x) max(0,0.45-x+0.05)*2))
  adjust<-0.15
  
  if(new) plot(0,0,xlim=xlim,ylim=ylim,type="n",xlab=NA,ylab=NA,xaxt="n",yaxt="n",bty="n")
  for(i in 1:nrow(segs)){
    xy1<-c(boxes$xs[boxes$name==segs$from[i]],boxes$ys[boxes$name==segs$from[i]])
    xy2<-c(boxes$xs[boxes$name==segs$to[i]],boxes$ys[boxes$name==segs$to[i]])
    #unless horizontal or curved move higher down and lower up
    vertmlt<-(((xy1[2]>xy2[2])-0.5)/0.5)*(xy1[2]!=xy2[2])*(segs$curve[i]==0)
    xy1[2]<-xy1[2]+adjust*vertmlt*-1
    xy2[2]<-xy2[2]+adjust*vertmlt
    #unless vertical move toward center
    horizmlt<-(((xy1[1]>xy2[1])-0.5)/0.5)*(xy1[1]!=xy2[1])
    xy1[1]<-xy1[1]-adjust*horizmlt
    xy2[1]<-xy2[1]-adjust*horizmlt*-1
    #if bidirectional, move up or down
    xy1[2]<-xy1[2]+adjust*0.25*segs$bidir[i]*horizmlt
    xy2[2]<-xy2[2]+adjust*0.25*segs$bidir[i]*horizmlt
    #if curved, move left or right
    xy1[1]<-xy1[1]-sign(segs$curve[i])*adjust
    xy2[1]<-xy2[1]-sign(segs$curve[i])*adjust
    #label up or down if horizontal
    if(xy1[2]==xy2[2] & segs$bidir[i]==TRUE) {
      if(xy1[1]>xy2[1]) thislab<-paste0(labs[i],"\n\n") else thislab<-paste0("\n\n",labs[i])
    } else thislab<-labs[i]
    #label adjustments
    tpos<-ifelse((xy1[1]!=xy2[1]),3,ifelse((xy1[1]<0) & (segs$curve[i]==0) | (segs$curve[i]<0),2,4))
    curvedarrow(xy1,xy2,curve=segs$curve[i],arr.pos=0.5,endhead=FALSE,lcol=segs$lcol[i],lwd=segs$lwd[i],lty=segs$lty[i],lab=thislab,lab.pos=segs$labpos[i],tpos=tpos)
    #text(mean(xy1[1],xy2[2]),mean(xy1[2],xy2[2]),labs[i])
  }
  text(boxes$xs,boxes$ys,boxes$labels)
  if(!is.na(letter)) text(xlim[1],ylim[2],letter,adj=c(0,1),cex=2)
}

boxes1<-data.frame(
  names=c("Inflam","HHInflam","ViInflam","PDSCont","HHPDSCont","ViPDSCont"),
  xs=c(0.4,0.4,0.4,-0.4,-0.4,-0.4),
  ys=c(-0.8,0,0.8,-0.8,0,0.8),
  labpos=c(0.5,0.5,0.5,0.5,0.5,0.5),
  labels=c("Inflam","Household\nInflam","Community\nInflam","Contagion\nDisgust","Household\nContagion\nDisgust","Community\nContagion\nDisgust"),stringsAsFactors = FALSE)

boxes2<-data.frame(
  names=c("Inflam","HHInflam","ViInflam","PDSFood","HHPDSFood","ViPDSFood"),
  xs=c(0.4,0.4,0.4,-0.4,-0.4,-0.4),
  ys=c(-0.8,0,0.8,-0.8,0,0.8),
  labpos=c(0.5,0.5,0.5,0.5,0.5,0.5),
  labels=c("Inflam","Household\nInflam","Community\nInflam","Food\nDisgust","Household\nFood\nDisgust","Community\nFood\nDisgust"),stringsAsFactors = FALSE)

boxes3<-data.frame(
  names=c("Parasites","HHParasites","ViParasites","PDSCont","HHPDSCont","ViPDSCont"),
  xs=c(0.4,0.4,0.4,-0.4,-0.4,-0.4),
  ys=c(-0.8,0,0.8,-0.8,0,0.8),
  labpos=c(0.5,0.5,0.5,0.5,0.5,0.5),
  labels=c("Parasites","Household\nParasites","Community\nParasites","Contagion\nDisgust","Household\nContagion\nDisgust","Community\nContagion\nDisgust"),stringsAsFactors = FALSE)

boxes4<-data.frame(
  names=c("Parasites","HHParasites","ViParasites","PDSFood","HHPDSFood","ViPDSFood"),
  xs=c(0.4,0.4,0.4,-0.4,-0.4,-0.4),
  ys=c(-0.8,0,0.8,-0.8,0,0.8),
  labpos=c(0.5,0.5,0.5,0.5,0.5,0.5),
  labels=c("Parasites","Household\nParasites","Community\nParasites","Food\nDisgust","Household\nFood\nDisgust","Community\nFood\nDisgust"),stringsAsFactors = FALSE)

#Complete model plots
load("mod11.Rdata")
tiff("complete Models.tif",width=3200,height=2500,res=400,pointsize=11,compression="lzw")
par(mfrow=c(2,2),mar=c(0,0,0,0))
pathplot(mod11,boxes1,"A")
pathplot(modF11,boxes2,"B")
pathplot(modP11,boxes3,"C")
pathplot(modPF11,boxes4,"D")
dev.off()

load("mod11b.Rdata")
tiff("complete Modelsb.tif",width=3200,height=2500,res=400,pointsize=11,compression="lzw")
par(mfrow=c(2,2),mar=c(0,0,0,0))
pathplot(mod11b,boxes1,"A")
pathplot(modF11b,boxes2,"B")
pathplot(modP11b,boxes3,"C")
pathplot(modPF11b,boxes4,"D")
dev.off()

#large mediation models
load("modCI.Rdata")
load("modF.Rdata")
load("modP.Rdata")
load("modPF.Rdata")

tx <- -0.8
tiff("Mediation.tif",width=3500,height=3000,res=350,pointsize=11,compression="lzw")
par(mfrow=c(1,4),mar=c(0,0,0,0))
#Inflam-Cont
boxesT<-boxes1
boxesT$ys<-boxesT$ys+2
boxesT$xs<-boxesT$xs-0.15
pathplot(mod1,boxesT,ylim=c(-4.8,2.4))

boxesT$ys<-boxesT$ys+0.2
pathplot(mod3,boxesT,new=FALSE)
boxesT$xs<-boxesT$xs+0.35
pathplot(mod2,boxesT,new=FALSE)
boxesT$xs<-boxesT$xs+0.55

boxesT<-boxes1
pathplot(mod6,boxesT,new=FALSE)
boxesT$ys<-boxesT$ys+0.4
pathplot(mod7,boxesT,new=FALSE)

boxesT<-boxes1
boxesT$ys<-boxesT$ys-1.6
pathplot(mod8,boxesT,new=FALSE)

boxesT<-boxes1
boxesT$ys<-boxesT$ys-4.0
pathplot(mod11,boxesT,new=FALSE)

l<-1
text(tx,2.4,LETTERS[l],cex=2)
text(tx,0.6,LETTERS[l+4],cex=2)
text(tx,-1.4,LETTERS[l+8],cex=2)
text(tx,-3.0,LETTERS[l+12],cex=2)
l<-l+1

#Inflam-Food
boxesT<-boxes2
boxesT$ys<-boxesT$ys+2
boxesT$xs<-boxesT$xs-0.15
pathplot(modF1,boxesT,ylim=c(-4.8,2.4))

boxesT$ys<-boxesT$ys+0.2
pathplot(modF3,boxesT,new=FALSE)
boxesT$xs<-boxesT$xs+0.35
pathplot(modF2,boxesT,new=FALSE)
boxesT$xs<-boxesT$xs+0.55

boxesT<-boxes2
pathplot(modF6,boxesT,new=FALSE)
boxesT$ys<-boxesT$ys+0.4
pathplot(modF7,boxesT,new=FALSE)

boxesT<-boxes2
boxesT$ys<-boxesT$ys-1.6
pathplot(modF8,boxesT,new=FALSE)

boxesT<-boxes2
boxesT$ys<-boxesT$ys-4.0
pathplot(modF11,boxesT,new=FALSE)

text(tx,2.4,LETTERS[l],cex=2)
text(tx,0.6,LETTERS[l+4],cex=2)
text(tx,-1.4,LETTERS[l+8],cex=2)
text(tx,-3.0,LETTERS[l+12],cex=2)
l<-l+1

#Parasites-Cont
boxesT<-boxes3
boxesT$ys<-boxesT$ys+2
boxesT$xs<-boxesT$xs-0.15
pathplot(modP1,boxesT,ylim=c(-4.8,2.4))

boxesT$ys<-boxesT$ys+0.2
pathplot(modP3,boxesT,new=FALSE)
boxesT$xs<-boxesT$xs+0.35
pathplot(modP2,boxesT,new=FALSE)
boxesT$xs<-boxesT$xs+0.55

boxesT<-boxes3
pathplot(modP6,boxesT,new=FALSE)
boxesT$ys<-boxesT$ys+0.4
pathplot(modP7,boxesT,new=FALSE)

boxesT<-boxes3
boxesT$ys<-boxesT$ys-1.6
pathplot(modP8,boxesT,new=FALSE)

boxesT<-boxes3
boxesT$ys<-boxesT$ys-4.0
pathplot(modP11,boxesT,new=FALSE)

text(tx,2.4,LETTERS[l],cex=2)
text(tx,0.6,LETTERS[l+4],cex=2)
text(tx,-1.4,LETTERS[l+8],cex=2)
text(tx,-3.0,LETTERS[l+12],cex=2)
l<-l+1

#Parasites-Food
boxesT<-boxes4
boxesT$ys<-boxesT$ys+2
boxesT$xs<-boxesT$xs-0.15
pathplot(modPF1,boxesT,ylim=c(-4.8,2.4))

boxesT$ys<-boxesT$ys+0.2
pathplot(modPF3,boxesT,new=FALSE)
boxesT$xs<-boxesT$xs+0.35
pathplot(modPF2,boxesT,new=FALSE)
boxesT$xs<-boxesT$xs+0.55

boxesT<-boxes4
pathplot(modPF6,boxesT,new=FALSE)
boxesT$ys<-boxesT$ys+0.4
pathplot(modPF7,boxesT,new=FALSE)

boxesT<-boxes4
boxesT$ys<-boxesT$ys-1.6
pathplot(modPF8,boxesT,new=FALSE)

boxesT<-boxes4
boxesT$ys<-boxesT$ys-4.0
pathplot(modPF11,boxesT,new=FALSE)

text(tx,2.4,LETTERS[l],cex=2)
text(tx,0.6,LETTERS[l+4],cex=2)
text(tx,-1.4,LETTERS[l+8],cex=2)
text(tx,-3.0,LETTERS[l+12],cex=2)
l<-l+1
dev.off()


tx <- -0.8
tiff("Mediationb.tif",width=3500,height=3000,res=300,pointsize=11,compression="lzw")
par(mfrow=c(1,4),mar=c(0,0,0,0))
#Inflam-Cont
boxesT<-boxes1
boxesT$ys<-boxesT$ys+2
boxesT$xs<-boxesT$xs-0.15
pathplot(mod1b,boxesT,ylim=c(-4.8,2.4))
boxesT$ys<-boxesT$ys+0.2
pathplot(mod2b,boxesT,new=FALSE)
boxesT$xs<-boxesT$xs+0.3
pathplot(mod3b,boxesT,new=FALSE)

boxesT<-boxes1
pathplot(mod6b,boxesT,new=FALSE)
boxesT$ys<-boxesT$ys+0.4
pathplot(mod7b,boxesT,new=FALSE)

boxesT<-boxes1
boxesT$ys<-boxesT$ys-1.6
pathplot(mod8b,boxesT,new=FALSE)

boxesT<-boxes1
boxesT$ys<-boxesT$ys-4.0
pathplot(mod11b,boxesT,new=FALSE)

l<-1
text(tx,2.4,LETTERS[l],cex=3)
text(tx,0.6,LETTERS[l+4],cex=3)
text(tx,-1.4,LETTERS[l+8],cex=3)
text(tx,-3.0,LETTERS[l+12],cex=3)
l<-l+1

#Inflam-Food
boxesT<-boxes2
boxesT$ys<-boxesT$ys+2
boxesT$xs<-boxesT$xs-0.15
pathplot(modF1b,boxesT,ylim=c(-4.8,2.4))

boxesT$ys<-boxesT$ys+0.2
pathplot(modF2b,boxesT,new=FALSE)
boxesT$xs<-boxesT$xs+0.3
pathplot(modF3b,boxesT,new=FALSE)

boxesT$xs<-boxesT$xs+0.6

boxesT<-boxes2
pathplot(modF6b,boxesT,new=FALSE)
boxesT$ys<-boxesT$ys+0.4
pathplot(modF7b,boxesT,new=FALSE)

boxesT<-boxes2
boxesT$ys<-boxesT$ys-1.6
pathplot(modF8b,boxesT,new=FALSE)

boxesT<-boxes2
boxesT$ys<-boxesT$ys-4.0
pathplot(modF11b,boxesT,new=FALSE)

text(tx,2.4,LETTERS[l],cex=3)
text(tx,0.6,LETTERS[l+4],cex=3)
text(tx,-1.4,LETTERS[l+8],cex=3)
text(tx,-3.0,LETTERS[l+12],cex=3)
l<-l+1

#Parasites-Cont
boxesT<-boxes3
boxesT$ys<-boxesT$ys+2
boxesT$xs<-boxesT$xs-0.15
pathplot(modP1b,boxesT,ylim=c(-4.8,2.4))
boxesT$ys<-boxesT$ys+0.2
pathplot(modP2b,boxesT,new=FALSE)
boxesT$xs<-boxesT$xs+0.3
pathplot(modP3b,boxesT,new=FALSE)

boxesT<-boxes3
pathplot(modP6b,boxesT,new=FALSE)
boxesT$ys<-boxesT$ys+0.4
pathplot(modP7b,boxesT,new=FALSE)

boxesT<-boxes3
boxesT$ys<-boxesT$ys-1.6
pathplot(modP8b,boxesT,new=FALSE)

boxesT<-boxes3
boxesT$ys<-boxesT$ys-4.0
pathplot(modP11b,boxesT,new=FALSE)

text(tx,2.4,LETTERS[l],cex=3)
text(tx,0.6,LETTERS[l+4],cex=3)
text(tx,-1.4,LETTERS[l+8],cex=3)
text(tx,-3.0,LETTERS[l+12],cex=3)
l<-l+1

#Parasites-Food
boxesT<-boxes4
boxesT$ys<-boxesT$ys+2
boxesT$xs<-boxesT$xs+0.15
pathplot(modPF1b,boxesT,ylim=c(-4.8,2.4))

boxesT$ys<-boxesT$ys+0.2
pathplot(modPF3b,boxesT,new=FALSE)
boxesT$xs<-boxesT$xs-0.3
pathplot(modPF2b,boxesT,new=FALSE)
boxesT$xs<-boxesT$xs+0.6

boxesT<-boxes4
pathplot(modPF6b,boxesT,new=FALSE)
boxesT$ys<-boxesT$ys+0.4
pathplot(modPF7b,boxesT,new=FALSE)

boxesT<-boxes4
boxesT$ys<-boxesT$ys-1.6
pathplot(modPF8b,boxesT,new=FALSE)

boxesT<-boxes4
boxesT$ys<-boxesT$ys-4.0
pathplot(modPF11b,boxesT,new=FALSE)

text(tx,2.4,LETTERS[l],cex=3)
text(tx,0.6,LETTERS[l+4],cex=3)
text(tx,-1.4,LETTERS[l+8],cex=3)
text(tx,-3.0,LETTERS[l+12],cex=3)
l<-l+1
dev.off()

#Market Integration Path Plots
boxes1MI<-data.frame(
  names=c("Inflam","HHInflam","ViInflam","PDSCont","HHPDSCont","ViPDSCont","MSOL","HOUSE","TSOL","MSOL2","HOUSE2","TSOL2"),
  xs=c(0.5,0.5,0.5,-0.5,-0.5,-0.5,1.5,1.5,1.5,-1.5,-1.5,-1.5),
  ys=c(-0.8,0,0.8,-0.8,0,0.8,0.4,-0.4,-1.2,0.4,-0.4,-1.2),
  labpos=c(0.5,0.6,0.4,0.5,0.6,0.4,0.3,0.3,0.3,0.3,0.3,0.3),
  labels=c("Inflam","Household\nInflam","Community\nInflam","Contagion\nDisgust","Household\nContagion\nDisgust","Community\nContagion\nDisgust","MSOL","HSOL","TSOL","MSOL","HSOL","TSOL"),stringsAsFactors = FALSE)

boxes2MI<-data.frame(
  names=c("Inflam","HHInflam","ViInflam","PDSFood","HHPDSFood","ViPDSFood","MSOL","HOUSE","TSOL","MSOL2","HOUSE2","TSOL2"),
  xs=c(0.5,0.5,0.5,-0.5,-0.5,-0.5,1.5,1.5,1.5,-1.5,-1.5,-1.5),
  ys=c(-0.8,0,0.8,-0.8,0,0.8,0.4,-0.4,-1.2,0.4,-0.4,-1.2),
  labpos=c(0.5,0.6,0.4,0.5,0.6,0.4,0.3,0.3,0.3,0.3,0.3,0.3),
  labels=c("Inflam","Household\nInflam","Community\nInflam","Food\nDisgust","Household\nFood\nDisgust","Community\nFood\nDisgust","MSOL","HSOL","TSOL","MSOL","HSOL","TSOL"),stringsAsFactors = FALSE)

boxes3MI<-data.frame(
  names=c("Parasites","HHParasites","ViParasites","PDSCont","HHPDSCont","ViPDSCont","MSOL","HOUSE","TSOL","MSOL2","HOUSE2","TSOL2"),
  xs=c(0.5,0.5,0.5,-0.5,-0.5,-0.5,1.5,1.5,1.5,-1.5,-1.5,-1.5),
  ys=c(-0.8,0,0.8,-0.8,0,0.8,0.4,-0.4,-1.2,0.4,-0.4,-1.2),
  labpos=c(0.5,0.6,0.4,0.5,0.6,0.4,0.3,0.3,0.3,0.3,0.3,0.3),
  labels=c("Parasites","Household\nParasites","Community\nParasites","Contagion\nDisgust","Household\nContagion\nDisgust","Community\nContagion\nDisgust","MSOL","HSOL","TSOL","MSOL","HSOL","TSOL"),stringsAsFactors = FALSE)

boxes4MI<-data.frame(
  names=c("Parasites","HHParasites","ViParasites","PDSFood","HHPDSFood","ViPDSFood","MSOL","HOUSE","TSOL","MSOL2","HOUSE2","TSOL2"),
  xs=c(0.5,0.5,0.5,-0.5,-0.5,-0.5,1.5,1.5,1.5,-1.5,-1.5,-1.5),
  ys=c(-0.8,0,0.8,-0.8,0,0.8,0.4,-0.4,-1.2,0.4,-0.4,-1.2),
  labpos=c(0.5,0.6,0.4,0.5,0.6,0.4,0.3,0.3,0.3,0.3,0.3,0.3),
  labels=c("Parasites","Household\nParasites","Community\nParasites","Food\nDisgust","Household\nFood\nDisgust","Community\nFood\nDisgust","MSOL","HSOL","TSOL","MSOL","HSOL","TSOL"),stringsAsFactors = FALSE)


load("mod11MI.Rdata")
tiff("MI Models.tif",width=3700,height=3000,res=350,pointsize=11,compression="lzw")
par(mfrow=c(2,2),mar=c(0,0,0,0))
tx <- -1.5
ty <- 0.8
#Inflam-Cont
pathplot(mod11MI,boxes1MI,xlim=c(-1.8,1.8),ylim=c(-1.5,1),split=TRUE)
l<-1
text(tx,ty,LETTERS[l],cex=3); l <- l +1
pathplot(modF11MI,boxes2MI,xlim=c(-1.8,1.8),ylim=c(-1.5,1),split=TRUE)
text(tx,ty,LETTERS[l],cex=3); l <- l +1
pathplot(modP11MI,boxes3MI,xlim=c(-1.8,1.8),ylim=c(-1.5,1),split=TRUE)
text(tx,ty,LETTERS[l],cex=3); l <- l +1
pathplot(modPF11MI,boxes4MI,xlim=c(-1.8,1.8),ylim=c(-1.5,1),split=TRUE)
text(tx,ty,LETTERS[l],cex=3); l <- l +1
dev.off()
