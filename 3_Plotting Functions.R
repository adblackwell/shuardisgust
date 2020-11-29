##For path plots
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

## Note: curvedarrow source code is borrow from the 'diagram' package
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
    segments$from[(segments$from %in% c("MSOL","TSOL","HOUSE")) & (segments$to %in% c("PDSTotal","HHPDSTotal","PDSCont","HHPDSCont","PDSFood","HHPDSFood"))]<-paste0(segments$from[(segments$from %in% c("MSOL","TSOL","HOUSE")) & (segments$to %in% c("PDSTotal","HHPDSTotal","PDSCont","HHPDSCont","PDSFood","HHPDSFood"))],"2")
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

##For MI density plots. When there are two models with the same pathway, these two posteriors are combined.
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

denspoly<-function(den,hdi,col,shift){
  y<-den$y[den$x >= hdi[1] & den$x <= hdi[2]]
  y<-y/max(y) + shift
  x<-den$x[den$x >= hdi[1] & den$x <= hdi[2]]
  lines(x,y)
  polygon(c(x[1],x,x[length(x)]),c(shift,y,shift),col=col.alpha(col,0.7),border=NA)
}

#Similarity Matrix Plot
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

