library(rethinking)
#For 'rethinking' see https://xcelab.net/rm/software/ 
#and Statistical Rethinking: A Bayesian Course with Examples in R and STAN, 2nd Edition, by Richard McElreath
#Used here for convenience and for fast Bayesian estimates. Could easily replace that code with something else.

#Function to simulate a biomarker sequence based on given parameters
simbiomarker<-function(DP=0.80, #disgust protection (per sd of disgust)
                       freq=0.015, #Infection frequency, average risk per day
                       duration=7, #infection duration in days
                       durationsd=1.5, #infectionsd
                       durationlog=TRUE, #whether infection duration is log normal
                       hl=0.95, #biomarker reduction per day
                       stack=FALSE, #infections stack (i.e. for helminths where load can increase)
                       days=730,#length to simulate)
                       n=75){ #sample size   
  disgust<-rnorm(n,0,1)
  dailyrisk <- 1-((1-freq)^exp(log(DP)*disgust))
    exp(log(DP)*disgust) * freq
  newinfections<-t(sapply(dailyrisk,function(risk) sample(c(0,1),days,replace=TRUE,prob=c(1,risk))))
  
  infections<-newinfections
  for(i in 1:nrow(infections)){
    for(j in 1:(ncol(infections)-1)){
      if(newinfections[i,j]==1){
        if(durationlog) dur<-round(exp(rnorm(1,log(duration),log(durationsd)))) else dur<-round(rnorm(1,duration,durationsd))
        dur<-min(ncol(newinfections)-j,dur)
        infections[i,(j+1):(j+dur-1)]<-infections[i,(j+1):(j+dur-1)]+1
    }
  }}

  if(!stack) infections[infections>0]<-1
  
  biomarker<-infections
  for(i in 2:ncol(biomarker)){
    previous<-biomarker[,i-1]
    now<-biomarker[,i-1]*hl
    biomarker[,i][biomarker[,i]<previous]<-round(now[biomarker[,i]<previous],3)
  }
  return(list(disgust=disgust,newinfections=newinfections,infections=infections,biomarker=biomarker))  
}

#Function to sample biomarkers from a simulated sequence and check association with disgust variable
#Function also tests whether continuous variables provide more power than converting variables to binary infection variables
biosampmap<-function(simbiom,n=75,cutoff=quantile(simbiom$biomarker,0.75)){
  samp<-sample(nrow(simbiom$biomarker),n)
  biosamp<-apply(simbiom$biomarker[samp,(ncol(simbiom$biomarker)-365):ncol(simbiom$biomarker)],1,function(x) x[sample(length(x),1)])
  infectsamp<-ifelse(biosamp>=cutoff,1,0)
  d<-data.frame(biosamp=scale(log(biosamp+0.1)),infectsamp=infectsamp,disgust=scale(simbiom$disgust[samp]))

  m1 <- quap( alist(
    biosamp ~ dnorm( mu , sigma ) ,
    mu <- a+b*disgust ,
    c(a,b) ~ dnorm(0,1) , 
    sigma ~ dcauchy(0,1)
  ) , data=d,start = list(a=0,b=0,sigma=1))
  
  m2 <- quap( alist(
    infectsamp ~ dbinom( 1 , p ),
    logit(p) <-  a+b*disgust,
    c(a,b) ~ dnorm(0,1)), data=d)
  
  post1<-extract.samples(m1,1000)[,2]
  post2<-extract.samples(m2,1000)[,2]
  return(cbind(post1,post2))
}

#single simulations for illustration
inflammation<-simbiomarker(DP=0.80,freq=0.01,duration=7,durationsd=1.5,durationlog=TRUE,hl=0.95,stack=FALSE,days=1000,n=750)
  
helminth<-simbiomarker(DP=0.80,freq=0.005,duration=2000,durationsd=300,durationlog=FALSE,hl=0.99,stack=TRUE,days=7000,n=750)
  
tiff("infection.tif",width=2800,height=3000,res=400,pointsize=12,compression="lzw")
par(mfrow=c(3,2),mar=c(4,5,2,1))
plot(inflammation$biomarker[1,-(1:365)],type="l",ylim=c(0,1),xlab=NA,ylab="Biomarker",main="Inflammation")
plot(helminth$biomarker[1,-(1:3500)],type="l",ylim=c(0,20),xlab=NA,ylab="Biomarker",main="Parasite")
plot(inflammation$biomarker[2,-(1:365)],type="l",ylim=c(0,1),xlab=NA,ylab="Biomarker")
plot(helminth$biomarker[2,-(1:3500)],type="l",ylim=c(0,20),xlab=NA,ylab="Biomarker")
plot(inflammation$biomarker[3,-(1:365)],type="l",ylim=c(0,1),xlab="Days",ylab="Biomarker")
plot(helminth$biomarker[3,-(1:3500)],type="l",ylim=c(0,20),xlab="Days",ylab="Biomarker")
dev.off()

##False Positive test, with no effect 
inflammation<-simbiomarker(DP=1,freq=0.025,duration=7,durationsd=1.5,durationlog=TRUE,hl=0.95,stack=FALSE,days=1000,n=750)
inflamtest<-foreach(i=1:300, .combine=rbind, .packages = "rethinking") %dopar% {
  biosampmap(inflammation,75)
}
PerBPs<-t(sapply(PerBP,function(y) apply(inflamtest,2,function(x) c(mean(x[(y-l+1):y]), sum(x[(y-l+1):y]<0)/l))))
#Trials with effects < -0.25
sum(PerBPs[,1] < -0.25)/nrow(PerBPs)
#Trials with effects < -0.25 and posteriors that 95% exclude 0
sum((PerBPs[,1] < -0.25) & (PerBPs[,2] > 0.95))/nrow(PerBPs)


library(foreach)
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

#power test
DPs<-seq(0.50,1.0,0.05)
freqs<-seq(0.001,0.04,0.002)
tests<-expand.grid(DP=DPs, freq=freqs)
#sample size of 75
inflampowermap<-foreach(DP=tests$DP,freq=tests$freq, .combine=rbind, .packages = "foreach") %dopar% {
  inflammation<-simbiomarker(DP=DP,freq=freq,duration=7,durationsd=1.5,durationlog=TRUE,hl=0.95,stack=FALSE,days=1000,n=750)
  inflamtest<-foreach(i=1:50, .combine=rbind, .packages = "rethinking") %dopar% {
    biosampmap(inflammation,75)
  }
  means<-apply(inflamtest,2,mean)
  hpdis<-apply(inflamtest,2,HPDI,prob=0.95)
  Bps<-apply(inflamtest,2,function(x) sum(x<0)/length(x))
  l<-1000
  PerBP<-seq(l,nrow(inflamtest),l)
  PerBPs<-t(sapply(PerBP,function(y) apply(inflamtest,2,function(x)(sum(x[(y-l+1):y]<0)/l) > 0.90)))
  c(means[1],hpdis[,1],Bps[1],means[2],hpdis[,2],Bps[2],colSums(PerBPs)/nrow(PerBPs))
}

#helminths
DPs<-seq(0.50,1.0,0.05)
freqs<-seq(0.001,0.04,0.002)
tests<-expand.grid(DP=DPs, freq=freqs)

helminthpowermap<-foreach(DP=tests$DP,freq=tests$freq, .combine=rbind, .packages = "foreach") %dopar% {
  helminth<-simbiomarker(DP=DP,freq=freq,duration=2000,durationsd=300,durationlog=FALSE,hl=0.99,stack=TRUE,days=7000,n=750)
  inflamtest<-foreach(i=1:50, .combine=rbind, .packages = "rethinking") %dopar% {
    biosampmap(helminth,75)
  }
  means<-apply(inflamtest,2,mean)
  hpdis<-apply(inflamtest,2,HPDI,prob=0.95)
  Bps<-apply(inflamtest,2,function(x) sum(x<0)/length(x))
  l<-1000
  PerBP<-seq(l,nrow(inflamtest),l)
  PerBPs<-t(sapply(PerBP,function(y) apply(inflamtest,2,function(x)(sum(x[(y-l+1):y]<0)/l) > 0.90)))
  c(means[1],hpdis[,1],Bps[1],means[2],hpdis[,2],Bps[2],colSums(PerBPs)/nrow(PerBPs))
}

save(inflampowermap,helminthpowermap,file="simulations.Rdata")


#False positive test



tiff("parameter.tif",width=2000,height=1000,res=300,pointsize=9,compression="lzw")

par(mfrow=c(1,2))
freq<-0.015
plot(inflampowermap[tests$freq==freq,1],tests$DP[tests$freq==freq],xlab="Biomarker Parameter Estimate",ylab="Odds-ratio for infection for 1 SD Disgust",xlim=c(-1.1,0.5),pch=19)
for(i in 1:length(tests$DP[tests$freq==freq])){
  lines(inflampowermap[tests$freq==freq,][i,2:3],tests$DP[tests$freq==freq][c(i,i)])
}
abline(v=0)
p<-par("usr")
mtext("A", side=3, line=1, at=p[1]-(p[2]-p[1])/10,cex=1.5)

plot(helminthpowermap[tests$freq==freq,1],tests$DP[tests$freq==freq],xlab="Biomarker Parameter Estimate",ylab="Odds-ratio for infection for 1 SD Disgust",xlim=c(-1.1,0.5),pch=19)
for(i in 1:length(tests$DP[tests$freq==freq])){
  lines(helminthpowermap[tests$freq==freq,][i,2:3],tests$DP[tests$freq==freq][c(i,i)])
}
abline(v=0)
p<-par("usr")
mtext("B", side=3, line=1, at=p[1]-(p[2]-p[1])/10,cex=1.5)

dev.off()

#length and frequency
durations<-exp(seq(1,8,0.2))
freqs<-seq(0.001,0.04,0.002)
tests<-expand.grid(duration=durations, freq=freqs)
lfpower<-foreach(duration=tests$duration,freq=tests$freq, .combine=rbind, .packages = "foreach") %dopar% {
  inflammation<-simbiomarker(DP=0.7,freq=freq,duration=duration,durationsd=(duration/6),durationlog=FALSE,hl=0.95,stack=TRUE,days=7000,n=750)
  inflamtest<-foreach(i=1:50, .combine=rbind) %dopar% {
    biosampmap(inflammation,75)
  }
  means<-apply(inflamtest,2,mean)
  hpdis<-apply(inflamtest,2,HPDI,prob=0.95)
  Bps<-apply(inflamtest,2,function(x) sum(x<0)/length(x))
  l<-1000
  PerBP<-seq(l,nrow(inflamtest),l)
  PerBPs<-t(sapply(PerBP,function(y) apply(inflamtest,2,function(x)(sum(x[(y-l+1):y]<0)/l) > 0.90)))
  c(means[1],hpdis[,1],Bps[1],means[2],hpdis[,2],Bps[2],colSums(PerBPs)/nrow(PerBPs))
}

#modified version of filled contour to move layoutr commands out of function (so plots can be combined)
filled.contour2<-function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
                                                                        length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
                           ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
                           levels = pretty(zlim, nlevels), nlevels = 20, color.palette = function(n) hcl.colors(n,"YlOrRd", rev = TRUE), col = color.palette(length(levels) - 1), plot.title, plot.axes, key.title, key.axes, asp = NA, 
                           xaxs = "i", yaxs = "i", las = 1, axes = TRUE, cont=NA, 
                           frame.plot = axes, ...) 
{
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  # mar.orig <- (par.orig <- par(c("mar", "las", 
  #                                "mfrow")))$mar
  # on.exit(par(par.orig))
  # w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  # layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  par(las = las)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
              yaxs = "i")
  rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
  if (missing(key.axes)) {
    if (axes) 
      axis(4)
  }
  else key.axes
  box()
  if (!missing(key.title)) 
    key.title
  mar <- mar.orig
  mar[4L] <- 1
  par(mar = mar)
  plot.new()
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, 
              asp = asp)
  .filled.contour(x, y, z, levels, col)
  if(!is.na(cont)) contour(x,y,z,levels=cont,drawlabels = FALSE, axes = FALSE, frame.plot = FALSE, add = TRUE)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}

tiff("powerplot.tif",width=3500,height=2600,res=400,pointsize=12,compression="lzw")

par(mar=c(4,4,3,1))
mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
w <- (3 + mar.orig[2L]) * par("csi") * 1.5
layout(matrix(c(2, 1, 4, 3, 6, 5, 8, 7), ncol = 4,byrow=TRUE), widths = c(1, lcm(w),1, lcm(w)))

z<-matrix(inflampowermap[,9],ncol=length(DPs),byrow=TRUE)
filled.contour2(1-(1-freqs)^30,DPs,z,xlab="Monthly Probability of Infection",ylab="Odds-ratio for infection for 1 SD Disgust",key.title=mtext("Power",3,1,cex=0.7),zlim=c(0,1),cont=0.8)
p<-par("usr")
mtext("A", side=3, line=1, at=p[1]-(p[2]-p[1])/10,cex=1.5)

z<-matrix(inflampowermap[,3],ncol=length(DPs),byrow=TRUE) - matrix(inflampowermap[,2],ncol=length(DPs),byrow=TRUE)
filled.contour2(1-(1-freqs)^30,DPs,z,xlab="Monthly Probability of Infection",ylab="Odds-ratio for infection for 1 SD Disgust",key.title=mtext("HPDI",3,1,cex=0.7))
p<-par("usr")
mtext("B", side=3, line=1, at=p[1]-(p[2]-p[1])/10,cex=1.5)

z<-matrix(helminthpowermap[,9],ncol=length(DPs),byrow=TRUE)
filled.contour2(1-(1-freqs)^30,DPs,z,xlab="Monthly Probability of Infection",ylab="Odds-ratio for infection for 1 SD Disgust",key.title=mtext("Power",3,1,cex=0.7),zlim=c(0,1),cont=0.8)
p<-par("usr")
mtext("C", side=3, line=1, at=p[1]-(p[2]-p[1])/10,cex=1.5)

z<-abs(matrix(helminthpowermap[,3],ncol=length(DPs),byrow=TRUE) - matrix(inflampowermap[,2],ncol=length(DPs),byrow=TRUE))
filled.contour2(1-(1-freqs)^30,DPs,z,xlab="Monthly Probability of Infection",ylab="Odds-ratio for infection for 1 SD Disgust",key.title=mtext("HPDI",3,1,cex=0.7))
p<-par("usr")
mtext("D", side=3, line=1, at=p[1]-(p[2]-p[1])/10,cex=1.5)

dev.off()

tiff("durationplot.tif",width=1750,height=1300,res=300,pointsize=12,compression="lzw")

mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
w <- (3 + mar.orig[2]) * par("csi") * 2.54
layout(matrix(c(2, 1), nc = 2), widths = c(1, lcm(w)))

z<-matrix(lfpower[,9],ncol=length(durations),byrow=TRUE)
labs<-c(5,10,20,40,80,160,320,640,1280,2560)
xlabs <- 1-(1-seq(0,0.04,1))^30
filled.contour2(1-(1-freqs)^30,log(durations),z,xlab="Monthly Probability of Infection",ylab="Infection Duration (Days)",key.title="Power",zlim=c(0.5,1),plot.axes = {axis(2,log(labs),labs); axis(1,seq(0,0.7,0.1))},cont=0.8)

dev.off()