library(flextable)
library(officer)

##Function for a less compact version of the tables. Ultimately not used in the paper but left here for posterity.
# modelflex<-function(model,cor=TRUE){
#   mt<-as.data.frame(summary(model)$fixed[,1:4])
#   mr<-as.data.frame(summary(model)$random$Family[,1:4])
#   row.names(mr)<-paste0("h",row.names(mr))
#   if(cor==FALSE) mr<-mr[!grepl("[hv]cor\\(",row.names(mr)),]
#   mr2<-as.data.frame(summary(model)$random$Village[,1:4])
#   if(nrow(mr2)>0) row.names(mr2)<-paste0("v",row.names(mr2))
#   if(cor==FALSE | nrow(mr2)>0) mr2<-mr2[!grepl("[hv]cor\\(",row.names(mr2)),]
#   mt<-rbind(mt,mr)
#   mt<-rbind(mt,mr2)
#   vnames<-row.names(mt)
#   vnames<-gsub("PDSCont","Contagion Disgust",vnames)
#   vnames<-gsub("PDSFood","Food Disgust",vnames)
#   vnames<-gsub("PDSPest","Component 3",vnames)
#   vnames<-gsub("PDSTotal","Total Disgust",vnames)
#   vars<-Reduce(rbind,strsplit(vnames,"_"))
#   mt$Dependent<-vars[,1]
#   mt$Independent<-vars[,2]
#   vars2<-Reduce(rbind,strsplit(vnames,"[(_]"))
#   mt$Dependent[grepl("[hv]sd\\(",vnames)]<-vars2[grepl("[hv]sd\\(",vnames),2]
#   mt$Independent[grepl("hsd\\(",vnames)]<-"sd(Household)"
#   mt$Independent[grepl("vsd\\(",vnames)]<-"sd(Community)"
#   mt$Dependent[grepl("[hv]cor\\(",vnames)]<-NA
#   mt$Independent[grepl("hcor\\(",vnames)]<-"cor(Household)"
#   mt$Independent[grepl("vcor\\(",vnames)]<-"cor(Community)"
#   
#   mf<-capture.output(summary(model)$formula)
#   mf<-c("Model Formula:",mf)
#   mf<-gsub("PDSCont","Contagion",mf)
#   mf<-gsub("PDSFood","Food",mf)
#   mf<-gsub("PDSPest","Component3",mf)
#   mf<-gsub("Family","Household",mf)
#   mf<-c(mf,"HH = Household mean, excluding target individual. Vi = Village mean, excluding target household","Items below the grey bar are group level effects for Household")
#   
#   ft<-flextable(mt,col_keys = c("Dependent","Independent","Estimate","l-95% CI","u-95% CI"))
#   ft<-colformat_num(ft,j=c("Estimate","l-95% CI","u-95% CI"))
#   ft<-autofit(ft)
#   ft<-hline(ft,i=nrow(mt)-(nrow(mr)+nrow(mr2)), border = fp_border(color="gray"))
#   ft<-add_footer_lines(ft,values=mf)
#   ft<-font(ft,fontname="Calibri",part="all")
#   ft<-fontsize(ft,size=11)
#   ft<-fontsize(ft, part="footer", size=10)
#   ft<-padding(ft, part="footer", padding=0)
#   ft<-bold(ft,part="header")
#   ft
# }
# 
# 
# load("modT11.Rdata")
# ft11<-modelflex(modT11)
# ftp11<-modelflex(modTP11)
# save_as_docx("S1" = ft11, "S2" = ftp11, path = "Model Tables Single Component.docx")
# 
# load("mod11.Rdata")
# f11<-modelflex(mod11)
# ff11<-modelflex(modF11)
# fp11<-modelflex(modP11)
# fpf11<-modelflex(modPF11)
# fv11<-modelflex(modV11)
# fpv11<-modelflex(modPV11)
# 
# save_as_docx("S8" = f11, "S9" = ff11, "S10" = fv11, "S11" = fp11, "S12" = fpf11, "S13" = fpv11, path = "Model Tables.docx")

modelflexC<-function(models,cor=TRUE,columns=rep(NA,length(models))){
  mt2<-NA
  for(model in models){
    mt<-as.data.frame(summary(model)$fixed[,1:4])
    mr<-as.data.frame(summary(model)$random$Family[,1:4])
    row.names(mr)<-paste0("h",row.names(mr))
    if(cor==FALSE) mr<-mr[!grepl("[hv]cor\\(",row.names(mr)),]
    mr2<-as.data.frame(summary(model)$random$Village[,1:4])
    if(nrow(mr2)>0) row.names(mr2)<-paste0("v",row.names(mr2))
    if(cor==FALSE & nrow(mr2)>0) mr2<-mr2[!grepl("[hv]cor\\(",row.names(mr2)),]
    mt<-rbind(mt,mr)
    mt<-rbind(mt,mr2)
    mt$Estimate<-paste0(format(round(mt$Estimate,2),nsmall=2,digits=2)," (",format(round(mt$"l-95% CI",2),nsmall=2,digits=2),",",format(round(mt$"u-95% CI",2),nsmall=2,digits=2) ,")")
    vnames<-row.names(mt)
    vnames<-gsub("PDSCont","Disgust",vnames)
    vnames<-gsub("PDSFood","Disgust",vnames)
    vnames<-gsub("PDSPest","Disgust",vnames)
    vnames<-gsub("PDSTotal","Disgust",vnames)
    vars<-Reduce(rbind,strsplit(vnames,"_"))
    mt$Dependent<-vars[,1]
    mt$Independent<-vars[,2]
    vars2<-Reduce(rbind,strsplit(vnames,"[(_]"))
    mt$Dependent[grepl("[hv]sd\\(",vnames)]<-vars2[grepl("[hv]sd\\(",vnames),2]
    mt$Independent[grepl("hsd\\(",vnames)]<-"sd(Household)"
    mt$Independent[grepl("vsd\\(",vnames)]<-"sd(Community)"
    mt$Dependent[grepl("[hv]cor\\(",vnames)]<-NA
    mt$Independent[grepl("hcor\\(",vnames)]<-"cor(Household)"
    mt$Independent[grepl("vcor\\(",vnames)]<-"cor(Community)"
    
    mf<-capture.output(summary(model)$formula)
    mf<-c("Model Formula:",mf)
    mf<-gsub("PDSCont","Disgust",mf)
    mf<-gsub("PDSFood","Disgust",mf)
    mf<-gsub("PDSPest","Disgust",mf)
    mf<-gsub("Family","Household",mf)
    mf<-c(mf,"HH = Household mean, excluding target individual. Vi = Village mean, excluding target household","Items below the grey bar are group level effects for Household")
    if(is.na(mt2)) mt2<-mt[,c("Dependent","Independent","Estimate")] 
    else {
      if(mt2$Dependent[1]!=mt$Dependent[1]){
        mt2$Dependent<-gsub("Inflam","Infection",mt2$Dependent)
        mt2$Dependent<-gsub("Parasites","Infection",mt2$Dependent)
        mt2$Independent<-gsub("Inflam","Infection",mt2$Independent)
        mt2$Independent<-gsub("Parasites","Infection",mt2$Independent)
      }
      mt2<-cbind(mt2,mt$Estimate)
    }
  }
  names(mt2)<-c("Dependent","Independent",columns)
  ft<-flextable(mt2)
  #ft<-colformat_num(ft,j=c("Estimate","l-95% CI","u-95% CI"))
  ft<-autofit(ft)
  ft<-hline(ft,i=nrow(mt)-(nrow(mr)+nrow(mr2)), border = fp_border(color="gray"))
  ft<-add_footer_lines(ft,values=mf)
  ft<-font(ft,fontname="Calibri",part="all")
  ft<-fontsize(ft,size=11)
  ft<-fontsize(ft, part="footer", size=10)
  ft<-padding(ft, part="footer", padding=0)
  ft<-bold(ft,part="header")
  ft
}

fcT11<-modelflexC(list(modT11,modTP11),TRUE,c("Inflammation","Parasites"))
fc3I<-modelflexC(list(mod11,modF11,modV11),TRUE,c("C1:Contagion","C2:Food","C3:Other"))
fc3P<-modelflexC(list(modP11,modPF11,modPV11),TRUE,c("C1:Contagion","C2:Food","C3:Other"))

save_as_docx("S6" = fcT11, "S7" = fc3I, "S8" = fc3P, path = "Compact Model Tables.docx")


fcMI1<-modelflexC(list(modT11MI,modTP11MI),FALSE,c("Inflammation","Parasites"))
fcMI3I<-modelflexC(list(mod11MI,modF11MI,modV11MI),FALSE,c("C1:Contagion","C2:Food","C3:Other"))
fcMI3P<-modelflexC(list(modP11MI,modPF11MI,modPV11MI),FALSE,c("C1:Contagion","C2:Food","C3:Other"))

save_as_docx("S9" = fcMI1, "S10" = fcMI3I, "S11" = fcMI3P, path = "Compact Model MI Tables.docx")


modelflex2<-function(models){
  mts<-NA
  mfs<-NA
  for(model in models){
    mt<-as.data.frame(summary(model)$fixed[,1:4])
    mr<-as.data.frame(Reduce(rbind,summary(model)$random)[,1:4])
    row.names(mr)<-c("sd(Household)","sd(Community)")
    mt<-rbind(mt,mr)
    mt$Independent<-row.names(mt)
    mt$Independent<-gsub("PDSCont","Contagion Disgust",mt$Independent)
    mt$Independent<-gsub("PDSFood","Food Disgust",mt$Independent)
    mt$Independent<-gsub("PDSTotal","Total Disgust",mt$Independent)
    mt$Independent<-gsub("PDSPest","Component 3",mt$Independent)
    mf<-capture.output(summary(model)$formula)
    mt$Dependent<-NA
    mt$Dependent[1]<-c(strsplit(mf," ")[[1]][1])
    mf<-gsub("PDSCont","Contagion",mf)
    mf<-gsub("PDSFood","Food",mf)
    mf<-gsub("PDSTotal","Disgust",mf)
    mf<-gsub("PDSPest","Comp3",mf)
    mf<-gsub("Family","Household",mf)
    if(is.na(mts[1])) mts<-mt else mts<-rbind(mts,mt)
    if(is.na(mfs[1])) mfs<-mf else mfs<-c(mfs,mf)
  }
  
  mfs<-c("Model Formula:",mfs)
  ft<-flextable(mts,col_keys = c("Dependent","Independent","Estimate","l-95% CI","u-95% CI"))
  ft<-colformat_num(ft,j=c("Estimate","l-95% CI","u-95% CI"))
  ft<-autofit(ft)
  nr<-nrow(mt)
  ft<-hline(ft,i=seq(nr,nr*length(models),nr), border = fp_border(color="black"))
  ft<-font(ft,fontname="Calibri",part="all")
  ft<-fontsize(ft,size=11)
  ft<-bold(ft,part="header")
  ft
}
load("simpleREmodels.Rdata")
reC2<-modelflex2(list(modTI1x,modTP1x))
save_as_docx(reC2, path = "RE Model Table Single Component.docx")



