library(flextable)
library(officer)

modelflex<-function(model){
  mt<-as.data.frame(summary(model)$fixed[,1:4])
  mr<-as.data.frame(summary(model)$random$Family[,1:4])
  mt<-rbind(mt,mr)
  vnames<-row.names(mt)
  vnames<-gsub("PDSCont","Contagion",vnames)
  vnames<-gsub("PDSFood","Food",vnames)
  vars<-Reduce(rbind,strsplit(vnames,"_"))
  mt$Dependent<-vars[,1]
  mt$Independent<-vars[,2]
  vars2<-Reduce(rbind,strsplit(vnames,"[(_]"))
  mt$Dependent[grepl("sd\\(",vnames)]<-vars2[grepl("sd\\(",vnames),2]
  mt$Independent[grepl("sd\\(",vnames)]<-"sd(Household)"
  mt$Dependent[grepl("cor\\(",vnames)]<-NA
  mt$Independent[grepl("cor\\(",vnames)]<-"cor(Household)"
  
  mf<-capture.output(summary(model)$formula)
  mf<-c("Model Formula:",mf)
  mf<-gsub("PDSCont","Contagion",mf)
  mf<-gsub("PDSFood","Food",mf)
  mf<-gsub("Family","Household",mf)
  mf<-c(mf,"HH = Household mean, excluding target individual. Vi = Village mean, excluding target household","Items below the grey bar are group level effects for Household")
  
  ft<-flextable(mt,col_keys = c("Dependent","Independent","Estimate","l-95% CI","u-95% CI"))
  ft<-colformat_num(ft,j=c("Estimate","l-95% CI","u-95% CI"))
  ft<-autofit(ft)
  ft<-hline(ft,i=nrow(mt)-nrow(mr), border = fp_border(color="gray"))
  ft<-add_footer_lines(ft,values=mf)
  ft<-font(ft,fontname="Calibri",part="all")
  ft<-fontsize(ft,size=11)
  ft<-fontsize(ft, part="footer", size=10)
  ft<-padding(ft, part="footer", padding=0)
  ft<-bold(ft,part="header")
  ft
}

f11<-modelflex(mod11)
ff11<-modelflex(modF11)
fp11<-modelflex(modP11)
fpf11<-modelflex(modPF11)

save_as_docx("S2" = f11, "S3" = ff11, "S4" = fp11, "S5" = fpf11, path = "Model Tables.docx")

modelflex2<-function(models){
  mts<-NA
  mfs<-NA
  for(model in models){
    mt<-as.data.frame(summary(model)$fixed[,1:4])
    mr<-as.data.frame(Reduce(rbind,summary(model)$random)[,1:4])
    row.names(mr)<-c("sd(Household)","sd(Community)")
    mt<-rbind(mt,mr)
    mt$Independent<-row.names(mt)
    mt$Independent<-gsub("PDSCont","Contagion",mt$Independent)
    mt$Independent<-gsub("PDSFood","Food",mt$Independent)
  
    mf<-capture.output(summary(model)$formula)
    mt$Dependent<-NA
    mt$Dependent[1]<-c(strsplit(mf," ")[[1]][1])
    mf<-gsub("PDSCont","Contagion",mf)
    mf<-gsub("PDSFood","Food",mf)
    mf<-gsub("Family","Household",mf)
    if(is.na(mts[1])) mts<-mt else mts<-rbind(mts,mt)
    if(is.na(mfs[1])) mfs<-mf else mfs<-c(mfs,mf)
  }
  
  mfs<-c("Model Formula:",mfs)
  ft<-flextable(mts,col_keys = c("Dependent","Independent","Estimate","l-95% CI","u-95% CI"))
  ft<-colformat_num(ft,j=c("Estimate","l-95% CI","u-95% CI"))
  ft<-autofit(ft)
  ft<-hline(ft,i=c(6,12,18,24,30), border = fp_border(color="black"))
  ft<-font(ft,fontname="Calibri",part="all")
  ft<-fontsize(ft,size=11)
  ft<-bold(ft,part="header")
  ft
}

reC<-modelflex2(list(modC1x,modF1x,modP1x,modPF1x,mod1pestx,modP1pestx))
save_as_docx(reC, path = "RE Model Table.docx")
