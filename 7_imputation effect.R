#Check correlations between imputed datasets
Inflams<-Reduce(cbind,lapply(dsets,function(x) x[,"Inflam"]))

cc<-cor(Inflams)
mean(cc[upper.tri(cc)])
sd(cc[upper.tri(cc)])

cc<-cor(Inflams[is.na(d2$IL6) | is.na(d2$CRP),])
mean(cc[upper.tri(cc)])
sd(cc[upper.tri(cc)])

Parasites<-Reduce(cbind,lapply(dsets,function(x) x[,"Parasites"]))
cc<-cor(Parasites)
mean(cc[upper.tri(cc)])
sd(cc[upper.tri(cc)])

cc<-cor(Parasites[is.na(d2$IgE),])
mean(cc[upper.tri(cc)])
sd(cc[upper.tri(cc)])
