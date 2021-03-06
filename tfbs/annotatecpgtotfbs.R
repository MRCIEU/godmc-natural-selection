#generatetfbsbedfiles.R
#meancpgbymqtl.R

#collapse overlaps
library(dplyr)
path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/encode_tfbs"
l<-list.files(path)
l2<-gsub(".bed","",l)
l2<-gsub("tfbs_","",l2)

id<-read.table("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation-tfbs/link_file_tfbs2.txt",he=T)
id$Type<-gsub("_.*","",id$Type)
m<-match(id$Annotation,l2)
id<-id[m,]

w<-which(id$Tissue=="Blood")
id<-id[w,]
l.bl<-l2[w]

load("~/repo/godmc_phase2_analysis/07_enrichments/tfbsbycpg.Robj")
names(ov.all)[6:ncol(ov.all)]<-l2
w<-which(names(ov.all)%in%l.bl)

tfbs<-data.frame()

ov.all<-data.frame(ov.all)
rs<-rowSums(ov.all[,w])
tfbs<-data.frame(cpg=ov.all$cpg,tfbs_rs=rs,tfbs=0)

tfbs2 <- tfbs %>% group_by(cpg) %>% summarize(max.tfbs.chr=max(tfbs_rs))
w<-which(tfbs2$max.tfbs.chr>0)
tfbs2$tfbs<-0
tfbs2$tfbs[w]<-"tfbs"
tfbs2$tfbs[-w]<-"no tfbs"

save(tfbs2,file="tfbsbycpg_summarised.Robj")

