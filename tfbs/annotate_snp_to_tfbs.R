library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(GenomicRanges)
library(readr)

#ov.all<-data.frame()
#for (i in 1:23){
#cat(i,"\n")
#load(paste0("chr",i,".Robj"))
#ov.all<-rbind(ov.all,ov)
#}
#save(ov.all,file="tfbsbycpg.Robj")

#load("../07_enrichments/tfbsbycpg.Robj")
#path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/encode_tfbs"
#l<-list.files(path)
#l2<-gsub("tfbs_","",l)
#names(ov.all)[6:ncol(ov.all)]<-l2
#save(ov.all,file="tfbsbycpg.Robj")

#nvars<-names(ov.all)

#ov.all<-data.frame(ov.all)
#rs<-rowSums(ov.all[,-1:-6])
#w<-which(rs>0)
#rs[w]<-"tfbs"
#rs[-w]<-"no tfbs"

#tfbs<-data.frame(cpg=ov.all$cpg,tfbs=rs)

load("~/repo/godmc_phase2_analysis/results/enrichments/snpcontrolsets_selection_se.rdata")
f.all<-data.table(f.all)
f.all$snpchr<-gsub("chrX","chr23",f.all$snpchr)
temp<-data.frame()

for (chr in 1:23) {
snp_cis<-read_tsv(paste0("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/snpcpgpval.chr",chr,".cistrans3.txt.gz"))
snp_cis<-data.frame(snp_cis)
snpchr<-paste0("chr",chr)
cat(snpchr,"\t")
w<-which(f.all$snpchr==snpchr)
f.all.chr<-f.all[w,]
m<-match(f.all.chr$snp_tested,snp_cis$snp)
f.all.chr<-data.frame(f.all.chr,id=snp_cis$id[m],cpg=snp_cis$cpg[m])
temp<-rbind(temp,f.all.chr)
write.csv(f.all.chr,file=paste0(snpchr,".csv"),col.names=F,quote=F)
}

table(f.all$snpchr)
table(temp$snpchr)
#write.csv(t(names(f.all.chr)),"header.csv",quote=F)

#cat header.csv > allchrs.csv
#cat chr*.csv >> allchrs.csv
f.all<-temp
f.all$snpchr<-gsub("chr23","chrX",f.all$snpchr)

rm(temp,f.all.chr,snp_cis)
#f.all<-read_csv("allchrs.csv")

save(f.all,file="~/repo/godmc_phase2_analysis/results/enrichments/snpcontrolsets_selection_se.rdata")
##

gr_snp = with(f.all,GRanges(seqnames=snpchr,ranges=IRanges(min,max),strand=Rle("+")))

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
l<-l[w]

tfbs<-data.frame(f.all$SNP)
for (i in 1:length(l)){
#for (i in 1:10){
cat(i,"\n")
bed<-read.table(paste0(path,"/",l[i]))
bed$V1<-gsub("chrX","chr23",bed$V1)
#bed<-bed[which(bed$V1==chr),]

if(nrow(bed)>0){
bed$V2<-as.numeric(bed$V2)
bed$V3<-as.numeric(bed$V3)
bed<-unique(data.frame(chr=bed$V1,start=bed$V2,end=bed$V3,strand="*",antibody=bed$V4,celltype=bed$V5))
gr_range<-makeGRangesFromDataFrame(bed, keep.extra.columns=TRUE,starts.in.df.are.0based=TRUE) 
overlap<-countOverlaps(gr_snp,gr_range)
print(length(which(overlap==1)))
w<-which(overlap>0)
if(length(w)>0){overlap[w]<-1}
}

if(nrow(bed)==0){overlap<-rep(0,nrow(r))}

tfbs<-data.frame(tfbs,overlap)
}

names(tfbs)<-c("SNP",paste0(id$Celltype,"_",id$Type))

tfbs.sum<-rowSums(tfbs[,-1])
w<-which(tfbs.sum>0)
if(length(w)>0){tfbs.sum[w]<-1}
tfbs$any_tfbs<-tfbs.sum

g<-grep("CTCF",names(tfbs))
tfbs.ctcf<-rowSums(tfbs[,g])
w<-which(tfbs.ctcf>0)
if(length(w)>0){tfbs.ctcf[w]<-1}
tfbs$any_ctcf<-tfbs.ctcf


save(tfbs,file="tfbsbysnp.Robj")

m<-match(f.all$SNP,tfbs$SNP)
f.all<-data.frame(f.all,any_ctcf=tfbs$any_ctcf[m],any_tfbs=tfbs$any_tfbs[m])
save(f.all,file="snpcontrolsets_selection_se_tfbs.rdata")

library(data.table)
load("snpcontrolsets_selection_se_tfbs.rdata")
r<-fread("~/gerp/ex1.hg19_gerp++gt2_dropped",sep="\t")
spl<-strsplit(r$V3,split=" ")
spl<-do.call("rbind",spl)
r<-data.frame(SNP=as.character(spl[,6]),gerp_gt2=r$V2)

m<-match(f.all$SNP,as.character(r$SNP))
f.all<-data.frame(f.all,gerp_gt2=r[m,"gerp_gt2"])
save(f.all,file="snpcontrolsets_selection_se_tfbs.rdata")



