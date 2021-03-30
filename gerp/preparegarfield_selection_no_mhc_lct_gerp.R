arguments<-commandArgs(T)
i<-as.numeric(arguments[1])
 
###
load("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/13_selection/snpcontrolsets_selection_se_tfbs.rdata")
length(unique(f.all$SNP))

#[1] 10085072
 
table(f.all$mqtl_clumped)

#  FALSE    TRUE 
#5980960  224648 

w1<-which(f.all$snpchr=="chr6"&f.all$snppos>24570005&f.all$snppos<38377657)
w2<-which(f.all$snpchr=="chr2"&f.all$snppos>129608646&f.all$snppos<143608646)

w<-c(w1,w2)
if(length(w)>0){
f.all<-f.all[-w,]
}

w<-which(f.all$snptype%in%"INDEL")
indels<-f.all[w,"snppos"]
f.all<-f.all[-w,]

length(which(!is.na(f.all$gerp_gt2)))
#[1] 548495

df.out<-read.table(paste("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/snpcpgpval.chr",i,".cistrans.txt.gz",sep=""),he=T)

w1<-which(df.out$snpchr=="chr6"&df.out$snppos>24570005&df.out$snppos<38377657)
w2<-which(df.out$snpchr=="chr2"&df.out$min>129608646&df.out$max<143608646)

w<-c(w1,w2)
if(length(w)>0){
df.out<-df.out[-w,]
}

w<-which(df.out$type%in%"INDEL")
if(length(w)>0){
df.out<-df.out[-w,]
}


bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered.bim")
g<-grep("INDEL",bim$V2)
bim<-bim[-g,]

#flip<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/ref/flipped_snps.txt",he=F)
#w<-which(bim$V2%in%flip[,1])
#bim<-bim[-w,]

#indels<-read.table("/panfs/panasas01/shared-godmc/INDELs/indels_equal_seq_length.txt")
#w<-which(bim$V2%in%indels[,1])
#bim<-bim[-w,]

 #ALL-cis+trans
bim<-bim[bim$V1==i,]
m<-match(bim$V4,df.out$snppos)
df.out2<-data.frame(snp=bim$V4,pval=df.out$min_pval[m])

w<-which(is.na(df.out2$pval))
df.out2$pval[w]<-1
 
length(which(!is.na(df.out2$pval)))

chr<-paste("chr",i,sep="")
f.all.chr<-f.all[f.all$snpchr==chr,]
m<-match(df.out2$snp,f.all.chr$snppos)
f.all2<-f.all.chr[m,c("snppos","gerp_gt2")]

w<-which(is.na(f.all2$gerp_gt2))
f.all2$gerp<-f.all2$gerp_gt2
f.all2$gerp[w]<-0
f.all2$gerp[-w]<-1

df.out3<-df.out2
o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/mqtl_gerp_no_mhc_lct/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)


#f.all2$annot<-paste(f.all2$ihs,f.all2$fst,f.all2$xpehhchb,f.all2$xpehhyri,f.all2$sds,sep="")
f.all2$annot<-f.all2$gerp

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation/"
p<-paste("chr",i,sep="")
r<-read.table(paste(path,"chr",i,sep=""))

if(i==6){
w<-which(r$V1>24570005&r$V1<38377657)}

if(i==2){
w<-which(r$V1>129608646&r$V1<143608646)}

if(length(w)>0){
r<-r[-w,]
}

w<-which(r$V1%in%indels)
r<-r[-w,]

m<-match(r$V1,f.all2$snppos)
f.all2<-f.all2[m,]
f.all2$snppos<-r$V1
w<-which(is.na(m))
f.all2$annot[w]<-0

#y<-nchar(f.all2$annot)
#m<-max(nchar(f.all2$annot)) #
#m<-which(y==m)
#f.all2$annot[w]<-f.all2$annot[m[1]] #00000

write.table(data.frame(f.all2$snppos,f.all2$annot),paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation_gerp_no_lct_mhc/chr",i,sep=""),sep=" ",col.names=F,row.names=F,quote=F)

##
#trans
#set cispval on 1

w<-which(df.out$snp_cis!="FALSE")
#table(df.out[w,"cis"],df.out[w,"trans"])
df.out2<-df.out
df.out2[w,"pval"]<-1

o<-order(df.out2$pval)
df.out2<-df.out2[o,]

bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered.bim")
g<-grep("INDEL",bim$V2)
bim<-bim[-g,]

bim<-bim[bim$V1==i,]

m<-match(bim$V4,df.out2$snppos)
df.out2<-data.frame(snp=bim$V4,pval=df.out2$pval[m])
w<-which(is.na(df.out2$pval))
df.out2$pval[w]<-1
 
length(which(!is.na(df.out2$pval)))

chr<-paste("chr",i,sep="")
f.all.chr<-f.all[f.all$snpchr==chr,]
m<-match(df.out2$snp,f.all.chr$snppos)
f.all2<-f.all.chr[m,c("snppos","gerp_gt2")]

w<-which(is.na(f.all2$gerp_gt2))
f.all2$gerp<-f.all2$gerp_gt2
f.all2$gerp[w]<-0
f.all2$gerp[-w]<-1

df.out3<-df.out2
o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/transmqtl_gerp_no_lct_mhc/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)


f.all2$annot<-f.all2$gerp

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation/"
p<-paste("chr",i,sep="")
r<-read.table(paste(path,"chr",i,sep=""))

if(i==6){
w<-which(r$V1>24570005&r$V1<38377657)}

if(i==2){
w<-which(r$V1>129608646&r$V1<143608646)}

if(length(w)>0){
r<-r[-w,]
}

w<-which(r$V1%in%indels)
r<-r[-w,]

m<-match(r$V1,f.all2$snppos)
f.all2<-f.all2[m,]
f.all2$snppos<-r$V1
w<-which(is.na(m))
f.all2$annot[w]<-0


#y<-nchar(f.all2$annot)
#m<-max(nchar(f.all2$annot))
#m<-which(y==m)
#f.all2$annot[w]<-f.all2$annot[m[1]]

write.table(data.frame(f.all2$snppos,f.all2$annot),paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation_gerp_trans_no_lct_mhc/chr",i,sep=""),sep=" ",col.names=F,row.names=F,quote=F)


#
#cis
#set trans and ambivalent pval on 1

w<-which(df.out$snp_cis!="TRUE")
#table(df.out[w,"cis"],df.out[w,"trans"])
df.out2<-df.out
df.out2[w,"pval"]<-1

o<-order(df.out2$pval)
df.out2<-df.out2[o,]

bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered.bim")
g<-grep("INDEL",bim$V2)
bim<-bim[-g,]

bim<-bim[bim$V1==i,]
m<-match(bim$V4,df.out2$snppos)
df.out2<-data.frame(snp=bim$V4,pval=df.out2$pval[m])
w<-which(is.na(df.out2$pval))
df.out2$pval[w]<-1
 
length(which(!is.na(df.out2$pval)))

chr<-paste("chr",i,sep="")
f.all.chr<-f.all[f.all$snpchr==chr,]
m<-match(df.out2$snp,f.all.chr$snppos)
f.all2<-f.all.chr[m,c("snppos","gerp_gt2")]

w<-which(is.na(f.all2$gerp_gt2))
f.all2$gerp<-f.all2$gerp_gt2
f.all2$gerp[w]<-0
f.all2$gerp[-w]<-1

df.out3<-df.out2
o<-order(df.out3$snp)
df.out3<-df.out3[o,]

write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/cismqtl_gerp_no_lct_mhc/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)


f.all2$annot<-f.all2$gerp

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation/"
p<-paste("chr",i,sep="")
r<-read.table(paste(path,"chr",i,sep=""))
if(i==6){
w<-which(r$V1>24570005&r$V1<38377657)}

if(i==3){
w<-which(r$V1>129608646&r$V1<143608646)}

if(length(w)>0){
r<-r[-w,]
}

w<-which(r$V1%in%indels)
r<-r[-w,]

m<-match(r$V1,f.all2$snppos)
f.all2<-f.all2[m,]
f.all2$snppos<-r$V1
w<-which(is.na(m))
f.all2$annot[w]<-0

#y<-nchar(f.all2$annot)
#m<-max(nchar(f.all2$annot))
#m<-which(y==m)
#f.all2$annot[w]<-f.all2$annot[m[1]]

write.table(data.frame(f.all2$snppos,f.all2$annot),paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation_gerp_cis_no_lct_mhc/chr",i,sep=""),sep=" ",col.names=F,row.names=F,quote=F)

###

#
#ambivalent
#set trans and cis pval on 1

w<-which(df.out$snp_cis!="ambivalent")
#table(df.out[w,"cis"],df.out[w,"trans"])
df.out2<-df.out
df.out2[w,"pval"]<-1

o<-order(df.out2$pval)
df.out2<-df.out2[o,]

bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered.bim")
g<-grep("INDEL",bim$V2)
bim<-bim[-g,]

bim<-bim[bim$V1==i,]
m<-match(bim$V4,df.out2$snppos)
df.out2<-data.frame(snp=bim$V4,pval=df.out2$pval[m])
w<-which(is.na(df.out2$pval))
df.out2$pval[w]<-1
 
length(which(!is.na(df.out2$pval)))

chr<-paste("chr",i,sep="")
f.all.chr<-f.all[f.all$snpchr==chr,]
m<-match(df.out2$snp,f.all.chr$snppos)
f.all2<-f.all.chr[m,c("snppos","gerp_gt2")]

w<-which(is.na(f.all2$gerp_gt2))
f.all2$gerp<-f.all2$gerp_gt2
f.all2$gerp[w]<-0
f.all2$gerp[-w]<-1

df.out3<-df.out2
o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/ambivalentmqtl_gerp_no_lct_mhc/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)


f.all2$annot<-f.all2$gerp

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation/"
p<-paste("chr",i,sep="")
r<-read.table(paste(path,"chr",i,sep=""))
if(i==6){
w<-which(r$V1>24570005&r$V1<38377657)}

if(i==3){
w<-which(r$V1>129608646&r$V1<143608646)}

if(length(w)>0){
r<-r[-w,]
}
w<-which(r$V1%in%indels)
r<-r[-w,]

m<-match(r$V1,f.all2$snppos)
f.all2<-f.all2[m,]
f.all2$snppos<-r$V1
w<-which(is.na(m))
f.all2$annot[w]<-0


#y<-nchar(f.all2$annot)
#m<-max(nchar(f.all2$annot))
#m<-which(y==m)
#f.all2$annot[w]<-f.all2$annot[m[1]]

write.table(data.frame(f.all2$snppos,f.all2$annot),paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation_gerp_ambivalent_no_lct_mhc/chr",i,sep=""),sep=" ",col.names=F,row.names=F,quote=F)
