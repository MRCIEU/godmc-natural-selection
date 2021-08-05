#zcat GoDMC_browser1.gz |wc -l
#136927575
#mqtl<-fread("./data/GoDMC_browser1.gz")
#Avoidable 1.287 seconds. This file is very unusual: it ends abruptly without a final newline, and also its size is a multiple of 4096 bytes. Please properly end the last row with a newline using for example 'echo >> file' to avoid this  time to copy.
#|--------------------------------------------------|
#|==================================================|
#There were 50 or more warnings (use warnings() to see the first 50)

#nrow(mqtl) #8571321

#awk '{print $0 >> $1".txt"}' GoDMC_browser1
#wc -l [0-9]*txt # 136927575

library(data.table)
n<-c("snpchr","snppos1","snppos2","ProbeID","SNP","SNP_Chr","SNP_bp","A1","A2","Freq","b_GWAS","se_GWAS","p_GWAS","b_eQTL","se_eQTL","p_eQTL","b_SMR","se_SMR","p_SMR","p_SMR_multi","p_HEIDI","nsnp_HEIDI","GWA_id")
#cls<-cols(snpchr = col_double(),snppos1 = col_double(),snppos2 = col_double(),ProbeID = col_character(),SNP = col_character(),SNP_Chr = col_double(),SNP_bp = col_double(), A1 = col_character(),A2 = col_character(),Freq = col_double(),b_GWAS = col_double(),se_GWAS = col_double(),p_GWAS = col_double(),b_eQTL = col_double(),se_eQTL = col_double(),p_eQTL = col_double(),b_SMR = col_double(),se_SMR = col_double(),p_SMR = col_double(),p_SMR_multi = col_double(),p_HEIDI = col_double(),nsnp_HEIDI = col_double(),GWA_id = col_character())

#mqtl<-read_delim("./data/GoDMC_browser1.gz",col_names=n,col_types = cls)
#nrow(mqtl) #8568524

mqtl<-data.frame()
for (i in 1:22){
cat(i,"\n")
chr<-fread(paste0("./data/",i,".txt"))
cat(nrow(chr),"\n")
mqtl<-rbind(mqtl,chr)
}
nrow(mqtl)

names(mqtl)<-n
save(mqtl,file="GoDMC.Robj")
