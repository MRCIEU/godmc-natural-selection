#GoDMC_browser1.gz/tbi - tabixed results file from running the SMR method for all CpGs against 723 GWAS traits

#eQTLGen_browser1.gz/tbi - tabixed results file from running the SMR method for all genes against 723 GWAS traits

#GoDMC_eQTLGen_multi.msmr.gz - raw output file from the SMR method (multi-SNP) to investigate all cis-pairs for CpGs -> eGenes between the 2 consortia.

#Additionally there are 3 .RData files - genes.RData, CpGs.RData & traits.RData. The web app mainly used these as mapping files (e.g. each trait had a random UKB id for example, ENSG ids to gene symbols from Ensembl).


```{r}
library(data.table)
library(tidyverse)
datadir="/projects/MRC-IEU/research/projects/ieu2/p5/091/working/data/godmc_eqtlgen"
#load(paste0(datadir,"/","traits.RData"))

eqtl_mqtl_trait<-read.table("./data/eQTL_mQTL_trait.txt",sep="\t",he=T)
nrow(eqtl_mqtl_trait) #2962
length(unique(eqtl_mqtl_trait$cpg)) #1313
min(eqtl_mqtl_trait$p_both_topsnp)
#[1] "1.00E-07
length(which(is.na(eqtl_mqtl_trait$p_gene_topsnp)&is.na(eqtl_mqtl_trait$p_cpg_topsnp)))
length(which(is.na(eqtl_mqtl_trait$p_cpg_topsnp)))
```

Sometimes there are multiple SNPs for a CpG. Count how many heterogeneity Q stats are significant (FDR < 0.05) per CpG

```{r}
outcpg_3way <- eqtl_mqtl_trait %>% 
  group_by(cpg) %>%
  summarise(
    n=n(), 
    nsig=sum(p_both_topsnp < 1e-7),
    min=min(p_both_topsnp),
    max=max(p_both_topsnp)
  )
str(outcpg_3way)

subset(outcpg_3way, n != 1) %>%
  {hist(pmin(-log10(.$min) - -log10(.$max), 10), breaks=100)}
```

```{r}
subset(outcpg_3way, n != 1) %>%
  {hist(.$nsig / .$n)}
```

Create a relatively lenient variable for eqtl_mqtl_trait specificity where a CPG is cell-specific if any SNP has significant heterogeneity

```{r}
outcpg_3way$threeway <- outcpg_3way$nsig >= 1
table(outcpg_3way$threeway)
save(outcpg_3way, file="outcpg_3way.rdata")
```

eqtl_mqtl<-fread("./data/GoDMC_eQTLGen_multi.msmr.gz")
nrow(eqtl_mqtl) #1661343

length(which(is.na(eqtl_mqtl$p_Outco)&is.na(eqtl_mqtl$p_Expo)))
length(which(is.na(eqtl_mqtl$p_Outco)))
length(which(is.na(eqtl_mqtl$p_Expo)))
#72543

head(eqtl_mqtl[which(is.na(eqtl_mqtl$p_SMR)&is.na(eqtl_mqtl$p_Expo)),])
#      Expo_ID Expo_Chr  Expo_Gene Expo_bp        Outco_ID Outco_Chr
#1: cg23635910        1 cg23635910  912055 ENSG00000225880         1
#2: cg23635910        1 cg23635910  912055 ENSG00000228794         1
#3: cg05151709        1 cg05151709  846155 ENSG00000228794         1
#4: cg08546707        1 cg08546707  846195 ENSG00000228794         1
#5: cg03648020        1 cg03648020  931327 ENSG00000188976         1
#6: cg04891761        1 cg04891761 1001973 ENSG00000188976         1
#        Outco_Gene Outco_bp    topSNP topSNP_chr topSNP_bp A1 A2     Freq
#1: ENSG00000225880   762244 rs3131972          1    752721  A  G 0.161034
#2: ENSG00000228794   778907 rs3131972          1    752721  A  G 0.161034
#3: ENSG00000228794   778907 rs4970385          1    831489  C  T 0.296223
#4: ENSG00000228794   778907 rs4970385          1    831489  C  T 0.296223
#5: ENSG00000188976   887136 rs9442394          1   1006223  G  A 0.425447
#6: ENSG00000188976   887136 rs9442394          1   1006223  G  A 0.425447
#      b_Outco   se_Outco      p_Outco   b_Expo se_Expo p_Expo      b_SMR se_SMR
#1:  0.1504220 0.00953211 4.231491e-56 0.100712     NaN    NaN  1.4935900    NaN
#2: -0.0511553 0.00957391 9.131993e-08 0.100712     NaN    NaN -0.5079380    NaN
#3: -0.0559109 0.00933996 2.148002e-09 0.727730     NaN    NaN -0.0768292    NaN
#4: -0.0559109 0.00933996 2.148002e-09 0.727034     NaN    NaN -0.0769027    NaN
#5: -0.0386755 0.00901107 1.770710e-05 0.469885     NaN    NaN -0.0823084    NaN
#6: -0.0386755 0.00901107 1.770710e-05 0.165608     NaN    NaN -0.2335370    NaN
#   p_SMR  p_SMR_multi      p_HEIDI nsnp_HEIDI
#1:   NaN 2.033875e-05 7.087603e-05          5
#2:   NaN 2.762255e-05 1.592127e-04         10
#3:   NaN 1.286210e-05           NA         NA
#4:   NaN 1.286210e-05           NA         NA
#5:   NaN 1.422482e-05 1.335608e-15          8
#6:   NaN 2.693449e-05 3.664995e-09         14

outcpg_gen <- eqtl_mqtl %>% 
  group_by(Expo_ID) %>%
  summarise(
    n=n(), 
    nsig=sum(p_SMR < 1e-7),
    min=min(p_SMR),
    max=max(p_SMR)
  )
str(outcpg_gen)
nrow(outcpg_gen) #158730

subset(outcpg_gen, n != 1) %>%
  {hist(pmin(-log10(.$min) - -log10(.$max), 10), breaks=100)}
```

```{r}
subset(outcpg_gen, n != 1) %>%
  {hist(.$nsig / .$n)}
```
Create a relatively lenient variable for eQTL-mQTL specificity where a CPG is qtl-specific if any SNP has significant eqtl and mqtl

```{r}
outcpg_gen$specific <- outcpg_gen$nsig >= 1
#table(outcpg_gen$cpg_gen)
```

```{r}
save(outcpg_gen, file="outcpg_gen.rdata")
```

```{r}
load("GoDMC.Robj")

nrow(mqtl)
#[1] 136927575
length(which(is.na(mqtl$p_eQTL)&is.na(mqtl$p_SMR)))
#[1] 3519538

table(mqtl$nsnp_HEIDI)


#       3        4        5        6        7        8        9       10 
# 5322837  4751511  4369298  4016679  3883502  3506462  3482194  3337956 
#      11       12       13       14       15       16       17       18 
# 3145854  3090694  3078492  2924344  2932471  2762926  2590398  2636680 
#      19       20 
# 2548875 64787755

max(mqtl$p_eQTL,na.rm=T)
#[1] 4.99601e-08

mqtl<-mqtl[which(mqtl$p_eQTL<1e-8),]
outcpg <- mqtl %>% 
  group_by(ProbeID) %>%
  summarise(
    n=n(), 
    nsig=sum(p_SMR < 1e-7),
    min=min(p_SMR),
    max=max(p_SMR)
  )
str(outcpg)
nrow(outcpg) #178810

subset(outcpg, n != 1) %>%
  {hist(pmin(-log10(.$min) - -log10(.$max), 10), breaks=100)}
```

```{r}
subset(outcpg, n != 1) %>%
  {hist(.$nsig / .$n)}
```
Create a relatively lenient variable for eQTL-mQTL specificity where a CPG is qtl-specific if any SNP has significant eqtl and mqtl

```{r}
outcpg$specific <- outcpg$nsig >= 1
#table(outcpg$cpg_mqtl)
```

```{r}
save(outcpg, file="outcpg.rdata")
```

```
#zcat /projects/MRC-IEU/research/projects/ieu2/p5/091/working/data/godmc_eqtlgen/eQTLGen_browser1.gz |wc -l
#11767689
eqtl<-fread(paste0(datadir,"/eQTLGen_browser1.gz"))
n<-c("snpchr","snppos1","snppos2","ProbeID","SNP","SNP_Chr","SNP_bp","A1","A2","Freq","b_GWAS","se_GWAS","p_GWAS","b_eQTL","se_eQTL","p_eQTL","b_SMR","se_SMR","p_SMR","p_SMR_multi","p_HEIDI","nsnp_HEIDI","GWA_id")
names(eqtl)<-n
nrow(eqtl)
#[1] 11767689
length(which(is.na(eqtl$p_eQTL)))
#0
length(which(is.na(eqtl$p_SMR)))
#0
length(which(is.na(eqtl$p_GWA)))
#0

```
load(paste0(datadir,"/CpGs.RData"))
nrow(CpGs)
#[1] 181293