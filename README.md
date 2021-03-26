# Analysis of natural selection on DNAm using GoDMC

Work in progress.

Project folder (https://ieuportal.epi.bris.ac.uk/project_detail/627d97b1-d9d8-417a-ae28-832e866ecc6c/):

```
/projects/MRC-IEU/research/projects/ieu2/p5/091 
```


## Todo


1. Do separate slopes for SNPs overlapping TFBS vs no overlap (Josine) (DONE)
2. Move analysis scripts to https://github.com/MRCIEU/godmc-natural-selection and data to project folder - https://ieuportal.epi.bris.ac.uk/project_detail/627d97b1-d9d8-417a-ae28-832e866ecc6c/
3. Do separate S fitting for CpGs overlapping TFBS vs no overlap – remove SNPs that associate with multiple DNAm sites that fall in the two categories: 
    - Get data (Josine) 
    - Run analysis (Dan)  
    - Make github repo for method (Dan) 
4. S search to be either more granular or a true heuristic search (Dan) 
5. Flip signs on S to be intuitive (-S = negative selection) 
6. Create figure 1 from red/blue/black with confidence intervals around slope (Dan) 
7. Add GERP score to Figure 1 and Sfig1 (Josine) 
8. Estimate the enrichment of high SDS/Fst SNPs in complex traits - are they more likely to be disease traits - use file here on bc3: /newshared/godmc/1kg_reference_ph3/mQTLSNP_selection_se_tfbs.Robj.  (Gib) 
