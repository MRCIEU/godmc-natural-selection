# Analysis of natural selection on DNAm using GoDMC

Work in progress.

Project folder (https://ieuportal.epi.bris.ac.uk/project_detail/627d97b1-d9d8-417a-ae28-832e866ecc6c/):

```
/projects/MRC-IEU/research/projects/ieu2/p5/091 
```

Paper: https://outlook.office.com/mail/group/groups.bristol.ac.uk/grp-godmc-natural_selection/files

Github repo for method: https://github.com/danjlawson/robustarchitecture


## Todo

1. Move analysis scripts to https://github.com/MRCIEU/godmc-natural-selection and data to project folder - https://ieuportal.epi.bris.ac.uk/project_detail/627d97b1-d9d8-417a-ae28-832e866ecc6c/
7. Add GERP score to Figure 1 and Sfig1 (Josine) 
8. Josine to send Gib the cross-tissue effect estimates. Gib to estimate heterogeneity between tissues per mQTL. Josine to regression heterogeneity against selection scores
9. Compare enrichments for selection scores for:

-40k DNAm-trait pairs

-5.5k eQTL trait pairs

-3000 mQTL-eQTL trait pairs

-non colocalizing mQTL

-non-localizing eQTL

Data can be found here:
```
/projects/MRC-IEU/research/projects/ieu2/p5/091/working/data/godmc_eqtlgen
```

## Done

1. Do separate slopes for SNPs overlapping TFBS vs no overlap (Josine) (DONE)
3. Do separate S fitting for CpGs overlapping TFBS vs no overlap – remove SNPs that associate with multiple DNAm sites that fall in the two categories: 
    - Get data (DONE, Josine) 
    - Run analysis (DONE, Dan)  
    - Make github repo for method (DONE, Dan) 
4. S search to be either more granular or a true heuristic search (DONE, Dan) 
5. Flip signs on S to be intuitive (DONE, -S = negative selection) 
6. Create figure 1 from red/blue/black with confidence intervals around slope (Dan) 
