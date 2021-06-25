# Heterogeneity analysis of per-tissue mQTL associations

To run

```
cp /projects/MRC-IEU/research/projects/ieu2/p5/091/working/data/tissue_specificity/tissue.RData .
Rscript heterogeneity.r
Rscript -e "rmarkdown::render('heterogeneity.rmd')"
```

The `heterogeneity.rdata` file contains an `out` object which is a data frame of the heterogeneity test statistics per DNAm

The `outcpg.rdata` file contains a per-CpG summary of whether it is tissue specific. See `heterogeneity.nb.html` for details.
