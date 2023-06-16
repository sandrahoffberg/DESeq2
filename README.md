## DESeq2 analysis in Rmarkdown

Differential expression analysis is used to identify differences in the transcriptome (gene expression) across a cohort of samples. Often, it will be used to define the differences between multiple biological conditions (e.g. drug treated vs. untreated samples). 


This capsule can be run either in an RStudio cloud workstation, or as a Reproducible Run. 


### Input data

This capsule downloads data during the R session. RNAseq data are obtained from the EBI Expression Atlas (GXA). Specifically data set E-GEOD-50760 which corresponds to PMID: 25049118. This data consists of 54 samples from 18 individuals. Each individual has a primary colorectal cancer sample, a metastatic liver sample, and a normal sample of the surrounding colonic epithilium. The quantification data required to run differential expression analysis using DEseq2 are raw readcounts for either genes or transcripts. 

Input data for DEseq2 consists of non-normalized sequence read counts at either the gene or transcript level. No preliminary normalization of this data is needed. DEseq2 will internally corrects for differences in library size, using the raw counts. The tool HTseq was used to obtain this information here.

### Output

- MA plot: MA plots display a log ratio (M) vs an average (A) in order to visualize the differences between two groups. In general we would expect the expression of genes to remain consistent between conditions and so the MA plot should be similar to the shape of a trumpet with most points residing on a y intercept of 0.

- PCA: Look at how samples group by treatment

- Volcano Plot: The volcano plot enables to simultaneously capture the effect size and significance of each tested gene.

- Plots by gene: plot the normalized counts for a single gene in order to get an idea of what is occurring for that gene across the sample cohort.


### Online Documentation

This DESeq2 demo follows tutorials available online

[Obtaining data and running DESeq2](https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/) can be found here. <br>
[Plotting](https://lashlock.github.io/compbio/R_presentation.html) can be found here. 