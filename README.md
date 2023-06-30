[![Code Ocean Logo](images/CO_logo_135x72.png)](http://codeocean.com/product)

<hr>

# DESeq2

Differential expression analysis is used to identify differences in the transcriptome (gene expression) across a cohort of samples. Often, it will be used to define the differences between multiple biological conditions (e.g. drug treated vs. untreated samples). 


This capsule can be run either in an RStudio cloud workstation, or as a Reproducible Run. 


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "../results/")
```

## DESEQ2 R Tutorial


#### Load libraries that are pre-installed on the Environment UI page.

### Input data

This capsule can accommodate 4 different types of input data: 

1. Transcript abundance data: 
from Salmon, Sailfish, kallisto, RSEM, etc. 

Some advantages of using the above methods for transcript abundance estimation are: (i) this approach corrects for potential changes in gene length across samples (e.g. from differential isoform usage) (Trapnell et al. 2013), (ii) some of these methods (Salmon, Sailfish, kallisto) are substantially faster and require less memory and disk usage compared to alignment-based methods that require creation and storage of BAM files, and (iii) it is possible to avoid discarding those fragments that can align to multiple genes with homologous sequence, thus increasing sensitivity (Robert and Watson 2015).


2. Counts data: 
the user should provide the counts matrix, the information about the samples (the columns of the count matrix) as a DataFrame or data.frame, and the design formula.
from featurecounts, Seurat, etc


It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. DESeq2 will not make guesses as to which column of the count matrix belongs to which row of the column data, these must be provided to DESeq2 already in consistent order.


3. htseq-count data:

4. 









4This capsule downloads data during the R session. RNAseq data are obtained from the EBI Expression Atlas (GXA). Specifically data set E-GEOD-50760 which corresponds to PMID: 25049118. This data consists of 54 samples from 18 individuals. Each individual has a primary colorectal cancer sample, a metastatic liver sample, and a normal sample of the surrounding colonic epithilium. The quantification data required to run differential expression analysis using DEseq2 are raw readcounts for either genes or transcripts. 

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



<hr>

[Code Ocean](https://codeocean.com/) is a cloud-based computational platform that aims to make it easy for researchers to share, discover, and run code.<br /><br />
[![Code Ocean Logo](images/CO_logo_68x36.png)](https://www.codeocean.com)