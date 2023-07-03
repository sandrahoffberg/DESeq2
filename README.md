[![Code Ocean Logo](images/CO_logo_135x72.png)](http://codeocean.com/product)

<hr>

# DESeq2 analysis and plotting

Differential expression analysis is used to identify differences in the transcriptome (gene expression) across a cohort of samples. Often, it will be used to define the differences between multiple biological conditions (e.g. drug treated vs. untreated samples). 

**The design formula must be edited as appropriate for the input data**:

The design formula expresses the condition we want to measure the effect of, controlling for batch differences. The formula should be a tilde (~) followed by the variables (columns of the metadata file) with plus signs between them (it will be coerced into an formula if it is not already). The design can be changed later, however then all differential analysis steps should be repeated, as the design formula is used to estimate the dispersions and to estimate the log2 fold changes of the model. DESeq2 can analyze any possible experimental design that can be expressed with fixed effects terms (multiple factors, designs with interactions, designs with continuous variables, splines, and so on are all possible).


For more information, see the [DESeq2 Vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).


## Input data

Input data for DEseq2 consists of non-normalized sequence read counts at either the gene or transcript level. No preliminary normalization of this data is needed. DEseq2 will internally corrects for differences in library size, using the raw counts. This capsule can accommodate 4 different types of input data: 

1. Transcript abundance data: <br>
from Salmon, Sailfish, kallisto, RSEM, etc. 

Some advantages of using the above methods for transcript abundance estimation are: (i) this approach corrects for potential changes in gene length across samples (e.g. from differential isoform usage) (Trapnell et al. 2013), (ii) some of these methods (Salmon, Sailfish, kallisto) are substantially faster and require less memory and disk usage compared to alignment-based methods that require creation and storage of BAM files, and (iii) it is possible to avoid discarding those fragments that can align to multiple genes with homologous sequence, thus increasing sensitivity (Robert and Watson 2015).  Multiple data files will be read in by giving the capsule a string found in the file names (i.e., "treated" can file treated.txt and untreated.txt). 


2. Counts data: <br>
from featurecounts, Seurat, etc.

The user should provide the counts matrix, the information about the samples (the columns of the count matrix) as a DataFrame or data.frame, and the design formula. It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. DESeq2 will not make guesses as to which column of the count matrix belongs to which row of the column data, these must be provided to DESeq2 already in consistent order.


3. htseq-count data: <br>
from HTSeq Python package

Multiple data files will be read in by giving the capsule a string found in the file names (i.e., "treated" can file treated.txt and untreated.txt). 


4. SummarizedExperiment

If one has already created or obtained a SummarizedExperiment, the single file with data and metadata can be read in. 


## Output

- CSV file with the results table

- MA plot: MA plots display a log ratio (M) vs an average (A) in order to visualize the differences between two groups. In general we would expect the expression of genes to remain consistent between conditions and so the MA plot should be similar to the shape of a trumpet with most points residing on a y intercept of 0.

- PCA: Visualize how the samples group by treatment

- Volcano Plot: The volcano plot enables to simultaneously capture the effect size and significance of each tested gene.

- Plots by gene: plot the normalized counts for a single gene in order to get an idea of what is occurring for that gene across the sample cohort.

## App Panel Parameters

1. type of input data: Transcipt abundance, Count matrix, htseq count, or summarized experiment.
2. File path to a single data file, for Count matrix or summarized experiment input data types. 
3. String to locate multiple input data files, for Transcript abundance or htseq count data. 
4. Table that correlated trabscipts to genes. 
5. Type of transcript analysis performed upstreamm e.g., kallisto, RSEM, Salmon, sailfish
6. Metadata file
7. Name of column in the metadata file that indicates which treatment the sample was in, e.g., condition, tissueType, dex
8. Minimum number of reads across all samples to include that sample in the data matrix
9. Alpha value for significance
10. Number of genes to save individual plots for
11. Name of output CSV file


## Cite 

Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi:10.1186/s13059-014-0550-8.



<hr>

[Code Ocean](https://codeocean.com/) is a cloud-based computational platform that aims to make it easy for researchers to share, discover, and run code.<br /><br />
[![Code Ocean Logo](images/CO_logo_68x36.png)](https://www.codeocean.com)