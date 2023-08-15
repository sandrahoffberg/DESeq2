[![Code Ocean Logo](images/CO_logo_135x72.png)](http://codeocean.com/product)

<hr>

# DESeq2 analysis and plotting

Differential expression analysis is used to identify differences in the transcriptome (gene expression) across a cohort of samples. Often, it will be used to define the differences between multiple biological conditions (e.g. drug treated vs. untreated samples). 

**The design formula must be edited as appropriate for the input data**:

The design formula expresses the condition we want to measure the effect of, controlling for batch differences. The formula should be a tilde (~) followed by the variables (columns of the metadata file) with plus signs between them (it will be coerced into an formula if it is not already). The design can be changed later, however then all differential analysis steps should be repeated, as the design formula is used to estimate the dispersions and to estimate the log2 fold changes of the model. DESeq2 can analyze any possible experimental design that can be expressed with fixed effects terms (multiple factors, designs with interactions, designs with continuous variables, splines, and so on are all possible).


For more information, see the [DESeq2 Vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).


## Input data

Input data for DEseq2 consists of non-normalized sequence read counts at either the gene or transcript level. No preliminary normalization of this data is needed. DEseq2 internally corrects for differences in library size, using the raw counts. This capsule can accommodate 4 different types of input data: 

1. Transcript abundance data: <br>
from Salmon, Sailfish, kallisto, RSEM, etc. 

Some advantages of using the above methods for transcript abundance estimation are: (i) this approach corrects for potential changes in gene length across samples (e.g. from differential isoform usage) (Trapnell et al. 2013), (ii) some of these methods (Salmon, Sailfish, kallisto) are substantially faster and require less memory and disk usage compared to alignment-based methods that require creation and storage of BAM files, and (iii) it is possible to avoid discarding those fragments that can align to multiple genes with homologous sequence, thus increasing sensitivity (Robert and Watson 2015).  Multiple data files will be read in by giving the capsule a string found in the file names (i.e., "treated" can file treated.txt and untreated.txt). 

If the transcript files are autodetected, they must contain "abundance" in the name. 


2. Counts data: <br>
from featurecounts, Seurat, etc.

The user should provide (1) the counts matrix , (2) information about the samples (the columns of the count matrix) as a DataFrame or data.frame, and (3) the design formula. It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. DESeq2 will not make guesses as to which column of the count matrix belongs to which row of the column data, these must be provided to DESeq2 already in consistent order.

If the counts matrix is autodetected, it must contain "raw" in the name. 

3. htseq-count data: <br>
from HTSeq Python package

Multiple data files will be read in by giving the capsule a string found in the file names (i.e., "treated" can file treated.txt and untreated.txt). 


4. SummarizedExperiment

If one has already created or obtained a SummarizedExperiment, the single file with data and metadata can be read in. 


## Output

- **DESeq2_results.csv**: CSV file with the results table

- **MA_plot.png**: MA plots display a log ratio (M) vs an average (A) in order to visualize the differences between two groups. In general we would expect the expression of genes to remain consistent between conditions and so the MA plot should be similar to the shape of a trumpet with most points residing on a y intercept of 0.

- **PCA.png**: Visualize how the samples group by treatment

- **volcano_plot.png**: The volcano plot enables to simultaneously capture the effect size and significance of each tested gene.

- **/plots_by_gene**: A folder containing a file for each gene that plots the normalized counts for a single gene in order to get an idea of what is occurring for that gene across the sample cohort.

## App Panel Parameters

Type of input data
- Transcipt abundance, Count matrix, htseq count, or summarized experiment.

File path to a single data file
- For Count matrix or summarized experiment input data types. 

Key word contained in multiple input files
- for Transcript abundance or htseq count data. 

Table that matches transcripts to genes.
- for Transcipt abundance data. In CSV format, column 1 should be the transcript name and column 2 should be the Gene ID.  If the table is automatically detected, it must contain "gene" in the name. See the example below

TXNAME,GENEID <br>
ENST00000456328.2,ENSG00000223972.5 <br>
ENST00000450305.2,ENSG00000223972.5 <br> 
ENST00000473358.1,ENSG00000243485.5 <br>

Type of transcript analysis
- What anlaysis was performed updstream e.g., kallisto, RSEM, Salmon, Sailfish

Path to Sample Metadata 
- Metadata file, required for Transcript Abundance, HTseq-data, and Counts data. The format of the metadata file is flexible, but it should be a CSV file that lists the condition corresponding to each sample. The design formula reflects the condition, and must be edited as appropriate. If the metadata file is automatically detected, it must contain "Metadata" or "metadata" in the name.

Name of treatment/control column
- Column in the metadata file that indicates which treatment the sample was in, e.g., condition, tissueType, dex.  [Default: condition]

Minimum number of reads
- Minimum reads across all samples to include that sample in the data matrix [Default: 5]

Alpha value for significance
- the significance cutoff used for optimizing the independent filtering [Default: 0.1]. If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.

Number of genes to save individual plots for
- The # most significant genes to save plots for in the **/plots_by_gene** directory. [Default: 6]

Name of output file
- Name of CSV file without file extension. [Default: DESeq2_results]


## Cite 

Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi:10.1186/s13059-014-0550-8.

<hr>

[Code Ocean](https://codeocean.com/) is a cloud-based computational platform that aims to make it easy for researchers to share, discover, and run code.<br /><br />
[![Code Ocean Logo](images/CO_logo_68x36.png)](https://www.codeocean.com)