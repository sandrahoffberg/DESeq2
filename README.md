[![Code Ocean Logo](images/CO_logo_135x72.png)](http://codeocean.com/product)

<hr>

# DESeq2 analysis and plotting

Differential expression analysis is used to identify differences in the transcriptome (gene expression) across a cohort of samples. Often, it will be used to define the differences between multiple biological conditions (e.g. drug treated vs. untreated samples). 

For more information, see the [DESeq2 Vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

For example parameters that work with the sample data included in this capsule, please reference [deseq_sample_arguments](deseq_sample_arguments.csv).

## Input data

Input data for DEseq2 consists of non-normalized sequence read counts at either the gene or transcript level. No preliminary normalization of this data is needed. DEseq2 internally corrects for differences in library size using the raw counts. This capsule can accommodate 4 different types of input data: 

1. Transcript abundance data: <br>
from Salmon, Sailfish, kallisto, RSEM, etc. 

Some advantages of using the above methods for transcript abundance estimation are: (i) this approach corrects for potential changes in gene length across samples (e.g. from differential isoform usage) (Trapnell et al. 2013), (ii) some of these methods (Salmon, Sailfish, kallisto) are substantially faster and require less memory and disk usage compared to alignment-based methods that require creation and storage of BAM files, and (iii) it is possible to avoid discarding those fragments that can align to multiple genes with homologous sequence, thus increasing sensitivity (Robert and Watson 2015).  Multiple data files will be read in by giving the capsule a string found in the file names (i.e., "treated" can find treated.txt and untreated.txt). 

If the transcript files are to be autodetected, they must contain "abundance" in the name. 


2. Counts data: <br>
from featurecounts, Seurat, etc.

The user should provide (1) the counts matrix , (2) information about the samples (the columns of the count matrix) as a DataFrame or data.frame, and (3) the design formula. It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. DESeq2 will not make guesses as to which column of the count matrix belongs to which row of the column data, these must be provided to DESeq2 already in consistent order.

If the counts matrix is to be autodetected, it must contain "raw" in the name. 

3. htseq-count data: <br>
from HTSeq Python package

Multiple data files will be read in by giving the capsule a string found in the file names (i.e., "treated" can find treated.txt and untreated.txt). 


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
- **Transcipt abundance**, **Counts data**, **HTseq data**, or **summarized experiment**. [Default: Transcript_abundance]

Input data directory
- Directory to search for input data, e.g., /data/HTseq_data. If not provided, this will default to "../data/{input_data_type}". For example, ../data/Transcript_abundance when the **Type of input data** is "Transcript_abundance" [Default: ../data/Transcript_abundance]

File path to a single data file
- For **Counts_data** or **Summarized_experiment** input data types. If not specified, first file containing "raw" will be chosen. [Default: ../data/Counts_data/raw-counts.csv]

Key word contained in multiple input files
- for **Transcript abundance** or **HTseq data**. For the sample data included, this should be "abundance" for Transcript abundance and "treated" for HTseq_data (includes treated and untreated). If not specified, all files containing "abundance.tsv" will be used [Default: abundance.tsv]

Table that matches transcripts to genes.
- for **Transcript abundance** data. In CSV format, column 1 should be the transcript name and column 2 should be the Gene ID.  See the example below 

TXNAME,GENEID <br>
ENST00000456328.2,ENSG00000223972.5 <br>
ENST00000450305.2,ENSG00000223972.5 <br> 
ENST00000473358.1,ENSG00000243485.5 <br>

- If the file is not specified, The first file containing "gen" is used. [Default: ../data/Transcript_abundance/tx2gencode.v27.csv]

Type of **Transcript abundance** analysis
- What anlaysis was performed upstream e.g., kallisto, RSEM, Salmon, Sailfish. This only applies for **Transcript abundance** inputs. [Deafult: kallisto]

Path to Sample Metadata 
- Metadata file, required for **Transcript Abundance**, **Counts data**, **HTseq data**. The format of the metadata file is flexible, but it should be a CSV file that lists the condition corresponding to each sample. The design formula should use a column from this file to separate the control and test conditions. If the metadata file is not specified, the first file containing "Metadata" or "metadata" is used. [Default: ../data/Transcript_abundance/sample_metadata.csv]

Design formula
- A design formula tells the software which sources of variation to test for. This includes both the factor of interest as well as any additional covariates that are sources of variation. The formula should be a tilde (~) followed by the variables (columns of the metadata file) with plus signs between them. The design formula should have all of the factors in the metadata that account for major sources of variation in the data. DESeq2 can analyze any possible experimental design that can be expressed with fixed effects terms (multiple factors, designs with interactions, designs with continuous variables, splines, and so on are all possible). [Default: ~ condition]

- For the sample data included

    | Type of input data | Design | 
    | :---  | :--- | 
    | Counts_data | ~ individualID + tissueType | 
    | HTseq_data | ~ condition | 
    | Summarized_experiment | ~ cell + dex | 
    | Transcript_abundance | ~ condition |

Name of treatment/control column
- Column in the metadata file that indicates which treatment the sample was in, e.g., condition, tissueType, dex. [Default: condition]

- For the sample data included

    | Type of input data | Control Column | 
    | :---  | :--- | 
    | Counts_data | tissueType | 
    | HTseq_data | condition | 
    | Summarized_experiment | dex | 
    | Transcript_abundance | condition |

Minimum number of reads
- Minimum reads across all samples to include that sample in the data matrix [Default: 5]

Alpha value for significance
- the significance cutoff used for optimizing the independent filtering [Default: 0.1]. If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.

Number of genes to save individual plots for
- The # most significant genes to save plots for in the **/plots_by_gene** directory. [Default: 6]

Name of output file
- Name of CSV file without file extension. [Default: {type of input data}_DESeq2_results]


## Cite 

Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi:10.1186/s13059-014-0550-8.

<hr>

[Code Ocean](https://codeocean.com/) is a cloud-based computational platform that aims to make it easy for researchers to share, discover, and run code.<br /><br />
[![Code Ocean Logo](images/CO_logo_68x36.png)](https://www.codeocean.com)