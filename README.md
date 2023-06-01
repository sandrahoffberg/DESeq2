# DESeq2

### How this capsule works: 

This is the fifth capsule in the pipeline. Running this capsule is optional. <br>
DESeq2 estimates variance-mean dependence in count data from high-throughput sequencing assays and tests for differential expression based on a model using the negative binomial distribution.

### Inputs: 

- List of control genes (.csv)
- List of target genes (.csv)
- {sample}.raw.exp.csv from Seurat output 

### Outputs: 

- Many csv files named with {sample}_{condition}

### App Panel Parameters: 
- Run this capsule? (yes/no)
- Raw expression file from Seurat output (.raw.exp.csv)
- The gene target list (.csv)
- The metadata file (.csv)
- Control condition (e.g., DMSO_0_144)