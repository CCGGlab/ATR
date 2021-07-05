This notebook contains high level information on the bioinformatics pipeline that was used for the analysis reported in Szydzik J et al. ATR inhibition enables complete tumour regression in ALK-driven NB mouse models, Nature Communications 2021

# Environment
  
Analysis was performed in a Conda environment. See **ATR.yml** for details. **scripts/Rpacks** describes R packages that were installed independently.

# Experimental conditions

RNA-Seq data cell lines after ATR inhibition with 50nM BAY 1895344 

- CLB-BAR: 0h (ctrl), 24h, 48h
- CLB-GE: 0h (ctrl), 24h, 48h

Proteomics data cell lines after ATR inhibition with 50nM BAY 1895344 

- CLB-GE: 0h (ctrl), 24h, 48h (Exp 3336 & 3337)
- CLB-BAR: 0h (ctrl), 24h, 48h (Exp 3403)
- CLB-BAR: 0h (ctrl), 1h, 3h, 6h, 24h (Exp 4060)
- CLB-BAR: SYNCHRONIZED, 6h (DMSO), 6h (Treatment) (Exp 4082)

Phosphoproteomics data cell lines after ATR inhibition with 50nM BAY 1895344 

- CLB-BAR: SYNCHRONIZED, 6h (DMSO), 6h (Treatment) (Exp 4082)

RNA-Seq data ALK and ALKAL mice after 3d ATR inhibition with 50 nM BAY 1895344. 

# Data

- Raw RNA-Seq data have been deposited in Arrayexpress 
  - Cell lines: E-MTAB-10603
  - Mice: E-MTAB-10616 (BAY treatment) & E-MTAB-9600 (Control samples)
- Processed data are available in *data/ATR_data.RData*. This file containes the following objects:
  - normalized_counts: DESeq2-normalized counts from cell line RNA-Seq data
  - res_diff_expr: DESeq2 output from cell line RNA-Seq data
  - normalized_counts_prot: Normalized counts from cell line proteomics data
  - res_diff_expr_prot: ROTS output from cell line proteomics data
  - normalized_counts_PP: Normalized counts from cell line phosphoproteomics data
  - res_diff_expr_PP: ROTS output from cell line phosphoproteomics data
  - sample_info: sample information cell line RNA-Seq data
  - sample_info_prot: sample information cell line (phospho-)proteomics data
  - normalized_counts_mice: DESeq2-normalized counts from mice RNA-Seq data
  - res_diff_expr_mice: DESeq2 output from mice RNA-Seq data
  - sample_info_mice: sample information mice RNA-Seq data
  - SQ_sites: Information on S/TQ sites of phosphoproteomics data
  - geneset_ls: genesets downloaded from MSigDB v7.2
  - MGI_to_HGNC: MGI to HGNC mapping table
  - ppi_db: High confidence interactions from STRING

# Manuscript

The main analysis, as reported in the manuscript

## RNA-Seq cell lines
```{r}
source("scripts/manuscript_RNA.R")
```

## Proteomics cell lines
```{r}
source("scripts/manuscript_prot.R")
```

## PP cell lines

Basic analysis: volcano plot, ST/Q sites, GSEA (+ table)
```{r}
source("scripts/manuscript_PP.R")
```

Heatmap with main results
```{r}
source("scripts/manuscript_PP_hm.R")
```

Network analysis
```{r}
source("scripts/manuscript_PP_network.R")
```

Motif analysis
```{r}
source("scripts/manuscript_PP_motifs.R")
```

## RNA-Seq mice
```{r}
source("scripts/manuscript_mice.R")
```

## Summary table
```{r}
source("scripts/manuscript_summary.R")
```

