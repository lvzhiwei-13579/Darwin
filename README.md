Project Name: Molecular Characteristics and Core Gene Analysis of High-Grade Gliomas (WHO III vs IV)

Project Goal: Through bioinformatics analysis, screen differentially expressed genes, enriched pathways, and core regulatory genes between WHO grade III and IV gliomas, and elucidate their molecular mechanisms.
I. Environment Preparation
Operating Environment: R 4.2.0 or higher
Installation of Required R Packages:
install.packages(c("tidyverse", "limma", "pheatmap", "ggplot2", "dplyr", "igraph", "readr"))# Bioconductor packagesif (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "biomaRt"))
II. Analysis Workflow and Script Description
1. Data Preprocessing (data_preprocessing.R)
Function: Perform cleaning, standardization, and sample matching on gene expression data and clinical data.
Input Files:
Gene expression data: GENE.txt (rows: genes, columns: samples, expression matrix)
Clinical data: clinical.txt (including sample ID, WHO grade, survival information, IDH mutation status, etc.
Core Steps:
Read data and check format;
Filter low-expression genes (mean ≥ 1) and perform log2 normalization on high-expression data;
Extract key clinical features (survival time, survival status, WHO grade, IDH status) and fill missing values;
Match common samples between gene data and clinical data, and save preprocessing results.
Output Files:
Preprocessed data: preprocessed_data/preprocessed_data.RData
Matched gene/clinical data: preprocessed_data/gene_matched.csv, preprocessed_data/clinical_matched.csv
2. Differentially Expressed Genes (DEGs) Analysis (DEGs_analysis.R)
Function: Screen significantly differentially expressed genes between WHO grade III and IV gliomas and visualize the results.
Input File: Preprocessed gene and clinical data (preprocessed_data.RData)
Core Steps:
Screen WHO grade III and IV samples and construct a grouping matrix;
Perform differential analysis using the limma package, with thresholds: |log2 (fold change)| > 1 and adjusted P-value < 0.05;
Visualization: volcano plot (showing distribution of differential genes), heatmap (showing expression patterns of top 50 differential genes, with separate plots for up-regulated/down-regulated genes).
3. Functional Enrichment Analysis (enrichment_analysis.R)
Function: Perform GO (Gene Ontology) and KEGG (Kyoto Encyclopedia of Genes and Genomes) enrichment analyses on differential genes to interpret biological functions.
Input File: List of differential genes (DEGs_results/filtered_differential_gene_list.csv)
Core Steps:
Gene ID conversion: Convert gene symbols to Entrez IDs;
GO enrichment: Analyze biological processes (BP), cellular components (CC), and molecular functions (MF) involved in differential genes;
KEGG enrichment: Analyze signaling pathways enriched in differential genes (supporting the use of offline cache when online database connection fails);
Visualization: GO enrichment dot plot, KEGG enrichment bar plot.
Output Files:
Enrichment results: Enrichment_results/GO_enrichment_results.csv, Enrichment_results/KEGG_enrichment_results.csv
Visualization results: Enrichment_results/GO_enrichment_plot.png, Enrichment_results/KEGG_enrichment_plot.png
4. Core Gene Screening (core_genes_analysis.R)
Function: Screen core regulatory genes based on protein-protein interaction (PPI) network.
Input Files:
KEGG enrichment results (Enrichment_results/KEGG_enrichment_results.csv)
List of differential genes (DEGs_results/filtered_differential_gene_list.csv)
STRING database protein interaction data: string_interactions.txt (to be downloaded from STRING in advance, including protein interaction relationships and confidence scores)
Core Steps:
Convert protein IDs in STRING data to gene symbols;
Construct PPI network based on genes from significantly enriched KEGG pathways (confidence threshold > 400);
Calculate the degree of genes in the network, and screen genes with degree values in the top 20% as core genes;
Visualization: PPI network with core genes labeled.
III. Running Instructions
Run the scripts in order: data_preprocessing.R → DEGs_analysis.R → enrichment_analysis.R → core_genes_analysis.R;
All input file paths need to be modified according to their actual storage locations;
If KEGG online connection fails, the script will automatically download offline cache data (approximately 100MB), so ensure network connectivity.
