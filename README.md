# Bulk-RNA-seq-Differential-Expression-Analysis-using-DESeq2

## Project overview
This project follows a standard bulk RNA-seq differential expression analysis pipeline, aiming to identify genes that change expression between wildtype and silenced condition. And check if the differentially expressed genes are enriched in specific pathways, helping to interpret functional consequences of the condition difference.

## Project Aim
The main aim of this analysis is to:

- Identify **genes that are differentially expressed** between wildtype and silenced samples.
- Determine the **direction and magnitude of expression changes**
- Interpret these changes in a biological context by identifying **enriched pathways**

## General workflow
The analysis follows a standard RNA-seq pipeline:

Data loading  

⬇️

Preprocessing and filtering  

⬇️

Exploratory data analysis (PCA, sample distance) 

⬇️  

Differential expression analysis (DESeq2)  

⬇️  

Visualization (volcano plot)  

⬇️  

Pathway enrichment analysis (ORA, KEGG)

## Key steps

**1. Input data and experimental design**

The analysis starts with two key inputs: a count matrix and a metadata table. The count matrix contains raw gene expression values, where each row is a gene and each column is a sample. The metadata describes each sample, particularly which condition it belongs to. A crucial step here is ensuring that the order of samples in the count matrix matches the metadata, so that each sample is correctly labeled. The condition variable is then defined as a factor, with silenced set as the reference, which determines how fold changes will be interpreted later.

**2. Building the DESeq2 model**

The count data and metadata are combined into a DESeq2 object, which defines the statistical model used to test for differential expression. Before running the analysis, genes with very low counts are filtered out, because they provide little useful information and increase noise. At this stage, the dataset is prepared for both exploration and formal statistical testing.

**3. Exploratory data analysis**

Before testing for differential expression, the data is explored using variance stabilizing transformation. This transformation makes the data more suitable for visualization by stabilizing variance across genes. Using the transformed data, PCA and sample distance heatmaps are generated to assess whether samples cluster according to condition and to detect potential outliers. This step ensures that the dataset behaves as expected before proceeding to statistical testing.

**4. Differential expression analysis**

The core step is running DESeq2, which models the count data using a negative binomial distribution and tests whether gene expression differs between wildtype and silenced samples. The output includes, for each gene, a log2 fold change (indicating direction and magnitude of change) and an adjusted p-value (padj), which accounts for multiple testing. Genes are then classified as up-regulated, down-regulated, or not significant based on thresholds for both statistical significance and fold change.

**5. Visualization of results**

The results are visualized using a volcano plot, where the x-axis represents the log2 fold change and the y-axis represents statistical significance. This allows quick identification of genes that are both strongly changed and statistically significant. Genes on the right side are up-regulated in wildtype, while genes on the left are higher in silenced samples.

**6. Functional interpretation through ORA**

To move from gene lists to biological meaning, over-representation analysis is performed using KEGG pathways. This step tests whether the differentially expressed genes are enriched in specific biological pathways more than expected by chance. The results are summarized in a bubble plot, where larger points indicate a higher proportion of pathway genes represented (RichFactor), and color reflects statistical significance.

## Results

The analysis identified a set of genes significantly differentially expressed between wildtype and silenced samples, including both up-regulated and down-regulated genes.

Exploratory analysis showed that samples cluster according to condition, indicating a strong transcriptional effect.

Pathway enrichment analysis revealed that differentially expressed genes are associated with specific biological pathways, providing insight into the functional impact of the silencing condition.



