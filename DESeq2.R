# ------------------------
## 1. Load RNA-seq dataset:
# ------------------------

#Load expression object
FOXA1 <- readRDS("FOXA1_PC3_exprs.rds")

# Inspect object structure
str(FOXA1)

#Remove Ensembl version suffix from gene IDs
rownames(FOXA1) = sub("\\..*", "", rownames(FOXA1))

# Inspect sample names
colnames(FOXA1)


# -------------------------
## 2. Extract count matrix:
# -------------------------

#Select RNA-seq count columns
counts = FOXA1[, 5:12]

#Inspect counts
counts


# ------------------------
## 3. Build metadata table:
# -------------------------

#Create metadata using sample names
metadata = data.frame(
  barcode = colnames(counts)
)

#Assign experimental condition based on sample name, "psi" = silenced, otherwise = wildtype
metadata$condition = ifelse(grepl("psi", metadata$barcode),"silenced", "wildtype")


# -------------------------
## 4. Prepare DESeq2 input:
# -------------------------

library(DESeq2)

#Convert counts into integer matrix (DESeq2 requires raw integer counts)
counts = floor(as.matrix(counts))

#Convert condition to factor, reference = silenced, comparison = wildtype
metadata$condition = factor(metadata$condition,
                            levels = c("silenced", "wildtype"))

#Create DESeq2 object
dds = DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)


# ---------------------------------
## 5. Filter lowly expressed genes:
#----------------------------------

#Dimensions before filtering
dim(dds)

#Remove genes with extremely low counts (Helps reduce noise and improve statistical power)
dds = dds[rowSums(counts(dds)) > 1, ]

#Dimensions after filtering
dim(dds)


# ---------------------------------------------
## 6.Variance stabilizing transformation (VST):
# ---------------------------------------------

#Transform counts for visualization and PCA (VST is used for exploratory analysis,not for the DESeq2 statistical test itself)
vsd = vst(dds, blind = FALSE)


# ---------------------------------------
## 7. Principal Component Analysis (PCA):
#----------------------------------------

library(ggplot2)

#Extract PCA coordinates
pca.data = plotPCA(vsd, returnData = TRUE)

#Percentage of variance explained by PCs
percentVar = round(100 * attr(pca.data, "percentVar"))

#PCA plot (Samples clustering together are transcriptionally similar)
p_pca <- ggplot(pca.data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_bw()

ggsave("pca.png", plot = p_pca, width = 6, height = 5)

p_pca


# --------------------------------------
## 8. Sample-to-sample distance heatmap:
# --------------------------------------

library(ComplexHeatmap)
library(RColorBrewer)

#Compute pairwise sample distances
sampleDists = dist(t(assay(vsd)))

#Convert distance object into matrix
mat = as.matrix(sampleDists)

#Define heatmap color palette
col_fun = colorRampPalette(brewer.pal(9, "Blues"))(255)

png("distance_heatmap.png", width = 800, height = 800)

Heatmap(mat, name = "Distance", col = col_fun)

dev.off()


# ------------------------------------------------
## 9. Run DESeq2 differential expression analysis:
# ------------------------------------------------

#Run DESeq2 pipeline (normalization/ dispersion estimation/ model fitting/ statistical testing)
dds = DESeq(dds)

#Extract results(wildtype vs silenced)
res = results(dds, name = "condition_wildtype_vs_silenced")

#Sort genes by adjusted p-value
res = res[order(res$padj), ]

# --------------------------
## 10. Result distributions:
# --------------------------


#Distribution of average normalized expression
ggplot(res, aes(x = log10(baseMean + 1))) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  labs(title = "baseMean distribution", x = "log10(baseMean + 1)") +
  theme_bw()

#Distribution of fold changes
ggplot(res, aes(x = log2FoldChange)) +
  geom_density(fill = "lightgreen", alpha = 0.5, na.rm = TRUE) +
  labs(title = "log2FC distribution", x = "log2FoldChange") +
  theme_bw()

#Define DE thresholds:
#Statistical significance threshold
padj.th = 0.1
#Fold change threshold
l2fc.th = 0.5

#Initialize gene classification
res$status = "none"

#Up-regulated in wildtype
res$status[res$padj <= padj.th & res$log2FoldChange >= l2fc.th] = "up-regulated"

#Down-regulated in wildtype (higher in silenced)
res$status[res$padj <= padj.th & res$log2FoldChange <= -l2fc.th] = "down-regulated"

#Count genes in each category
table(res$status)


# --------------------
## 11. Volcano plot:
# --------------------

#Volcano plot(x-axis = fold change, y-axis = statistical significance)
p_volcano <- ggplot(res, aes(log2FoldChange, -log10(padj), color = status)) +
  geom_point(alpha = 0.6) +
  #Fold-change thresholds
  geom_vline(xintercept = c(-l2fc.th, l2fc.th), linetype = "dashed") +
  #Adjusted p-value threshold
  geom_hline(yintercept = -log10(padj.th), linetype = "dashed") +
  scale_color_manual(values = c("grey", "red", "blue")) +
  theme_bw()

ggsave("volcano.png", plot = p_volcano, width = 6, height = 5)

p_volcano

#More stringent threshold
padj.th2 = 0.05
l2fc.th2 = 1

#Count highly significant DE genes
sum(res$padj <= padj.th2 & abs(res$log2FoldChange) >= l2fc.th2, na.rm = TRUE)


# -----------------------------------
## 12.  Over-representation analysis:
# -----------------------------------

library(msigdbr)
library(clusterProfiler)

#Add gene symbols
res$gene_name = FOXA1$gene_name[match(rownames(res), rownames(FOXA1))]

#Load Hallmark gene sets
hallmark = msigdbr(species = "Homo sapiens", category = "H")
hallmark = hallmark[, c("gs_name", "gene_symbol")]

#Separate up/down-regulated genes
up.genes = na.omit(res$gene_name[res$status == "up-regulated"])
down.genes = na.omit(res$gene_name[res$status == "down-regulated"])


# --------------------------------
## 13. Run enrichment analysis:
# --------------------------------

#ORA for up-regulated genes
ora_up = enricher(up.genes, TERM2GENE = hallmark)

#ORA for down-regulated genes
ora_down = enricher(down.genes, TERM2GENE = hallmark)

#Convert results into data frames
ora_up = as.data.frame(ora_up)
ora_down = as.data.frame(ora_down)

#Add category labels
ora_up$status = "up"
ora_down$status = "down"

#Merge results
ora = rbind(ora_up, ora_down)


# -------------------------
## 14. Compute RichFactor:
# -------------------------

# RichFactor = proportion of pathway genes represented in DE list
ora$RichFactor = sapply(strsplit(ora$GeneRatio, "/"),
                        function(x) as.numeric(x[1]) / as.numeric(x[2]))

#Keep top 15 pathways
ora = ora[order(ora$p.adjust), ]
ora_top = ora[1:15, ]


# ----------------------
## 15. ORA bubble plot:
# ----------------------

#Bubble plot(size = RichFactor, color = adjusted p-value, shape = up/down category)
p_ora <- ggplot(ora_top, aes(RichFactor, reorder(Description, RichFactor))) +
  geom_point(aes(size = RichFactor,
                 fill = p.adjust,
                 shape = status),
             color = "black") +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_viridis_c() +
  theme_bw()

ggsave("ora.png", plot = p_ora, width = 7, height = 5)

p_ora

