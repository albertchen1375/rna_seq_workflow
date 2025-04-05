# ======================
# STEP 0: Load Libraries
# ======================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "pheatmap", "ggplot2", "RColorBrewer", "apeglm", "EnhancedVolcano"))

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(apeglm)
library(EnhancedVolcano)

# =========================
# STEP 1: Load Input Files
# =========================
# Replace with your file paths
count_file <- "counts.csv"       # Raw counts file (genes x samples)
metadata_file <- "metadata.csv"  # Metadata (samples x condition)

count_data <- read.csv(count_file, row.names = 1)
col_data <- read.csv(metadata_file, row.names = 1)

# Ensure sample names match
stopifnot(all(colnames(count_data) %in% rownames(col_data)))

# =========================
# STEP 2: Create DESeq2 Object
# =========================
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~ condition)  # Replace with your variable name

# ============================
# STEP 3: Pre-Filtering
# ============================
dds <- dds[rowSums(counts(dds)) > 10, ]  # Remove low-count genes

# ===================================
# STEP 4: Run DESeq2 Normalization
# ===================================
dds <- DESeq(dds)

# =======================================
# STEP 5: Transformation for Visualization
# =======================================
vsd <- vst(dds, blind = FALSE)

# PCA
pca_plot <- plotPCA(vsd, intgroup = "condition") + ggtitle("PCA of Samples")
print(pca_plot)

# Sample distance heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample Distance Heatmap")

# ==============================
# STEP 6: Differential Expression
# ==============================
res <- results(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")  # Use 'resultsNames(dds)' if unsure of coef

# Order by adjusted p-value
resOrdered <- res[order(res$padj), ]
write.csv(as.data.frame(resOrdered), file = "deseq2_results.csv")

# =========================
# STEP 7: MA Plot
# =========================
plotMA(res, ylim = c(-5, 5), main = "MA Plot")

# =========================
# STEP 8: Volcano Plot
# =========================
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "Volcano Plot")

# =========================
# STEP 9: Heatmap of Top Genes
# =========================
select <- order(res$padj)[1:20]  # Top 20 DE genes
mat <- assay(vsd)[select, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = col_data, main = "Top 20 DE Genes")

# =========================
# STEP 10: Save Normalized Counts
# =========================
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(as.data.frame(normalized_counts), file = "normalized_counts.csv")
