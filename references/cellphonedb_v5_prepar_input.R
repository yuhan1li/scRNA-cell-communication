# Cellphonedb V5 数据准备与运行
# 来源: KS-codes(2025) / 视频教程7：单细胞通讯之Cellphonedb V5

library(Seurat)
library(tidyverse)

# ============================================================================
# 1. 数据准备 - Seurat → V5 格式
# ============================================================================
# 注意: V5 格式是 Cell x Gene（与 V3 的 Gene x Cell 相反！）

seu <- readRDS("D:/.../your_scRNA_seurat.rds")

# 表达矩阵 (Cell x Gene)
counts_matrix <- GetAssayData(seu, slot = "counts")
counts_df <- as.data.frame(t(as.matrix(counts_matrix)))
counts_df$cell_id <- rownames(counts_df)
write.csv(counts_df, "cpdb_counts.csv", row.names = FALSE)

# 元数据
meta_df <- data.frame(cell_id = colnames(seu), cell_type = as.character(Idents(seu)))
write.csv(meta_df, "cpdb_meta.csv", row.names = FALSE)

# ============================================================================
# 2. Python 端运行
# ============================================================================
# pip install cellphonedb
# cellphonedb statics cpdb_counts.csv cpdb_meta.csv \
#   --counts-data gene_name --output-path output_v5 --threads 8 \
#   --iterations 1000 --pvalues-method permtest

# ============================================================================
# 3. 结果可视化
# ============================================================================
library(pheatmap)

pvals <- read.csv("output_v5/statistical_analysis_pvalues_*.csv", row.names = 1, check.names = FALSE)
means <- read.csv("output_v5/statistical_analysis_significant_means_*.csv", row.names = 1, check.names = FALSE)

sig_means <- means[apply(pvals < 0.05, 1, sum) > 0, ]

pheatmap(sig_means,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE,
         main = "Cellphonedb V5: Significant Interactions")
