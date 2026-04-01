# Cellphonedb V5 数据准备代码
# 来源: KS-codes(2025) / 视频教程7：单细胞通讯之Cellphonedb V5（2024更新版）

library(Seurat)

# ============================================================================
# 数据准备
# ============================================================================
# 加载 Seurat 对象
scRNA <- readRDS("your_annotated_scRNA.rds")
DefaultAssay(scRNA) <- "RNA"

# ============================================================================
# 1. 提取归一化表达矩阵
# ============================================================================
Normalized_counts <- GetAssayData(scRNA, layer = 'data')
Normalized_counts <- as.data.frame(Normalized_counts)

# 基因名列（第一列）
Normalized_counts <- cbind(rownames(Normalized_counts), Normalized_counts)
colnames(Normalized_counts)[1] <- 'Gene'

# ============================================================================
# 2. 准备元数据
# ============================================================================
metadata <- data.frame(
  Cell = rownames(scRNA@meta.data),
  cell_type = scRNA$celltype
)

# 重要：空格替换为下划线
metadata$cell_type <- gsub(' ', '_', metadata$cell_type)

# ============================================================================
# 3. 保存文件（制表符分隔，无引号）
# ============================================================================
write.table(Normalized_counts, "Normalized_counts.txt",
            row.names = F, sep = '\t', quote = F)

write.table(metadata, "cellphonedb_meta.txt",
            row.names = F, sep = '\t', quote = F)

# ============================================================================
# 4. Python 端运行命令
# ============================================================================
# 有生物学重复：
# cellphonedb statics Normalized_counts.txt cellphonedb_meta.txt \
#   --counts-data hgnc_symbol \
#   --output-path output_v5 \
#   --threads 8 \
#   --iterations 1000 \
#   --pvalues-method permtest

# 无生物学重复：
# cellphonedb statics-single Normalized_counts.txt cellphonedb_meta.txt \
#   --counts-data hgnc_symbol \
#   --output-path output_single \
#   --threads 8
