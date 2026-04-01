# MultiNichenet 无重复样本多组受配体差异分析
# 来源: KS-codes(2025) / 1-scRNA 细胞通讯分析 / multinichentr（视频教程）
# 适用: 无样本重复或样本较少的临床数据分析

library(MultiNichenet)
library(Seurat)
library(tidyverse)

# 网络文件（百度云: hqp3）
load("lr_network_human_allInfo.rds")
load("ligand_target_matrix.rds")

seu <- readRDS("scRNA_nichenet.rds")
# 元数据需包含: sample_id, group_id, celltype_id

# ============================================================================
# 三、定义 contrasts（比较组）
# ============================================================================
contrasts <- list(contrast_1 = c("SD", "HC"))

# ============================================================================
# 四、定义 sender / receiver
# ============================================================================
sender_celltypes <- c("Tcell", "Fibroblast", "Macrophage", "Dendritic")
receiver_celltypes <- c("Tcell", "Fibroblast", "Macrophage", "Cancer_cell")

# ============================================================================
# 五、细胞类型过滤
# ============================================================================
abundance_data <- get_abundance_info(
  seurat_obj = seu, sample_id = "sample_id",
  celltype_id = "celltype_id", min_cells = 10)

# ============================================================================
# 六、基因过滤（无重复方法）
# ============================================================================
gene_filtering <- get_frac_exprs_sampleAgnostic(
  abundance_data = abundance_data, fraction_cutoff = 0.05)

# ============================================================================
# 七、处理表达数据
# ============================================================================
expression_data <- process_abundance_expression_info(
  seurat_obj = seu, abundance_data = abundance_data,
  gene_filtering = gene_filtering, sample_id = "sample_id")

# ============================================================================
# 八、差异分析（无重复）
# ============================================================================
de_info <- get_DE_info_sampleAgnostic(
  seurat_obj = seu, expression_data = expression_data,
  sample_id = "sample_id", group_id = "group_id",
  celltype_id = "celltype_id", logfc_threshold = 0.25, padj_cutoff = 0.05)

# ============================================================================
# 九、配体活性预测
# ============================================================================
ligand_activities <- get_ligand_activities_targets_DEgenes(
  expression_data = expression_data, de_info = de_info,
  contrasts = contrasts, ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network_human_allInfo)

# ============================================================================
# 十、生成优先级排序表
# ============================================================================
prioritization <- generate_prioritization_tables(
  ligand_activities = ligand_activities, de_info = de_info,
  expression_data = expression_data, abundance_data = abundance_data,
  contrasts = contrasts, celltype_id = "celltype_id",
  sender_celltypes = sender_celltypes, receiver_celltypes = receiver_celltypes)

# ============================================================================
# 十一、可视化
# ============================================================================
# Circos 比较图
make_circos_group_comparison(
  seurat_obj = seu, prioritization = prioritization,
  contrast = "contrast_1", celltype_id = "celltype_id", group_id = "group_id")

# 配体-受体-活性-靶基因综合图
make_sample_lr_prod_activity_plots(
  seurat_obj = seu, prioritization = prioritization,
  contrast = "contrast_1", top_n = 10)

# ============================================================================
# 十二、保存结果
# ============================================================================
saveRDS(prioritization, "multinichenet_prioritization.rds")
saveRDS(ligand_activities, "multinichenet_ligand_activities.rds")
write.csv(prioritization$contrast_1, "multinichenet_prioritization_contrast1.csv", row.names = FALSE)
