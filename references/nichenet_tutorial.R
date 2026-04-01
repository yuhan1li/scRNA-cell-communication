# Nichenet 单样本/组间分析流程
# 来源: KS-codes(2025) / Nichenet V2（视频教程）
# 适用: 配体→受体→信号→靶基因完整调控链条分析

library(nichenetr)
library(Seurat)
library(tidyverse)

seu <- readRDS("nichenet_seurat.rds")
# 元数据需包含: sample_id, group_id, celltype_id

load("lr_network_human_allInfo.rds")
load("ligand_target_matrix.rds")

sender_celltypes <- c("Fibroblast", "Macrophage", "DC")
receiver_celltypes <- c("T_cell", "Cancer_cell")

# 感兴趣的靶基因（差异表达基因）
interesting_genes <- rownames(FindMarkers(seu, ident.1 = "treated", ident.2 = "control",
  group.by = "group_id", subset.ident = "Cancer_cell", assay = "RNA"))
interesting_genes <- interesting_genes[interesting_genes %in% rownames(ligand_target_matrix)]

# 配体活性分析
ligand_activities <- predict_ligand_activities(
  seurat_obj = seu, geneset = interesting_genes,
  ligands = unique(lr_network_human_allInfo$from),
  lr_network = lr_network_human_allInfo,
  ligand_target_matrix = ligand_target_matrix)

top_ligands <- ligand_activities %>% arrange(-pearson) %>% pull(ligand) %>% head(10)

# 受体分析
receptors_expressed <- find_receptors_expressed(
  seurat_obj = seu, receivers = receiver_celltypes,
  lr_network = lr_network_human_allInfo)

# Circos 图
make_circos_group_comparison(seu, ligand_activities, ligand_target_matrix,
  group1 = "treated", group2 = "control")

# 配体-受体-靶基因综合图
make_lr_prod_activity_plots(seu, ligand_activities, top_ligands = 10)

# 保存
write.csv(ligand_activities, "nichenet_ligand_activities.csv")
saveRDS(ligand_activities, "nichenet_ligand_activities.rds")
