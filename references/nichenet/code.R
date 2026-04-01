# Nichenet 调控网络分析代码
# 来源: KS-codes(2025) / Nichenet V2（视频教程）

library(nichenetr)
library(Seurat)

seu <- readRDS("nichenet_seurat.rds")
load("lr_network_human_allInfo.rds")
load("ligand_target_matrix.rds")

# ============================================================================
# 1. 定义 sender 和 receiver 细胞类型
# ============================================================================
sender_celltypes <- c("Fibroblast", "Macrophage", "DC")
receiver_celltypes <- c("T_cell", "Cancer_cell")

# ============================================================================
# 2. 差异基因分析（Receiver 细胞中）
# ============================================================================
# 感兴趣的靶基因：疾病组 vs 对照组在特定细胞类型中的差异基因
interesting_genes <- rownames(
  FindMarkers(seu, ident.1 = "treated", ident.2 = "control",
              group.by = "group_id", subset.ident = "Cancer_cell", assay = "RNA")
)

# 只保留在配体-靶基因矩阵中的基因
interesting_genes <- interesting_genes[interesting_genes %in% rownames(ligand_target_matrix)]

# ============================================================================
# 3. 配体活性预测
# ============================================================================
ligand_activities <- predict_ligand_activities(
  seurat_obj = seu,
  geneset = interesting_genes,
  ligands = unique(lr_network_human_allInfo$from),
  lr_network = lr_network_human_allInfo,
  ligand_target_matrix = ligand_target_matrix)

# 按活性排序，取 top 配体
top_ligands <- ligand_activities %>% arrange(-pearson) %>% pull(ligand) %>% head(10)

# ============================================================================
# 4. 受体分析（receiver 细胞表达的受体）
# ============================================================================
receptors_expressed <- find_receptors_expressed(
  seurat_obj = seu, receivers = receiver_celltypes,
  lr_network = lr_network_human_allInfo)

# ============================================================================
# 5. Circos 可视化
# ============================================================================
make_circos_group_comparison(
  seurat_obj = seu,
  ligand_activities = ligand_activities,
  ligand_target_matrix = ligand_target_matrix,
  group1 = "treated", group2 = "control"
)

# ============================================================================
# 6. 配体-受体-活性-靶基因综合图
# ============================================================================
make_lr_prod_activity_plots(seu, ligand_activities, top_ligands = 10)

# ============================================================================
# 7. 保存结果
# ============================================================================
write.csv(ligand_activities, "nichenet_ligand_activities.csv", row.names = FALSE)
saveRDS(ligand_activities, "nichenet_ligand_activities.rds")
