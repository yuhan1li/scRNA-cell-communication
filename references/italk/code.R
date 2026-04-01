# iTALK 细胞通讯分析代码
# 来源: KS-codes(2025) / iTALK
# 注意: iTALK 不提供统计显著性检验，仅适合快速初步探索

library(iTALK)
library(Seurat)
library(igraph)

seu <- readRDS("your_scRNA.rds")

# ============================================================================
# 1. 差异基因分析（找各细胞类型的 marker）
# ============================================================================
de_genes <- FindAllMarkers(seu, only.pos = TRUE, logfc.threshold = 0.25)

# ============================================================================
# 2. 提取配体和受体基因
# ============================================================================
# 从差异基因中筛选 sender（配体）和 receiver（受体）细胞类型的 marker
ligands <- de_genes %>%
  filter(cluster %in% c("Macrophage", "Fibroblast")) %>%
  pull(gene)

receptors <- de_genes %>%
  filter(cluster %in% c("T_cell", "Cancer_cell")) %>%
  pull(gene)

# ============================================================================
# 3. 构建配体-受体网络数据框
# ============================================================================
# 随机采样生成互作对（演示用）
set.seed(123)
edges <- data.frame(
  from = sample(ligands, 50, replace = TRUE),
  to = sample(receptors, 50, replace = TRUE),
  weight = runif(50, 0.1, 1)
)

# ============================================================================
# 4. 可视化网络图
# ============================================================================
g <- graph_from_data_frame(edges, directed = TRUE)

# 基础网络图
plot(g, vertex.size = 10, vertex.color = "lightblue",
     edge.width = E(g)$weight * 3, layout = layout_with_kk)

# ============================================================================
# 5. ggraph 美化版本
# ============================================================================
library(ggraph)

ggraph(g, layout = "fr") +
  geom_edge_fan(aes(width = weight), alpha = 0.6) +
  geom_node_point(aes(size = 10, color = name)) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void()

# ============================================================================
# 重要提示:
# iTALK 的网络是随机生成的，不代表真实通讯强度
# 建议使用 CellChat 或 Cellphonedb 进行正式分析
# ============================================================================
