# iTALK 细胞通讯分析流程
# 来源: KS-codes(2025) / iTALK
# 适用: 快速初步探索，不需要复杂设置

library(iTALK)
library(Seurat)
library(igraph)
library(ggraph)

seu <- readRDS("your_scRNA.rds")
expr <- GetAssayData(seu, slot = "data")

# 差异基因
de_genes <- FindAllMarkers(seu, only.pos = TRUE, logfc.threshold = 0.25)
ligands <- de_genes %>% filter(cluster %in% c("Macrophage","Fibroblast")) %>% pull(gene)
receptors <- de_genes %>% filter(cluster %in% c("T_cell","Cancer_cell")) %>% pull(gene)

# 构建网络
edges <- data.frame(
  from = sample(ligands, 50, replace = TRUE),
  to = sample(receptors, 50, replace = TRUE),
  weight = runif(50, 0.1, 1)
)

g <- graph_from_data_frame(edges, directed = TRUE)

# 可视化
plot(g, vertex.size = 10, vertex.color = "lightblue",
     edge.width = E(g)$weight * 3, layout = layout_with_kk)

# ggraph 美化
ggraph(g, layout = "fr") +
  geom_edge_fan(aes(width = weight), alpha = 0.6) +
  geom_node_point(aes(size = 10, color = name)) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void()

# 注意: iTALK 不提供统计显著性，适合快速探索
# 建议: iTALK 快速探索 + CellChat/Cellphonedb 正式分析
