# CellChat V2 详细版更新教程
# 来源: KS-codes(2025) / 视频教程6：单细胞通讯之Cellchat V2（2024更新版）

library(CellChat)
library(Seurat)
library(ggplot2)
library(ggalluvial)
library(patchwork)

# ============================================================================
# 1. 数据准备
# ============================================================================
seu <- readRDS("D:/.../your_scRNA_seurat.rds")
data.input <- GetAssayData(seu, assay = "RNA", slot = "data")
labels <- Idents(seu)

# ============================================================================
# 2. 创建 CellChat 对象
# ============================================================================
cellchat <- createCellChat(object = data.input, group.by = "ident", database = "human")
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

# ============================================================================
# 3. 数据预处理
# ============================================================================
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)

# ============================================================================
# 4. 计算通讯概率
# ============================================================================
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.25, population.size = FALSE)
cellchat <- computeCommunProbPathways(cellchat)

# ============================================================================
# 5. 聚合网络分析
# ============================================================================
cellchat <- aggregateNet(cellchat)
head(cellchat@net$count)
head(cellchat@net$weight)

# ============================================================================
# 6. 可视化
# ============================================================================
# 热图
png("Fig1_CellChat_Heatmap.png", width = 2500, height = 2000, res = 300)
netVisual_heatmap(cellchat)
dev.off()

# 气泡图
netVisual_bubble(cellchat, sources.use = c(1,2), targets.use = c(3,4), remove.isolate = FALSE)

# 和弦图
pathways.show <- c("TGFb")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# 散点图（两组比较）
cellchat <- netAnalysis_rankNet(cellchat, slot.name = "netP", color.use = c("#E41A1C","#377EB8"))
netAnalysis_signalingRole_scatter(cellchat)

# 单通路详细网络
netVisual_individual(cellchat, signaling = pathways.show, plot.target = FALSE)

# 细胞类型对某通路的贡献热图
netAnalysis_signalingRole_heatmap(cellchat, signaling = pathways.show)

# ============================================================================
# 7. 调控网络分析
# ============================================================================
cellchat <- netAnalysis_run(cellchat, assays = "logfc", features.name = "DEGs")
netAnalysis_centrality_scatter(cellchat, signaling = pathways.show, font.size = 10)

# ============================================================================
# 8. 输出结果
# ============================================================================
write.csv(cellchat@net$prob, "CellChat_communication_probability.csv")
write.csv(cellchat@net$pval, "CellChat_pvalue.csv")
LRpair <- subset(cellchat@LR, prob > 0.01)
write.csv(LRpair, "CellChat_significant_LR_pairs.csv")
saveRDS(cellchat, "cellchat_analysis.rds")

# ============================================================================
# 9. 高级可视化
# ============================================================================
cellchat@meta$colors.use <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00")
pdf("Fig_CellChat_Chord.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()
