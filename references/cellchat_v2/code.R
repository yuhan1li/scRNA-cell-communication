# CellChat V2 完整分析代码
# 来源: KS-codes(2025) / 视频教程6：单细胞通讯之Cellchat V2（2024更新版）

library(CellChat)
library(Seurat)

# ============================================================================
# 1. 数据准备
# ============================================================================
# 加载 Seurat 对象
sceV5 <- load("sceV5.RData")  # 示例数据

# 分组提取（示例：HD vs MDA5 两组）
HD <- subset(sceV5, group1 == 'HD')
MDA <- subset(sceV5, group1 == 'MDA5')

DefaultAssay(HD) <- 'SCT'

# 提取表达矩阵和元数据
HD_input <- GetAssayData(HD, layer = 'data')
HD_meta <- HD@meta.data[, c("group1", "cell_type")]
colnames(HD_meta) <- c("group", "labels")

# 检查一致性
identical(colnames(HD_input), rownames(HD_meta))  # 必须返回 TRUE

# ============================================================================
# 2. 创建 CellChat 对象
# ============================================================================
HD.cellchat <- createCellChat(object = HD_input, meta = HD_meta, group.by = "labels")
levels(HD.cellchat@idents)
groupSize <- as.numeric(table(HD.cellchat@idents))

# ============================================================================
# 3. 设置数据库
# ============================================================================
# 人源数据
CellChatDB <- CellChatDB.human
# 鼠源数据用 CellChatDB.mouse
HD.cellchat@DB <- CellChatDB

# ============================================================================
# 4. 预处理
# ============================================================================
HD.cellchat <- subsetData(HD.cellchat)
future::plan("multisession", workers = 2)  # 并行计算
HD.cellchat <- identifyOverExpressedGenes(HD.cellchat)
HD.cellchat <- identifyOverExpressedInteractions(HD.cellchat)

# ============================================================================
# 5. 计算通讯概率
# ============================================================================
HD.cellchat <- computeCommunProb(HD.cellchat, type = "triMean")
HD.cellchat <- filterCommunication(HD.cellchat, min.cells = 10)
HD.cellchat <- computeCommunProbPathway(HD.cellchat)
HD.cellchat <- aggregateNet(HD.cellchat)

# ============================================================================
# 6. 可视化
# ============================================================================
# 6.1 互作数目 Circle 图
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(HD.cellchat@net$count, vertex.weight = groupSize,
                weight.scale = T, label.edge = F,
                title.name = "Number of interactions")
netVisual_circle(HD.cellchat@net$weight, vertex.weight = groupSize,
                weight.scale = T, label.edge = F,
                title.name = "Interaction weights/strength")

# 6.2 热图
pheatmap::pheatmap(HD.cellchat@net$count, cluster_cols = F, cluster_rows = F,
                  fontsize = 10, display_numbers = T, number_color = "black")

# 6.3 指定通路可视化
pathways.show <- c("MHC-II")
netVisual_aggregate(HD.cellchat, signaling = pathways.show)  # Circle
netVisual_aggregate(HD.cellchat, signaling = pathways.show, layout = "chord")  # 和弦图
netVisual_heatmap(HD.cellchat, signaling = pathways.show, color.heatmap = "Reds")  # 热图

# 6.4 分组和弦图
group.cellType <- c(rep("immune", 5), rep("other", 4))
names(group.cellType) <- levels(HD.cellchat@idents)
netVisual_chord_cell(HD.cellchat, signaling = pathways.show, group = group.cellType)

# 6.5 气泡图
netVisual_bubble(HD.cellchat, sources.use = "Tcell", remove.isolate = FALSE)

# 6.6 单个配体-受体对可视化
pairLR <- extractEnrichedLR(HD.cellchat, signaling = pathways.show, geneLR.return = FALSE)
netVisual_individual(HD.cellchat, signaling = pathways.show,
                    pairLR.use = "HLA-DPA1_CD4", layout = "circle")

# ============================================================================
# 7. 系统分析
# ============================================================================
# 7.1 网络中心性分析
HD.cellchat <- netAnalysis_computeCentrality(HD.cellchat, slot.name = "netP")

# 7.2 信号角色散点图
netAnalysis_signalingRole_scatter(HD.cellchat)

# 7.3 信号角色网络
netAnalysis_signalingRole_network(HD.cellchat, signaling = pathways.show,
                                width = 8, height = 2.5, font.size = 10)

# 7.4 进出信号热图
ht1 <- netAnalysis_signalingRole_heatmap(HD.cellchat, pattern = "outgoing", width = 8, height = 12)
ht2 <- netAnalysis_signalingRole_heatmap(HD.cellchat, pattern = "incoming", width = 8, height = 12)
ht1 + ht2

# ============================================================================
# 8. 多组比较（合并分析）
# ============================================================================
MDA.cellchat <- createCellChat(object = GetAssayData(MDA, layer = 'data'),
                               meta = data.frame(group = "MDA", labels = MDA@meta.data$cell_type),
                               group.by = "labels")
MDA.cellchat@DB <- CellChatDB
MDA.cellchat <- subsetData(MDA.cellchat)
MDA.cellchat <- identifyOverExpressedGenes(MDA.cellchat)
MDA.cellchat <- identifyOverExpressedInteractions(MDA.cellchat)
MDA.cellchat <- computeCommunProb(MDA.cellchat, type = "triMean")
MDA.cellchat <- filterCommunication(MDA.cellchat, min.cells = 10)
MDA.cellchat <- computeCommunProbPathway(MDA.cellchat)
MDA.cellchat <- aggregateNet(MDA.cellchat)

# 合并
object.list <- list(HD = HD.cellchat, MDA = MDA.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# 8.1 互作数目比较
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2), measure = "weight")
gg1 + gg2

# 8.2 差异互作网络
netVisual_diffInteraction(cellchat, weight.scale = T)

# 8.3 热图比较
netVisual_heatmap(cellchat)
netVisual_heatmap(cellchat, measure = "weight")

# 8.4 排名网络
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

# 8.5 气泡图比较
netVisual_bubble(cellchat, sources.use = "Tcell", comparison = c(2, 1), angle.x = 45)

# 8.6 差异表达分析（识别上下调信号）
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets",
                                     pos.dataset = "HD", features.name = "HD.merged",
                                     only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,
                                     thresh.p = 0.05, group.DE.combined = T)
net <- netMappingDEG(cellchat, features.name = "HD.merged", variable.all = TRUE)
net.up <- subsetCommunication(cellchat, net = net, datasets = "MDA",
                              ligand.logFC = 0.05, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "HD",
                              ligand.logFC = -0.05, receptor.logFC = NULL)

# ============================================================================
# 9. 保存结果
# ============================================================================
saveRDS(HD.cellchat, "HD.cellchat.rds")
saveRDS(cellchat, "cellchat_merged.rds")
