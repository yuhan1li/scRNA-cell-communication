# scRNA Cell-Cell Communication Analysis

单细胞转录组细胞间通讯（Cell-Cell Communication）分析工具集，支持 6 种主流分析方法。

> 代码来源：KS-codes(2025) 代码库

---

## 🎯 方法选择指南

| 方法 | 推荐场景 | 特点 |
|------|---------|------|
| **CellChat V2** | 文章发表首选，综合分析 | 图形美观，组间比较，多种可视化 |
| **Cellphonedb V5** | 严格统计推断 | 数据库大，统计严格 |
| **Nichenet** | 配体→靶基因机制 | 调控链条分析 |
| **MultiNichenet** | 少样本/无重复 | 多组比较专用，无需生物学重复 |
| **iTALK** | 快速初步探索 | 简单直观 |

## 🚀 快速决策

```
发表文章 → CellChat V2
少样本/无重复 → MultiNichenet
严格统计推断 → Cellphonedb V5
配体→靶基因调控 → Nichenet
多组比较分析 → MultiNichenet
快速初步探索 → iTALK
```

---

## 📦 依赖环境

### R 环境（推荐 R >= 4.2）

```r
# 基础包
install.packages(c("Seurat", "tidyverse", "patchwork", "SingleCellExperiment"))

# CellChat V2
devtools::install_github("immunogenomics/presto")  # 可能需要的依赖
devtools::install_github("jinworks/CellChat")

# MultiNichenet
install.packages("multinichenetr")

# Nichenet
devtools::install_github("jinworks/nichenetr")
```

### Python 环境（Cellphonedb）

```bash
pip install cellphonedb
```

---

## 📁 文件结构

```
scRNA-cell-communication/
├── README.md
└── references/
    ├── cellchat_v2_tutorial.R          # CellChat V2 完整分析流程
    ├── cellphonedb_v5_prepar_input.R  # Cellphonedb V5 数据准备
    ├── multinichenet_tutorial.R       # MultiNichenet 无重复分析
    ├── nichenet_tutorial.R             # Nichenet 调控网络
    └── italk_tutorial.R               # iTALK 快速可视化
```

---

## 🔬 方法一：CellChat V2（推荐首选）

### 数据准备

```r
library(CellChat)
library(Seurat)

# 加载数据
sceV5 <- load("sceV5.RData")

# 分组提取（示例：HD vs MDA5）
HD <- subset(sceV5, group1 == 'HD')
MDA <- subset(sceV5, group1 == 'MDA5')

# 设置 assay
DefaultAssay(HD) <- 'SCT'

# 获取表达矩阵和元数据
HD_input <- GetAssayData(HD, layer = 'data')
HD_meta <- HD@meta.data[, c("group1", "cell_type")]
colnames(HD_meta) <- c("group", "labels")

# 检查一致性
identical(colnames(HD_input), rownames(HD_meta))  # 应返回 TRUE
```

### 创建 CellChat 对象

```r
HD.cellchat <- createCellChat(object = HD_input, meta = HD_meta, group.by = "labels")

# 设置数据库（人源数据）
CellChatDB <- CellChatDB.human
HD.cellchat@DB <- CellChatDB

# 查看细胞类型
levels(HD.cellchat@idents)
groupSize <- as.numeric(table(HD.cellchat@idents))
```

### 预处理与计算

```r
# 预处理
HD.cellchat <- subsetData(HD.cellchat)
future::plan("multisession", workers = 2)  # 并行计算

# 识别过表达基因和互作
HD.cellchat <- identifyOverExpressedGenes(HD.cellchat)
HD.cellchat <- identifyOverExpressedInteractions(HD.cellchat)

# 计算通讯概率（triMean 方法）
HD.cellchat <- computeCommunProb(HD.cellchat, type = "triMean")

# 过滤低可信度通讯
HD.cellchat <- filterCommunication(HD.cellchat, min.cells = 10)

# 计算通路级别通讯
HD.cellchat <- computeCommunProbPathway(HD.cellchat)

# 聚合网络
HD.cellchat <- aggregateNet(HD.cellchat)
```

### 可视化

```r
# 1. 互作数目 Circle 图
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(HD.cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F,
                 title.name = "Number of interactions")
netVisual_circle(HD.cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F,
                 title.name = "Interaction weights/strength")

# 2. 热图
pheatmap::pheatmap(HD.cellchat@net$count, cluster_cols = F, cluster_rows = F,
                   fontsize = 10, display_numbers = T, number_color = "black")

# 3. 指定通路可视化（Circle / Chord / Heatmap）
pathways.show <- c("MHC-II")

# Circle 图
netVisual_aggregate(HD.cellchat, signaling = pathways.show)

# 和弦图
netVisual_aggregate(HD.cellchat, signaling = pathways.show, layout = "chord")

# 热图
netVisual_heatmap(HD.cellchat, signaling = pathways.show, color.heatmap = "Reds")

# 4. 分组和弦图
group.cellType <- c(rep("immune", 5), rep("other", 4))
names(group.cellType) <- c("Kers", "Mast", "Men", "Mon", "Tcell", "ECs", "Fibs", "lang", "SMCs")
netVisual_chord_cell(HD.cellchat, signaling = pathways.show, group = group.cellType)

# 5. 受配体气泡图
netVisual_bubble(HD.cellchat, sources.use = "Tcell", remove.isolate = FALSE)

# 6. 提取特定通路的受配体对
pairLR <- extractEnrichedLR(HD.cellchat, signaling = pathways.show, geneLR.return = FALSE)
netVisual_individual(HD.cellchat, signaling = pathways.show,
                     pairLR.use = "HLA-DPA1_CD4", layout = "circle")
```

### 系统分析

```r
# 网络中心性分析
HD.cellchat <- netAnalysis_computeCentrality(HD.cellchat, slot.name = "netP")

# 信号角色散点图
netAnalysis_signalingRole_scatter(HD.cellchat)

# 特定通路信号角色网络
netAnalysis_signalingRole_network(HD.cellchat, signaling = pathways.show,
                                  width = 8, height = 2.5, font.size = 10)

# 进出信号热图
ht1 <- netAnalysis_signalingRole_heatmap(HD.cellchat, pattern = "outgoing", width = 8, height = 12)
ht2 <- netAnalysis_signalingRole_heatmap(HD.cellchat, pattern = "incoming", width = 8, height = 12)
ht1 + ht2

# NMF 模式识别
library(NMF)
selectK(HD.cellchat, pattern = "outgoing")  # 选择合适的 k 值

nPatterns = 2
HD.cellchat <- identifyCommunicationPatterns(HD.cellchat, pattern = "outgoing", k = nPatterns)

# 冲击图
netAnalysis_river(HD.cellchat, pattern = "outgoing")

# 气泡图
netAnalysis_dot(HD.cellchat, pattern = "outgoing")

# 功能/结构相似性聚类
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
```

### 多组比较

```r
# 合并多组 CellChat 对象
object.list <- list(HD = HD.cellchat, MDA = MDA.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# 互作数目比较
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2), measure = "weight")
gg1 + gg2

# 差异互作网络
netVisual_diffInteraction(cellchat, weight.scale = T)

# 热图比较
netVisual_heatmap(cellchat)
netVisual_heatmap(cellchat, measure = "weight")

# 排名网络
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

# 气泡图比较
netVisual_bubble(cellchat, sources.use = "ECs", comparison = c(2, 1), angle.x = 45)

# 差异表达分析
pos.dataset = "HD"
features.name = paste0(pos.dataset, ".merged")
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets",
                                       pos.dataset = pos.dataset, features.name = features.name,
                                       only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,
                                       thresh.p = 0.05, group.DE.combined = T)

# 获取上下调信号
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
net.up <- subsetCommunication(cellchat, net = net, datasets = "MDA",
                              ligand.logFC = 0.05, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "HD",
                                ligand.logFC = -0.05, receptor.logFC = NULL)

# 可视化上下调
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = "ECs",
                        comparison = c(1, 2), angle.x = 90,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = "ECs",
                        comparison = c(1, 2), angle.x = 90,
                        title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

# 保存
saveRDS(HD.cellchat, "HD.cellchat.rds")
save(cellchat, file = "cellchat_merge.RData")
```

---

## 🔬 方法二：Cellphonedb V5

### R 端数据准备

```r
library(Seurat)

# 归一化数据
Normalized_counts <- GetAssayData(scRNA, layer = 'data')
Normalized_counts <- as.data.frame(Normalized_counts)
Normalized_counts <- cbind(rownames(Normalized_counts), Normalized_counts)
colnames(Normalized_counts)[1] <- 'Gene'

# 元数据
metadata <- data.frame(Cell = rownames(scRNA@meta.data),
                       cell_type = scRNA$celltype)
metadata$cell_type <- gsub(' ', '_', metadata$cell_type)

# 保存（制表符分隔）
write.table(Normalized_counts, "Normalized_counts.txt",
            row.names = F, sep = '\t', quote = F)
write.table(metadata, "cellphonedb_meta.txt",
            row.names = F, sep = '\t', quote = F)
```

### Python 端运行

```bash
# 有生物学重复
cellphonedb statics Normalized_counts.txt cellphonedb_meta.txt \
  --counts-data hgnc_symbol \
  --output-path output_v5 \
  --threads 8 \
  --iterations 1000 \
  --pvalues-method permtest

# 无生物学重复
cellphonedb statics-single Normalized_counts.txt cellphonedb_meta.txt \
  --counts-data hgnc_symbol \
  --output-path output_single \
  --threads 8
```

### 结果解读

Cellphonedb V5 输出文件：
- `statistical_analysis_significant_means_*.txt` — 显著互作的平均表达值
- `statistical_analysis_pvalues_*.txt` — 配对置换检验 p 值
- 列名为 `Celltype_A|Celltype_B` 格式的行表示 A→B 的通讯

---

## 🔬 方法三：MultiNichenet（无重复/少样本）

> 适用：每组样本数 < 3 或无生物学重复的临床数据

### 数据准备

```r
library(SingleCellExperiment)
library(multinichenetr)
library(Seurat)

# 加载示例数据
load("scRNA_nichenet.RData")

# Seurat → SingleCellExperiment 转换
sce <- Seurat::as.SingleCellExperiment(sce, assay = "RNA")
sce <- alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()

# 加载网络文件（需从 MultiNichenet 官网下载）
lr_network_all <- readRDS("lr_network_human_allInfo.rds") %>%
  mutate(ligand = convert_alias_to_symbols(ligand, organism = "human"),
         receptor = convert_alias_to_symbols(receptor, organism = "human")) %>%
  mutate(ligand = make.names(ligand), receptor = make.names(receptor))
lr_network <- lr_network_all %>% distinct(ligand, receptor)

ligand_target_matrix <- readRDS("ligand_target_matrix.rds")
colnames(ligand_target_matrix) <- colnames(ligand_target_matrix) %>%
  convert_alias_to_symbols(organism = "human") %>% make.names()
rownames(ligand_target_matrix) <- rownames(ligand_target_matrix) %>%
  convert_alias_to_symbols(organism = "human") %>% make.names()
lr_network <- lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
ligand_target_matrix <- ligand_target_matrix[, lr_network$ligand %>% unique()]
```

### 设置比较参数

```r
# 定义分组信息
sample_id <- "orig.ident"
group_id <- "group"
celltype_id <- "celltype"
batches <- NA

# 整理 colData
SummarizedExperiment::colData(sce)$orig.ident <- SummarizedExperiment::colData(sce)$orig.ident %>% make.names()
SummarizedExperiment::colData(sce)$group <- SummarizedExperiment::colData(sce)$group %>% make.names()
SummarizedExperiment::colData(sce)$celltype <- SummarizedExperiment::colData(sce)$celltype %>% make.names()

# 定义比较组（如 SD vs HC）
contrasts_oi <- c("'SD-HC','HC-SD'")
contrast_tbl <- tibble(contrast = c("SD-HC", "HC-SD"), group = c("SD", "HC"))

# 定义 sender / receiver
senders_oi <- SummarizedExperiment::colData(sce)[, celltype_id] %>% unique()
receivers_oi <- SummarizedExperiment::colData(sce)[, celltype_id] %>% unique()
```

### 分析流程

```r
# Step 1: 细胞类型过滤
min_cells <- 10
abundance_info <- get_abundance_info(
  sce = sce, sample_id = sample_id, group_id = group_id,
  celltype_id = celltype_id, min_cells = min_cells,
  senders_oi = senders_oi, receivers_oi = receivers_oi, batches = batches)

# Step 2: 基因过滤
fraction_cutoff <- 0.05
min_sample_prop <- 1
frq_list <- get_frac_exprs_sampleAgnostic(
  sce = sce, sample_id = sample_id, celltype_id = celltype_id,
  group_id = group_id, batches = batches, min_cells = min_cells,
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)

genes_oi <- frq_list$expressed_df %>% filter(expressed == TRUE) %>% pull(gene) %>% unique()
sce <- sce[genes_oi, ]

# Step 3: 表达量处理
abundance_expression_info <- process_abundance_expression_info(
  sce = sce, sample_id = group_id, group_id = group_id,
  celltype_id = celltype_id, min_cells = min_cells,
  senders_oi = senders_oi, receivers_oi = receivers_oi,
  lr_network = lr_network, batches = batches, frq_list = frq_list,
  abundance_info = abundance_info)

# Step 4: 差异分析（无重复）
DE_info <- get_DE_info_sampleAgnostic(
  sce = sce, group_id = group_id, celltype_id = celltype_id,
  contrasts_oi = contrasts_oi, expressed_df = frq_list$expressed_df,
  min_cells = min_cells, contrast_tbl = contrast_tbl)

celltype_de <- DE_info$celltype_de_findmarkers

sender_receiver_de <- multinichenetr::combine_sender_receiver_de(
  sender_de = celltype_de, receiver_de = celltype_de,
  senders_oi = senders_oi, receivers_oi = receivers_oi, lr_network = lr_network)

# Step 5: 配体活性预测
logFC_threshold <- 0.3
p_val_threshold <- 0.05
p_val_adj <- TRUE
top_n_target <- 250

ligand_activities_targets_DEgenes <- suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold, p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj, top_n_target = top_n_target, verbose = F, n.cores = 2
  )
))

# Step 6: 优先级排序
metadata_combined <- SummarizedExperiment::colData(sce) %>% tibble::as_tibble()
grouping_tbl <- metadata_combined[, c(group_id)] %>% tibble::as_tibble() %>%
  distinct() %>% rename(group = group_id) %>% mutate(sample = group)

prioritization_tables <- suppressMessages(multinichenetr::generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_de %>% distinct(sender, receiver),
  grouping_tbl = grouping_tbl, scenario = "no_frac_LR_expr",
  fraction_cutoff = fraction_cutoff,
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = FALSE
))

# Step 7: 保存
multinichenet_output <- list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl
)
saveRDS(multinichenet_output, "multinichenet_output.rds")
```

### 可视化

```r
# Circos 比较图
prioritized_tbl_oi_all <- get_top_n_lr_pairs(prioritization_tables, top_n = 50, rank_per_group = F)
prioritized_tbl_oi <- prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>%
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] <- 0

senders_receivers <- union(prioritized_tbl_oi$sender, prioritized_tbl_oi$receiver) %>% sort()
colors_sender <- RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>%
  magrittr::set_names(senders_receivers)
colors_receiver <- RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>%
  magrittr::set_names(senders_receivers)

circos_list <- make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)

# 配体-受体-活性-靶基因综合图
group_oi <- "SD"
prioritized_tbl_oi_M_30 <- get_top_n_lr_pairs(prioritization_tables, top_n = 30, groups_oi = group_oi)
plot_oi <- make_sample_lr_prod_activity_plots(prioritization_tables, prioritized_tbl_oi_M_30)
plot_oi
```

---

## 🔬 方法四：Nichenet

```r
library(nichenetr)
library(Seurat)

seu <- readRDS("nichenet_seurat.rds")
load("lr_network_human_allInfo.rds")
load("ligand_target_matrix.rds")

# 感兴趣的靶基因（来自差异分析）
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

# Circos 图
make_circos_group_comparison(seu, ligand_activities, ligand_target_matrix,
  group1 = "treated", group2 = "control")

# 配体-受体-靶基因综合图
make_lr_prod_activity_plots(seu, ligand_activities, top_ligands = 10)

# 保存
saveRDS(ligand_activities, "nichenet_ligand_activities.rds")
```

---

## 🔬 方法五：iTALK（快速探索）

```r
library(iTALK)
library(Seurat)
library(igraph)

seu <- readRDS("your_scRNA.rds")
expr <- GetAssayData(seu, slot = "data")

# 差异基因
de_genes <- FindAllMarkers(seu, only.pos = TRUE, logfc.threshold = 0.25)
ligands <- de_genes %>% filter(cluster %in% c("Macrophage", "Fibroblast")) %>% pull(gene)
receptors <- de_genes %>% filter(cluster %in% c("T_cell", "Cancer_cell")) %>% pull(gene)

# 构建网络
edges <- data.frame(
  from = sample(ligands, 50, replace = TRUE),
  to = sample(receptors, 50, replace = TRUE),
  weight = runif(50, 0.1, 1)
)

g <- graph_from_data_frame(edges, directed = TRUE)
plot(g, vertex.size = 10, vertex.color = "lightblue",
     edge.width = E(g)$weight * 3, layout = layout_with_kk)

# 注意：iTALK 不提供统计显著性检验，仅适合快速初步探索
```

---

## ⚠️ 常见问题

| 问题 | 解决方案 |
|------|---------|
| 找不到配体-受体对 | 检查基因名是否为 Symbol（不是 Ensemble ID） |
| CellChat 结果太少 | 尝试降低 `computeCommunProb` 的 `trim` 参数，或换用 `truncatedMean` |
| Cellphonedb 运行报错 | 确保 Python 端 cellphonedb 版本 >= 5.0 |
| MultiNichenet 基因过滤过严 | 降低 `fraction_cutoff`（默认 0.05）|
| 内存不足 | 使用 `future::plan("multisession")` 并行计算，减少 `workers` 数量 |

---

## 📚 方法对比

| 特性 | CellChat V2 | Cellphonedb V5 | MultiNichenet | Nichenet | iTALK |
|------|------------|----------------|----------------|----------|-------|
| **统计学方法** | 概率模型 | 置换检验 | N/A（差异分析）| 调控网络推断 | 无统计 |
| **重复需求** | 需要 | 需要 | **无需重复** | 需要 | 需要 |
| **多组比较** | ✅ | ❌ | ✅ | ❌ | ❌ |
| **配体→靶基因** | ❌ | ❌ | ✅ | ✅ | ❌ |
| **图形美观度** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐ |
| **发表友好** | ✅ | ✅ | ✅ | ✅ | ❌ |

---

## License

MIT
