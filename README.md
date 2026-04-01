# scRNA Cell-Cell Communication Analysis

单细胞转录组细胞间通讯（Cell-Cell Communication）分析工具集，支持 6 种主流分析方法。

> 代码来源：KS-codes(2025) 代码库
>
> 完整代码文件位于：`~/github/scRNA-cell-communication/references/`

---

## 📂 完整目录结构

```
1-scRNA 细胞通讯分析/
├── Cellchat V1/
│   └── Cellchat(详细注释版).R              # 单样本 + 多组 CellChat V1 完整流程
├── Cellchat V2 （视频教程）/
│   ├── 【视频教程】cellchat V2详细分析教程及可视化/
│   │   ├── cellchat V2详细版更新教程.R       # CellChat V2 完整分析 + 可视化（含 Figure 1-34）
│   │   ├── cellchat V2 示例数据及视频教程.txt  # 视频链接 + 示例数据下载
│   │   └── Figure/figure1-34.pdf           # 34 张示例图
│   ├── cellchat数据库更新及自定义/
│   │   ├── cellchat数据库更新及自定义.R         # CellChat 自定义数据库构建（整合 CPDB 等）
│   │   └── cpdb数据库文件/                    # CPDB 数据库文件（interaction/gene/complex/protein）
│   └── （VIP）2024-9-4 cellchat多组受配体比较气泡图函数.zip
├── Cellphonedb V3/
│   ├── cellphonedb V3 分析及可视化大全（分析终端操作）.R  # Cellphonedb V3 完整 R 端教程
│   └── cellphonedb受配体可视化气泡图（ggplot2）.zip      # ggplot2 气泡图可视化
├── Cellphonedb V5 （视频教程）/
│   ├── 【视频教程】cpdb V5分析教程/
│   │   ├── prepar_input_for_cpdb5.R             # Cellphonedb V5 数据准备脚本
│   │   ├── 单细胞cellphonedb_v5教程.html       # 详细 HTML 教程
│   │   └── output/                            # V5 统计结果文件
│   ├── cellphonedb数据库更新及自定义/
│   │   ├── ccdb_cpdb_merge/                   # CellChatDB + Cellphonedb 合并
│   │   └── cellphonedb.zip                    # V5 数据库文件
│   └── （VIP）2024-8-17 cellphonedb V5多组受配体分析可视化函数.zip
├── iTALK/
│   └── iTALK受配体差异分析.R                   # iTALK 差异分析 + 环形和弦图
├── multinichentr（视频教程）/
│   ├── 1-multinechenetr适用于有样本重复的数据/
│   │   ├── multinichenetr分步基础分析及可视化教程.R  # 有重复： pseudobulk + EdgeR
│   │   └── 【multinechenetr视频教程】.txt
│   └── 2-multinechenetr适用于无样本重复或者样本较少的数据/
│       └── multinichenetr无重复样本多组受配体差异分析.R  # 无重复：FindMarkers
└── Nichenet V2（视频教程）/
    ├── 1-nichenetr基本分析流程【视频演示】/
    ├── 2-nichenetr分析流程及基本可视化修饰【视频演示】/
    ├── 3-nichenetr配体靶基因信号推断及网络可视化/
    └── 4-nichenetr分析之基于ligand-receptor的ligands排序【视频演示】/
```

---

## 🎯 方法选择指南

| 方法 | 推荐场景 | 特点 |
|------|---------|------|
| **CellChat V2** | 文章发表首选，综合分析 | 图形美观，组间比较，多种可视化 |
| **Cellphonedb V5** | 严格统计推断 | 数据库大，统计严格 |
| **Nichenet V2** | 配体→靶基因机制 | 调控链条分析 |
| **MultiNichenet** | 少样本/无重复 | 多组比较专用，无需生物学重复 |
| **iTALK** | 快速初步探索 | 简单直观 |

## 🚀 快速决策

```
发表文章 → CellChat V2
少样本/无重复 → MultiNichenet
严格统计推断 → Cellphonedb V5
配体→靶基因调控 → Nichenet V2
多组比较分析 → MultiNichenet
快速初步探索 → iTALK
数据库自定义 → CellChat V2 / Cellphonedb V5
```

---

## 📦 依赖环境

### R 环境（推荐 R >= 4.2）

```r
# 基础包
install.packages(c("Seurat", "tidyverse", "patchwork", "SingleCellExperiment"))

# CellChat V2（推荐）
devtools::install_github("immunogenomics/presto")
devtools::install_github("jinworks/CellChat")  # V2 最新

# MultiNichenet
install.packages("multinichenetr")

# Nichenet V2
devtools::install_github("jinworks/nichenetr")
```

### Python 环境（Cellphonedb）

```bash
# 推荐创建独立 conda 环境
conda create -n cellphonedb python=3.9
conda activate cellphonedb
pip install cellphonedb -i https://pypi.douban.com/simple --trusted-host pypi.douban.com

# 验证安装
cellphonedb --help
```

---

## 🔬 方法一：CellChat V2（推荐首选）

> CellChat V2 vs V1 主要区别：
> - V2 新增 `updateCellChatDB()` 支持自定义数据库（整合 CPDB 等）
> - V2 `computeCommunProb()` 支持 `type = "triMean"`（默认）和 `"truncatedMean"`
> - V2 `mergeCellChat()` 支持多组合并比较
> - V2 网络分析更加丰富（netAnalysis 系列函数）

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
# V2 推荐写法
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

# 并行计算加速
future::plan("multisession", workers = 2)

# 识别过表达基因和互作
HD.cellchat <- identifyOverExpressedGenes(HD.cellchat)
HD.cellchat <- identifyOverExpressedInteractions(HD.cellchat)

# 计算通讯概率（V2 推荐 triMean）
HD.cellchat <- computeCommunProb(HD.cellchat, type = "triMean")

# 过滤低可信度通讯
HD.cellchat <- filterCommunication(HD.cellchat, min.cells = 10)

# 计算通路级别通讯
HD.cellchat <- computeCommunProbPathway(HD.cellchat)

# 聚合网络
HD.cellchat <- aggregateNet(HD.cellchat)
```

### 可视化（完整示例）

```r
# 1. 互作数目 Circle 图
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(HD.cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F,
                 title.name = "Number of interactions")
netVisual_circle(HD.cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F,
                 title.name = "Interaction weights/strength")

# 2. 热图展示互作数据
pheatmap::pheatmap(HD.cellchat@net$count, border_color = "black",
                   cluster_cols = F, fontsize = 10, cluster_rows = F,
                   display_numbers = T, number_color = "black",
                   number_format = "%.0f")

# 3. 指定通路可视化
pathways.show <- c("MHC-II")

# Circle 图
netVisual_aggregate(HD.cellchat, signaling = pathways.show)

# 和弦图
netVisual_aggregate(HD.cellchat, signaling = pathways.show, layout = "chord")

# 热图
netVisual_heatmap(HD.cellchat, signaling = pathways.show, color.heatmap = "Reds")

# 4. CellPhoneDB style 受配体气泡图
netVisual_bubble(HD.cellchat, sources.use = "Tcell", remove.isolate = FALSE)

# 5. 分组和弦图
group.cellType <- c(rep("immune", 5), rep("other", 4))
names(group.cellType) <- c("Kers", "Mast", "Men", "Mon", "Tcell", "ECs", "Fibs", "lang", "SMCs")
netVisual_chord_cell(HD.cellchat, signaling = pathways.show, group = group.cellType)

# 6. 提取特定通路的受配体对
pairLR <- extractEnrichedLR(HD.cellchat, signaling = pathways.show, geneLR.return = FALSE)
netVisual_individual(HD.cellchat, signaling = pathways.show,
                     pairLR.use = "HLA-DPA1_CD4", layout = "circle")
```

### 系统分析（Signal Role + Pattern）

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

# 保存
saveRDS(HD.cellchat, "HD.cellchat.rds")
save(cellchat, file = "cellchat_merge.RData")
```

### CellChat 数据库更新与自定义

```r
# 查看 updateCellChatDB 函数参数
help(updateCellChatDB)
# updateCellChatDB(
#   db,                    # 数据框，至少需要 ligand 和 receptor 两列
#   gene_info = NULL,      # 基因信息（至少一列 Symbol），设为 NULL 时需指定物种
#   other_info = NULL,      # 其他信息（如 complex_input）
#   trim.pathway = FALSE,   # 是否删除无 pathway 注释的互作对
#   merged = FALSE,         # 是否与现有 CellChat 数据库合并
#   species_target = NULL    # 物种： "human" 或 "mouse"
# )

# 整合 CPDB 数据库（人源）
library(CellChat)

# 读取 CPDB 数据库文件
interaction_input <- read.csv('./cpdb数据库文件/interaction_input.csv')
complex_input <- read.csv('./cpdb数据库文件/complex_input.csv', row.names = 1)
geneInfo <- read.csv('./cpdb数据库文件/gene_input.csv')

# 整理格式
geneInfo$Symbol <- geneInfo$hgnc_symbol
geneInfo <- unique(geneInfo)

# 受配体信息
idx_partnerA <- match(interaction_input$partner_a, geneInfo$uniprot)
interaction_input$ligand <- interaction_input$partner_a
interaction_input$ligand[!is.na(idx_partnerA)] <- geneInfo$hgnc_symbol[idx_partnerA[!is.na(idx_partnerA)]]

idx_partnerB <- match(interaction_input$partner_b, geneInfo$uniprot)
interaction_input$receptor <- interaction_input$partner_b
interaction_input$receptor[!is.na(idx_partnerB)] <- geneInfo$hgnc_symbol[idx_partnerB[!is.na(idx_partnerB)]]

interaction_input$interaction_name <- interaction_input$interactors
interaction_input$pathway_name <- gsub(".*by ", "", interaction_input$classification)

# 仅用 CPDB 数据库
db.new.cpdb <- updateCellChatDB(db = interaction_input, gene_info = geneInfo,
                                other_info = list(complex = complex_input),
                                trim.pathway = T, species_target = "human", merged = F)

# CPDB + CellChatDB 合并
db.new.cpdb_ccdb <- updateCellChatDB(db = interaction_input, gene_info = geneInfo,
                                      other_info = list(complex = complex_input),
                                      trim.pathway = T, species_target = "human", merged = T)

# 保存
save(db.new.cpdb, file = 'db_new_cpdb.RData')
save(db.new.cpdb_ccdb, file = 'db.new.cpdb_ccdb.RData')

# 分析时选择自定义数据库
object <- createCellChat(...)
object@DB <- db.new.cpdb  # 使用 CPDB 数据库
```

---

## 🔬 方法二：Cellphonedb V5

> Cellphonedb V5 vs V3 主要区别：
> - V5 支持 `cellphonedb statics`（有重复）和 `cellphonedb statics-single`（无重复）
> - V5 数据库更大，包含更多受配体对
> - V5 使用 Python API，`prepar_input_for_cpdb5.R` 准备数据

### R 端数据准备

```r
library(Seurat)

# 方法 1：手动准备
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

# 方法 2：使用官方准备脚本
source("prepar_input_for_cpdb5.R")
```

### Python 端运行

```bash
# 有生物学重复（推荐 iterations=1000）
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

# 输出文件说明
# statistical_analysis_deconvoluted_*.txt      # 去卷积丰度
# statistical_analysis_means_*.txt             # 平均表达值
# statistical_analysis_pvalues_*.txt            # 置换检验 p 值
# statistical_analysis_significant_means_*.txt  # 显著互作（p < 0.05）的平均表达值
```

### 可视化（ggplot2 气泡图）

```r
library(tidyverse)
library(ggplot2)

# 读取结果文件
all_pval <- read.table("pvalues.txt", header = T, stringsAsFactors = F,
                       sep = '\t', comment.char = '', check.names = F)
all_means <- read.table("means.txt", header = T, stringsAsFactors = F,
                        sep = '\t', comment.char = '', check.names = F)

intr_pairs <- all_pval$interacting_pair
all_pval <- all_pval[, -c(1:11)]
all_means <- all_means[, -c(1:11)]

# 选择特定细胞互作（如 T_cell → 其他）
selected_celltype <- c("T_cell|DC", "T_cell|Endo", "T_cell|Fib",
                       "T_cell|Mon", "T_cell|NK")

# 筛选显著互作（p < 0.05）
sig_pairs <- all_pval
sig_pairs <- sig_pairs[, -c(1, 3:11)]
sig_pairs <- sig_pairs[rowSums(sig_pairs <= 0.05) != 0, ]

# 选择展示的受配体对
selected_pairs <- sig_pairs$interacting_pair[1:40]

# 合并 mean 与 pvalue
sel_pval <- all_pval[match(selected_pairs, intr_pairs), selected_celltype]
sel_means <- all_means[match(selected_pairs, intr_pairs), selected_celltype]
df_names <- expand.grid(selected_pairs, selected_celltype)
pval <- unlist(sel_pval)
pval[pval == 0] <- 0.00001
plot.data <- cbind(df_names, pval)

plot.data$mean <- unlist(sel_means)
colnames(plot.data) <- c("interacting_pair", "celltype_pair", "pval", "mean")

# 绘制气泡图
ggplot(plot.data, aes(x = celltype_pair, y = interacting_pair)) +
  geom_point(aes(size = -log10(pval), color = mean)) +
  scale_color_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

---

## 🔬 方法三：MultiNichenet（无重复/少样本）

> **核心优势：无需生物学重复，适合临床数据（每组 < 3 个样本）**
>
> 两种模式：
> - **有重复**：使用 pseudobulk + EdgeR（`get_abundance_info` + pseudobulk）
> - **无重复**：使用 FindMarkers 差异分析

### 数据准备

```r
library(SingleCellExperiment)
library(multinichenetr)
library(Seurat)

# Seurat → SingleCellExperiment 转换
sce <- Seurat::as.SingleCellExperiment(sce, assay = "RNA")
sce <- alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()

# 加载网络文件（从 MultiNichenet 官网下载）
lr_network_all <- readRDS("lr_network_human_allInfo.rds") %>%
  mutate(ligand = convert_alias_to_symbols(ligand, organism = "human"),
         receptor = convert_alias_to_symbols(receptor, organism = "human")) %>%
  mutate(ligand = make.names(ligand),
         receptor = make.names(receptor))
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

# make.names 确保有效变量名
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

### 分析流程（有重复：pseudobulk + EdgeR）

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

# Step 4: Pseudobulk 差异分析
DE_info <- get_DE_info(
  sce = sce, sample_id = sample_id, group_id = group_id,
  celltype_id = celltype_id, contrasts_oi = contrasts_oi,
  min_cells = min_cells, fraction_cutoff = fraction_cutoff,
  min_sample_prop = min_sample_prop, contrast_tbl = contrast_tbl,
  batches = batches, per_celltype = FALSE)

celltype_de <- DE_info$celltype_de

sender_receiver_de <- multinichenetr::combine_sender_receiver_de(
  sender_de = celltype_de, receiver_de = celltype_de,
  senders_oi = senders_oi, receivers_oi = receivers_oi,
  lr_network = lr_network)

# Step 5: 配体活性预测
logFC_threshold <- 0.3
p_val_threshold <- 0.05
top_n_target <- 250

ligand_activities_targets_DEgenes <- suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold, p_val_threshold = p_val_threshold,
    p_val_adj = TRUE, top_n_target = top_n_target, verbose = F, n.cores = 2
  )
))

# Step 6: 优先级排序
prioritization_tables <- suppressMessages(multinichenetr::generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_de %>% distinct(sender, receiver),
  grouping_tbl = grouping_tbl, scenario = "no_frac_LR_expr",
  fraction_cutoff = fraction_cutoff
))

# Step 7: 保存
saveRDS(prioritization_tables, "multinichenet_output.rds")
```

### 分析流程（无重复：FindMarkers）

```r
# 与有重复流程类似，Step 4 替换为：
DE_info <- get_DE_info_sampleAgnostic(
  sce = sce, group_id = group_id, celltype_id = celltype_id,
  contrasts_oi = contrasts_oi, expressed_df = frq_list$expressed_df,
  min_cells = min_cells, contrast_tbl = contrast_tbl)

celltype_de <- DE_info$celltype_de_findmarkers
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

## 🔬 方法四：Nichenet V2

### 基本分析流程

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

> **注意：iTALK 不提供统计显著性检验，仅适合快速初步探索**

```r
library(iTALK)
library(Seurat)
library(igraph)

seu <- readRDS("your_scRNA.rds")
expr <- GetAssayData(seu, slot = "data")
expr$cell_type <- seu@meta.data$celltype
expr$compare_group <- seu@meta.data$orig.ident

# 两组之间差异基因分析
deg_endo <- DEG(expr %>% filter(cell_type == 'Endothelial'),
                method = 'Wilcox', contrast = c("AEH", "HC"))
deg_endo$cell_type <- 'Endothelial'  # 必须有，否则报错

deg_ly <- DEG(expr %>% filter(cell_type == 'Lymphocytes'),
              method = 'Wilcox', contrast = c("AEH", "HC"))
deg_ly$cell_type <- 'Lymphocyte'

# 寻找受配体对
res <- FindLR(deg_ly, deg_endo, datatype = 'DEG', comm_type = 'growth factor')
res <- res[order(res$cell_from_logFC * res$cell_to_logFC, decreasing = TRUE), ]

# 环形图
LRPlot(res[1:60, ],
       datatype = 'DEG',
       cell_col = cell_col,
       link.arr.lwd = 1,
       print.cell = TRUE,
       track.height_1 = uh(1, "mm"),
       track.height_2 = uh(15, "mm"),
       text.vjust = "0.5cm")
```

---

## ⚠️ 常见问题

| 问题 | 解决方案 |
|------|---------|
| 找不到配体-受体对 | 检查基因名是否为 Symbol（不是 Ensemble ID） |
| CellChat 结果太少 | 尝试降低 `computeCommunProb` 的 `trim` 参数，或换用 `truncatedMean` |
| Cellphonedb 运行报错 | 确保 Python 端 cellphonedb 版本 >= 5.0 |
| MultiNichenet 基因过滤过严 | 降低 `fraction_cutoff`（默认 0.05），或 `min_sample_prop` |
| 内存不足 | 使用 `future::plan("multisession")` 并行计算，减少 `workers` 数量 |
| Cellphonedb V5 无重复报错 | 使用 `cellphonedb statics-single` 命令替代 `statics` |
| iTALK 作图报错 | 确保 `cell_type` 列名正确：celltype 必须重命名以符合 iTALK 要求 |

---

## 📚 方法对比

| 特性 | CellChat V2 | Cellphonedb V5 | MultiNichenet | Nichenet V2 | iTALK |
|------|------------|----------------|----------------|-------------|-------|
| **统计学方法** | 概率模型 | 置换检验 | Pseudobulk/FindMarkers | 调控网络推断 | 无统计 |
| **重复需求** | 需要 | 需要 | **无需重复** | 需要 | 需要 |
| **多组比较** | ✅ | ❌ | ✅ | ❌ | ❌ |
| **配体→靶基因** | ❌ | ❌ | ✅ | ✅ | ❌ |
| **图形美观度** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐ |
| **发表友好** | ✅ | ✅ | ✅ | ✅ | ❌ |
| **数据库自定义** | ✅ | ✅ | ❌ | ❌ | ❌ |

---

## 📌 数据库下载链接

| 工具 | 下载说明 |
|------|---------|
| **CellChatDB** | 内置，直接加载 `CellChatDB.human` / `CellChatDB.mouse` |
| **Cellphonedb** | `pip install cellphonedb`；数据库文件：https://github.com/ventolab/cellphonedb-data |
| **MultiNichenet** | 链接见视频教程文件夹 `【multinechenetr数据库模型】.txt` |
| **Nichenet** | 链接见 `nichenet-modle-data-for human -and - mouse.txt` |

---

## License

MIT
