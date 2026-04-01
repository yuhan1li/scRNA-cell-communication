# scRNA Cell-Cell Communication Analysis

单细胞转录组细胞间通讯（Cell-Cell Communication）分析工具集，支持 6 种主流分析方法。

## 🎯 方法选择

| 方法 | 推荐场景 | 特点 |
|------|---------|------|
| **CellChat V2** | 文章发表首选 | 图形美观，调控网络完整 |
| **Cellphonedb V5** | 严格统计推断 | 数据库 >2000 对 |
| **Nichenet** | 配体→靶基因机制 | 调控链条分析 |
| **MultiNichenet** | 少样本/无重复 | 多组比较专用 |
| **iTALK** | 快速初步探索 | 简单直观 |

## 🚀 快速选择

```
发表文章 → CellChat V2
少样本无重复 → MultiNichenet
严格统计 → Cellphonedb V5
配体→靶基因机制 → Nichenet
初步探索 → iTALK
```

## 📦 安装依赖

### R 环境

```r
# 推荐 R >= 4.2
install.packages("Seurat")
install.packages("tidyverse")
install.packages("patchwork")
```

### Python 环境（Cellphonedb）

```bash
pip install cellphonedb
```

## 📁 文件结构

```
scRNA-cell-communication/
├── README.md
└── references/
    ├── cellchat_v2_tutorial.R          # CellChat V2 完整分析流程
    ├── cellphonedb_v5_prepar_input.R  # Cellphonedb V5 数据准备
    ├── nichenet_tutorial.R             # Nichenet 调控网络
    ├── multinichenet_tutorial.R       # MultiNichenet 少样本分析
    └── italk_tutorial.R               # iTALK 快速可视化
```

## 🔬 快速开始

### CellChat V2

```r
library(CellChat)

seu <- readRDS("your_seurat.rds")
data.input <- GetAssayData(seu, slot = "data")
cellchat <- createCellChat(data.input, group.by = "ident", database = "human")

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.25)
cellchat <- computeCommunProbPathways(cellchat)

# 可视化
netVisual_heatmap(cellchat)
netVisual_bubble(cellchat)
```

### Cellphonedb V5

```r
# 数据导出
counts_matrix <- GetAssayData(seu, slot = "counts")
counts_df <- as.data.frame(t(as.matrix(counts_matrix)))
counts_df$cell_id <- rownames(counts_df)
write.csv(counts_df, "cpdb_counts.csv", row.names = FALSE)

meta_df <- data.frame(cell_id = colnames(seu), cell_type = as.character(Idents(seu)))
write.csv(meta_df, "cpdb_meta.csv", row.names = FALSE)
```

```bash
# Python 端运行
cellphonedb statics cpdb_counts.csv cpdb_meta.csv \
  --counts-data gene_name --output-path output_v5 --threads 8 \
  --iterations 1000 --pvalues-method permtest
```

## ⚠️ 常见问题

- **找不到配体-受体对？** 检查基因名是否为 Symbol（不是 Ensemble ID）
- **CellChat vs Cellphonedb？** 建议两者都用，互为验证
- **可视化太复杂？** 筛选 top N 或特定细胞类型对

## 📚 代码来源

KS-codes(2025) 代码库

## License

MIT
