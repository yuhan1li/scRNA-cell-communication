# scRNA 细胞通讯分析 — 方法选择指南

## 🔬 6 种方法快速选择

| 方法 | 推荐场景 | 参考代码 |
|------|---------|---------|
| **CellChat V2** | 文章发表首选，图形美观 | `references/cellchat_v2_tutorial.R` |
| **Cellphonedb V5** | 严格统计推断 | `references/cellphonedb_v5_prepar_input.R` |
| **Nichenet** | 配体→靶基因调控机制 | `references/nichenet_tutorial.R` |
| **MultiNichenet** | 少样本/无重复/多组比较 | `references/multinichenet_tutorial.R` |
| **iTALK** | 快速初步探索 | `references/italk_tutorial.R` |

## 📊 决策树

```
发表文章 → CellChat V2
少样本无重复 → MultiNichenet
严格统计 → Cellphonedb V5
配体→靶基因 → Nichenet
初步探索 → iTALK
```

## ⚠️ 常见问题

- **找不到配体-受体对？** 检查基因名是否为 Symbol（不是 Ensemble ID）
- **可视化太复杂？** 筛选 top N 或特定细胞类型对
- **CellChat vs Cellphonedb？** 两者都用，互为验证

## 📚 快速启动

```r
# CellChat
library(CellChat)
seu <- readRDS("your_seurat.rds")
data.input <- GetAssayData(seu, slot = "data")
cellchat <- createCellChat(data.input, group.by = "ident", database = "human")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.25)
cellchat <- computeCommunProbPathways(cellchat)
netVisual_heatmap(cellchat)  # 最基础热图
```

## 📂 参考文件

代码来源：KS-codes(2025)
详细参数说明见各 `.R` 文件内注释
