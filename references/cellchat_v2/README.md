# CellChat V2 分析指南

## 方法简介

CellChat V2 是目前最流行的细胞通讯分析方法，通过概率模型推断细胞间的配体-受体通讯强度。图形美观，自带多种可视化方式，支持多组比较，是文章发表的首选工具。

**核心原理**：基于配体和受体在细胞组中的平均表达量，结合 CellChatDB 数据库（~3300 对经验证的相互作用），计算通讯概率。

## 输入文件

| 文件 | 格式 | 说明 |
|------|------|------|
| Seurat 对象 | `.rds` | 包含 `RNA` assay 的 SCT 降维结果，细胞注释在 `Idents()` 中 |

**示例数据**：

```r
# 加载 Seurat 对象
seu <- readRDS("scRNA_annotated.rds")
table(Idents(seu))  # 查看细胞类型和数量

# 分组提取示例（疾病 vs 对照）
disease <- subset(seu, group == "Disease")
control <- subset(seu, group == "Normal")
```

## 分析流程

1. **数据准备** → GetAssayData 提取表达矩阵
2. **创建 CellChat 对象** → `createCellChat()`
3. **设置数据库** → `CellChatDB.human` 或 `.mouse`
4. **预处理** → `subsetData()` + `identifyOverExpressedGenes()`
5. **计算通讯概率** → `computeCommunProb()`
6. **聚合网络** → `aggregateNet()`
7. **可视化** → 6 种图形
8. **多组比较**（可选）→ `mergeCellChat()` + 差异分析

## 输出结果

| 文件 | 内容 |
|------|------|
| `net@count` | 细胞类型对之间的互作数目 |
| `net@weight` | 互作强度（权重）|
| `netP@pathways` | 信号通路列表 |
| `netP@prob` | 配体-受体对通讯概率矩阵 |
| `.rds` | 完整的 CellChat 对象（可重复使用）|

## 适用场景

- ✅ 多组比较（疾病 vs 对照）
- ✅ 细胞类型较多的数据（> 5 种细胞类型）
- ✅ 需要美观图形发表文章
- ✅ 需要识别主要信号发送者/接收者
- ✅ 需要通路富集分析

## 注意事项

- CellChatDB v2 包含约 3300 个验证过的相互作用
- 默认使用 `triMean` 方法，产生的互作较少但置信度高
- 如需更多互作，可尝试 `truncatedMean` 方法
- 建议每组至少 3 个生物学重复
