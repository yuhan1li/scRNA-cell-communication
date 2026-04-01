# Nichenet 分析指南

## 方法简介

Nichenet 通过配体-受体网络和靶基因调控矩阵，推断**从配体到下游靶基因的完整调控链条**。与 CellChat/Cellphonedb 不同，Nichenet 关注的是"哪个配体影响了哪些靶基因"，适合机制研究。

**核心原理**：结合差异基因分析，预测哪些配体具有调控活性，并推断其下游靶基因。

## 输入文件

| 文件 | 格式 | 说明 |
|------|------|------|
| Seurat 对象 | `.rds` | 包含差异基因分析结果 |
| `lr_network_human_allInfo.rds` | `.rds` | 配体-受体网络 |
| `ligand_target_matrix.rds` | `.rds` | 配体-靶基因调控矩阵 |

### 差异基因准备

```r
# 在 receiver 细胞中做差异分析
interesting_genes <- rownames(
  FindMarkers(seu, ident.1 = "treated", ident.2 = "control",
               group.by = "group_id", subset.ident = "Cancer_cell", assay = "RNA")
)
```

## 分析流程

1. **差异基因分析** → 在 receiver 细胞中找差异基因
2. **配体活性预测** → `predict_ligand_activities()`
3. **Circos 图** → 配体-受体-靶基因可视化
4. **综合图** → 配体-受体-活性-靶基因气泡图

## 输出结果

| 文件 | 内容 |
|------|------|
| `ligand_activities.csv` | 配体活性评分表 |
| `ligand_activities.rds` | 完整配体活性结果 |

## 适用场景

- ✅ 探究下游靶基因机制
- ✅ 配体→受体→靶基因完整链条
- ✅ 需要识别关键调控配体

## 注意事项

- 与 MultiNichenet 不同，Nichenet **需要生物学重复**
- 靶基因集应来自差异分析
- 适合 Receiver 细胞（如肿瘤细胞）受到哪些配体调控
