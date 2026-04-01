# MultiNichenet 分析指南

## 方法简介

MultiNichenet 是专为**少样本或无生物学重复**设计的细胞通讯分析方法。适用于临床样本（每组样本数 < 3 个）。通过差异表达分析结合配体-受体网络，推断优先级最高的互作对。

**核心原理**：对 sender 和 receiver 细胞分别做差异分析，结合配体活性预测和靶基因调控网络，计算每个配体-受体对的综合评分。

## 输入文件

| 文件 | 格式 | 说明 |
|------|------|------|
| SingleCellExperiment 对象 | `.rds` | 包含 `RNA` assay，需转换自 Seurat |
| `lr_network_human_allInfo.rds` | `.rds` | 配体-受体网络（需从 MultiNichenet 官网下载）|
| `ligand_target_matrix_nsga2r_final.rds` | `.rds` | 配体-靶基因调控矩阵（需从 MultiNichenet 官网下载）|

### Seurat → SingleCellExperiment 转换

```r
library(SingleCellExperiment)
library(Seurat)

seu <- readRDS("your_scRNA.rds")
sce <- as.SingleCellExperiment(seu, assay = "RNA")
```

### 元数据要求

colData 中需要包含：
- `orig.ident`：样本 ID
- `group`：分组（如 Disease / Normal）
- `celltype`：细胞类型

## 网络文件下载

需从 MultiNichenet GitHub 下载 human 或 mouse 网络文件：
```r
# human network files
lr_network_human_allInfo.rds
ligand_target_matrix_nsga2r_final.rds
```

## 分析流程

1. **数据转换** → Seurat → SingleCellExperiment
2. **过滤** → 细胞类型过滤 + 基因过滤
3. **表达量处理** → `process_abundance_expression_info()`
4. **差异分析** → `get_DE_info_sampleAgnostic()`（无重复方法）
5. **配体活性预测** → `get_ligand_activities_targets_DEgenes()`
6. **优先级排序** → `generate_prioritization_tables()`
7. **可视化** → Circos 图 + 综合图

## 输出结果

| 文件 | 内容 |
|------|------|
| `multinichenet_output.rds` | 完整分析结果列表 |
| `group_prioritization_tbl` | 每组的优先级排序表 |
| `celltype_de` | 细胞类型差异基因 |

**优先级表列说明**：
- `sender` / `receiver`：通讯细胞类型
- `ligand` / `receptor`：配体-受体对
- `prioritization_score`：综合评分（越高越显著）

## 适用场景

- ✅ **无生物学重复**（最大优势）
- ✅ 每组样本数 < 3 个
- ✅ 多组比较（如疾病 vs 对照 vs 处理）
- ✅ 需要同时考虑配体和受体差异

## 注意事项

- 需要准备网络文件（从官网下载）
- 基因过滤 `fraction_cutoff` 默认 0.05，过严时可下调
- sender 和 receiver 细胞类型相同时也能分析
- `scenario = "no_frac_LR_expr"` 用于无重复分析
