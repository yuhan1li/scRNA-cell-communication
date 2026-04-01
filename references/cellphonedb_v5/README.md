# Cellphonedb V5 分析指南

## 方法简介

Cellphonedb V5 是 Cellphonedb 的升级版本，使用更严格的统计方法（置换检验）来判断配体-受体对是否显著通讯。数据库涵盖人/小鼠的蛋白相互作用，支持复合物（complex）形式的受配体。

**核心原理**：对配体和受体在细胞类型对之间进行置换检验，判断互作是否显著。

## 输入文件

Cellphonedb 需要 **Python 端运行**，R 端负责数据准备。

### R 端：导出数据

| 文件 | 格式 | 说明 |
|------|------|------|
| `counts.txt` | 制表符分隔 | 基因名 × 细胞名，表达值矩阵（归一化后）|
| `metadata.txt` | 制表符分隔 | 两列：`Cell`（细胞名）+ `cell_type`（细胞类型）|

**注意事项**：
- 细胞类型名称中的**空格**会替换为**下划线**
- 行名不包含在矩阵中，第一列为基因名
- 使用 `layer = 'data'`（归一化后数据）而非原始 counts

### Python 端：运行分析

```bash
# 有生物学重复
cellphonedb statics counts.txt metadata.txt \
  --counts-data hgnc_symbol \
  --output-path output_v5 \
  --threads 8 \
  --iterations 1000 \
  --pvalues-method permtest

# 无生物学重复（单样本）
cellphonedb statics-single counts.txt metadata.txt \
  --counts-data hgnc_symbol \
  --output-path output_single \
  --threads 8
```

## 输出文件

| 文件 | 内容 |
|------|------|
| `statistical_analysis_significant_means_*.txt` | 显著互作的平均表达值 |
| `statistical_analysis_pvalues_*.txt` | 置换检验 p 值 |

**结果解读**：
- 列名格式 `Celltype_A|Celltype_B` 表示 A→B 的通讯
- `rank` 列表示显著性排名
- `classification` 列表示信号通路分类

## 适用场景

- ✅ 需要严格统计检验
- ✅ 数据库覆盖全面（包含复合物）
- ✅ 细胞类型较多
- ✅ 需要与 CellChat 互相验证

## 注意事项

- Cellphonedb **不支持多组比较**，只能分析单样本或合并样本
- 元数据中的细胞类型名称需与表达矩阵列名一致
- 建议在 Python 虚拟环境中安装 cellphonedb，避免包冲突
