# scRNA Cell-Cell Communication Analysis

单细胞转录组细胞间通讯（Cell-Cell Communication, CCC）分析工具集。

## 📂 目录结构

```
scRNA-cell-communication/
├── README.md
└── references/
    ├── cellchat_v2/          # CellChat V2 分析
    ├── cellphonedb_v5/      # Cellphonedb V5 分析
    ├── nichenet/             # Nichenet 调控网络
    ├── multinichenet/         # MultiNichenet 少样本分析
    └── italk/                # iTALK 快速探索
```

每个子文件夹包含：
- `README.md` — 方法介绍、输入输出、适用场景
- `code.R` — 完整分析代码

## 🔬 六种方法概览

| 方法 | 输入数据 | 输出结果 | 适用场景 |
|------|---------|---------|---------|
| **CellChat V2** | Seurat 对象 / 表达矩阵 + 细胞注释 | 通讯网络、通路图、气泡图、弦图 | 文章发表首选，图形美观 |
| **Cellphonedb V5** | 归一化矩阵 + 元数据 | 显著互作对及 p 值 | 严格统计推断 |
| **Nichenet** | Seurat 对象 + 差异基因 | 配体-受体-靶基因调控链条 | 机制研究 |
| **MultiNichenet** | SingleCellExperiment 对象 | 优先级排序表、Circos 图 | **少样本/无重复**、多组比较 |
| **iTALK** | Seurat 对象 + 差异基因 | 网络可视化图 | 快速初步探索 |

## 🚀 快速选择

```
有生物学重复，文章发表 → CellChat V2
每组样本数 < 3，无重复 → MultiNichenet
严格统计检验，数据库全面 → Cellphonedb V5
探究配体下游靶基因机制 → Nichenet
初步快速探索，无需统计 → iTALK
```

## 📥 输入文件准备

### Seurat 对象

所有基于 R 的方法都需要一个注释好的 Seurat 对象：

```r
seu <- readRDS("your_annotated_scRNA.rds")
table(Idents(seu))  # 查看细胞类型
```

### MultiNichenet 特殊要求

需要转换为 SingleCellExperiment 对象：

```r
sce <- as.SingleCellExperiment(seu, assay = "RNA")
```

### Cellphonedb Python 要求

需要两个文件（制表符分隔）：
- **counts 文件**：基因名 × 细胞名矩阵（归一化后）
- **metadata 文件**：cell_id × cell_type 两列

## 📊 各方法对比

| 特性 | CellChat V2 | Cellphonedb V5 | MultiNichenet | Nichenet | iTALK |
|------|------------|----------------|----------------|----------|-------|
| **统计学方法** | 概率模型 | 置换检验 | 差异分析 | 调控网络推断 | 无统计 |
| **需要生物学重复** | ✅ | ✅ | ❌ | ✅ | ✅ |
| **多组比较** | ✅ | ❌ | ✅ | ❌ | ❌ |
| **配体→靶基因链条** | ❌ | ❌ | ✅ | ✅ | ❌ |
| **发表级别图形** | ✅ | ✅ | ✅ | ✅ | ❌ |

## ⚠️ 常见问题

**Q: CellChat 和 Cellphonedb 哪个更好？**
A: 两者原理不同，建议同时使用互为验证。CellChat 图形更美观，Cellphonedb 统计更严格。

**Q: 我的样本没有生物学重复怎么办？**
A: 使用 MultiNichenet，专为无重复或少样本（每组 < 3 个）设计。

**Q: 为什么找不到配体-受体对？**
A: 检查基因名是否为 Symbol 格式（不是 Ensemble ID）。
