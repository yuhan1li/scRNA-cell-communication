# MultiNichenet 无重复样本分析完整代码
# 来源: KS-codes(2025) / 1-scRNA 细胞通讯分析 / multinichentr（视频教程）
# 适用: 没有样本重复或每组样本数少于3个的数据

library(SingleCellExperiment)
library(multinichenetr)
library(Seurat)
library(tidyverse)

# ============================================================================
# 1. 数据加载与转换
# ============================================================================
load("scRNA_nichenet.RData")  # 示例 Seurat 对象

# Seurat → SingleCellExperiment 转换
sce <- Seurat::as.SingleCellExperiment(sce, assay = "RNA")
sce <- alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()

# ============================================================================
# 2. 加载网络文件（需从 MultiNichenet 官网下载）
# ============================================================================
lr_network_all <- readRDS("lr_network_human_allInfo.rds") %>%
  mutate(ligand = convert_alias_to_symbols(ligand, organism = "human"),
         receptor = convert_alias_to_symbols(receptor, organism = "human")) %>%
  mutate(ligand = make.names(ligand), receptor = make.names(receptor))
lr_network <- lr_network_all %>% distinct(ligand, receptor)

ligand_target_matrix <- readRDS("ligand_target_matrix_nsga2r_final.rds")
colnames(ligand_target_matrix) <- colnames(ligand_target_matrix) %>%
  convert_alias_to_symbols(organism = "human") %>% make.names()
rownames(ligand_target_matrix) <- rownames(ligand_target_matrix) %>%
  convert_alias_to_symbols(organism = "human") %>% make.names()

# 只保留在配体-靶基因矩阵中的配体
lr_network <- lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
ligand_target_matrix <- ligand_target_matrix[, lr_network$ligand %>% unique()]

# ============================================================================
# 3. 设置分组参数
# ============================================================================
sample_id <- "orig.ident"
group_id <- "group"
celltype_id <- "celltype"
batches <- NA

# 整理 colData
SummarizedExperiment::colData(sce)$orig.ident <-
  SummarizedExperiment::colData(sce)$orig.ident %>% make.names()
SummarizedExperiment::colData(sce)$group <-
  SummarizedExperiment::colData(sce)$group %>% make.names()
SummarizedExperiment::colData(sce)$celltype <-
  SummarizedExperiment::colData(sce)$celltype %>% make.names()

# 定义比较组（如 SD vs HC）
contrasts_oi <- c("'SD-HC','HC-SD'")
contrast_tbl <- tibble(contrast = c("SD-HC", "HC-SD"), group = c("SD", "HC"))

# 定义 sender 和 receiver（通常是全部细胞类型）
senders_oi <- SummarizedExperiment::colData(sce)[, celltype_id] %>% unique()
receivers_oi <- SummarizedExperiment::colData(sce)[, celltype_id] %>% unique()

# ============================================================================
# 4. Step 1: 细胞类型过滤
# ============================================================================
min_cells <- 10
abundance_info <- get_abundance_info(
  sce = sce, sample_id = sample_id, group_id = group_id,
  celltype_id = celltype_id, min_cells = min_cells,
  senders_oi = senders_oi, receivers_oi = receivers_oi, batches = batches)

# ============================================================================
# 5. Step 2: 基因过滤
# ============================================================================
fraction_cutoff <- 0.05  # 基因在多少比例的样本中表达
min_sample_prop <- 1

frq_list <- get_frac_exprs_sampleAgnostic(
  sce = sce, sample_id = sample_id, celltype_id = celltype_id,
  group_id = group_id, batches = batches, min_cells = min_cells,
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)

genes_oi <- frq_list$expressed_df %>% filter(expressed == TRUE) %>% pull(gene) %>% unique()
sce <- sce[genes_oi, ]

# ============================================================================
# 6. Step 3: 表达量处理
# ============================================================================
abundance_expression_info <- process_abundance_expression_info(
  sce = sce, sample_id = group_id, group_id = group_id,
  celltype_id = celltype_id, min_cells = min_cells,
  senders_oi = senders_oi, receivers_oi = receivers_oi,
  lr_network = lr_network, batches = batches, frq_list = frq_list,
  abundance_info = abundance_info)

# ============================================================================
# 7. Step 4: 差异分析（无重复）
# ============================================================================
DE_info <- get_DE_info_sampleAgnostic(
  sce = sce, group_id = group_id, celltype_id = celltype_id,
  contrasts_oi = contrasts_oi, expressed_df = frq_list$expressed_df,
  min_cells = min_cells, contrast_tbl = contrast_tbl)

celltype_de <- DE_info$celltype_de_findmarkers

# 合并 sender 和 receiver 差异基因
sender_receiver_de <- multinichenetr::combine_sender_receiver_de(
  sender_de = celltype_de, receiver_de = celltype_de,
  senders_oi = senders_oi, receivers_oi = receivers_oi, lr_network = lr_network)

# ============================================================================
# 8. Step 5: 配体活性预测
# ============================================================================
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
    p_val_adj = p_val_adj, top_n_target = top_n_target,
    verbose = F, n.cores = 2
  )
))

# ============================================================================
# 9. Step 6: 优先级排序
# ============================================================================
metadata_combined <- SummarizedExperiment::colData(sce) %>% tibble::as_tibble()
grouping_tbl <- metadata_combined[, c(group_id)] %>% tibble::as_tibble() %>%
  distinct() %>% rename(group = group_id) %>% mutate(sample = group)

prioritization_tables <- suppressMessages(multinichenetr::generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_de %>% distinct(sender, receiver),
  grouping_tbl = grouping_tbl,
  scenario = "no_frac_LR_expr",
  fraction_cutoff = fraction_cutoff,
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = FALSE
))

# ============================================================================
# 10. 保存结果
# ============================================================================
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

# ============================================================================
# 11. 可视化
# ============================================================================
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
