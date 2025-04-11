# Rutils/find_all_markers.R
#-------------------------------------------------------------------------------
# scFlowKit: Find Marker Genes for All Clusters
#-------------------------------------------------------------------------------
#
# find_all_markers: 对所有聚类执行差异表达分析，识别标志基因。
# - 自动切换 RNA assay，进行归一化、变量基因选择与缩放。
# - 使用 Seurat::FindAllMarkers，输出完整结果和 topN 基因列表。
# - 结果保存为 CSV 文件（完整结果和 top 基因）
#
# 参数：
#   seu：Seurat 对象，包含 seurat_clusters 聚类信息。
#   output_dir：输出目录，用于保存 CSV 文件。
#   top_n：每个聚类保留的 top 基因数量，默认 5。
#   test.use：差异表达检验方法，默认 "MAST"。
#   only.pos：是否仅返回上调基因，默认 TRUE。
#   min.pct：基因在至少一组细胞中表达的最小比例，默认 0.25。
#   logfc.threshold：最小 log2FoldChange 阈值，默认 0.5。
# 返回：
#   返回一个 list，包括：
#     - all_markers：全部标志基因结果数据框
#     - top_markers：top N 基因子集
#-------------------------------------------------------------------------------

find_all_markers <- function(
  seu,
  output_dir = "results",
  top_n = 5,
  test.use = "MAST",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5
) {

  # -------------------- 参数检查 --------------------
  if (!inherits(seu, "Seurat")) {
    stop("seu 必须为 Seurat 对象！", call. = FALSE)
  }
  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    stop("Seurat 对象缺少 seurat_clusters 聚类信息！", call. = FALSE)
  }
  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir 必须为单一字符！", call. = FALSE)
  }
  if (!is.numeric(top_n) || top_n <= 0 || top_n != as.integer(top_n)) {
    stop("top_n 必须为正整数！", call. = FALSE)
  }
  if (!is.character(test.use)) {
    stop("test.use 必须为字符！", call. = FALSE)
  }
  if (!is.logical(only.pos)) {
    stop("only.pos 必须为逻辑值！", call. = FALSE)
  }
  if (!is.numeric(min.pct) || min.pct < 0 || min.pct > 1) {
    stop("min.pct 应为 0 到 1 之间的数值！", call. = FALSE)
  }
  if (!is.numeric(logfc.threshold)) {
    stop("logfc.threshold 应为数值型！", call. = FALSE)
  }

  # -------------------- 处理 RNA assay --------------------
  DefaultAssay(seu) <- "RNA"

  # -------------------- 查找标志基因 --------------------
  all_markers <- FindAllMarkers(
    object = seu,
    test.use = test.use,
    only.pos = only.pos,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    verbose = TRUE
  )

  # -------------------- 保存结果 --------------------
  out_all_file <- file.path(output_dir, "scFlowKit_cluster_markers.csv")
  out_top_file <- file.path(output_dir, sprintf("scFlowKit_top%d_markers.csv", top_n))

  write.csv(all_markers, out_all_file, row.names = FALSE)

  top_markers <- all_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()

  write.csv(top_markers, out_top_file, row.names = FALSE)

  return(list(
    all_markers = all_markers,
    top_markers = top_markers
  ))
}

#-------------------------------------------------------------------------------
