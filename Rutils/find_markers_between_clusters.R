# Rutils/find_markers_between_clusters.R
#-------------------------------------------------------------------------------
# scFlowKit: Compare Marker Genes Between Custom Cluster Groups
#-------------------------------------------------------------------------------
#
# find_markers_between_clusters: 比较任意两个聚类（或聚类组）之间的差异表达基因。
# - 支持多聚类分组（如 ident.1 = c("0", "1") vs ident.2 = c("2", "3")）。
# - 自动执行 Seurat::FindMarkers，并输出完整结果和 Top N 基因。
# - 结果保存为 CSV 文件。
#
# 参数：
#   seu：Seurat 对象，包含 seurat_clusters。
#   ident.1, ident.2：要比较的两个聚类（或聚类组），向量。
#   output_dir：输出目录。
#   top_n：保留的 top 差异基因数，默认 5。
#   test.use：差异分析方法（默认 MAST）。
#   only.pos：是否仅返回上调基因。
#   min.pct：最低表达比例。
#   logfc.threshold：最小 log2FC 阈值。
# 返回：
#   返回一个 list，包括：
#     - markers：全部差异基因
#     - top_markers：top N 基因
#-------------------------------------------------------------------------------

find_markers_between_clusters <- function(
  seu,
  ident.1,
  ident.2,
  output_dir = "results",
  top_n = 5,
  test.use = "MAST",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5
) {
  #-------------------- 参数检查 --------------------
  if (!inherits(seu, "Seurat")) {
    stop("seu 必须是 Seurat 对象！", call. = FALSE)
  }
  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    stop("Seurat 对象缺少 seurat_clusters！", call. = FALSE)
  }
  if (!is.character(ident.1) || length(ident.1) < 1) {
    stop("ident.1 必须为字符向量且至少包含一个聚类标签！", call. = FALSE)
  }
  if (!is.character(ident.2) || length(ident.2) < 1) {
    stop("ident.2 必须为字符向量且至少包含一个聚类标签！", call. = FALSE)
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
  if (!is.character(output_dir)) {
    stop("output_dir 应为字符型路径。", call. = FALSE)
  }
  if (!is.numeric(top_n) || top_n <= 0 || top_n != as.integer(top_n)) {
    stop("top_n 应为正整数。", call. = FALSE)
  }
  
  #-------------------- 差异分析 --------------------
  cli::cli_alert_info("比较聚类组：{paste(ident.1, collapse = ", ")} vs {paste(ident.2, collapse = ", ")}")

  markers <- Seurat::FindMarkers(
    object = seu,
    ident.1 = ident.1,
    ident.2 = ident.2,
    subset.ident = c(ident.1, ident.2),  # 限制分析的聚类
    test.use = test.use,
    only.pos = only.pos,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    verbose = TRUE
  )
  # 将行名（基因名）转换为 gene 列
  markers <- tibble::rownames_to_column(markers, var = "gene")

  #-------------------- 提取 Top N --------------------
  top_markers <- markers %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    dplyr::slice_head(n = top_n)

  #-------------------- 保存结果 --------------------
  group_label <- paste0("cluster_", 
                        paste(ident.1, collapse = "_"), 
                        "_vs_", 
                        paste(ident.2, collapse = "_"))

  all_file <- file.path(output_dir, 
                        glue::glue("scFlowKit_markers_{group_label}.csv"))
  top_file <- file.path(output_dir, 
                        glue::glue("scFlowKit_top{top_n}_markers_{group_label}.csv"))

  write.csv(markers, all_file, row.names = FALSE)
  write.csv(top_markers, top_file, row.names = FALSE)

  cli::cli_alert_success("差异分析结果已保存至：{all_file}")
  cli::cli_alert_success("Top {top_n} 基因结果已保存至：{top_file}")

  return(list(
    markers = markers,
    top_markers = top_markers
  ))
}

#-------------------------------------------------------------------------------
