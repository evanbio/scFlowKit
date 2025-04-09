# Rutils/find_conserved_markers.R
#-------------------------------------------------------------------------------
# scFlowKit: Find Conserved Marker Genes across Conditions
#-------------------------------------------------------------------------------
#
# find_conserved_markers: 基于 grouping.var 对所有聚类执行保守标志基因分析。
# - 使用 Seurat::FindConservedMarkers 分析每个聚类在不同分组间的一致表达标志基因。
# - 可指定输出目录和每个聚类保留的 top N 个基因。
# - 分析结果包含完整结果和 Top 标志基因列表，保存为 CSV 文件。
#
# 参数：
#   seu：Seurat 对象，包含 seurat_clusters 聚类信息。
#   grouping.var：元数据中的条件分组变量（例如 condition）。
#   output_dir：输出目录。
#   top_n：每个聚类保留的 top 基因数量，默认 5。
#   test.use：差异表达检验方法，默认 "MAST"。
#   only.pos：是否仅返回上调基因，默认 TRUE。
#   min.pct：基因在至少一组细胞中表达的最小比例，默认 0.25。
#   logfc.threshold：最小 log2FoldChange 阈值，默认 0.5。
# 返回：
#   返回一个 list，包括：
#     - conserved_markers：所有聚类的保守标志基因
#     - top_conserved_markers：每个聚类 top N 保守标志基因
#-------------------------------------------------------------------------------

find_conserved_markers <- function(
  seu,
  grouping.var,
  output_dir = "results",
  top_n = 5,
  test.use = "MAST",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5
) {
  #-------------------- 参数检查 --------------------
  cli::cli_h2("🧬 Step 4.3: 查找保守标志基因（FindConservedMarkers）")

  if (!inherits(seu, "Seurat")) {
    stop("seu 必须是 Seurat 对象！", call. = FALSE)
  }
  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    stop("meta.data 中缺少 seurat_clusters 列！", call. = FALSE)
  }
  if (!grouping.var %in% colnames(seu@meta.data)) {
    stop("meta.data 中缺少分组变量：", grouping.var, call. = FALSE)
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

  #-------------------- 展示分组情况 --------------------
  cli::cli_text("分组变量 ({grouping.var}) 的分布情况：")
  print(table(seu@meta.data[[grouping.var]]))

  #-------------------- 查找保守标志基因 --------------------
  clusters <- unique(seu@meta.data$seurat_clusters)
  cli::cli_alert_info("检测到 {length(clusters)} 个聚类：{paste(clusters, collapse = ', ')}")

  # 查看每个分组中聚类的组成
  split_tab <- table(seu@meta.data[[grouping.var]], seu@meta.data$seurat_clusters)
  cli::cli_text("各分组中聚类组成：")
  print(split_tab)

  # 遍历聚类进行保守基因分析
  markers_list <- list()

  for (cluster in clusters) {
    cli::cli_alert("分析聚类 {cluster} 的保守标志基因...")

    result <- FindConservedMarkers(
      seu,
      ident.1 = cluster, # 目标聚类
      grouping.var = grouping.var,
      test.use = test.use,
      only.pos = only.pos,
      min.pct = min.pct,
      logfc.threshold = logfc.threshold,
      verbose = TRUE
    )

    if (nrow(result) > 0) {
      result <- tibble::rownames_to_column(result, var = "gene")
      result$cluster <- cluster
      markers_list[[as.character(cluster)]] <- result
    }
  }

  #-------------------- 合并结果并提取 topN --------------------
  conserved_markers <- dplyr::bind_rows(markers_list)
  log2fc_cols <- grep("_avg_log2FC$", names(conserved_markers), value = TRUE)

  conserved_markers <- conserved_markers %>%
    dplyr::mutate(mean_log2FC = rowMeans(dplyr::select(., dplyr::all_of(log2fc_cols))))

  top_conserved_markers <- conserved_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc(mean_log2FC)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()

  #-------------------- 保存结果 --------------------
  all_file <- file.path(output_dir, "scFlowKit_conserved_markers.csv")
  top_file <- file.path(output_dir, sprintf("scFlowKit_top%d_conserved_markers.csv", top_n))

  write.csv(conserved_markers, all_file, row.names = FALSE)
  write.csv(top_conserved_markers, top_file, row.names = FALSE)

  cli::cli_alert_success("已保存完整保守标志基因结果至：{all_file}")
  cli::cli_alert_success("已保存 Top {top_n} 保守标志基因至：{top_file}")

  return(list(
    conserved_markers = conserved_markers,
    top_conserved_markers = top_conserved_markers
  ))
}

#-------------------------------------------------------------------------------
