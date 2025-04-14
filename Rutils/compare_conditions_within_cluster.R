# Rutils/compare_conditions_within_cluster.R
#-------------------------------------------------------------------------------
# scFlowKit: Compare Differential Gene Expression Between Conditions Within Each Cluster
#-------------------------------------------------------------------------------
#
# 背景说明：
# 本函数在每个指定的聚类（cluster）中，比较不同条件（例如 "WT" vs "KO"）之间的
# 差异表达基因（DEG），用于分析条件对单细胞聚类的特异性影响。
#
# 策略说明：
# - 对每个聚类子集 Seurat 对象；
# - 使用 FindMarkers 比较指定条件（ident.1 vs ident.2）；
# - 返回完整差异基因表（full_table）和 Top N 基因表（top_table）；
# - 输出 CSV 文件，并记录日志。
#
# 参数说明：
# - seu           : Seurat 对象，包含表达矩阵和元数据
# - group.by      : 条件分组变量（meta.data 中的列名，例如 "condition"）
# - ident.1       : 条件 1（字符向量，例如 "WT"）
# - ident.2       : 条件 2（字符向量，例如 "KO"）
# - cluster_col   : 聚类列（默认 "seurat_clusters"）
# - cluster_id    : 指定的聚类 ID 向量（默认 NULL，表示使用所有聚类）
# - output_dir    : 输出目录（默认 "results/tables"）
# - test.use      : 差异分析方法（默认 a"MAST"）
# - only.pos      : 是否只返回上调基因（默认 TRUE）
# - min.pct       : 最小表达比例（默认 0.25）
# - logfc.threshold: 最小 log2FC（默认 0.5）
# - top_n         : Top N 基因输出（默认 5）
# - log           : 是否保存日志（默认 TRUE）
#
# 返回值：
# - list 包含：
#   - full_table : 完整差异基因表，包含 cluster, gene, avg_log2FC, p_val_adj, comparison
#   - top_table  : Top N 基因表，按 avg_log2FC 排序
#
# 额外输出：
# - CSV 文件：
#   - results/tables/compare_{ident.1}_vs_{ident.2}_bycluster_full.csv
#   - results/tables/compare_{ident.1}_vs_{ident.2}_bycluster_top{top_n}.csv
# - 日志文件：logs/find_markers/compare_conditions_within_cluster.log
#
#-------------------------------------------------------------------------------

library(tidyverse)

compare_conditions_within_cluster <- function(
  seu,
  group.by,
  ident.1,
  ident.2,
  cluster_col = "seurat_clusters",
  cluster_id = NULL,
  output_dir = "results/tables",
  test.use = "MAST",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5,
  top_n = 5,
  log = TRUE
) {
  #-------------------- 依赖检查 --------------------
  for (pkg in c("tidyverse", "cli", "Seurat")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("请先安装 ", pkg, " 包", call. = FALSE)
    }
  }

  #-------------------- 日志设置 --------------------
  log_file <- "logs/find_markers/compare_conditions_within_cluster.log"
  if (log) {
    dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
    sink(log_file, append = TRUE, split = TRUE)
    cat("\n", strrep("-", 70), "\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Start of compare_conditions_within_cluster\n")
    cat("• Group by         : ", group.by, "\n")
    cat("• Ident 1          : ", ident.1, "\n")
    cat("• Ident 2          : ", ident.2, "\n")
    cat("• Cluster column   : ", cluster_col, "\n")
    cat("• Cluster IDs      : ", if (is.null(cluster_id)) "All" else paste(cluster_id, collapse = ", "), "\n")
    cat("• Test method      : ", test.use, "\n")
    cat("• Only positive    : ", only.pos, "\n")
    cat("• Min pct          : ", min.pct, "\n")
    cat("• logFC threshold  : ", logfc.threshold, "\n")
    cat("• Top N            : ", top_n, "\n")
    cat(strrep("-", 70), "\n")
  }

  #-------------------- 参数校验 --------------------
  if (!inherits(seu, "Seurat")) {
    stop("seu 必须为 Seurat 对象", call. = FALSE)
  }
  if (!group.by %in% colnames(seu@meta.data)) {
    stop("未找到分组列：", group.by, call. = FALSE)
  }
  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop("未找到聚类列：", cluster_col, call. = FALSE)
  }
  if (!is.character(ident.1) || length(ident.1) != 1) {
    stop("ident.1 必须为单一字符值", call. = FALSE)
  }
  if (!is.character(ident.2) || length(ident.2) != 1) {
    stop("ident.2 必须为单一字符值", call. = FALSE)
  }
  if (!is.character(test.use) || length(test.use) != 1) {
    stop("test.use 必须为单一字符值", call. = FALSE)
  }
  if (!is.logical(only.pos) || length(only.pos) != 1) {
    stop("only.pos 必须为单一逻辑值", call. = FALSE)
  }
  if (!is.numeric(min.pct) || min.pct < 0 || min.pct > 1) {
    stop("min.pct 必须为 0 到 1 之间的数值", call. = FALSE)
  }
  if (!is.numeric(logfc.threshold) || logfc.threshold < 0) {
    stop("logfc.threshold 必须为非负数", call. = FALSE)
  }
  if (!is.numeric(top_n) || top_n < 1 || top_n != round(top_n)) {
    stop("top_n 必须为正整数", call. = FALSE)
  }
  if (!is.logical(log) || length(log) != 1) {
    stop("log 必须为单一逻辑值", call. = FALSE)
  }

  #-------------------- 设置默认 assay --------------------
  DefaultAssay(seu) <- "RNA"
  cli::cli_h1("开始比较 {ident.1} vs {ident.2} 在每个聚类的差异表达基因")

  #-------------------- 获取聚类 ID --------------------
  if (is.null(cluster_id)) {
    cluster_id <- unique(as.character(seu@meta.data[[cluster_col]]))
  } else {
    cluster_id <- as.character(cluster_id)
    invalid_clusters <- setdiff(cluster_id, unique(seu@meta.data[[cluster_col]]))
    if (length(invalid_clusters) > 0) {
      stop("以下 cluster_id 不在 ", cluster_col, " 中：", paste(invalid_clusters, collapse = ", "), call. = FALSE)
    }
  }
  cli::cli_alert_info("处理 {length(cluster_id)} 个聚类：{paste(cluster_id, collapse = ', ')}")

  #-------------------- 差异分析 --------------------
  all_results <- list()
  top_results <- list()

  for (cl in cluster_id) {
    cli::cli_h2("处理聚类 {cl}")

    # 子集细胞
    subset_cells <- seu@meta.data[[cluster_col]] == cl
    if (sum(subset_cells) == 0) {
      cli::cli_alert_warning("聚类 {cl} 中无细胞，跳过")
      if (log) cat("警告：聚类 ", cl, " 中无细胞，跳过\n")
      next
    }

    seu_sub <- subset(seu, cells = WhichCells(seu)[subset_cells])
    Idents(seu_sub) <- factor(seu_sub@meta.data[[group.by]], levels = unique(seu_sub@meta.data[[group.by]]))

    # 检查条件是否存在
    if (!all(c(ident.1, ident.2) %in% levels(Idents(seu_sub)))) {
      cli::cli_alert_warning("聚类 {cl} 中不包含组别 {ident.1} 或 {ident.2}，跳过")
      if (log) cat("警告：聚类 ", cl, " 中不包含组别 ", ident.1, " 或 ", ident.2, "\n")
      next
    }

    # 差异分析
    markers <- tryCatch(
      {
        Seurat::FindMarkers(
          object = seu_sub,
          ident.1 = ident.1,
          ident.2 = ident.2,
          test.use = test.use,
          only.pos = only.pos,
          min.pct = min.pct,
          logfc.threshold = logfc.threshold,
          verbose = FALSE
        )
      },
      error = function(e) {
        cli::cli_alert_warning("聚类 {cl} 差异分析失败：{e$message}")
        if (log) cat("错误：聚类 ", cl, " 差异分析失败 - ", e$message, "\n")
        return(NULL)
      }
    )

    if (is.null(markers) || nrow(markers) == 0) {
      cli::cli_alert_warning("聚类 {cl} 无差异基因，跳过")
      if (log) cat("警告：聚类 ", cl, " 无差异基因\n")
      next
    }

    # 整理结果
    markers <- markers %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::mutate(
        cluster = cl,
        comparison = paste0(ident.1, "_vs_", ident.2)
      ) %>%
      dplyr::select(cluster, comparison, gene, avg_log2FC, p_val_adj, dplyr::everything())

    all_results[[cl]] <- markers
    top_results[[cl]] <- markers %>%
      dplyr::arrange(desc(avg_log2FC)) %>%
      dplyr::slice_head(n = top_n)
  }

  #-------------------- 合并结果 --------------------
  full_table <- dplyr::bind_rows(all_results)
  top_table <- dplyr::bind_rows(top_results)

  if (nrow(full_table) == 0) {
    cli::cli_alert_warning("无差异基因，请检查数据或参数")
    if (log) {
      cat("错误：无差异基因\n")
      sink()
    }
    stop("无差异基因，请调整参数", call. = FALSE)
  }

  #-------------------- 保存输出 --------------------
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  filename_prefix <- paste0("compare_", ident.1, "_vs_", ident.2, "_bycluster")
  full_file <- file.path(output_dir, paste0(filename_prefix, "_full.csv"))
  top_file <- file.path(output_dir, paste0(filename_prefix, "_top", top_n, ".csv"))

  write.csv(full_table, full_file, row.names = FALSE)
  write.csv(top_table, top_file, row.names = FALSE)
  cli::cli_alert_success("已保存结果至：")
  cli::cli_ul(c(full_file, top_file))

  #-------------------- 打印结果 --------------------
  cli::cli_h1("差异分析结果")
  cli::cli_alert_info("共发现 {nrow(full_table)} 个差异基因，覆盖 {length(unique(full_table$cluster))} 个聚类")
  cli::cli_alert_info("Top {top_n} 基因（前几行）：")
  print(head(top_table))

  if (log) {
    cat("共处理聚类：", length(unique(full_table$cluster)), "\n")
    cat("差异基因总数：", nrow(full_table), "\n")
    cat("Top ", top_n, " 基因（前几行）：\n")
    print(head(top_table))
    cat(strrep("-", 70), "\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - End of compare_conditions_within_cluster\n")
    sink()
    cli::cli_alert_success("日志已保存至：{log_file}")
  }

  return(list(
    full_table = full_table,
    top_table = top_table
  ))
}