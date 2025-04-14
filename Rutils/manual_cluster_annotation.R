# Rutils/manual_cluster_annotation.R
#-------------------------------------------------------------------------------
# scFlowKit: Suggest Cluster Annotations by Matching Top DEGs to Marker Set
#-------------------------------------------------------------------------------
#
# 背景说明：
# 本函数用于将每个 cluster 的 Top N 差异表达基因与标准 marker set 匹配，
# 根据交集数量推荐最可能的细胞类型标签，用于手动注释参考。
#
# 策略说明：
# - 对 DEG 数据筛选 avg_log2FC、pct.1、pct.2 满足阈值的基因；
# - 在筛选后，取每个 cluster 的 top N 基因（按 avg_log2FC 排序）；
# - 统计这些 Top N 与 marker set 的匹配数，得出推荐标签；
# - 返回标签向量、命中矩阵、汇总表，包含 second_best 和比例等辅助列。
#
# 参数说明：
# - deg            : FindAllMarkers 风格数据框，需包含 cluster, gene, avg_log2FC, pct.1, pct.2
# - marker_set     : 命名 list，形如 list("B cells" = c("CD79A", "MS4A1"))
# - set_name       : marker set 的名称，用于命名 CSV 输出（默认 "default"）
# - top_n          : 每个 cluster 取 Top N 基因（默认 5）
# - logfc_cutoff   : 最小 avg_log2FC（默认 0.25）
# - pct1_cutoff    : 最小 pct.1（默认 0.5）
# - pct2_cutoff    : 最大 pct.2（默认 0.2）
# - default        : 若无命中，填充默认标签（默认 "ambiguous"）
# - log            : 是否记录日志（默认 TRUE）
#
# 返回值：
# - list 包含：
#   - result_table  : cluster, best_match, second_best, best_hit_count, marker_set_size, all_hit_count, best_hit_ratio
#   - result_vector : 推荐标签（命名向量）
#   - anno_matrix   : 矩阵 [cluster × celltype]，值为交集基因数
#
# 额外输出：
# - 全局变量 marker_match_result 保存结果
# - CSV 文件保存 result_table，路径为 results/tables/cell_annotation/manual_cluster_annotation_{set_name}.csv
# - 日志记录运行信息，路径为 logs/cell_annotation/manual_cluster_annotation.log
#
#-------------------------------------------------------------------------------

library(tidyverse)

manual_cluster_annotation <- function(
    deg,
    marker_set,
    set_name = "default",
    top_n = 5,
    logfc_cutoff = 0.25,
    pct1_cutoff = 0.5,
    pct2_cutoff = 0.2,
    default = "ambiguous",
    log = TRUE
) {
  #-------------------- 依赖检查 --------------------
  for (pkg in c("tidyverse", "cli")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("请先安装 ", pkg, " 包", call. = FALSE)
    }
  }
  
  #-------------------- 日志设置 --------------------
  log_file <- "logs/cell_annotation/manual_cluster_annotation.log"
  if (log) {
    dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
    sink(log_file, append = TRUE, split = TRUE)
    cat("\n", strrep("-", 70), "\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Start of manual_cluster_annotation\n")
    cat("• Input clusters    : ", length(unique(deg$cluster)), "\n")
    cat("• Input genes       : ", length(unique(deg$gene)), "\n")
    cat("• Marker celltypes  : ", length(marker_set), "\n")
    cat("• Marker set name   : ", set_name, "\n")
    cat("• Top N             : ", top_n, "\n")
    cat("• logFC cutoff      : ", logfc_cutoff, "\n")
    cat("• pct.1 cutoff      : ", pct1_cutoff, "\n")
    cat("• pct.2 cutoff      : ", pct2_cutoff, "\n")
    cat("• Default label     : ", default, "\n")
    cat(strrep("-", 70), "\n")
  }
  
  #-------------------- 参数校验 --------------------
  required_cols <- c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2")
  if (!all(required_cols %in% colnames(deg))) {
    stop("DEG 数据必须包含以下列：", paste(required_cols, collapse = ", "), call. = FALSE)
  }
  if (!is.list(marker_set) || is.null(names(marker_set)) || length(marker_set) == 0) {
    stop("marker_set 必须为非空的命名 list", call. = FALSE)
  }
  for (ct in names(marker_set)) {
    if (!is.character(marker_set[[ct]]) || length(marker_set[[ct]]) == 0) {
      stop("marker_set 中的 ", ct, " 必须为非空字符向量", call. = FALSE)
    }
  }
  if (!is.character(set_name) || length(set_name) != 1) {
    stop("set_name 必须为单一字符值", call. = FALSE)
  }
  if (!is.numeric(top_n) || top_n < 1 || top_n != round(top_n)) {
    stop("top_n 必须为正整数", call. = FALSE)
  }
  if (!is.numeric(logfc_cutoff) || logfc_cutoff < 0) {
    stop("logfc_cutoff 必须为非负数", call. = FALSE)
  }
  if (!is.numeric(pct1_cutoff) || pct1_cutoff < 0 || pct1_cutoff > 1) {
    stop("pct1_cutoff 必须为 0 到 1 之间的数值", call. = FALSE)
  }
  if (!is.numeric(pct2_cutoff) || pct2_cutoff < 0 || pct2_cutoff > 1) {
    stop("pct2_cutoff 必须为 0 到 1 之间的数值", call. = FALSE)
  }
  if (!is.character(default) || length(default) != 1) {
    stop("default 必须为单一字符值", call. = FALSE)
  }
  if (!is.logical(log) || length(log) != 1) {
    stop("log 必须为单一逻辑值", call. = FALSE)
  }
  
  #-------------------- 筛选 DEG 数据 --------------------
  cli::cli_h1("筛选 Top {top_n} 差异表达基因")
  
  deg_filtered <- deg %>%
    dplyr::filter(
      avg_log2FC >= logfc_cutoff,
      pct.1 >= pct1_cutoff,
      pct.2 <= pct2_cutoff
    ) %>%
    dplyr::arrange(cluster, desc(avg_log2FC)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()
  
  if (nrow(deg_filtered) == 0) {
    cli::cli_alert_warning("过滤后无基因满足条件，请检查 logfc_cutoff, pct1_cutoff, pct2_cutoff")
    if (log) {
      cat("错误：过滤后无基因满足条件\n")
      sink()
    }
    stop("过滤后无基因，请调整阈值", call. = FALSE)
  }
  
  deg$cluster <- as.character(deg$cluster)
  
  clusters <- unique(deg_filtered$cluster)
  if (length(clusters) == 0) {
    cli::cli_alert_warning("无聚类包含满足条件的基因")
    if (log) {
      cat("错误：无聚类包含满足条件的基因\n")
      sink()
    }
    stop("无聚类包含基因，请调整阈值", call. = FALSE)
  }
  celltypes <- names(marker_set)
  
  cli::cli_alert_info("过滤后保留 {nrow(deg_filtered)} 个基因，覆盖 {length(clusters)} 个聚类")
  
  #-------------------- 交集匹配计算 --------------------
  cli::cli_h1("匹配 Top 基因与 Marker Set")
  
  match_matrix <- matrix(0, nrow = length(clusters), ncol = length(celltypes),
                         dimnames = list(clusters, celltypes))
  
  for (cl in clusters) {
    genes <- deg_filtered$gene[deg_filtered$cluster == cl]
    for (ct in celltypes) {
      hits <- intersect(genes, marker_set[[ct]])
      match_matrix[cl, ct] <- length(hits)
    }
  }
  
  #-------------------- 推荐标签计算 --------------------
  result_vector <- apply(match_matrix, 1, function(row) {
    if (max(row) == 0) return(default)
    names(row)[which.max(row)]
  })
  names(result_vector) <- clusters
  
  #-------------------- 汇总结果生成 --------------------
  result_table <- tibble(
    cluster = rownames(match_matrix),
    best_match = result_vector,
    best_hit_count = apply(match_matrix, 1, max),
    second_best = apply(match_matrix, 1, function(x) {
      sorted <- sort(x, decreasing = TRUE)
      if (length(sorted) >= 2 && sorted[2] > 0) names(sorted)[2] else default
    }),
    marker_set_size = sapply(result_vector, function(ct) {
      if (ct == default) 0 else length(marker_set[[ct]])
    }),
    all_hit_count = apply(match_matrix, 1, sum),
    best_hit_ratio = ifelse(marker_set_size > 0, best_hit_count / marker_set_size, 0)
  )
  
  # 保存全局变量
  marker_match_result <<- list(
    result_table = result_table,
    result_vector = result_vector,
    anno_matrix = match_matrix
  )
  cli::cli_alert_info("结果已保存至全局变量 'marker_match_result'")
  
  # 保存 CSV
  output_dir <- "results/tables/cell_annotation"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  outfile <- file.path(output_dir, paste0("manual_cluster_annotation_", set_name, ".csv"))
  write.csv(result_table, outfile, row.names = FALSE)
  cli::cli_alert_success("结果已保存至：{outfile}")
  
  # 打印输出
  cli::cli_h2("匹配结果")
  cli::cli_alert_info("以下为 Top {top_n} DEG 与 {length(celltypes)} 个细胞类型的匹配结果（marker set: {set_name}）：")
  print(result_table)
  
  if (log) {
    cat("共处理 cluster 数量：", nrow(result_table), "\n")
    cat("匹配结果（marker set: ", set_name, "）：\n")
    print(result_table)
    cat(strrep("-", 70), "\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - End of manual_cluster_annotation\n")
    sink()
    cli::cli_alert_success("日志已保存至：{log_file}")
  }
  
  cli::cli_alert_success("完成 Top {top_n} DEG 与 {length(celltypes)} 个细胞类型的匹配（marker set: {set_name})")
  
  return(list(
    result_table = result_table,
    result_vector = result_vector,
    anno_matrix = match_matrix
  ))
}