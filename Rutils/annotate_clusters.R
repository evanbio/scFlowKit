# Rutils/annotate_clusters.R
#-------------------------------------------------------------------------------
# scFlowKit: Assign Cluster Labels Based on Cell Labels
#-------------------------------------------------------------------------------
#
# 背景说明：
# 本函数基于Seurat对象的细胞标签和聚类标签，为每个聚类分配主导标签。
# 如果某细胞标签在聚类中的比例达到指定阈值，则采用该标签，否则标记为默认值。
#
# 本函数主要执行以下步骤：
# - 校验输入参数
# - 计算每个聚类的细胞标签比例
# - 为聚类分配标签（满足比例阈值或默认值）
# - 将结果写入Seurat对象的meta.data
# - 保存额外结果到全局变量 cluster_anno_result
#
# 参数说明：
# - seu           : Seurat 对象
# - cell_label_col: 细胞标签列名（默认 "consensus_label"）
# - cluster_col   : 聚类标签列名（默认 "seurat_clusters"）
# - min_ratio     : 最低比例阈值（默认 0.7）
# - default       : 未达阈值时的默认标签（默认 "ambiguous"）
# - output_col    : 输出列名（默认 "consensus_anno_cluster"）
# - log           : 是否记录日志
#
# 返回值：
# - 更新后的Seurat对象（包含新列 output_col）
# - 额外结果自动保存到全局变量 cluster_anno_result，包含：
#   - result_table：聚类和分配标签的数据框
#   - result_vector：聚类标签向量
#   - anno_vector：细胞级的映射标签向量
#
#-------------------------------------------------------------------------------

library(tidyverse)

annotate_clusters <- function(
    seu,
    cell_label_col = "consensus_label",
    cluster_col = "seurat_clusters",
    min_ratio = 0.7,
    default = "ambiguous",
    output_col = "consensus_anno_cluster",
    log = TRUE
) {
  #-------------------- 依赖包检查 --------------------
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("请先安装 Seurat 包！", call. = FALSE)
  }
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("请先安装 cli 包！", call. = FALSE)
  }
  if (!requireNamespace("tidyverse", quietly = TRUE)) {
    stop("请先安装 tidyverse 包！", call. = FALSE)
  }
  
  #-------------------- 参数校验 --------------------
  if (!inherits(seu, "Seurat")) {
    stop("参数 'seu' 必须为 Seurat 对象。", call. = FALSE)
  }
  if (!is.character(cell_label_col) || length(cell_label_col) != 1) {
    stop("参数 'cell_label_col' 必须为单一字符值。", call. = FALSE)
  }
  if (!is.character(cluster_col) || length(cluster_col) != 1) {
    stop("参数 'cluster_col' 必须为单一字符值。", call. = FALSE)
  }
  if (!cell_label_col %in% colnames(seu@meta.data)) {
    stop("列 '", cell_label_col, "' 在 seu@meta.data 中不存在！", call. = FALSE)
  }
  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop("列 '", cluster_col, "' 在 seu@meta.data 中不存在！", call. = FALSE)
  }
  if (!is.numeric(min_ratio) || min_ratio < 0 || min_ratio > 1) {
    stop("参数 'min_ratio' 必须为 0 到 1 之间的数值。", call. = FALSE)
  }
  if (!is.character(default) || length(default) != 1) {
    stop("参数 'default' 必须为单一字符值。", call. = FALSE)
  }
  if (!is.character(output_col) || length(output_col) != 1) {
    stop("参数 'output_col' 必须为单一字符值。", call. = FALSE)
  }
  if (!is.logical(log) || length(log) != 1) {
    stop("参数 'log' 必须为单一逻辑值（TRUE 或 FALSE）。", call. = FALSE)
  }
  
  #-------------------- 日志设置 --------------------
  log_file <- file.path("logs/cell_annotation", paste0(output_col, ".log"))
  
  if (log) {
    dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
    sink(log_file, append = TRUE, split = TRUE)
    
    cat("\n", strrep("-", 70), "\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Start of cluster annotation run\n")
    cat("• Input cells      : ", ncol(seu), "\n")
    cat("• Cell label col  : ", cell_label_col, "\n")
    cat("• Cluster col     : ", cluster_col, "\n")
    cat("• Min ratio       : ", min_ratio, "\n")
    cat("• Default label   : ", default, "\n")
    cat("• Output column   : ", output_col, "\n")
    cat(strrep("-", 70), "\n")
  }
  
  #-------------------------------------------------------------------------------
  # 计算聚类标签
  #-------------------------------------------------------------------------------
  
  cli::cli_h1("执行聚类标签分配")
  
  cli::cli_alert_info("基于 {cell_label_col} 为聚类分配标签...")
  
  # 提取数据
  df <- tibble(
    cluster = as.character(seu@meta.data[[cluster_col]]),
    cell_label = as.character(seu@meta.data[[cell_label_col]]),
    cell = colnames(seu)
  )
  
  # 计算每个聚类的标签比例并分配主导标签
  result_table <- df %>%
    group_by(cluster, cell_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(cluster) %>%
    mutate(prop = n / sum(n)) %>%
    arrange(cluster, desc(n)) %>%
    slice_head(n = 1) %>%
    mutate(
      !!output_col := if_else(prop >= min_ratio, cell_label, default)
    ) %>%
    select(cluster, !!output_col)
  
  # 映射回细胞
  cluster_to_label <- deframe(result_table %>% select(cluster, !!output_col))
  anno_vector <- cluster_to_label[df$cluster]
  names(anno_vector) <- df$cell
  
  # 创建 result_vector（聚类标签向量）
  result_vector <- cluster_to_label
  names(result_vector) <- result_table$cluster
  
  # 写入 Seurat 对象
  seu[[output_col]] <- anno_vector
  
  # 保存额外结果到全局变量
  cluster_anno_result <<- list(
    result_table = result_table,
    result_vector = result_vector,
    anno_vector = anno_vector
  )
  cli::cli_alert_info("额外结果（result_table, result_vector, anno_vector）已保存至全局变量 'cluster_anno_result'")
  
  # 打印 result_table
  cli::cli_h2("聚类标签结果")
  cli::cli_alert_info("以下为 {output_col} 的分配结果：")
  print(result_table)
  
  #-------------------------------------------------------------------------------
  # 保存结果
  #-------------------------------------------------------------------------------
  
  output_dir <- "results/tables/cell_annotation"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  outfile <- file.path(output_dir, paste0(output_col, ".csv"))
  write.csv(result_table, outfile, row.names = FALSE)
  cli::cli_alert_success("聚类标签结果已保存至：{outfile}")
  
  # 如果 log = TRUE，关闭 sink
  if (log) {
    cat(strrep("-", 70), "\n")
    cat("Cluster annotation complete\n")
    cat("聚类标签已统计，详见上方 R 控制台输出。\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - End of cluster annotation run\n")
    sink()
    cli::cli_alert_success("日志已保存至：{log_file}")
  }
  
  #-------------------------------------------------------------------------------
  # 返回结果
  #-------------------------------------------------------------------------------
  
  cli::cli_alert_success("聚类标签分配完成！")
  return(seu)
}