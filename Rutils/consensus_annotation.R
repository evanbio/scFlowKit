# Rutils/consensus_annotation.R
#-------------------------------------------------------------------------------
# scFlowKit: Perform Consensus Annotation Using Majority Voting
#-------------------------------------------------------------------------------
#
# 背景说明：
# 本函数基于Seurat对象的多个注释列，通过多数投票法生成统一的共识标签列。
# 适用于整合多种细胞类型注释结果，生成可靠的最终标签。
#
# 本函数主要执行以下步骤：
# - 参数与输入校验
# - 计算共识标签（多数投票法）
# - 将共识标签写入Seurat对象的meta.data
# - 保存额外结果到全局变量 consensus_result
#
# 返回值：
# - 更新后的Seurat对象（包含共识标签列）
# - 额外结果自动保存到全局变量 consensus_result，包含：
#   - anno_table：包含 cell 和 consensus_label 的注释表
#   - anno_vector：以细胞名命名的共识标签向量
#   - label_counts：各共识标签的计数统计
#
#-------------------------------------------------------------------------------

consensus_annotation <- function(
  seu,
  label_cols,
  output_col = "consensus_label",
  min_agreement = 2,
  default = "ambiguous",
  log = TRUE
) {
  #-------------------- 依赖包检查 --------------------
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("请先安装 Seurat 包！", call. = FALSE)
  }
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("请先安装 cli 包！", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("请先安装 tibble 包！", call. = FALSE)
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("请先安装 stringr 包！", call. = FALSE)
  }

  #-------------------- 参数校验 --------------------
  if (!inherits(seu, "Seurat")) {
    stop("参数 'seu' 必须为 Seurat 对象。", call. = FALSE)
  }
  if (!is.character(label_cols) || length(label_cols) == 0) {
    stop("参数 'label_cols' 必须为非空字符向量。", call. = FALSE)
  }
  if (!all(label_cols %in% colnames(seu@meta.data))) {
    missing_cols <- label_cols[!label_cols %in% colnames(seu@meta.data)]
    stop("以下列在 seu@meta.data 中不存在：", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (!is.character(output_col) || length(output_col) != 1) {
    stop("参数 'output_col' 必须为单一字符值。", call. = FALSE)
  }
  if (!is.numeric(min_agreement) || min_agreement <= 0 || min_agreement != round(min_agreement)) {
    stop("参数 'min_agreement' 必须为正整数。", call. = FALSE)
  }
  if (!is.character(default) || length(default) != 1) {
    stop("参数 'default' 必须为单一字符值。", call. = FALSE)
  }
  if (!is.logical(log) || length(log) != 1) {
    stop("参数 'log' 必须为单一逻辑值（TRUE 或 FALSE）。", call. = FALSE)
  }

  #-------------------- 日志设置 --------------------
  log_file <- "logs/cell_annotation/consensus_annotation.log"

  if (log) {
    dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
    sink(log_file, append = TRUE, split = TRUE)

    cat("\n", strrep("-", 70), "\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Start of consensus annotation run\n")
    cat("• Input cells      : ", ncol(seu), "\n")
    cat("• Label columns   : ", paste(label_cols, collapse = ", "), "\n")
    cat("• Output column   : ", output_col, "\n")
    cat("• Min agreement   : ", min_agreement, "\n")
    cat("• Default label   : ", default, "\n")
    cat(strrep("-", 70), "\n")
  }

  #-------------------------------------------------------------------------------
  # 执行共识计算
  #-------------------------------------------------------------------------------

  cli::cli_h1("执行共识注释")

  cli::cli_alert_info("基于 {length(label_cols)} 个注释列计算共识标签...")
  labels_df <- seu@meta.data[, label_cols, drop = FALSE]

  consensus_vector <- apply(labels_df, 1, function(x) {
    x <- x[!is.na(x) & x != "" & x != "ambiguous"]  # 过滤掉 "ambiguous" 标签
    if (length(x) == 0) return(default)
    label_counts <- table(x)
    top_label <- names(label_counts)[which.max(label_counts)]
    if (label_counts[[top_label]] >= min_agreement) {
      return(top_label)
    } else {
      return(default)
    }
  })

  # 将共识标签写入 Seurat 对象
  seu[[output_col]] <- consensus_vector

  # 创建注释表
  anno_table <- tibble::tibble(
    cell = colnames(seu),
    !!output_col := consensus_vector
  )

  # 创建命名向量
  anno_vector <- consensus_vector
  names(anno_vector) <- colnames(seu)

  # 计算标签分布
  label_counts <- table(consensus_vector)

  # 保存额外结果到全局变量
  consensus_result <<- list(
    anno_table = anno_table,
    anno_vector = anno_vector,
    label_counts = label_counts
  )
  cli::cli_alert_info("额外结果已保存至全局变量 'consensus_result'")

  # 打印分布统计
  cli::cli_h2("共识标签分布")
  cli::cli_alert_info("以下为 {output_col} 的计数统计：")
  print(label_counts)

  #-------------------------------------------------------------------------------
  # 保存结果
  #-------------------------------------------------------------------------------

  output_dir <- "results/tables/cell_annotation"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  outfile <- file.path(output_dir, "consensus_annotation.csv")
  write.csv(anno_table, outfile, row.names = FALSE)
  cli::cli_alert_success("共识注释结果已保存至：{outfile}")

  # 如果 log = TRUE，关闭 sink
  if (log) {
    cat(strrep("-", 70), "\n")
    cat("Consensus annotation complete\n")
    cat("共识标签已统计，详见上方 R 控制台输出。\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - End of consensus annotation run\n")
    sink()
    cli::cli_alert_success("日志已保存至：{log_file}")
  }

  #-------------------------------------------------------------------------------
  # 返回结果
  #-------------------------------------------------------------------------------

  cli::cli_alert_success("共识注释完成！")
  return(seu)
}