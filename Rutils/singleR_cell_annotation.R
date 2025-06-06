# Rutils/singleR_cell_annotation.R
#-------------------------------------------------------------------------------
# scFlowKit: Perform Cell-Level Cell Type Annotation Using SingleR
#-------------------------------------------------------------------------------
#
# 背景说明：
# SingleR 是一种基于参考表达谱的自动细胞类型注释工具。本函数使用 SingleR 的 cell mode，
# 通过将目标数据中的每个细胞与参考数据集的表达谱进行比较，为每个细胞分配一个独立的
# 细胞类型标签。
#
# 本函数的主要步骤：
# - 检查输入的 SingleCellExperiment 对象（target_sce 和 ref_sce）及其 logcounts 层。
# - 检查 labels 的有效性。
# - 执行 SingleR cell-level 注释。
# - 处理结果并返回包含原始结果和细胞标签的列表。
#
# 返回值：
# 本函数返回一个列表，包含以下三个元素，提供了从原始输出到最终注释的完整结果：
# - singleR_cell_annotation: SingleR::SingleR 的原始结果，是一个 DFrame 对象，包含每个细胞的
#   得分、初步标签和剪枝后的最终标签等详细信息。
# - anno_table: 数据框，列为 `cell`（细胞 ID）和 `singleR_cell_label`（预测标签），展示每个
#   细胞的注释结果。
# - anno_vector: 命名向量，以细胞 ID 为名，值为预测标签，用于直接整合到 Seurat 对象。
#
#-------------------------------------------------------------------------------

singleR_cell_annotation <- function(
  target_sce,           # 查询 SingleCellExperiment 对象
  ref_sce,              # 参考 SingleCellExperiment 对象
  labels,               # 参考数据的标签向量（字符向量或因子）
  genes = "de",         # 特征基因选择方法，默认 "de"
  sd.thresh = 1,        # 标准差过滤阈值，默认 1
  de.method = "classic",# 差异分析方法，默认 "classic"
  quantile = 0.8,       # 相似性过滤的分位数，默认 0.8
  fine.tune = TRUE,     # 是否进行微调，默认 TRUE
  prune = TRUE,         # 是否进行标签剪枝，默认 TRUE
  assay.type.test = "logcounts", # 测试集 assay 层，默认 "logcounts"
  assay.type.ref = "logcounts",  # 参考集 assay 层，默认 "logcounts"
  log = TRUE            # 是否将日志写入文件
) {
  #-------------------------------------------------------------------------------
  # 参数校验：确保输入为合法参数
  #-------------------------------------------------------------------------------

  # 检查必要包
  if (!requireNamespace("SingleR", quietly = TRUE)) {
    stop("请先安装 SingleR 包！", call. = FALSE)
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("请先安装 SingleCellExperiment 包！", call. = FALSE)
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

  # 检查输入对象类型
  if (!inherits(target_sce, "SingleCellExperiment")) {
    stop("参数 'target_sce' 必须为 SingleCellExperiment 对象。", call. = FALSE)
  }
  if (!inherits(ref_sce, "SingleCellExperiment")) {
    stop("参数 'ref_sce' 必须为 SingleCellExperiment 对象。", call. = FALSE)
  }

  # 检查 logcounts 层
  if (!"logcounts" %in% SummarizedExperiment::assayNames(target_sce)) {
    stop("参数 'target_sce' 缺少 logcounts 层！", call. = FALSE)
  }
  if (!"logcounts" %in% SummarizedExperiment::assayNames(ref_sce)) {
    stop("参数 'ref_sce' 缺少 logcounts 层！", call. = FALSE)
  }

  # 检查 labels
  if (!is.character(labels) && !is.factor(labels)) {
    stop("参数 'labels' 必须为字符向量或因子。", call. = FALSE)
  }
  if (length(labels) != ncol(ref_sce)) {
    stop("参数 'labels' 的长度 (", length(labels), ") 与 ref_sce 的细胞数 (", 
         ncol(ref_sce), ") 不匹配！", call. = FALSE)
  }

  # 检查 log
  if (!is.logical(log) || length(log) != 1) {
    stop("参数 'log' 必须为单一逻辑值（TRUE 或 FALSE）。", call. = FALSE)
  }

  # 检查其他参数
  if (!genes %in% c("de", "sd", "all")) {
    stop("参数 'genes' 必须为 'de', 'sd' 或 'all'。", call. = FALSE)
  }
  if (!is.numeric(sd.thresh) || sd.thresh < 0) {
    stop("参数 'sd.thresh' 必须为非负数值。", call. = FALSE)
  }
  if (!de.method %in% c("classic", "wilcox", "t")) {
    stop("参数 'de.method' 必须为 'classic', 'wilcox' 或 't'。", call. = FALSE)
  }
  if (!is.numeric(quantile) || quantile < 0 || quantile > 1) {
    stop("参数 'quantile' 必须为 0 到 1 之间的数值。", call. = FALSE)
  }
  if (!is.logical(fine.tune) || length(fine.tune) != 1) {
    stop("参数 'fine.tune' 必须为单一逻辑值。", call. = FALSE)
  }
  if (!is.logical(prune) || length(prune) != 1) {
    stop("参数 'prune' 必须为单一逻辑值。", call. = FALSE)
  }

  # 设置日志文件路径
  log_file <- "logs/cell_annotation/singleR_cell_annotation.log"
  if (log) {
    dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
    sink(log_file, append = TRUE, split = TRUE)
    cat("\n", strrep("-", 70), "\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Start of SingleR cell run\n")
    cat("• Reference cells: ", ncol(ref_sce), "\n")
    cat("• Target cells   : ", ncol(target_sce), "\n")
    cat("• Genes method   : ", genes, "\n")
    cat("• sd.thresh      : ", sd.thresh, "\n")
    cat("• de.method      : ", de.method, "\n")
    cat("• quantile       : ", quantile, "\n")
    cat(strrep("-", 70), "\n")
  }

  #-------------------------------------------------------------------------------
  # 数据预处理：准备 labels
  #-------------------------------------------------------------------------------

  cli::cli_h1("执行 SingleR cell-level 细胞类型注释")

  #-------------------------------------------------------------------------------
  # 执行注释：使用 SingleR 进行 cell-level 预测
  #-------------------------------------------------------------------------------

  cli::cli_alert_info("执行 SingleR cell-level 注释...")
  singleR_cell_res <- SingleR::SingleR(
    test = target_sce,
    ref = ref_sce,
    labels = labels,
    genes = genes,
    sd.thresh = sd.thresh,
    de.method = de.method,
    quantile = quantile,
    fine.tune = fine.tune,
    prune = prune,
    assay.type.test = assay.type.test,
    assay.type.ref = assay.type.ref
  )

  #-------------------------------------------------------------------------------
  # 结果处理：提取注释表和向量
  #-------------------------------------------------------------------------------

  # 提取细胞标签
  anno_vector <- singleR_cell_res$pruned.labels
  names(anno_vector) <- rownames(singleR_cell_res)
  na_count <- sum(is.na(anno_vector))
  if (na_count > 0) {
    cli::cli_alert_warning("{na_count} 个细胞标签为 NA，已替换为 'ambiguous'")
  }
  anno_vector[is.na(anno_vector)] <- "ambiguous"

  # 创建细胞级注释表
  anno_table <- tibble::tibble(
    cell = names(anno_vector),
    singleR_cell_label = anno_vector
  )

  # 打印注释分布
  cli::cli_h2("注释结果分布")
  cli::cli_alert_info("以下为 singleR_cell_label 的计数统计：")
  print(table(anno_vector))

  # 如果 log = TRUE，关闭 sink
  if (log) {
    cat(strrep("-", 70), "\n")
    cat("SingleR cell annotation complete\n")
    cat("注释标签已统计，详见上方 R 控制台输出。\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - End of SingleR cell run\n")
    sink()
    cli::cli_alert_success("日志已保存至：{log_file}")
  }

  #-------------------------------------------------------------------------------
  # 保存结果：保存两列数据框为 CSV
  #-------------------------------------------------------------------------------

  output_dir <- "results/tables/cell_annotation"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  outfile <- file.path(output_dir, "singleR_cell_annotation.csv")
  write.csv(anno_table, outfile, row.names = FALSE)
  cli::cli_alert_success("注释结果已保存至：{outfile}")

  #-------------------------------------------------------------------------------
  # 返回结果列表
  #-------------------------------------------------------------------------------

  cli::cli_alert_success("细胞类型注释完成！")
  return(list(
    singleR_cell_annotation = singleR_cell_res,
    anno_table = anno_table,
    anno_vector = anno_vector
  ))
}

#-------------------------------------------------------------------------------