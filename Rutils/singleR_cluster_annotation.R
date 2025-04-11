# Rutils/singleR_cluster_annotation.R
#-------------------------------------------------------------------------------
# scFlowKit: Perform Cluster-Level Cell Type Annotation Using SingleR
#-------------------------------------------------------------------------------
#
# 背景说明：
# SingleR 是一种基于参考表达谱的自动细胞类型注释工具。本函数使用 SingleR 的 cluster mode，
# 通过将目标数据中的聚类与参考数据集的表达谱进行比较，为每个聚类分配一个细胞类型标签。
#
# 本函数的主要步骤：
# - 检查输入的 SingleCellExperiment 对象（target_sce 和 ref_sce）及其 logcounts 层。
# - 检查 clusters 和 labels 的有效性。
# - 执行 SingleR cluster-level 注释。
# - 处理结果并返回包含原始结果、聚类注释表和细胞标签的列表。
#
# 返回值：
# - 一个列表，包含以下元素：
#   - singleR_cluster_annotation：SingleR::SingleR 的原始结果（DFrame 对象）。
#   - cluster_anno_table：两列数据框（cluster 和 singleR_cluster_label）。
#   - cluster_anno_vector：以聚类名为名的预测标签向量。
#   - anno_table：两列数据框（cell 和 singleR_cluster_label）。
#   - anno_vector：以细胞名为名的预测标签向量。
#
#-------------------------------------------------------------------------------

singleR_cluster_annotation <- function(
  target_sce,           # 查询 SingleCellExperiment 对象
  ref_sce,              # 参考 SingleCellExperiment 对象
  labels,               # 参考数据的标签向量（字符向量或因子）
  clusters,             # 目标数据的聚类向量（字符向量或因子）
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

  # 检查 clusters
  if (!is.character(clusters) && !is.factor(clusters)) {
    stop("参数 'labels' 必须为字符向量或因子。", call. = FALSE)
  }
  if (length(clusters) != ncol(target_sce)) {
    stop("参数 'clusters' 的长度 (", length(clusters), ") 与 target_sce 的细胞数 (", 
         ncol(target_sce), ") 不匹配！", call. = FALSE)
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
  log_file <- "logs/cell_annotation/singleR_cluster_annotation.log"
  if (log) {
    dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
    sink(log_file, append = TRUE, split = TRUE)
    cat("\n", strrep("-", 70), "\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Start of SingleR cluster run\n")
    cat("• Reference cells: ", ncol(ref_sce), "\n")
    cat("• Target cells   : ", ncol(target_sce), "\n")
    cat("• Genes method   : ", genes, "\n")
    cat("• sd.thresh      : ", sd.thresh, "\n")
    cat("• de.method      : ", de.method, "\n")
    cat("• quantile       : ", quantile, "\n")
    cat(strrep("-", 70), "\n")
  }

  #-------------------------------------------------------------------------------
  # 数据预处理：准备 labels 和 clusters
  #-------------------------------------------------------------------------------

  cli::cli_h1("执行 SingleR cluster-level 细胞类型注释")

  #-------------------------------------------------------------------------------
  # 执行注释：使用 SingleR 进行 cluster-level 预测
  #-------------------------------------------------------------------------------

  cli::cli_alert_info("执行 SingleR cluster-level 注释...")
  singleR_cluster_res <- SingleR::SingleR(
    test = target_sce,
    ref = ref_sce,
    labels = labels,
    clusters = clusters,
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

  # 提取聚类标签
  cluster_anno_vector <- singleR_cluster_res$pruned.labels
  names(cluster_anno_vector) <- rownames(singleR_cluster_res)
  na_count <- sum(is.na(cluster_anno_vector))
  if (na_count > 0) {
    cli::cli_alert_warning("{na_count} 个聚类标签为 NA，已替换为 'ambiguous'")
  }
  cluster_anno_vector[is.na(cluster_anno_vector)] <- "ambiguous"

  # 创建聚类注释表
  cluster_anno_table <- tibble::tibble(
    cluster = names(cluster_anno_vector),
    singleR_cluster_label = cluster_anno_vector
  )

  # 创建细胞级注释表
  anno_table <- tibble::tibble(
    cell = colnames(target_sce),
    singleR_cluster_label = cluster_anno_vector[as.character(clusters)]
  )

  # 创建命名向量
  anno_vector <- anno_table$singleR_cluster_label
  names(anno_vector) <- anno_table$cell

  # 打印注释分布
  cli::cli_h2("注释结果分布")
  cli::cli_alert_info("以下为 cluster_anno_table 的结果：")
  print(cluster_anno_table)
  cli::cli_alert_info("以下为 singleR_cluster_label 的计数统计：")
  print(table(anno_vector))

  # 如果 log = TRUE，关闭 sink
  if (log) {
    cat(strrep("-", 70), "\n")
    cat("SingleR cluster annotation complete\n")
    cat("注释标签已统计，详见上方 R 控制台输出。\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - End of SingleR cluster run\n")
    sink()
    cli::cli_alert_success("日志已保存至：{log_file}")
  }

  #-------------------------------------------------------------------------------
  # 保存结果：保存两列数据框为 CSV
  #-------------------------------------------------------------------------------

  output_dir <- "results/tables/cell_annotation"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  outfile <- file.path(output_dir, "singleR_cluster_annotation.csv")
  write.csv(anno_table, outfile, row.names = FALSE)
  cli::cli_alert_success("注释结果已保存至：{outfile}")

  #-------------------------------------------------------------------------------
  # 返回结果列表
  #-------------------------------------------------------------------------------

  cli::cli_alert_success("细胞类型注释完成！")
  return(list(
    singleR_cluster_annotation = singleR_cluster_res,
    cluster_anno_table = cluster_anno_table,
    cluster_anno_vector = cluster_anno_vector,
    anno_table = anno_table,
    anno_vector = anno_vector
  ))
}

#-------------------------------------------------------------------------------