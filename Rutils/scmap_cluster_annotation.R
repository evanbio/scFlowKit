# Rutils/scmap_cluster_annotation.R
#-------------------------------------------------------------------------------
# scFlowKit: Perform Cell Type Annotation Using scmap-cluster
#-------------------------------------------------------------------------------
#
# 背景说明：
# 在单细胞 RNA 测序分析中，scmap-cluster 是一种高效的细胞类型注释方法，通过将查询数据集
#（target SCE）投影到参考数据集（reference SCE）的聚类索引上，预测细胞的类型。
#
# 本函数的主要步骤：
# - 检查输入的 SingleCellExperiment 对象及其 logcounts 层。
# - 使用 scmap::selectFeatures 选择 Top N 个特征基因。
# - 根据 include_genes 强制包含基因，exclude_genes 强制排除基因。
# - 使用 scmap::indexCluster 构建索引并执行注释。
# - 返回包含特征基因、索引、原始结果和注释表的列表。
#
# 返回值：
# - 一个列表，包含以下元素：
#   - feature_genes：选择的特征基因名称向量。
#   - cluster_index：在参考数据 metadata 中生成的 scmap_cluster_index 数据框。
#   - scmap_cluster_annotation：scmap::scmapCluster 的原始结果。
#   - anno_table：两列数据框（细胞名和预测标签）。
#   - anno_vector：以细胞名为名的预测标签向量。
#
#-------------------------------------------------------------------------------

scmap_cluster_annotation <- function(
  target_sce,           # 查询 SingleCellExperiment 对象
  ref_sce,              # 参考 SingleCellExperiment 对象
  label_col,            # 参考数据中用于注释的列名
  include_genes = NULL, # 可选：强制包含的基因集合（大小写不敏感匹配）
  exclude_genes = NULL, # 可选：强制排除的基因集合（大小写不敏感匹配）
  n_features = 500,     # 特征基因数量，默认 500
  threshold = 0.1,      # 相似性阈值，默认 0.1
  log = TRUE            # 是否将日志写入文件
) {
  #-------------------------------------------------------------------------------
  # 参数校验：确保输入为合法参数
  #-------------------------------------------------------------------------------

  # 检查必要包
  if (!requireNamespace("scmap", quietly = TRUE)) {
    stop("请先安装 scmap 包！", call. = FALSE)
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("请先安装 Seurat 包！", call. = FALSE)
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
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("请先安装 dplyr 包！", call. = FALSE)
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

  # 检查 label_col
  if (!is.character(label_col) || length(label_col) != 1 ||
      !label_col %in% colnames(SummarizedExperiment::colData(ref_sce))) {
    stop("参数 'label_col' 必须为单一字符值且存在于 ref_sce 的 colData 中！",
         call. = FALSE)
  }

  # 检查 include_genes 和 exclude_genes
  if (!is.null(include_genes) && (!is.character(include_genes) ||
                                  length(include_genes) == 0)) {
    stop("参数 'include_genes' 必须为非空字符向量或 NULL。", call. = FALSE)
  }
  if (!is.null(exclude_genes) && (!is.character(exclude_genes) ||
                                  length(exclude_genes) == 0)) {
    stop("参数 'exclude_genes' 必须为非空字符向量或 NULL。", call. = FALSE)
  }

  # 检查 n_features
  if (!is.numeric(n_features) || n_features <= 0 || n_features != round(n_features)) {
    stop("参数 'n_features' 必须为正整数。", call. = FALSE)
  }

  # 检查 threshold
  if (!is.numeric(threshold) || threshold < 0 || threshold > 1) {
    stop("参数 'threshold' 必须为 0 到 1 之间的数值。", call. = FALSE)
  }

  # 检查 log
  if (!is.logical(log) || length(log) != 1) {
    stop("参数 'log' 必须为单一逻辑值（TRUE 或 FALSE）。", call. = FALSE)
  }

  # 设置日志文件路径
  log_file <- "logs/cell_annotation/scmap_cluster_annotation.log"

  # 如果 log = TRUE，启动 sink 捕获控制台输出
  if (log) {
    dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
    sink(log_file, append = TRUE, split = TRUE)
    cat("\n", strrep("-", 70), "\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Start of scmap-cluster run\n")
    cat(strrep("-", 70), "\n")
  }

  #-------------------------------------------------------------------------------
  # 数据预处理：设置细胞类型列并创建 feature_symbol
  #-------------------------------------------------------------------------------

  cli::cli_h1("执行 scmap-cluster 细胞类型注释")

  # 将用户指定的 label_col 转换为 cell_type1 并转为因子
  ref_sce$cell_type1 <- as.factor(ref_sce[[label_col]])

  # 为参考和查询数据创建 feature_symbol
  SummarizedExperiment::rowData(ref_sce)$feature_symbol <- rownames(ref_sce)
  SummarizedExperiment::rowData(target_sce)$feature_symbol <- rownames(target_sce)

  #-------------------------------------------------------------------------------
  # 特征选择：选择 Top N 个特征基因并根据 include/exclude_genes 调整
  #-------------------------------------------------------------------------------

  cli::cli_alert_info("选择 Top {n_features} 个特征基因...")

  # 选择特征基因
  ref_sce <- scmap::selectFeatures(ref_sce, n_features = n_features,
                                   suppress_plot = TRUE)

  # 获取所有基因
  all_genes <- rownames(ref_sce)

  # 强制包含 include_genes（设置为 TRUE）
  if (!is.null(include_genes)) {
    matched_include <- Seurat::CaseMatch(include_genes, all_genes)
    if (length(matched_include) == 0) {
      cli::cli_alert_warning("未匹配到任何指定的 include_genes，跳过强制包含。")
    } else {
      SummarizedExperiment::rowData(ref_sce)$scmap_features[
        all_genes %in% matched_include] <- TRUE
    }
  }

  # 强制排除 exclude_genes（设置为 FALSE）
  if (!is.null(exclude_genes)) {
    matched_exclude <- Seurat::CaseMatch(exclude_genes, all_genes)
    if (length(matched_exclude) == 0) {
      cli::cli_alert_warning("未匹配到任何指定的 exclude_genes，跳过强制排除。")
    } else {
      SummarizedExperiment::rowData(ref_sce)$scmap_features[
        all_genes %in% matched_exclude] <- FALSE
    }
  }

  feature_genes <- rownames(ref_sce)[
    SummarizedExperiment::rowData(ref_sce)$scmap_features]
  cli::cli_alert_success("共选择 {length(feature_genes)} 个特征基因")

  #-------------------------------------------------------------------------------
  # 构建索引：在 metadata 中生成 scmap_cluster_index
  #-------------------------------------------------------------------------------

  cli::cli_alert_info("构建 scmap-cluster 索引...")
  ref_sce <- scmap::indexCluster(ref_sce)
  cluster_index <- metadata(ref_sce)$scmap_cluster_index

  cli::cli_alert_success(
    paste0("索引构建完成：", nrow(cluster_index), " 行 × ", ncol(cluster_index), " 列")
  )
  #-------------------------------------------------------------------------------
  # 执行注释：使用 scmapCluster 预测细胞类型
  #-------------------------------------------------------------------------------

  cli::cli_alert_info("执行细胞类型注释（threshold = {threshold}）...")
  scmap_cluster_res <- scmap::scmapCluster(
    projection = target_sce,
    index_list = list(reference = cluster_index),
    threshold = threshold
  )

  #-------------------------------------------------------------------------------
  # 结果处理：提取注释表和向量
  #-------------------------------------------------------------------------------

  # 创建两列数据框
  anno_table <- tibble::tibble(
    cell = colnames(target_sce),
    scmap_cluster_label = as.character(scmap_cluster_res$combined_labs)
  )

  # 创建命名向量
  anno_vector <- scmap_cluster_res$combined_labs
  names(anno_vector) <- colnames(target_sce)

  # 打印注释分布
  cli::cli_h2("注释结果分布")
  cli::cli_alert_info("以下为 scmap_label 的计数统计：")
  print(table(scmap_cluster_res$combined_labs))

  # 如果 log = TRUE，关闭 sink
  if (log) {
    cat(strrep("-", 70), "\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - End of scmap-cluster run\n")
    sink()
    cli::cli_alert_success("日志已保存至：{log_file}")
  }

  #-------------------------------------------------------------------------------
  # 保存结果：仅保存两列数据框为 CSV
  #-------------------------------------------------------------------------------

  output_dir <- "results/tables/cell_annotation"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  outfile <- file.path(output_dir, "scmap_cluster_annotation.csv")
  write.csv(anno_table, outfile, row.names = FALSE)
  cli::cli_alert_success("注释结果已保存至：{outfile}")

  #-------------------------------------------------------------------------------
  # 返回结果列表
  #-------------------------------------------------------------------------------

  cli::cli_alert_success("细胞类型注释完成！")
  return(list(
    feature_genes = feature_genes,
    cluster_index = cluster_index,
    scmap_cluster_annotation = scmap_cluster_res,
    anno_table = anno_table,
    anno_vector = anno_vector
  ))
}

#-------------------------------------------------------------------------------