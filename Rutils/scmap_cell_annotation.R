# Rutils/scmap_cell_annotation.R
#-------------------------------------------------------------------------------
# scFlowKit: Perform Cell Type Annotation Using scmap-cell
#-------------------------------------------------------------------------------
#
# 背景说明：
# scmap-cell 是 scmap 包提供的另一种细胞类型注释方法。它通过将目标细胞的表达谱
# 与参考数据集中每个单细胞进行相似性计算，从而预测细胞类型，适用于更加细粒度的注释。
#
# 本函数主要执行以下步骤：
# - 参数与输入校验
# - 使用 scmap::selectFeatures 选择特征基因（支持 include/exclude 自定义）
# - 构建 cell-level 索引（scmap::indexCell）
# - 进行注释（scmap::scmapCell）
# - 返回包含特征基因、索引、注释表和向量的结果列表
#
# 返回值：
# - 一个列表，包括：
#   - feature_genes：特征基因名称
#   - cell_index：参考数据中构建的索引
#   - scmap_cell_annotation：scmap::scmapCell 的原始结果
#   - anno_table：包含 cell 与 scmap_cell_label 的注释表
#   - anno_vector：以细胞名命名的标签向量
#
#-------------------------------------------------------------------------------

scmap_cell_annotation <- function(
  target_sce,
  ref_sce,
  label_col,
  include_genes = NULL,
  exclude_genes = NULL,
  n_features = 500,
  w = 10,
  log = TRUE
) {
  #-------------------- 依赖包检查 --------------------
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
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("请先安装 stringr 包！", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("请先安装 tibble 包！", call. = FALSE)
  }
  #-------------------- 参数校验 --------------------
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

  # 检查 w
  if (!is.numeric(w) || w <= 0 || w != round(w)) {
    stop("参数 'w' 必须为正整数。", call. = FALSE)
  }

  # 检查 log
  if (!is.logical(log) || length(log) != 1) {
    stop("参数 'log' 必须为单一逻辑值（TRUE 或 FALSE）。", call. = FALSE)
  }

  #-------------------- 日志设置 --------------------
  log_file <- "logs/cell_annotation/scmap_cell_annotation.log"

  if (log) {
    dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
    sink(log_file, append = TRUE, split = TRUE)

    cat("\n", strrep("-", 70), "\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Start of scmap-cell run\n")
    cat("• Reference cells: ", ncol(ref_sce), "\n")
    cat("• Target cells   : ", ncol(target_sce), "\n")
    cat("• Label column   : ", label_col, "\n")
    cat("• n_features     : ", n_features, "\n")
    cat("• w              : ", w, "\n")
    cat(strrep("-", 70), "\n")
  }

  #-------------------------------------------------------------------------------
  # 数据预处理：设置细胞类型列并创建 feature_symbol
  #-------------------------------------------------------------------------------

  cli::cli_h1("执行 scmap-cell 细胞类型注释")

  ref_sce$cell_type1 <- as.factor(ref_sce[[label_col]])
  SummarizedExperiment::rowData(ref_sce)$feature_symbol <- rownames(ref_sce)
  SummarizedExperiment::rowData(target_sce)$feature_symbol <- rownames(target_sce)

  #-------------------------------------------------------------------------------
  # 特征选择：选择 Top N 个特征基因并根据 include/exclude_genes 调整
  #-------------------------------------------------------------------------------

  cli::cli_alert_info("选择 Top {n_features} 个特征基因...")
  
  # 选择特征基因
  ref_sce <- scmap::selectFeatures(ref_sce, n_features = n_features, suppress_plot = TRUE)
  
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
  # 构建索引：在 metadata 中生成 scmap_cell_index
  # - M：将表达矩阵分成 M 个子集（用于加速计算）
  #      默认行为：若特征基因数量 ≤ 1000，则 M = n_features / 10，否则 M = 100
  # - k：对每个子集执行 k-means 聚类（用于索引压缩）
  #      默认行为：k = sqrt(参考细胞数量)
  # - 可通过 scmap::indexCell(..., M = , k = ) 手动指定
  #-------------------------------------------------------------------------------

  cli::cli_alert_info("构建 scmap-cell 索引...")
  
  ref_sce <- scmap::indexCell(ref_sce)
  cell_index <- metadata(ref_sce)$scmap_cell_index

  #-------------------------------------------------------------------------------
  # 执行注释：使用 scmapCell 预测细胞类型
  #-------------------------------------------------------------------------------

  cli::cli_alert_info("执行 scmap-cell 注释（w = {w}）...")
  scmap_cell_res <- scmap::scmapCell(
    projection = target_sce,
    index_list = list(ref = cell_index),
    w = w
  )

  #-------------------------------------------------------------------------------
  # 结果处理：提取注释表和向量
  #-------------------------------------------------------------------------------

  # 最近邻矩阵（列为 target cells）
  nn_cells <- scmap_cell_res$ref$cells  # matrix: w x N
  ref_labels <- as.character(ref_sce[[label_col]])

  # 定义众数提取函数
  get_majority_label <- function(indices) {
    labels <- ref_labels[indices]
    tbl <- table(labels)
    top <- which(tbl == max(tbl))  # 找到出现次数最多的标签（可能多个）

    if (length(top) == 1) {
      return(names(top))           # 唯一众数
    } else {
      return("ambiguous")          # 多个并列最多
    }
  } 

  # 应用于每一列（即每个 target cell）
  predicted_labels <- apply(nn_cells, 2, get_majority_label)

  # 创建两列数据框
  anno_table <- tibble::tibble(
    cell = colnames(target_sce),
    scmap_cell_label = predicted_labels
  )

  # 创建命名向量
  anno_vector <- predicted_labels
  names(anno_vector) <- colnames(target_sce)

  # 打印注释分布
  cli::cli_h2("注释结果分布")
  cli::cli_alert_info("以下为 scmap_cell_label 的计数统计：")
  print(table(predicted_labels))

  # 如果 log = TRUE，关闭 sink
  if (log) {
    cat(strrep("-", 70), "\n")
    cat("scmap-cell annotation complete\n")
    cat("注释标签已统计，详见上方 R 控制台输出。\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - End of scmap-cell run\n")
    sink()
    cli::cli_alert_success("日志已保存至：{log_file}")
  }
  #-------------------------------------------------------------------------------
  # 保存结果：仅保存两列数据框为 CSV
  #-------------------------------------------------------------------------------

  output_dir <- "results/tables/cell_annotation"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  outfile <- file.path(output_dir, "scmap_cell_annotation.csv")
  write.csv(anno_table, outfile, row.names = FALSE)
  cli::cli_alert_success("注释结果已保存至：{outfile}")

  #-------------------------------------------------------------------------------
  # 返回结果列表
  #-------------------------------------------------------------------------------

  cli::cli_alert_success("细胞类型注释完成！")
  return(list(
    feature_genes = feature_genes,
    cell_index = cell_index,
    scmap_cell_annotation = scmap_cell_res,
    anno_table = anno_table,
    anno_vector = anno_vector
  ))
}
