# Rutils/filter_cells.R
#-------------------------------------------------------------------------------

# scFlowKit: Filter Low-Quality Cells for Single-Cell RNA-seq Data
#-------------------------------------------------------------------------------

# filter_cells: 过滤低质量细胞
# 参数:
#   sce: Seurat 对象，包含单细胞 RNA-seq 数据
#   min_umi: 最小 UMI 计数，默认 500
#   max_umi: 最大 UMI 计数，默认 Inf（不限制，用户可指定）
#   min_genes: 最小基因数，默认 200
#   max_genes: 最大基因数，默认 Inf（不限制，用户可指定）
#   max_mito: 最大线粒体基因比例（百分比），默认 10
#   min_ratio: 最小 log10_ratio_features_to_umi，默认 0.8
#   filter_hb: 是否过滤红细胞基因比例高的细胞，默认 FALSE
#   max_hb: 最大红细胞基因比例（百分比），默认 5（仅当 filter_hb = TRUE 时生效）
#   filter_ribo: 是否过滤核糖体基因比例高的细胞，默认 FALSE
#   max_ribo: 最大核糖体基因比例（百分比），默认 50（仅当 filter_ribo = TRUE 时生效）
filter_cells <- function(sce, 
                         min_umi = 500, 
                         max_umi = Inf, 
                         min_genes = 200, 
                         max_genes = Inf, 
                         max_mito = 10,
                         min_ratio = 0.8,
                         filter_hb = FALSE, 
                         max_hb = 5,
                         filter_ribo = FALSE, 
                         max_ribo = 50) {
  # 验证输入参数是否为 Seurat 对象
  if (!inherits(sce, "Seurat")) {
    stop("参数 'sce' 必须为 Seurat 对象！", call. = FALSE)
  }

  # 验证数值参数是否为数值类型且合理
  if (!is.numeric(min_umi) || min_umi < 0) {
    stop("参数 'min_umi' 必须为非负数值！", call. = FALSE)
  }
  if (!is.numeric(max_umi) || max_umi <= min_umi) {
    stop("参数 'max_umi' 必须为数值且大于 min_umi！", call. = FALSE)
  }
  if (!is.numeric(min_genes) || min_genes < 0) {
    stop("参数 'min_genes' 必须为非负数值！", call. = FALSE)
  }
  if (!is.numeric(max_genes) || max_genes <= min_genes) {
    stop("参数 'max_genes' 必须为数值且大于 min_genes！", call. = FALSE)
  }
  if (!is.numeric(max_mito) || max_mito < 0) {
    stop("参数 'max_mito' 必须为非负数值！", call. = FALSE)
  }
  if (!is.numeric(min_ratio) || min_ratio < 0) {
    stop("参数 'min_ratio' 必须为非负数值！", call. = FALSE)
  }

  # 验证 filter_hb 和 filter_ribo 是否为逻辑值
  if (!is.logical(filter_hb)) {
    stop("参数 'filter_hb' 必须为逻辑值！", call. = FALSE)
  }
  if (!is.logical(filter_ribo)) {
    stop("参数 'filter_ribo' 必须为逻辑值！", call. = FALSE)
  }

  # 验证 max_hb 和 max_ribo 是否为非负数值
  if (!is.numeric(max_hb) || max_hb < 0) {
    stop("参数 'max_hb' 必须为非负数值！", call. = FALSE)
  }
  if (!is.numeric(max_ribo) || max_ribo < 0) {
    stop("参数 'max_ribo' 必须为非负数值！", call. = FALSE)
  }

  # 验证元数据是否包含必要的质控指标
  required_metrics <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "log10_ratio_features_to_umi")
  missing_metrics <- setdiff(required_metrics, colnames(sce@meta.data))
  if (length(missing_metrics) > 0) {
    stop("元数据缺少必要的质控指标：", paste(missing_metrics, collapse = ", "), call. = FALSE)
  }

  # 如果 filter_hb = TRUE，验证元数据是否包含 percent_hb
  if (filter_hb && !"percent_hb" %in% colnames(sce@meta.data)) {
    stop("元数据缺少 'percent_hb' 指标，请在 calculate_qc_metrics 中设置 calculate_hb = TRUE！", call. = FALSE)
  }

  # 如果 filter_ribo = TRUE，验证元数据是否包含 percent_ribo
  if (filter_ribo && !"percent_ribo" %in% colnames(sce@meta.data)) {
    stop("元数据缺少 'percent_ribo' 指标，请在 calculate_qc_metrics 中设置 calculate_ribo = TRUE！", call. = FALSE)
  }

  # 提示用户正在过滤数据
  message("正在过滤低质量细胞...")

  # 构建过滤条件，使用 sprintf 格式化字符串
  subset_conditions <- list(
    nCount_RNA = sprintf("nCount_RNA > %f & nCount_RNA < %f", min_umi, max_umi),
    nFeature_RNA = sprintf("nFeature_RNA > %f & nFeature_RNA < %f", min_genes, max_genes),
    percent_mito = sprintf("percent_mito < %f", max_mito),
    log10_ratio_features_to_umi = sprintf("log10_ratio_features_to_umi > %f", min_ratio)
  )

  # 如果 filter_hb = TRUE，添加红细胞基因比例过滤条件
  if (filter_hb) {
    subset_conditions$percent_hb <- sprintf("percent_hb < %f", max_hb)
  }

  # 如果 filter_ribo = TRUE，添加核糖体基因比例过滤条件
  if (filter_ribo) {
    subset_conditions$percent_ribo <- sprintf("percent_ribo < %f", max_ribo)
  }

  # 合并过滤条件
  subset_expr <- paste(unlist(subset_conditions), collapse = " & ")

  # 过滤低质量细胞
  sce <- subset(sce, subset = !!rlang::parse_expr(subset_expr))

  # 提示用户过滤完成
  message("低质量细胞过滤完成！")

  return(sce)
}

#-------------------------------------------------------------------------------