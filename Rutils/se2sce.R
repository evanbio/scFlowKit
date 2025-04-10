# Rutils/se2sce.R
#-------------------------------------------------------------------------------

# 转换工具：将 SummarizedExperiment 对象转换为 SingleCellExperiment 对象
#-------------------------------------------------------------------------------
#
# 背景介绍:
#   - SummarizedExperiment 是 Bioconductor 中用于存储基因组数据的通用数据结构，包含表达矩阵、行数据、列数据和元数据。
#   - SingleCellExperiment 是 SummarizedExperiment 的扩展，专为单细胞数据设计，支持存储降维结果、细胞元数据等。
#   - 本函数提供以下功能：
#     - 将 SummarizedExperiment 对象转换为 SingleCellExperiment 对象。
#     - 保留所有原始数据（表达矩阵、行数据、列数据、元数据）。
#
# 参数说明:
#   - se: SummarizedExperiment 对象（输入）
#
# 返回值:
#   - SingleCellExperiment 对象
#
# 依赖包:
#   - SummarizedExperiment (处理 SummarizedExperiment 对象)
#   - SingleCellExperiment (创建 SingleCellExperiment 对象)
#   - cli（格式化输出）

se2sce <- function(se) {

  # 加载必要的包
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("请先安装 SummarizedExperiment 包：BiocManager::install('SummarizedExperiment')", call. = FALSE)
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("请先安装 SingleCellExperiment 包：BiocManager::install('SingleCellExperiment')", call. = FALSE)
  }
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("请先安装 cli 包：install.packages('cli')", call. = FALSE)
  }
  # 检查输入参数
  if (!inherits(se, "SummarizedExperiment")) {
    stop("输入对象 se 必须是 SummarizedExperiment 类！", call. = FALSE)
  }

  if (length(SummarizedExperiment::assays(se)) == 0) {
    stop("输入的 SummarizedExperiment 对象 assays 为空，请检查输入。", call. = FALSE)
  }

  # 打印 assays 名称信息
  cli::cli_h3("检测到的表达矩阵 (assays):")
  assay_names <- SummarizedExperiment::assayNames(se)
  cli::cli_text("{.field 共 {length(assay_names)} 个 assay：}")
  cli::cli_li(assay_names)

  # 提取组件并转换
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = SummarizedExperiment::assays(se),
    rowData = SummarizedExperiment::rowData(se),
    colData = SummarizedExperiment::colData(se),
    metadata = SummarizedExperiment::metadata(se)
  )

  cli::cli_alert_success("已成功转换为 SingleCellExperiment 对象。")

  # 返回 SingleCellExperiment 对象
  return(sce)
}

#-------------------------------------------------------------------------------
# 示例用法: 将 SummarizedExperiment 对象转换为 SingleCellExperiment 对象
#-------------------------------------------------------------------------------

# # 创建一个简单的 SummarizedExperiment 对象
# counts <- matrix(rnorm(200), nrow = 10)
# col_data <- DataFrame(cell_id = paste0("cell_", 1:20))
# row_data <- DataFrame(gene_id = paste0("gene_", 1:10))
# se <- SummarizedExperiment(assays = list(counts = counts), colData = col_data, rowData = row_data)
#
# # 转换为 SingleCellExperiment 对象
# sce <- se2sce(se)
# print(sce)

#-------------------------------------------------------------------------------

