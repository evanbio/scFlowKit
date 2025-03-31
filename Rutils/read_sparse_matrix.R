# scRNAutils/read_sparse_matrix.R
#-------------------------------------------------------------------------------

# 读取单细胞 RNA 测序的特征-条形码矩阵
#-------------------------------------------------------------------------------
#
# 背景介绍:
#   - 本函数用于读取 CellRanger 输出的特征-条形码矩阵文件夹，生成带有行名（特征）和列名（条形码）的稀疏矩阵。
#   - 输入文件夹需包含以下文件：
#     - barcodes.tsv.gz：条形码列表
#     - features.tsv.gz：特征信息（通常包括基因 ID 和名称）
#     - matrix.mtx.gz：稀疏矩阵数据
#   - 函数会验证文件存在性、数据完整性，并提供详细的日志反馈。
#
# 参数说明:
#   - folder_path: 包含特征-条形码矩阵文件的文件夹路径
#   - barcode_file: 条形码文件名（默认 "barcodes.tsv.gz"）
#   - feature_file: 特征文件名（默认 "features.tsv.gz"）
#   - matrix_file: 矩阵文件名（默认 "matrix.mtx.gz"）
#
# 返回值:
#   - 一个稀疏矩阵，行名为特征 ID（ensembl），列名为条形码
#
# 依赖包:
#   - Matrix (稀疏矩阵处理)
#   - data.table (高效文件读取)
#   - cli (命令行交互提示)

read_sparse_matrix <- function(folder_path,
                               barcode_file = "barcodes.tsv.gz",
                               feature_file = "features.tsv.gz",
                               matrix_file = "matrix.mtx.gz") {
  
  # 加载必要的包
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("请先安装 Matrix 包：install.packages('Matrix')", call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("请先安装 data.table 包：install.packages('data.table')", call. = FALSE)
  }
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("请先安装 cli 包：install.packages('cli')", call. = FALSE)
  }
  library(Matrix)
  library(data.table)
  library(cli)
  
  cli_h1("读取特征-条形码矩阵：{.path {folder_path}}")
  
  # 定义文件路径
  paths <- list(
    barcodes = file.path(folder_path, barcode_file),
    features = file.path(folder_path, feature_file),
    matrix = file.path(folder_path, matrix_file)
  )
  
  # 检查文件是否存在
  missing_files <- names(paths)[!vapply(paths, file.exists, logical(1))]
  if (length(missing_files) > 0) {
    cli_alert_danger("以下文件缺失：{paste(missing_files, collapse = ', ')}。请检查路径：{.path {folder_path}}")
    return(invisible(NULL))
  }
  
  # 读取条形码文件
  cli_text("读取条形码文件：{.path {paths$barcodes}}")
  barcodes <- tryCatch(
    data.table::fread(paths$barcodes, header = FALSE, col.names = "barcode"),
    error = function(e) {
      cli_alert_danger("读取条形码文件失败：{e$message}")
      return(NULL)
    }
  )
  if (is.null(barcodes) || nrow(barcodes) == 0) {
    cli_alert_danger("条形码文件为空或读取失败。")
    return(invisible(NULL))
  }
  
  # 读取特征文件
  cli_text("读取特征文件：{.path {paths$features}}")
  features <- tryCatch(
    data.table::fread(paths$features, header = FALSE, 
                      col.names = c("ensembl", "symbol", "type")),
    error = function(e) {
      cli_alert_danger("读取特征文件失败：{e$message}")
      return(NULL)
    }
  )
  if (is.null(features) || nrow(features) == 0) {
    cli_alert_danger("特征文件为空或读取失败。")
    return(invisible(NULL))
  }
  
  # 读取矩阵文件
  cli_text("读取矩阵文件：{.path {paths$matrix}}")
  mat <- tryCatch(
    Matrix::readMM(paths$matrix),
    error = function(e) {
      cli_alert_danger("读取矩阵文件失败：{e$message}")
      return(NULL)
    }
  )
  if (is.null(mat)) {
    cli_alert_danger("矩阵文件读取失败。")
    return(invisible(NULL))
  }
  
  # 验证矩阵维度
  if (nrow(mat) != nrow(features) || ncol(mat) != nrow(barcodes)) {
    cli_alert_danger("矩阵维度 ({nrow(mat)} x {ncol(mat)}) 与特征 ({nrow(features)}) 或条形码 ({nrow(barcodes)}) 不匹配。")
    return(invisible(NULL))
  }
  
  # 设置行名和列名
  rownames(mat) <- features$ensembl
  colnames(mat) <- barcodes$barcode
  
  cli_alert_success(sprintf("成功导入稀疏矩阵：%d 列（条形码） x %d 行（特征）。", 
                            ncol(mat), nrow(mat)))
  
  return(mat)
}

#-------------------------------------------------------------------------------
# 示例用法: 读取 CellRanger 输出的特征-条形码矩阵
#-------------------------------------------------------------------------------

# # 示例 1: 使用默认文件名读取矩阵
# mat <- read_sparse_matrix("/path/to/cellranger/output")
# if (!is.null(mat)) {
#   cat("矩阵维度：", dim(mat), "\n")
#   cat("前几个条形码：", head(colnames(mat)), "\n")
#   cat("前几个特征：", head(rownames(mat)), "\n")
# }
#
# # 示例 2: 指定自定义文件名
# mat <- read_sparse_matrix(
#   folder_path = "/path/to/custom/output",
#   barcode_file = "custom_barcodes.tsv.gz",
#   feature_file = "custom_features.tsv.gz",
#   matrix_file = "custom_matrix.mtx.gz"
# )
# if (!is.null(mat)) {
#   cat("成功读取自定义矩阵，维度：", dim(mat), "\n")
# }
