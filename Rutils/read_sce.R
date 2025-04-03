# scRNAutils/read_sce.R
#-------------------------------------------------------------------------------

# 读取单细胞 RNA 测序数据为 SingleCellExperiment 对象
#-------------------------------------------------------------------------------
#
# 背景介绍:
#   - 本函数读取单细胞 RNA 测序数据，支持 CellRanger 的 .mtx 文件夹、10x Genomics 的 .h5 文件、AnnData 的 .h5ad 文件和 Loom 的 .loom 文件。
#   - 对于 .mtx 文件夹，可选择手动读取（type = "mtx"）或 Seurat::Read10X（type = "10x"）。
#
# 参数说明:
#   - input_path: 输入路径，可以是文件夹（.mtx 模式）或文件（.h5/.h5ad/.loom 模式）
#   - type: 数据类型，默认根据 input_path 推断，可选 "mtx"、"10x"、"h5"、"h5ad"、"loom"
#   - barcode_file: 条形码文件名（默认 "barcodes.tsv.gz"，仅 type = "mtx"）
#   - feature_file: 特征文件名（默认 "features.tsv.gz"，仅 type = "mtx"）
#   - matrix_file: 矩阵文件名（默认 "matrix.mtx.gz"，仅 type = "mtx"）
#   - gene.column: Read10X 中特征文件的基因列（默认 2，即 Symbol，仅 type = "10x"）
#   - cell.column: Read10X 中条形码文件的细胞列（默认 1，仅 type = "10x"）
#   - unique.features: 是否确保基因名唯一（默认 TRUE，仅 type = "10x"）
#   - strip.suffix: 是否移除条形码后缀（如 -1，默认 FALSE，仅 type = "10x"）
#
# 返回值:
#   - SingleCellExperiment 对象，若失败则返回 NULL
#
# 依赖包:
#   - SingleCellExperiment (SCE 对象)
#   - Matrix (稀疏矩阵处理)
#   - data.table (高效文件读取，type = "mtx")
#   - rhdf5 (HDF5 文件读取，type = "h5")
#   - zellkonverter (AnnData 文件读取，type = "h5ad")
#   - LoomExperiment (Loom 文件读取，type = "loom")
#   - Seurat (Read10X 函数，type = "10x")
#   - cli (命令行交互提示)

read_sce <- function(input_path,
                     type = NULL,
                     barcode_file = "barcodes.tsv.gz",
                     feature_file = "features.tsv.gz",
                     matrix_file = "matrix.mtx.gz",
                     gene.column = 2,
                     cell.column = 1,
                     unique.features = TRUE,
                     strip.suffix = FALSE) {
  
  # -------------------- 加载依赖包 --------------------
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("请先安装 SingleCellExperiment 包：BiocManager::install('SingleCellExperiment')", call. = FALSE)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("请先安装 Matrix 包：install.packages('Matrix')", call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("请先安装 data.table 包：install.packages('data.table')", call. = FALSE)
  }
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("请先安装 rhdf5 包：BiocManager::install('rhdf5')", call. = FALSE)
  }
  if (!requireNamespace("zellkonverter", quietly = TRUE)) {
    stop("请先安装 zellkonverter 包：BiocManager::install('zellkonverter')", call. = FALSE)
  }
  if (!requireNamespace("LoomExperiment", quietly = TRUE)) {
    stop("请先安装 LoomExperiment 包：BiocManager::install('LoomExperiment')", call. = FALSE)
  }
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("请先安装 cli 包：install.packages('cli')", call. = FALSE)
  }

  library(SingleCellExperiment)
  library(Matrix)
  library(data.table)
  library(rhdf5)
  library(zellkonverter)
  library(LoomExperiment)
  library(cli)
  
  cli_h1("读取单细胞数据：{.path {input_path}}")
  
  # -------------------- 自动推断数据类型 --------------------
  if (is.null(type)) {
    if (grepl("\\.h5ad$", input_path)) {
      type <- "h5ad"
    } else if (grepl("\\.h5$", input_path)) {
      type <- "h5"
    } else if (grepl("\\.loom$", input_path)) {
      type <- "loom"
    } else {
      type <- "mtx"  # 默认文件夹为手动读取
    }
  }
  cli_text("数据类型：{type}")
  
  # 根据 type 读取数据
  # -------------------- 读取 .h5ad 文件 --------------------
  if (type == "h5ad") {
    cli_text("使用 zellkonverter 读取 .h5ad 文件")
    sce <- tryCatch(
      zellkonverter::readH5AD(input_path),
      error = function(e) {
        cli_alert_danger("读取 .h5ad 文件失败：{e$message}")
        return(NULL)
      }
    )
    if (is.null(sce)) return(invisible(NULL))

  # -------------------- 读取 .h5 文件（10X） --------------------  
  } else if (type == "h5") {
    cli_text("使用 rhdf5 读取 .h5 文件")
    data <- tryCatch(h5read(input_path, "matrix/data"), error = function(e) NULL)
    indices <- tryCatch(h5read(input_path, "matrix/indices"), error = function(e) NULL)
    indptr <- tryCatch(h5read(input_path, "matrix/indptr"), error = function(e) NULL)
    shape <- tryCatch(h5read(input_path, "matrix/shape"), error = function(e) NULL)
    
    if (is.null(data) || is.null(indices) || is.null(indptr) || is.null(shape)) {
      cli_alert_danger("读取 .h5 矩阵组件失败，请检查文件结构。")
      return(invisible(NULL))
    }
    
    mat <- sparseMatrix(
      i = rep(seq_along(indptr[-length(indptr)]), diff(indptr)),
      j = indices + 1,
      x = data,
      dims = shape
    )
    
    barcodes <- tryCatch(h5read(input_path, "matrix/barcodes"), error = function(e) NULL)
    features <- tryCatch({
      data.frame(
        ensembl = h5read(input_path, "matrix/features/id"),
        symbol = h5read(input_path, "matrix/features/name"),
        type = h5read(input_path, "matrix/features/feature_type")
      )
    }, error = function(e) NULL)
    
    if (is.null(barcodes) || is.null(features)) {
      cli_alert_danger("读取 .h5 元数据失败，请检查文件结构。")
      return(invisible(NULL))
    }
    
    # 设置矩阵行名和列名
    rownames(mat) <- features$ensembl
    colnames(mat) <- barcodes
    
    sce <- SingleCellExperiment(
      assays = list(counts = mat),
      colData = data.frame(barcode = barcodes),
      rowData = features
    )

    # ✅ 检查行名完整性
    if (any(is.na(rownames(mat)) | rownames(mat) == "")) {
      cli_alert_warning("⚠️ 检测到部分基因行名为空，请检查基因注释信息是否正确")
    }

  # -------------------- 读取 .loom 文件 --------------------  
  } else if (type == "loom") {
    cli_text("使用 LoomExperiment 读取 .loom 文件")
    sce <- tryCatch(
      import(input_path, type = "LoomExperiment"),
      error = function(e) {
        cli_alert_danger("读取 .loom 文件失败：{e$message}")
        return(NULL)
      }
    )
    if (is.null(sce)) return(invisible(NULL))
    
    if (!"counts" %in% SummarizedExperiment::assayNames(sce)) {
      SummarizedExperiment::assayNames(sce)[1] <- "counts"
    }

  # -------------------- 读取 10X 文件夹（type = '10x'） --------------------  
  } else if (type == "10x") {
    cli_text("使用 Seurat::Read10X 读取 CellRanger 数据")
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("请先安装 Seurat 包：install.packages('Seurat')", call. = FALSE)
    }

    library(Seurat)
    
    mat <- tryCatch(
      Read10X(
        data.dir = input_path,
        gene.column = gene.column,
        cell.column = cell.column,
        unique.features = unique.features,
        strip.suffix = strip.suffix
      ),
      error = function(e) {
        cli_alert_danger("Read10X 读取失败：{e$message}")
        return(NULL)
      }
    )
    
    if (is.null(mat)) return(invisible(NULL))
    
    barcodes <- colnames(mat)
    features <- data.frame(rowname = rownames(mat))
    if (gene.column == 1) {
      colnames(features) <- "ensembl"
    } else {
      colnames(features) <- "symbol"
    }
    
    feature_path <- file.path(input_path, feature_file)
    if (file.exists(feature_path)) {
      features_full <- tryCatch(
        data.table::fread(feature_path, header = FALSE, 
                          col.names = c("ensembl", "symbol", "type")),
        error = function(e) NULL
      )
      if (!is.null(features_full)) features <- features_full
    }
    
    # 设置矩阵行名和列名（Read10X 已提供，但确保一致）
    rownames(mat) <- features$ensembl
    colnames(mat) <- barcodes

    # ✅ 检查行名完整性
    if (any(is.na(rownames(mat)) | rownames(mat) == "")) {
      cli_alert_warning("⚠️ 检测到部分基因行名为空，请检查基因注释信息是否正确")
    }

    sce <- SingleCellExperiment(
      assays = list(counts = mat),
      colData = data.frame(barcode = barcodes),
      rowData = features
    )

  # -------------------- 读取手动 mtx 文件夹（type = 'mtx'） --------------------  
  } else if (type == "mtx") {
    cli_text("手动读取 CellRanger 数据")
    paths <- list(
      barcodes = file.path(input_path, barcode_file),
      features = file.path(input_path, feature_file),
      matrix = file.path(input_path, matrix_file)
    )
    
    missing_files <- names(paths)[!vapply(paths, file.exists, logical(1))]
    if (length(missing_files) > 0) {
      cli_alert_danger("以下文件缺失：{paste(missing_files, collapse = ', ')}。")
      return(invisible(NULL))
    }
    
    cli_text("读取条形码文件：{.path {paths$barcodes}}")
    barcodes <- tryCatch(
      data.table::fread(paths$barcodes, header = FALSE, col.names = "barcode"),
      error = function(e) { cli_alert_danger("读取失败：{e$message}"); NULL }
    )
    if (is.null(barcodes) || nrow(barcodes) == 0) return(invisible(NULL))
    
    cli_text("读取特征文件：{.path {paths$features}}")
    features <- tryCatch(
      data.table::fread(paths$features, header = FALSE, 
                        col.names = c("ensembl", "symbol", "type")),
      error = function(e) { cli_alert_danger("读取失败：{e$message}"); NULL }
    )
    if (is.null(features) || nrow(features) == 0) return(invisible(NULL))
    
    cli_text("读取矩阵文件：{.path {paths$matrix}}")
    mat <- tryCatch(
      Matrix::readMM(paths$matrix),
      error = function(e) { cli_alert_danger("读取失败：{e$message}"); NULL }
    )
    if (is.null(mat)) return(invisible(NULL))
    
    if (nrow(mat) != nrow(features) || ncol(mat) != nrow(barcodes)) {
      cli_alert_danger("矩阵维度不匹配：{nrow(mat)} x {ncol(mat)} vs {nrow(features)} x {nrow(barcodes)}")
      return(invisible(NULL))
    }
    
    # 转换为 dgCMatrix
    mat <- as(mat, "dgCMatrix")
    
    # 设置矩阵行名和列名
    rownames(mat) <- features$ensembl
    colnames(mat) <- barcodes$barcode

    # ✅ 检查行名完整性
    if (any(is.na(rownames(mat)) | rownames(mat) == "")) {
      cli_alert_warning("⚠️ 检测到部分基因行名为空，请检查基因注释信息是否正确")
    }

    sce <- SingleCellExperiment(
      assays = list(counts = mat),
      colData = data.frame(barcode = barcodes$barcode),
      rowData = features
    )
    
  } else {
    cli_alert_danger("不支持的 type 参数：{type}，可选值为 'mtx', '10x', 'h5', 'h5ad', 'loom'")
    return(invisible(NULL))
  }
  
  # -------------------- 返回结果 --------------------
  cli_alert_success(sprintf("成功创建 SingleCellExperiment 对象：%d 细胞 x %d 特征", 
                            ncol(sce), nrow(sce)))
  
  return(sce)
}

#-------------------------------------------------------------------------------
# 示例用法
#-------------------------------------------------------------------------------

# # .mtx 文件夹（手动读取）
# sce <- read_sce("/path/to/cellranger/output")
#
# # .mtx 文件夹（用 Read10X，默认 Symbol）
# sce <- read_sce("/path/to/cellranger/output", type = "10x")
#
# # .mtx 文件夹（用 Read10X，用 Ensembl，去除后缀）
# sce <- read_sce("/path/to/cellranger/output", type = "10x", gene.column = 1, strip.suffix = TRUE)
#
# # .h5 文件
# sce <- read_sce("/path/to/filtered_feature_bc_matrix.h5")
#
# # .h5ad 文件
# sce <- read_sce("/path/to/data.h5ad")
#
# # .loom 文件
# sce <- read_sce("/path/to/data.loom")
#
# # 检查结果
# if (!is.null(sce)) {
#   print(sce)              # 查看 SCE 结构
#   head(rownames(sce))     # 检查行名
#   head(colnames(sce))     # 检查列名
#   head(colData(sce))      # 检查 colData
#   head(rowData(sce))      # 检查 rowData
# }
#-------------------------------------------------------------------------------