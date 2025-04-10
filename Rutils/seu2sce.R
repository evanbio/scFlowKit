# Rutils/seu2sce.R
#-------------------------------------------------------------------------------

# 转换工具：将 Seurat 对象转换为 SingleCellExperiment 对象
#-------------------------------------------------------------------------------
#
# 背景介绍:
#   - Seurat 是单细胞 RNA-seq 分析的常用 R 包，其核心数据结构是 Seurat 对象，包含表达矩阵、细胞元数据、降维结果等。
#   - SingleCellExperiment 是 Bioconductor 中用于单细胞数据的通用数据结构，支持存储表达矩阵、元数据、降维结果等。
#   - 本函数提供以下功能：
#     - 基于 as.SingleCellExperiment 将 Seurat 对象转换为 SingleCellExperiment 对象。
#     - 可选择是否提取降维结果和高变基因标记。
#
# 参数说明:
#   - seu: Seurat 对象（输入）
#   - reductions: 逻辑值，是否提取降维结果（默认 FALSE）
#   - assay：设置 DefaultAssay，默认 "RNA"
#   - gene_metadata: 逻辑值，是否提取基因注释信息（默认 FALSE）
#
# 返回值:
#   - SingleCellExperiment 对象
#
# 依赖包:
#   - Seurat (处理 Seurat 对象)
#   - SingleCellExperiment (创建 SingleCellExperiment 对象)

seu2sce <- function(seu,
                    assay = "RNA",
                    reductions = FALSE,
                    gene_metadata = FALSE) {

  #-------------------- 依赖检查 --------------------
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("请先安装 Seurat 包：install.packages('Seurat')", call. = FALSE)
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("请先安装 SingleCellExperiment 包：BiocManager::install('SingleCellExperiment')", call. = FALSE)
  }
  if (!requireNamespace("cli", quietly = TRUE)) {
    install.packages("cli", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  }


  #-------------------- 参数检查 --------------------
  if (!inherits(seu, "Seurat")) {
    stop("输入对象 seu 必须是 Seurat 类！", call. = FALSE)
  }

  if (!assay %in% names(seu@assays)) {
    stop("Seurat 对象中不包含指定的 assay: ", assay, call. = FALSE)
  }

  #-------------------- 设置 assay --------------------
  Seurat::DefaultAssay(seu) <- assay
  cli::cli_alert_info("DefaultAssay 已设置为：'{assay}'")

  #-------------------- 主转换过程 --------------------
  # 使用 as.SingleCellExperiment 进行转换
  sce <- Seurat::as.SingleCellExperiment(seu)

  # 降维结果（如果需要）
  if (reductions) {
    for (reduction in names(seu@reductions)) {
      red_name <- paste0("seurat_", reduction)
      reducedDim(sce, type = red_name) <- Seurat::Embeddings(seu, reduction = reduction)
    }
  }

  # 提取高变基因标记（如果需要）
  if (gene_metadata) {
    gene_info <- data.frame(
      gene_id = rownames(sce),
      variable_feature = rownames(sce) %in% Seurat::VariableFeatures(seu),
      row.names = rownames(sce)
    )
    SingleCellExperiment::rowData(sce) <- gene_info
  }

  # 返回 SingleCellExperiment 对象
  return(sce)
}

#-------------------------------------------------------------------------------
# 示例用法: 将 Seurat 对象转换为 SingleCellExperiment 对象
#-------------------------------------------------------------------------------

# # 示例 Seurat 对象
# sce <- seu2sce(seu)
#
# # 转换为 SingleCellExperiment 对象（提取降维和基因注释信息）
# sce_with_extras <- seu2sce(seu, reductions = TRUE, gene_metadata = TRUE)
# print(sce_with_extras)

#-------------------------------------------------------------------------------
