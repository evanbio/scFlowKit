# Rutils/sce2seu.R
#-------------------------------------------------------------------------------
# 转换工具：将 SingleCellExperiment 对象转换为 Seurat 对象
#-------------------------------------------------------------------------------
#
# 背景介绍:
#   - Bioconductor 提供的单细胞数据往往是 SingleCellExperiment 格式（SCE）。
#   - 若要在 Seurat 流程中分析这些数据，需先转换为 Seurat 对象。
#
# 参数说明:
#   - sce: SCE 对象（SingleCellExperiment）
#   - counts_assay: 要作为 counts 的层名称（默认 "counts"）
#   - data_assay: 要作为 data 的层（可选，默认 NULL 表示不设 logcounts）
#   - project: Seurat 对象的 project.name（默认 "sce_import"）
#
# 返回值:
#   - Seurat 对象（seu）

sce2seu <- function(sce,
                    counts_assay = "counts",
                    project = "sce_import") {

  library(Seurat)
  library(SingleCellExperiment)

  # 检查对象类型
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("输入对象 sce 必须是 SingleCellExperiment 类！", call. = FALSE)
  }

  # 检查 counts 层是否存在
  if (!counts_assay %in% assayNames(sce)) {
    stop("指定的 counts_assay 不存在：", counts_assay, call. = FALSE)
  }

  # 提取 counts 数据
  counts <- assay(sce, counts_assay)

  # 检查行名是否存在
  if (is.null(rownames(counts)) || any(rownames(counts) == "")) {
    stop("counts matrix 缺少有效的行名（基因名），请检查！", call. = FALSE)
  }

  # 创建 Seurat 对象
  seu <- CreateSeuratObject(counts = counts, project = project)

  # 添加细胞元数据
  seu <- AddMetaData(seu, metadata = as.data.frame(colData(sce)))

  cli::cli_alert_success("✅ SCE 对象已成功转换为 Seurat 对象！")
  cli::cli_text("📦 包含 {ncol(seu)} 个细胞，{nrow(seu)} 个基因。")

  return(seu)
}

#-------------------------------------------------------------------------------
# 示例用法
#-------------------------------------------------------------------------------

# sce <- readRDS("your_sce_data.rds")  # 读取 SCE 对象
# seu <- sce2seu(sce)