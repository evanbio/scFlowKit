# Rutils/explore_seurat.R
#-------------------------------------------------------------------------------
# 工具函数：探索 Seurat 对象结构，支持多种访问方式（v4+ 兼容）
#-------------------------------------------------------------------------------
#
# 背景介绍:
#   - 在单细胞分析中，了解 Seurat 对象内部结构有助于数据理解与调试。
#   - 本函数支持结构概览、表达矩阵提取、元数据查看、默认 assay、特征名和数据层查看等操作。
#   - 同时保留不同方法（推荐 / 兼容 / 旧版）供学习参考。
#
# 参数说明:
#   - seu          : Seurat 对象
#   - explore_mode : 是否启用探索模式（默认 TRUE）
#
# 返回值:
#   - 无返回值，打印探索信息（仅作交互输出使用）

explore_seurat <- function(seu, explore_mode = TRUE) {
  #---------------------------
  # 显式加载必要包
  #---------------------------
  if (!requireNamespace("cli", quietly = TRUE)) stop("请安装 cli 包")
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("请安装 Seurat 包")
  if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("请安装 SeuratObject 包")

  if (!explore_mode || !exists("seu")) return(invisible(NULL))

  cli::cli_h1("🔍 Seurat 对象结构探索")

  #------------------------------------------------------------
  # 1. 查看整体结构
  #------------------------------------------------------------
  cli::cli_h2("📦 Seurat 对象结构")
  utils::str(seu)

  #------------------------------------------------------------
  # 2. 表达矩阵（counts 层）
  #------------------------------------------------------------
  cli::cli_h2("🧬 表达矩阵（counts）")

  # 方法 1：底层语法（不推荐）
  # counts1 <- seu[["RNA"]]@counts
  # print(counts1[1:5, 1:5])

  # 方法 2：推荐方式
  counts <- Seurat::GetAssayData(seu, assay = "RNA", slot = "counts")
  print(counts[1:5, 1:5])

  #------------------------------------------------------------
  # 3. 提取 metadata（细胞元数据）
  #------------------------------------------------------------
  cli::cli_h2("📋 细胞元数据")

  # 方法 1：旧版本方式
  # metadata1 <- seu@meta.data

  # 方法 2：colData 接口（兼容 SCE）
  # metadata2 <- colData(seu)

  # 方法 3：推荐写法
  metadata <- seu[[]]
  utils::head(metadata)
  cli::cli_text("📌 列名：{paste(colnames(metadata), collapse = ', ')}")

  #------------------------------------------------------------
  # 4. 默认激活的 assay
  #------------------------------------------------------------
  cli::cli_h2("🔧 当前默认 assay")
  active_assay <- Seurat::DefaultAssay(seu)  # nolint
  cli::cli_text("默认 assay：{active_assay}")

  #------------------------------------------------------------
  # 5. 查看基因名与细胞名
  #------------------------------------------------------------
  cli::cli_h2("🧬 基因名与细胞名")

  # 方法 1：底层访问
  # genes1 <- rownames(seu@assays$RNA)
  # cells1 <- colnames(seu@assays$RNA)

  # 方法 2：推荐方式
  genes <- rownames(seu) # nolint
  cells <- colnames(seu) # nolint

  cli::cli_text("🔹 Genes: {paste(head(genes), collapse = ', ')}")
  cli::cli_text("🔹 Cells: {paste(head(cells), collapse = ', ')}")

  #------------------------------------------------------------
  # 6. 查看数据层（slots / layers）
  #------------------------------------------------------------
  cli::cli_h2("📚 当前数据层（slots/layers）")
  cli::cli_text("层：{paste(SeuratObject::Layers(seu), collapse = ', ')}")
}
