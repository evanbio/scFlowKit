# Rutils/remove_doublets.R
#-------------------------------------------------------------------------------

# scFlowKit: Detect and Remove Doublets Using DoubletFinder
#-------------------------------------------------------------------------------

# 双细胞检测背景介绍
# - 在单细胞 RNA-seq 数据中，双细胞（doublets）是指两个或更多细胞被错误识别为一个细胞的情况。
# - 双细胞会干扰下游分析（如聚类、差异表达分析），因此需要在数据预处理阶段检测和去除。
# - DoubletFinder 是一种常用的双细胞检测工具，基于以下原理：
#   - 通过 PCA 降维后的数据，计算每个细胞的邻居关系。
#   - 模拟生成人工双细胞（artificial doublets），比较真实细胞与人工双细胞的相似性。
#   - 使用 pANN（proportion of artificial nearest neighbors）指标，判断细胞是否为双细胞。
# - DoubletFinder 的关键参数：
#   - pN：人工双细胞比例（通常 0.25）。
#   - pK：邻居数量比例（通过参数扫描确定最佳值）。
#   - nExp：预期双细胞数量（基于双细胞比例和细胞总数计算）。
# - 本函数的主要步骤：
#   1. 参数扫描：确定最佳 pK 值。
#   2. 双细胞检测：运行 DoubletFinder，生成双细胞标签。
#   3. 双细胞去除：根据标签过滤掉双细胞。
# - 注意事项：
#   - 需要在 PCA 降维后运行（RunPCA）。
#   - 双细胞比例（doublet_rate）需根据实验设计调整（例如 10x Genomics 数据通常为 0.05-0.08）。
#   - 本函数默认使用 RNA assay 数据（sct = FALSE），如使用 SCTransform 数据需设置 sct = TRUE。

# remove_doublets: 使用 DoubletFinder 检测和去除双细胞
# 参数:
#   seu: Seurat 对象，包含单细胞 RNA-seq 数据（需已完成 PCA 降维）
#   PCs: 用于双细胞检测的主成分范围（默认 1:20）
#   doublet_rate: 假设的双细胞比例（默认 0.08，可根据实验设计调整）
#   pN: DoubletFinder 参数，人工双细胞比例（默认 0.25）
#   sct: 是否使用 SCTransform 数据（默认 FALSE，使用 RNA assay）
# 
# 返回:
#   过滤双细胞后的 Seurat 对象
remove_doublets <- function(seu, PCs = 1:20, doublet_rate = 0.08, pN = 0.25, sct = FALSE) {

  cli::cli_h2("检测并去除双细胞")

  # ------------------------- 参数检查 -----------------------
  # 验证输入参数是否为 Seurat 对象
  if (!inherits(seu, "Seurat")) {
    stop("参数 'seu' 必须为 Seurat 对象！", call. = FALSE)
  }

  # 验证是否已完成 PCA 降维
  if (!"pca" %in% names(seu@reductions)) {
    stop("Seurat 对象尚未完成 PCA 降维，请先运行 RunPCA！", call. = FALSE)
  }

  # 验证 PCs 参数是否为数值向量
  if (!is.numeric(PCs) || length(PCs) < 2) {
    stop("参数 'PCs' 必须为数值向量，且长度至少为 2！", call. = FALSE)
  }

  # 验证 PCs 是否在有效范围内
  max_pc <- ncol(seu@reductions$pca@cell.embeddings)
  if (max(PCs) > max_pc) {
    stop("PCs 参数超出 PCA 主成分数量（最大为 ", max_pc, "）！", call. = FALSE)
  }

  # 验证 doublet_rate 参数
  if (!is.numeric(doublet_rate) || doublet_rate <= 0 || doublet_rate >= 1) {
    stop("参数 'doublet_rate' 必须为 0 到 1 之间的数值！", call. = FALSE)
  }

  # 验证 pN 参数
  if (!is.numeric(pN) || pN <= 0 || pN >= 1) {
    stop("参数 'pN' 必须为 0 到 1 之间的数值！", call. = FALSE)
  }

  # 验证 sct 参数
  if (!is.logical(sct)) {
    stop("参数 'sct' 必须为逻辑值！", call. = FALSE)
  }

  # 加载 DoubletFinder 包
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    stop("请先安装 DoubletFinder 包：install.packages('DoubletFinder')", call. = FALSE)
  }
  suppressPackageStartupMessages(library(DoubletFinder))

  # ---------------- 参数扫描确定最佳 pK ----------------
  cli::cli_text("正在扫描 pK 参数...")
  sweep.res <- paramSweep(seu, PCs = PCs, sct = sct, num.cores = 1)  # 单核运行，兼容 Windows
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)

  # 选择最佳 pK 值（通常是最大 BCMVN 对应的 pK）
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  cli::cli_text("最佳 pK 值为：{pK}")

  # ---------------- 计算预期双细胞数量 ----------------
  nExp_poi <- round(doublet_rate * ncol(seu) * ncol(seu) / 10000)
  cli::cli_text("预期双细胞数量：{nExp_poi}")

  # ---------------- 运行 DoubletFinder ----------------
  cli::cli_text("运行 DoubletFinder 检测双细胞...")
  seu <- doubletFinder(seu, 
                       PCs = PCs, 
                       pN = pN, 
                       pK = pK, 
                       nExp = nExp_poi, 
                       sct = sct)

  # 关闭所有未使用的连接，避免警告
  closeAllConnections()

  # ---------------- 统计并去除双细胞 ----------------

  # DoubletFinder 添加的列名通常为 "DF.classifications_pN_pK_nExp"
  df_col <- paste0("DF.classifications_", pN, "_", pK, "_", nExp_poi)

  # 统计实际双细胞数量
  doublet_counts <- table(seu@meta.data[[df_col]])
  cli::cli_text("实际双细胞检测统计：")
  print(doublet_counts)

  actual_doublets <- if ("Doublet" %in% names(doublet_counts)) doublet_counts["Doublet"] else 0
  cli::cli_text("实际检测到的双细胞数量：{actual_doublets}（预期：{nExp_poi}）")

  # 过滤掉双细胞
  cli::cli_text("正在去除双细胞...")
  seu <- subset(seu, subset = !!sym(df_col) == "Singlet")

  cli::cli_alert_success("双细胞去除完成！")

  return(seu)
}

#-------------------------------------------------------------------------------