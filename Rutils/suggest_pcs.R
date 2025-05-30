# Rutils/suggest_pcs.R
#-------------------------------------------------------------------------------

# scFlowKit: Suggest Number of Principal Components for Downstream Analysis
#-------------------------------------------------------------------------------

# 主成分选择背景说明：
# - 本函数用于基于 PCA 降维结果推荐用于聚类与可视化的主成分数量。
# - 在 PCA 过程中，每个主成分（PC）解释原始数据的一部分变异；
# - 通过组合两种常用指标选择合理数量的 PCs，避免信息损失或噪声过多：
#   - 累计方差 > 90%，且当前 PC 贡献 < 5%
#   - 方差变化差值 > 0.1% 的最后一个 PC（即“肘部拐点”）

#-------------------------------------------------------------------------------
# suggest_pcs: 推荐主成分数量
#-------------------------------------------------------------------------------
# 参数:
#   seu: Seurat 对象，包含 PCA 降维结果
#   reduction: 使用的降维结果名称，默认为 "pca"
#   verbose: 是否输出建议过程信息，默认 TRUE
#
# 返回值:
#   - 建议使用的 PC 数量（整数值）
#-------------------------------------------------------------------------------

suggest_pcs <- function(seu, reduction = "pca", verbose = TRUE) {
  # ---------------- 参数检查 ----------------
  if (!inherits(seu, "Seurat")) {
    stop("参数 'seu' 必须是 Seurat 对象", call. = FALSE)
  }
  if (!reduction %in% names(seu@reductions)) {
    stop(glue::glue("Seurat 对象中未找到降维结果 '{reduction}'"), call. = FALSE)
  }

  # ---------------- 提取 PCA 信息 ----------------
  stdev <- seu[[reduction]]@stdev # stdev：每个 PC 的标准差，表示变异大小
  total_var <- sum(stdev)
  var_ratio <- (stdev / total_var) * 100 # var_ratio：每个 PC 的方差贡献比例（%）
  cum_var <- cumsum(var_ratio) # 前 k 个 PCs 的累计方差贡献（%）

  # ---------------- 指标 1 ----------------
  # - 累计解释度 > 90%，且当前 PC 解释度 < 5%
  pc_cutoff_1 <- which(cum_var > 90 & var_ratio < 5)[1] # 确保涵盖大部分变异，同时避免噪声

  # ---------------- 指标 2 ----------------
  # - 解释度下降 > 0.1% 的最后一个主成分（肘部）
  var_diff <- head(var_ratio, -1) - tail(var_ratio, -1) # 找到方差贡献下降速度显著减缓的 PC（肘部位置）
  pc_cutoff_2 <- if (any(var_diff > 0.1)) {
    max(which(var_diff > 0.1)) + 1
  } else {
    NA
  }

  # ---------------- 最终推荐 ----------------
  pcs_suggested <- min(pc_cutoff_1, pc_cutoff_2, na.rm = TRUE)

  # ---------------- 输出信息 ----------------
  if (verbose) {
    cli::cli_h2("📊 建议主成分数")
    cli::cli_text("🔹 指标 1：累计方差 > 90%，且单 PC < 5% → {pc_cutoff_1}")
    cli::cli_text("🔹 指标 2：解释度下降 > 0.1% 的最后一个 PC → {pc_cutoff_2}")
    cli::cli_alert_info("🎯 推荐使用前 {pcs_suggested} 个主成分用于后续分析")
  }

  return(pcs_suggested)
}

#-------------------------------------------------------------------------------
