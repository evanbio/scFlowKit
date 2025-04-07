# Rutils/calculate_qc_metrics.R
#-------------------------------------------------------------------------------
# scFlowKit: Calculate QC Metrics for Single-Cell RNA-seq Data
#-------------------------------------------------------------------------------
# 
# 📌 背景说明：
# 在单细胞 RNA 测序分析中，质控指标（QC metrics）是初步评估细胞质量的核心步骤，
# 可用于识别低质量细胞、污染、双细胞等问题，确保下游分析的可靠性。
# 
# 常见质控指标（QC Metrics）介绍
# - 在单细胞 RNA-seq 分析中，质控指标用于评估细胞和数据的质量，帮助过滤低质量细胞。
# - 常见的质控指标包括：
#   - nCount_RNA：每个细胞的总 UMI 计数（unique molecular identifiers），反映测序深度。
#     - 意义：UMI 计数过低可能表示细胞质量差或测序深度不足；过高可能表示双细胞（doublets）。
#   - nFeature_RNA：每个细胞的基因数（unique genes detected），反映细胞的基因表达多样性。
#     - 意义：基因数过低可能表示细胞死亡或低质量；过高可能表示双细胞或污染。
#   - percent_mito：线粒体基因比例，反映细胞的代谢状态和存活状态。
#     - 意义：线粒体基因比例过高通常表示细胞凋亡或应激（例如 >10% 可能需要过滤）。
#   - percent_hb：红细胞基因比例（可选），反映红细胞污染。
#     - 意义：红细胞基因比例过高可能表示样本污染（例如 PBMC 数据中不应有大量红细胞）。
#   - percent_ribo：核糖体基因比例（可选），反映细胞的翻译活性。
#     - 意义：核糖体基因比例过高可能表示细胞处于高代谢状态，需结合实验背景判断。
#   - log10_ratio_features_to_umi：log10(nFeature_RNA) / log10(nCount_RNA)，反映基因表达复杂性。
#     - 意义：值过低可能表示 RNA 降解或低质量细胞（参考 HBC 设置）。
#------------------------------------------------------------------------------- 
# - Seurat 对象创建时的默认指标：
#   - 在创建 Seurat 对象时（CreateSeuratObject），会自动计算 nCount_RNA 和 nFeature_RNA。
#   - nCount_RNA：基于 counts 矩阵的列和（colSums(counts)）。
#   - nFeature_RNA：基于 counts 矩阵中非零基因的列和（colSums(counts > 0)）。

# - 本函数计算的指标：
#   - nCount_RNA 和 nFeature_RNA（如果未计算或者误删除，则重新计算）。
#   - percent_mito：线粒体基因比例（基于 "^MT-" 开头的基因）。
#   - percent_hb：红细胞基因比例（可选，基于用户提供或默认的红细胞基因列表）。
#   - percent_ribo：核糖体基因比例（可选，基于用户提供或默认的核糖体基因列表）。
#   - log10_ratio_features_to_umi：log10(nFeature_RNA) / log10(nCount_RNA)（参考 HBC 设置）。

# 主函数：计算质控指标
# 参数:
#   seu: Seurat 对象，包含单细胞 RNA-seq 数据
#   calculate_hb: 是否计算红细胞基因比例，默认 FALSE
#   calculate_ribo: 是否计算核糖体基因比例，默认 FALSE
#   hb_genes: 红细胞基因列表，默认 NULL（使用内置默认基因集）
#   ribo_genes: 核糖体基因列表，默认 NULL（使用内置默认基因集）

# 返回值
# - 添加了质控指标列的 Seurat 对象（修改原始 meta.data）

calculate_qc_metrics <- function(sce, 
                                 calculate_hb = FALSE, 
                                 calculate_ribo = FALSE,
                                 hb_genes = NULL,
                                 ribo_genes = NULL) {
  #-------------------------------------------------------------------------------
  # ✅ 参数校验：确保输入为合法 Seurat 对象，并检查选项参数
  #-------------------------------------------------------------------------------

  # 核查对象类型（必须为 Seurat 对象）
  if (!inherits(seu, "Seurat")) {
    stop("参数 'seu' 必须为 Seurat 对象。", call. = FALSE)
  }

  # 验证是否为逻辑值（布尔）
  if (!is.logical(calculate_hb)) {
    stop("参数 'calculate_hb' 必须为逻辑值（TRUE 或 FALSE）。", call. = FALSE)
  }
  if (!is.logical(calculate_ribo)) {
    stop("参数 'calculate_ribo' 必须为逻辑值（TRUE 或 FALSE）。", call. = FALSE)
  }

  # 自定义基因集校验（如指定，必须为非空字符向量）
  if (!is.null(hb_genes) && (!is.character(hb_genes) || length(hb_genes) == 0)) {
    stop("参数 'hb_genes' 必须为非空字符向量（或使用默认值 NULL）。", call. = FALSE)
  }
  if (!is.null(ribo_genes) && (!is.character(ribo_genes) || length(ribo_genes) == 0)) {
    stop("参数 'ribo_genes' 必须为非空字符向量（或使用默认值 NULL）。", call. = FALSE)
  }

  cli::cli_h1("🧪 计算质控指标")
  
  #-------------------------------------------------------------------------------
  # 🚧 检查并补充默认的计数指标（如丢失）
  #-------------------------------------------------------------------------------
  if (!"nCount_RNA" %in% colnames(seu@meta.data)) {
    cli::cli_alert_info("重新计算 nCount_RNA（每个细胞的总计数）...")
    seu@meta.data$nCount_RNA <- Matrix::colSums(seu[["RNA"]]$counts)
  } else {
    cli::cli_alert_success("已检测到 nCount_RNA，跳过重算。")
  }

  if (!"nFeature_RNA" %in% colnames(seu@meta.data)) {
    cli::cli_alert_info("重新计算 nFeature_RNA（每个细胞的表达基因数）...")
    seu@meta.data$nFeature_RNA <- Matrix::colSums(seu[["RNA"]]$counts > 0)
  } else {
    cli::cli_alert_success("已检测到 nFeature_RNA，跳过重算。")
  }

  #-------------------------------------------------------------------------------
  # 🔬 计算线粒体基因比例（^MT-）
  #-------------------------------------------------------------------------------
  cli::cli_alert_info("计算 percent_mito（线粒体基因比例）...")
  seu <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent_mito")

  #-------------------------------------------------------------------------------
  # 🔢 计算表达复杂度指标：log10(nFeature_RNA) / log10(nCount_RNA)
  #-------------------------------------------------------------------------------
  cli::cli_alert_info("计算表达复杂度指标 log10_ratio_features_to_umi ...")
  seu$log10_ratio_features_to_umi <- log10(seu$nFeature_RNA) / log10(seu$nCount_RNA)

  #-------------------------------------------------------------------------------
  # 🩸 计算红细胞基因比例（可选）
  #-------------------------------------------------------------------------------
  default_hb_genes <- c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "HBE1", "HBZ",
                        "HBM", "HBQ1", "HBA", "HBD1", "HBG", "HBE", "HBZP1", "HBDP1",
                        "HBG1P1", "HBG2P1", "HBA1P1", "HBA2P1")

  if (calculate_hb) {
    cli::cli_alert_info("计算 percent_hb（红细胞基因比例）...")
    hb_use <- if (is.null(hb_genes)) default_hb_genes else hb_genes
    matched_hb <- Seurat::CaseMatch(hb_use, rownames(seu))

    if (length(matched_hb) == 0) {
      cli::cli_alert_warning("未匹配到任何红细胞基因，percent_hb 将为 0。")
      seu@meta.data$percent_hb <- 0
    } else {
      seu <- Seurat::PercentageFeatureSet(seu, features = matched_hb, col.name = "percent_hb")
    }
  }
  
  #-------------------------------------------------------------------------------
  # ⚙️ 计算核糖体基因比例（可选）
  #-------------------------------------------------------------------------------
  default_ribo_genes <- c("RPL5", "RPS6", "RPL10", "RPS19", "RPL11", "RPS3", "RPL13A",
                          "RPS4X", "RPL7", "RPS27A", "RPL23A", "RPS18", "RPL32", "RPS15A",
                          "RPL37A", "RPS11", "RPL39", "RPS29", "RPLP1", "RPS24")

  if (calculate_ribo) {
    cli::cli_alert_info("计算 percent_ribo（核糖体基因比例）...")
    ribo_use <- if (is.null(ribo_genes)) default_ribo_genes else ribo_genes
    matched_ribo <- Seurat::CaseMatch(ribo_use, rownames(seu))

    if (length(matched_ribo) == 0) {
      cli::cli_alert_warning("未匹配到任何核糖体基因，percent_ribo 将为 0。")
      seu@meta.data$percent_ribo <- 0
    } else {
      seu <- Seurat::PercentageFeatureSet(seu, features = matched_ribo, col.name = "percent_ribo")
    }
  }
  
  #-------------------------------------------------------------------------------
  # 📌 打印随机基因名供检查
  #-------------------------------------------------------------------------------
  cli::cli_alert_info("示例基因名（随机 10 个）：{paste(sample(rownames(seu), 10), collapse = ', ')}")

  #-------------------------------------------------------------------------------
  # ✅ 返回计算完成的 Seurat 对象
  #-------------------------------------------------------------------------------
  cli::cli_alert_success("质控指标计算完成！")
  return(seu)
}

#-------------------------------------------------------------------------------