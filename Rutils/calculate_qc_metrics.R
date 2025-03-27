# Rutils/calculate_qc_metrics.R
#-------------------------------------------------------------------------------

# scFlowKit: Calculate QC Metrics for Single-Cell RNA-seq Data
#-------------------------------------------------------------------------------

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
#   sce: Seurat 对象，包含单细胞 RNA-seq 数据
#   calculate_hb: 是否计算红细胞基因比例，默认 FALSE
#   calculate_ribo: 是否计算核糖体基因比例，默认 FALSE
#   hb_genes: 红细胞基因列表，默认 NULL（使用内置默认基因集）
#   ribo_genes: 核糖体基因列表，默认 NULL（使用内置默认基因集）
calculate_qc_metrics <- function(sce, 
                                 calculate_hb = FALSE, 
                                 calculate_ribo = FALSE,
                                 hb_genes = NULL,
                                 ribo_genes = NULL) {
  # 验证输入参数是否为 Seurat 对象
  if (!inherits(sce, "Seurat")) {
    stop("参数 'sce' 必须为 Seurat 对象！", call. = FALSE)
  }
  
  # 验证 calculate_hb 和 calculate_ribo 是否为逻辑值
  if (!is.logical(calculate_hb)) {
    stop("参数 'calculate_hb' 必须为逻辑值！", call. = FALSE)
  }
  if (!is.logical(calculate_ribo)) {
    stop("参数 'calculate_ribo' 必须为逻辑值！", call. = FALSE)
  }
  
  # 验证 hb_genes 和 ribo_genes（如果非 NULL）
  if (!is.null(hb_genes) && (!is.character(hb_genes) || length(hb_genes) == 0)) {
    stop("参数 'hb_genes' 必须为非空字符向量！", call. = FALSE)
  }
  if (!is.null(ribo_genes) && (!is.character(ribo_genes) || length(ribo_genes) == 0)) {
    stop("参数 'ribo_genes' 必须为非空字符向量！", call. = FALSE)
  }
  
  # 提示用户正在计算质控指标
  message("正在计算质控指标...")
  
  # 检查是否需要计算 nCount_RNA 和 nFeature_RNA
  # - 如果 Seurat 对象已包含这些指标，则直接使用
  if (!"nCount_RNA" %in% colnames(sce@meta.data)) {
    message("计算 nCount_RNA（总 UMI 计数）...")
    sce@meta.data$nCount_RNA <- Matrix::colSums(sce[["RNA"]]$counts)
  } else {
    message("使用已有的 nCount_RNA（总 UMI 计数）...")
  }
  
  if (!"nFeature_RNA" %in% colnames(sce@meta.data)) {
    message("计算 nFeature_RNA（基因数）...")
    sce@meta.data$nFeature_RNA <- Matrix::colSums(sce[["RNA"]]$counts > 0)
  } else {
    message("使用已有的 nFeature_RNA（基因数）...")
  }
  
  # 计算线粒体基因比例
  # - 线粒体基因通常以 "MT-" 开头（人类数据）
  message("计算线粒体基因比例...")
  sce <- Seurat::PercentageFeatureSet(sce, pattern = "^MT-", col.name = "percent_mito")

  # 计算 log10_ratio_features_to_umi
  # - log10(nFeature_RNA) / log10(nCount_RNA)
  message("计算 log10_ratio_features_to_umi（log10(nFeature_RNA) / log10(nCount_RNA)）...")
  sce$log10_ratio_features_to_umi <- log10(sce$nFeature_RNA) / log10(sce$nCount_RNA)
  
  # 定义默认红细胞基因集
  # - 常见的红细胞基因
  default_hb_genes <- c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "HBE1", "HBZ", 
                        "HBM", "HBQ1", "HBA", "HBD1", "HBG", "HBE", "HBZP1", "HBDP1", 
                        "HBG1P1", "HBG2P1", "HBA1P1", "HBA2P1")
  
  # 定义默认核糖体基因集
  # - 常见的核糖体基因（示例，可扩展）
  default_ribo_genes <- c("RPL5", "RPS6", "RPL10", "RPS19", "RPL11", "RPS3", "RPL13A", 
                          "RPS4X", "RPL7", "RPS27A", "RPL23A", "RPS18", "RPL32", "RPS15A", 
                          "RPL37A", "RPS11", "RPL39", "RPS29", "RPLP1", "RPS24")
  
  # 获取 Seurat 对象的基因列表（原始大小写）
  sce_genes <- rownames(sce)
  
  # 打印随机取样的 10 个基因名，检查大小写
  message("Seurat 对象的基因名（随机取样 10 个，检查大小写）：", paste(sample(sce_genes, 10), collapse = ", "))
  
  # 计算红细胞基因比例（可选）
  # - 如果用户未提供基因列表，使用默认基因集
  if (calculate_hb) {
    message("计算红细胞基因比例...")
    hb_genes_to_use <- if (is.null(hb_genes)) default_hb_genes else hb_genes
    
    # 使用 CaseMatch 进行大小写不敏感匹配
    # - search: 目标基因集 (hb_genes_to_use)
    # - match: Seurat 对象的基因名 (sce_genes)
    matched_hb_genes <- Seurat::CaseMatch(search = hb_genes_to_use, match = sce_genes)
    
    # 检查基因是否存在
    missing_hb_genes <- setdiff(hb_genes_to_use, matched_hb_genes)
    if (length(missing_hb_genes) > 0) {
      warning("以下红细胞基因未在 Seurat 对象中找到，将被忽略：", paste(missing_hb_genes, collapse = ", "), call. = FALSE)
    }
    
    # 如果所有基因都不存在，发出警告
    if (length(matched_hb_genes) == 0) {
      warning("没有找到任何红细胞基因，percent_hb 将为 0！", call. = FALSE)
      sce@meta.data$percent_hb <- 0
    } else {
      # 使用匹配后的基因名（已转换为 sce_genes 中的格式）计算比例
      sce <- Seurat::PercentageFeatureSet(sce, features = matched_hb_genes, col.name = "percent_hb")
    }
  }
  
  # 计算核糖体基因比例（可选）
  # - 如果用户未提供基因列表，使用默认基因集
  if (calculate_ribo) {
    message("计算核糖体基因比例...")
    ribo_genes_to_use <- if (is.null(ribo_genes)) default_ribo_genes else ribo_genes
    
    # 使用 CaseMatch 进行大小写不敏感匹配
    # - search: 目标基因集 (ribo_genes_to_use)
    # - match: Seurat 对象的基因名 (sce_genes)
    matched_ribo_genes <- Seurat::CaseMatch(search = ribo_genes_to_use, match = sce_genes)
    
    # 检查基因是否存在
    missing_ribo_genes <- setdiff(ribo_genes_to_use, matched_ribo_genes)
    if (length(missing_ribo_genes) > 0) {
      warning("以下核糖体基因未在 Seurat 对象中找到，将被忽略：", paste(missing_ribo_genes, collapse = ", "), call. = FALSE)
    }
    
    # 如果所有基因都不存在，发出警告
    if (length(matched_ribo_genes) == 0) {
      warning("没有找到任何核糖体基因，percent_ribo 将为 0！", call. = FALSE)
      sce@meta.data$percent_ribo <- 0
    } else {
      # 使用匹配后的基因名（已转换为 sce_genes 中的格式）计算比例
      sce <- Seurat::PercentageFeatureSet(sce, features = matched_ribo_genes, col.name = "percent_ribo")
    }
  }
  
  # 提示用户计算完成
  message("质控指标计算完成！")
  
  return(sce)
}

#-------------------------------------------------------------------------------