# Rutils/preprocess.R
#-------------------------------------------------------------------------------

# scFlowKit: Preprocess Single-Cell RNA-seq Data
#-------------------------------------------------------------------------------

# 自定义函数：预处理单细胞 RNA-seq 数据
#-------------------------------------------------------------------------------

# calculate_qc_metrics: 计算质控指标
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


#-------------------------------------------------------------------------------

# plot_qc_metrics: 可视化质控指标
# 参数:
#   sce: Seurat 对象，包含质控指标
#   output_dir: 输出目录，用于保存质控图
#   pt.size: VlnPlot 中点的显示大小，默认 0（不显示点）
plot_qc_metrics <- function(sce, output_dir, pt.size = 0) {
  # 验证输入参数是否为 Seurat 对象
  if (!inherits(sce, "Seurat")) {
    stop("参数 'sce' 必须为 Seurat 对象！", call. = FALSE)
  }

  # 验证 output_dir 是否为字符类型
  if (!is.character(output_dir)) {
    stop("参数 'output_dir' 必须为字符类型！", call. = FALSE)
  }

  # 验证 pt.size 是否为非负数值
  if (!is.numeric(pt.size) || pt.size < 0) {
    stop("参数 'pt.size' 必须为非负数值！", call. = FALSE)
  }

  # 验证元数据是否包含必要的质控指标
  required_metrics <- c("nCount_RNA", "nFeature_RNA", "percent_mito")
  missing_metrics <- setdiff(required_metrics, colnames(sce@meta.data))
  if (length(missing_metrics) > 0) {
    stop("元数据缺少必要的质控指标：", paste(missing_metrics, collapse = ", "), call. = FALSE)
  }

  # 提示用户正在可视化质控指标
  message("正在可视化质控指标...")

  # 确保输出目录存在
  figures_dir <- file.path(output_dir, "figures")
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  # 加载 Seurat、patchwork 和 ggplot2 包
  # - 已在 main.R 中加载，此处为模块独立性考虑
  # - 用户可移除此部分以避免重复加载
  library(Seurat)
  library(patchwork)
  library(ggplot2)

  # 使用 Seurat 的 VlnPlot 绘制质控指标分布
  # UMI 计数分布
  p1 <- VlnPlot(sce, features = "nCount_RNA", pt.size = pt.size, layer = "counts") +
        labs(title = "UMI Counts per Cell")
  
  # 基因数分布
  p2 <- VlnPlot(sce, features = "nFeature_RNA", pt.size = pt.size, layer = "counts") +
        labs(title = "Genes per Cell")
  
  # 线粒体基因比例分布
  p3 <- VlnPlot(sce, features = "percent_mito", pt.size = pt.size, layer = "counts") +
        labs(title = "Mitochondrial Gene %")
  
  # 使用 patchwork 将三张 VlnPlot 图拼成一张
  combined_vln_plot <- p1 + p2 + p3 + plot_layout(ncol = 3)
  
  # 保存 VlnPlot 拼图
  ggsave(file.path(figures_dir, "qc_metrics_combined.png"), combined_vln_plot, width = 15, height = 5)

  # 使用 Seurat 的 FeatureScatter 绘制质控指标散点图
  # UMI 计数 vs 基因数
  s1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
        labs(title = "UMI Counts vs Genes")
  
  # UMI 计数 vs 线粒体基因比例
  s2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent_mito") +
        labs(title = "UMI Counts vs Mito %")
  
  # 基因数 vs 线粒体基因比例
  s3 <- FeatureScatter(sce, feature1 = "nFeature_RNA", feature2 = "percent_mito") +
        labs(title = "Genes vs Mito %")
  
  # 使用 patchwork 将三张 FeatureScatter 图拼成一张
  combined_scatter_plot <- s1 + s2 + s3 + plot_layout(ncol = 3)
  
  # 保存 FeatureScatter 拼图
  ggsave(file.path(figures_dir, "qc_metrics_scatter_combined.png"), combined_scatter_plot, width = 15, height = 5)

  # 提示用户可视化完成
  message("质控指标可视化完成！")
}

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------

# filter_cells_genes: 过滤低质量细胞
# 参数:
#   sce: Seurat 对象，包含单细胞 RNA-seq 数据
#   min_umi: 最小 UMI 计数，默认 500
#   max_umi: 最大 UMI 计数，默认 15000（去除异常高值）
#   min_genes: 最小基因数，默认 200
#   max_genes: 最大基因数，默认 5000
#   max_mito: 最大线粒体基因比例（百分比），默认 10
#   filter_hb: 是否过滤红细胞基因比例高的细胞，默认 FALSE
#   max_hb: 最大红细胞基因比例（百分比），默认 5（仅当 filter_hb = TRUE 时生效）
#   filter_ribo: 是否过滤核糖体基因比例高的细胞，默认 FALSE
#   max_ribo: 最大核糖体基因比例（百分比），默认 50（仅当 filter_ribo = TRUE 时生效）
filter_cells_genes <- function(sce, 
                               min_umi = 500, 
                               max_umi = 15000, 
                               min_genes = 200, 
                               max_genes = 5000, 
                               max_mito = 10,
                               filter_hb = FALSE, 
                               max_hb = 5,
                               filter_ribo = FALSE, 
                               max_ribo = 50) {
  # 验证输入参数是否为 Seurat 对象
  if (!inherits(sce, "Seurat")) {
    stop("参数 'sce' 必须为 Seurat 对象！", call. = FALSE)
  }

  # 验证数值参数是否为数值类型且合理
  if (!is.numeric(min_umi) || min_umi < 0) {
    stop("参数 'min_umi' 必须为非负数值！", call. = FALSE)
  }
  if (!is.numeric(max_umi) || max_umi <= min_umi) {
    stop("参数 'max_umi' 必须为数值且大于 min_umi！", call. = FALSE)
  }
  if (!is.numeric(min_genes) || min_genes < 0) {
    stop("参数 'min_genes' 必须为非负数值！", call. = FALSE)
  }
  if (!is.numeric(max_genes) || max_genes <= min_genes) {
    stop("参数 'max_genes' 必须为数值且大于 min_genes！", call. = FALSE)
  }
  if (!is.numeric(max_mito) || max_mito < 0) {
    stop("参数 'max_mito' 必须为非负数值！", call. = FALSE)
  }

  # 验证 filter_hb 和 filter_ribo 是否为逻辑值
  if (!is.logical(filter_hb)) {
    stop("参数 'filter_hb' 必须为逻辑值！", call. = FALSE)
  }
  if (!is.logical(filter_ribo)) {
    stop("参数 'filter_ribo' 必须为逻辑值！", call. = FALSE)
  }

  # 验证 max_hb 和 max_ribo 是否为非负数值
  if (!is.numeric(max_hb) || max_hb < 0) {
    stop("参数 'max_hb' 必须为非负数值！", call. = FALSE)
  }
  if (!is.numeric(max_ribo) || max_ribo < 0) {
    stop("参数 'max_ribo' 必须为非负数值！", call. = FALSE)
  }

  # 验证元数据是否包含必要的质控指标
  required_metrics <- c("nCount_RNA", "nFeature_RNA", "percent_mito")
  missing_metrics <- setdiff(required_metrics, colnames(sce@meta.data))
  if (length(missing_metrics) > 0) {
    stop("元数据缺少必要的质控指标：", paste(missing_metrics, collapse = ", "), call. = FALSE)
  }

  # 如果 filter_hb = TRUE，验证元数据是否包含 percent_hb
  if (filter_hb && !"percent_hb" %in% colnames(sce@meta.data)) {
    stop("元数据缺少 'percent_hb' 指标，请在 calculate_qc_metrics 中设置 calculate_hb = TRUE！", call. = FALSE)
  }

  # 如果 filter_ribo = TRUE，验证元数据是否包含 percent_ribo
  if (filter_ribo && !"percent_ribo" %in% colnames(sce@meta.data)) {
    stop("元数据缺少 'percent_ribo' 指标，请在 calculate_qc_metrics 中设置 calculate_ribo = TRUE！", call. = FALSE)
  }

  # 提示用户正在过滤数据
  message("正在过滤低质量细胞和基因...")

  # 构建过滤条件，使用sprintf格式化字符串,%d是整数占位符，会被替换为后面的参数
  subset_conditions <- list(
    nCount_RNA = sprintf("nCount_RNA > %f & nCount_RNA < %f", min_umi, max_umi),
    nFeature_RNA = sprintf("nFeature_RNA > %f & nFeature_RNA < %f", min_genes, max_genes),
    percent_mito = sprintf("percent_mito < %f", max_mito)
  )

  # 如果 filter_hb = TRUE，添加红细胞基因比例过滤条件
  if (filter_hb) {
    subset_conditions$percent_hb <- sprintf("percent_hb < %f", max_hb)
  }

  # 如果 filter_ribo = TRUE，添加核糖体基因比例过滤条件
  if (filter_ribo) {
    subset_conditions$percent_ribo <- sprintf("percent_ribo < %f", max_ribo)
  }

  # 合并过滤条件
  subset_expr <- paste(unlist(subset_conditions), collapse = " & ")

  # 过滤低质量细胞
  sce <- subset(sce, subset = !!rlang::parse_expr(subset_expr))

  # 提示用户过滤完成
  message("低质量细胞和基因过滤完成！")

  return(sce)
}

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# remove_doublets: 使用 DoubletFinder 检测和去除双细胞
# 参数:
#   sce: Seurat 对象，包含单细胞 RNA-seq 数据（需已完成 PCA 降维）
#   PCs: 用于双细胞检测的主成分范围（默认 1:20）
#   doublet_rate: 假设的双细胞比例（默认 0.08，可根据实验设计调整）
#   pN: DoubletFinder 参数，pN 值（默认 0.25）
#   sct: 是否使用 SCTransform 数据（默认 FALSE）
# 返回:
#   过滤双细胞后的 Seurat 对象
#-------------------------------------------------------------------------------

# remove_doublets: 使用 DoubletFinder 检测和去除双细胞
# 参数:
#   sce: Seurat 对象，包含单细胞 RNA-seq 数据（需已完成 PCA 降维）
#   PCs: 用于双细胞检测的主成分范围（默认 1:20）
#   doublet_rate: 假设的双细胞比例（默认 0.08，可根据实验设计调整）
#   pN: DoubletFinder 参数，pN 值（默认 0.25）
#   sct: 是否使用 SCTransform 数据（默认 FALSE）
# 返回:
#   过滤双细胞后的 Seurat 对象
remove_doublets <- function(sce, PCs = 1:20, doublet_rate = 0.08, pN = 0.25, sct = FALSE) {
  # 验证输入参数是否为 Seurat 对象
  if (!inherits(sce, "Seurat")) {
    stop("参数 'sce' 必须为 Seurat 对象！", call. = FALSE)
  }

  # 验证是否已完成 PCA 降维
  if (!"pca" %in% names(sce@reductions)) {
    stop("Seurat 对象尚未完成 PCA 降维，请先运行 RunPCA！", call. = FALSE)
  }

  # 验证 PCs 参数是否为数值向量
  if (!is.numeric(PCs) || length(PCs) < 2) {
    stop("参数 'PCs' 必须为数值向量，且长度至少为 2！", call. = FALSE)
  }

  # 验证 PCs 是否在有效范围内
  max_pc <- ncol(sce@reductions$pca@cell.embeddings)
  if (max(PCs) > max_pc) {
    stop("PCs 参数超出 PCA 降维的主成分数量（最大为 ", max_pc, "）！", call. = FALSE)
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
  library(DoubletFinder)

  # 扫描参数，找到最佳 pK 值
  message("扫描参数，找到最佳 pK 值...")
  sweep.res <- paramSweep(sce, PCs = PCs, sct = sct, num.cores = 1)  # 单核运行，兼容 Windows
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)

  # 选择最佳 pK 值（通常是最大 BCMVN 对应的 pK）
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  message("最佳 pK 值：", pK)

  # 计算预期双细胞数量
  nExp_poi <- round(doublet_rate * ncol(sce) * ncol(sce) / 10000)
  message("预期双细胞数量：", nExp_poi)

  # 运行 DoubletFinder 预测双细胞
  message("运行 DoubletFinder 预测双细胞...")
  sce <- doubletFinder(sce, 
                       PCs = PCs, 
                       pN = pN, 
                       pK = pK, 
                       nExp = nExp_poi, 
                       sct = sct)

  # 关闭所有未使用的连接，避免警告
  closeAllConnections()

  # DoubletFinder 添加的列名通常为 "DF.classifications_pN_pK_nExp"
  df_col <- paste0("DF.classifications_", pN, "_", pK, "_", nExp_poi)

  # 过滤掉双细胞
  message("过滤双细胞...")
  sce <- subset(sce, subset = !!sym(df_col) == "Singlet")

  return(sce)
}

#-------------------------------------------------------------------------------


