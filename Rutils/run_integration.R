# Rutils/run_integration.R
#-------------------------------------------------------------------------------

# scFlowKit: Run Integration for Single-Cell RNA-seq Data
#-------------------------------------------------------------------------------

# Integration 背景介绍
# - 在单细胞 RNA-seq 数据中，批次效应（batch effect）是常见的技术噪声，可能掩盖生物学差异。
# - Seurat 提供多种整合方法，包括 CCA（Canonical Correlation Analysis）和 Harmony。
# - 整合方法：
#   - 不整合（method = "none"）：不进行显式批次校正（需要所有样本一起 SCTransform）。
#   - CCA（method = "cca"）：使用 anchors 整合（需要分组后的 Seurat 对象列表）。
#   - Harmony（method = "harmony"）：在 PCA 空间中校正批次效应。
# - 注意事项：
#   - 需要在 SCTransform 或 NormalizeData 后运行。
#   - 整合后的数据存储在 integrated assay（CCA）或 harmony 降维结果（Harmony）中。
#   - 计算量和内存需求因方法而异。

# run_integration: 运行 Seurat 整合
# 参数:
#   sce_list: Seurat 对象列表（可以是单个 Seurat 对象 list(sce) 或分组后的列表 list(donor1, donor2, ...)）
#   method: 整合方法（"none", "cca", "harmony"，默认 "cca"）
#   assay: 输入的 assay 名称（默认 "SCT"，也可为 "RNA")
#   k.anchor: 寻找 anchors 时的 k 参数（默认 5，仅在 method = "cca" 时有效）
#   k.filter: 过滤 anchors 时的 k 参数（默认 200，仅在 method = "cca" 时有效）
#   k.score: 评分 anchors 时的 k 参数（默认 30，仅在 method = "cca" 时有效）
#   new.assay.name: 整合后新 assay 的名称（默认 "integrated"，仅在 method = "cca" 时有效）
#   dims: 使用的维度（默认 1:30）
#   npcs: PCA 的主成分数量（默认 50，仅在 method = "harmony" 时有效）
#   harmony.group.by: Harmony 的分组变量（默认 "sample"，仅在 method = "harmony" 时有效）
#   harmony.theta: Harmony 的校正强度（默认 2，仅在 method = "harmony" 时有效）
#   harmony.sigma: Harmony 的软聚类分散度（默认 0.1，仅在 method = "harmony" 时有效）
#   variable.features.n: 选择的高变基因数量（默认 2000，仅在 method = "harmony" 或 method = "cca" 时有效）
#   verbose: 是否显示进度信息（默认 TRUE）
# 返回:
#   运行整合后的 Seurat 对象（根据 method 不同，包含 integrated assay 或 harmony 降维结果）
run_integration <- function(sce_list,
                            method = "cca",
                            assay = "SCT",
                            k.anchor = 5,
                            k.filter = 200,
                            k.score = 30,
                            new.assay.name = "integrated",
                            dims = 1:30,
                            npcs = 50,
                            harmony.group.by = "sample",
                            harmony.theta = 2,
                            harmony.sigma = 0.1,
                            variable.features.n = 2000,
                            verbose = TRUE) {
  # 验证输入参数：sce_list 必须提供
  if (is.null(sce_list)) {
    stop("必须提供 'sce_list' 参数！", call. = FALSE)
  }
  
  # 如果 sce_list 是单个 Seurat 对象，转换为列表
  if (inherits(sce_list, "Seurat")) {
    sce_list <- list(sce_list)
  }
  
  # 验证 sce_list 的每个元素是否为 Seurat 对象
  if (!all(sapply(sce_list, inherits, "Seurat"))) {
    stop("参数 'sce_list' 的每个元素必须为 Seurat 对象！", call. = FALSE)
  }
  
  # 检测输入类型并打印提示信息
  total_cells <- sum(sapply(sce_list, ncol))
  if (length(sce_list) == 1) {
    message("检测到输入为单个 Seurat 对象（包含 ", total_cells, " 个细胞）。")
  } else {
    message("检测到输入为分组的 Seurat 对象列表（包含 ", length(sce_list), " 个样本，", total_cells, " 个细胞）。")
  }
  
  # 验证 assay 是否存在
  if (!all(sapply(sce_list, function(x) assay %in% names(x@assays)))) {
    stop("参数 'assay' 必须存在于所有 Seurat 对象的 assay 中！", call. = FALSE)
  }
  
  # 验证 method 是否有效
  if (!method %in% c("none", "cca", "harmony")) {
    stop("参数 'method' 必须为 'none', 'cca' 或 'harmony'！", call. = FALSE)
  }
  
  # 验证 k.anchor 是否为正整数（仅在 method = "cca" 时）
  if (method == "cca" && (!is.numeric(k.anchor) || k.anchor <= 0)) {
    stop("参数 'k.anchor' 必须为正整数！", call. = FALSE)
  }
  
  # 验证 k.filter 是否为正整数（仅在 method = "cca" 时）
  if (method == "cca" && (!is.numeric(k.filter) || k.filter <= 0)) {
    stop("参数 'k.filter' 必须为正整数！", call. = FALSE)
  }
  
  # 验证 k.score 是否为正整数（仅在 method = "cca" 时）
  if (method == "cca" && (!is.numeric(k.score) || k.score <= 0)) {
    stop("参数 'k.score' 必须为正整数！", call. = FALSE)
  }
  
  # 验证 new.assay.name 是否为字符（仅在 method = "cca" 时）
  if (method == "cca" && (!is.character(new.assay.name) || length(new.assay.name) != 1)) {
    stop("参数 'new.assay.name' 必须为单一字符！", call. = FALSE)
  }
  
  # 验证 dims 是否为数值向量
  if (!is.numeric(dims) || length(dims) < 1) {
    stop("参数 'dims' 必须为数值向量，且长度至少为 1！", call. = FALSE)
  }
  
  # 验证 npcs 是否为正整数（仅在 method = "harmony" 时）
  if (method == "harmony" && (!is.numeric(npcs) || npcs <= 0)) {
    stop("参数 'npcs' 必须为正整数！", call. = FALSE)
  }
  
  # 验证 harmony.group.by 是否为字符（仅在 method = "harmony" 时）
  if (method == "harmony" && (!is.character(harmony.group.by) || length(harmony.group.by) != 1)) {
    stop("参数 'harmony.group.by' 必须为单一字符！", call. = FALSE)
  }
  
  # 验证 harmony.theta 是否为数值（仅在 method = "harmony" 时）
  if (method == "harmony" && !is.numeric(harmony.theta)) {
    stop("参数 'harmony.theta' 必须为数值！", call. = FALSE)
  }
  
  # 验证 harmony.sigma 是否为正数值（仅在 method = "harmony" 时）
  if (method == "harmony" && (!is.numeric(harmony.sigma) || harmony.sigma <= 0)) {
    stop("参数 'harmony.sigma' 必须为正数值！", call. = FALSE)
  }
  
  # 验证 variable.features.n 是否为正整数（仅在 method = "harmony" 或 method = "cca" 时）
  if ((method == "harmony" || method == "cca") && (!is.numeric(variable.features.n) || variable.features.n <= 0)) {
    stop("参数 'variable.features.n' 必须为正整数！", call. = FALSE)
  }
  
  # 验证 verbose 是否为逻辑值
  if (!is.logical(verbose)) {
    stop("参数 'verbose' 必须为逻辑值！", call. = FALSE)
  }
  
  # 加载 Seurat 包
  library(Seurat)
  
  # 根据 method 选择整合方式
  if (method == "none") {
    # 不整合：检查输入是否为单个 Seurat 对象
    if (length(sce_list) > 1) {
      stop("不整合（method = 'none'）需要所有样本一起运行 SCTransform，输入为分组的 Seurat 对象列表（包含 ", length(sce_list), " 个样本，", total_cells, " 个细胞）。请重新运行 SCTransform，不设置 split.by 参数。")
    }
    message("输入为单个 Seurat 对象（包含 ", total_cells, " 个细胞），直接返回...")
    sce_list <- sce_list[[1]]
    message("不整合完成！")
    message("Seurat 对象基本信息：")
    print(sce_list)
    return(sce_list)
  } else if (method == "cca") {
    # CCA 整合：需要分组的 Seurat 对象列表
    if (length(sce_list) == 1) {
      stop("CCA 整合（method = 'cca'）需要分组的 Seurat 对象列表，输入为单个 Seurat 对象（包含 ", total_cells, " 个细胞）。请重新运行 SCTransform，设置 split.by 参数（例如 split.by = 'sample'）。")
    }
    message("使用 CCA 方法进行整合...")
    # 选择高变基因（anchor features）
    message("选择高变基因（anchor features）...")
    anchor_features <- SelectIntegrationFeatures(object.list = sce_list,
                                                 nfeatures = variable.features.n,
                                                 # assay = rep(assay, length(sce_list)),
                                                 verbose = verbose)
    # 准备 SCTransform 数据
    message("准备 SCTransform 数据...")
    sce_list <- PrepSCTIntegration(object.list = sce_list,
                                   anchor.features = anchor_features,
                                   # assay = rep(assay, length(sce_list)),
                                   verbose = verbose)
    # 寻找整合 anchors
    message("寻找整合 anchors...")
    message("过程需要花费一段时间...")
    anchors <- FindIntegrationAnchors(object.list = sce_list,
                                      # assay = rep(assay, length(sce_list)),
                                      normalization.method = "SCT",
                                      reduction = method,
                                      anchor.features = anchor_features,
                                      k.anchor = k.anchor,
                                      k.filter = k.filter,
                                      k.score = k.score,
                                      dims = dims,
                                      verbose = verbose)
    # 运行整合
    message("运行整合...")
    sce_integrated <- IntegrateData(anchorset = anchors,
                                    new.assay.name = new.assay.name,
                                    normalization.method = "SCT",
                                    dims = dims,
                                    verbose = verbose)
    message("整合完成！")
    message("整合后的 Seurat 对象基本信息：")
    print(sce_integrated)
    return(sce_integrated)
  } else if (method == "harmony") {
    # Harmony 整合：需要单个 Seurat 对象
    message("使用 Harmony 方法进行整合...")
    if (length(sce_list) == 1) {
      message("输入为单个 Seurat 对象（包含 ", total_cells, " 个细胞），直接进行 PCA 和 Harmony 整合...")
      seurat_obj <- sce_list[[1]]
    } else {
      message("输入为分组的 Seurat 对象列表（包含 ", length(sce_list), " 个样本，", total_cells, " 个细胞），将合并所有样本...")
      # 选择高变基因
      message("选择高变基因...")
      var_features <- SelectIntegrationFeatures(object.list = sce_list,
                                                nfeatures = variable.features.n,
                                                assay = rep(assay, length(sce_list)),
                                                verbose = verbose)
      # 合并所有样本
      message("合并所有样本...")
      seurat_obj <- merge(sce_list[[1]], sce_list[-1], merge.data = TRUE)
      # 设置默认 assay 为 SCT
      DefaultAssay(seurat_obj) <- assay
      # 设置高变基因
      VariableFeatures(seurat_obj, assay = assay) <- var_features
    }
    # 运行 PCA
    message("运行 PCA...")
    seurat_obj <- RunPCA(seurat_obj,
                         assay = assay,
                         npcs = npcs,
                         verbose = verbose)
    # 运行 Harmony
    message("运行 Harmony 整合...")
    if (!requireNamespace("harmony", quietly = TRUE)) {
      stop("Harmony 整合需要安装 harmony 包，请安装：\n",
           "if (!requireNamespace('BiocManager', quietly = TRUE))\n",
           "    install.packages('BiocManager')\n",
           "BiocManager::install('harmony')")
    }
    library(harmony)
    sce_integrated <- RunHarmony(seurat_obj,
                                 group.by.vars = harmony.group.by,
                                 theta = harmony.theta,
                                 sigma = harmony.sigma,
                                 dims.use = dims,
                                 verbose = verbose)
    message("整合完成！")
    message("整合后的 Seurat 对象基本信息：")
    print(sce_integrated)
    return(sce_integrated)
  }
}

#-------------------------------------------------------------------------------