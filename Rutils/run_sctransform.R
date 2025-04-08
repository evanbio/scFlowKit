# Rutils/run_sctransform.R
#-------------------------------------------------------------------------------

# scFlowKit: Run SCTransform for Single-Cell RNA-seq Data
#-------------------------------------------------------------------------------

# SCTransform 背景介绍
# - SCTransform 是 Seurat 提供的一种标准化方法，替代传统的 NormalizeData、FindVariableFeatures 和 ScaleData。
# - 核心原理：
#   - 使用正则化负二项模型（regularized negative binomial regression）对基因表达数据进行标准化。
#   - 自动处理测序深度（UMI 总数）和批次效应（通过 vars.to.regress 参数）。
#   - 同时生成标准化后的数据（data 层）、高变基因（variable features）和缩放后的数据（scale.data 层）。
# - 优点：
#   - 更适合多样本数据，能够减少批次效应。
#   - 简化流程，减少手动参数调整。
# - 注意事项：
#   - 计算量较大，可能需要更多内存和时间。
#   - 默认不回归细胞周期效应，可通过 vars.to.regress 参数指定（例如 c("S.Score", "G2M.Score")）。
#   - 输出为一个新的 SCT assay，原始 RNA assay 保持不变。
#   - 当 vst.flavor = "v2" 时，推荐安装 glmGamPoi 包以加速计算。

# run_sctransform: 运行 SCTransform 标准化
# 参数:
#   seu: Seurat 对象，包含单细胞 RNA-seq 数据
#   vars.to.regress: 回归掉的变量（默认 NULL，可指定为 c("S.Score", "G2M.Score") 或 "percent_mito")
#   variable.features.n: 选择的高变基因数量（默认 3000）
#   assay: 输入的 assay 名称（默认 "RNA")
#   split.by: 分组变量，用于多样本分别运行 SCTransform（默认 NULL，例如 "sample")
#   method: SCTransform 的方法（默认 "glmGamPoi"，可选 "poisson" 或 "nbinom"）
#   vst.flavor: SCTransform 的变体（默认 "v2"，可选 "v1"）
#   ncells: 用于计算的细胞数量（默认 NULL，动态设置）
#   seed.use: 随机种子（默认 42）
#   verbose: 是否显示进度信息（默认 TRUE）
# 返回:
#   如果 split.by 不为 NULL，返回分组后的 Seurat 对象列表；否则返回单个 Seurat 对象

run_sctransform <- function(seu,
                            vars.to.regress = NULL,
                            variable.features.n = 3000,
                            assay = "RNA",
                            split.by = NULL,
                            method = "glmGamPoi",
                            vst.flavor = "v2",
                            ncells = NULL,
                            seed.use = 42,
                            verbose = TRUE) {
  
  cli::cli_h2("SCTransform 标准化")

  # ---------------- 参数检查 ----------------

  # 验证输入参数是否为 Seurat 对象
  if (!inherits(seu, "Seurat")) {
    stop("参数 'seu' 必须为 Seurat 对象！", call. = FALSE)
  }

  # 验证 assay 是否存在
  if (!assay %in% names(seu@assays)) {
    stop("参数 'assay' 必须为 Seurat 对象中的一个 assay！", call. = FALSE)
  }

  # 验证 vars.to.regress 是否为 NULL 或元数据中的列
  if (!is.null(vars.to.regress) && !all(vars.to.regress %in% colnames(seu@meta.data))) {
    stop("参数 'vars.to.regress' 中的变量必须存在于元数据中！", call. = FALSE)
  }

  # 验证 split.by 是否为 NULL 或元数据中的列
  if (!is.null(split.by) && (!is.character(split.by) || length(split.by) != 1 || !split.by %in% colnames(seu@meta.data))) {
    stop("参数 'split.by' 必须为单一字符，且存在于元数据中！", call. = FALSE)
  }

  # 验证 variable.features.n 是否为正整数
  if (!is.numeric(variable.features.n) || variable.features.n <= 0) {
    stop("参数 'variable.features.n' 必须为正整数！", call. = FALSE)
  }

  # 验证 method 是否有效
  if (!method %in% c("glmGamPoi", "poisson", "nbinom")) {
    stop("参数 'method' 必须为 'glmGamPoi', 'poisson' 或 'nbinom'！", call. = FALSE)
  }

  # 验证 vst.flavor 是否有效
  if (!vst.flavor %in% c("v1", "v2")) {
    stop("参数 'vst.flavor' 必须为 'v1' 或 'v2'！", call. = FALSE)
  }

  # 验证 ncells 是否为 NULL 或正整数
  if (!is.null(ncells) && (!is.numeric(ncells) || ncells <= 0)) {
    stop("参数 'ncells' 必须为 NULL 或正整数！", call. = FALSE)
  }

  # 验证 seed.use 是否为数值
  if (!is.numeric(seed.use)) {
    stop("参数 'seed.use' 必须为数值！", call. = FALSE)
  }

  # 验证 verbose 是否为逻辑值
  if (!is.logical(verbose)) {
    stop("参数 'verbose' 必须为逻辑值！", call. = FALSE)
  }

  # ---------------- 环境准备 ----------------
  suppressPackageStartupMessages(library(Seurat))

  # 检查 glmGamPoi 包是否安装（当 vst.flavor = "v2" 且 method = "glmGamPoi" 时）
  if (vst.flavor == "v2" && method == "glmGamPoi") {
    if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
      warning("glmGamPoi 包未安装，SCTransform 将回退到较慢的估计方法！\n",
              "请安装 glmGamPoi 包以加速计算：\n",
              "if (!requireNamespace('BiocManager', quietly = TRUE))\n",
              "    install.packages('BiocManager')\n",
              "BiocManager::install('glmGamPoi')")
      method <- "poisson"  # 回退到 poisson 方法
    }
  }

  # ---------------- 执行 SCTransform ----------------
  # 如果 split.by 不为 NULL，则按分组运行 SCTransform
  if (!is.null(split.by)) {
    cli::cli_text("按 {split.by} 分组运行 SCTransform...", .envir = environment())
    # 分割 Seurat 对象
    seu_list <- SplitObject(seu, split.by = split.by)

    # 对每个分组运行 SCTransform
    seu_list <- lapply(seu_list, function(obj) {
      # 动态设置 ncells：如果细胞数量小于 5000，则使用所有细胞；否则使用 5000
      n_cells <- ncol(obj)
      dynamic_ncells <- if (is.null(ncells)) {
        if (n_cells < 5000) n_cells else 5000
      } else {
        ncells
      }
      cli::cli_text("当前样本细胞数：{n_cells}，设置 ncells = {dynamic_ncells}")
      SCTransform(obj,
                  assay = assay,
                  vars.to.regress = vars.to.regress,
                  variable.features.n = variable.features.n,
                  method = method,
                  vst.flavor = vst.flavor,
                  ncells = dynamic_ncells,
                  seed.use = seed.use,
                  verbose = verbose)
    })

    # 返回分组结果（不合并）
    cli::cli_alert_success("SCTransform 标准化完成，返回分组后的 Seurat 对象列表！")
    return(seu_list)
    
  } else {
    # 动态设置 ncells：如果细胞数量小于 5000，则使用所有细胞；否则使用 5000
    n_cells <- ncol(seu)
    dynamic_ncells <- if (is.null(ncells)) {
      if (n_cells < 5000) n_cells else 5000
    } else {
      ncells
    }
    cli::cli_text("总细胞数：{n_cells}，设置 ncells = {dynamic_ncells}")

    # 直接运行 SCTransform
    seu <- SCTransform(seu,
                       assay = assay,
                       vars.to.regress = vars.to.regress,
                       variable.features.n = variable.features.n,
                       method = method,
                       vst.flavor = vst.flavor,
                       ncells = dynamic_ncells,
                       seed.use = seed.use,
                       verbose = verbose)
    cli::cli_alert_success("SCTransform 标准化完成！")
    return(seu)
  }
}
#-------------------------------------------------------------------------------