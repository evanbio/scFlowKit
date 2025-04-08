#-------------------------------------------------------------------------------
# scFlowKit: A Modular Single-Cell RNA-seq Analysis Pipeline
#-------------------------------------------------------------------------------

# 项目信息
# Name       : scFlowKit
# Description: 一个模块化的单细胞 RNA-seq 分析流程，支持 R 与 Python 互通
# GitHub     : https://github.com/evanbio/scFlowKit
# License    : MIT License
# Version    : [TODO: 后续填写版本号]

# 作者信息
# Author     : Evan Zhou (Yibin Zhou)
# Email      : evanzhou.bio@gmail.com
# Website    : https://academic.evanzhou.org
# GitHub     : https://github.com/evanbio
# LinkedIn   : https://www.linkedin.com/in/yibin-zhou

# 当前版本功能：
# - ✅ 数据加载：支持多格式（10X .mtx/.h5, AnnData .h5ad, Loom .loom, SCE .rds）
# - ✅ 预处理模块：QC、过滤、标准化（LogNorm/VST）、高变基因识别、细胞周期评分
# - ✅ 降维与聚类：PCA、t-SNE、UMAP、Louvain/Leiden 聚类
# - ✅ 差异表达分析：Wilcoxon/MAST 等方法识别 marker genes
# - ✅ 细胞类型注释：支持 marker-based 自动注释 + 自定义 marker 集合

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 安装依赖包
#-------------------------------------------------------------------------------

# 请先运行 Rutils/install.R 安装必要的依赖包 (初次运行main.R时)
# 运行以下命令：
# source("Rutils/install.R")

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 载入依赖包
#-------------------------------------------------------------------------------

# 自动载入必要的依赖包
# 运行以下命令：
source("Rutils/load.R")

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 全局参数设置
#-------------------------------------------------------------------------------

# 通过 source 方式加载项目全局配置（路径设置、环境参数等）
source("Rutils/global_config.R")

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 主流程：单细胞 RNA-seq 分析
# 主分析代码从这里开始
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 步骤 1：数据导入
#-------------------------------------------------------------------------------

# 导入数据加载模块
source("Rutils/load_data.R")
source("Rutils/sce2seu.R")
source("scRNAutils/read_sce.R") 

# 模式控制（后续可参数化）

use_sce <- FALSE  # ✅ 使用 SCE 模式（.rds/.h5ad/.loom/.h5/.mtx）

if (use_sce) {
  #-----------------------------------------------
  # 加载数据（SCE 模式）
  #-----------------------------------------------
  cli_h1("🧬 步骤 1：加载 SingleCellExperiment 数据")

    # ✅ 支持多个 SCE 文件（路径向量）
  sce_input_paths <- c(
    "data/raw/sample1.rds",
    "data/raw/sample2.h5ad"
    # 可继续添加更多样本
  )

  # 自动提取样本名（去掉路径和扩展名）
  sample_names <- basename(tools::file_path_sans_ext(sce_input_paths))

  # 检查文件是否存在
  missing_paths <- sce_input_paths[!file.exists(sce_input_paths)]
  if (length(missing_paths) > 0) {
    cli::cli_alert_danger("❌ 以下 SCE 文件不存在：{paste(missing_paths, collapse = ', ')}")
    stop()
  }

  # 初始化结果列表
  seu_list <- list()

  # 遍历读取每个样本
  for (i in seq_along(sce_input_paths)) {
    path <- sce_input_paths[i]
    name <- sample_names[i]

    cli::cli_h2("📦 读取样本：{name}")

    sce <- if (grepl("\\.rds$", path)) {
      cli::cli_text("📄 使用 readRDS 读取 {name}")
      readRDS(path)
    } else {
      cli::cli_text("📄 使用 read_sce 读取 {name}")
      read_sce(path)
    }

    if (is.null(sce)) {
      cli::cli_alert_warning("⚠️ 读取失败：{name}，跳过该样本。")
      next
    }

    # 转换为 Seurat 对象
    seu_tmp <- sce2seu(sce, counts_assay = "counts", project = name)
    seu_tmp$sample <- name
    seu_list[[name]] <- seu_tmp
  }

  # 合并 Seurat 对象
  if (length(seu_list) == 0) {
    cli::cli_alert_danger("❌ 没有任何 SCE 文件成功加载，终止。")
    stop()
  } else if (length(seu_list) == 1) {
    seu <- seu_list[[1]]
  } else {
    cli::cli_alert_info("🧪 合并多个 Seurat 对象 ...")
    seu <- merge(seu_list[[1]], y = seu_list[-1], add.cell.id = names(seu_list))
    seu <- JoinLayers(seu)
  }

} else {
  #-----------------------------------------------
  # 默认路径：使用 10X 格式（原始数据）
  #-----------------------------------------------
  
  # 指定数据集名称（用户可根据实际数据集调整）
  # - 支持单个整合样本（例如 "5k_pbmc_combined"）
  # - 支持多个样本（例如 c("5k_pbmc_donor1", "5k_pbmc_donor2", "5k_pbmc_donor3", "5k_pbmc_donor4")）
  # - 或者使用 list.files 自动读取数据目录下的所有样本：dataset_name <- list.files(path = data_path, pattern = "5k_pbmc_donor[0-9]+", full.names = FALSE)

  # 加载数据
  cli_h1("🧬 步骤 1：加载单细胞 RNA-seq 数据（10X）")

  dataset_name <- c("5k_pbmc_donor1", "5k_pbmc_donor2", "5k_pbmc_donor3", "5k_pbmc_donor4")
  if (length(dataset_name) == 1) {
    # 单个整合样本，直接加载
    seu <- load_data(base_path = data_path, 
                    dataset_name = dataset_name,
                    min_cells = 10,         # 基因至少在 10 个细胞中表达
                    min_features = 40,      # 细胞至少表达 40 个基因
                    project = dataset_name,  # Seurat 对象项目名称
                    assay = "RNA")          # 测序类型
  } else {
    # 多个样本，分别加载并整合
    seu_list <- list()
    for (ds in dataset_name) {
      cli_alert_info("📂 加载样本：{ds} ...")
      seu_tmp <- load_data(base_path = data_path, 
                          dataset_name = ds,
                          min_cells = 10,
                          min_features = 40,
                          project = ds,  # 使用当前样本名作为 Seurat 对象的 project 名称
                          assay = "RNA")
      # 添加样本信息（sample information）
      seu_tmp$sample <- ds
      seu_list[[ds]] <- seu_tmp
    }
    
    # 整合多个样本
    cli_alert_info("🧪 整合多个样本中 ...")
    seu <- merge(seu_list[[1]], y = seu_list[-1], add.cell.id = names(seu_list))
    
    # 合并多个count层
    # - add.cell.id 生成了多个层（counts.5k_pbmc_donor1, counts.5k_pbmc_donor2 等）
    # - JoinLayers 将这些层合并为一个 counts 层
    cli_alert_info("🧼 合并 assays 下的层（JoinLayers）...")
    seu <- JoinLayers(seu)
  }
}

# 输出基本信息，确认加载成功
cli_alert_success("✅ Seurat 对象已成功创建，基本信息如下：")
print(seu)

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 1.5：可选探索 Seurat 对象结构（开发调试模式）
#-------------------------------------------------------------------------------

source("Rutils/explore_seurat.R")  # 加载结构探索函数

# 控制是否启用探索模式
explore_mode <- TRUE  # ✅ 可改为 FALSE 关闭结构探索

# 调用探索函数
explore_seurat(seu, explore_mode = explore_mode)

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 2：数据预处理和质控
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 步骤 2.1：计算质控指标
#-------------------------------------------------------------------------------

# 导入质控指标计算模块
source("Rutils/calculate_qc_metrics.R")

#-------------------------------------------------------------------------------
# 说明：
# - 自动计算以下质控指标：
#     - nCount_RNA（每个细胞的总 UMI 计数）
#     - nFeature_RNA（每个细胞表达的基因数）
#     - percent_mito（线粒体基因比例，基于 "^MT-"）
#     - log10_ratio_features_to_umi（表达复杂度）
# - 可选启用：
#     - percent_hb：红细胞基因比例（calculate_hb = TRUE）
#     - percent_ribo：核糖体基因比例（calculate_ribo = TRUE）
# - 支持自定义基因集（hb_genes, ribo_genes）
# - 会随机打印 10 个基因名，供大小写格式检查
#-------------------------------------------------------------------------------

cli::cli_h1("🧪 步骤 2.1：计算质控指标")

# 计算质控指标（默认不计算红细胞 / 核糖体比例）
seu <- calculate_qc_metrics(seu,
                            calculate_hb = FALSE,     # ❌ 不计算 percent_hb
                            calculate_ribo = FALSE)   # ❌ 不计算 percent_ribo

# 查看计算结果（Seurat 元数据前几行）
cli::cli_h2("📋 质控指标计算结果（前几行）")
print(head(seu@meta.data))

#-------------------------------------------------------------------------------
# 步骤 2.2：可视化质控指标
#-------------------------------------------------------------------------------

# 导入质控指标可视化模块
source("Rutils/plot_qc_metrics.R")

#-------------------------------------------------------------------------------
# 说明：
# - 可视化 QC 指标分布与相关性：
#     - 分布图：小提琴图（Violin）
#     - 相关性：散点图（Scatter）+ 综合视图（Comprehensive）
# - 可视指标：
#     - nCount_RNA、nFeature_RNA、percent_mito、log10_ratio_features_to_umi
# - 每行显示一个指标（共 4 行），每列为不同样本（按 sample 字段分组）
# - 阈值线设置：
#     - nCount_RNA ≥ 500
#     - nFeature_RNA ≥ 300
#     - percent_mito ≤ 10
#     - ratio ≥ 0.8
# - 点大小设为 0.1，展示单细胞分布
# - 输出路径：
#     - 小提琴图：results/figures/qc_metrics_combined.png
#     - 散点图：results/figures/qc_metrics_scatter_combined.png
#     - 综合散点图：results/figures/qc_metrics_comprehensive.png
#-------------------------------------------------------------------------------

cli::cli_h1("📊 步骤 2.2：可视化质控指标")

plot_qc_metrics(seu,
                output_dir = output_dir,
                pt.size = 0.1,                # 设置点大小为 0.1
                umi_threshold = 500,
                feature_threshold = 300,
                mito_threshold = 10,
                ratio_threshold = 0.8)

# 打印输出文件路径
cli::cli_alert_success("✅ 质控图表已保存至：")
cli::cli_text("📌 小提琴图：{file.path(output_dir, 'figures', 'qc_metrics_combined.png')}")
cli::cli_text("📌 散点图：{file.path(output_dir, 'figures', 'qc_metrics_scatter_combined.png')}")
cli::cli_text("📌 综合图：{file.path(output_dir, 'figures', 'qc_metrics_comprehensive.png')}")

#-------------------------------------------------------------------------------
# 步骤 2.3：过滤低质量细胞
#-------------------------------------------------------------------------------

# 导入过滤低质量细胞模块
source("Rutils/filter_cells.R")

#-------------------------------------------------------------------------------
# 说明：
# - 基于质控指标过滤低质量细胞
# - 过滤条件：
#   - nCount_RNA > 500（最小 UMI 计数）
#   - nFeature_RNA > 300（最小基因数）
#   - percent_mito < 10（最大线粒体基因比例）
#   - log10_ratio_features_to_umi > 0.8（最小 log10 比值）
# - 不设置最大 UMI 计数和最大基因数（默认 Inf）
# - 不过滤红细胞和核糖体基因比例（默认 FALSE）

cli::cli_h2("🧹 步骤 2.3：过滤低质量细胞")

# 应用过滤函数
seu <- filter_cells(
  seu,
  min_umi    = 500,
  max_umi    = Inf,
  min_genes  = 300,
  max_genes  = Inf,
  max_mito   = 10,
  min_ratio  = 0.8,
  filter_hb  = FALSE,
  filter_ribo = FALSE
)

# 查看过滤后 Seurat 对象信息
cli::cli_text("📦 过滤后对象概览：")
print(seu)

# 保存过滤后的数据（中间点）
cli::cli_text("💾 保存过滤后的 Seurat 对象至文件：")
saveRDS(seu, file = file.path(processed_data_dir, "scFlowKit_filtered.rds"))
cli::cli_alert_success("✅ 已保存：{file.path(processed_data_dir, 'scFlowKit_filtered.rds')}")

#-------------------------------------------------------------------------------
# 步骤 2.4：标准化和对数化数据
#-------------------------------------------------------------------------------

# 可选：从 .rds 文件加载预处理后的 Seurat 对象（跳过步骤 2.1 到 2.3）
# - 加载路径：processed_data_dir/scFlowKit_filtered.rds
# - 确保 processed_data_dir 已定义
# sce <- readRDS(file.path(processed_data_dir, "scFlowKit_filtered.rds"))
# 
# 📊 标准化（Normalization）+ 对数化（Log-Transformation）
# - 标准化：将每个细胞的总表达量缩放到 10,000
# - 对数化：log1p 变换，稳定数据分布
# - 使用 Seurat 的 NormalizeData 函数，参数：
#   - normalization.method = "LogNormalize" 指定标准化方法，可选 "LogNormalize" 或 "CLR"
#   - scale_factor = 10000 指定缩放因子

cli::cli_h2("🧮 步骤 2.4：标准化与对数化数据")

seu <- Seurat::NormalizeData(
  object = seu,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

cli::cli_alert_success("✅ 数据标准化与对数化完成！")

# 查看标准化后的 Seurat 对象基本信息
cli::cli_text("📦 Seurat 对象结构（标准化后）:")
print(seu) # 标准化后的 Seurat 对象会增加一个 data 层

# 确认数据层结构，检查是否生成 data 层
cli::cli_text("📚 当前数据层（slots/layers）：")
print(SeuratObject::Layers(seu))

#-------------------------------------------------------------------------------
# 步骤 2.5：寻找可变基因
#-------------------------------------------------------------------------------

# - 识别高变异基因（highly variable genes），用于后续降维和聚类
# - 使用 Seurat 的 FindVariableFeatures 函数，参数：
#   - selection.method = "vst" 使用 variance stabilizing transformation 方法
#   - nfeatures = 2000 选择 2000 个高变异基因
#   - verbose = TRUE 显示进度条
# - 可变基因存储在 sce@assays$RNA@var.features 中
# - 基因的均值和方差存储在 sce@assays$RNA@meta.data 中（Seurat 5.0 及以上版本）

cli::cli_h2("🚀 步骤 2.5：寻找可变基因")

seu <- Seurat::FindVariableFeatures(
  object = seu,
  selection.method = "vst",  # 使用 VST 方法
  nfeatures = 2000,          # 选择 2000 个高变基因
  verbose = TRUE             # 显示进度条
)

# 提取全部可变基因名称（使用 VariableFeatures）
all_variable_genes <- Seurat::VariableFeatures(seu)
cli::cli_alert_info("共检测到 {length(all_variable_genes)} 个可变基因")

# 提取可变基因名称，打印变化最大的 10 个基因
top_variable_genes <- head(all_variable_genes, 10)
cli::cli_text("Top 10 可变基因：{paste(top_variable_genes, collapse = ', ')}")

# 📊 可视化：均值-方差散点图 + top 10 基因标注
cli::cli_text("🎨 绘制可变基因散点图...")

variable_feature_plot <- Seurat::VariableFeaturePlot(
  object = seu,
  log = NULL,                      # 默认根据数据决定
  col = c("black", "red"),         # 黑色 = 普通基因，红色 = 高变基因
  pt.size = 1                      # 点大小
)

variable_feature_plot <- Seurat::LabelPoints(
  plot = variable_feature_plot,
  points = top_variable_genes,
  repel = TRUE,                   # 避免标签重叠
  xnudge = 0.3,
  ynudge = 0.05
)

# 保存散点图
ggsave(
  filename = file.path(output_dir, "figures", "variable_features_plot.png"),
  plot = variable_feature_plot,
  width = 8,
  height = 6
)

cli::cli_alert_success("✅ 可变基因散点图保存成功：figures/variable_features_plot.png")



#-------------------------------------------------------------------------------
# 步骤 2.6：细胞周期评分
#-------------------------------------------------------------------------------

# - 使用 Seurat 的 CellCycleScoring 函数为细胞分配周期阶段
# - 使用 Seurat 内置的 cc.genes 数据集（包含 S 期和 G2/M 期基因）
# - cc.genes 包含 43 个 S 期基因和 54 个 G2/M 期基因
# - 使用标准化后的 data 层数据
# - 输出：添加 S.Score、G2M.Score 和 Phase 列到 sce@meta.data
# - 细胞周期阶段推断规则：
#   - S.Score > 0 且高于 G2M.Score：S 期
#   - G2M.Score > 0 且高于 S.Score：G2/M 期
#   - 两者均低或接近：G1 期
# - 放在 ScaleData 之前，以便 ScaleData 可以回归掉细胞周期影响

cli::cli_h2("🔬 步骤 2.6：细胞周期评分")

# 打印基因集数量
cli::cli_text("S 期基因数：{length(Seurat::cc.genes$s.genes)}")
cli::cli_text("G2/M 期基因数：{length(Seurat::cc.genes$g2m.genes)}")

# 进行评分
seu <- Seurat::CellCycleScoring(
  object = seu,
  s.features = Seurat::cc.genes$s.genes,
  g2m.features = Seurat::cc.genes$g2m.genes,
  set.ident = FALSE
)

# 查看评分结果（前 6 个细胞）
cli::cli_alert_info("细胞周期评分结果（前 6 个细胞）：")
print(head(seu@meta.data[, c("S.Score", "G2M.Score", "Phase")]))

cli::cli_alert_success("✅ 细胞周期评分完成！")

#-------------------------------------------------------------------------------
# 步骤 2.7：数据缩放（为双细胞检测准备）
#-------------------------------------------------------------------------------

# - 对基因进行中心化和缩放（零均值、单位方差）
# - 使用 Seurat 的 ScaleData 函数，常用参数：
#   - features 要缩放的基因，默认 NULL（使用高变异基因 VariableFeatures(sce)）
#   - vars.to.regress 回归掉的变量（比如 "percent_mito" 去除线粒体比例影响）
#   - scale.max 缩放后表达量的最大值（默认 10）
#   - do.scale 是否进行缩放（默认 TRUE）
#   - do.center 是否进行中心化（默认 TRUE）
#   - verbose = TRUE 显示进度信息
# - 这里选择使用高变异基因（features = VariableFeatures(sce)），以减少计算量
# - 结果存储在 sce[["RNA"]]$scale.data 中

cli::cli_h2("📐 步骤 2.7：数据缩放（为双细胞检测准备）")

seu <- Seurat::ScaleData(
  object = seu,
  features = VariableFeatures(sce),  # 使用高变异基因
  vars.to.regress = NULL,  # 不回归任何变量（可设置为 c("S.Score", "G2M.Score")）
  scale.max = 10,  # 缩放后表达量最大值
  do.scale = TRUE,  # 进行缩放
  do.center = TRUE,  # 进行中心化
  verbose = TRUE)  # 显示进度信息  


# 打印 Seurat 对象基本信息
cli::cli_text("✅ 数据缩放完成，Seurat 对象信息如下：")
print(seu)

# 展示缩放数据示例
cli::cli_text("缩放后表达矩阵（前 5 个基因 × 前 5 个细胞）：")
print(seu[["RNA"]]@scale.data[1:5, 1:5])

#-------------------------------------------------------------------------------
# 步骤 2.8：PCA 降维及可视化（为双细胞检测准备）
#-------------------------------------------------------------------------------

# - 使用 PCA 进行降维，基于高变异基因
# - 使用 Seurat 的 RunPCA 函数，常用参数：
#   - features 使用高变异基因（默认 VariableFeatures(sce)）
#   - npcs = 50 选择前 50 个主成分（可调整）
#   - verbose = TRUE 显示进度信息
# - 结果存储在 sce@reductions$pca 中

cli::cli_h2("📉 步骤 2.8：PCA 降维及可视化")

# ----------------- PCA 降维 -----------------
cli::cli_text("🚀 运行 PCA（基于高变异基因）...")
seu <- Seurat::RunPCA(
  object = seu,
  features = Seurat::VariableFeatures(seu), # 默认使用高变异基因
  npcs = 50,
  verbose = TRUE
)

# 输出降维后的 Seurat 对象信息
cli::cli_text("✅ PCA 结果已添加至 reductions$pca 中")
print(seu)


# 可视化 PCA 结果
source("Rutils/plot_sc_pca.R")

# ----------------- 可视化 PCA（细胞周期相关） -----------------
# - 第一行：按 sample 和 Phase 分组
# - 第二行：按 Phase 分组，按 Phase 分面
# - 保存为 output_dir/figures/preliminary_phase_pca_dimplot.png
cli::cli_text("🎨 可视化 PCA：按细胞周期分组...")
plot_sc_pca(
  seu,
  output_dir = output_dir,
  reduction = "pca",
  group.by = "Phase",
  split.by = "Phase",
  prefix = "preliminary_phase",       # <--- 设置前缀避免冲突
  plot_elbow = TRUE,      # <--- 是否绘制 ElbowPlot
  plot_heatmap = TRUE,    # <--- 是否绘制 Heatmap
  width = 10,
  height = 10
)

# 对连续变量进行分段
# ----------------- 分段处理 percent_mito -----------------
cli::cli_text("🔢 分段处理线粒体比例 percent_mito...")
if (is.numeric(head(seu@meta.data$percent_mito, 3))) {
  quartiles <- quantile(seu@meta.data$percent_mito, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  seu@meta.data$percent_mito_binned <- cut(
    seu@meta.data$percent_mito,
    breaks = quartiles,
    labels = c("Q1", "Q2", "Q3", "Q4"),
    include.lowest = TRUE
  )
}

# ----------------- 可视化 PCA（线粒体比例分段） -----------------
# - 第一行：按 sample 和 percent_mito 分组
# - 第二行：按 percent_mito 分组，按 percent_mito 分面
# - 保存为 output_dir/figures/preliminary_mito_binned_pca_dimplot.png
cli::cli_text("🎨 可视化 PCA：按线粒体比例分段分组...")
plot_sc_pca(
  seu,
  output_dir = output_dir,
  reduction = "pca",
  group.by = "percent_mito_binned",
  split.by = "percent_mito_binned",
  prefix = "preliminary_mito_binned",     # <--- 避免覆盖
  plot_elbow = FALSE,
  plot_heatmap = FALSE,
  width = 10,
  height = 10
)

# ----------------- 输出路径提示 -----------------
cli::cli_alert_success("🎯 PCA 图表已保存：")
cli::cli_text("- 细胞周期 PCA 图：{file.path(output_dir, 'figures', 'preliminary_phase_pca_dimplot.png')}")
cli::cli_text("- 线粒体 PCA 图：{file.path(output_dir, 'figures', 'preliminary_mito_binned_pca_dimplot.png')}")
cli::cli_text("📁 图表目录：{file.path(output_dir, 'figures')}")
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 步骤 2.9：双细胞检测和去除
#-------------------------------------------------------------------------------

# 导入双细胞检测模块
source("Rutils/remove_doublets.R")

# - 使用 DoubletFinder 检测和去除双细胞
# - DoubletFinder 是一种基于 PCA 降维的双细胞检测工具，通过模拟人工双细胞并比较邻居关系来识别双细胞。
# - 参数说明：
#   - PCs = 1:20：使用前 20 个主成分（之前 PCA 已计算 50 个主成分）
#   - doublet_rate = 0.08：假设双细胞比例为 8%（根据实验设计调整，10x Genomics 数据通常为 0.05-0.08）
#   - pN = 0.25：DoubletFinder 默认 pN 值（人工双细胞比例）
#   - sct = FALSE：不使用 SCTransform 数据（当前使用 RNA assay）
# - 封装在 remove_doublets 函数中（位于 Rutils/remove_doublets.R）
# - 预期输出：
#   - DoubletFinder 会在元数据中添加双细胞标签（DF.classifications）
#   - 过滤后的 Seurat 对象仅保留 Singlet 细胞

# 记录过滤前的细胞数量
pre_cell_count <- ncol(seu)
cli::cli_text("过滤前细胞数量：{pre_cell_count}")

cli::cli_h2("步骤 2.9：检测并去除双细胞")
sce <- remove_doublets(sce,
                       PCs = 1:20,  # 使用前 20 个主成分
                       doublet_rate = 0.08,  # 假设双细胞比例为 8%
                       pN = 0.25,  # DoubletFinder 参数 pN
                       sct = FALSE)  # 不使用 SCTransform 数据

# 记录过滤后的细胞数量
post_cell_count <- ncol(sce)
removed_count <- pre_cell_count - post_cell_count

# 输出过滤双细胞后的 Seurat 对象信息
cli::cli_text("过滤双细胞后的细胞数量：{post_cell_count}")
cli::cli_text("去除的细胞数量：{removed_count}")

cli::cli_text("过滤双细胞后的 Seurat 对象基本信息：")
print(seu)

# 保存去除双细胞后的数据（中间点）
filtered_rds_path <- file.path(processed_data_dir, "scFlowKit_doublet_removed.rds")
saveRDS(seu, file = filtered_rds_path)
cli::cli_alert_success("已保存 Seurat 对象（去除双细胞）：{filtered_rds_path}")

#-------------------------------------------------------------------------------
# 步骤 2.10：SCTransform 标准化
# 这个步骤需要消耗较多时间和内存
#-------------------------------------------------------------------------------

# 导入 SCTransform 模块
source("Rutils/run_sctransform.R")

# 可选：从 .rds 文件加载去除双细胞后的 Seurat 对象（跳过步骤 2.1 到 2.9）
# - 加载路径：processed_data_dir/scFlowKit_doublet_removed.rds
# - 确保 processed_data_dir 已定义
# seu <- readRDS(file.path(processed_data_dir, "scFlowKit_doublet_removed.rds"))

# - 使用 SCTransform 进行标准化，作用涵盖了 NormalizeData、FindVariableFeatures 和 ScaleData
# - 参数说明：
#   - vars.to.regress = NULL：不回归任何变量（可设置为 c("S.Score", "G2M.Score") 或 "percent_mito")
#   - variable.features.n = 3000：选择 3000 个高变基因
#   - assay = "RNA"：使用 RNA assay 作为输入
#   - split.by = "sample"：按 sample 分组运行 SCTransform（多样本数据）
#   - method = "glmGamPoi"：使用 glmGamPoi 方法（推荐）
#   - vst.flavor = "v2"：使用 SCTransform v2 变体
#   - ncells = NULL：动态设置（小于 5000 使用所有细胞，否则 5000）
#   - seed.use = 42：设置随机种子，确保可重复性
# - 预期输出：
#   - 返回分组后的 Seurat 对象列表（donor1 到 donor4），每个对象包含 SCT assay
#   - 后续整合步骤将处理分组对象列表
#   - 原始 RNA assay 保持不变

# ----------------- 执行 SCTransform -----------------
cli::cli_h2("步骤 2.10：SCTransform 标准化")
seu_list <- run_sctransform(seu,
                            vars.to.regress = NULL,
                            variable.features.n = 3000,
                            assay = "RNA",
                            split.by = "sample",  # 按 sample 分组运行
                            method = "glmGamPoi",
                            vst.flavor = "v2",
                            ncells = NULL, # 动态设置
                            seed.use = 42,
                            verbose = TRUE)

# 打印分组结果信息
cli::cli_text("标准化完成，输出分组的 Seurat 对象列表：")
print(seu_list)

# 保存SCTransform 标准化后的数据（中间点）
save_path <- file.path(processed_data_dir, "scFlowKit_run_sctransform.rds")
saveRDS(seu_list, file = save_path)
cli::cli_alert_success("✅ 保存SCTransform 标准化后的 Seurat 对象：{save_path}")

#-------------------------------------------------------------------------------
# 步骤 2.11：Integration（整合多样本）
#-------------------------------------------------------------------------------

# 导入 Integration 模块
source("Rutils/run_integration.R")

# 可选：从 .rds 文件加载 SCTransform 标准化后的 Seurat 对象列表（跳过步骤 2.10）
# - 加载路径：processed_data_dir/scFlowKit_run_sctransform.rds
# - 确保 processed_data_dir 已定义
# seu_list <- readRDS(file.path(processed_data_dir, "scFlowKit_run_sctransform.rds"))

# - 使用 IntegrateData 或 Harmony 整合多样本数据，去除批次效应
# - 参数说明：
#   - seu_list：分组后的 Seurat 对象列表（donor1 到 donor4）
#   - method = "cca"：使用 CCA 方法（Seurat 默认）
#   - assay = "SCT"：使用 SCT assay 作为输入
#   - k.anchor = 5：寻找 anchors 时的 k 参数（仅 CCA）
#   - k.filter = 200：过滤 anchors 时的 k 参数（仅 CCA）
#   - k.score = 30：评分 anchors 时的 k 参数（仅 CCA）
#   - new.assay.name = "integrated"：整合后新 assay 的名称（仅 CCA）
#   - dims = 1:30：使用的维度
#   - npcs = 50：PCA 的主成分数量（仅 Harmony）
#   - variable.features.n = 2000：选择的高变基因数量（仅 Harmony）
# - 预期输出：
#   - method = "cca"：创建 integrated assay，包含整合后的数据
#   - method = "harmony"：创建 harmony 降维结果，包含校正后的 PCA 嵌入
#   - method = "none"：合并后的 Seurat 对象，包含 SCT assay

cli::cli_h2("步骤 2.11：整合多样本数据（Integration）")

# 执行整合
cli::cli_text("运行整合函数 run_integration()，方法为 'cca'...")
seu_integrated <- run_integration(seu_list,
                                  method = "cca",  # 使用 CCA 方法
                                  assay = "SCT",
                                  k.anchor = 5,
                                  k.filter = 200,
                                  k.score = 30,
                                  new.assay.name = "integrated",
                                  dims = 1:30,
                                  npcs = 50,
                                  variable.features.n = 2000,
                                  verbose = TRUE)

# 保存整合后的 Seurat 对象（中间点）
cli::cli_text("保存整合后的 Seurat 对象...")
saveRDS(seu, file = file.path(processed_data_dir, "scFlowKit_integrated.rds"))
cli::cli_alert_success("已保存至：{file.path(processed_data_dir, 'scFlowKit_integrated.rds')}")


#-------------------------------------------------------------------------------
# 步骤 2.12：PCA（再次运行，基于整合后的数据）
#-------------------------------------------------------------------------------

# - 在整合后的数据（或未整合的 SCT assay）上运行 PCA 降维，用于后续聚类和可视化
# - 参数说明：
#   - sce_integrated：整合后的 Seurat 对象（包含 integrated assay 或 harmony 降维结果，或未整合的 SCT assay）
#   - method：整合方法（"none", "cca", "harmony"，默认 "cca"）
#   - assay：输入的 assay 名称（默认 "integrated"，仅在 method = "cca" 时有效；method = "none" 时使用 "SCT"）
#   - reduction：输入的降维结果（默认 "harmony"，仅在 method = "harmony" 时有效）
#   - npcs = 50：PCA 的主成分数量
#   - seed.use = 42：设置随机种子，确保可重复性
# - 预期输出：
#   - 包含 PCA 降维结果的 Seurat 对象（reductions 槽中）

cli::cli_h2("步骤 2.12：PCA（基于整合数据）")

# 根据整合方法选择 PCA 输入
method <- "cca"  # 根据实际使用的整合方法设置（"none", "cca", "harmony"）

if (method == "none") {
  # 未整合：使用 SCT assay 运行 PCA
  cli::cli_text("未整合，使用 SCT assay 运行 PCA...")
  seu <- RunPCA(seu,
                assay = "SCT",
                npcs = 50,
                seed.use = 42,
                verbose = TRUE)

} else if (method == "cca") {
  # 使用 integrated assay 运行 PCA
  cli::cli_text("未整合，使用 SCT assay 运行 PCA...")
  seu <- RunPCA(seu,
                assay = "SCT",
                npcs = 50,
                seed.use = 42,
                verbose = TRUE)
                
} else if (method == "harmony") {
  # 使用 harmony 降维结果运行 PCA
  cli::cli_text("使用 harmony 降维结果运行 PCA...")
  seu <- RunPCA(seu,
                reduction = "harmony",
                npcs = 50,
                seed.use = 42,
                verbose = TRUE)
}

# 打印 PCA 结果信息
cli::cli_alert_success("PCA 降维完成！")
cli::cli_text("Seurat 对象降维信息：")
print(seu)

# 保存 PCA 降维后的 Seurat 对象（中间点）
save_path <- file.path(processed_data_dir, "scFlowKit_pca.rds")
cli::cli_text("保存降维结果到：{save_path}")
saveRDS(seu, file = save_path)
cli::cli_alert_success("保存完成！")


#-------------------------------------------------------------------------------
# 步骤 2.13：可视化和探索 PCA 结果
#-------------------------------------------------------------------------------
# - 使用 visualize_pca 函数可视化 PCA 降维结果，检查细胞周期相关分布
# - 参数说明：
#   - sce_integrated：整合后的 Seurat 对象（包含 PCA 降维结果）
#   - output_dir：保存图形的目录
#   - reduction = "pca"：使用的降维结果
#   - dims = c(1, 2)：DimPlot 使用的 PCA 维度（PC1 和 PC2）
#   - group.by = "Phase"：按细胞周期阶段分组（G1, S, G2M）
#   - split.by = "Phase"：按细胞周期阶段分面
#   - ndims = 50：ElbowPlot 显示的主成分数量
#   - width = 10：图形宽度
#   - height = 10：图形高度
# - 预期输出：
#   - DimPlot：PCA 的二维散点图，按 Phase 分组和分面
#   - ElbowPlot：主成分的方差贡献图
#   - Heatmap：主成分的热图


# - 使用 visualize_pca 函数可视化 PCA 降维结果，检查细胞周期相关分布
message("步骤 2.13：可视化和探索 PCA 结果...")

# 导入 visualize_pca 模块
source("Rutils/plot_sc_pca.R")

# 运行 plot_sc_pca
plot_sc_pca(seu_integrated,
            output_dir = output_dir,
            reduction = "pca",
            dims = c(1, 2),  # 使用 PC1 和 PC2
            group.by = "Phase",  # 按细胞周期阶段分组
            split.by = "Phase",  # 按细胞周期阶段分面
            ndims = 50,  # 显示前 50 个主成分
            width = 10,
            height = 10,
            dpi = 300)

#-------------------------------------------------------------------------------
# 探索 PCA 结果：选择主成分（PCs）数量
#-------------------------------------------------------------------------------

# - 本步骤将根据 PCA 结果，使用 suggest_pcs 函数自动推荐合适的主成分数量
# - 该函数结合两种指标判断：
#     - 指标 1：累计方差贡献 > 90%，且当前 PC 贡献 < 5%
#     - 指标 2：主成分解释度下降显著（> 0.1%）的最后一个主成分（肘部位置）
# - 推荐值为以上两个指标中较小者，更稳健、避免噪声干扰

# 导入辅助函数
source("Rutils/suggest_pcs.R")

# 运行推荐函数
pcs_to_use <- suggest_pcs(seu_integrated, reduction = "pca", verbose = TRUE)

# 可选：打印前 10 个 PC 的 top 5 驱动基因（用于手动审阅）
message("前 10 个主成分的 Top 5 驱动基因如下：")
print(seu_integrated[["pca"]], dims = 1:10, nfeatures = 5)

# # 最终选择：基于 SCTransform 的经验值，使用前 40 个 PCs
# # - SCTransform 更准确，40 个 PCs 是合理的折中（保留足够变异，控制计算复杂度）
# pcs_to_use <- 40
# message("最终选择的 PC 数量（基于 SCTransform 经验值）：", pcs_to_use)
#---------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 3.1：构建邻居图（FindNeighbors）
#-------------------------------------------------------------------------------

# 可选：从 .rds 文件加载 PCA 降维后的 Seurat 对象（跳过步骤 2.1 到 2.12）
# - 加载路径：processed_data_dir/scFlowKit_pca.rds
# - 确保 processed_data_dir 已定义
# sce_integrated <- readRDS(file = file.path(processed_data_dir, "scFlowKit_pca.rds"))

# - 基于 PCA 空间构建细胞间的邻居图
# - 使用 Seurat 的 FindNeighbors 函数，常用参数：
#   - reduction = "pca"：使用 PCA 降维结果
#   - dims = 1:20：使用 PCA 的前 20 个主成分
#   - k.param = 20：邻居数量（默认 20）
#   - verbose = TRUE：显示进度信息
# - 结果存储在 sce_integrated@graphs 中（包括 integrated_nn 和 integrated_snn）
message("步骤 3.1：寻找邻居...")
sce_integrated <- FindNeighbors(sce_integrated,
                                reduction = "pca",
                                dims = 1:20,  # 使用前 20 个主成分
                                k.param = 20,  # 邻居数量
                                verbose = TRUE)

# 输出邻居图信息
message("邻居图信息：")
print(names(sce_integrated@graphs))

#-------------------------------------------------------------------------------
# 步骤 3.2：聚类
#-------------------------------------------------------------------------------

# - 基于邻居图（SNN）对细胞进行聚类
# - 使用 Seurat 的 FindClusters 函数，常用参数：
#   - resolution = 0.8：分辨率（控制聚类数量）
#     - 高分辨率（> 1.0）：生成更多、更小的聚类，适合发现细粒度的细胞亚群
#     - 低分辨率（< 0.5）：生成更少、较大的聚类，适合发现大类细胞群
#   - algorithm = 1：使用原始 Louvain 算法
#   - verbose = TRUE：显示进度信息
# - 结果存储在 sce_integrated@meta.data$seurat_clusters 中

# 定义多个 resolution 值
resolutions <- c(0.4, 0.6, 0.8, 1.0, 1.4)
message("测试的 resolution 值：", paste(resolutions, collapse = ", "))

message("步骤 3.2：聚类...")
# 运行 FindClusters，测试多个 resolution 值
message("运行 FindClusters（测试多个 resolution 值）...")
sce_integrated <- FindClusters(sce_integrated,
                               resolution = resolutions,  # 分辨率，控制聚类数量
                               algorithm = 1,  # 使用原始 Louvain 算法
                               verbose = TRUE)

# 输出每个 resolution 的聚类数量和分布
for (res in resolutions) {
  col_name <- paste0("integrated_snn_res.", res)  # 修正列名
  message("Resolution ", res, " 聚类数量：")
  print(length(unique(sce_integrated@meta.data[[col_name]])))
  message("Resolution ", res, " 聚类分布：")
  print(table(sce_integrated@meta.data[[col_name]]))
}

# 最终选择 resolution = 0.8（默认值）
message("最终选择 resolution = 0.8 进行聚类...")
sce_integrated@meta.data$seurat_clusters <- sce_integrated@meta.data[["integrated_snn_res.0.8"]]

# 输出最终聚类数量和分布
message("最终聚类数量（resolution = 0.8）：")
print(length(unique(sce_integrated@meta.data$seurat_clusters)))
message("最终聚类分布（resolution = 0.8）：")
print(table(sce_integrated@meta.data$seurat_clusters))

# 保存聚类后的 Seurat 对象（中间点）
message("保存聚类后的 Seurat 对象...")
saveRDS(sce_integrated, file = file.path(processed_data_dir, "scFlowKit_clustered.rds"))
message("已保存至：", file.path(processed_data_dir, "scFlowKit_clustered.rds"))


#-------------------------------------------------------------------------------
# 步骤 3.3：运行 t-SNE 降维和可视化
#-------------------------------------------------------------------------------

# 可选：从 .rds 文件加载聚类后的 Seurat 对象（跳过步骤 2.1 到 3.2）
# - 加载路径：processed_data_dir/scFlowKit_clustered.rds
# - 确保 processed_data_dir 已定义
# sce_integrated <- readRDS(file = file.path(processed_data_dir, "scFlowKit_clustered.rds"))

# - 使用 t-SNE 进行降维，基于 PCA 空间
# - 使用 Seurat 的 RunTSNE 函数，常用参数：
#   - reduction = "pca"：使用 PCA 降维结果
#   - dims = 1:20：使用 PCA 的前 20 个主成分（与 FindNeighbors 一致）
#   - seed.use = 1：设置随机种子，确保结果可重复
#   - dim.embed = 2：降维后的维度（默认 2D）
#   - verbose = TRUE：显示进度信息
# - 结果存储在 sce_integrated@reductions$tsne 中
message("步骤 3.3：运行 t-SNE 降维...")
sce_integrated <- RunTSNE(sce_integrated,
               reduction = "pca",  # 使用 PCA 降维结果
               dims = 1:5,  # 使用前 5 个主成分
               seed.use = 1,  # 设置随机种子
               dim.embed = 2,  # 降维到 2D
               verbose = TRUE)  # 显示进度信息

# 输出 t-SNE 降维后的 Seurat 对象信息
message("t-SNE 降维后的 Seurat 对象基本信息：")
print(sce_integrated)
#  2 dimensional reductions calculated: pca, tsne

# 可视化 t-SNE 结果
# - 使用 DimPlot 绘制 t-SNE 散点图，展示 TSNE_1 和 TSNE_2 的分布
# - 按聚类结果（seurat_clusters）分组，观察聚类效果
message("可视化 t-SNE 结果...")

tsne_plot_clusters <- DimPlot(sce_integrated, 
                     reduction = "tsne",  # 使用 t-SNE 降维结果
                     group.by = "seurat_clusters",  # 按聚类结果分组
                     label = TRUE,  # 显示分组标签
                     repel = TRUE) +  # 避免标签重叠
  labs(title = "t-SNE Plot by Clusters") 

ggsave(file.path(output_dir, "figures/tsne_plot_clusters.png"), tsne_plot_clusters, width = 8, height = 6)
message("t-SNE 聚类图已保存至：", file.path(output_dir, "figures/tsne_plot_clusters.png"))

# 可视化 t-SNE 结果（按细胞周期阶段分组）
tsne_phase_plot <- DimPlot(sce_integrated, 
                           reduction = "tsne",  # 使用 t-SNE 降维结果
                           group.by = "Phase",  # 按细胞周期阶段分组
                           label = TRUE,  # 显示分组标签
                           repel = TRUE) +  # 避免标签重叠
  labs(title = "t-SNE Plot by Phase")

ggsave(file.path(output_dir, "figures/tsne_plot_phase.png"), tsne_phase_plot, width = 8, height = 6)
message("t-SNE 细胞周期图已保存至：", file.path(output_dir, "figures/tsne_plot_phase.png"))

# 可视化 t-SNE 结果（按样本分组）
tsne_plot_sample <- DimPlot(sce_integrated,
                            reduction = "tsne",  # 使用 t-SNE 降维结果
                            group.by = "sample",  # 按样本分组
                            label = TRUE,  # 显示分组标签
                            repel = TRUE) +  # 避免标签重叠
  labs(title = "t-SNE Plot by Sample")

ggsave(file.path(output_dir, "figures/tsne_plot_sample.png"), plot = tsne_plot_sample,
       width = 8, height = 6, dpi = 300)
message("t-SNE 样本图已保存至：", file.path(output_dir, "figures/tsne_plot_sample.png"))

# 按聚类分组，按样本分面
tsne_plot_clusters_split <- DimPlot(sce_integrated,
                                    reduction = "tsne",  # 使用 t-SNE 降维结果
                                    group.by = "seurat_clusters",  # 按聚类分组
                                    split.by = "sample",  # 按样本分面
                                    label = TRUE,  # 显示分组标签
                                    repel = TRUE) +  # 避免标签重叠
  labs(title = "t-SNE Plot by Clusters, Split by Sample")

# 保存 t-SNE 聚类分面图
ggsave(file.path(output_dir, "figures/tsne_plot_clusters_split_by_sample.png"),
       plot = tsne_plot_clusters_split,
       width = 12,  # 增加宽度以适应分面
       height = 6,
       dpi = 300)
message("t-SNE 聚类分面图已保存至：", file.path(output_dir, "figures/tsne_plot_clusters_split_by_sample.png"))  
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 3.4：运行 UMAP 降维和可视化
#-------------------------------------------------------------------------------

# - 使用 UMAP 进行降维，基于 PCA 空间和邻居图（SNN）
# - 使用 Seurat 的 RunUMAP 函数，常用参数：
#   - reduction = "pca" 使用 PCA 降维结果
#   - dims 使用 PCA 的主成分（默认 1:10）
#   - n.neighbors = 30 邻居数量（默认 30）
#   - min.dist = 0.3 最小距离（默认 0.3，控制点之间的距离）
#   - n.components = 2 降维后的维度（默认 2D）
#   - seed.use = 1 设置随机种子，确保结果可重复
#   - verbose = TRUE 显示进度信息
# - 结果存储在 sce@reductions$umap 中
message("步骤 3.4：运行 UMAP 降维...")
sce_integrated <- RunUMAP(sce_integrated,
               reduction = "pca",
               dims = 1:10,  # 使用前 10 个主成分
               n.neighbors = 30,  # 邻居数量
               min.dist = 0.3,  # 最小距离
               n.components = 2,  # 降维到 2D
               seed.use = 1,  # 设置随机种子
               verbose = TRUE)  # 显示进度信息

# 输出 UMAP 降维后的 Seurat 对象信息
message("UMAP 降维后的 Seurat 对象基本信息：")
print(sce_integrated)
#  3 dimensional reductions calculated: pca, tsne, umap

# 可视化 UMAP 结果
# - 使用 DimPlot 绘制 UMAP 散点图，展示 UMAP_1 和 UMAP_2 的分布
# - 按聚类结果（seurat_clusters）分组，观察聚类效果
message("可视化 UMAP 结果...")

umap_plot_clusters <- DimPlot(sce_integrated,
                     reduction = "umap",  # 使用 UMAP 降维结果
                     group.by = "seurat_clusters",  # 按聚类结果分组
                     label = TRUE,  # 显示分组标签
                     repel = TRUE) +  # 避免标签重叠
  labs(title = "UMAP Plot by Clusters")
ggsave(file.path(output_dir, "figures/umap_plot_clusters.png"), umap_plot_clusters, width = 8, height = 6)
message("UMAP 聚类图已保存至：", file.path(output_dir, "figures/umap_plot_clusters.png"))

# 可视化 UMAP 结果（按细胞周期阶段分组）
umap_plot_phase <- DimPlot(sce_integrated,
                           reduction = "umap",  # 使用 UMAP 降维结果
                           group.by = "Phase",  # 按细胞周期阶段分组
                           label = TRUE,  # 显示分组标签
                           repel = TRUE) +  # 避免标签重叠
  labs(title = "UMAP Plot by Phase")

ggsave(file.path(output_dir, "figures/umap_plot_phase.png"), umap_plot_phase, width = 8, height = 6)

# 按样本分组
umap_plot_sample <- DimPlot(sce_integrated,
                            reduction = "umap",  # 使用 UMAP 降维结果
                            group.by = "sample",  # 按样本分组
                            label = TRUE,  # 显示分组标签
                            repel = TRUE) +   # 避免标签重叠
  labs(title = "UMAP Plot by Sample")

ggsave(file.path(output_dir, "figures/umap_plot_sample.png"), plot = umap_plot_sample,
       width = 8, height = 6, dpi = 300)
message("UMAP 样本图已保存至：", file.path(output_dir, "figures/umap_plot_sample.png"))

# 按聚类分组，按样本分面
umap_plot_clusters_split <- DimPlot(sce_integrated,
                                    reduction = "umap",  # 使用 UMAP 降维结果
                                    group.by = "seurat_clusters",  # 按聚类分组
                                    split.by = "sample",  # 按样本分面
                                    label = TRUE,  # 显示分组标签
                                    repel = TRUE) +  # 避免标签重叠
  labs(title = "UMAP Plot by Clusters, Split by Sample")

# 保存 UMAP 聚类分面图
ggsave(file.path(output_dir, "figures/umap_plot_clusters_split_by_sample.png"),
       plot = umap_plot_clusters_split,
       width = 12,  # 增加宽度以适应分面
       height = 6,
       dpi = 300)
message("UMAP 聚类分面图已保存至：", file.path(output_dir, "figures/umap_plot_clusters_split_by_sample.png"))
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 3.5：保存聚类结果
#-------------------------------------------------------------------------------

# - 保存包含聚类和降维结果的 Seurat 对象为 Rds 文件
# - 文件路径：processed_data_dir/scFlowKit_umap.rds
# - 包含质控、过滤、标准化、可变基因选择、细胞周期评分、缩放、降维（PCA、t-SNE、UMAP）、双细胞去除和聚类的结果
message("步骤 3.5：保存聚类结果...")
dir.create(processed_data_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(sce_integrated, file = file.path(processed_data_dir, "scFlowKit_umap.rds"))
message("聚类结果已保存至：", file.path(processed_data_dir, "scFlowKit_umap.rds"))

# - 使用 Embeddings 提取降维结果并保存为 CSV 文件
# - 提取 PCA、t-SNE 和 UMAP 的降维坐标
# - 保存路径：processed_data_dir/
message("保存降维结果（PCA、t-SNE、UMAP）...")

# 提取 PCA 降维坐标
pca_embeddings <- Embeddings(sce_integrated, reduction = "pca")
write.csv(pca_embeddings, file = file.path(processed_data_dir, "scFlowKit_pca_embeddings.csv"), row.names = TRUE)

# 提取 t-SNE 降维坐标
tsne_embeddings <- Embeddings(sce_integrated, reduction = "tsne")
write.csv(tsne_embeddings, file = file.path(processed_data_dir, "scFlowKit_tsne_embeddings.csv"), row.names = TRUE)

# 提取 UMAP 降维坐标
umap_embeddings <- Embeddings(sce_integrated, reduction = "umap")
write.csv(umap_embeddings, file = file.path(processed_data_dir, "scFlowKit_umap_embeddings.csv"), row.names = TRUE)

message("降维结果已保存至：")
message("PCA 降维坐标：", file.path(processed_data_dir, "scFlowKit_pca_embeddings.csv"))
message("t-SNE 降维坐标：", file.path(processed_data_dir, "scFlowKit_tsne_embeddings.csv"))
message("UMAP 降维坐标：", file.path(processed_data_dir, "scFlowKit_umap_embeddings.csv"))


#-------------------------------------------------------------------------------
# 步骤 4.1：差异表达分析
#-------------------------------------------------------------------------------

# - SCTransform 标准化仅针对 3000 个变异最大的基因进行，主要用于降维和聚类
# - 我们感兴趣的许多基因可能并不存在于这些数据中，因此差异分析时需要切换回RNA assay
# - 寻找每个聚类的标志基因（marker genes）
# - 使用 Seurat 的 FindAllMarkers 函数，常用参数：
#   - test.use = "MAST" 使用 MAST 检验
#   - only.pos = TRUE 仅返回上调的基因
#   - min.pct = 0.25 基因在至少一组细胞中的最低表达比例
#   - logfc.threshold = 0.5 基因的最小 log2 折叠变化阈值
#   - verbose = TRUE 显示进度信息
# - 样本量过大时，可以通过抽样减少计算量
# - 例如：sce_sub = subset(sce, downsample = 100)，并在 FindAllMarkers 中使用 sce_sub
# - 结果为差异表达基因列表，包含以下字段：
#   - p_val：原始 p 值，表示基因在当前聚类与其他聚类之间的表达差异的显著性（未经多重检验校正）
#   - avg_log2FC：平均 log2 折叠变化，表示基因在当前聚类（ident.1）相对于其他聚类（ident.2）的表达差异（正值表示上调，负值表示下调）
#   - pct.1：基因在当前聚类（ident.1）中的表达比例（即表达该基因的细胞占当前聚类总细胞的比例）
#   - pct.2：基因在其他聚类（ident.2）中的表达比例（即表达该基因的细胞占其他聚类总细胞的比例）
#   - p_val_adj：调整后的 p 值（Bonferroni 校正），用于多重检验校正，控制假阳性率
#   - cluster：基因所属的聚类（例如 0, 1, 2, ...）
#   - gene：基因名称


# 可选：从 .rds 文件加载聚类后的 Seurat 对象（跳过步骤 2.1 到 3.5）
# - 加载路径：processed_data_dir/scFlowKit_umap.rds
# - 确保 processed_data_dir 已定义
sce_integrated <- readRDS(file = file.path(processed_data_dir, "scFlowKit_umap.rds"))

message("步骤 4.1：差异表达分析...")

# 切换到 RNA assay
message("切换到 RNA assay...")
DefaultAssay(sce_integrated) <- "RNA"
message("切换后的 Seurat 对象信息：")
print(sce_integrated)

# 合并 RNA assay 的 layers
message("合并 RNA assay 的 layers...")
sce_integrated <- JoinLayers(sce_integrated, assay = "RNA")
message("合并后的 Seurat 对象信息：")
print(sce_integrated)

# 重新标准化、选择高变基因和缩放数据（使用管道符）
message("重新标准化、选择高变基因和缩放 RNA assay...")
sce_integrated <- sce_integrated %>%
  NormalizeData(assay = "RNA", 
                normalization.method = "LogNormalize", 
                scale.factor = 10000,  # 默认值
                verbose = TRUE) %>%
  FindVariableFeatures(assay = "RNA",
                       selection.method = "vst",  # 使用 vst 方法
                       nfeatures = 2000,  # 选择 2000 个高变基因
                       verbose = TRUE) %>%
  ScaleData(assay = "RNA",
            features = rownames(sce_integrated),  # 使用所有基因（46517 个）
            verbose = TRUE)
message("标准化、高变基因选择和缩放后的 Seurat 对象信息：")
print(sce_integrated)

# 获取所有聚类标签
clusters <- unique(sce_integrated@meta.data$seurat_clusters)
message("聚类数量：", length(clusters))

# 使用 FindAllMarkers 分析所有聚类的标志基因
all_markers_df <- FindAllMarkers(sce_integrated,
                                 test.use = "MAST",  # 使用 MAST 检验
                                 only.pos = TRUE,  # 仅返回上调的基因
                                 min.pct = 0.25,  # 最低表达比例
                                 logfc.threshold = 0.5,  # 最小 log2 折叠变化
                                 verbose = TRUE)  # 显示进度信息

# 保存标志基因结果为 CSV 文件
message("保存标志基因结果...")
write.csv(all_markers_df, 
          file = file.path(output_dir, "tables", "scFlowKit_cluster_markers.csv"), 
          row.names = FALSE)
message("标志基因结果已保存至：", file.path(output_dir, "table", "scFlowKit_cluster_markers.csv"))

# 输出每个聚类的 top 5 标志基因
message("每个聚类的 top 5 标志基因：")
# - top 5 定义：按 log2 折叠变化（avg_log2FC）降序排序
top_markers <- all_markers_df %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()

# 按聚类分组打印 top 5 标志基因，并保存到文件
write.csv(top_markers, 
          file = file.path(output_dir, "tables", "scFlowKit_top5_markers.csv"), 
          row.names = FALSE)
message("Top 5 标志基因已保存至：", file.path(output_dir, "table", "scFlowKit_top5_markers.csv"))


# 按聚类分组打印 top 5 标志基因
for (cluster in unique(top_markers$cluster)) {
  message("聚类 ", cluster, "：")
  cluster_top <- top_markers[top_markers$cluster == cluster, ]
  print(cluster_top[, c("gene", "p_val_adj", "avg_log2FC", "pct.1", "pct.2")])
}

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 4.2：可视化标志基因
#-------------------------------------------------------------------------------

# - 可视化每个聚类的 top 标志基因，验证聚类结果
# - 使用 Seurat 的 FeaturePlot、VlnPlot 和 DotPlot 函数
# - FeaturePlot 在 UMAP 空间上绘制基因表达（散点图）
# - VlnPlot 绘制基因在不同聚类中的表达分布（小提琴图）
# - DotPlot 绘制基因在不同聚类中的表达比例和表达量（点图）
# - 可视化结果保存到子目录 results/figures/marker_visualization/
message("步骤 4.2：可视化标志基因...")

# 确保活跃 assay 是 RNA（差异表达分析基于 RNA assay）
DefaultAssay(sce_integrated) <- "RNA"

# 读取 top 5 标志基因文件（如果从头运行，可以直接使用 top_markers）
# top_markers <- read.csv(file.path(output_dir, "table", paste0(dataset_name, "_top5_markers.csv")))

# 创建子目录用于存储可视化结果
marker_viz_dir <- file.path(output_dir, "figures/marker_visualization")
dir.create(marker_viz_dir, recursive = TRUE, showWarnings = FALSE)

# 提取 top 基因列表（从 top_markers 中获取）
top_genes <- unique(top_markers$gene)
message("Top 标志基因数量：", length(top_genes))

# 使用 FeaturePlot 可视化 top 基因在 UMAP 空间的表达
message("绘制 FeaturePlot（UMAP 空间）...")
for (gene in top_genes) {
  p <- FeaturePlot(sce_integrated,
                   features = gene,
                   reduction = "umap",  # 使用 UMAP 降维结果
                   pt.size = 0.5,  # 调整点的大小，减少重叠
                   alpha = 0.7,  # 调整透明度，提高清晰度
                   label = TRUE,  # 显示基因名
                   repel = TRUE)  # 避免标签重叠
  ggsave(file.path(marker_viz_dir, paste0("umap_feature_", gene, ".png")), p, width = 8, height = 6, dpi = 300)
}

# 使用 VlnPlot 可视化 top 基因在不同聚类中的表达分布
message("绘制 VlnPlot（按聚类分组）...")
for (gene in top_genes) {
  p <- VlnPlot(sce_integrated,
               features = gene,
               group.by = "seurat_clusters",  # 按聚类分组
               pt.size = 0)  # 不显示单个细胞的点
  ggsave(file.path(marker_viz_dir, paste0("vlnplot_", gene, ".png")), p, width = 10, height = 6, dpi = 300)
}

# 使用 DotPlot 可视化 top 基因在不同聚类中的表达比例和表达量
message("绘制 DotPlot（按聚类分组）...")
# 使用 top_genes（78 个基因），每次 10 个基因一组
gene_groups <- split(top_genes, ceiling(seq_along(top_genes) / 10))
message("DotPlot 分组数量：", length(gene_groups))

# 按基因分组绘制 DotPlot，横坐标为簇，纵坐标为基因
for (i in seq_along(gene_groups)) {
  group_genes <- gene_groups[[i]]
  message("绘制 DotPlot 分组 ", i, "（基因：", paste(group_genes, collapse = ", "), "）...")
  p <- DotPlot(sce_integrated,
               features = group_genes,
               group.by = "seurat_clusters",  # 横坐标为簇
               dot.scale = 6) +  # 调整点的大小
    coord_flip() +  # 翻转坐标轴，横坐标为簇，纵坐标为基因
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 旋转横坐标标签（簇）
  ggsave(file.path(marker_viz_dir, paste0("dotplot_top_markers_group_", i, ".png")), p, width = 12, height = 6, dpi = 300)
}

message("标志基因可视化已完成，图表保存至：", marker_viz_dir)

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 4.3：寻找保守的标志基因（考虑样本条件）
#-------------------------------------------------------------------------------

# - 引入样本条件分组（例如 tumor vs normal），寻找保守的标志基因
# - 使用 Seurat 的 FindConservedMarkers 函数，常用参数：
#   - ident.1：目标聚类（例如 "0"）
#   - grouping.var：条件分组变量（例如 "condition"）
#   - test.use = "MAST"：使用 MAST 检验
#   - only.pos = TRUE：仅返回上调的基因
#   - min.pct = 0.25：基因在至少一组细胞中的最低表达比例
#   - logfc.threshold = 0.5：基因的最小 log2 折叠变化阈值
#   - verbose = TRUE：显示进度信息
# - 结果为保守的标志基因列表，包含以下字段：
#   - gene：基因名称
#   - p_val：原始 p 值（每个条件分别计算）
#   - avg_log2FC：平均 log2 折叠变化（每个条件分别计算）
#   - pct.1：基因在目标聚类中的表达比例（每个条件分别计算）
#   - pct.2：基因在其他聚类中的表达比例（每个条件分别计算）
#   - p_val_adj：调整后的 p 值（Bonferroni 校正）
#   - max_pval：所有条件中的最大 p 值
#   - minimump_p_val：所有条件中的最小 p 值（综合 p 值）
#   - cluster：基因所属的聚类
message("步骤 4.3：寻找保守的标志基因（考虑样本条件）...")

# 确保活跃 assay 是 RNA（差异表达分析基于 RNA assay）
DefaultAssay(sce_integrated) <- "RNA"

# 引入条件分组（假设 donor1 和 donor2 是 tumor，donor3 和 donor4 是 normal）
message("引入条件分组（tumor vs normal）...")
sce_integrated@meta.data$condition <- ifelse(sce_integrated@meta.data$sample %in% c("5k_pbmc_donor1", "5k_pbmc_donor2"), "condition1", "condition2")
message("条件分组分布：")
print(table(sce_integrated@meta.data$condition))

# 获取所有聚类标签
clusters <- unique(sce_integrated@meta.data$seurat_clusters)
message("聚类数量：", length(clusters))

# 使用 FindConservedMarkers 寻找每个聚类的保守标志基因
message("运行 FindConservedMarkers 寻找保守标志基因...")
conserved_markers_list <- list()
for (cluster in clusters) {
  message("寻找聚类 ", cluster, " 的保守标志基因...")
  markers <- FindConservedMarkers(sce_integrated,
                                  ident.1 = cluster,  # 目标聚类
                                  grouping.var = "condition",  # 按条件分组（tumor vs normal）
                                  test.use = "MAST",  # 使用 MAST 检验
                                  only.pos = TRUE,  # 仅返回上调的基因
                                  min.pct = 0.25,  # 最低表达比例
                                  logfc.threshold = 0.5,  # 最小 log2 折叠变化
                                  verbose = TRUE)
  if (nrow(markers) > 0) {
    # 将行名（基因名）转换为 gene 列
    markers <- markers %>%
      tibble::rownames_to_column(var = "gene")
    markers$cluster <- cluster
    # 重置行名，避免合并时添加前缀
    rownames(markers) <- NULL
    conserved_markers_list[[as.character(cluster)]] <- markers
  }
}

# 合并所有聚类的保守标志基因
message("合并所有聚类的保守标志基因...")
conserved_markers_df <- do.call(rbind, conserved_markers_list)

# 保存保守标志基因结果为 CSV 文件
message("保存保守标志基因结果...")
conserved_markers_file <- file.path(output_dir, "tables", "scFlowKit_conserved_markers.csv")
write.csv(conserved_markers_df, file = conserved_markers_file, row.names = FALSE)
message("保守标志基因结果已保存至：", conserved_markers_file)

# 输出每个聚类的 top 5 保守标志基因
message("提取每个聚类的 top 5 保守标志基因...")
# 动态选择所有以 _avg_log2FC 结尾的列
log2fc_cols <- grep("_avg_log2FC$", names(conserved_markers_df), value = TRUE)

# 计算所有条件的 avg_log2FC 均值
conserved_markers_df <- conserved_markers_df %>%
  dplyr::mutate(mean_log2FC = rowMeans(dplyr::select(., all_of(log2fc_cols))))

top_conserved_markers <- conserved_markers_df %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(desc(mean_log2FC)) %>% # 按 avg_log2FC 均值降序排序
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()

# 保存 top 5 保守标志基因到文件
top_conserved_markers_file <- file.path(output_dir, "tables", "scFlowKit_top5_conserved_markers.csv")
write.csv(top_conserved_markers, file = top_conserved_markers_file, row.names = FALSE)
message("Top 5 保守标志基因已保存至：", top_conserved_markers_file)

# 按聚类分组打印 top 5 保守标志基因
message("打印每个聚类的 top 5 保守标志基因：")
# 动态选择打印字段：gene, minimump_p_val, 以及所有 _avg_log2FC 列
print_cols <- c("gene", "minimump_p_val", log2fc_cols)
for (cluster in unique(top_conserved_markers$cluster)) {
  message("聚类 ", cluster, "：")
  cluster_top <- top_conserved_markers[top_conserved_markers$cluster == cluster, ]
  print(cluster_top[, print_cols])
}


#-------------------------------------------------------------------------------
# 步骤 4.4：比较任意聚类的差异表达基因
#-------------------------------------------------------------------------------

# - 允许用户指定一个或多个聚类，比较它们之间的差异表达基因
# - 比如聚类0，2，3均被确定为T细胞，可以通过进一步比较用于区分，找到更精细的亚类
# - 或许能进一步区分，比如找到聚类0或许还主要表达初始T细胞的marker，而2，3更接近成熟T细胞
# - 使用 Seurat 的 FindMarkers 函数，常用参数：
#   - ident.1：第一个分组（例如 c("0", "1")）
#   - ident.2：第二个分组（例如 c("2", "3")）
#   - test.use = "MAST"：使用 MAST 检验
#   - only.pos = TRUE：仅返回上调的基因
#   - min.pct = 0.25：基因在至少一组细胞中的最低表达比例
#   - logfc.threshold = 0.5：基因的最小 log2 折叠变化阈值
#   - verbose = TRUE：显示进度信息
# - 结果为差异表达基因列表，包含以下字段：
#   - gene：基因名称
#   - p_val：原始 p 值
#   - avg_log2FC：平均 log2 折叠变化
#   - pct.1：基因在 ident.1 中的表达比例
#   - pct.2：基因在 ident.2 中的表达比例
#   - p_val_adj：调整后的 p 值（Bonferroni 校正）

message("步骤 4.4：比较任意聚类的差异表达基因...")

# 确保活跃 assay 是 RNA（差异表达分析基于 RNA assay）
DefaultAssay(sce_integrated) <- "RNA"

# 获取所有聚类标签
clusters <- unique(sce_integrated@meta.data$seurat_clusters)
message("可用聚类：", paste(clusters, collapse = ", "))

# 用户指定要比较的两个分组（支持一个或多个聚类）
clusters1 <- c("0")  # 第一个分组（例如聚类 0 ）
clusters2 <- c("2", "3")  # 第二个分组（例如聚类 2 和 3）
message("比较的分组：", paste(clusters1, collapse = ","), " vs ", paste(clusters2, collapse = ","))

# 使用 FindMarkers 比较两个分组
message("运行 FindMarkers 比较分组 ", paste(clusters1, collapse = ","), " vs ", paste(clusters2, collapse = ","), "...")
markers <- FindMarkers(sce_integrated,
                       ident.1 = clusters1,  # 第一个分组
                       ident.2 = clusters2,  # 第二个分组
                       subset.ident = c(clusters1, clusters2),  # 限制分析的聚类
                       test.use = "MAST",  # 使用 MAST 检验
                       only.pos = TRUE,  # 仅返回上调的基因
                       min.pct = 0.25,  # 最低表达比例
                       logfc.threshold = 0.5,  # 最小 log2 折叠变化
                       verbose = TRUE)

# 将行名（基因名）转换为 gene 列
markers <- markers %>%
  tibble::rownames_to_column(var = "gene")

# 保存差异表达基因结果为 CSV 文件
message("保存差异表达基因结果...")
comparison_name <- paste0("cluster_", paste(clusters1, collapse = "_"), "_vs_", paste(clusters2, collapse = "_"))
markers_file <- file.path(output_dir, "tables", paste0("scFlowKit_markers_", comparison_name, ".csv"))
write.csv(markers, file = markers_file, row.names = FALSE)
message("差异表达基因结果已保存至：", markers_file)

# 输出 top 5 差异表达基因（按 avg_log2FC 排序）
message("提取 top 5 差异表达基因（按 avg_log2FC 排序）...")
top_markers <- markers %>%
  dplyr::arrange(desc(avg_log2FC)) %>%  # 按 avg_log2FC 降序排序
  dplyr::slice_head(n = 5)

# 保存 top 5 差异表达基因到文件
top_markers_file <- file.path(output_dir, "tables", paste0("scFlowKit_top5_markers_", comparison_name, ".csv"))
write.csv(top_markers, file = top_markers_file, row.names = FALSE)
message("Top 5 差异表达基因已保存至：", top_markers_file)

# 打印 top 5 差异表达基因
message("打印 top 5 差异表达基因：")
print_cols <- c("gene", "p_val_adj", "avg_log2FC", "pct.1", "pct.2")
print(top_markers[, print_cols])

# 步骤 4.5：细胞注释
#-------------------------------------------------------------------------------
# 继续后续分析

# SingleR注释

# 差异基因注释

# 已知marker + 点图注释




# 步骤 4.6：可视化
#-------------------------------------------------------------------------------
# 继续后续分析

# 基因，metadata

# ridgeplot
# vlnplot
# featureplot
# dotplot
# doheatmap


# 步骤 4.7：Condition之间比较
#-------------------------------------------------------------------------------
# 继续后续分析

# 不同条件下同一细胞类型之间的差异


# 步骤 4.7：差异基因富集分析 & 通路打分
#-------------------------------------------------------------------------------
# 继续后续分析

# 步骤 5.1：轨迹分析
#-------------------------------------------------------------------------------
# 继续后续分析

# 步骤 5.2：snRNA
#-------------------------------------------------------------------------------
# 继续后续分析

# 步骤 5.2：整合蛋白丰度
#-------------------------------------------------------------------------------
# 继续后续分析

