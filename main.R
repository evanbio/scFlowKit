#-------------------------------------------------------------------------------
# scFlowKit: A Modular Single-Cell RNA-seq Analysis Pipeline
#-------------------------------------------------------------------------------

# 项目信息
# Name: scFlowKit
# Description: 一个模块化的单细胞 RNA-seq 分析流程，支持 R 和 Python
# GitHub: https://github.com/chriswang001121/scFlowKit
# License: MIT License
# Version: [占位，待后续补充]

# 作者信息
# Name: chriswang001121
# Email: chriswang001121@gmail.com
# Homepage: https://scholar.pulppoetry.org
# LinkedIn: https://www.linkedin.com/in/yibin-zhou


# 当前版本功能
# - 数据加载：支持 10X Genomics 数据格式（.h5 或 .mtx）

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

# 数据路径：指定输入数据的基础路径，告诉程序从哪里加载单细胞数据
# - 应指向 data/raw/ 文件夹，具体数据集子路径（如 pbmc3k）由模块函数指定
# - 示例：如果你的数据在 data/raw/pbmc3k/ 文件夹，保持此值为 "data/raw/"
data_path <- "data/raw/"

# 输出目录：指定分析结果的保存路径，告诉程序把结果保存在哪里
# - 分析结果（比如图表、表格）会保存在此目录下的 figures/ 和 tables/ 子文件夹
# - 示例：如果想保存在 results/run1/ 文件夹，设置为 "results/run1/"
output_dir <- "results/"

# 预处理数据保存路径：指定预处理后数据的保存路径
# - 预处理数据会保存在此目录下
# - 示例：data/processed/
processed_data_dir <- "data/processed/"

# 设置系统报错语言为英文，便于调试和阅读
# - 在非英文环境下运行 R 时，确保错误信息以英文显示
Sys.setenv(LANGUAGE = "en")

# 禁止字符向量自动转换为因子，提升数据处理灵活性
# - 处理单细胞数据时，避免字符型元数据（如细胞名）被转为因子
options(stringsAsFactors = FALSE)


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

# 指定数据集名称（用户可根据实际数据集调整）
dataset_name <- "5k_pbmc_donor4"

# 加载数据
message("步骤 1：正在加载单细胞 RNA-seq 数据...")
sce <- load_data(base_path = data_path, 
                 dataset_name = dataset_name,
                 min_cells = 3,          # 基因至少在 3 个细胞中表达
                 min_features = 40,      # 细胞至少表达 40 个基因
                 project = "scFlowKit",  # Seurat 对象项目名称
                 assay = "RNA")          # 测序类型

# 输出基本信息，确认加载成功
message("Seurat 对象基本信息：")
print(sce)
# 加载结果（2025年3月24日）：
# An object of class Seurat 
# 25554 features across 5709 samples within 1 assay 
# Active assay: RNA (25554 features, 0 variable features)
#   1 layer present: counts


#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 数据集介绍
#-------------------------------------------------------------------------------

# 数据集信息
# - 数据集：5k Human PBMCs (Donor 4)
# - 来源：10x Genomics 公开数据集
# - 描述：Universal 3' Gene Expression dataset analyzed using Cell Ranger 9.0.0
# - 下载链接：https://www.10xgenomics.com/datasets/5k-human-pbmcs-donor-4-universal-3-prime-gene-expression-dataset-analyzed-using-cell-ranger-9-0-0
# - 实验背景：外周血单核细胞（PBMC）取自一名健康的女性捐献者（年龄 36-50 岁），由 10x Genomics 从 Cellular Technologies Limited 获取。
#            使用 Chromium GEM-X Single Cell 3' Reagent Kits v4 和 Illumina NovaSeq 6000 测序，平均每个细胞约 36,000 个读对。
# - 测序参数：
#   - Read 1：28 个周期
#   - i7 索引：10 个周期
#   - i5 索引：10 个周期
#   - Read 2：90 个周期
# - 细胞数量：约 5,000 个外周血单核细胞（PBMC）
# - 数据格式：10x Genomics 标准格式（matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz）
# - 数据路径：data/raw/5k_pbmc_donor4/
# - 许可：Creative Commons Attribution 4.0 International (CC BY 4.0)
# - 使用说明：
#   - 数据为过滤后的基因表达矩阵（filtered feature barcode matrix），已去除低质量条码，可直接用于下游分析。
#   - 文件为压缩格式（.gz），Seurat 的 Read10X() 函数可直接加载，无需解压。
# - 补充文件（存放在 docs/5k_pbmc_donor4/ 文件夹）

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 1.5：探索 Seurat 对象结构（注释掉，仅供学习参考）
#-------------------------------------------------------------------------------

# # 查看 Seurat 对象的整体结构
# str(sce)

# # 提取 assay 数据（RNA 表达矩阵）
# # - sce[["RNA"]] 提取 RNA assay
# # - sce[["RNA"]]$counts 提取原始计数矩阵
# # - sce[["RNA"]]$data 提取标准化后的数据（当前为空，因为尚未标准化）
# rna_counts <- sce[["RNA"]]$counts                    # 方法 1：通过层级索引提取
# rna_counts <- GetAssayData(sce,                      # 方法 2：通过函数提取
#                            assay = "RNA",
#                            layer = "counts")
# # 查看前几个基因和细胞的表达量
# rna_counts[1:5, 1:5]

# # 提取 metadata（元数据）
# # - sce@meta.data 包含细胞的元数据，如细胞 ID、每个细胞的基因数等
# metadata <- sce@meta.data                            # 方法 1：通过层级索引提取
# metadata <- colData(sce)                             # 方法 2：通过函数提取
# metadata <- sce[[]]                                  # 方法 3：Seurat 4.0.0 之后的版本可以直接提取，默认提取 metadata
# # 查看元数据的列名
# colnames(metadata)
# # 查看前几行元数据
# head(metadata)

# # 查看当前激活的 assay
# active_assay <- DefaultAssay(sce)
# message("当前激活的 assay：", active_assay)

# # 查看 Seurat 对象的基因名（features）
# features <- rownames(sce@assays$RNA)                 # 方法 1：通过层级索引提取
# features <- rownames(sce)                            # 方法 2：通过函数提取
# # 查看前几个基因名
# head(features)

# # 查看 Seurat 对象的细胞名（samples）
# cells <- colnames(sce@assays$RNA)                    # 方法 1：通过层级索引提取
# cells <- colnames(sce)                               # 方法 2：通过函数提取
# # 查看前几个细胞名
# head(cells)

# # 查看 Seurat 对象的层（layers）
# # - 当前只有 counts 层，标准化后会有 data 层
# layers <- SeuratObject::Layers(sce)
# message("当前数据层：", paste(layers, collapse = ", "))

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 2：数据预处理和质控
#-------------------------------------------------------------------------------

# 导入预处理模块
source("Rutils/preprocess.R")

# 步骤 2.1：计算质控指标
#-------------------------------------------------------------------------------
# - 计算线粒体基因比例、总 UMI 计数、基因数
# - 可选计算红细胞和核糖体基因比例（默认不计算）
# - 可通过 hb_genes 和 ribo_genes 参数传入自定义基因集，例如：
#   - hb_genes = c("HBB", "HBA1", "HBA2") 指定红细胞基因
#   - ribo_genes = c("RPL5", "RPS6", "RPL10") 指定核糖体基因
# - 注意：函数会检查 sce 的基因名大小写模式（全大写、全小写、首字母大写），并转换基因集以匹配
# - 如果基因名大小写不匹配或基因不存在，会发出警告，建议用户检查基因名并手动调整

message("步骤 2.1：计算质控指标...")
sce <- calculate_qc_metrics(sce,
                            calculate_hb = FALSE,    # 不计算红细胞基因比例
                            calculate_ribo = FALSE)  # 不计算核糖体基因比例

# 查看计算后的元数据
# - 包含 nCount_RNA、nFeature_RNA、percent_mito 等指标
message("质控指标计算结果（前几行元数据）：")
head(sce@meta.data)
#                    orig.ident nCount_RNA nFeature_RNA percent_mito
# AAACCAAAGGCGCTTG-1  scFlowKit      10084         2847     9.044030
# AAACCAAAGTAGGACG-1  scFlowKit       9336         3108     3.309769
# AAACCATTCCAGCTAA-1  scFlowKit       7089         2751     3.879250
# AAACCATTCCATCCGC-1  scFlowKit       6664         2421     6.692677
# AAACCATTCGACCAGT-1  scFlowKit      11082         3244     2.752211
# AAACCCTGTGATGAAT-1  scFlowKit      13111         3590     4.050034

#-------------------------------------------------------------------------------

# 步骤 2.2：可视化质控指标，帮助用户检查数据质量
# - 绘制基因数、总 UMI 计数、线粒体比例的分布图（拼成一张图）
# - 使用 VlnPlot 绘制小提琴图，保存为 output_dir/figures/qc_metrics_combined.png
# - 使用 FeatureScatter 绘制散点图，保存为 output_dir/figures/qc_metrics_scatter_combined.png
# - 可通过 pt.size 参数控制 VlnPlot 中点的大小（默认 0，不显示点）
message("步骤 2.2：可视化质控指标...")
plot_qc_metrics(sce,
                output_dir = output_dir,
                pt.size = 0)  # 默认不显示点，可调整为 0.1 等

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 步骤 2.3：过滤低质量细胞和基因
#-------------------------------------------------------------------------------
# - 过滤标准：细胞基因数、总 UMI 计数、线粒体比例
# - 可选过滤红细胞和核糖体基因比例（默认不启用）
# - 用户可通过参数调整过滤阈值，例如：
#   - min_umi = 500, max_umi = 15000 指定 UMI 计数范围
#   - min_genes = 200, max_genes = 5000 指定基因数范围
#   - max_mito = 10 指定最大线粒体基因比例
message("步骤 2.3：过滤低质量细胞和基因...")
sce <- filter_cells_genes(sce,
                          min_umi = 500,           # 最小 UMI 计数
                          max_umi = 15000,         # 最大 UMI 计数
                          min_genes = 200,         # 最小基因数
                          max_genes = 5000,        # 最大基因数
                          max_mito = 10,           # 最大线粒体基因比例（百分比）
                          filter_hb = FALSE,       # 不过滤红细胞基因比例
                          max_hb = 5,              # 最大红细胞基因比例（未启用）
                          filter_ribo = FALSE,     # 不过滤核糖体基因比例
                          max_ribo = 50)           # 最大核糖体基因比例（未启用）

# 输出过滤后的 Seurat 对象信息
message("过滤后的 Seurat 对象基本信息：")
print(sce)

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 2.4：标准化和对数化数据
#-------------------------------------------------------------------------------
# - 标准化：将每个细胞的总表达量缩放到 10,000
# - 对数化：log1p 变换，稳定数据分布
# - 使用 Seurat 的 NormalizeData 函数，参数：
#   - normalization.method = "LogNormalize" 指定标准化方法，可选 "LogNormalize" 或 "CLR"
#   - scale_factor = 10000 指定缩放因子
message("步骤 2.4：标准化和对数化数据...")
sce <- NormalizeData(sce,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000)  # 默认缩放因子为 10000

# 提示用户标准化完成
message("数据标准化和对数化完成！")

# 输出标准化后的 Seurat 对象信息
message("标准化后的 Seurat 对象基本信息：")
print(sce) # 标准化后的 Seurat 对象会增加一个 data 层

# 检查数据层，确认标准化结果
message("标准化后的数据层：")
print(SeuratObject::Layers(sce))

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 步骤 2.5：保存预处理数据（标准化后）
#-------------------------------------------------------------------------------
# - 保存预处理后的 Seurat 对象为 Rds 文件，路径：processed_data_dir/
# - 文件名格式：dataset_name_normalized.rds
# - 包含质控、过滤和标准化的结果
# - 后续步骤（可变基因选择、缩放、降维、聚类等）会继续处理此数据
message("步骤 2.5：保存预处理数据（标准化后）...")
dir.create(processed_data_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(sce, file = paste0(processed_data_dir, "/", dataset_name, "_normalized.rds"))
message("预处理数据已保存至：", paste0(processed_data_dir, "/", dataset_name, "_normalized.rds"))

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 2.6：寻找可变基因（可选：从 .rds 文件加载数据）
#-------------------------------------------------------------------------------

# 可选：从 .rds 文件加载预处理后的 Seurat 对象（跳过步骤 2.1 到 2.5）
# - 加载路径：processed_data_dir/dataset_name_normalized.rds
# - 确保 processed_data_dir 和 dataset_name 已定义
# sce <- readRDS(file = paste0(processed_data_dir, "/", dataset_name, "_normalized.rds"))

# - 识别高变异基因（highly variable genes），用于后续降维和聚类
# - 使用 Seurat 的 FindVariableFeatures 函数，参数：
#   - selection.method = "vst" 使用 variance stabilizing transformation 方法
#   - nfeatures = 2000 选择 2000 个高变异基因
#   - verbose = TRUE 显示进度条
# - 可变基因存储在 sce@assays$RNA@var.features 中
# - 基因的均值和方差存储在 sce@assays$RNA@meta.data 中（Seurat 5.0 及以上版本）
message("步骤 2.6：寻找可变基因...")
sce <- FindVariableFeatures(sce,
                            selection.method = "vst",
                            nfeatures = 2000,  # 选择 2000 个高变异基因
                            verbose = TRUE)  # 显示进度条

# 提取全部可变基因名称（使用 VariableFeatures）
all_variable_genes <- VariableFeatures(sce)
message("找到的可变基因数量：", length(all_variable_genes))

# 提取可变基因名称，打印变化最大的 10 个基因
top_variable_genes <- head(all_variable_genes, 10)
message("变化最大的 10 个可变基因：", paste(top_variable_genes, collapse = ", "))

# 可视化可变基因
# - 使用 VariableFeaturePlot 绘制均值-方差散点图，常用参数：
#   - log = NULL,  # 默认根据数据决定
#   - col = c("black", "red") 指定非高变异基因和高变异基因的颜色
#   - pt.size = 1 指定点的大小
# - 使用 LabelPoints 标注变化最大的 10 个基因，常用参数：
#   - points 要标注的基因名称，字符向量
#   - repel = TRUE 避免标签重叠
#   - xnudge/ynudge 调整标签位置（默认 0）
# - 保存为 output_dir/figures/variable_features_plot.png
message("绘制可变基因散点图...")
variable_feature_plot <- VariableFeaturePlot(sce,
                                             log = NULL,  # 默认根据数据决定
                                             col = c("black", "red"),  # 非高变异基因黑色，高变异基因红色
                                             pt.size = 1)  # 点的大小
variable_feature_plot <- LabelPoints(plot = variable_feature_plot,
                                     points = top_variable_genes,
                                     repel = TRUE,  # 避免标签重叠
                                     xnudge = 0.3,  # 标签水平偏移
                                     ynudge = 0.05)  # 标签垂直偏移
ggsave(file.path(output_dir, "figures/variable_features_plot.png"), variable_feature_plot, width = 8, height = 6)

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 2.7：细胞周期评分
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
message("步骤 2.7：细胞周期评分...")
message("S 期基因数量：", length(cc.genes$s.genes))
message("G2/M 期基因数量：", length(cc.genes$g2m.genes))

sce <- CellCycleScoring(sce,
                        s.features = cc.genes$s.genes,
                        g2m.features = cc.genes$g2m.genes,
                        set.ident = FALSE)

# 输出细胞周期评分结果（前几行元数据）
message("细胞周期评分结果（前几行元数据）：")
print(head(sce@meta.data[, c("S.Score", "G2M.Score", "Phase")]))

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 2.8：数据缩放
#-------------------------------------------------------------------------------
# - 对基因进行中心化和缩放（零均值、单位方差）
# - 使用 Seurat 的 ScaleData 函数，常用参数：
#   - features 要缩放的基因，默认 NULL（使用高变异基因 VariableFeatures(sce)）
#   - vars.to.regress 回归掉的变量（比如 "percent_mito" 去除线粒体比例影响）
#   - scale.max 缩放后表达量的最大值（默认 10）
#   - do.scale 是否进行缩放（默认 TRUE）
#   - do.center 是否进行中心化（默认 TRUE）
#   - verbose = TRUE 显示进度信息
# - 内存允许时，建议使用全部基因（features = rownames(sce)），以保留所有基因信息
# - 内存有限时，使用高变异基因（features = VariableFeatures(sce)），以减少计算量
# - 函数默认：features = NULL 时，ScaleData 使用高变异基因
# - 这里选择使用全部基因
# - 结果存储在 sce[["RNA"]]$scale.data 中
message("步骤 2.8：数据缩放...")
sce <- ScaleData(sce,
                 features = rownames(sce),  # 使用全部基因（内存允许时推荐）
                 vars.to.regress = NULL,  # 不回归任何变量（可设置为 c("S.Score", "G2M.Score")）
                 scale.max = 10,  # 缩放后表达量最大值
                 do.scale = TRUE,  # 进行缩放
                 do.center = TRUE,  # 进行中心化
                 verbose = TRUE)  # 显示进度信息

# 输出缩放后的 Seurat 对象信息
message("缩放后的 Seurat 对象基本信息：")
print(sce)
# Active assay: RNA (25554 features, 2000 variable features)
# 3 layers present: counts, data, scale.data

# 打印 scale.data 的前 5 个基因和前 5 个细胞
message("scale.data 前 5 个基因和前 5 个细胞：")
print(sce[["RNA"]]$scale.data[1:5, 1:5])

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 2.9：降维（PCA）
#-------------------------------------------------------------------------------

# - 使用 PCA 进行降维，基于高变异基因
# - 使用 Seurat 的 RunPCA 函数，常用参数：
#   - features 使用高变异基因（默认 VariableFeatures(sce)）
#   - npcs = 50 选择前 50 个主成分（可调整）
#   - verbose = TRUE 显示进度信息
# - 结果存储在 sce@reductions$pca 中
message("步骤 2.9：降维（PCA）...")
sce <- RunPCA(sce,
              features = VariableFeatures(sce),  # 使用高变异基因（默认）
              npcs = 50,  # 计算前 50 个主成分
              verbose = TRUE)  # 显示进度信息

# 输出降维后的 Seurat 对象信息
message("降维后的 Seurat 对象基本信息：")
print(sce)
# An object of class Seurat 
# 25554 features across 4413 samples within 1 assay 
# Active assay: RNA (25554 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 1 dimensional reduction calculated: pca

# 可视化 PCA 结果
# - 使用 DimPlot 绘制 PCA 散点图，展示 PC1 和 PC2 的分布
# - 按细胞周期阶段（Phase）分组，观察细胞周期对 PCA 分布的影响
# - 使用 ElbowPlot 展示主成分的方差解释比例
message("可视化 PCA 结果...")
pca_plot <- DimPlot(sce, 
                    reduction = "pca",  # 使用 PCA 降维结果
                    dims = c(1, 2),  # 绘制 PC1 和 PC2
                    group.by = "Phase",  # 按细胞周期阶段分组
                    label = TRUE,  # 显示分组标签
                    repel = TRUE)  # 避免标签重叠
ggsave(file.path(output_dir, "figures/pca_plot.png"), pca_plot, width = 8, height = 6)

# 绘制 ElbowPlot 展示主成分的方差解释比例
elbow_plot <- ElbowPlot(sce, 
                        reduction = "pca",  # 使用 PCA 降维结果
                        ndims = 50)  # 显示前 50 个主成分
ggsave(file.path(output_dir, "figures/pca_elbow_plot.png"), elbow_plot, width = 8, height = 6)

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 2.10：去除双细胞（使用 DoubletFinder）
#-------------------------------------------------------------------------------

# - 使用 DoubletFinder 检测和去除双细胞
# - 假设双细胞比例为 8%（可根据实验设计调整）
# - 使用 PCA 结果（PCs）进行双细胞检测
# - 封装在 remove_doublets 函数中（位于 Rutils/preprocess.R）
message("步骤 2.10：去除双细胞...")
sce <- remove_doublets(sce,
                       PCs = 1:20,  # 使用前 20 个主成分
                       doublet_rate = 0.08,  # 假设双细胞比例为 8%
                       pN = 0.25,  # DoubletFinder 参数 pN
                       sct = FALSE)  # 不使用 SCTransform 数据

# 输出过滤双细胞后的 Seurat 对象信息
message("过滤双细胞后的 Seurat 对象基本信息：")
print(sce)

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 2.11：重新寻找可变基因（可选：从 .rds 文件加载数据）
#-------------------------------------------------------------------------------

# 可选：从 .rds 文件加载预处理后的 Seurat 对象（跳过步骤 2.1 到 2.5）
# - 加载路径：processed_data_dir/dataset_name_normalized.rds
# - 确保 processed_data_dir 和 dataset_name 已定义
# sce <- readRDS(file = paste0(processed_data_dir, "/", dataset_name, "_normalized.rds"))

# - 识别高变异基因（highly variable genes），用于后续降维和聚类
# - 使用 Seurat 的 FindVariableFeatures 函数，参数：
#   - selection.method = "vst" 使用 variance stabilizing transformation 方法
#   - nfeatures = 2000 选择 2000 个高变异基因
#   - verbose = TRUE 显示进度条
# - 可变基因存储在 sce@assays$RNA@var.features 中
# - 基因的均值和方差存储在 sce@assays$RNA@meta.data 中（Seurat 5.0 及以上版本）
message("步骤 2.11：重新寻找可变基因...")
sce <- FindVariableFeatures(sce,
                            selection.method = "vst",
                            nfeatures = 2000,  # 选择 2000 个高变异基因
                            verbose = TRUE)  # 显示进度条

# 提取全部可变基因名称（使用 VariableFeatures）
all_variable_genes <- VariableFeatures(sce)
message("找到的可变基因数量：", length(all_variable_genes))

# 提取可变基因名称，打印变化最大的 10 个基因
top_variable_genes <- head(all_variable_genes, 10)
message("变化最大的 10 个可变基因：", paste(top_variable_genes, collapse = ", "))

# 可视化可变基因
# - 使用 VariableFeaturePlot 绘制均值-方差散点图，常用参数：
#   - log = NULL,  # 默认根据数据决定
#   - col = c("black", "red") 指定非高变异基因和高变异基因的颜色
#   - pt.size = 1 指定点的大小
# - 使用 LabelPoints 标注变化最大的 10 个基因，常用参数：
#   - points 要标注的基因名称，字符向量
#   - repel = TRUE 避免标签重叠
#   - xnudge/ynudge 调整标签位置（默认 0）
# - 保存为 output_dir/figures/variable_features_plot_after_doublet_removal.png
message("绘制可变基因散点图...")
variable_feature_plot <- VariableFeaturePlot(sce,
                                             log = NULL,  # 默认根据数据决定
                                             col = c("black", "red"),  # 非高变异基因黑色，高变异基因红色
                                             pt.size = 1)  # 点的大小
variable_feature_plot <- LabelPoints(plot = variable_feature_plot,
                                     points = top_variable_genes,
                                     repel = TRUE,  # 避免标签重叠
                                     xnudge = 0.3,  # 标签水平偏移
                                     ynudge = 0.05)  # 标签垂直偏移
ggsave(file.path(output_dir, "figures/variable_features_plot_after_doublet_removal.png"), variable_feature_plot, width = 8, height = 6)

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 2.12：重新细胞周期评分
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
message("步骤 2.12：重新细胞周期评分...")
message("S 期基因数量：", length(cc.genes$s.genes))
message("G2/M 期基因数量：", length(cc.genes$g2m.genes))

sce <- CellCycleScoring(sce,
                        s.features = cc.genes$s.genes,
                        g2m.features = cc.genes$g2m.genes,
                        set.ident = FALSE)

# 输出细胞周期评分结果（前几行元数据）
message("细胞周期评分结果（前几行元数据）：")
print(head(sce@meta.data[, c("S.Score", "G2M.Score", "Phase")]))

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 2.13：重新数据缩放
#-------------------------------------------------------------------------------
# - 对基因进行中心化和缩放（零均值、单位方差）
# - 使用 Seurat 的 ScaleData 函数，常用参数：
#   - features 要缩放的基因，默认 NULL（使用高变异基因 VariableFeatures(sce)）
#   - vars.to.regress 回归掉的变量（比如 "percent_mito" 去除线粒体比例影响）
#   - scale.max 缩放后表达量的最大值（默认 10）
#   - do.scale 是否进行缩放（默认 TRUE）
#   - do.center 是否进行中心化（默认 TRUE）
#   - verbose = TRUE 显示进度信息
# - 内存允许时，建议使用全部基因（features = rownames(sce)），以保留所有基因信息
# - 内存有限时，使用高变异基因（features = VariableFeatures(sce)），以减少计算量
# - 函数默认：features = NULL 时，ScaleData 使用高变异基因
# - 这里选择使用全部基因
# - 结果存储在 sce[["RNA"]]$scale.data 中
message("步骤 2.13：重新数据缩放...")
sce <- ScaleData(sce,
                 features = rownames(sce),  # 使用全部基因（内存允许时推荐）
                 vars.to.regress = NULL,  # 不回归任何变量（可设置为 c("S.Score", "G2M.Score")）
                 scale.max = 10,  # 缩放后表达量最大值
                 do.scale = TRUE,  # 进行缩放
                 do.center = TRUE,  # 进行中心化
                 verbose = TRUE)  # 显示进度信息

# 输出缩放后的 Seurat 对象信息
message("缩放后的 Seurat 对象基本信息：")
print(sce)
# Active assay: RNA (25554 features, 2000 variable features)
# 3 layers present: counts, data, scale.data

# 打印 scale.data 的前 5 个基因和前 5 个细胞
message("scale.data 前 5 个基因和前 5 个细胞：")
print(sce[["RNA"]]$scale.data[1:5, 1:5])

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 2.14：重新降维（PCA）
#-------------------------------------------------------------------------------

# - 使用 PCA 进行降维，基于高变异基因
# - 使用 Seurat 的 RunPCA 函数，常用参数：
#   - features 使用高变异基因（默认 VariableFeatures(sce)）
#   - npcs = 50 选择前 50 个主成分（可调整）
#   - verbose = TRUE 显示进度信息
# - 结果存储在 sce@reductions$pca 中
message("步骤 2.14：重新降维（PCA）...")
sce <- RunPCA(sce,
              features = VariableFeatures(sce),  # 使用高变异基因（默认）
              npcs = 50,  # 计算前 50 个主成分
              verbose = TRUE)  # 显示进度信息

# 输出降维后的 Seurat 对象信息
message("降维后的 Seurat 对象基本信息：")
print(sce)
# An object of class Seurat 
# 25554 features across 4257 samples within 1 assay 
# Active assay: RNA (25554 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 1 dimensional reduction calculated: pca

# 可视化 PCA 结果
# - 使用 DimPlot 绘制 PCA 散点图，展示 PC1 和 PC2 的分布
# - 按细胞周期阶段（Phase）分组，观察细胞周期对 PCA 分布的影响
# - 使用 ElbowPlot 展示主成分的方差解释比例
message("可视化 PCA 结果...")
pca_plot <- DimPlot(sce, 
                    reduction = "pca",  # 使用 PCA 降维结果
                    dims = c(1, 2),  # 绘制 PC1 和 PC2
                    group.by = "Phase",  # 按细胞周期阶段分组
                    label = TRUE,  # 显示分组标签
                    repel = TRUE)  # 避免标签重叠
ggsave(file.path(output_dir, "figures/pca_plot_after_doublet_removal.png"), pca_plot, width = 8, height = 6)

# 绘制 ElbowPlot 展示主成分的方差解释比例
elbow_plot <- ElbowPlot(sce, 
                        reduction = "pca",  # 使用 PCA 降维结果
                        ndims = 50)  # 显示前 50 个主成分
ggsave(file.path(output_dir, "figures/pca_elbow_plot_after_doublet_removal.png"), elbow_plot, width = 8, height = 6)

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 步骤 2.15：保存预处理数据（最终预处理结果）
#-------------------------------------------------------------------------------

# - 保存预处理后的 Seurat 对象为 Rds 文件，路径：processed_data_dir/
# - 文件名格式：dataset_name_processed.rds
# - 包含质控、过滤、标准化、可变基因选择、细胞周期评分、缩放、降维和双细胞去除的结果
message("步骤 2.15：保存预处理数据（最终预处理结果）...")
dir.create(processed_data_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(sce, file = paste0(processed_data_dir, "/", dataset_name, "_pca.rds"))
message("预处理数据已保存至：", paste0(processed_data_dir, "/", dataset_name, "_pca.rds"))

#-------------------------------------------------------------------------------