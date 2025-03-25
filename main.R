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

