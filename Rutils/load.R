#-------------------------------------------------------------------------------

# scFlowKit: Load Script for Dependencies
#-------------------------------------------------------------------------------

# 载入依赖包
#-------------------------------------------------------------------------------

library(Seurat)        # 加载 Seurat 包，用于单细胞 RNA-seq 数据分析
library(patchwork)     # 加载 patchwork 包，用于拼图（例如组合多个 ggplot 图）
library(ggplot2)       # 加载 ggplot2 包，用于绘图（例如 UMAP、t-SNE 可视化）
library(rlang)         # 加载 rlang 包，用于解析过滤条件表达式（例如动态条件）
library(stringr)       # 加载 stringr 包，字符串处理（例如首字母大写、字符串操作）
library(DoubletFinder) # 加载 DoubletFinder 包，用于双细胞检测
library(MAST)          # 加载 MAST 包，用于差异表达分析
library(tidyverse)     # 加载 tidyverse 包，包含 dplyr, tidyr 等，用于数据处理和可视化
library(glmGamPoi)     # 加载 glmGamPoi 包，用于加速 SCTransform 的计算
library(harmony)       # 加载 harmony 包，用于批次效应校正（Harmony 整合方法）

# 提示用户载入完成
cat("scFlowKit 依赖包载入完成！\n")

#-------------------------------------------------------------------------------
