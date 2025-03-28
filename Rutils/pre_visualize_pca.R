# Rutils/pre_visualize_pca.R
#-------------------------------------------------------------------------------

# scFlowKit: Visualize PCA Results for Single-Cell RNA-seq Data
#-------------------------------------------------------------------------------

# pre_visualize_pca: 可视化 PCA 结果
# 参数:
#   sce: Seurat 对象，包含 PCA 降维结果
#   output_dir: 输出目录，用于保存 PCA 图
#   reduction: 降维方法，默认 "pca"
#   group.by: 分组变量，字符形式，例如 "sample", "Phase", "percent_mito_binned"
#   split.by: 分面变量，字符形式，例如 "sample", "Phase", "percent_mito_binned"
#   width: 保存图片的宽度，默认 10
#   height: 保存图片的高度，默认 10
pre_visualize_pca <- function(sce, 
                              output_dir, 
                              reduction = "pca", 
                              group.by = "sample", 
                              split.by = "sample", 
                              width = 10, 
                              height = 10) {
  # 验证输入参数是否为 Seurat 对象
  if (!inherits(sce, "Seurat")) {
    stop("参数 'sce' 必须为 Seurat 对象！", call. = FALSE)
  }
  
  # 验证 output_dir 是否为字符类型
  if (!is.character(output_dir)) {
    stop("参数 'output_dir' 必须为字符类型！", call. = FALSE)
  }
  
  # 验证 group.by 是否为字符且存在于元数据中
  if (!is.character(group.by) || length(group.by) != 1 || !group.by %in% colnames(sce@meta.data)) {
    stop("参数 'group.by' 必须为单一字符，且存在于元数据中！", call. = FALSE)
  }
  
  # 验证 split.by 是否为字符且存在于元数据中
  if (!is.character(split.by) || length(split.by) != 1 || !split.by %in% colnames(sce@meta.data)) {
    stop("参数 'split.by' 必须为单一字符，且存在于元数据中！", call. = FALSE)
  }
  
  # 确保输出目录存在
  figures_dir <- file.path(output_dir, "figures")
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 加载必要的包
  library(Seurat)
  library(patchwork)
  
  # 第一行：按 sample 分组
  p1 <- DimPlot(sce,
                reduction = reduction,
                group.by = "sample",
                label = FALSE) +
    labs(title = "PCA by Sample")
  
  # 第一行：按 group.by 分组
  p2 <- DimPlot(sce,
                reduction = reduction,
                group.by = group.by,
                label = FALSE) +
    labs(title = paste0("PCA by ", group.by))
  
  # 第二行：按 group.by 分组，按 split.by 分面
  p3 <- DimPlot(sce,
                reduction = reduction,
                group.by = group.by,
                split.by = split.by,
                label = FALSE) +
    labs(title = paste0("PCA by ", group.by, ", Split by ", split.by))
  
  # 使用 patchwork 组合图表（2 行布局）
  combined_plot <- (p1 | p2) / p3 + plot_layout(heights = c(1, 1))
  
  # 动态生成文件名
  filename <- paste0("preliminary_pca_", group.by, ".png")
  
  # 保存组合图
  ggsave(file.path(figures_dir, filename), 
         plot = combined_plot, 
         width = width, 
         height = height,
         dpi = 300)  
}

#-------------------------------------------------------------------------------