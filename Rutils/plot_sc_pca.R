# Rutils/plot_sc_pca.R
#-------------------------------------------------------------------------------
# scFlowKit: Visualize PCA Results for Single-Cell RNA-seq Data
#-------------------------------------------------------------------------------
# 
# plot_sc_pca: 可视化单细胞 PCA 降维结果
# 参数:
#   sce: Seurat 对象，包含 PCA 降维结果
#   output_dir: 输出目录，用于保存 PCA 图
#   reduction: 降维方法，默认 "pca"
#   dims: DimPlot 使用的 PCA 维度，默认 c(1, 2)，即 PC1 和 PC2
#   group.by: 分组变量，字符形式，例如 "sample", "Phase", "percent_mito_binned"
#   split.by: 分面变量，字符形式，例如 "sample", "Phase", "percent_mito_binned"
#   ndims: ElbowPlot 显示的主成分数量，默认 50
#   width: 保存图片的宽度，默认 10
#   height: 保存图片的高度，默认 10
#   dpi: 图形分辨率，默认 300
#
# 返回:
#   无返回值，直接保存 DimPlot、ElbowPlot 和 DimHeatmap到指定目录
plot_sc_pca <- function(sce,
                        output_dir,
                        reduction = "pca",
                        dims = c(1, 2),
                        group.by = "sample",
                        split.by = "sample",
                        ndims = 50,
                        prefix = NULL,
                        plot_elbow = FALSE,
                        plot_heatmap = FALSE,
                        width = 10,
                        height = 10,
                        dpi = 300) {

  # ------------------------- 参数检查 -------------------------                          
  # 验证输入参数是否为 Seurat 对象
  if (!inherits(sce, "Seurat")) {
    stop("参数 'sce' 必须为 Seurat 对象！", call. = FALSE)
  }

  # 验证 output_dir 是否为字符类型
  if (!is.character(output_dir)) {
    stop("参数 'output_dir' 必须为字符类型！", call. = FALSE)
  }

  # 验证 reduction 是否存在
  if (!reduction %in% names(sce@reductions)) {
    stop("参数 'reduction' 必须为 Seurat 对象中的一个降维结果！", call. = FALSE)
  }

  # 验证 dims 是否为数值向量
  if (!is.numeric(dims) || length(dims) != 2) {
    stop("参数 'dims' 必须为长度为 2 的数值向量！", call. = FALSE)
  }

  # 验证 group.by 是否为字符且存在于元数据中
  if (!is.character(group.by) || length(group.by) != 1 || !group.by %in% colnames(sce@meta.data)) {
    stop("参数 'group.by' 必须为单一字符，且存在于元数据中！", call. = FALSE)
  }

  # 验证 split.by 是否为字符且存在于元数据中
  if (!is.character(split.by) || length(split.by) != 1 || !split.by %in% colnames(sce@meta.data)) {
    stop("参数 'split.by' 必须为单一字符，且存在于元数据中！", call. = FALSE)
  }

  # 验证 ndims 是否为正整数
  if (!is.numeric(ndims) || ndims <= 0) {
    stop("参数 'ndims' 必须为正整数！", call. = FALSE)
  }

  if (!is.null(prefix) && (!is.character(prefix) || length(prefix) != 1)) {
    stop("参数 'prefix' 必须为长度为 1 的字符向量，或设为 NULL！", call. = FALSE)
  }

  if (!is.logical(plot_elbow) || length(plot_elbow) != 1) {
    stop("参数 'plot_elbow' 必须为单一逻辑值（TRUE 或 FALSE）！", call. = FALSE)
  }

  if (!is.logical(plot_heatmap) || length(plot_heatmap) != 1) {
    stop("参数 'plot_heatmap' 必须为单一逻辑值（TRUE 或 FALSE）！", call. = FALSE)
  }
  
  # ------------------------- 准备环境 -------------------------
  cli::cli_h2("🎯 开始绘制 PCA 可视化图")

  suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(patchwork)
  })

  if (is.null(prefix)) {
    prefix <- group.by
  }

  # 确保输出目录存在
  figures_dir <- file.path(output_dir, "figures")
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  # ------------------------- DimPlot 图 -------------------------

  # 第一行左：按 sample 分组
  cli::cli_text("绘制 DimPlot（按 sample 分组）")
  p1 <- DimPlot(sce,
                reduction = reduction,
                dims = dims,
                group.by = "sample",
                label = FALSE,
                pt.size = 0.5) +
    labs(title = "PCA by Sample") 

  # 第一行右：按 group.by 分组
  cli::cli_text("绘制 DimPlot（按 {group.by} 分组）", .envir = environment())
  p2 <- DimPlot(sce,
                reduction = reduction,
                dims = dims,
                group.by = group.by,
                label = FALSE,
                pt.size = 0.5) +
    labs(title = paste0("PCA by ", group.by)) 

  # 第二行：按 group.by 分组，按 split.by 分面
  cli::cli_text("绘制 DimPlot（按 {group.by} 分组，按 {split.by} 分面）...", .envir = environment())
  p3 <- DimPlot(sce,
                reduction = reduction,
                dims = dims,
                group.by = group.by,
                split.by = split.by,
                label = FALSE,
                pt.size = 0.5) +
    labs(title = paste0("PCA by ", group.by, ", Split by ", split.by)) 

  # 使用 patchwork 组合图表（2 行布局）
  combined_plot <- (p1 | p2) / p3 + plot_layout(heights = c(1, 1))

  # 动态生成文件名
  pca_file <- file.path(figures_dir, paste0(prefix, "_pca_dimplot.png"))

  # 保存组合图
  ggsave(file.path(figures_dir, filename),
         plot = combined_plot,
         width = width,
         height = height,
         dpi = dpi)
  cli::cli_alert_success("✅ DimPlot 已保存：{pca_file}")

  # ------------------------- ElbowPlot 图 -------------------------
  if (plot_elbow) {
    cli::cli_text("📐 绘制 ElbowPlot（主成分解释度）")
    p_elbow <- ElbowPlot(seu, reduction = reduction, ndims = ndims) +
              labs(title = "Elbow Plot of PCA")

    # 保存 ElbowPlot
    elbow_file <- file.path(figures_dir, paste0(prefix, "_pca_elbowplot.png"))
    ggsave(elbow_file, 
          plot = p_elbow, 
          width = width, 
          height = height, 
          dpi = dpi)
    cli::cli_alert_success("✅ ElbowPlot 已保存：{elbow_file}")
  }

  # ------------------------- DimHeatmap 图 -------------------------
  # 绘制 DimHeatmap（前 9 个 PCs）
  if (plot_heatmap) {
    cli::cli_text("🔥 绘制 DimHeatmap（前 9 个主成分）")
    p_heatmap <- DimHeatmap(sce,
                            dims = 1:9,
                            cells = 500,
                            balanced = TRUE)

    # 保存 DimHeatmap
    heatmap_file <- file.path(figures_dir, paste0(prefix, "_pca_dimheatmap.png"))
    ggsave(heatmap_file,
          plot = p_heatmap,
          width = 12,
          height = 10,
          dpi = dpi)
    cli::cli_alert_success("✅ DimHeatmap 已保存：{heatmap_file}")
  }

  cli::cli_h2("🎉 PCA 可视化完成！")
}

#-------------------------------------------------------------------------------