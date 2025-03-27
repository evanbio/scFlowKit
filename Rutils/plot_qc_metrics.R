# Rutils/plot_qc_metrics.R
#-------------------------------------------------------------------------------

# scFlowKit: Visualize QC Metrics for Single-Cell RNA-seq Data
#-------------------------------------------------------------------------------

# 质控指标可视化背景介绍
# - 在单细胞 RNA-seq 分析中，可视化质控指标是评估数据质量的重要步骤，帮助识别低质量细胞和异常模式。
# - 常见的质控指标可视化方法：
#   - 小提琴图（VlnPlot）：展示每个细胞的指标分布（如 UMI 计数、基因数、线粒体基因比例）。
#     - 意义：直观显示指标的分布范围和集中趋势，帮助识别异常值（例如 UMI 计数过低或线粒体比例过高的细胞）。
#   - 散点图（FeatureScatter）：展示两个指标之间的关系（如 UMI 计数 vs 基因数、UMI 计数 vs 线粒体基因比例）。
#     - 意义：揭示指标之间的相关性，例如 UMI 计数和基因数通常呈正相关，线粒体比例过高的细胞可能质量较差。
# - 本函数可视化的图表：
#   - 小提琴图（VlnPlot）：绘制 nCount_RNA（UMI 计数）、nFeature_RNA（基因数）、percent_mito（线粒体基因比例）、log10_ratio_features_to_umi 的分布，按样本分组。
#   - 散点图（FeatureScatter）：绘制 UMI 计数 vs 基因数、UMI 计数 vs 线粒体基因比例、基因数 vs 线粒体基因比例、UMI 计数 vs log10_ratio_features_to_umi，按样本分组。
#   - 综合散点图（ggplot2）：绘制 UMI 计数 vs 基因数（按样本分面），以线粒体基因比例为颜色，添加拟合线和阈值线。
# - 输出：
#   - 所有小提琴图拼成一张图（qc_metrics_combined.png）。
#   - 所有散点图拼成一张图（qc_metrics_scatter_combined.png）。
#   - 综合散点图（qc_metrics_comprehensive.png）。

# plot_qc_metrics: 可视化质控指标
# 参数:
#   sce: Seurat 对象，包含质控指标
#   output_dir: 输出目录，用于保存质控图
#   pt.size: VlnPlot 和 FeatureScatter 中点的显示大小，默认 0.1（显示点）
#   umi_threshold: UMI 计数阈值，默认 500
#   feature_threshold: 基因数阈值，默认 300
#   mito_threshold: 线粒体基因比例阈值，默认 10
#   ratio_threshold: log10_ratio_features_to_umi 阈值，默认 0.8
plot_qc_metrics <- function(sce, output_dir, pt.size = 0.1, 
                            umi_threshold = 500, feature_threshold = 300, 
                            mito_threshold = 10, ratio_threshold = 0.8) {
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

  # 验证阈值参数是否为正数值
  if (!is.numeric(umi_threshold) || umi_threshold <= 0) {
    stop("参数 'umi_threshold' 必须为正数值！", call. = FALSE)
  }
  if (!is.numeric(feature_threshold) || feature_threshold <= 0) {
    stop("参数 'feature_threshold' 必须为正数值！", call. = FALSE)
  }
  if (!is.numeric(mito_threshold) || mito_threshold <= 0) {
    stop("参数 'mito_threshold' 必须为正数值！", call. = FALSE)
  }
  if (!is.numeric(ratio_threshold) || ratio_threshold <= 0) {
    stop("参数 'ratio_threshold' 必须为正数值！", call. = FALSE)
  }

  # 验证元数据是否包含必要的质控指标
  required_metrics <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "log10_ratio_features_to_umi")
  missing_metrics <- setdiff(required_metrics, colnames(sce@meta.data))
  if (length(missing_metrics) > 0) {
    stop("元数据缺少必要的质控指标：", paste(missing_metrics, collapse = ", "), call. = FALSE)
  }

  # 验证元数据是否包含 sample 字段（用于分组）
  if (!"sample" %in% colnames(sce@meta.data)) {
    stop("元数据缺少 'sample' 字段，无法按样本分组！", call. = FALSE)
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
  p1 <- VlnPlot(sce, features = "nCount_RNA", pt.size = pt.size, layer = "counts", group.by = "sample", log = TRUE) +
        labs(title = "UMI Counts per Cell") +
        geom_hline(yintercept = umi_threshold, linetype = "dashed", color = "red")  # 阈值线
  
  # 基因数分布
  p2 <- VlnPlot(sce, features = "nFeature_RNA", pt.size = pt.size, layer = "counts", group.by = "sample") +
        labs(title = "Genes per Cell") +
        geom_hline(yintercept = feature_threshold, linetype = "dashed", color = "red")  # 阈值线
  
  # 线粒体基因比例分布
  p3 <- VlnPlot(sce, features = "percent_mito", pt.size = pt.size, layer = "counts", group.by = "sample") +
        labs(title = "Mitochondrial Gene %") +
        geom_hline(yintercept = mito_threshold, linetype = "dashed", color = "red")  # 阈值线
  
  # log10_ratio_features_to_umi 分布
  p4 <- VlnPlot(sce, features = "log10_ratio_features_to_umi", pt.size = pt.size, layer = "counts", group.by = "sample") +
        labs(title = "Log10 Ratio Features to UMI") +
        geom_hline(yintercept = ratio_threshold, linetype = "dashed", color = "red")  # 阈值线
  
  # 使用 patchwork 将四张 VlnPlot 图拼成一张（4 行）
  combined_vln_plot <- p1 / p2 / p3 / p4 + plot_layout(ncol = 1)
  
  # 保存 VlnPlot 拼图
  ggsave(file.path(figures_dir, "qc_metrics_combined.png"), combined_vln_plot, width = 15, height = 20)

  # 使用 Seurat 的 FeatureScatter 绘制质控指标散点图
  # UMI 计数 vs 基因数
  s1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample") +
        labs(title = "UMI Counts vs Genes") +
        geom_vline(xintercept = umi_threshold, linetype = "dashed", color = "red") +
        geom_hline(yintercept = feature_threshold, linetype = "dashed", color = "red") +
        theme(plot.title = element_text(size = 14),  # 增大标题字体
              axis.title = element_text(size = 12),  # 增大轴标签字体
              axis.text = element_text(size = 10))    # 增大轴刻度字体
  
  # UMI 计数 vs 线粒体基因比例
  s2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent_mito", group.by = "sample") +
        labs(title = "UMI Counts vs Mito %") +
        geom_vline(xintercept = umi_threshold, linetype = "dashed", color = "red") +
        geom_hline(yintercept = mito_threshold, linetype = "dashed", color = "red") +
        theme(plot.title = element_text(size = 14),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 10))
  
  # 基因数 vs 线粒体基因比例
  s3 <- FeatureScatter(sce, feature1 = "nFeature_RNA", feature2 = "percent_mito", group.by = "sample") +
        labs(title = "Genes vs Mito %") +
        geom_vline(xintercept = feature_threshold, linetype = "dashed", color = "red") +
        geom_hline(yintercept = mito_threshold, linetype = "dashed", color = "red") +
        theme(plot.title = element_text(size = 14),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 10))
  
  # UMI 计数 vs log10_ratio_features_to_umi
  s4 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "log10_ratio_features_to_umi", group.by = "sample") +
        labs(title = "UMI Counts vs Log10 Ratio Features to UMI") +
        geom_vline(xintercept = umi_threshold, linetype = "dashed", color = "red") +
        geom_hline(yintercept = ratio_threshold, linetype = "dashed", color = "red") +
        theme(plot.title = element_text(size = 14),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 10))
  
  # 使用 patchwork 将四张 FeatureScatter 图拼成一张（4 行）
  combined_scatter_plot <- s1 / s2 / s3 / s4 + plot_layout(ncol = 1)
  
  # 保存 FeatureScatter 拼图
  ggsave(file.path(figures_dir, "qc_metrics_scatter_combined.png"), combined_scatter_plot, width = 10, height = 20)

  # 使用 ggplot2 绘制综合散点图（参考 HBC）
  # - 绘制 nCount_RNA vs nFeature_RNA，按样本分面，以 percent_mito 为颜色
  # - 添加拟合线和阈值线
  metadata <- sce@meta.data
  p_comprehensive <- metadata %>% 
    ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = percent_mito)) + 
    geom_point(size = 1, alpha = 0.8) +  # 优化点大小和透明度
    scale_colour_gradient(low = "gray90", high = "#8856a7") +  # 调整颜色梯度
    stat_smooth(method = "lm", color = "darkblue", linewidth = 1, linetype = "solid") +  # 调整拟合线
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    labs(title = "UMI Counts vs Genes by Sample", 
         x = "UMI Counts (log10)", 
         y = "Genes Detected (log10)") +
    geom_vline(xintercept = umi_threshold, linetype = "dashed", color = "red") +
    geom_hline(yintercept = feature_threshold, linetype = "dashed", color = "red") +
    facet_wrap(~sample, labeller = label_both) +  # 优化分面标签
    theme(strip.text = element_text(size = 8))  # 调整分面标签字体大小
  
  # 保存综合散点图
  ggsave(file.path(figures_dir, "qc_metrics_comprehensive.png"), p_comprehensive, width = 15, height = 10)

  # 提示用户可视化完成
  message("质控指标可视化完成！")
}

#-------------------------------------------------------------------------------