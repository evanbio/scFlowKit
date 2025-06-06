# Rutils/plot_marker.R
#-------------------------------------------------------------------------------
# scFlowKit: 可视化 Marker 基因（FeaturePlot / VlnPlot / DotPlot）
#-------------------------------------------------------------------------------
#
# plot_marker: 可视化每个 cluster 的 marker 基因表达情况。
# - 输入为 top_marker_df，必须包含 gene 和 cluster 两列。
# - 输出 FeaturePlot、VlnPlot、DotPlot 三类图像，保存至指定目录。
#
# 参数：
#   seu：Seurat 对象。
#   top_markers_df：marker 基因表，需包含 gene 和 cluster 两列。
#   outdir：输出目录（文件夹路径）。
#   group.by：分组变量，用于 VlnPlot 和 DotPlot，默认 "seurat_clusters"。
# 返回：
#   无（直接保存图像文件）。
#-------------------------------------------------------------------------------

plot_marker <- function(
  seu,
  top_markers_df,
  outdir = "results/figures/markers",
  group.by = "seurat_clusters"
) {
  #-------------------- 参数检查 --------------------

  if (!inherits(seu, "Seurat")) {
    stop("seu 必须为 Seurat 对象！", call. = FALSE)
  }
  if (!all(c("gene", "cluster") %in% colnames(top_markers_df))) {
    stop("top_markers_df 必须包含 gene 和 cluster 列！", call. = FALSE)
  }
  if (!group.by %in% colnames(seu@meta.data)) {
    stop("group.by 字段在 Seurat 对象 meta.data 中不存在：", group.by, call. = FALSE)
  }

  #-------------------- 准备绘图 --------------------
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  cli::cli_alert_info("绘图总数：{nrow(top_markers_df)} 个（每个 cluster topN 个基因）")

  #-------------------- FeaturePlot --------------------
  # 使用 FeaturePlot 可视化 top 基因在 UMAP 空间的表达
  cli::cli_h3("绘制 FeaturePlot...")
  for (i in seq_len(nrow(top_markers_df))) {
    gene <- top_markers_df$gene[i]
    cluster <- top_markers_df$cluster[i]

    p <- Seurat::FeaturePlot(
      object = seu,
      features = gene,
      reduction = "umap",
      pt.size = 0.5,
      alpha = 0.7,
      label = TRUE,  # 显示基因名
      repel = TRUE   # 避免标签重叠
    )

    filename <- paste0("feature_", gene, "_cluster", cluster, ".png")
    ggsave(file.path(outdir, filename), p, width = 8, height = 6, dpi = 300)
  }

  #-------------------- VlnPlot --------------------
  # 使用 VlnPlot 可视化 top 基因在不同聚类中的表达分布
  cli::cli_h3("绘制 VlnPlot...")
  for (i in seq_len(nrow(top_markers_df))) {
    gene <- top_markers_df$gene[i]
    cluster <- top_markers_df$cluster[i]

    p <- Seurat::VlnPlot(
      object = seu,
      features = gene,
      group.by = group.by, # 按聚类分组
      pt.size = 0
    )
    filename <- paste0("vln_", gene, "_cluster", cluster, ".png")
    ggsave(file.path(outdir, filename), p, width = 10, height = 6, dpi = 300)
  }
  
  #-------------------- DotPlot --------------------
  # 使用 DotPlot 可视化 top 基因在不同聚类中的表达比例和表达量
  cli::cli_h3("绘制 DotPlot...")
  dot_groups <- split(top_markers_df, top_markers_df$cluster)

  for (cluster in names(dot_groups)) {
    genes <- dot_groups[[cluster]]$gene

    p <- Seurat::DotPlot(
      object = seu,
      features = genes,
      group.by = group.by,  # 横坐标为簇
      dot.scale = 6
    ) +
      coord_flip() +  # 翻转坐标轴，横坐标为簇，纵坐标为基因
      ggtitle(paste0("DotPlot - Cluster ", cluster)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    filename <- paste0("dotplot_cluster", cluster, ".png")
    ggsave(file.path(outdir, filename), p, width = 12, height = 6, dpi = 300)
  }

  cli::cli_alert_success("全部图像保存至：{outdir}")
}
