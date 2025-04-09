
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
DefaultAssay(seu_integrated) <- "RNA"

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
  p <- FeaturePlot(seu_integrated,
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
  p <- VlnPlot(seu_integrated,
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
  p <- DotPlot(seu_integrated,
               features = group_genes,
               group.by = "seurat_clusters",  # 横坐标为簇
               dot.scale = 6) +  # 调整点的大小
    coord_flip() +  # 翻转坐标轴，横坐标为簇，纵坐标为基因
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 旋转横坐标标签（簇）
  ggsave(file.path(marker_viz_dir, paste0("dotplot_top_markers_group_", i, ".png")), p, width = 12, height = 6, dpi = 300)
}

message("标志基因可视化已完成，图表保存至：", marker_viz_dir)

#-------------------------------------------------------------------------------