# Rutils/preview_palette.R
#-------------------------------------------------------------------------------

# 颜色方案预览工具：从 RDS 文件提取颜色并展示为图形
#-------------------------------------------------------------------------------
#
# 背景介绍:
#   - 本函数用于预览颜色方案的视觉效果，帮助用户选择合适的配色。
#   - 直接从 RDS 文件中读取颜色方案，支持多种图形展示。
#   - 支持的图形类型：
#     - bar：条形图，展示颜色顺序和对比。
#     - pie：饼图，展示颜色分布。
#     - point：点图，模拟数据点的颜色效果。
#     - rect：矩形图，平铺展示颜色块。
#     - circle：圆形图，圆点排列展示颜色。
#
# 参数说明:
#   - name: 颜色方案名称（字符串，例如 "vividset"）
#   - type: 颜色类型，"sequential"、"diverging" 或 "qualitative"（默认 "sequential"）
#   - n: 返回的颜色数量（正整数，默认 NULL 表示使用全部颜色）
#   - plot_type: 图形类型，"bar"、"pie"、"point"、"rect" 或 "circle"（默认 "bar"）
#   - title: 图形标题（默认使用 name）
#   - palette_rds: RDS 文件路径（默认 "colors/color_palettes.rds"）
#
# 返回值:
#   - 无（直接生成图形展示颜色方案）
#
# 依赖包:
#   - cli (命令行交互提示)

preview_palette <- function(name, 
                            type = c("sequential", "diverging", "qualitative"),
                            n = NULL,
                            plot_type = c("bar", "pie", "point", "rect", "circle"),
                            title = name,
                            palette_rds = "colors/color_palettes.rds") {
  
  # 加载必要的包
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("请先安装 cli 包：install.packages('cli')", call. = FALSE)
  }
  library(cli)
  
  # 检查 plot_type
  plot_type <- match.arg(plot_type)
  
  # 检查 RDS 文件是否存在
  if (!file.exists(palette_rds)) {
    cli_alert_danger("颜色方案文件不存在：{.path {palette_rds}}")
    stop("请检查路径")
  }
  
  # 读取 RDS 文件
  palettes <- readRDS(palette_rds)
  
  # 检查 type 是否合法
  valid_types <- names(palettes)
  type <- match.arg(type)
  if (!type %in% valid_types) {
    cli_alert_danger("无效的类型 '{type}'，可选类型为：{.val {valid_types}}")
    stop("类型错误")
  }
  
  # 检查 name 的合法性
  if (!is.character(name) || length(name) != 1) {
    stop("颜色方案名称 (name) 必须是长度为 1 的字符串！")
  }
  
  # 检查指定 type 下是否存在 name
  if (!name %in% names(palettes[[type]])) {
    found_type <- NULL
    for (t in valid_types) {
      if (name %in% names(palettes[[t]])) {
        found_type <- t
        break
      }
    }
    if (is.null(found_type)) {
      cli_alert_danger("未找到颜色方案 '{name}' 在任何类型中")
      stop("颜色方案不存在")
    } else {
      cli_alert_warning("未在类型 '{type}' 中找到 '{name}'，但在类型 '{found_type}' 中找到")
      cli_alert_info("建议更改 type 参数为 '{found_type}'")
      stop("类型与名称不匹配")
    }
  }
  
  # 提取颜色
  colors <- palettes[[type]][[name]]
  max_len <- length(colors)
  
  # 如果指定 n，调整颜色数量
  if (!is.null(n)) {
    if (!is.numeric(n) || n != round(n) || n <= 0) {
      cli_alert_danger("参数 'n' 必须是正整数，当前值为：{.val {n}}")
      stop("n 必须是正整数")
    }
    if (n > max_len) {
      cli_alert_danger("请求的颜色数量 ({n}) 超过方案 '{name}' 的最大长度 ({max_len})")
      stop("颜色数量超出范围")
    }
    colors <- colors[seq_len(n)]
  }
  
  num_colors <- length(colors)
  cli_alert_success("提取 '{name}' 成功，颜色数：{.val {num_colors}}")
  
  # 根据 plot_type 绘制图形并显示 HEX 码
  switch(plot_type,
         "bar" = {
           barplot(rep(1, num_colors), 
                   col = colors, 
                   border = NA, 
                   space = 0, 
                   axes = FALSE, 
                   main = title, 
                   names.arg = colors, 
                   las = 2,  # 垂直显示 HEX 码
                   cex.names = 0.8)  # 调整字体大小
         },
         "pie" = {
           pie(rep(1, num_colors), 
               col = colors, 
               labels = colors,  # 显示 HEX 码
               border = "white", 
               main = title, 
               cex = 0.8)  # 调整字体大小
         },
         "point" = {
           plot(seq_len(num_colors), 
                rep(1, num_colors), 
                pch = 19, 
                cex = 5, 
                col = colors, 
                axes = FALSE, 
                xlab = "", 
                ylab = "", 
                main = title)
           text(seq_len(num_colors), 
                rep(1.2, num_colors), 
                labels = colors, 
                pos = 3,  # HEX 码在上方
                cex = 0.8)
         },
         "rect" = {
           plot(0, 0, type = "n", 
                xlim = c(0, num_colors), 
                ylim = c(0, 1), 
                axes = FALSE, 
                xlab = "", 
                ylab = "", 
                main = title)
           rect(0:(num_colors-1), 0, 1:num_colors, 1, 
                col = colors, 
                border = NA)
           text((0:(num_colors-1) + 1:num_colors) / 2, 0.5, 
                labels = colors, 
                col = "white", 
                cex = 0.8)
         },
         "circle" = {
           plot(0, 0, type = "n", 
                xlim = c(0, num_colors), 
                ylim = c(0, 1), 
                axes = FALSE, 
                xlab = "", 
                ylab = "", 
                main = title)
           symbols(seq_len(num_colors) - 0.5, rep(0.5, num_colors), 
                   circles = rep(0.4, num_colors), 
                   inches = FALSE, 
                   bg = colors, 
                   add = TRUE)
           text(seq_len(num_colors) - 0.5, 0.5, 
                labels = colors, 
                col = "white", 
                cex = 0.8)
         },
         stop("不支持的图形类型，目前仅支持：'bar', 'pie', 'point', 'rect', 'circle'")
  )
  
  cli_alert_info("预览 '{name}' 完成，展示类型：{.val {plot_type}}，颜色数：{.val {num_colors}}")
}

#-------------------------------------------------------------------------------
# 示例用法: 预览颜色方案
#-------------------------------------------------------------------------------

# 前提：假设已通过 compile_palettes() 生成了 colors/color_palettes.rds
# 包含以下颜色方案：
# - qualitative/vividset (9 个颜色)
# - qualitative/softtrio (3 个颜色)
# - qualitative/harmonysix (6 个颜色)

# # 示例 1: 用条形图预览 vividset（全部颜色）
# preview_palette("vividset", type = "qualitative", plot_type = "bar")
# 
# # 示例 2: 用饼图预览 softtrio（全部颜色）
# preview_palette("softtrio", type = "qualitative", plot_type = "pie")
# 
# # 示例 3: 用点图预览 harmonysix（前 4 个颜色）
# preview_palette("harmonysix", type = "qualitative", n = 4, plot_type = "point")
# 
# # 示例 4: 用矩形图预览 vividset（前 5 个颜色）
# preview_palette("vividset", type = "qualitative", n = 5, plot_type = "rect", 
#                 title = "VividSet Preview (5 Colors)")
# 
# # 示例 5: 用圆形图预览 softtrio（全部颜色）
# preview_palette("softtrio", type = "qualitative", plot_type = "circle")
# 
# # 示例 6: 测试不支持的图形类型
# preview_palette("softtrio", type = "qualitative", plot_type = "line")
# 
# # 示例 7: 测试非正整数 n (softtrio，n = 0)
# preview_palette("softtrio", type = "qualitative", n = 0, plot_type = "bar")

#-------------------------------------------------------------------------------