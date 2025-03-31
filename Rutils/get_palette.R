# Rutils/get_palette.R
#-------------------------------------------------------------------------------

# 颜色方案提取工具：从 RDS 文件中提取指定颜色方案
#-------------------------------------------------------------------------------
#
# 背景介绍:
#   - 本函数用于从编译好的 RDS 文件中提取颜色方案，提供灵活的颜色数量控制。
#   - RDS 文件由 compile_palettes() 生成，包含 sequential、diverging 和 qualitative 三类颜色方案。
#   - 支持以下功能：
#     - 根据名称和类型提取颜色方案。
#     - 可指定返回的颜色数量（n）。
#     - 当类型错误时，自动查找并提示正确的类型。
#     - 返回 HEX 颜色值向量，直接用于可视化或其他操作。
#
# 参数说明:
#   - name: 颜色方案名称（字符串，例如 "vividset"）
#   - type: 颜色类型，"sequential"、"diverging" 或 "qualitative"（默认 "sequential"）
#   - n: 返回的颜色数量（正整数，默认 NULL 表示返回全部颜色）
#   - palette_rds: RDS 文件路径（默认 "colors/color_palettes.rds"）
#
# 返回值:
#   - 字符向量，包含指定颜色方案的 HEX 颜色值
#
# 依赖包:
#   - cli (命令行交互提示)

get_palette <- function(name, 
                        type = c("sequential", "diverging", "qualitative"),
                        n = NULL,
                        palette_rds = "colors/color_palettes.rds") {
  
  # 加载必要的包
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("请先安装 cli 包：install.packages('cli')", call. = FALSE)
  }
  library(cli)
  
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
    # 遍历所有类型查找 name
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
  
  cli_alert_success("提取 '{name}' 成功，颜色数：{.val {max_len}}")
  
  # 如果未指定 n，返回全部颜色
  if (is.null(n)) {
    return(colors)
  }
  
  # 检查 n 是否为正整数
  if (!is.numeric(n) || n != round(n) || n <= 0) {
    cli_alert_danger("参数 'n' 必须是正整数，当前值为：{.val {n}}")
    stop("n 必须是正整数")
  }
  
  # 检查请求数量是否超过最大长度
  if (n > max_len) {
    stop(sprintf("请求的颜色数量 (%d) 超过方案 '%s' 的最大长度 (%d)！", 
                 n, name, max_len))
  }
  
  # 返回指定数量的颜色
  return(colors[seq_len(n)])
}

#-------------------------------------------------------------------------------
# 示例用法: 从 RDS 文件提取颜色方案
#-------------------------------------------------------------------------------

# 前提：假设已通过 compile_palettes() 生成了 colors/color_palettes.rds
# 包含以下颜色方案：
# - qualitative/vividset (9 个颜色)
# - qualitative/softtrio (3 个颜色)
# - qualitative/harmonysix (6 个颜色)

# # 示例 1: 提取全部颜色 (vividset)
# colors_vividset <- get_palette("vividset", type = "qualitative")
# cat("提取 'vividset' 全部颜色 (", length(colors_vividset), "个):", 
#     paste(colors_vividset, collapse = ", "), "\n")
# 
# # 示例 2: 提取指定数量颜色 (softtrio，前 2 个)
# colors_softtrio <- get_palette("softtrio", type = "qualitative", n = 2)
# cat("提取 'softtrio' 前 2 个颜色:", paste(colors_softtrio, collapse = ", "), "\n")
# 
# # 示例 3: 提取指定数量颜色 (harmonysix，6 个颜色)
# colors_harmonysix <- get_palette("harmonysix", type = "qualitative", n = 6)
# cat("提取 'harmonysix' 6 个颜色:", paste(colors_harmonysix, collapse = ", "), "\n")
# 
# # 示例 4: 测试类型错误 (vividset 在 qualitative，但指定 sequential)
# # 这会报错并提示正确类型
# get_palette("vividset", type = "sequential")
# 
# # 示例 5: 测试非正整数 n (softtrio，n = 0)
# # 这会报错并提示 n 必须是正整数
# get_palette("softtrio", type = "qualitative", n = 0)
# 
# # 示例 6: 测试非整数 n (softtrio，n = 1.5)
# # 这会报错并提示 n 必须是正整数
# get_palette("softtrio", type = "qualitative", n = 1.5)
# 
# # 示例 7: 用于可视化 (vividset 前 5 个颜色)
# colors_vividset_5 <- get_palette("vividset", type = "qualitative", n = 5)
# barplot(rep(1, 5), col = colors_vividset_5, main = "vividset (5 colors)")

#-------------------------------------------------------------------------------