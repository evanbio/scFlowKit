# Rutils/create_palette_json.R
#-------------------------------------------------------------------------------

# 颜色方案创建工具：使用 jsonlite 保存自定义配色方案到 JSON 文件
#-------------------------------------------------------------------------------
#
# 背景介绍:
#   - 颜色方案（palette）在数据可视化中起到关键作用，用于区分类别、展示渐变或对比。
#   - JSON 文件格式轻量且跨平台，便于存储和共享自定义颜色方案。
#   - 本函数提供以下功能：
#     - 创建并保存自定义颜色方案为 JSON 文件。
#     - 支持三种颜色类型：
#       - Sequential（连续型）：适用于渐进式数据。
#       - Diverging（发散型）：适用于突出中间值的数据。
#       - Qualitative（分类型）：适用于离散类别区分。
#     - 自动记录创建日志，便于追溯和管理。
#
# 参数说明:
#   - name: 颜色方案名称（字符串，例如 "Blues"）
#   - type: 颜色类型，"sequential"、"diverging" 或 "qualitative"（默认 "sequential"）
#   - colors: 十六进制颜色值向量（支持透明度，例如 "#E64B35B2"）
#   - color_dir: 颜色方案保存的根目录（默认 "colors"）
#   - log_file: 日志文件路径（默认 "colors/palette_creation.log"）
#
# 返回值:
#   - 无（隐式返回 JSON 文件路径），结果保存为 JSON 文件并记录日志
#
# 依赖包:
#   - jsonlite (JSON 文件读写)
#   - cli (命令行交互提示)

create_palette_json <- function(name, 
                                type = c("sequential", "diverging", "qualitative"), 
                                colors, 
                                color_dir = "colors",
                                log_file = file.path(color_dir, "palette_creation.log")) {
  
  # 加载必要的包
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("请先安装 jsonlite 包：install.packages('jsonlite')", call. = FALSE)
  }
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("请先安装 cli 包：install.packages('cli')", call. = FALSE)
  }
  library(jsonlite)
  library(cli)
  
  type <- match.arg(type)
  
  # 检查颜色值
  if (!all(grepl("^#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?$", colors))) {
    stop("颜色值必须是有效的 HEX 码，例如 '#FF5733' 或 '#FF5733B2'")
  }
  
  # 检查 name 的合法性
  if (!is.character(name) || length(name) != 1) {
    stop("颜色方案名称 (name) 必须是长度为 1 的字符串！")
  }
  
  # 检查并创建目录
  palette_dir <- file.path(color_dir, type)
  if (!dir.exists(palette_dir)) {
    dir.create(palette_dir, recursive = TRUE)
    cli_alert_info("已自动创建目录：{.path {palette_dir}}")
  }
  
  # 构建颜色方案列表
  palette_info <- list(
    name = name,
    type = type,
    colors = colors
  )
  
  # JSON 文件路径
  json_file <- file.path(palette_dir, paste0(name, ".json"))
  
  # 检查文件是否存在
  if (file.exists(json_file)) {
    cli_alert_warning("文件已存在：{.path {json_file}}，将被覆盖")
  }
  
  # 保存为 JSON
  tryCatch({
    write_json(palette_info, path = json_file, pretty = TRUE, auto_unbox = TRUE)
    cli_alert_success("已成功创建配色 JSON 文件：{.path {json_file}}")
  }, error = function(e) {
    cli_alert_danger("保存 JSON 文件失败：{e$message}")
    stop(e)
  })
  
  # 记录日志
  log_entry <- paste(Sys.time(), 
                     "| 类型:", type, 
                     "| 名称:", name, 
                     "| 颜色数量:", length(colors), 
                     "| 路径:", json_file, 
                     "\n")
  
  tryCatch({
    cat(log_entry, file = log_file, append = TRUE)
    cli_alert_info("日志已追加到文件：{.path {log_file}}")
  }, error = function(e) {
    cli_alert_danger("日志记录失败：{e$message}")
  })
  
  # 返回文件路径（隐式）
  invisible(json_file)
}

#-------------------------------------------------------------------------------
# 示例用法: 创建并保存颜色方案到 JSON 文件
#-------------------------------------------------------------------------------
# 
# # 示例 1: 创建 Sequential 类型颜色方案 "blues"
# create_palette_json(
#   name = "blues",
#   type = "sequential",
#   colors = c("#deebf7", "#9ecae1", "#3182bd")
# )
# cat("已创建 Sequential 类型颜色方案: colors/sequential/blues.json\n")
# 
# # 示例 2: 创建 Diverging 类型颜色方案 "piyg"（带透明度）
# create_palette_json(
#   name = "piyg",
#   type = "diverging",
#   colors = c("#E64B35B2", "#00A087B2", "#3C5488B2")
# )
# cat("已创建 Diverging 类型颜色方案: colors/diverging/piyg.json\n")
# 
# # 示例 3: 创建 Qualitative 类型颜色方案 "vividset"
# create_palette_json(
#   name = "vividset",
#   type = "qualitative",
#   colors = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")
# )
# cat("已创建 Qualitative 类型颜色方案: colors/qualitative/vividset.json\n")

#-------------------------------------------------------------------------------