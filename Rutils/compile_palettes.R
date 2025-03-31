# Rutils/compile_palettes.R
#-------------------------------------------------------------------------------

# 颜色方案编译工具：从 JSON 文件编译颜色方案并保存为 RDS 文件
#-------------------------------------------------------------------------------
#
# 背景介绍:
#   - 颜色方案（palette）在数据可视化中用于区分类别、显示渐变或突出对比。
#   - 本工具从指定目录下的 JSON 文件中读取颜色方案，编译为统一的结构，并保存为 RDS 文件。
#   - 支持三种颜色类型：
#     - Sequential（连续型）：适合渐进式数据。
#     - Diverging（发散型）：适合突出中间值的数据。
#     - Qualitative（分类型）：适合离散类别区分。
#   - 编译过程会记录详细日志，包括成功、警告和错误信息。
#
# 参数说明:
#   - color_dir: 颜色方案 JSON 文件的根目录（默认 "colors"）
#   - output_rds: 编译后的 RDS 文件保存路径（默认 "colors/color_palettes.rds"）
#   - log_file: 日志文件保存路径（默认 "colors/compile_palettes.log"）
#
# 返回值:
#   - 无（隐式返回 output_rds 文件路径），结果保存为 RDS 文件
#
# 依赖包:
#   - jsonlite (JSON 文件读写)
#   - cli (命令行交互提示)

compile_palettes <- function(color_dir = "colors",
                             output_rds = file.path(color_dir, "color_palettes.rds"),
                             log_file = file.path(color_dir, "compile_palettes.log")) {
  
  # 加载必要的包
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("请先安装 jsonlite 包：install.packages('jsonlite')", call. = FALSE)
  }
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("请先安装 cli 包：install.packages('cli')", call. = FALSE)
  }
  library(jsonlite)
  library(cli)
  
  cli_h1("开始编译颜色方案 (JSON → RDS)")
  
  # 只读取特定子文件夹的 JSON 文件
  json_files <- unlist(lapply(c("sequential", "diverging", "qualitative"), function(subdir) {
    list.files(file.path(color_dir, subdir), pattern = "\\.json$", full.names = TRUE)
  }))
  
  if (length(json_files) == 0) {
    cli_alert_warning("未找到任何 JSON 文件，请检查路径：{.path {color_dir}}")
    return(invisible(NULL))
  }
  
  # 初始化颜色集合
  palettes <- list(
    sequential = list(),
    diverging = list(),
    qualitative = list()
  )
  
  # 编译过程日志（初始化）
  compile_log <- c(sprintf("\n=== [%s] 编译开始 ===", Sys.time()))
  
  for (json_file in json_files) {
    palette_info <- fromJSON(json_file)
    
    required_fields <- c("name", "type", "colors")
    if (!all(required_fields %in% names(palette_info))) {
      msg <- sprintf("JSON 缺少必要字段，跳过：%s", json_file)
      cli_alert_warning(msg)
      compile_log <- c(compile_log, sprintf("[%s] 警告: %s", Sys.time(), msg))
      next
    }
    
    type <- palette_info$type
    name <- palette_info$name
    colors <- palette_info$colors
    
    if (!type %in% names(palettes)) {
      msg <- sprintf("未知颜色类型 '%s'，跳过：%s", type, json_file)
      cli_alert_warning(msg)
      compile_log <- c(compile_log, sprintf("[%s] 警告: %s", Sys.time(), msg))
      next
    }
    
    palettes[[type]][[name]] <- colors
    
    msg <- sprintf("成功合并方案 '%s' (类型: %s, 颜色数: %d)", name, type, length(colors))
    cli_alert_success(msg)
    compile_log <- c(compile_log, sprintf("[%s] 成功: %s", Sys.time(), msg))
  }
  
  # 保存 RDS 文件
  tryCatch({
    saveRDS(palettes, file = output_rds)
    msg <- sprintf("已保存所有颜色方案到 RDS 文件：%s", output_rds)
    cli_alert_success(msg)
    compile_log <- c(compile_log, sprintf("[%s] 完成: %s", Sys.time(), msg))
  }, error = function(e) {
    msg <- sprintf("保存 RDS 失败：%s", e$message)
    cli_alert_danger(msg)
    compile_log <- c(compile_log, sprintf("[%s] 错误: %s", Sys.time(), msg))
  })
  
  # 追加编译日志
  cat(paste(compile_log, collapse = "\n"), file = log_file, append = TRUE, sep = "\n")
  cli_alert_info("编译日志已追加到：{.path {log_file}}")
  
  invisible(output_rds)
}

#-------------------------------------------------------------------------------
# 示例用法: 编译颜色方案并保存 RDS 文件
#-------------------------------------------------------------------------------

# 前提：假设已通过 create_palette_json() 创建了以下颜色方案
# colors/qualitative/vividset.json
# colors/qualitative/softtrio.json
# colors/qualitative/harmonysix.json

# # 示例 1: 使用默认路径编译颜色方案
# compile_palettes()
# cat("颜色方案已编译并保存至: colors/color_palettes.rds\n")
# 
# # 示例 2: 指定自定义输出路径
# compile_palettes(
#   color_dir = "colors",
#   output_rds = "data/custom_palettes.rds",
#   log_file = "data/compile_log.txt"
# )
# cat("颜色方案已编译并保存至: data/custom_palettes.rds\n")
# 
# # 示例 3: 检查编译结果
# palettes <- readRDS("colors/color_palettes.rds")
# print("已编译的颜色方案数量:")
# print(length(unlist(palettes, recursive = FALSE)))
# print("示例颜色方案 (vividset):")
# print(palettes$qualitative$vividset)

#-------------------------------------------------------------------------------