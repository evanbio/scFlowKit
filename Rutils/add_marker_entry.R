# Rutils/add_marker_entry.R
#-------------------------------------------------------------------------------
# 向现有 marker set 添加新条目并更新原文件
#-------------------------------------------------------------------------------
#
# 背景说明：
# 本函数用于向已有的 marker set 添加新条目，并通过调用 create_marker_set 更新原文件。
# 版本号由用户通过 version 参数指定，文件名保持不变。
#
# 参数说明：
# - set_name   : marker set 的基础名称（单一字符值，不带 .json 后缀）
# - cell_type  : 新增的细胞类型名称（单一字符值）
# - marker_genes : 新增的 marker genes（字符向量或逗号分隔的字符串）
# - version    : 新版本号（单一字符值，例如 "v2"，由用户指定）
# - source     : 数据来源（单一字符值）
# - log_file   : 日志保存路径，默认 "logs/marker_set/add_marker_entry.log"
#
# 返回值：
# - 无返回值，直接更新 JSON 文件，并记录日志
#-------------------------------------------------------------------------------

add_marker_entry <- function(
  set_name,
  cell_type,
  marker_genes,
  version,
  source,
  log_file = "logs/marker_set/add_marker_entry.log"
) {
  #-------------------- 依赖检查 --------------------
  for (pkg in c("tibble", "stringr", "jsonlite", "cli", "dplyr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("请先安装 ", pkg, " 包！", call. = FALSE)
    }
  }

  #-------------------- 参数校验 --------------------
  if (!is.character(set_name) || length(set_name) != 1) {
    stop("参数 'set_name' 必须为单一字符值！", call. = FALSE)
  }
  if (grepl("\\.json$", set_name)) {
    cli::cli_alert_warning("参数 'set_name' 不需包含 .json 后缀，已自动移除")
    set_name <- sub("\\.json$", "", set_name)
  }
  if (!grepl("^[a-zA-Z0-9_]+$", set_name)) {
    stop("参数 'set_name' 只能包含字母、数字和下划线！", call. = FALSE)
  }
  if (!is.character(cell_type) || length(cell_type) != 1 || cell_type == "") {
    stop("参数 'cell_type' 必须为单一非空字符值！", call. = FALSE)
  }
  if (!is.character(marker_genes) || length(marker_genes) == 0) {
    stop("参数 'marker_genes' 必须为非空字符向量！", call. = FALSE)
  }
  if (!is.character(version) || length(version) != 1) {
    stop("参数 'version' 必须为单一字符值！", call. = FALSE)
  }
  if (!is.character(source) || length(source) != 1) {
    stop("参数 'source' 必须为单一字符值！", call. = FALSE)
  }
  if (!is.character(log_file) || length(log_file) != 1) {
    stop("参数 'log_file' 必须为单一字符值！", call. = FALSE)
  }

  #-------------------- 确定文件路径 --------------------
  json_file <- file.path("marker_sets", "human", paste0(set_name, ".json"))
  if (!file.exists(json_file)) {
    json_file <- file.path("marker_sets", "mouse", paste0(set_name, ".json"))
    if (!file.exists(json_file)) {
      stop("marker set 文件 '", set_name, "' 不存在，请先用 create_marker_set 创建！", call. = FALSE)
    }
  }

  #-------------------- 读取现有 JSON 文件 --------------------
  cli::cli_alert_info("读取现有 marker set：{json_file}")
  json_data <- jsonlite::fromJSON(json_file, simplifyVector = TRUE)

  # 检查 JSON 结构
  if (!all(c("meta", "markers", "entries") %in% names(json_data))) {
    stop("JSON 文件格式错误，缺少 meta, markers 或 entries！", call. = FALSE)
  }
  species_std <- json_data$meta$species

  # 保存原 meta 信息
  old_meta <- json_data$meta

  # 打印原 meta 信息到窗口（使用 cli 美化）
  cli::cli_h2("原 meta 信息")
  cli::cli_ul()
  cli::cli_li("物种: {json_data$meta$species}")
  cli::cli_li("版本: {json_data$meta$version}")
  cli::cli_li("来源: {paste(json_data$meta$source, collapse = ', ')}")
  cli::cli_li("创建时间: {json_data$meta$created}")
  cli::cli_end()

  #-------------------- 提取现有数据 --------------------
  df_existing <- tibble::as_tibble(json_data$entries)

  #-------------------- 处理新条目 --------------------
  if (length(marker_genes) == 1 && grepl(",", marker_genes)) {
    new_genes <- strsplit(trimws(marker_genes), ",\\s*")[[1]]
  } else {
    new_genes <- as.character(marker_genes)
  }
  new_df <- tibble::tibble(
    cell_type = cell_type,
    marker_genes = list(new_genes),
    source = source
  )

  # 合并数据（覆盖同名 cell_type）
  if (cell_type %in% df_existing$cell_type) {
    cli::cli_alert_warning("细胞类型 '{cell_type}' 已存在，将覆盖原有记录")
    df_existing <- df_existing[df_existing$cell_type != cell_type, ]
  }
  df_combined <- dplyr::bind_rows(df_existing, new_df)

  #-------------------- 调用 create_marker_set 更新文件 --------------------
  source("Rutils/create_marker_set.R")
  create_marker_set(
    data = df_combined,
    species = species_std,
    set_name = set_name,
    version = version,
    log_file = "logs/marker_set/create_marker_set.log"  # 保持默认值
  )

  # 重新读取更新后的 JSON 文件以获取新 meta
  updated_json_data <- jsonlite::fromJSON(json_file, simplifyVector = TRUE)
  new_meta <- updated_json_data$meta

  # 打印新 meta 信息到窗口（使用 cli 美化）
  cli::cli_h2("新 meta 信息")
  cli::cli_ul()
  cli::cli_li("物种: {updated_json_data$meta$species}")
  cli::cli_li("版本: {updated_json_data$meta$version}")
  cli::cli_li("来源: {paste(updated_json_data$meta$source, collapse = ', ')}")
  cli::cli_li("创建时间: {updated_json_data$meta$created}")
  cli::cli_end()



  #-------------------- 写入日志 --------------------
  dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)

  # 原 meta 信息
  old_meta_str <- paste(
    "  物种             ：", old_meta$species, "\n",
    "  版本             ：", old_meta$version, "\n",
    "  来源             ：", paste(old_meta$source, collapse = ", "), "\n",
    "  创建时间         ：", old_meta$created, "\n"
  )

  # 新 meta 信息
  new_meta_str <- paste(
    "  物种             ：", new_meta$species, "\n",
    "  版本             ：", new_meta$version, "\n",
    "  来源             ：", paste(new_meta$source, collapse = ", "), "\n",
    "  创建时间         ：", new_meta$created, "\n"
  )

  cat(
    "\n", paste0(rep("-", 70), collapse = ""), "\n",
    "🧬 添加 marker entry 并更新：", basename(json_file), "\n",
    "📁 文件路径         ：", json_file, "\n",
    "🔬 原 meta 信息：\n", old_meta_str,
    "🔬 新 meta 信息：\n", new_meta_str,
    "🧠 添加的细胞类型   ：", cell_type, "\n",
    "🧬 marker genes     ：", paste(new_genes, collapse = ", "), "\n",
    "📌 来源             ：", source, "\n",
    "🕒 更新时间         ：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    paste0(rep("-", 70), collapse = ""), "\n",
    file = log_file, append = TRUE
  )

  cli::cli_alert_success("成功添加新条目并更新 marker set：{json_file}")
}

# 使用示例
# 先创建一个初始 marker set
df <- data.frame(
  cell_type = c("CD4+ T cells", "NK cells"),
  marker_genes = c("CD3D, IL7R, CCR7", "NKG7, GNLY"),
  source = "Manual"
)
source("Rutils/create_marker_set.R")
create_marker_set(data = df, species = "human", set_name = "pan_immune", version = "v1")

# # 添加新条目并更新
# source("Rutils/add_marker_entry.R")
add_marker_entry(
  set_name = "pan_immune",
  cell_type = "B cells",
  marker_genes = c("CD19", "CD79A", "MS4A1"),
  version = "v2",
  source = "Literature"
)
