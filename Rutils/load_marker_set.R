# Rutils/load_marker_set.R
#-------------------------------------------------------------------------------
# 加载 marker set 数据集
#-------------------------------------------------------------------------------
#
# 背景说明：
# 本函数用于从 marker_sets 目录加载指定的 marker set JSON 文件，
# 并返回整个数据集的结构，包括 meta, markers 和 entries。
# entries 部分被转换为数据框格式，包含 cell_type, marker_genes, source 三列。
#
# 参数说明：
# - species   : 物种名称（"human" 或 "mouse"）
# - set_name  : marker set 名称（不带 .json 后缀）
# - log_file  : 日志保存路径，默认 "logs/marker_set/load_marker_set.log"
#
# 返回值：
# - 一个列表，包含以下元素：
#   - meta: 元信息（物种、版本等）
#   - markers: 命名列表，细胞类型到 marker genes 的映射（例如 list("T_cells" = c("CD3D", "CD3E")))
#   - entries: 数据框，包含 cell_type, marker_genes, source 三列
#-------------------------------------------------------------------------------

load_marker_set <- function(
  species,
  set_name,
  log_file = "logs/marker_set/load_marker_set.log"
) {
  #-------------------- 依赖检查 --------------------
  for (pkg in c("tibble", "jsonlite", "cli", "dplyr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("请先安装 ", pkg, " 包！", call. = FALSE)
    }
  }

  # 显式加载 dplyr 包，确保 %>% 可用
  library(dplyr)

  #-------------------- 参数校验 --------------------
  if (!is.character(species) || length(species) != 1 || !tolower(species) %in% c("human", "mouse")) {
    stop("参数 'species' 必须为 'human' 或 'mouse'！", call. = FALSE)
  }
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
  if (!is.character(log_file) || length(log_file) != 1) {
    stop("参数 'log_file' 必须为单一字符值！", call. = FALSE)
  }

  #-------------------- 确定文件路径 --------------------
  json_file <- file.path("marker_sets", tolower(species), paste0(set_name, ".json"))
  if (!file.exists(json_file)) {
    stop("marker set 文件 '", json_file, "' 不存在！", call. = FALSE)
  }

  #-------------------- 读取 JSON 文件 --------------------
  cli::cli_alert_info("加载 marker set：{json_file}")
  json_data <- jsonlite::fromJSON(json_file, simplifyVector = FALSE)

  # 检查 JSON 结构
  if (!all(c("meta", "markers", "entries") %in% names(json_data))) {
    stop("JSON 文件格式错误，缺少 meta, markers 或 entries！", call. = FALSE)
  }

  #-------------------- 提取 meta 和 markers --------------------
  meta <- json_data$meta
  markers <- json_data$markers

  # 校验 markers 结构：必须是列表（允许为空）
  if (!is.list(markers)) {
    stop("markers 必须是列表！", call. = FALSE)
  }

  # 如果 markers 不为空，转换每个元素为字符向量
  if (length(markers) > 0) {
    # 检查每个元素是否是 list 或 vector
    is_valid <- sapply(markers, function(x) is.list(x) || is.vector(x))
    if (!all(is_valid)) {
      stop("markers 元素必须是 list 或 vector！问题元素：{names(markers)[!is_valid]}", call. = FALSE)
    }
    # 将每个元素展平为字符向量
    markers <- lapply(markers, function(x) as.character(unlist(x)))
    # 再次检查是否成功转换为字符向量
    is_character <- sapply(markers, is.character)
    if (!all(is_character)) {
      stop("markers 元素无法转换为字符向量！问题元素：{names(markers)[!is_character]}", call. = FALSE)
    }
  }

  # 打印 meta 信息到窗口（使用 cli 美化）
  cli::cli_h2("meta 信息")
  cli::cli_ul()
  cli::cli_li("物种: {meta$species}")
  cli::cli_li("版本: {meta$version}")
  cli::cli_li("来源: {paste(meta$source, collapse = ', ')}")
  cli::cli_li("创建时间: {meta$created}")
  cli::cli_end()

  #-------------------- 转换 entries 为数据框 --------------------
  entries_tbl <- tibble::tibble(
    cell_type = sapply(json_data$entries, function(x) x$cell_type),
    marker_genes = lapply(json_data$entries, function(x) x$marker_genes),
    source = sapply(json_data$entries, function(x) x$source)
  )

  # 确保 cell_type 和 source 是字符向量
  entries_tbl <- entries_tbl %>%
    dplyr::mutate(
      cell_type = as.character(cell_type),
      source = as.character(source)
    )

  #-------------------- 构造返回列表 --------------------
  result <- list(
    meta = meta,
    markers = markers,
    entries = entries_tbl
  )

  #-------------------- 写入日志 --------------------
  dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
  cat(
    "\n", paste0(rep("-", 70), collapse = ""), "\n",
    "🧬 加载 marker set：", basename(json_file), "\n",
    "📁 文件路径         ：", json_file, "\n",
    "🔬 物种             ：", meta$species, "\n",
    "🧠 细胞类型         ：", paste(entries_tbl$cell_type, collapse = ", "), "\n",
    "📌 来源             ：", paste(unique(entries_tbl$source), collapse = ", "), "\n",
    "📌 版本             ：", meta$version, "\n",
    "🕒 加载时间         ：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    paste0(rep("-", 70), collapse = ""), "\n",
    file = log_file, append = TRUE
  )

  cli::cli_alert_success("成功加载 marker set：{json_file}")
  return(result)
}

# # 使用示例：加载 marker set
# source("Rutils/load_marker_set.R")
# marker_set <- load_marker_set(
#   species = "human",
#   set_name = "pbmc_22_10x"
# )

# # 查看 markers
# print(marker_set$markers)
