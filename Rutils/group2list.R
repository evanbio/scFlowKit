# Rutils/group2list.R
#-------------------------------------------------------------------------------
# 将数据框按键列分组，将值列聚合成列表
#-------------------------------------------------------------------------------
#
# 背景说明：
# 本函数用于将数据框按指定键列（key_col）分组，将值列（value_col）聚合成列表。
# 特别适用于为 create_marker_set 准备输入数据，将一对多关系转换为 cell_type 和 marker_genes 格式。
#
# 参数说明：
# - data      : 输入数据框
# - key_col   : 分组列名（用作键，例如细胞类型）
# - value_col : 值列名（聚合成列表，例如基因）
#
# 返回值：
# - 一个数据框，包含 key_col 和 value_col 对应的两列
#-------------------------------------------------------------------------------

group2list <- function(
  data,
  key_col,
  value_col
) {
  #-------------------- 依赖检查 --------------------
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("请先安装 dplyr 包！", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("请先安装 tibble 包！", call. = FALSE)
  }
  
  # 显式加载 dplyr 包，确保 %>% 可用
  library(dplyr)

  #-------------------- 参数校验 --------------------
  if (!is.data.frame(data)) {
    stop("参数 'data' 必须为数据框！", call. = FALSE)
  }
  if (!is.character(key_col) || length(key_col) != 1 || !key_col %in% names(data)) {
    stop("参数 'key_col' 必须为单一字符值且存在于数据框中！", call. = FALSE)
  }
  if (!is.character(value_col) || length(value_col) != 1 || !value_col %in% names(data)) {
    stop("参数 'value_col' 必须为单一字符值且存在于数据框中！", call. = FALSE)
  }

  #-------------------- 数据转换 --------------------
  # 按 key_col 分组，将 value_col 聚合成列表
  result <- data %>%
    dplyr::group_by(.data[[key_col]]) %>%
    dplyr::summarise(
      !!value_col := list(.data[[value_col]]),
      .groups = "drop"
    )

  # 转换为 tibble
  result <- tibble::as_tibble(result)

  return(result)
}

# # 使用示例：转换为 create_marker_set 输入
# df <- data.frame(
#   cell_type = rep("B_CELLS_MEMORY", 15),
#   marker_genes = c("AIM2", "BANK1", "BLK", "CD19", "CD27", "CD37", "CD69", "CD79A",
#                    "CD79B", "FAIM3", "FAM65B", "FCGR2B", "FCRL2", "GUSB", "HLA-DOB")
# )

# source("Rutils/group2list.R")
# converted_df <- group2list(
#   data = df,
#   key_col = "cell_type",
#   value_col = "marker_genes"
# )