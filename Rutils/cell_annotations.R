# cell_annotations.R
#-------------------------------------------------------------------------------

# scFlowKit: Load Cell Type Annotations
#-------------------------------------------------------------------------------

# 自定义函数：加载和处理细胞类型注释数据
#-------------------------------------------------------------------------------

# load_cell_annotations: 加载 JSON 格式的细胞类型注释数据
# 参数:
#   file_path: JSON 文件路径（默认 "docs/cell_annotation.json"）

load_cell_annotations <- function(file_path = "docs/cell_annotation.json") {
  # 验证输入参数是否为字符类型
  if (!is.character(file_path)) {
    stop("参数 'file_path' 必须为字符类型！", call. = FALSE)
  }
  
  # 检查文件是否存在
  if (!file.exists(file_path)) {
    stop("文件不存在: ", file_path, call. = FALSE)
  }
  
  # 提示用户正在加载数据
  message("正在加载细胞类型注释数据: [", file_path, "]")
  
  # 加载所需包
  library(jsonlite)
  
  # 读取 JSON 文件
  annotations <- fromJSON(file_path, simplifyDataFrame = FALSE)
  
  # 提取所有 source
  sources <- unique(sapply(annotations, function(x) x$source))
  
  # 按 source 分组
  result <- lapply(sources, function(src) {
    src_data <- annotations[sapply(annotations, function(x) x$source == src)]
    list(
      cell_type = sapply(src_data, function(x) x$cell_type),
      marker_genes = lapply(src_data, function(x) x$marker_genes),
      source = sapply(src_data, function(x) x$source)
    )
  })
  
  # 设置顶级键名为 source
  names(result) <- sources
  
  # 提示用户加载完成
  message("细胞类型注释数据加载完成！")
  
  return(result)
}

#-------------------------------------------------------------------------------

# extract_markers: 从注释数据中提取特定 source 的标记基因
# 参数:
#   cell_anno: load_cell_annotations 返回的注释对象
#   source: 要提取的 source 名称（比如 "HBC Tutorial"）

extract_markers <- function(cell_anno, source) {
  # 验证输入参数
  if (!is.list(cell_anno)) {
    stop("参数 'cell_anno' 必须为列表类型！", call. = FALSE)
  }
  if (!is.character(source)) {
    stop("参数 'source' 必须为字符类型！", call. = FALSE)
  }
  
  # 检查 source 是否存在
  if (!source %in% names(cell_anno)) {
    stop("Source 不存在: ", source, call. = FALSE)
  }
  
  # 提取该 source 的数据
  src_data <- cell_anno[[source]]
  
  # 检查 cell_type 是否有重复
  cell_types <- src_data$cell_type
  if (any(duplicated(cell_types))) {
    warning("在 source ", source, " 中发现重复的 cell_type: ",
            paste(cell_types[duplicated(cell_types)], collapse = ", "),
            call. = FALSE)
  }
  
  # 创建以 cell_type 为键的 markers 列表
  cell_markers <- setNames(src_data$marker_genes, cell_types)
  
  return(cell_markers)
}

#-------------------------------------------------------------------------------
# 示例用法（注释掉）
#-------------------------------------------------------------------------------
# # 加载脚本
# source("cell_annotations.R")
#
# # 加载默认路径的注释数据
# cell_anno <- load_cell_annotations()
#
# # 提取 "HBC Tutorial" 的标记基因
# cell_markers <- extract_markers(cell_anno, "HBC Tutorial")
#
# # 查看结果
# str(cell_markers)
# print(cell_markers[["CD8+ T cells"]])  # 输出: "CD3D" "CD8A"
#-------------------------------------------------------------------------------