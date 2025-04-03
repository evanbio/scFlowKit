# Rutils/load_data.R
#-------------------------------------------------------------------------------

# scFlowKit: Load Single-Cell RNA-seq Data
#-------------------------------------------------------------------------------


# load_data: 加载 10X Genomics 格式的单细胞 RNA-seq 数据
# 参数:
#   base_path: 数据基础路径（比如 "data/raw/"）
#   dataset_name: 数据集名称（比如 "pbmc3k"）
#   gene_column: 基因列编号（genes.tsv 中基因标识符列），默认 2
#   cell_column: 细胞列编号（barcodes.tsv 中条形码列），默认 1
#   unique_features: 是否保留唯一特征（去除重复基因名），默认 TRUE
#   strip_suffix: 是否去除基因名后缀，默认 FALSE
#   min_cells: 基因至少在多少细胞中表达，默认 3
#   min_features: 细胞至少表达多少基因，默认 40
#   project: Seurat 对象项目名称，默认 "scFlowKit"
#   assay: 测序类型，默认 "RNA"（表示 RNA 数据）

load_data <- function(base_path, dataset_name, 
                      gene_column = 2, cell_column = 1, 
                      unique_features = TRUE, strip_suffix = FALSE,
                      min_cells = 3, min_features = 40, 
                      project = "scFlowKit", assay = "RNA") {
  
  # 使用 cli 美化输出
  library(cli)

  # ---------------- 参数检查 ----------------
  # 验证输入参数是否为字符类型
  if (!is.character(base_path) || !is.character(dataset_name)) {
    stop("参数 'base_path' 和 'dataset_name' 必须为字符类型！", call. = FALSE)
  }

  # 验证 gene_column 和 cell_column 是否为正整数
  if (!is.numeric(gene_column) || gene_column < 1 || gene_column != as.integer(gene_column)) {
    stop("参数 'gene_column' 必须为正整数！", call. = FALSE)
  }
  if (!is.numeric(cell_column) || cell_column < 1 || cell_column != as.integer(cell_column)) {
    stop("参数 'cell_column' 必须为正整数！", call. = FALSE)
  }

  # 验证 unique_features 和 strip_suffix 是否为逻辑值
  if (!is.logical(unique_features)) {
    stop("参数 'unique_features' 必须为逻辑值！", call. = FALSE)
  }
  if (!is.logical(strip_suffix)) {
    stop("参数 'strip_suffix' 必须为逻辑值！", call. = FALSE)
  }

  # 验证 min_cells 和 min_features 是否为非负整数
  if (!is.numeric(min_cells) || min_cells < 0 || min_cells != as.integer(min_cells)) {
    stop("参数 'min_cells' 必须为非负整数！", call. = FALSE)
  }
  if (!is.numeric(min_features) || min_features < 0 || min_features != as.integer(min_features)) {
    stop("参数 'min_features' 必须为非负整数！", call. = FALSE)
  }

  # 验证 project 和 assay 是否为字符类型
  if (!is.character(project)) {
    stop("参数 'project' 必须为字符类型！", call. = FALSE)
  }
  if (!is.character(assay)) {
    stop("参数 'assay' 必须为字符类型！", call. = FALSE)
  }

  # ---------------- 路径与数据读取 ----------------
  # 拼接完整数据路径
  data_path <- file.path(base_path, dataset_name)

  # 检查数据路径是否存在
  if (!dir.exists(data_path)) {
    cli_alert_danger("数据路径不存在：{data_path}")
    stop()
  }
  
  cli_h2("📥 加载数据集 {.strong {dataset_name}}")
  cli_text("📂 路径：{data_path}")
  cli_text("🧪 过滤参数：min.cells = {min_cells}, min.features = {min_features}")

  # 加载 10X Genomics 数据
  data <- Seurat::Read10X(data.dir = data_path,
                          gene.column = gene_column,
                          cell.column = cell_column,
                          unique.features = unique_features,
                          strip.suffix = strip_suffix)
  
  # ---------------- 创建 Seurat 对象 ----------------
  seurat_obj <- Seurat::CreateSeuratObject(counts = data,
                                           min.cells = min_cells,
                                           min.features = min_features,
                                           project = project,
                                           assay = assay)

  # 提示用户加载完成
  cli_alert_success("✅ 加载完成：共 {ncol(seurat_obj)} 个细胞，{nrow(seurat_obj)} 个基因")
  return(seurat_obj)
}

#-------------------------------------------------------------------------------