# Rutils/create_marker_set.R
#-------------------------------------------------------------------------------
# 创建标准化 marker set JSON 文件（支持从空或已有数据创建）
#-------------------------------------------------------------------------------
#
# 背景说明：
# 本函数用于创建标准化的标记基因集（marker set）JSON 文件，存储在 marker_sets/ 目录下，
# 支持从空数据或已有数据框创建。物种名称和基因大小写会根据输入自动标准化。
#
# 本函数主要执行以下步骤：
# - 检查依赖包和输入参数
# - 标准化物种名称和基因大小写
# - 处理输入数据并生成 markers 和 entries
# - 写入 JSON 文件并记录日志
#
# 参数说明：
# - data       : 可以为 NULL（创建空集），也可为包含 cell_type, marker_genes, source 的 data.frame
# - species    : 支持多种人鼠拼写形式，自动标准化为 "Human"/"Mouse"，决定基因名大小写
# - set_name   : marker set 文件名（不带 .json 后缀）
# - version    : marker set 的版本号，默认 "v1"
# - log_file   : 日志保存路径，默认写入 logs/marker_set/create_marker_set.log
#
# 返回值：
# - 无返回值，直接写入 JSON 文件，并记录日志
#-------------------------------------------------------------------------------

create_marker_set <- function(
    data = NULL,
    species = c("human", "mouse", "homo sapiens", "mus musculus", "hs", "mm", "Human", "Mouse"),
    set_name,
    version = "v1",
    log_file = "logs/marker_set/create_marker_set.log"
) {
  #-------------------- 依赖检查 --------------------
  for (pkg in c("tibble", "stringr", "jsonlite", "cli", "dplyr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("请先安装 ", pkg, " 包！", call. = FALSE)
    }
  }
  
  #-------------------- 参数校验 --------------------
  species <- species[1]
  if (!is.character(species) || length(species) != 1) {
    stop("参数 'species' 必须为单一字符值！", call. = FALSE)
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
  if (!is.character(version) || length(version) != 1) {
    stop("参数 'version' 必须为单一字符值！", call. = FALSE)
  }
  if (!is.character(log_file) || length(log_file) != 1) {
    stop("参数 'log_file' 必须为单一字符值！", call. = FALSE)
  }
  
  #-------------------- 物种匹配 + 基因大小写 --------------------
  species_std <- dplyr::case_match(
    tolower(species),
    "human" ~ "Human",
    "homo sapiens" ~ "Human",
    "hs" ~ "Human",
    "mouse" ~ "Mouse",
    "mus musculus" ~ "Mouse",
    "mm" ~ "Mouse",
    .default = NA
  )
  if (is.na(species_std)) stop("暂不支持该物种：", species)
  
  gene_case <- if (species_std == "Human") "upper" else "lower"
  
  #-------------------- 数据处理 --------------------
  if (is.null(data)) {
    df <- tibble::tibble(
      cell_type = character(),
      marker_genes = list(),
      source = character()
    )
  } else {
    if (!is.data.frame(data)) stop("参数 'data' 必须为 data.frame！", call. = FALSE)
    required_cols <- c("cell_type", "marker_genes", "source")
    if (!all(required_cols %in% colnames(data))) {
      stop("数据框必须包含列：cell_type, marker_genes, source！", call. = FALSE)
    }
    df <- tibble::as_tibble(data[, required_cols])
    
    if (!is.character(df$cell_type) || any(df$cell_type == "")) {
      stop("列 'cell_type' 必须为非空字符向量！", call. = FALSE)
    }
    if (!is.character(df$source)) {
      stop("列 'source' 必须为字符向量！", call. = FALSE)
    }
    if (is.character(df$marker_genes)) {
      df$marker_genes <- lapply(df$marker_genes, function(x) strsplit(trimws(x), ",\\s*")[[1]])
    } else if (!is.list(df$marker_genes)) {
      stop("列 'marker_genes' 必须为字符向量或 list！", call. = FALSE)
    }
    df$marker_genes <- lapply(df$marker_genes, function(g) {
      genes <- unique(trimws(as.character(g)))
      if (length(genes) == 0 || any(genes == "")) {
        stop("marker_genes 中存在空值或无效基因名！", call. = FALSE)
      }
      switch(gene_case,
             upper = toupper(genes),
             lower = tolower(genes),
             genes
      )
    })
  }
  
  #-------------------- 转换为 JSON 兼容结构 --------------------
  marker_list <- setNames(df$marker_genes, df$cell_type)
  entries_list <- lapply(seq_len(nrow(df)), function(i) {
    list(
      cell_type = df$cell_type[i],
      marker_genes = df$marker_genes[[i]],
      source = df$source[i]
    )
  })
  json_data <- list(
    meta = list(
      species = species_std,
      source = if (nrow(df) > 0) unique(df$source) else "Empty",
      version = version,
      created = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ),
    markers = marker_list,
    entries = entries_list
  )
  
  #-------------------- 写入 JSON 文件 --------------------
  out_dir <- file.path("marker_sets", tolower(species_std))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  outfile <- file.path(out_dir, paste0(set_name, ".json"))
  
  jsonlite::write_json(json_data, path = outfile, pretty = TRUE, auto_unbox = TRUE)
  
  #-------------------- 写入日志 --------------------
  dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
  cat(
    "\n", paste0(rep("-", 70), collapse = ""), "\n",
    "🧬 创建 marker set：", set_name, "\n",
    "📁 文件路径        ：", outfile, "\n",
    "🔬 物种            ：", species_std, "\n",
    "🧠 基因格式        ：", gene_case, "\n",
    "🧾 条目数量        ：", nrow(df), "\n",
    "📌 版本            ：", version, "\n",
    "🕒 创建时间        ：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    paste0(rep("-", 70), collapse = ""), "\n",
    file = log_file, append = TRUE
  )
  
  #-------------------- 提示 --------------------
  cli::cli_alert_success("成功写入 marker set：{outfile}")
}

#-------------------- 示例 ----------------------

# # 示例 1：创建简单的免疫细胞 marker set
# df_simple <- data.frame(
#   cell_type = c("CD4+ T cells", "NK cells", "B cells"),
#   marker_genes = c("CD3D, IL7R, CCR7", "NKG7, GNLY", "CD19, CD79A, MS4A1"),
#   source = "Manual curation from literature"
# )

# create_marker_set(
#   data = df_simple,
#   species = "human",
#   set_name = "pan_immune_v2",
#   version = "v2"
# )

# # 示例 2：创建空的 marker set，用于后续手动填充
# create_marker_set(
#   data = NULL,
#   species = "mouse",
#   set_name = "brain_cells_v1",
#   version = "v1"
# )

# # 示例 3：创建复杂的 marker set，包含多种细胞类型和来源
# df_complex <- data.frame(
#   cell_type = c("CD8+ T cells", "Macrophages", "Neurons", "Astrocytes"),
#   marker_genes = c("CD8A, GZMB, PRF1", "CD68, CD163", "RBFOX3, SNAP25", "GFAP, S100B"),
#   source = c("10x Genomics PBMC dataset", "Manual", "Allen Brain Atlas", "Manual")
# )

# create_marker_set(
#   data = df_complex,
#   species = "human",
#   set_name = "mixed_tissues_v1",
#   version = "v1",
#   log_file = "logs/marker_set/custom_marker_set.log"
# )

# # 示例 4：使用 list 格式的 marker_genes 输入
# df_list <- data.frame(
#   cell_type = c("Monocytes", "Dendritic cells"),
#   marker_genes = I(list(c("CD14", "FCGR3A"), c("CD1C", "HLA-DR"))),
#   source = "Published scRNA-seq study"
# )

# create_marker_set(
#   data = df_list,
#   species = "human",
#   set_name = "myeloid_cells_v1",
#   version = "v1"
# )

marker_genes = I(list(c("CD14", "FCGR3A"), c("CD1C", "HLA-DR")))
