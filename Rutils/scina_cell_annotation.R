# Rutils/scina_cell_annotation.R
#-------------------------------------------------------------------------------
# scFlowKit: Perform Cell Type Annotation Using SCINA
#-------------------------------------------------------------------------------
#
# 背景说明：
# SCINA（Semi-supervised Category Identification and Assignment）是一种
# 基于 marker 基因的半监督注释方法。通过已知 marker gene 集合与细胞表达谱的
# 拟合，自动为每个细胞赋予细胞类型标签。
#
# 本函数主要执行以下步骤：
# - 参数与输入校验
# - 运行 SCINA 注释（支持自定义迭代与收敛参数）
# - 保存日志与注释表
# - 返回包含原始对象、标签表、标签向量的结果列表
#
# 输入要求：
# - expr         : 标准化后的表达矩阵（基因 × 细胞），建议使用 logcounts 或 Seurat data slot
# - marker_list  : 命名 list，元素为每种细胞类型对应的 marker 基因向量，示例：
#                  list(
#                    "T_cells" = c("CD3D", "CD3E"),
#                    "B_cells" = c("CD79A", "MS4A1")
#                  )
#
# 参数说明：
# - expr                 : 表达矩阵，必须为标准化数据（不建议使用 raw counts）
# - marker_list          : 命名 list，每个元素为一类细胞的 marker 基因名（字符向量）
# - max_iter             : 最大迭代次数（默认 100）
# - convergence_n        : 收敛判据（默认 10）
# - convergence_rate     : 收敛速率（默认 0.99）
# - sensitivity_cutoff   : 灵敏度阈值（默认 1）
# - rm_overlap           : 是否自动移除 marker 重叠（默认 FALSE）
# - allow_unknown        : 是否允许部分细胞注释为 unknown（默认 TRUE）
# - log                  : 是否保存日志至 logs/cell_annotation/scina_cell_annotation.log（默认 TRUE）
#
# 返回结果：
# 返回一个列表，包含以下内容：
# - scina_cell_annotation : SCINA::SCINA 原始返回结果
# - anno_table             : 数据框，包含细胞 ID 与对应标签（列名：cell, scina_cell_label）
# - anno_vector            : 命名向量，cell_name → label，可用于直接写入 Seurat meta.data
#-------------------------------------------------------------------------------
scina_cell_annotation <- function(
  expr,
  marker_list,
  max_iter = 100,
  convergence_n = 10,
  convergence_rate = 0.99,
  sensitivity_cutoff = 1,
  rm_overlap = FALSE,
  allow_unknown = TRUE,
  log = TRUE
) {
  #-------------------- 依赖包检查 --------------------
  # 检查必要包
  if (!requireNamespace("SCINA", quietly = TRUE)) {
    stop("请先安装 SCINA 包！", call. = FALSE)
  }
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("请先安装 cli 包！", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("请先安装 tibble 包！", call. = FALSE)
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("请先安装 stringr 包！", call. = FALSE)
  }

  #-------------------- 参数校验 --------------------
  # 检查 expr
  if (!is.matrix(expr)) {
    stop("参数 'expr' 必须为矩阵类型！", call. = FALSE)
  }
  if (nrow(expr) == 0 || ncol(expr) == 0) {
    stop("参数 'expr' 不能为空矩阵！", call. = FALSE)
  }
  if (is.null(rownames(expr)) || is.null(colnames(expr))) {
    stop("参数 'expr' 必须包含行名（基因名）和列名（细胞名）！", call. = FALSE)
  }

  # 检查 marker_list
  if (!is.list(marker_list)) {
    stop("参数 'marker_list' 必须为 list 类型！", call. = FALSE)
  }
  if (length(marker_list) == 0) {
    stop("参数 'marker_list' 不能为空！", call. = FALSE)
  }
  if (is.null(names(marker_list)) || any(names(marker_list) == "")) {
    stop("参数 'marker_list' 必须为命名 list，且每个元素需有细胞类型名称！", call. = FALSE)
  }
  for (markers in marker_list) {
    if (!is.character(markers) || length(markers) == 0) {
      stop("参数 'marker_list' 的每个元素必须为非空字符向量！", call. = FALSE)
    }
  }

  # 检查数值参数
  if (!is.numeric(max_iter) || max_iter <= 0 || max_iter != round(max_iter)) {
    stop("参数 'max_iter' 必须为正整数！", call. = FALSE)
  }
  if (!is.numeric(convergence_n) || convergence_n <= 0 || convergence_n != round(convergence_n)) {
    stop("参数 'convergence_n' 必须为正整数！", call. = FALSE)
  }
  if (!is.numeric(convergence_rate) || convergence_rate <= 0 || convergence_rate > 1) {
    stop("参数 'convergence_rate' 必须为 0 到 1 之间的数值！", call. = FALSE)
  }
  if (!is.numeric(sensitivity_cutoff) || sensitivity_cutoff < 0) {
    stop("参数 'sensitivity_cutoff' 必须为非负数值！", call. = FALSE)
  }

  # 检查逻辑参数
  if (!is.logical(rm_overlap) || length(rm_overlap) != 1) {
    stop("参数 'rm_overlap' 必须为单一逻辑值！", call. = FALSE)
  }
  if (!is.logical(allow_unknown) || length(allow_unknown) != 1) {
    stop("参数 'allow_unknown' 必须为单一逻辑值！", call. = FALSE)
  }
  if (!is.logical(log) || length(log) != 1) {
    stop("参数 'log' 必须为单一逻辑值！", call. = FALSE)
  }

  #-------------------- 日志设置 --------------------
  log_file <- "logs/cell_annotation/scina_cell_annotation.log"
  if (log) {
    dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
    sink(log_file, append = TRUE, split = TRUE)
    cat("\n", strrep("-", 70), "\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Start of SCINA run\n")
    cat("• 细胞数       ：", ncol(expr), "\n")
    cat("• 基因数       ：", nrow(expr), "\n")
    cat("• 细胞类型数   ：", length(marker_list), "\n")
    cat("• max_iter     ：", max_iter, "\n")
    cat("• convergence_n：", convergence_n, "\n")
    cat("• convergence_rate ：", convergence_rate, "\n")
    cat("• sensitivity_cutoff：", sensitivity_cutoff, "\n")
    cat("• allow_unknown     ：", allow_unknown, "\n")
    cat("• rm_overlap        ：", rm_overlap, "\n")
    cat(strrep("-", 70), "\n")
  }

  #-------------------- 执行注释 --------------------
  cli::cli_h1("执行 SCINA 细胞类型注释")

  # 执行注释
  scina_cell_res <- SCINA::SCINA(
    expr = as.matrix(expr),
    signatures = marker_list,
    max_iter = max_iter,
    convergence_n = convergence_n,
    convergence_rate = convergence_rate,
    sensitivity_cutoff = sensitivity_cutoff,
    rm_overlap = rm_overlap,
    allow_unknown = allow_unknown,
    log_file = NULL  # 使用自定义日志
  )

  #-------------------- 构建注释表 --------------------
  anno_vector <- scina_cell_res$cell_labels
  names(anno_vector) <- colnames(expr)
  na_count <- sum(is.na(anno_vector))
  if (na_count > 0) {
    cli::cli_alert_warning("{na_count} 个细胞标签为 NA，已替换为 'ambiguous'")
  }
  anno_vector[is.na(anno_vector)] <- "ambiguous"

  anno_table <- tibble::tibble(
    cell = colnames(expr),
    scina_cell_label = anno_vector
  )

  #-------------------- 打印注释分布 --------------------
  cli::cli_h2("注释结果分布")
  cli::cli_alert_info("以下为 scina_cell_label 的计数统计：")
  print(table(anno_vector))

  #-------------------- 结束日志 --------------------
  if (log) {
    cat(strrep("-", 70), "\n")
    cat("SCINA annotation complete\n")
    cat("注释标签已统计，详见上方 R 控制台输出。\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - End of SCINA run\n")
    sink()
    cli::cli_alert_success("日志已保存至：{log_file}")
  }

  #-------------------- 保存注释结果 --------------------
  output_dir <- "results/tables/cell_annotation"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  outfile <- file.path(output_dir, "scina_cell_annotation.csv")
  write.csv(anno_table, outfile, row.names = FALSE)
  cli::cli_alert_success("注释结果已保存至：{outfile}")

  #-------------------- 返回结果 --------------------
  cli::cli_alert_success("细胞类型注释完成！")
  return(list(
    scina_cell_annotation = scina_cell_res,
    anno_table = anno_table,
    anno_vector = anno_vector
  ))
}