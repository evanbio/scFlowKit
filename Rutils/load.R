# Rutils/load.R
#-------------------------------------------------------------------------------
# scFlowKit: Load Script for Dependencies (with cli styling)
#-------------------------------------------------------------------------------

# 使用 cli 风格输出
if (!requireNamespace("cli", quietly = TRUE)) {
  stop("请先安装 [cli] 包后再运行此脚本：install.packages('cli')", call. = FALSE)
}
library(cli)

cli_h1("scFlowKit · 依赖包加载器")

# 包列表及说明
required_packages <- list(
  Seurat         = "单细胞 RNA-seq 主工具包",
  patchwork      = "拼图组合图形",
  ggplot2        = "绘图基础",
  rlang          = "表达式解析支持",
  stringr        = "字符串操作",
  DoubletFinder  = "双细胞检测",
  MAST           = "差异表达分析",
  tidyverse      = "数据处理套件",
  glmGamPoi      = "SCTransform 加速器",
  harmony        = "批次整合",
  jsonlite       = "JSON 数据处理",
  EnhancedVolcano= "火山图绘制",
  SingleR        = "自动细胞注释",
  celldex        = "参考数据库",
  cli            = "美化控制台信息"
)

# 加载包并提示状态
for (pkg in names(required_packages)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cli_alert_danger("缺少包 {.pkg {pkg}} —— {required_packages[[pkg]]}")
    cli_alert_info("请运行 {.code source('Rutils/install.R')} 安装缺失的包。\n")
  } else {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    cli_alert_success("已加载包 {.pkg {pkg}} —— {required_packages[[pkg]]}")
  }
}

cli_h2("所有依赖包加载流程已完成 ✅")
cli_text("如有缺失包，请根据提示补充安装后重新运行。")

#-------------------------------------------------------------------------------
