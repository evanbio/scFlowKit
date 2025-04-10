# Rutils/download_celldex_refs.R
#-------------------------------------------------------------------------------
# scFlowKit: 下载并保存 celldex 提供的 7 个参考数据集
#-------------------------------------------------------------------------------
#
# 本脚本用于批量下载 celldex 提供的免疫与组织特异性参考数据集，并以 `.rds` 
# 格式保存在 `data/external/` 文件夹下。同时生成一个 `celldex_reference.rda` 文件，
# 其中包含所有参考对象，供自动注释流程调用。
#
# 参考数据集来源：celldex 包（https://bioconductor.org/packages/celldex）
#
# 数据集包括：
#   1. BlueprintEncodeData                - 人类免疫系统中的主要细胞亚群（Blueprint & ENCODE）
#   2. DatabaseImmuneCellExpressionData   - 多个数据库整合的免疫细胞表达谱
#   3. HumanPrimaryCellAtlasData          - 人类原代细胞表达谱（来自 Mabbott et al.）
#   4. ImmGenData                         - 小鼠免疫细胞图谱（ImmGen 项目）
#   5. MonacoImmuneData                   - 人类外周血免疫亚群表达谱（Monaco et al.）
#   6. MouseRNAseqData                    - 小鼠组织来源免疫细胞的 RNA-seq 表达数据
#   7. NovershternHematopoieticData       - 人类造血干/祖细胞谱系表达数据（Novershtern et al.）
#
# 输出文件：
#   - data/external/*.rds（共 7 个，单独保存每个数据集）
#   - data/external/celldex_reference.rda（一个文件包含全部数据）
#-------------------------------------------------------------------------------

#-------------------- 环境准备 --------------------
if (!requireNamespace("celldex", quietly = TRUE)) {
  BiocManager::install("celldex")
}
if (!requireNamespace("cli", quietly = TRUE)) {
  install.packages("cli")
}

library(celldex)
library(cli)

output_dir <- "data/external"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

#-------------------- 下载并保存 --------------------
cli::cli_h1("📦 下载 celldex 参考数据集")

cli::cli_alert_info("开始下载 BlueprintEncodeData ...")
BlueprintEncodeData <- celldex::BlueprintEncodeData()
saveRDS(BlueprintEncodeData, file.path(output_dir, "BlueprintEncodeData.rds"))

cli::cli_alert_info("开始下载 DatabaseImmuneCellExpressionData ...")
DatabaseImmuneCellExpressionData <- celldex::DatabaseImmuneCellExpressionData()
saveRDS(DatabaseImmuneCellExpressionData, file.path(output_dir, "DatabaseImmuneCellExpressionData.rds"))

cli::cli_alert_info("开始下载 HumanPrimaryCellAtlasData ...")
HumanPrimaryCellAtlasData <- celldex::HumanPrimaryCellAtlasData()
saveRDS(HumanPrimaryCellAtlasData, file.path(output_dir, "HumanPrimaryCellAtlasData.rds"))

cli::cli_alert_info("开始下载 ImmGenData ...")
ImmGenData <- celldex::ImmGenData()
saveRDS(ImmGenData, file.path(output_dir, "ImmGenData.rds"))

cli::cli_alert_info("开始下载 MonacoImmuneData ...")
MonacoImmuneData <- celldex::MonacoImmuneData()
saveRDS(MonacoImmuneData, file.path(output_dir, "MonacoImmuneData.rds"))

cli::cli_alert_info("开始下载 MouseRNAseqData ...")
MouseRNAseqData <- celldex::MouseRNAseqData()
saveRDS(MouseRNAseqData, file.path(output_dir, "MouseRNAseqData.rds"))

cli::cli_alert_info("开始下载 NovershternHematopoieticData ...")
NovershternHematopoieticData <- celldex::NovershternHematopoieticData()
saveRDS(NovershternHematopoieticData, file.path(output_dir, "NovershternHematopoieticData.rds"))

#-------------------- 保存为合集 RDA --------------------
cli::cli_alert_info("保存所有对象为 celldex_reference.rda ...")
save(
  BlueprintEncodeData,
  DatabaseImmuneCellExpressionData,
  HumanPrimaryCellAtlasData,
  ImmGenData,
  MonacoImmuneData,
  MouseRNAseqData,
  NovershternHematopoieticData,
  file = file.path(output_dir, "celldex_reference.rda")
)

cli::cli_h2("✅ 所有 celldex 数据集保存完毕！")
