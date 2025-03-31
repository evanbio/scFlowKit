# Rutils/get_gene_annotation.R
#-------------------------------------------------------------------------------

# 基因注释工具：使用 biomaRt 从 Ensembl 获取人类或小鼠基因注释
#-------------------------------------------------------------------------------

# 背景介绍
# - 基因注释是生物信息分析的重要环节，提供基因的标准化标识及功能信息。
# - biomaRt 提供实时访问 Ensembl 数据库的接口，确保信息准确且及时。
# - 本函数从 Ensembl 提取如下信息：
#   - 基因Symbol（人类：hgnc_symbol，小鼠：mgi_symbol）
#   - Ensembl ID
#   - Entrez ID
#   - 基因描述（Description）
#   - 基因类型（gene_biotype，如 protein-coding, lncRNA 等）
#   - 染色体位置（chromosome、start、end、strand）
#   - Ensembl 数据库版本号和数据下载日期
#
# - 通过 remove_empty_symbol 和 remove_na_entrez 参数，可灵活控制是否去除Symbol或Entrez ID为空的记录。
#
# get_gene_annotation: 获取基因注释并存储为RDS文件
# 参数:
#   species: 物种选择，"human" 或 "mouse"（默认 "human"）
#   remove_empty_symbol: 是否去除Symbol为空的记录（默认 FALSE）
#   remove_na_entrez: 是否去除Entrez ID为NA的记录（默认 FALSE）
# 返回:
#   基因注释数据框

get_gene_annotation <- function(species = c("human", "mouse"), 
                                remove_empty_symbol = FALSE,
                                remove_na_entrez = FALSE) {
  
  # 加载必要的包
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("请先安装 biomaRt 包：BiocManager::install('biomaRt')", call. = FALSE)
  }
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("请先安装 cli 包：install.packages('cli')", call. = FALSE)
  }
  library(biomaRt)
  library(dplyr)
  library(cli)
  
  species <- match.arg(species)
  
  # 选择物种对应数据集与Symbol
  if (species == "human") {
    dataset <- "hsapiens_gene_ensembl"
    symbol_attr <- "hgnc_symbol"
  } else if (species == "mouse") {
    dataset <- "mmusculus_gene_ensembl"
    symbol_attr <- "mgi_symbol"
  }
  
  mart <- useMart("ensembl", dataset = dataset)
  
  # 获取Ensembl数据库版本
  ensembl_version <- biomaRt::listEnsemblArchives() %>% 
    filter(current_release == "*") %>% 
    pull(version)
  
  # 下载基因注释信息
  cli_alert_info("[{Sys.time()}] 正在从Ensembl获取 {.field {species}} 基因注释 (版本 {.val {ensembl_version}})...")
  annotation <- getBM(attributes = c("ensembl_gene_id",
                                     symbol_attr,
                                     "entrezgene_id",
                                     "gene_biotype",
                                     "chromosome_name",
                                     "start_position",
                                     "end_position",
                                     "strand",
                                     "description"),
                      mart = mart)
  
  # 标准化列名
  colnames(annotation) <- c("ensembl_id",
                            "symbol",
                            "entrez_id",
                            "gene_type",
                            "chromosome",
                            "start",
                            "end",
                            "strand",
                            "description")
  
  # 添加物种、版本、下载日期
  annotation <- annotation %>% 
    mutate(species = species,
           ensembl_version = ensembl_version,
           download_date = Sys.Date())
  
  # 根据参数去除空值
  if (remove_empty_symbol) {
    annotation <- annotation %>% filter(!is.na(symbol) & symbol != "")
  }
  
  if (remove_na_entrez) {
    annotation <- annotation %>% filter(!is.na(entrez_id))
  }
  
  cli_alert_success("[{Sys.time()}] 获取完成，共 {.val {nrow(annotation)}} 条记录")
  
  return(annotation)
}

#-------------------------------------------------------------------------------
# 示例脚本: 下载并保存四个基因注释RDS文件到 data/external
#-------------------------------------------------------------------------------
# dir.create("data/external", recursive = TRUE, showWarnings = FALSE)

# 人类基因注释，保留全部记录
# human_all <- get_gene_annotation("human")
# saveRDS(human_all, "data/external/human_gene_annotation_all.rds")

# 人类基因注释，去除Symbol空值和Entrez ID为NA
# human_filtered <- get_gene_annotation("human", TRUE, TRUE)
# saveRDS(human_filtered, "data/external/human_gene_annotation_filtered.rds")

# 小鼠基因注释，保留全部记录
# mouse_all <- get_gene_annotation("mouse")
# saveRDS(mouse_all, "data/external/mouse_gene_annotation_all.rds")

# 小鼠基因注释，去除Symbol空值和Entrez ID为NA
# mouse_filtered <- get_gene_annotation("mouse", TRUE, TRUE)
# saveRDS(mouse_filtered, "data/external/mouse_gene_annotation_filtered.rds")

#-------------------------------------------------------------------------------