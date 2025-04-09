# Rutils/install.R
#-------------------------------------------------------------------------------

# scFlowKit: Installation Script for Dependencies
#-------------------------------------------------------------------------------

# inst_pkg: 安装 R 包（支持 CRAN、GitHub、Bioconductor、本地）
# 参数:
#   pkg: 需要安装的包名向量（CRAN/Bioconductor）或 GitHub 仓库名向量（格式 "user/repo"）
#   source: 安装来源，可选 "CRAN", "GitHub", "Bioconductor", "local"，默认 "CRAN"
#   path: 本地安装时指定路径（.tar.gz 或包目录），默认 NULL
#   ...: 其他参数（如 lib、dependencies 等）

inst_pkg <- function(pkg = NULL,
                     source = c("CRAN", "GitHub", "Bioconductor", "local"),
                     path = NULL, ...) {
  # 1.1 标准化 source 参数：支持简写 + 大小写不敏感
  source_input <- tolower(source[1])
  source_matched <- switch(source_input,
                           "cran" = "CRAN",
                           "gh"   = "GitHub", "github" = "GitHub",
                           "bio"  = "Bioconductor", "bioc" = "Bioconductor", "bioconductor" = "Bioconductor",
                           "local" = "local",
                           stop("无效的 source 参数：", source_input,
                                "。可用选项：CRAN, GitHub, Bioconductor, local", call. = FALSE)
  )
  
  # 1.2 检查 pkg 是否为空（非 local 情况下）
  if (is.null(pkg) && source_matched != "local") {
    stop("请提供要安装的包名或 GitHub 仓库名！", call. = FALSE)
  }
  
  # 1.3 跳过已安装包（local 不检查）
  if (!is.null(pkg) && source_matched != "local") {
    pkg_name <- if (source_matched == "GitHub") basename(pkg) else pkg
    if (requireNamespace(pkg_name, quietly = TRUE)) {
      message("包 [", pkg_name, "] 已安装，跳过安装。")
      return(invisible(NULL))
    }
  }
  
  # 1.4 根据来源进行安装
  if (source_matched == "CRAN") {
    message("从 CRAN 安装包: [", paste(pkg, collapse = ", "), "]")
    install.packages(pkg,
                     repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
                     ...)
    
  } else if (source_matched == "GitHub") {
    message("从 GitHub 安装包: [", paste(pkg, collapse = ", "), "]")
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools",
                       repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
    }
    for (p in pkg) devtools::install_github(p, ...)
    
  } else if (source_matched == "Bioconductor") {
    message("从 Bioconductor 安装包: [", paste(pkg, collapse = ", "), "]")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager",
                       repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
    }
    old_mirror <- getOption("BioC_mirror")
    options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
    on.exit(options(BioC_mirror = old_mirror), add = TRUE)
    BiocManager::install(pkg, ...)
    
  } else if (source_matched == "local") {
    if (is.null(path)) stop("请提供本地包路径！", call. = FALSE)
    message("从本地路径安装包: [", path, "]")
    install.packages(path, repos = NULL, type = "source", ...)
  }
  
  message("安装完成！")
  invisible(NULL)
}

#-------------------------------------------------------------------------------

# 安装 scFlowKit 所需的依赖包
#-------------------------------------------------------------------------------

# 安装 Seurat 包（从 CRAN 安装）
inst_pkg("Seurat", source = "CRAN")
inst_pkg("patchwork", source = "CRAN")
inst_pkg("ggplot2", source = "CRAN")
inst_pkg("rlang", source = "CRAN")
inst_pkg("stringr", source = "CRAN")
inst_pkg("chris-mcginnis-ucsf/DoubletFinder", source = "GitHub")
inst_pkg("MAST", source = "Bioconductor")
inst_pkg("tidyverse",source = "CRAN")
inst_pkg("glmGamPoi",source = "Bioconductor")
inst_pkg("harmony", source = "CRAN")
inst_pkg("jsonlite", source = "CRAN")
inst_pkg("EnhancedVolcano", source = "Bioconductor")
inst_pkg("SingleR",source = "Bioconductor")
inst_pkg("celldex",source = "Bioconductor")
inst_pkg("cli", source = "CRAN")

# 提示用户安装完成
cat("scFlowKit 依赖包安装完成！")
cat("请在运行 main.R 前确保所有包已正确安装！")

#-------------------------------------------------------------------------------