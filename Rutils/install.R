#-------------------------------------------------------------------------------
# scFlowKit: Installation Script for Dependencies
#-------------------------------------------------------------------------------

# 自定义函数：安装 R 包
#-------------------------------------------------------------------------------

# inst_pkg: 安装 R 包（支持 CRAN、GitHub、Bioconductor、本地）
# 参数:
#   pkg: 需要安装的包名向量（CRAN/Bioconductor）或 GitHub 仓库名向量（格式 "user/repo"）
#   source: 安装来源，可选 "CRAN", "GitHub", "Bioconductor", "local"，默认 "CRAN"
#   path: 本地安装时指定路径（.tar.gz 或包目录），默认 NULL
#   ...: 其他参数（如 lib、dependencies 等）

inst_pkg <- function(pkg = NULL,
                     source = "CRAN",
                     path = NULL, ...) {
  # 定义可用的安装来源
  valid_sources <- c("CRAN", "GitHub", "Bioconductor", "local")
  
  # 如果来源不在可选范围内，给出自定义提示并停止
  if (!source %in% valid_sources) {
    stop("无效的安装来源: ", source,
         "。可选值: CRAN, GitHub, Bioconductor, local", call. = FALSE)
  }
  
  # 如果需要安装 CRAN/GitHub/Bioconductor 包，但没有提供 pkg 名，给出提示并停止
  if (is.null(pkg) && source != "local") {
    stop("请提供要安装的包名或 GitHub 仓库名！", call. = FALSE)
  }
  
  # 对于 CRAN/GitHub/Bioconductor，若包已安装则跳过
  # （本地安装可能不知道包名，故不做此判断或仅在 pkg 不为空时判断）
  if (!is.null(pkg) && source != "local") {
    for (p in pkg) {
      if (requireNamespace(p, quietly = TRUE)) {
        message("包 [", p, "] 已安装，跳过安装。")
        return(invisible(NULL))
      }
    }
  }
  
  # 根据不同来源安装
  if (source == "CRAN") {
    message("从 CRAN (清华大学镜像) 安装包: [", paste(pkg, collapse = ", "), "]")
    install.packages(pkg,
                     repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
                     ...)
    
  } else if (source == "GitHub") {
    message("从 GitHub 安装包: [", paste(pkg, collapse = ", "), "]")
    # 如果 devtools 没安装则先安装 devtools
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools",
                       repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
    }
    for (p in pkg) {
      devtools::install_github(p, ...)
    }
    
  } else if (source == "Bioconductor") {
    message("从 Bioconductor (西湖大学镜像) 安装包: [", paste(pkg, collapse = ", "), "]")
    # 如果 BiocManager 没安装则先安装 BiocManager
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager",
                       repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
    }
    # 暂存原镜像设置，on.exit 在函数结束时恢复
    old_mirror <- getOption("BioC_mirror")
    options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
    on.exit(options(BioC_mirror = old_mirror), add = TRUE)
    
    BiocManager::install(pkg, ...)
    
  } else if (source == "local") {
    # 如果是本地安装，则需要 path
    if (is.null(path)) {
      stop("请提供本地包路径（.tar.gz 或包目录）！", call. = FALSE)
    }
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
inst_pkg("tidyversy",source = "CRAN")
inst_pkg("glmGamPoi",source = "Bioconductor")
inst_pkg("harmony", source = "CRAN")
inst_pkg("jsonlite", source = "CRAN")
inst_pkg("EnhancedVolcano", source = "Bioconductor")

# 提示用户安装完成
cat("scFlowKit 依赖包安装完成！")
cat("请在运行 main.R 前确保所有包已正确安装！")

#-------------------------------------------------------------------------------