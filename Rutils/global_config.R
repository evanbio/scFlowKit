# Rutils/global_config.R
#-------------------------------------------------------------------------------
# scFlowKit: 全局参数配置文件
#-------------------------------------------------------------------------------

# 📁 原始数据路径
# 示例：如果你的数据在 data/raw/5k_pbmc_donor1/ 目录下，则此路径应为 data/raw
data_path <- "data/raw"

# 📤 结果输出路径（图表/表格等将保存在此目录下）
# 示例：results/run1、results/test_batch 等
output_dir <- "results"

# 💾 预处理数据保存路径（用于存放 .rds 等）
processed_data_dir <- "data/processed"

# 📝 日志文件保存路径（可选，用于记录流程日志）
log_path <- "logs"

# 🌐 统一语言为英文（便于调试报错信息）
Sys.setenv(LANGUAGE = "en")

# 🧠 关闭字符转因子（避免生信数据类型混乱）
options(stringsAsFactors = FALSE)

# ✅ 控制台提示（加载配置）
cat("✅ 已加载全局配置文件。\n")
cat("📁 原始数据路径：", data_path, "\n")
cat("📤 结果输出路径：", output_dir, "\n")
cat("💾 预处理数据路径：", processed_data_dir, "\n")
cat("📝 日志保存路径：", log_path, "\n")
