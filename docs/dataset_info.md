# 数据集信息

---

## 基本信息

- **数据集**：5k Human PBMCs (Donor 1–4)
- **来源**：10x Genomics 公开数据集
- **描述**：Universal 3' Gene Expression datasets analyzed using Cell Ranger 9.0.0
- **下载链接**：
  - [Donor 1](https://www.10xgenomics.com/datasets/5k-human-pbmcs-donor-1-universal-3-prime-gene-expression-dataset-analyzed-using-cell-ranger-9-0-0)
  - [Donor 2](https://www.10xgenomics.com/datasets/5k-human-pbmcs-donor-2-universal-3-prime-gene-expression-dataset-analyzed-using-cell-ranger-9-0-0)
  - [Donor 3](https://www.10xgenomics.com/datasets/5k-human-pbmcs-donor-3-universal-3-prime-gene-expression-dataset-analyzed-using-cell-ranger-9-0-0)
  - [Donor 4](https://www.10xgenomics.com/datasets/5k-human-pbmcs-donor-4-universal-3-prime-gene-expression-dataset-analyzed-using-cell-ranger-9-0-0)

---

## 实验背景（适用于全部样本）

- 外周血单核细胞（PBMC）取自健康供体，由 10x Genomics 从 Cellular Technologies Limited 获取。
- 所有样本使用 Chromium GEM-X Single Cell 3' Reagent Kits v4 和 Illumina NovaSeq 6000 测序，平均每个细胞约 36,000 个读对。

---

## 测序参数

- **Read 1**：28 个周期
- **i7 索引**：10 个周期
- **i5 索引**：10 个周期
- **Read 2**：90 个周期

---

## 数据信息（按样本）

| Donor | 细胞数（估算） | 数据路径                      | 文件格式 | 存放文档目录         |
|-------|----------------|-------------------------------|-----------|----------------------|
| 1     | ~5,000         | `data/raw/5k_pbmc_donor1/`    | `.mtx.gz` 等 | `docs/5k_pbmc_donor1/` |
| 2     | ~5,000         | `data/raw/5k_pbmc_donor2/`    | `.mtx.gz` 等 | `docs/5k_pbmc_donor2/` |
| 3     | ~5,000         | `data/raw/5k_pbmc_donor3/`    | `.mtx.gz` 等 | `docs/5k_pbmc_donor3/` |
| 4     | ~5,000         | `data/raw/5k_pbmc_donor4/`    | `.mtx.gz` 等 | `docs/5k_pbmc_donor4/` |

---

## 使用说明

- 所有数据为过滤后的基因表达矩阵（filtered feature-barcode matrix），可直接用于下游分析。
- 压缩格式（`.gz`）可直接通过 Seurat 的 `Read10X()` 函数读取。

---

## 数据许可

Creative Commons Attribution 4.0 International (CC BY 4.0)

