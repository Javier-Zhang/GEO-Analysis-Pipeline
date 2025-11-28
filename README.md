# GEO Analysis Pipeline | GEO数据分析流程

<p align="center">
  <img src="https://www.ncbi.nlm.nih.gov/geo/img/geo_main.gif" alt="GEO Logo" width="200"/>
</p>

<p align="center">
  <strong>A comprehensive Shiny-based pipeline for GEO data analysis</strong><br>
  <strong>基于Shiny的完整GEO数据分析流程</strong>
</p>

<p align="center">
  <a href="#features">Features</a> •
  <a href="#installation">Installation</a> •
  <a href="#usage">Usage</a> •
  <a href="#example-datasets">Examples</a> •
  <a href="#license">License</a>
</p>

---

## 🎯 Overview | 概述

GEO Analysis Pipeline is an integrated R/Shiny application for analyzing Gene Expression Omnibus (GEO) datasets. It supports both **microarray** and **RNA-seq** data, providing a complete workflow from data download to differential expression analysis and report generation.

GEO Analysis Pipeline是一个集成的R/Shiny应用程序，用于分析基因表达综合数据库（GEO）数据集。它同时支持**微阵列**和**RNA-seq**数据，提供从数据下载到差异表达分析和报告生成的完整工作流程。

## ✨ Features | 功能特性

### Data Download | 数据下载
- 🔽 Download GEO datasets by GSE ID | 通过GSE ID下载GEO数据集
- 📁 Support for series matrix and supplementary files | 支持series matrix和supplementary文件
- 🗂️ Automatic directory structure creation | 自动创建目录结构

### Quality Control | 质量控制
- 📊 Expression distribution analysis | 表达分布分析
- 🔗 Sample correlation heatmap | 样本相关性热图
- 📈 PCA analysis and visualization | PCA分析和可视化
- ⚠️ Outlier detection and removal | 异常样本检测和移除

### Normalization | 标准化
- **Microarray | 微阵列**: RMA, GCRMA, MAS5, Quantile, VSN
- **RNA-seq**: VST, rlog (DESeq2), TMM, RLE (edgeR), voom (limma)
- 🔄 ComBat batch effect correction | ComBat批次效应校正
- 📉 Before/after comparison visualization | 标准化前后对比可视化

### Annotation | 注释
- 🏷️ Bioconductor org.db packages (org.Hs.eg.db, org.Mm.eg.db)
- 🌐 biomaRt online annotation | biomaRt在线注释
- 📋 GPL platform file parsing | GPL平台文件解析
- 🔀 Multiple probe collapsing methods | 多探针合并方法

### Differential Expression | 差异表达
- **Microarray | 微阵列**: limma
- **RNA-seq**: DESeq2, edgeR, limma-voom
- 🌋 Interactive volcano plot | 交互式火山图
- 📊 MA plot and heatmap | MA图和热图
- 📥 Result export (CSV, Excel, GCT) | 结果导出

### Reporting | 报告生成
- 📄 HTML/PDF report generation | HTML/PDF报告生成
- 📊 Interactive tables and plots | 交互式表格和图表
- 📋 Analysis parameter recording | 分析参数记录

## 📦 Installation | 安装

### Prerequisites | 前置要求

```r
# Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install core Bioconductor packages
BiocManager::install(c(
    "Biobase", "GEOquery", "limma", "affy",
    "DESeq2", "edgeR", "sva",
    "org.Hs.eg.db", "org.Mm.eg.db", "AnnotationDbi", "biomaRt"
))

# Install CRAN packages
install.packages(c(
    "shiny", "shinydashboard", "shinyWidgets",
    "ggplot2", "plotly", "pheatmap", "DT",
    "rmarkdown", "knitr", "openxlsx"
))
```

### Clone Repository | 克隆仓库

```bash
git clone https://github.com/Javier-Zhang/GEO-Analysis-Pipeline.git
cd GEO-Analysis-Pipeline
```

### Run Application | 运行应用

```r
# In R
shiny::runApp()
```

## 🚀 Usage | 使用方法

### 1. Data Download | 数据下载

1. Enter GSE ID (e.g., `GSE6791`) | 输入GSE ID
2. Select data type (Microarray/RNA-seq) | 选择数据类型
3. Click "Download" | 点击下载

### 2. Quality Control | 质量控制

1. Run QC analysis | 运行QC分析
2. Review outlier samples | 查看异常样本
3. Remove outliers if needed | 如需要，移除异常样本

### 3. Normalization | 标准化

1. Select normalization method | 选择标准化方法
2. Configure batch correction (optional) | 配置批次校正（可选）
3. Run normalization | 运行标准化

### 4. Annotation | 注释

1. Select annotation method | 选择注释方法
2. Choose organism and ID type | 选择物种和ID类型
3. Run annotation | 运行注释

### 5. Differential Expression | 差异表达

1. Select grouping variable | 选择分组变量
2. Set control and treatment groups | 设置对照组和实验组
3. Adjust thresholds | 调整阈值
4. Run analysis | 运行分析

### 6. Report | 报告

1. Configure report settings | 配置报告设置
2. Generate report | 生成报告
3. Download results | 下载结果

## 🔬 Example Datasets | 示例数据集

The pipeline includes example datasets for **cervical cancer research** | 包含**宫颈癌研究**示例数据集：

| GSE ID | Platform | Samples | Description |
|--------|----------|---------|-------------|
| GSE6791 | GPL570 | 42 | Cervical cancer vs Normal / 宫颈癌vs正常 |
| GSE7803 | GPL570 | 28 | Cervical cancer expression / 宫颈癌表达谱 |
| GSE9750 | GPL96 | 66 | Cervical cancer gene expression / 宫颈癌基因表达 |
| GSE63514 | GPL570 | 128 | Cervical cancer progression / 宫颈癌进展研究 |
| GSE39001 | GPL570 | 108 | HPV-associated cervical cancer / HPV相关宫颈癌 |

## 📁 Project Structure | 项目结构

```
GEO-Analysis-Pipeline/
├── app.R                          # Main Shiny application
├── R/
│   ├── microarray/
│   │   ├── download.R             # Data download
│   │   ├── qc.R                   # Quality control
│   │   ├── normalize.R            # Normalization
│   │   ├── annotation.R           # Annotation
│   │   ├── differential.R         # Differential analysis
│   │   └── report.R               # Report generation
│   ├── rnaseq/
│   │   ├── download.R             # RNA-seq download
│   │   ├── qc.R                   # RNA-seq QC
│   │   ├── quantification.R       # Quantification
│   │   ├── normalize.R            # RNA-seq normalization
│   │   ├── annotation.R           # RNA-seq annotation
│   │   ├── differential.R         # RNA-seq DEG
│   │   └── report.R               # RNA-seq report
│   ├── annotation/
│   │   ├── bioconductor.R         # org.db annotation
│   │   ├── biomart.R              # biomaRt annotation
│   │   └── gpl_parser.R           # GPL file parsing
│   └── utils/
│       ├── phenotype.R            # Phenotype processing
│       ├── correlation.R          # Correlation analysis
│       └── export.R               # Data export
├── modules/                       # Shiny UI/Server modules
│   ├── mod_download.R
│   ├── mod_qc.R
│   ├── mod_normalize.R
│   ├── mod_annotate.R
│   ├── mod_differential.R
│   └── mod_report.R
├── templates/                     # Report templates
│   ├── microarray_qc_report.Rmd
│   ├── rnaseq_qc_report.Rmd
│   └── analysis_report.Rmd
├── www/
│   └── styles.css                 # Custom styles
├── tests/
│   └── test_pipeline.R            # Unit tests
├── DESCRIPTION                    # Package description
├── NAMESPACE                      # Namespace exports
└── README.md                      # Documentation
```

## 🔧 Dependencies | 依赖

### Core | 核心

- R (>= 4.0.0)
- shiny, shinydashboard, shinyWidgets
- Biobase, GEOquery
- limma, affy
- ggplot2, plotly, pheatmap, DT

### Optional | 可选

- DESeq2, edgeR (RNA-seq analysis)
- sva (batch correction)
- org.Hs.eg.db, org.Mm.eg.db (annotation)
- biomaRt (online annotation)

## 📝 Citation | 引用

If you use this pipeline in your research, please cite | 如果您在研究中使用此流程，请引用:

```
GEO Analysis Pipeline: A Shiny-based workflow for GEO data analysis
https://github.com/Javier-Zhang/GEO-Analysis-Pipeline
```

## 📄 License | 许可证

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

本项目采用MIT许可证 - 详情请见[LICENSE](LICENSE)文件。

## 🤝 Contributing | 贡献

Contributions are welcome! Please feel free to submit a Pull Request.

欢迎贡献！请随时提交Pull Request。

## 📧 Contact | 联系

For questions or suggestions, please open an issue on GitHub.

如有问题或建议，请在GitHub上提交issue。

---

<p align="center">
Made with ❤️ for bioinformatics research<br>
为生物信息学研究而创建
</p>
