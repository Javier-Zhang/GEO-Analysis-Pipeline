#' ============================================================================
#' GEO 微阵列分析报告生成模块
#' Microarray Analysis Report Generation Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 生成HTML质控报告 / Generate HTML QC report
#' - 整合分析结果 / Integrate analysis results
#' - 参数记录 / Parameter recording
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

suppressPackageStartupMessages({
  library(rmarkdown)
  library(knitr)
})

#' 生成微阵列分析报告
#' Generate Microarray Analysis Report
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param qc_results QC结果（来自run_microarray_qc）/ QC results (from run_microarray_qc)
#' @param deg_results 差异分析结果（来自run_differential_analysis）/ DEG results
#' @param output_file 输出文件路径 / Output file path
#' @param title 报告标题 / Report title
#' @param template 模板路径 / Template path
#'
#' @return 字符串，生成的报告路径 / Character, path to generated report
#' @export
#'
#' @examples
#' \dontrun{
#' # 生成完整分析报告 / Generate complete analysis report
#' report_path <- generate_microarray_report(
#'   eset = eset,
#'   qc_results = qc_results,
#'   deg_results = deg_results,
#'   output_file = "results/analysis_report.html"
#' )
#' }
generate_microarray_report <- function(eset,
                                        qc_results = NULL,
                                        deg_results = NULL,
                                        output_file = "microarray_report.html",
                                        title = "Microarray Analysis Report / 微阵列分析报告",
                                        template = NULL) {
  
  message("生成分析报告... / Generating analysis report...")
  
  # 确保输出目录存在 / Ensure output directory exists
  output_dir <- dirname(output_file)
  if (output_dir != "." && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 准备报告数据 / Prepare report data
  report_data <- list(
    title = title,
    date = Sys.Date(),
    eset = eset,
    qc_results = qc_results,
    deg_results = deg_results,
    dataset_info = get_dataset_info(eset),
    parameters = get_analysis_parameters(qc_results, deg_results)
  )
  
  # 选择模板 / Choose template
  if (is.null(template)) {
    template <- system.file("templates/microarray_qc_report.Rmd", package = "GEOAnalysis")
    if (template == "") {
      # 使用内置模板 / Use built-in template
      template <- create_temp_template()
    }
  }
  
  tryCatch({
    # 渲染报告 / Render report
    rmarkdown::render(
      input = template,
      output_file = basename(output_file),
      output_dir = output_dir,
      params = report_data,
      quiet = TRUE
    )
    
    message(paste0("报告已生成 / Report generated: ", output_file))
    return(output_file)
    
  }, error = function(e) {
    warning(paste0("报告生成失败 / Report generation failed: ", conditionMessage(e)))
    # 尝试生成简化报告 / Try to generate simplified report
    return(generate_simple_report(report_data, output_file))
  })
}


#' 获取数据集信息
#' Get Dataset Information
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @return 列表，数据集信息 / List with dataset information
get_dataset_info <- function(eset) {
  
  info <- list(
    n_samples = ncol(exprs(eset)),
    n_features = nrow(exprs(eset)),
    sample_names = sampleNames(eset),
    platform = tryCatch(annotation(eset), error = function(e) "Unknown"),
    expression_range = range(exprs(eset), na.rm = TRUE),
    pheno_vars = varLabels(eset)
  )
  
  return(info)
}


#' 获取分析参数
#' Get Analysis Parameters
#'
#' @param qc_results QC结果 / QC results
#' @param deg_results 差异分析结果 / DEG results
#' @return 列表，分析参数 / List with analysis parameters
get_analysis_parameters <- function(qc_results = NULL, deg_results = NULL) {
  
  params <- list(
    qc = list(),
    deg = list()
  )
  
  if (!is.null(qc_results)) {
    params$qc <- list(
      outlier_threshold = 0.8,
      n_outliers = length(qc_results$outliers$all_outliers)
    )
  }
  
  if (!is.null(deg_results)) {
    params$deg <- list(
      lfc_threshold = deg_results$summary$lfc_threshold,
      p_threshold = deg_results$summary$p_threshold,
      adjust_method = deg_results$summary$adjust_method,
      n_up = deg_results$summary$n_up,
      n_down = deg_results$summary$n_down
    )
  }
  
  return(params)
}


#' 创建临时报告模板
#' Create Temporary Report Template
#'
#' @return 字符串，模板路径 / Character, template path
create_temp_template <- function() {
  
  template_content <- '---
title: "`r params$title`"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: flatly
    highlight: tango
params:
  title: "Microarray Analysis Report"
  date: !r Sys.Date()
  eset: NULL
  qc_results: NULL
  deg_results: NULL
  dataset_info: NULL
  parameters: NULL
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(ggplot2)
library(DT)
```

# 数据集信息 / Dataset Information

```{r dataset-info}
if (!is.null(params$dataset_info)) {
  info <- params$dataset_info
  cat(sprintf("- 样本数 / Number of samples: %d\\n", info$n_samples))
  cat(sprintf("- 特征数 / Number of features: %d\\n", info$n_features))
  cat(sprintf("- 平台 / Platform: %s\\n", info$platform))
}
```

# 质控结果 / QC Results

```{r qc-summary}
if (!is.null(params$qc_results)) {
  qc <- params$qc_results
  cat(sprintf("- 检测到异常样本 / Outliers detected: %d\\n", 
              length(qc$outliers$all_outliers)))
  cat(sprintf("- 平均样本相关性 / Mean sample correlation: %.3f\\n", 
              qc$summary$mean_correlation))
}
```

# 差异分析结果 / Differential Expression Results

```{r deg-summary}
if (!is.null(params$deg_results)) {
  deg <- params$deg_results
  cat(sprintf("- 上调基因 / Up-regulated genes: %d\\n", deg$summary$n_up))
  cat(sprintf("- 下调基因 / Down-regulated genes: %d\\n", deg$summary$n_down))
  cat(sprintf("- 显著差异基因 / Significant DEGs: %d\\n", deg$summary$n_significant))
}
```

```{r deg-table}
if (!is.null(params$deg_results)) {
  top_degs <- head(params$deg_results$deg_table[
    params$deg_results$deg_table$DEG_status != "Not Significant", ], 50)
  if (nrow(top_degs) > 0) {
    DT::datatable(top_degs, options = list(pageLength = 10))
  }
}
```

# 分析参数 / Analysis Parameters

```{r params}
if (!is.null(params$parameters$deg$lfc_threshold)) {
  cat(sprintf("- logFC阈值 / logFC threshold: %s\\n", params$parameters$deg$lfc_threshold))
  cat(sprintf("- P值阈值 / P-value threshold: %s\\n", params$parameters$deg$p_threshold))
  cat(sprintf("- 校正方法 / Adjustment method: %s\\n", params$parameters$deg$adjust_method))
}
```

---

*Report generated by GEO-Analysis-Pipeline*
'
  
  temp_file <- tempfile(fileext = ".Rmd")
  writeLines(template_content, temp_file)
  
  return(temp_file)
}


#' 生成简化报告
#' Generate Simple Report
#'
#' @param report_data 报告数据 / Report data
#' @param output_file 输出文件 / Output file
#' @return 字符串，报告路径 / Character, report path
generate_simple_report <- function(report_data, output_file) {
  
  # 生成简单的HTML报告 / Generate simple HTML report
  html_content <- paste0('
<!DOCTYPE html>
<html>
<head>
  <title>', report_data$title, '</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 40px; }
    h1 { color: #2c3e50; }
    h2 { color: #34495e; border-bottom: 1px solid #eee; padding-bottom: 10px; }
    .info-box { background: #f8f9fa; padding: 15px; border-radius: 5px; margin: 10px 0; }
    table { border-collapse: collapse; width: 100%; }
    th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
    th { background-color: #3498db; color: white; }
  </style>
</head>
<body>
  <h1>', report_data$title, '</h1>
  <p>生成日期 / Generated: ', as.character(Sys.Date()), '</p>
  
  <h2>数据集信息 / Dataset Information</h2>
  <div class="info-box">
    <p>样本数 / Samples: ', report_data$dataset_info$n_samples, '</p>
    <p>特征数 / Features: ', report_data$dataset_info$n_features, '</p>
    <p>平台 / Platform: ', report_data$dataset_info$platform, '</p>
  </div>
')
  
  if (!is.null(report_data$qc_results)) {
    html_content <- paste0(html_content, '
  <h2>质控结果 / QC Results</h2>
  <div class="info-box">
    <p>异常样本数 / Outliers: ', length(report_data$qc_results$outliers$all_outliers), '</p>
    <p>平均相关性 / Mean correlation: ', 
      round(report_data$qc_results$summary$mean_correlation, 3), '</p>
  </div>
')
  }
  
  if (!is.null(report_data$deg_results)) {
    html_content <- paste0(html_content, '
  <h2>差异分析结果 / DEG Results</h2>
  <div class="info-box">
    <p>上调基因 / Up-regulated: ', report_data$deg_results$summary$n_up, '</p>
    <p>下调基因 / Down-regulated: ', report_data$deg_results$summary$n_down, '</p>
    <p>显著差异基因 / Total significant: ', report_data$deg_results$summary$n_significant, '</p>
  </div>
')
  }
  
  html_content <- paste0(html_content, '
  <hr>
  <p><em>Generated by GEO-Analysis-Pipeline</em></p>
</body>
</html>
')
  
  writeLines(html_content, output_file)
  message(paste0("简化报告已生成 / Simple report generated: ", output_file))
  
  return(output_file)
}


#' 生成QC报告
#' Generate QC Report
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param qc_results QC结果 / QC results
#' @param output_file 输出文件 / Output file
#'
#' @return 字符串，报告路径 / Character, report path
#' @export
generate_qc_report <- function(eset, qc_results, output_file = "qc_report.html") {
  
  return(generate_microarray_report(
    eset = eset,
    qc_results = qc_results,
    output_file = output_file,
    title = "Microarray QC Report / 微阵列质控报告"
  ))
}


#' 生成差异分析报告
#' Generate Differential Analysis Report
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param deg_results 差异分析结果 / DEG results
#' @param output_file 输出文件 / Output file
#'
#' @return 字符串，报告路径 / Character, report path
#' @export
generate_deg_report <- function(eset, deg_results, output_file = "deg_report.html") {
  
  return(generate_microarray_report(
    eset = eset,
    deg_results = deg_results,
    output_file = output_file,
    title = "Differential Expression Report / 差异表达报告"
  ))
}


#' 保存分析会话
#' Save Analysis Session
#'
#' @param ... 要保存的对象 / Objects to save
#' @param output_file 输出文件 / Output file
#'
#' @return 无返回值 / No return value
#' @export
save_analysis_session <- function(..., output_file = "analysis_session.RData") {
  
  objects <- list(...)
  
  # 添加会话信息 / Add session information
  session_info <- list(
    date = Sys.time(),
    r_version = R.version.string,
    packages = sessionInfo()$otherPkgs
  )
  
  save(objects, session_info, file = output_file)
  message(paste0("分析会话已保存 / Analysis session saved: ", output_file))
}


#' 创建分析摘要表
#' Create Analysis Summary Table
#'
#' @param qc_results QC结果 / QC results
#' @param deg_results 差异分析结果列表 / List of DEG results
#'
#' @return 数据框，分析摘要 / Data frame with analysis summary
#' @export
create_summary_table <- function(qc_results = NULL, deg_results = NULL) {
  
  summary_df <- data.frame(
    Category = character(),
    Metric = character(),
    Value = character(),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(qc_results)) {
    qc_rows <- data.frame(
      Category = rep("QC", 4),
      Metric = c("样本数 / Samples", "异常样本 / Outliers", 
                 "平均相关性 / Mean Correlation", "PC1方差 / PC1 Variance"),
      Value = c(
        as.character(qc_results$summary$n_samples),
        as.character(qc_results$summary$n_outliers),
        round(qc_results$summary$mean_correlation, 3),
        paste0(round(qc_results$summary$pca_var_pc1, 1), "%")
      ),
      stringsAsFactors = FALSE
    )
    summary_df <- rbind(summary_df, qc_rows)
  }
  
  if (!is.null(deg_results)) {
    deg_rows <- data.frame(
      Category = rep("DEG", 4),
      Metric = c("总基因数 / Total Genes", "上调基因 / Up-regulated",
                 "下调基因 / Down-regulated", "显著基因 / Significant"),
      Value = c(
        as.character(deg_results$summary$n_total),
        as.character(deg_results$summary$n_up),
        as.character(deg_results$summary$n_down),
        as.character(deg_results$summary$n_significant)
      ),
      stringsAsFactors = FALSE
    )
    summary_df <- rbind(summary_df, deg_rows)
  }
  
  return(summary_df)
}
