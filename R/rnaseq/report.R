#' ============================================================================
#' GEO RNA-seq分析报告生成模块
#' RNA-seq Analysis Report Generation Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 生成HTML报告 / Generate HTML report
#' - 整合分析结果 / Integrate analysis results
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

suppressPackageStartupMessages({
  library(rmarkdown)
})

#' 生成RNA-seq分析报告
#' Generate RNA-seq Analysis Report
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param qc_results QC结果 / QC results
#' @param deg_results 差异分析结果 / DEG results
#' @param output_file 输出文件 / Output file
#' @param title 报告标题 / Report title
#'
#' @return 字符串，报告路径 / Character, report path
#' @export
generate_rnaseq_report <- function(counts,
                                    qc_results = NULL,
                                    deg_results = NULL,
                                    output_file = "rnaseq_report.html",
                                    title = "RNA-seq Analysis Report / RNA-seq分析报告") {
  
  message("生成RNA-seq报告... / Generating RNA-seq report...")
  
  # 确保输出目录存在 / Ensure output directory exists
  output_dir <- dirname(output_file)
  if (output_dir != "." && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 准备数据 / Prepare data
  report_data <- list(
    title = title,
    date = Sys.Date(),
    dataset_info = get_rnaseq_dataset_info(counts),
    qc_results = qc_results,
    deg_results = deg_results
  )
  
  # 生成简单HTML报告 / Generate simple HTML report
  html_content <- generate_rnaseq_html(report_data)
  
  writeLines(html_content, output_file)
  
  message(paste0("报告已生成 / Report generated: ", output_file))
  
  return(output_file)
}


#' 获取RNA-seq数据集信息
#' Get RNA-seq Dataset Information
#'
#' @param counts Counts矩阵 / Counts matrix
#' @return 列表，数据集信息 / List with dataset info
get_rnaseq_dataset_info <- function(counts) {
  
  list(
    n_genes = nrow(counts),
    n_samples = ncol(counts),
    total_counts = sum(counts),
    median_lib_size = median(colSums(counts)),
    sample_names = colnames(counts)
  )
}


#' 生成RNA-seq HTML报告内容
#' Generate RNA-seq HTML Report Content
#'
#' @param report_data 报告数据 / Report data
#' @return 字符串，HTML内容 / Character, HTML content
generate_rnaseq_html <- function(report_data) {
  
  html <- paste0('
<!DOCTYPE html>
<html>
<head>
  <title>', report_data$title, '</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 40px; max-width: 1200px; }
    h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
    h2 { color: #34495e; margin-top: 30px; }
    .info-box { background: #f8f9fa; padding: 20px; border-radius: 8px; margin: 15px 0; }
    .stat { display: inline-block; margin: 10px 20px; text-align: center; }
    .stat-value { font-size: 24px; font-weight: bold; color: #3498db; }
    .stat-label { font-size: 14px; color: #666; }
    table { border-collapse: collapse; width: 100%; margin: 15px 0; }
    th, td { border: 1px solid #ddd; padding: 10px; text-align: left; }
    th { background-color: #3498db; color: white; }
    tr:nth-child(even) { background-color: #f9f9f9; }
    .success { color: #27ae60; }
    .warning { color: #f39c12; }
  </style>
</head>
<body>
  <h1>', report_data$title, '</h1>
  <p>生成日期 / Generated: ', as.character(report_data$date), '</p>
  
  <h2>数据集概览 / Dataset Overview</h2>
  <div class="info-box">
    <div class="stat">
      <div class="stat-value">', format(report_data$dataset_info$n_genes, big.mark = ","), '</div>
      <div class="stat-label">基因数 / Genes</div>
    </div>
    <div class="stat">
      <div class="stat-value">', report_data$dataset_info$n_samples, '</div>
      <div class="stat-label">样本数 / Samples</div>
    </div>
    <div class="stat">
      <div class="stat-value">', format(round(report_data$dataset_info$median_lib_size / 1e6, 1)), 'M</div>
      <div class="stat-label">中位文库大小 / Median Library Size</div>
    </div>
  </div>
')
  
  # QC结果部分 / QC results section
  if (!is.null(report_data$qc_results)) {
    html <- paste0(html, '
  <h2>质控结果 / QC Results</h2>
  <div class="info-box">
    <div class="stat">
      <div class="stat-value">', report_data$qc_results$summary$n_outliers, '</div>
      <div class="stat-label">异常样本 / Outliers</div>
    </div>
    <div class="stat">
      <div class="stat-value">', round(report_data$qc_results$summary$median_detection_rate, 1), '%</div>
      <div class="stat-label">中位检测率 / Median Detection Rate</div>
    </div>
  </div>
')
  }
  
  # 差异分析结果部分 / DEG results section
  if (!is.null(report_data$deg_results)) {
    html <- paste0(html, '
  <h2>差异表达分析结果 / Differential Expression Results</h2>
  <div class="info-box">
    <div class="stat">
      <div class="stat-value success">', report_data$deg_results$summary$n_up, '</div>
      <div class="stat-label">上调基因 / Up-regulated</div>
    </div>
    <div class="stat">
      <div class="stat-value warning">', report_data$deg_results$summary$n_down, '</div>
      <div class="stat-label">下调基因 / Down-regulated</div>
    </div>
    <div class="stat">
      <div class="stat-value">', report_data$deg_results$summary$n_significant, '</div>
      <div class="stat-label">显著基因总数 / Total Significant</div>
    </div>
  </div>
  
  <h3>分析参数 / Analysis Parameters</h3>
  <table>
    <tr><th>参数 / Parameter</th><th>值 / Value</th></tr>
    <tr><td>分析方法 / Method</td><td>', report_data$deg_results$summary$method, '</td></tr>
    <tr><td>logFC阈值 / logFC Threshold</td><td>', report_data$deg_results$summary$lfc_threshold, '</td></tr>
    <tr><td>P值阈值 / P-value Threshold</td><td>', report_data$deg_results$summary$p_threshold, '</td></tr>
  </table>
')
    
    # Top差异基因表格 / Top DEG table
    top_degs <- head(report_data$deg_results$deg_table[
      report_data$deg_results$deg_table$DEG_status != "Not Significant", ], 20)
    
    if (nrow(top_degs) > 0) {
      html <- paste0(html, '
  <h3>Top差异基因 / Top Differentially Expressed Genes</h3>
  <table>
    <tr>
      <th>基因 / Gene</th>
      <th>log2FC</th>
      <th>adj. P-value</th>
      <th>状态 / Status</th>
    </tr>
')
      for (i in 1:nrow(top_degs)) {
        html <- paste0(html, '
    <tr>
      <td>', top_degs$gene[i], '</td>
      <td>', round(top_degs$log2FoldChange[i], 3), '</td>
      <td>', format(top_degs$padj[i], scientific = TRUE, digits = 3), '</td>
      <td>', top_degs$DEG_status[i], '</td>
    </tr>
')
      }
      html <- paste0(html, '
  </table>
')
    }
  }
  
  html <- paste0(html, '
  <hr>
  <p><em>Generated by GEO-Analysis-Pipeline</em></p>
</body>
</html>
')
  
  return(html)
}


#' 生成RNA-seq QC报告
#' Generate RNA-seq QC Report
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param qc_results QC结果 / QC results
#' @param output_file 输出文件 / Output file
#'
#' @return 字符串，报告路径 / Character, report path
#' @export
generate_rnaseq_qc_report <- function(counts, qc_results, output_file = "rnaseq_qc_report.html") {
  
  return(generate_rnaseq_report(
    counts = counts,
    qc_results = qc_results,
    output_file = output_file,
    title = "RNA-seq QC Report / RNA-seq质控报告"
  ))
}


#' 保存RNA-seq分析结果
#' Save RNA-seq Analysis Results
#'
#' @param results 分析结果 / Analysis results
#' @param output_dir 输出目录 / Output directory
#' @param prefix 文件前缀 / File prefix
#'
#' @return 无返回值 / No return value
#' @export
save_rnaseq_results <- function(results, output_dir = "results", prefix = "rnaseq") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 保存为RDS / Save as RDS
  rds_file <- file.path(output_dir, paste0(prefix, "_results.rds"))
  saveRDS(results, rds_file)
  
  # 如果有差异表达结果，保存为CSV / If DEG results, save as CSV
  if (!is.null(results$deg_table)) {
    csv_file <- file.path(output_dir, paste0(prefix, "_deg_table.csv"))
    write.csv(results$deg_table, csv_file, row.names = FALSE)
    message(paste0("差异基因表已保存 / DEG table saved: ", csv_file))
  }
  
  message(paste0("结果已保存 / Results saved to: ", output_dir))
}
