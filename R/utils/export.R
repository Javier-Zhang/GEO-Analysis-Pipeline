#' ============================================================================
#' 数据导出模块
#' Data Export Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - CSV/TSV/Excel导出 / CSV/TSV/Excel export
#' - GCT/CLS格式（GSEA用）/ GCT/CLS format (for GSEA)
#' - 标准化表达矩阵导出 / Normalized expression matrix export
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

#' 导出表达数据
#' Export Expression Data
#'
#' @param data 表达矩阵或ExpressionSet / Expression matrix or ExpressionSet
#' @param output_file 输出文件路径 / Output file path
#' @param format 输出格式 / Output format
#' @param ... 其他参数 / Other parameters
#'
#' @return 字符串，输出文件路径 / Character, output file path
#' @export
#'
#' @examples
#' \dontrun{
#' # 导出为CSV / Export as CSV
#' export_expression(expr_matrix, "expression.csv", format = "csv")
#' 
#' # 导出为GCT格式 / Export as GCT format
#' export_expression(expr_matrix, "expression.gct", format = "gct")
#' }
export_expression <- function(data,
                               output_file,
                               format = c("csv", "tsv", "xlsx", "gct", "txt"),
                               ...) {
  
  format <- match.arg(format)
  
  # 如果是ExpressionSet，提取表达矩阵 / If ExpressionSet, extract expression matrix
  if (inherits(data, "ExpressionSet")) {
    expr_matrix <- exprs(data)
  } else {
    expr_matrix <- as.matrix(data)
  }
  
  message(paste0("导出表达数据 / Exporting expression data..."))
  message(paste0("维度 / Dimensions: ", nrow(expr_matrix), " x ", ncol(expr_matrix)))
  
  # 确保输出目录存在 / Ensure output directory exists
  output_dir <- dirname(output_file)
  if (output_dir != "." && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  switch(format,
    "csv" = {
      write.csv(expr_matrix, output_file, row.names = TRUE)
    },
    "tsv" = {
      write.table(expr_matrix, output_file, sep = "\t", 
                  row.names = TRUE, col.names = NA, quote = FALSE)
    },
    "txt" = {
      write.table(expr_matrix, output_file, sep = "\t", 
                  row.names = TRUE, col.names = NA, quote = FALSE)
    },
    "xlsx" = {
      if (requireNamespace("openxlsx", quietly = TRUE)) {
        openxlsx::write.xlsx(as.data.frame(expr_matrix), output_file, rowNames = TRUE)
      } else {
        warning("openxlsx未安装，使用CSV格式 / openxlsx not installed, using CSV")
        write.csv(expr_matrix, gsub("\\.xlsx$", ".csv", output_file), row.names = TRUE)
      }
    },
    "gct" = {
      export_gct(expr_matrix, output_file, ...)
    }
  )
  
  message(paste0("导出完成 / Export complete: ", output_file))
  
  return(output_file)
}


#' 导出GCT格式文件（用于GSEA）
#' Export GCT Format File (for GSEA)
#'
#' @param expr_matrix 表达矩阵 / Expression matrix
#' @param output_file 输出文件 / Output file
#' @param description 基因描述向量（可选）/ Gene description vector (optional)
#'
#' @return 字符串，输出文件路径 / Character, output file path
#' @export
export_gct <- function(expr_matrix, output_file, description = NULL) {
  
  # GCT文件格式 / GCT file format
  # Line 1: #1.2
  # Line 2: n_genes n_samples
  # Line 3: NAME Description sample1 sample2 ...
  # Line 4+: gene_name description expr1 expr2 ...
  
  n_genes <- nrow(expr_matrix)
  n_samples <- ncol(expr_matrix)
  
  if (is.null(description)) {
    description <- rep("na", n_genes)
  }
  
  # 确保description长度匹配 / Ensure description length matches
  if (length(description) != n_genes) {
    description <- rep("na", n_genes)
  }
  
  # 打开文件 / Open file
  con <- file(output_file, "w")
  
  # 写入头部 / Write header
  writeLines("#1.2", con)
  writeLines(paste(n_genes, n_samples, sep = "\t"), con)
  writeLines(paste(c("NAME", "Description", colnames(expr_matrix)), collapse = "\t"), con)
  
  # 写入数据 / Write data
  for (i in seq_len(n_genes)) {
    line <- paste(c(rownames(expr_matrix)[i], description[i], expr_matrix[i, ]), collapse = "\t")
    writeLines(line, con)
  }
  
  close(con)
  
  message(paste0("GCT文件已导出 / GCT file exported: ", output_file))
  
  return(output_file)
}


#' 导出CLS格式文件（用于GSEA）
#' Export CLS Format File (for GSEA)
#'
#' @param group 分组向量 / Group vector
#' @param output_file 输出文件 / Output file
#'
#' @return 字符串，输出文件路径 / Character, output file path
#' @export
export_cls <- function(group, output_file) {
  
  # CLS文件格式 / CLS file format
  # Line 1: n_samples n_classes 1
  # Line 2: # class1 class2 ...
  # Line 3: group1 group2 group3 ...
  
  group <- as.factor(group)
  n_samples <- length(group)
  classes <- levels(group)
  n_classes <- length(classes)
  
  # 打开文件 / Open file
  con <- file(output_file, "w")
  
  # 写入头部 / Write header
  writeLines(paste(n_samples, n_classes, "1", sep = " "), con)
  writeLines(paste("#", paste(classes, collapse = " "), sep = " "), con)
  writeLines(paste(as.character(group), collapse = " "), con)
  
  close(con)
  
  message(paste0("CLS文件已导出 / CLS file exported: ", output_file))
  
  return(output_file)
}


#' 导出差异表达结果
#' Export Differential Expression Results
#'
#' @param deg_result 差异分析结果 / Differential analysis result
#' @param output_file 输出文件 / Output file
#' @param format 输出格式 / Output format
#' @param top_n 导出前N个基因（可选）/ Export top N genes (optional)
#'
#' @return 字符串，输出文件路径 / Character, output file path
#' @export
export_deg_results <- function(deg_result,
                                output_file,
                                format = c("csv", "tsv", "xlsx"),
                                top_n = NULL) {
  
  format <- match.arg(format)
  
  # 获取结果表 / Get result table
  if (is.list(deg_result) && "deg_table" %in% names(deg_result)) {
    results <- deg_result$deg_table
  } else {
    results <- deg_result
  }
  
  # 筛选top N / Filter top N
  if (!is.null(top_n) && top_n < nrow(results)) {
    results <- head(results, top_n)
  }
  
  # 导出 / Export
  switch(format,
    "csv" = {
      write.csv(results, output_file, row.names = TRUE)
    },
    "tsv" = {
      write.table(results, output_file, sep = "\t", 
                  row.names = TRUE, col.names = NA, quote = FALSE)
    },
    "xlsx" = {
      if (requireNamespace("openxlsx", quietly = TRUE)) {
        openxlsx::write.xlsx(results, output_file, rowNames = TRUE)
      } else {
        write.csv(results, gsub("\\.xlsx$", ".csv", output_file), row.names = TRUE)
      }
    }
  )
  
  message(paste0("差异分析结果已导出 / DEG results exported: ", output_file))
  
  return(output_file)
}


#' 批量导出分析结果
#' Batch Export Analysis Results
#'
#' @param results_list 结果列表 / List of results
#' @param output_dir 输出目录 / Output directory
#' @param prefix 文件前缀 / File prefix
#' @param format 输出格式 / Output format
#'
#' @return 字符向量，导出的文件路径 / Character vector of exported file paths
#' @export
batch_export <- function(results_list,
                          output_dir = "results",
                          prefix = "",
                          format = "csv") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  exported_files <- character()
  
  for (name in names(results_list)) {
    result <- results_list[[name]]
    
    if (is.null(result)) next
    
    file_name <- paste0(prefix, ifelse(prefix == "", "", "_"), name, ".", format)
    file_path <- file.path(output_dir, file_name)
    
    tryCatch({
      if (is.matrix(result) || is.data.frame(result)) {
        export_expression(result, file_path, format = format)
        exported_files <- c(exported_files, file_path)
      } else if (is.list(result) && "deg_table" %in% names(result)) {
        export_deg_results(result, file_path, format = format)
        exported_files <- c(exported_files, file_path)
      }
    }, error = function(e) {
      warning(paste0("导出 ", name, " 失败 / Export failed: ", conditionMessage(e)))
    })
  }
  
  message(paste0("批量导出完成，共 ", length(exported_files), " 个文件 / ",
                 "Batch export complete, ", length(exported_files), " files"))
  
  return(exported_files)
}


#' 导出为R对象
#' Export as R Object
#'
#' @param ... 要保存的对象 / Objects to save
#' @param output_file 输出文件 / Output file
#'
#' @return 字符串，输出文件路径 / Character, output file path
#' @export
export_rdata <- function(..., output_file = "results.RData") {
  
  save(..., file = output_file)
  message(paste0("R对象已保存 / R objects saved: ", output_file))
  
  return(output_file)
}


#' 导出为RDS对象
#' Export as RDS Object
#'
#' @param object 要保存的对象 / Object to save
#' @param output_file 输出文件 / Output file
#'
#' @return 字符串，输出文件路径 / Character, output file path
#' @export
export_rds <- function(object, output_file = "results.rds") {
  
  saveRDS(object, output_file)
  message(paste0("RDS对象已保存 / RDS object saved: ", output_file))
  
  return(output_file)
}


#' 创建结果汇总文件
#' Create Results Summary File
#'
#' @param qc_results QC结果 / QC results
#' @param deg_results 差异分析结果 / DEG results
#' @param output_file 输出文件 / Output file
#'
#' @return 字符串，输出文件路径 / Character, output file path
#' @export
create_results_summary <- function(qc_results = NULL,
                                    deg_results = NULL,
                                    output_file = "analysis_summary.txt") {
  
  lines <- c(
    "GEO数据分析结果汇总 / GEO Data Analysis Summary",
    paste0("生成时间 / Generated: ", Sys.time()),
    "=" %>% rep(50) %>% paste(collapse = ""),
    ""
  )
  
  # QC结果 / QC results
  if (!is.null(qc_results)) {
    lines <- c(lines,
      "质控结果 / QC Results",
      "-" %>% rep(30) %>% paste(collapse = ""),
      paste0("样本数 / Samples: ", qc_results$summary$n_samples),
      paste0("异常样本 / Outliers: ", qc_results$summary$n_outliers),
      paste0("平均相关性 / Mean Correlation: ", round(qc_results$summary$mean_correlation, 3)),
      ""
    )
  }
  
  # 差异分析结果 / DEG results
  if (!is.null(deg_results)) {
    lines <- c(lines,
      "差异分析结果 / Differential Expression Results",
      "-" %>% rep(30) %>% paste(collapse = ""),
      paste0("总基因数 / Total Genes: ", deg_results$summary$n_total),
      paste0("上调基因 / Up-regulated: ", deg_results$summary$n_up),
      paste0("下调基因 / Down-regulated: ", deg_results$summary$n_down),
      paste0("logFC阈值 / logFC Threshold: ", deg_results$summary$lfc_threshold),
      paste0("P值阈值 / P-value Threshold: ", deg_results$summary$p_threshold),
      ""
    )
  }
  
  writeLines(lines, output_file)
  message(paste0("汇总已保存 / Summary saved: ", output_file))
  
  return(output_file)
}
