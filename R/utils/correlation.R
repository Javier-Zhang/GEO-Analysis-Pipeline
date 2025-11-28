#' ============================================================================
#' 原始数据与Series Matrix相关性分析模块
#' Raw Data vs Series Matrix Correlation Analysis Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 原始数据vs series matrix相关性 / Raw data vs series matrix correlation
#' - 样本一致性验证 / Sample consistency validation
#' - 相关性可视化 / Correlation visualization
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

#' 比较原始数据和Series Matrix
#' Compare Raw Data and Series Matrix
#'
#' @param raw_data 原始表达数据（矩阵）/ Raw expression data (matrix)
#' @param matrix_data Series Matrix数据（矩阵）/ Series Matrix data (matrix)
#' @param method 相关性方法 / Correlation method
#'
#' @return 列表，比较结果 / List with comparison results
#' @export
#'
#' @examples
#' \dontrun{
#' result <- compare_raw_matrix(raw_expr, matrix_expr)
#' }
compare_raw_matrix <- function(raw_data,
                                matrix_data,
                                method = c("pearson", "spearman")) {
  
  method <- match.arg(method)
  
  message("比较原始数据和Series Matrix... / Comparing raw data and Series Matrix...")
  
  # 确保样本名匹配 / Ensure sample names match
  common_samples <- intersect(colnames(raw_data), colnames(matrix_data))
  
  if (length(common_samples) == 0) {
    warning("没有共同样本 / No common samples")
    return(NULL)
  }
  
  message(paste0("共同样本数 / Common samples: ", length(common_samples)))
  
  raw_data <- raw_data[, common_samples]
  matrix_data <- matrix_data[, common_samples]
  
  # 确保探针名匹配 / Ensure probe names match
  common_probes <- intersect(rownames(raw_data), rownames(matrix_data))
  
  if (length(common_probes) == 0) {
    warning("没有共同探针 / No common probes")
    return(NULL)
  }
  
  message(paste0("共同探针数 / Common probes: ", length(common_probes)))
  
  raw_data <- raw_data[common_probes, ]
  matrix_data <- matrix_data[common_probes, ]
  
  # 计算样本间相关性 / Calculate sample correlations
  sample_correlations <- calculate_sample_correlations(raw_data, matrix_data, method)
  
  # 计算总体相关性 / Calculate overall correlation
  overall_cor <- cor(as.vector(raw_data), as.vector(matrix_data), 
                     method = method, use = "pairwise.complete.obs")
  
  result <- list(
    n_common_samples = length(common_samples),
    n_common_probes = length(common_probes),
    sample_correlations = sample_correlations,
    overall_correlation = overall_cor,
    method = method,
    consistency_pass = all(sample_correlations > 0.9)
  )
  
  message(paste0("总体相关性 / Overall correlation: ", round(overall_cor, 4)))
  
  return(result)
}


#' 计算样本间相关性
#' Calculate Sample Correlations
#'
#' @param raw_data 原始数据 / Raw data
#' @param matrix_data Matrix数据 / Matrix data
#' @param method 相关性方法 / Correlation method
#' @return 数值向量，每个样本的相关性 / Numeric vector of correlations per sample
calculate_sample_correlations <- function(raw_data, matrix_data, method) {
  
  correlations <- sapply(colnames(raw_data), function(sample) {
    cor(raw_data[, sample], matrix_data[, sample], 
        method = method, use = "pairwise.complete.obs")
  })
  
  return(correlations)
}


#' 验证样本顺序一致性
#' Validate Sample Order Consistency
#'
#' @param data1 数据集1 / Dataset 1
#' @param data2 数据集2 / Dataset 2
#'
#' @return 逻辑值，TRUE表示一致 / Logical, TRUE if consistent
#' @export
validate_sample_order <- function(data1, data2) {
  
  samples1 <- colnames(data1)
  samples2 <- colnames(data2)
  
  if (length(samples1) != length(samples2)) {
    message("样本数不同 / Different number of samples")
    return(FALSE)
  }
  
  if (!all(samples1 == samples2)) {
    message("样本顺序不同 / Sample order differs")
    return(FALSE)
  }
  
  message("样本顺序一致 / Sample order consistent")
  return(TRUE)
}


#' 可视化相关性比较
#' Visualize Correlation Comparison
#'
#' @param comparison_result 比较结果 / Comparison result
#' @param output_dir 输出目录 / Output directory
#'
#' @return 列表，图文件路径 / List of plot file paths
#' @export
plot_correlation_comparison <- function(comparison_result, output_dir = NULL) {
  
  if (is.null(comparison_result)) {
    return(NULL)
  }
  
  plots <- list()
  
  # 1. 样本相关性条形图 / Sample correlation bar plot
  tryCatch({
    cor_df <- data.frame(
      Sample = names(comparison_result$sample_correlations),
      Correlation = comparison_result$sample_correlations
    )
    cor_df$Sample <- factor(cor_df$Sample, 
                            levels = cor_df$Sample[order(cor_df$Correlation)])
    
    p <- ggplot(cor_df, aes(x = Sample, y = Correlation)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
      coord_flip() +
      labs(title = "Sample Correlation: Raw vs Matrix / 原始数据vs矩阵相关性",
           x = "Sample / 样本",
           y = paste0(comparison_result$method, " correlation / 相关性")) +
      theme_minimal()
    
    plots$sample_correlation <- p
    
    if (!is.null(output_dir)) {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      ggsave(file.path(output_dir, "sample_correlation.png"), p, 
             width = 10, height = max(6, length(cor_df$Sample) * 0.2))
    }
  }, error = function(e) message("样本相关性图失败 / Sample correlation plot failed"))
  
  # 2. 相关性分布直方图 / Correlation distribution histogram
  tryCatch({
    cor_df <- data.frame(Correlation = comparison_result$sample_correlations)
    
    p <- ggplot(cor_df, aes(x = Correlation)) +
      geom_histogram(bins = 20, fill = "steelblue", color = "white") +
      geom_vline(xintercept = 0.9, linetype = "dashed", color = "red") +
      labs(title = "Correlation Distribution / 相关性分布",
           x = "Correlation / 相关性",
           y = "Count / 计数") +
      theme_minimal()
    
    plots$correlation_distribution <- p
    
    if (!is.null(output_dir)) {
      ggsave(file.path(output_dir, "correlation_distribution.png"), p, 
             width = 8, height = 6)
    }
  }, error = function(e) message("相关性分布图失败 / Correlation distribution plot failed"))
  
  return(plots)
}


#' 检测数据转换方式
#' Detect Data Transformation
#'
#' @param data1 数据集1 / Dataset 1
#' @param data2 数据集2 / Dataset 2
#'
#' @return 字符串，推测的转换方式 / Character, inferred transformation
#' @export
detect_transformation <- function(data1, data2) {
  
  # 比较数值范围 / Compare value ranges
  range1 <- range(data1, na.rm = TRUE)
  range2 <- range(data2, na.rm = TRUE)
  
  message(paste0("数据集1范围 / Dataset 1 range: ", round(range1[1], 2), " - ", round(range1[2], 2)))
  message(paste0("数据集2范围 / Dataset 2 range: ", round(range2[1], 2), " - ", round(range2[2], 2)))
  
  # 判断转换类型 / Determine transformation type
  if (range1[2] > 1000 && range2[2] < 30) {
    return("log2")
  } else if (range1[2] < 30 && range2[2] > 1000) {
    return("exp")
  } else if (abs(range1[2] - range2[2]) < range1[2] * 0.1) {
    return("none")
  } else {
    return("unknown")
  }
}


#' 计算表达矩阵相关性
#' Calculate Expression Matrix Correlation
#'
#' @param matrix1 表达矩阵1 / Expression matrix 1
#' @param matrix2 表达矩阵2 / Expression matrix 2
#' @param method 方法 / Method
#'
#' @return 数据框，相关性结果 / Data frame with correlation results
#' @export
correlate_matrices <- function(matrix1, matrix2, method = "pearson") {
  
  # 确保维度匹配 / Ensure dimensions match
  common_rows <- intersect(rownames(matrix1), rownames(matrix2))
  common_cols <- intersect(colnames(matrix1), colnames(matrix2))
  
  if (length(common_rows) == 0 || length(common_cols) == 0) {
    stop("没有共同的行或列 / No common rows or columns")
  }
  
  m1 <- matrix1[common_rows, common_cols]
  m2 <- matrix2[common_rows, common_cols]
  
  # 按行计算相关性（基因级）/ Calculate row-wise correlation (gene level)
  gene_cors <- sapply(seq_len(nrow(m1)), function(i) {
    cor(m1[i, ], m2[i, ], method = method, use = "pairwise.complete.obs")
  })
  names(gene_cors) <- common_rows
  
  # 按列计算相关性（样本级）/ Calculate column-wise correlation (sample level)
  sample_cors <- sapply(seq_len(ncol(m1)), function(i) {
    cor(m1[, i], m2[, i], method = method, use = "pairwise.complete.obs")
  })
  names(sample_cors) <- common_cols
  
  return(list(
    gene_correlations = gene_cors,
    sample_correlations = sample_cors,
    n_genes = length(common_rows),
    n_samples = length(common_cols),
    median_gene_cor = median(gene_cors, na.rm = TRUE),
    median_sample_cor = median(sample_cors, na.rm = TRUE)
  ))
}


#' 生成相关性报告
#' Generate Correlation Report
#'
#' @param comparison_result 比较结果 / Comparison result
#'
#' @return 字符串，报告文本 / Character, report text
#' @export
generate_correlation_report <- function(comparison_result) {
  
  if (is.null(comparison_result)) {
    return("无法生成报告：比较结果为空 / Cannot generate report: comparison result is NULL")
  }
  
  report <- paste0(
    "数据一致性检验报告 / Data Consistency Validation Report\n",
    "================================================\n\n",
    "共同样本数 / Common samples: ", comparison_result$n_common_samples, "\n",
    "共同探针数 / Common probes: ", comparison_result$n_common_probes, "\n",
    "相关性方法 / Correlation method: ", comparison_result$method, "\n",
    "总体相关性 / Overall correlation: ", round(comparison_result$overall_correlation, 4), "\n\n",
    "样本相关性统计 / Sample correlation statistics:\n",
    "  最小值 / Min: ", round(min(comparison_result$sample_correlations), 4), "\n",
    "  中位数 / Median: ", round(median(comparison_result$sample_correlations), 4), "\n",
    "  最大值 / Max: ", round(max(comparison_result$sample_correlations), 4), "\n\n",
    "一致性检验 / Consistency check: ", 
    if (comparison_result$consistency_pass) "通过 / PASS" else "未通过 / FAIL", "\n"
  )
  
  return(report)
}
