#' ============================================================================
#' GEO RNA-seq定量分析模块
#' RNA-seq Quantification Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - Counts数据处理 / Counts data processing
#' - 低表达基因过滤 / Low expression gene filtering
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

#' 处理RNA-seq counts数据
#' Process RNA-seq Counts Data
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param min_counts 最小counts数 / Minimum counts
#' @param min_samples 最小样本数 / Minimum number of samples
#' @param remove_zero_genes 是否移除全零基因 / Remove all-zero genes
#'
#' @return 列表，包含处理后的counts和过滤信息 / List with processed counts and filter info
#' @export
#'
#' @examples
#' \dontrun{
#' result <- process_counts(counts_matrix, min_counts = 10, min_samples = 3)
#' }
process_counts <- function(counts,
                           min_counts = 10,
                           min_samples = 2,
                           remove_zero_genes = TRUE) {
  
  message("处理counts数据... / Processing counts data...")
  
  original_genes <- nrow(counts)
  original_samples <- ncol(counts)
  
  # 记录过滤信息 / Record filter info
  filter_info <- list(
    original_genes = original_genes,
    original_samples = original_samples
  )
  
  # 1. 移除全零基因 / Remove all-zero genes
  if (remove_zero_genes) {
    zero_genes <- rowSums(counts) == 0
    counts <- counts[!zero_genes, ]
    filter_info$zero_genes_removed <- sum(zero_genes)
    message(paste0("移除全零基因 / Removed zero genes: ", sum(zero_genes)))
  }
  
  # 2. 过滤低表达基因 / Filter low expression genes
  keep <- rowSums(counts >= min_counts) >= min_samples
  counts_filtered <- counts[keep, ]
  
  filter_info$low_expr_removed <- sum(!keep)
  filter_info$final_genes <- nrow(counts_filtered)
  filter_info$final_samples <- ncol(counts_filtered)
  
  message(paste0("过滤低表达基因 / Filtered low expression genes: ", sum(!keep)))
  message(paste0("保留基因数 / Genes retained: ", nrow(counts_filtered)))
  
  return(list(
    counts = counts_filtered,
    filter_info = filter_info,
    removed_genes = rownames(counts)[!keep]
  ))
}


#' 计算CPM (Counts Per Million)
#' Calculate CPM (Counts Per Million)
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param log 是否进行log转换 / Whether to log transform
#' @param prior_count log转换时的先验值 / Prior count for log transform
#'
#' @return 矩阵，CPM值 / Matrix with CPM values
#' @export
calculate_cpm <- function(counts, log = FALSE, prior_count = 2) {
  
  lib_sizes <- colSums(counts)
  cpm <- t(t(counts) / lib_sizes * 1e6)
  
  if (log) {
    cpm <- log2(cpm + prior_count)
  }
  
  return(cpm)
}


#' 计算TPM (Transcripts Per Million)
#' Calculate TPM (Transcripts Per Million)
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param gene_lengths 基因长度向量 / Gene length vector
#'
#' @return 矩阵，TPM值 / Matrix with TPM values
#' @export
calculate_tpm <- function(counts, gene_lengths) {
  
  if (length(gene_lengths) != nrow(counts)) {
    stop("基因长度数量与counts行数不匹配 / Gene length count doesn't match counts rows")
  }
  
  # RPK (Reads Per Kilobase)
  rpk <- counts / (gene_lengths / 1000)
  
  # TPM
  scaling_factors <- colSums(rpk) / 1e6
  tpm <- t(t(rpk) / scaling_factors)
  
  return(tpm)
}


#' 计算FPKM (Fragments Per Kilobase Million)
#' Calculate FPKM (Fragments Per Kilobase Million)
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param gene_lengths 基因长度向量 / Gene length vector
#'
#' @return 矩阵，FPKM值 / Matrix with FPKM values
#' @export
calculate_fpkm <- function(counts, gene_lengths) {
  
  if (length(gene_lengths) != nrow(counts)) {
    stop("基因长度数量与counts行数不匹配 / Gene length count doesn't match counts rows")
  }
  
  lib_sizes <- colSums(counts)
  
  # FPKM = (counts * 10^9) / (gene_length * lib_size)
  fpkm <- t(t(counts) / lib_sizes) * 1e9 / gene_lengths
  
  return(fpkm)
}


#' 合并技术重复
#' Merge Technical Replicates
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param sample_groups 样本分组向量 / Sample group vector
#' @param method 合并方法 / Merge method
#'
#' @return 矩阵，合并后的counts / Merged counts matrix
#' @export
merge_technical_replicates <- function(counts, 
                                        sample_groups, 
                                        method = c("sum", "mean")) {
  
  method <- match.arg(method)
  
  unique_groups <- unique(sample_groups)
  merged <- matrix(0, nrow = nrow(counts), ncol = length(unique_groups))
  rownames(merged) <- rownames(counts)
  colnames(merged) <- unique_groups
  
  for (group in unique_groups) {
    cols <- which(sample_groups == group)
    
    if (length(cols) == 1) {
      merged[, group] <- counts[, cols]
    } else {
      if (method == "sum") {
        merged[, group] <- rowSums(counts[, cols])
      } else {
        merged[, group] <- rowMeans(counts[, cols])
      }
    }
  }
  
  message(paste0("合并技术重复：", ncol(counts), " -> ", ncol(merged), " 样本 / ",
                 "Merged replicates: ", ncol(counts), " -> ", ncol(merged), " samples"))
  
  return(merged)
}


#' 检测并移除批次效应
#' Detect and Remove Batch Effects
#'
#' @param counts Counts矩阵（已标准化）/ Normalized counts matrix
#' @param batch 批次变量 / Batch variable
#' @param covariates 协变量（可选）/ Covariates (optional)
#'
#' @return 矩阵，批次校正后的counts / Batch-corrected counts matrix
#' @export
remove_batch_effect_rnaseq <- function(counts, batch, covariates = NULL) {
  
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("需要limma包 / limma package required")
  }
  
  if (length(unique(batch)) < 2) {
    message("少于2个批次，无需校正 / Less than 2 batches, no correction needed")
    return(counts)
  }
  
  message("移除批次效应... / Removing batch effect...")
  
  # 使用limma::removeBatchEffect
  corrected <- limma::removeBatchEffect(counts, batch = batch, covariates = covariates)
  
  message("批次效应已移除 / Batch effect removed")
  
  return(corrected)
}


#' 获取高变异基因
#' Get Highly Variable Genes
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param n_top 返回top N基因 / Return top N genes
#' @param method 变异度计算方法 / Variability calculation method
#'
#' @return 字符向量，高变异基因ID / Character vector of highly variable gene IDs
#' @export
get_highly_variable_genes <- function(counts,
                                       n_top = 1000,
                                       method = c("var", "cv", "mad")) {
  
  method <- match.arg(method)
  
  # 对数转换 / Log transform
  log_counts <- log2(counts + 1)
  
  # 计算变异度 / Calculate variability
  variability <- switch(method,
    "var" = apply(log_counts, 1, var),
    "cv" = apply(log_counts, 1, sd) / (apply(log_counts, 1, mean) + 0.001),
    "mad" = apply(log_counts, 1, mad)
  )
  
  # 选择top基因 / Select top genes
  top_genes <- names(sort(variability, decreasing = TRUE))[1:min(n_top, length(variability))]
  
  return(top_genes)
}


#' 基因表达水平分类
#' Classify Gene Expression Levels
#'
#' @param counts Counts矩阵 / Counts matrix
#'
#' @return 数据框，基因表达分类 / Data frame with gene expression classification
#' @export
classify_gene_expression <- function(counts) {
  
  mean_expr <- rowMeans(counts)
  
  classification <- data.frame(
    gene = rownames(counts),
    mean_counts = mean_expr,
    level = cut(mean_expr,
                breaks = c(-Inf, 0, 1, 10, 100, 1000, Inf),
                labels = c("Zero", "Very Low", "Low", "Medium", "High", "Very High")),
    stringsAsFactors = FALSE
  )
  
  # 统计各级别数量 / Count each level
  level_summary <- table(classification$level)
  message("基因表达水平分布 / Gene expression level distribution:")
  print(level_summary)
  
  return(classification)
}
