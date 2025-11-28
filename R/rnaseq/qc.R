#' ============================================================================
#' GEO RNA-seq数据质控模块
#' RNA-seq Data Quality Control Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 读数分布统计 / Read distribution statistics
#' - 样本相关性和PCA / Sample correlation and PCA
#' - 基因检测率 / Gene detection rate
#' - 测序深度评估 / Sequencing depth assessment
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

#' 执行RNA-seq质控分析
#' Perform RNA-seq QC Analysis
#'
#' @param counts Counts矩阵（基因 x 样本）/ Counts matrix (genes x samples)
#' @param phenotype 表型数据框 / Phenotype data frame
#' @param output_dir 输出目录 / Output directory
#'
#' @return 列表，包含QC结果 / List containing QC results
#' @export
#'
#' @examples
#' \dontrun{
#' qc_results <- run_rnaseq_qc(counts_matrix, phenotype_data, output_dir = "qc")
#' }
run_rnaseq_qc <- function(counts, phenotype = NULL, output_dir = "qc_results") {
  
  message("开始RNA-seq质控分析... / Starting RNA-seq QC analysis...")
  
  # 创建输出目录 / Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  qc_results <- list(
    library_stats = NULL,
    gene_stats = NULL,
    correlation_matrix = NULL,
    pca_results = NULL,
    outliers = NULL,
    plots = list()
  )
  
  # 1. 文库统计 / Library statistics
  message("[1/5] 计算文库统计... / Calculating library statistics...")
  qc_results$library_stats <- calculate_library_stats(counts)
  
  # 2. 基因统计 / Gene statistics
  message("[2/5] 计算基因统计... / Calculating gene statistics...")
  qc_results$gene_stats <- calculate_gene_stats(counts)
  
  # 3. 样本相关性 / Sample correlation
  message("[3/5] 计算样本相关性... / Calculating sample correlation...")
  qc_results$correlation_matrix <- calculate_sample_correlation_rnaseq(counts)
  
  # 4. PCA分析 / PCA analysis
  message("[4/5] 执行PCA分析... / Performing PCA analysis...")
  qc_results$pca_results <- perform_rnaseq_pca(counts)
  
  # 5. 异常样本检测 / Outlier detection
  message("[5/5] 检测异常样本... / Detecting outlier samples...")
  qc_results$outliers <- detect_rnaseq_outliers(counts, qc_results$correlation_matrix)
  
  # 生成QC图 / Generate QC plots
  qc_results$plots <- generate_rnaseq_qc_plots(counts, qc_results, phenotype, output_dir)
  
  # 生成摘要 / Generate summary
  qc_results$summary <- generate_rnaseq_qc_summary(qc_results)
  
  message("RNA-seq质控完成！/ RNA-seq QC complete!")
  
  return(qc_results)
}


#' 计算文库统计
#' Calculate Library Statistics
#'
#' @param counts Counts矩阵 / Counts matrix
#' @return 数据框，文库统计 / Data frame with library statistics
calculate_library_stats <- function(counts) {
  
  stats <- data.frame(
    sample_id = colnames(counts),
    total_counts = colSums(counts),
    detected_genes = colSums(counts > 0),
    mean_counts_per_gene = colMeans(counts),
    median_counts_per_gene = apply(counts, 2, median),
    max_counts = apply(counts, 2, max),
    zero_count_genes = colSums(counts == 0),
    stringsAsFactors = FALSE
  )
  
  # 计算文库大小比例 / Calculate library size ratio
  stats$lib_size_ratio <- stats$total_counts / median(stats$total_counts)
  
  # 检测率 / Detection rate
  stats$detection_rate <- stats$detected_genes / nrow(counts) * 100
  
  return(stats)
}


#' 计算基因统计
#' Calculate Gene Statistics
#'
#' @param counts Counts矩阵 / Counts matrix
#' @return 数据框，基因统计 / Data frame with gene statistics
calculate_gene_stats <- function(counts) {
  
  stats <- data.frame(
    gene_id = rownames(counts),
    mean_counts = rowMeans(counts),
    median_counts = apply(counts, 1, median),
    sd_counts = apply(counts, 1, sd),
    max_counts = apply(counts, 1, max),
    n_samples_detected = rowSums(counts > 0),
    detection_rate = rowSums(counts > 0) / ncol(counts) * 100,
    stringsAsFactors = FALSE
  )
  
  # 分类基因表达水平 / Classify gene expression levels
  stats$expression_level <- cut(stats$mean_counts,
                                breaks = c(-Inf, 0, 1, 10, 100, 1000, Inf),
                                labels = c("Zero", "Very Low", "Low", "Medium", "High", "Very High"))
  
  return(stats)
}


#' 计算RNA-seq样本相关性
#' Calculate RNA-seq Sample Correlation
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param method 相关性方法 / Correlation method
#' @return 相关性矩阵 / Correlation matrix
calculate_sample_correlation_rnaseq <- function(counts, method = "spearman") {
  
  # 对数转换（+1防止log(0)）/ Log transform (+1 to avoid log(0))
  log_counts <- log2(counts + 1)
  
  # 过滤低表达基因 / Filter low expression genes
  keep <- rowMeans(counts) > 1
  log_counts_filtered <- log_counts[keep, ]
  
  cor_matrix <- cor(log_counts_filtered, method = method, use = "pairwise.complete.obs")
  
  return(cor_matrix)
}


#' 执行RNA-seq PCA分析
#' Perform RNA-seq PCA Analysis
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param n_top 使用表达最高的基因数 / Number of top expressed genes to use
#' @return 列表，PCA结果 / List with PCA results
perform_rnaseq_pca <- function(counts, n_top = 1000) {
  
  # 对数转换 / Log transform
  log_counts <- log2(counts + 1)
  
  # 选择高变异基因 / Select high variance genes
  vars <- apply(log_counts, 1, var)
  top_genes <- names(sort(vars, decreasing = TRUE)[1:min(n_top, length(vars))])
  
  log_counts_subset <- log_counts[top_genes, ]
  
  # PCA
  pca <- prcomp(t(log_counts_subset), center = TRUE, scale. = TRUE)
  
  var_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100
  
  result <- list(
    pca = pca,
    scores = pca$x,
    var_explained = var_explained,
    cumulative_var = cumsum(var_explained),
    top_genes = top_genes
  )
  
  return(result)
}


#' 检测RNA-seq异常样本
#' Detect RNA-seq Outlier Samples
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param cor_matrix 相关性矩阵 / Correlation matrix
#' @return 列表，异常样本信息 / List with outlier information
detect_rnaseq_outliers <- function(counts, cor_matrix) {
  
  # 方法1：基于相关性 / Method 1: Based on correlation
  mean_cors <- colMeans(cor_matrix)
  cor_threshold <- median(mean_cors) - 2 * mad(mean_cors)
  cor_outliers <- names(mean_cors)[mean_cors < cor_threshold]
  
  # 方法2：基于文库大小 / Method 2: Based on library size
  lib_sizes <- colSums(counts)
  lib_median <- median(lib_sizes)
  lib_mad <- mad(lib_sizes)
  lib_outliers <- colnames(counts)[abs(lib_sizes - lib_median) > 3 * lib_mad]
  
  # 方法3：基于检测基因数 / Method 3: Based on detected genes
  detected <- colSums(counts > 0)
  det_median <- median(detected)
  det_mad <- mad(detected)
  det_outliers <- colnames(counts)[abs(detected - det_median) > 3 * det_mad]
  
  result <- list(
    correlation_outliers = cor_outliers,
    library_size_outliers = lib_outliers,
    detection_outliers = det_outliers,
    all_outliers = unique(c(cor_outliers, lib_outliers, det_outliers)),
    mean_correlations = mean_cors,
    library_sizes = lib_sizes,
    detected_genes = detected
  )
  
  return(result)
}


#' 生成RNA-seq QC图
#' Generate RNA-seq QC Plots
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param qc_results QC结果 / QC results
#' @param phenotype 表型数据 / Phenotype data
#' @param output_dir 输出目录 / Output directory
#' @return 列表，图路径 / List of plot paths
generate_rnaseq_qc_plots <- function(counts, qc_results, phenotype, output_dir) {
  
  plots <- list()
  
  # 1. 文库大小图 / Library size plot
  tryCatch({
    lib_file <- file.path(output_dir, "library_size.png")
    png(lib_file, width = 1000, height = 600, res = 100)
    
    lib_df <- qc_results$library_stats
    lib_df$sample_id <- factor(lib_df$sample_id, levels = lib_df$sample_id[order(lib_df$total_counts)])
    
    p <- ggplot(lib_df, aes(x = sample_id, y = total_counts / 1e6)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(title = "Library Size / 文库大小",
           x = "Sample / 样本",
           y = "Total Counts (millions) / 总读数（百万）") +
      theme_minimal()
    
    print(p)
    dev.off()
    plots$library_size <- lib_file
  }, error = function(e) message("文库大小图失败 / Library size plot failed"))
  
  # 2. 基因检测率图 / Gene detection rate plot
  tryCatch({
    det_file <- file.path(output_dir, "detection_rate.png")
    png(det_file, width = 1000, height = 600, res = 100)
    
    det_df <- qc_results$library_stats
    det_df$sample_id <- factor(det_df$sample_id, levels = det_df$sample_id[order(det_df$detection_rate)])
    
    p <- ggplot(det_df, aes(x = sample_id, y = detection_rate)) +
      geom_bar(stat = "identity", fill = "darkgreen") +
      coord_flip() +
      labs(title = "Gene Detection Rate / 基因检测率",
           x = "Sample / 样本",
           y = "Detection Rate (%) / 检测率（%）") +
      theme_minimal()
    
    print(p)
    dev.off()
    plots$detection_rate <- det_file
  }, error = function(e) message("检测率图失败 / Detection rate plot failed"))
  
  # 3. 相关性热图 / Correlation heatmap
  tryCatch({
    cor_file <- file.path(output_dir, "correlation_heatmap.png")
    png(cor_file, width = 1000, height = 1000, res = 100)
    
    pheatmap(qc_results$correlation_matrix,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             main = "Sample Correlation / 样本相关性",
             show_rownames = ncol(counts) <= 50,
             show_colnames = ncol(counts) <= 50)
    
    dev.off()
    plots$correlation <- cor_file
  }, error = function(e) message("相关性热图失败 / Correlation heatmap failed"))
  
  # 4. PCA图 / PCA plot
  tryCatch({
    pca_file <- file.path(output_dir, "pca_plot.png")
    png(pca_file, width = 1000, height = 800, res = 100)
    
    pca_df <- data.frame(
      PC1 = qc_results$pca_results$scores[, 1],
      PC2 = qc_results$pca_results$scores[, 2],
      Sample = rownames(qc_results$pca_results$scores)
    )
    
    var_exp <- qc_results$pca_results$var_explained
    
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
      geom_point(size = 3, color = "steelblue") +
      theme_minimal() +
      labs(title = "PCA Plot / PCA图",
           x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
           y = paste0("PC2 (", round(var_exp[2], 1), "%)"))
    
    if (nrow(pca_df) <= 30) {
      p <- p + geom_text(vjust = -0.5, size = 3)
    }
    
    print(p)
    dev.off()
    plots$pca <- pca_file
  }, error = function(e) message("PCA图失败 / PCA plot failed"))
  
  # 5. 表达分布图 / Expression distribution plot
  tryCatch({
    dist_file <- file.path(output_dir, "expression_distribution.png")
    png(dist_file, width = 1000, height = 600, res = 100)
    
    log_counts <- log2(counts + 1)
    log_long <- reshape2::melt(log_counts)
    colnames(log_long) <- c("Gene", "Sample", "Expression")
    
    p <- ggplot(log_long, aes(x = Expression, color = Sample)) +
      geom_density() +
      labs(title = "Expression Distribution (log2) / 表达分布（log2）",
           x = "log2(counts + 1)",
           y = "Density / 密度") +
      theme_minimal() +
      theme(legend.position = "none")
    
    print(p)
    dev.off()
    plots$distribution <- dist_file
  }, error = function(e) message("分布图失败 / Distribution plot failed"))
  
  message(paste0("QC图已保存 / QC plots saved to: ", output_dir))
  
  return(plots)
}


#' 生成RNA-seq QC摘要
#' Generate RNA-seq QC Summary
#'
#' @param qc_results QC结果 / QC results
#' @return 列表，QC摘要 / List with QC summary
generate_rnaseq_qc_summary <- function(qc_results) {
  
  summary <- list(
    n_samples = nrow(qc_results$library_stats),
    total_reads = sum(qc_results$library_stats$total_counts),
    median_reads = median(qc_results$library_stats$total_counts),
    min_reads = min(qc_results$library_stats$total_counts),
    max_reads = max(qc_results$library_stats$total_counts),
    median_detection_rate = median(qc_results$library_stats$detection_rate),
    n_outliers = length(qc_results$outliers$all_outliers),
    outlier_samples = qc_results$outliers$all_outliers,
    pca_var_pc1 = qc_results$pca_results$var_explained[1],
    pca_var_pc2 = qc_results$pca_results$var_explained[2]
  )
  
  return(summary)
}


#' 过滤低质量样本
#' Filter Low Quality Samples
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param outliers 异常样本向量 / Vector of outlier samples
#' @return 矩阵，过滤后的counts / Filtered counts matrix
#' @export
filter_low_quality_samples <- function(counts, outliers) {
  
  if (length(outliers) == 0) {
    message("没有需要过滤的样本 / No samples to filter")
    return(counts)
  }
  
  keep <- !colnames(counts) %in% outliers
  counts_filtered <- counts[, keep, drop = FALSE]
  
  message(paste0("过滤了 ", length(outliers), " 个样本 / Filtered ", 
                 length(outliers), " samples"))
  
  return(counts_filtered)
}
