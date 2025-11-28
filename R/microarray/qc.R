#' ============================================================================
#' GEO 微阵列数据质控模块
#' Microarray Data Quality Control Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 样本质量评估 / Sample quality assessment
#' - 异常值检测 / Outlier detection
#' - 批次效应检测 / Batch effect detection
#' - 质控可视化 / QC visualization
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

suppressPackageStartupMessages({
  library(Biobase)
  library(affy)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

#' 执行完整的微阵列质控分析
#' Perform Complete Microarray QC Analysis
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param output_dir 输出目录 / Output directory
#' @param batch_var 批次变量名（可选）/ Batch variable name (optional)
#'
#' @return 列表，包含所有质控结果 / List containing all QC results
#' @export
run_microarray_qc <- function(eset, output_dir = "qc_results", batch_var = NULL) {
  
  # 创建输出目录 / Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  message("开始微阵列质控分析... / Starting microarray QC analysis...")
  
  qc_results <- list(
    sample_stats = NULL,
    correlation_matrix = NULL,
    pca_results = NULL,
    outliers = NULL,
    batch_effect = NULL,
    plots = list()
  )
  
  # 1. 基本统计 / Basic statistics
  message("[1/6] 计算样本统计量... / Calculating sample statistics...")
  qc_results$sample_stats <- calculate_sample_stats(eset)
  
  # 2. 相关性分析 / Correlation analysis
  message("[2/6] 计算样本间相关性... / Calculating sample correlations...")
  qc_results$correlation_matrix <- calculate_sample_correlation(eset)
  
  # 3. PCA分析 / PCA analysis
  message("[3/6] 执行PCA分析... / Performing PCA analysis...")
  qc_results$pca_results <- perform_pca_analysis(eset)
  
  # 4. 异常值检测 / Outlier detection
  message("[4/6] 检测异常样本... / Detecting outlier samples...")
  qc_results$outliers <- detect_outliers(eset, qc_results$correlation_matrix)
  
  # 5. 批次效应检测 / Batch effect detection
  if (!is.null(batch_var)) {
    message("[5/6] 检测批次效应... / Detecting batch effects...")
    qc_results$batch_effect <- detect_batch_effect(eset, batch_var)
  }
  
  # 6. 生成质控图 / Generate QC plots
  message("[6/6] 生成质控图... / Generating QC plots...")
  qc_results$plots <- generate_qc_plots(eset, qc_results, output_dir)
  
  # 生成质控摘要 / Generate QC summary
  qc_results$summary <- generate_qc_summary(qc_results)
  
  message("质控分析完成！/ QC analysis complete!")
  
  return(qc_results)
}

#' 计算样本统计量
#' Calculate Sample Statistics
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @return 数据框，样本统计信息 / Data frame with sample statistics
calculate_sample_stats <- function(eset) {
  
  expr_matrix <- exprs(eset)
  
  stats <- data.frame(
    sample_id = colnames(expr_matrix),
    mean_expression = colMeans(expr_matrix, na.rm = TRUE),
    median_expression = apply(expr_matrix, 2, median, na.rm = TRUE),
    sd_expression = apply(expr_matrix, 2, sd, na.rm = TRUE),
    min_expression = apply(expr_matrix, 2, min, na.rm = TRUE),
    max_expression = apply(expr_matrix, 2, max, na.rm = TRUE),
    iqr_expression = apply(expr_matrix, 2, IQR, na.rm = TRUE),
    n_missing = colSums(is.na(expr_matrix)),
    n_zeros = colSums(expr_matrix == 0, na.rm = TRUE),
    detection_rate = colSums(expr_matrix > 0, na.rm = TRUE) / nrow(expr_matrix) * 100,
    stringsAsFactors = FALSE
  )
  
  return(stats)
}

#' 计算样本间相关性
#' Calculate Sample Correlation
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param method 相关性方法 / Correlation method
#' @return 相关性矩阵 / Correlation matrix
calculate_sample_correlation <- function(eset, method = "pearson") {
  
  expr_matrix <- exprs(eset)
  cor_matrix <- cor(expr_matrix, method = method, use = "pairwise.complete.obs")
  
  return(cor_matrix)
}

#' 执行PCA分析
#' Perform PCA Analysis
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param n_components PCA组件数量 / Number of PCA components
#' @return 列表，PCA结果 / List with PCA results
perform_pca_analysis <- function(eset, n_components = 10) {
  
  expr_matrix <- exprs(eset)
  
  # 移除零方差基因 / Remove zero variance genes
  vars <- apply(expr_matrix, 1, var, na.rm = TRUE)
  expr_filtered <- expr_matrix[vars > 0, ]
  
  # 执行PCA / Perform PCA
  pca_result <- prcomp(t(expr_filtered), center = TRUE, scale. = TRUE)
  
  # 计算方差解释比例 / Calculate variance explained
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
  
  result <- list(
    pca = pca_result,
    scores = pca_result$x[, 1:min(n_components, ncol(pca_result$x))],
    loadings = pca_result$rotation[, 1:min(n_components, ncol(pca_result$rotation))],
    var_explained = var_explained[1:min(n_components, length(var_explained))],
    cumulative_var = cumsum(var_explained)[1:min(n_components, length(var_explained))]
  )
  
  return(result)
}

#' 检测异常样本
#' Detect Outlier Samples
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param cor_matrix 相关性矩阵 / Correlation matrix
#' @param threshold 异常阈值 / Outlier threshold
#' @return 列表，异常样本信息 / List with outlier information
detect_outliers <- function(eset, cor_matrix, threshold = 0.8) {
  
  # 方法1：基于平均相关性 / Method 1: Based on mean correlation
  mean_cors <- colMeans(cor_matrix, na.rm = TRUE)
  cor_outliers <- names(mean_cors)[mean_cors < threshold]
  
  # 方法2：基于表达量分布 / Method 2: Based on expression distribution
  expr_matrix <- exprs(eset)
  sample_means <- colMeans(expr_matrix, na.rm = TRUE)
  sample_sds <- apply(expr_matrix, 2, sd, na.rm = TRUE)
  
  # 使用MAD检测异常值 / Use MAD for outlier detection
  mean_mad <- mad(sample_means, na.rm = TRUE)
  mean_median <- median(sample_means, na.rm = TRUE)
  expr_outliers <- colnames(expr_matrix)[abs(sample_means - mean_median) > 3 * mean_mad]
  
  # 方法3：基于PCA / Method 3: Based on PCA
  pca <- prcomp(t(expr_matrix[apply(expr_matrix, 1, var, na.rm = TRUE) > 0, ]), 
                center = TRUE, scale. = TRUE)
  scores <- pca$x[, 1:2]
  
  # 计算Mahalanobis距离 / Calculate Mahalanobis distance
  center <- colMeans(scores)
  cov_matrix <- cov(scores)
  mahal_dist <- mahalanobis(scores, center, cov_matrix)
  pca_outliers <- rownames(scores)[mahal_dist > qchisq(0.975, df = 2)]
  
  result <- list(
    correlation_outliers = cor_outliers,
    expression_outliers = expr_outliers,
    pca_outliers = pca_outliers,
    all_outliers = unique(c(cor_outliers, expr_outliers, pca_outliers)),
    mean_correlations = mean_cors,
    mahalanobis_distances = mahal_dist,
    outlier_summary = data.frame(
      sample_id = colnames(expr_matrix),
      mean_correlation = mean_cors,
      is_cor_outlier = colnames(expr_matrix) %in% cor_outliers,
      is_expr_outlier = colnames(expr_matrix) %in% expr_outliers,
      is_pca_outlier = colnames(expr_matrix) %in% pca_outliers,
      stringsAsFactors = FALSE
    )
  )
  
  return(result)
}

#' 检测批次效应
#' Detect Batch Effect
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param batch_var 批次变量名 / Batch variable name
#' @return 列表，批次效应检测结果 / List with batch effect detection results
detect_batch_effect <- function(eset, batch_var) {
  
  if (!batch_var %in% varLabels(eset)) {
    warning(paste0("批次变量 '", batch_var, "' 不存在于表型数据中 / Batch variable not found in phenotype data"))
    return(NULL)
  }
  
  batch <- pData(eset)[[batch_var]]
  expr_matrix <- exprs(eset)
  
  # PVCA分析（简化版）/ PVCA analysis (simplified)
  # 计算批次解释的方差 / Calculate variance explained by batch
  pca <- prcomp(t(expr_matrix[apply(expr_matrix, 1, var, na.rm = TRUE) > 0, ]), 
                center = TRUE, scale. = TRUE)
  
  # 对前10个PC进行方差分解 / Variance decomposition for first 10 PCs
  n_pcs <- min(10, ncol(pca$x))
  batch_var_explained <- numeric(n_pcs)
  
  for (i in 1:n_pcs) {
    if (length(unique(batch)) > 1) {
      model <- aov(pca$x[, i] ~ batch)
      ss <- summary(model)[[1]]
      batch_var_explained[i] <- ss$`Sum Sq`[1] / sum(ss$`Sum Sq`) * 100
    }
  }
  
  result <- list(
    batch_variable = batch_var,
    n_batches = length(unique(batch)),
    batch_sizes = table(batch),
    batch_variance_per_pc = batch_var_explained,
    total_batch_effect = mean(batch_var_explained, na.rm = TRUE),
    has_significant_batch_effect = mean(batch_var_explained, na.rm = TRUE) > 10
  )
  
  return(result)
}

#' 生成质控图
#' Generate QC Plots
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param qc_results 质控结果 / QC results
#' @param output_dir 输出目录 / Output directory
#' @return 列表，图对象路径 / List of plot file paths
generate_qc_plots <- function(eset, qc_results, output_dir) {
  
  plots <- list()
  expr_matrix <- exprs(eset)
  
  # 1. 密度图 / Density plot
  tryCatch({
    density_file <- file.path(output_dir, "density_plot.png")
    png(density_file, width = 1200, height = 800, res = 100)
    
    expr_long <- reshape2::melt(expr_matrix)
    colnames(expr_long) <- c("Gene", "Sample", "Expression")
    
    p <- ggplot(expr_long, aes(x = Expression, color = Sample)) +
      geom_density(alpha = 0.5) +
      theme_minimal() +
      labs(title = "Expression Density Distribution / 表达密度分布",
           x = "Expression Value / 表达值",
           y = "Density / 密度") +
      theme(legend.position = "none")
    print(p)
    dev.off()
    plots$density <- density_file
  }, error = function(e) message("密度图生成失败 / Density plot failed"))
  
  # 2. 箱线图 / Box plot
  tryCatch({
    boxplot_file <- file.path(output_dir, "boxplot.png")
    png(boxplot_file, width = max(1200, ncol(expr_matrix) * 20), height = 800, res = 100)
    
    par(mar = c(10, 4, 4, 2))
    boxplot(expr_matrix, las = 2, col = rainbow(ncol(expr_matrix)),
            main = "Sample Expression Distribution / 样本表达分布",
            ylab = "Expression Value / 表达值")
    dev.off()
    plots$boxplot <- boxplot_file
  }, error = function(e) message("箱线图生成失败 / Boxplot failed"))
  
  # 3. 相关性热图 / Correlation heatmap
  tryCatch({
    heatmap_file <- file.path(output_dir, "correlation_heatmap.png")
    png(heatmap_file, width = 1000, height = 1000, res = 100)
    
    pheatmap(qc_results$correlation_matrix,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             main = "Sample Correlation Heatmap / 样本相关性热图",
             show_rownames = ncol(expr_matrix) <= 50,
             show_colnames = ncol(expr_matrix) <= 50)
    dev.off()
    plots$correlation_heatmap <- heatmap_file
  }, error = function(e) message("相关性热图生成失败 / Correlation heatmap failed"))
  
  # 4. PCA图 / PCA plot
  tryCatch({
    pca_file <- file.path(output_dir, "pca_plot.png")
    png(pca_file, width = 1000, height = 800, res = 100)
    
    pca_df <- data.frame(
      PC1 = qc_results$pca_results$scores[, 1],
      PC2 = qc_results$pca_results$scores[, 2],
      Sample = rownames(qc_results$pca_results$scores)
    )
    
    var_exp <- qc_results$pca_results$var_explained;
    
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
      geom_point(size = 3, color = "steelblue") +
      geom_text(vjust = -0.5, size = 3) +
      theme_minimal() +
      labs(title = "PCA Plot / PCA图",
           x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
           y = paste0("PC2 (", round(var_exp[2], 1), "%)")
    )
    print(p)
    dev.off()
    plots$pca <- pca_file
  }, error = function(e) message("PCA图生成失败 / PCA plot failed"))
  
  # 5. 方差解释图 / Variance explained plot
  tryCatch({
    scree_file <- file.path(output_dir, "scree_plot.png")
    png(scree_file, width = 800, height = 600, res = 100)
    
    var_df <- data.frame(
      PC = paste0("PC", 1:length(qc_results$pca_results$var_explained)),
      Variance = qc_results$pca_results$var_explained,
      Cumulative = qc_results$pca_results$cumulative_var
    )
    var_df$PC <- factor(var_df$PC, levels = var_df$PC)
    
    p <- ggplot(var_df, aes(x = PC, y = Variance)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_line(aes(y = Cumulative, group = 1), color = "red", size = 1) +
      geom_point(aes(y = Cumulative), color = "red", size = 2) +
      theme_minimal() +
      labs(title = "Variance Explained by PCs / PC方差解释",
           x = "Principal Component / 主成分",
           y = "Variance Explained (%) / 方差解释 (%)")
    print(p)
    dev.off()
    plots$scree <- scree_file
  }, error = function(e) message("方差解释图生成失败 / Scree plot failed"))
  
  # 6. 平均相关性图 / Mean correlation plot
  tryCatch({
    mean_cor_file <- file.path(output_dir, "mean_correlation.png")
    png(mean_cor_file, width = 1000, height = 600, res = 100)
    
    mean_cor_df <- data.frame(
      Sample = names(qc_results$outliers$mean_correlations),
      MeanCorrelation = qc_results$outliers$mean_correlations
    )
    mean_cor_df$IsOutlier <- mean_cor_df$Sample %in% qc_results$outliers$correlation_outliers
    mean_cor_df <- mean_cor_df[order(mean_cor_df$MeanCorrelation), ]
    mean_cor_df$Sample <- factor(mean_cor_df$Sample, levels = mean_cor_df$Sample)
    
    p <- ggplot(mean_cor_df, aes(x = Sample, y = MeanCorrelation, fill = IsOutlier)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "red")) +
      geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title = "Mean Sample Correlation / 平均样本相关性",
           x = "Sample / 样本",
           y = "Mean Correlation / 平均相关性")
    print(p)
    dev.off()
    plots$mean_correlation <- mean_cor_file
  }, error = function(e) message("平均相关性图生成失败 / Mean correlation plot failed"))
  
  return(plots)
}

#' 生成质控摘要
#' Generate QC Summary
#'
#' @param qc_results 质控结果 / QC results
#' @return 列表，质控摘要 / List with QC summary
generate_qc_summary <- function(qc_results) {
  
  summary <- list(
    n_samples = nrow(qc_results$sample_stats),
    n_outliers = length(qc_results$outliers$all_outliers),
    outlier_samples = qc_results$outliers$all_outliers,
    mean_correlation = mean(qc_results$outliers$mean_correlations, na.rm = TRUE),
    min_correlation = min(qc_results$outliers$mean_correlations, na.rm = TRUE),
    pca_var_pc1 = qc_results$pca_results$var_explained[1],
    pca_var_pc2 = qc_results$pca_results$var_explained[2],
    qc_pass = length(qc_results$outliers$all_outliers) == 0
  )
  
  # 质控建议 / QC recommendations
  summary$recommendations <- character()
  
  if (summary$n_outliers > 0) {
    summary$recommendations <- c(summary$recommendations,
      paste0("发现 ", summary$n_outliers, " 个异常样本，建议检查或移除 / ",
             summary$n_outliers, " outlier samples detected, consider removal"))
  }
  
  if (summary$min_correlation < 0.7) {
    summary$recommendations <- c(summary$recommendations,
      "部分样本相关性过低，可能存在质量问题 / Some samples have low correlation, potential quality issues")
  }
  
  if (summary$pca_var_pc1 > 50) {
    summary$recommendations <- c(summary$recommendations,
      "PC1解释方差过高，可能存在批次效应 / PC1 explains high variance, potential batch effect")
  }
  
  if (length(summary$recommendations) == 0) {
    summary$recommendations <- "所有样本通过质控 / All samples passed QC"
  }
  
  return(summary)
}

#' 移除异常样本
#' Remove Outlier Samples
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param outliers 要移除的样本名 / Sample names to remove
#' @return ExpressionSet，移除异常样本后 / ExpressionSet with outliers removed
#' @export
remove_outlier_samples <- function(eset, outliers) {
  
  if (length(outliers) == 0) {
    message("没有异常样本需要移除 / No outlier samples to remove")
    return(eset)
  }
  
  keep_samples <- !sampleNames(eset) %in% outliers
  eset_clean <- eset[, keep_samples]
  
  message(paste0("移除了 ", length(outliers), " 个异常样本 / Removed ", 
                 length(outliers), " outlier samples"))
  message(paste0("剩余 ", ncol(eset_clean), " 个样本 / ", 
                 ncol(eset_clean), " samples remaining"))
  
  return(eset_clean)
}