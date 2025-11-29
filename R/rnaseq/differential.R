#' ============================================================================
#' GEO RNA-seq差异表达分析模块
#' RNA-seq Differential Expression Analysis Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - DESeq2差异分析 / DESeq2 differential analysis
#' - edgeR差异分析 / edgeR differential analysis
#' - limma-voom差异分析 / limma-voom differential analysis
#' - 结果整合比较 / Result integration and comparison
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

#' 执行RNA-seq差异表达分析
#' Perform RNA-seq Differential Expression Analysis
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param phenotype 表型数据框 / Phenotype data frame
#' @param group_col 分组列名 / Group column name
#' @param contrast 对比 / Contrast
#' @param method 分析方法 / Analysis method
#' @param lfc_threshold logFC阈值 / logFC threshold
#' @param p_threshold P值阈值 / P-value threshold
#' @param output_dir 输出目录 / Output directory
#'
#' @return 列表，差异分析结果 / List with differential analysis results
#' @export
#'
#' @examples
#' \dontrun{
#' # DESeq2分析 / DESeq2 analysis
#' result <- run_rnaseq_differential(counts, phenotype, "condition", 
#'                                   contrast = c("condition", "cancer", "normal"))
#' }
run_rnaseq_differential <- function(counts,
                                     phenotype,
                                     group_col,
                                     contrast = NULL,
                                     method = c("deseq2", "edger", "limma_voom"),
                                     lfc_threshold = 1,
                                     p_threshold = 0.05,
                                     output_dir = NULL) {
  
  method <- match.arg(method)
  
  message(paste0("执行 ", toupper(method), " 差异分析... / Performing ", 
                 toupper(method), " differential analysis..."))
  
  # 确保样本匹配 / Ensure sample matching
  common_samples <- intersect(colnames(counts), rownames(phenotype))
  counts <- counts[, common_samples]
  phenotype <- phenotype[common_samples, , drop = FALSE]
  
  result <- switch(method,
    "deseq2" = run_deseq2_analysis(counts, phenotype, group_col, contrast),
    "edger" = run_edger_analysis(counts, phenotype, group_col, contrast),
    "limma_voom" = run_limma_voom_analysis(counts, phenotype, group_col, contrast)
  )
  
  # 添加差异状态 / Add differential status
  result$deg_table$DEG_status <- "Not Significant"
  result$deg_table$DEG_status[result$deg_table$padj < p_threshold & 
                               result$deg_table$log2FoldChange > lfc_threshold] <- "Up"
  result$deg_table$DEG_status[result$deg_table$padj < p_threshold & 
                               result$deg_table$log2FoldChange < -lfc_threshold] <- "Down"
  
  # 统计 / Summary
  result$summary <- list(
    n_total = nrow(result$deg_table),
    n_up = sum(result$deg_table$DEG_status == "Up"),
    n_down = sum(result$deg_table$DEG_status == "Down"),
    n_significant = sum(result$deg_table$DEG_status != "Not Significant"),
    method = method,
    lfc_threshold = lfc_threshold,
    p_threshold = p_threshold
  )
  
  message(paste0("上调基因 / Up-regulated: ", result$summary$n_up))
  message(paste0("下调基因 / Down-regulated: ", result$summary$n_down))
  
  # 生成图 / Generate plots
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    result$plots <- generate_rnaseq_deg_plots(result, output_dir, lfc_threshold, p_threshold)
  }
  
  message("差异分析完成！/ Differential analysis complete!")
  
  return(result)
}


#' DESeq2差异分析
#' DESeq2 Differential Analysis
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param phenotype 表型数据 / Phenotype data
#' @param group_col 分组列 / Group column
#' @param contrast 对比 / Contrast
#' @return 列表，分析结果 / List with analysis results
run_deseq2_analysis <- function(counts, phenotype, group_col, contrast) {
  
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("需要DESeq2包 / DESeq2 package required")
  }
  
  # 创建设计公式 / Create design formula
  design <- as.formula(paste0("~", group_col))
  
  # 确保分组是因子 / Ensure group is factor
  phenotype[[group_col]] <- factor(phenotype[[group_col]])
  
  # 创建DESeqDataSet / Create DESeqDataSet
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = phenotype,
    design = design
  )
  
  # 过滤低表达基因 / Filter low expression genes
  keep <- rowSums(DESeq2::counts(dds)) >= 10
  dds <- dds[keep, ]
  
  # 运行DESeq2 / Run DESeq2
  dds <- DESeq2::DESeq(dds)
  
  # 提取结果 / Extract results
  if (!is.null(contrast)) {
    res <- DESeq2::results(dds, contrast = contrast)
  } else {
    res <- DESeq2::results(dds)
  }
  
  # 转换为数据框 / Convert to data frame
  deg_table <- as.data.frame(res)
  deg_table$gene <- rownames(deg_table)
  deg_table <- deg_table[order(deg_table$padj), ]
  
  return(list(
    deg_table = deg_table,
    dds = dds,
    method = "deseq2"
  ))
}


#' edgeR差异分析
#' edgeR Differential Analysis
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param phenotype 表型数据 / Phenotype data
#' @param group_col 分组列 / Group column
#' @param contrast 对比 / Contrast
#' @return 列表，分析结果 / List with analysis results
run_edger_analysis <- function(counts, phenotype, group_col, contrast) {
  
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("需要edgeR包 / edgeR package required")
  }
  
  group <- factor(phenotype[[group_col]])
  
  # 创建DGEList / Create DGEList
  dge <- edgeR::DGEList(counts = counts, group = group)
  
  # 过滤低表达基因 / Filter low expression genes
  keep <- edgeR::filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # TMM标准化 / TMM normalization
  dge <- edgeR::calcNormFactors(dge)
  
  # 设计矩阵 / Design matrix
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  # 估计离散度 / Estimate dispersion
  dge <- edgeR::estimateDisp(dge, design)
  
  # 拟合模型 / Fit model
  fit <- edgeR::glmQLFit(dge, design)
  
  # 创建对比 / Create contrast
  if (!is.null(contrast) && length(contrast) == 3) {
    con <- limma::makeContrasts(
      contrasts = paste0(contrast[2], "-", contrast[3]),
      levels = design
    )
  } else {
    con <- limma::makeContrasts(
      contrasts = paste0(levels(group)[2], "-", levels(group)[1]),
      levels = design
    )
  }
  
  # 差异检验 / Differential test
  qlf <- edgeR::glmQLFTest(fit, contrast = con)
  
  # 提取结果 / Extract results
  deg_table <- edgeR::topTags(qlf, n = Inf)$table
  deg_table$gene <- rownames(deg_table)
  
  # 统一列名 / Standardize column names
  colnames(deg_table)[colnames(deg_table) == "logFC"] <- "log2FoldChange"
  colnames(deg_table)[colnames(deg_table) == "FDR"] <- "padj"
  colnames(deg_table)[colnames(deg_table) == "PValue"] <- "pvalue"
  
  return(list(
    deg_table = deg_table,
    dge = dge,
    fit = fit,
    method = "edger"
  ))
}


#' limma-voom差异分析
#' limma-voom Differential Analysis
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param phenotype 表型数据 / Phenotype data
#' @param group_col 分组列 / Group column
#' @param contrast 对比 / Contrast
#' @return 列表，分析结果 / List with analysis results
run_limma_voom_analysis <- function(counts, phenotype, group_col, contrast) {
  
  if (!requireNamespace("limma", quietly = TRUE) || 
      !requireNamespace("edgeR", quietly = TRUE)) {
    stop("需要limma和edgeR包 / limma and edgeR packages required")
  }
  
  group <- factor(phenotype[[group_col]])
  
  # 创建DGEList / Create DGEList
  dge <- edgeR::DGEList(counts = counts, group = group)
  
  # 过滤 / Filter
  keep <- edgeR::filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # TMM标准化 / TMM normalization
  dge <- edgeR::calcNormFactors(dge)
  
  # 设计矩阵 / Design matrix
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  # voom转换 / voom transformation
  v <- limma::voom(dge, design, plot = FALSE)
  
  # 拟合 / Fit
  fit <- limma::lmFit(v, design)
  
  # 创建对比 / Create contrast
  if (!is.null(contrast) && length(contrast) == 3) {
    con <- limma::makeContrasts(
      contrasts = paste0(contrast[2], "-", contrast[3]),
      levels = design
    )
  } else {
    con <- limma::makeContrasts(
      contrasts = paste0(levels(group)[2], "-", levels(group)[1]),
      levels = design
    )
  }
  
  fit2 <- limma::contrasts.fit(fit, con)
  fit2 <- limma::eBayes(fit2)
  
  # 提取结果 / Extract results
  deg_table <- limma::topTable(fit2, number = Inf)
  deg_table$gene <- rownames(deg_table)
  
  # 统一列名 / Standardize column names
  colnames(deg_table)[colnames(deg_table) == "logFC"] <- "log2FoldChange"
  colnames(deg_table)[colnames(deg_table) == "adj.P.Val"] <- "padj"
  colnames(deg_table)[colnames(deg_table) == "P.Value"] <- "pvalue"
  
  return(list(
    deg_table = deg_table,
    voom = v,
    fit = fit2,
    method = "limma_voom"
  ))
}


#' 生成RNA-seq差异分析图
#' Generate RNA-seq Differential Analysis Plots
#'
#' @param result 差异分析结果 / Differential analysis result
#' @param output_dir 输出目录 / Output directory
#' @param lfc_threshold logFC阈值 / logFC threshold
#' @param p_threshold P值阈值 / P-value threshold
#' @return 列表，图路径 / List of plot paths
generate_rnaseq_deg_plots <- function(result, output_dir, lfc_threshold, p_threshold) {
  
  plots <- list()
  deg_table <- result$deg_table
  
  # 1. 火山图 / Volcano plot
  tryCatch({
    volcano_file <- file.path(output_dir, "volcano_plot.png")
    png(volcano_file, width = 1000, height = 800, res = 100)
    
    volcano_df <- data.frame(
      logFC = deg_table$log2FoldChange,
      negLogP = -log10(deg_table$padj),
      Status = deg_table$DEG_status
    )
    volcano_df$negLogP[is.infinite(volcano_df$negLogP)] <- max(volcano_df$negLogP[is.finite(volcano_df$negLogP)]) + 1
    
    p <- ggplot(volcano_df, aes(x = logFC, y = negLogP, color = Status)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("Down" = "#3366CC", 
                                    "Not Significant" = "grey60", 
                                    "Up" = "#CC3366")) +
      geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed") +
      geom_hline(yintercept = -log10(p_threshold), linetype = "dashed") +
      labs(title = "Volcano Plot / 火山图",
           x = "log2 Fold Change",
           y = "-log10(adjusted p-value)") +
      theme_minimal()
    
    print(p)
    dev.off()
    plots$volcano <- volcano_file
  }, error = function(e) message("火山图失败 / Volcano plot failed"))
  
  # 2. MA图 / MA plot
  tryCatch({
    ma_file <- file.path(output_dir, "ma_plot.png")
    png(ma_file, width = 1000, height = 800, res = 100)
    
    if ("baseMean" %in% colnames(deg_table)) {
      ma_df <- data.frame(
        A = log2(deg_table$baseMean + 1),
        M = deg_table$log2FoldChange,
        Status = deg_table$DEG_status
      )
    } else {
      ma_df <- data.frame(
        A = deg_table$AveExpr,
        M = deg_table$log2FoldChange,
        Status = deg_table$DEG_status
      )
    }
    
    p <- ggplot(ma_df, aes(x = A, y = M, color = Status)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("Down" = "#3366CC", 
                                    "Not Significant" = "grey60", 
                                    "Up" = "#CC3366")) +
      geom_hline(yintercept = 0, color = "black") +
      labs(title = "MA Plot / MA图",
           x = "Average Expression",
           y = "log2 Fold Change") +
      theme_minimal()
    
    print(p)
    dev.off()
    plots$ma <- ma_file
  }, error = function(e) message("MA图失败 / MA plot failed"))
  
  message(paste0("图已保存 / Plots saved to: ", output_dir))
  
  return(plots)
}


#' 比较多种差异分析方法
#' Compare Multiple Differential Analysis Methods
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param phenotype 表型数据 / Phenotype data
#' @param group_col 分组列 / Group column
#' @param contrast 对比 / Contrast
#' @param methods 要比较的方法 / Methods to compare
#'
#' @return 列表，各方法结果 / List with results from each method
#' @export
compare_deg_methods <- function(counts,
                                 phenotype,
                                 group_col,
                                 contrast = NULL,
                                 methods = c("deseq2", "edger", "limma_voom")) {
  
  results <- list()
  
  for (method in methods) {
    message(paste0("\n运行 / Running: ", method))
    tryCatch({
      results[[method]] <- run_rnaseq_differential(
        counts, phenotype, group_col, contrast, method = method
      )
    }, error = function(e) {
      warning(paste0(method, " 失败 / failed: ", conditionMessage(e)))
    })
  }
  
  # 创建比较摘要 / Create comparison summary
  comparison <- create_method_comparison(results)
  
  return(list(
    results = results,
    comparison = comparison
  ))
}


#' 创建方法比较摘要
#' Create Method Comparison Summary
#'
#' @param results 各方法结果列表 / List of results from each method
#' @return 数据框，比较摘要 / Data frame with comparison summary
create_method_comparison <- function(results) {
  
  summary_df <- data.frame(
    Method = character(),
    Total_Tested = numeric(),
    Up_Regulated = numeric(),
    Down_Regulated = numeric(),
    Total_Significant = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (method in names(results)) {
    if (!is.null(results[[method]])) {
      summary_df <- rbind(summary_df, data.frame(
        Method = method,
        Total_Tested = results[[method]]$summary$n_total,
        Up_Regulated = results[[method]]$summary$n_up,
        Down_Regulated = results[[method]]$summary$n_down,
        Total_Significant = results[[method]]$summary$n_significant,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(summary_df)
}
