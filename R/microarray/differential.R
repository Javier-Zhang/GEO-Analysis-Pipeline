#' ============================================================================
#' GEO 微阵列差异表达分析模块
#' Microarray Differential Expression Analysis Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - limma差异表达分析 / limma differential expression analysis
#' - 支持两组/多组/配对设计 / Support two-group/multi-group/paired design
#' - 火山图、MA图、热图可视化 / Volcano plot, MA plot, heatmap visualization
#' - 结果筛选和导出 / Result filtering and export
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

#' 执行limma差异表达分析
#' Perform limma Differential Expression Analysis
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param group 分组变量（向量或表型列名）/ Group variable (vector or phenotype column name)
#' @param contrast 对比（如"Disease-Control"）/ Contrast (e.g., "Disease-Control")
#' @param design_type 设计类型 / Design type
#'   可选: "two_group", "multi_group", "paired", "continuous"
#' @param paired_var 配对变量（配对设计时使用）/ Paired variable (for paired design)
#' @param adjust_method 多重检验校正方法 / Multiple testing correction method
#' @param lfc_threshold logFC阈值 / logFC threshold
#' @param p_threshold P值阈值 / P-value threshold
#' @param output_dir 输出目录 / Output directory
#'
#' @return 列表，包含差异分析结果和可视化 / List containing DEG results and visualizations
#' @export
#'
#' @examples
#' \dontrun{
#' # 两组比较 / Two-group comparison
#' result <- run_differential_analysis(eset, group = "condition", 
#'                                     contrast = "Cancer-Normal")
#' 
#' # 宫颈癌数据分析示例 / Cervical cancer analysis example
#' # GSE6791: 宫颈癌 vs 正常
#' result <- run_differential_analysis(eset, group = pData(eset)$disease_state,
#'                                     contrast = "cervical_cancer-normal")
#' }
run_differential_analysis <- function(eset,
                                       group,
                                       contrast = NULL,
                                       design_type = c("two_group", "multi_group", "paired", "continuous"),
                                       paired_var = NULL,
                                       adjust_method = "BH",
                                       lfc_threshold = 1,
                                       p_threshold = 0.05,
                                       output_dir = NULL) {
  
  design_type <- match.arg(design_type)
  
  message("开始差异表达分析... / Starting differential expression analysis...")
  
  # 获取分组信息 / Get group information
  if (is.character(group) && length(group) == 1 && group %in% varLabels(eset)) {
    group_vector <- pData(eset)[[group]]
  } else {
    group_vector <- group
  }
  
  # 确保分组是因子 / Ensure group is factor
  group_factor <- factor(group_vector)
  
  # 构建设计矩阵 / Build design matrix
  design_result <- build_design_matrix(eset, group_factor, design_type, paired_var)
  design <- design_result$design
  
  message(paste0("设计类型 / Design type: ", design_type))
  message(paste0("分组水平 / Group levels: ", paste(levels(group_factor), collapse = ", ")))
  
  # 拟合线性模型 / Fit linear model
  message("拟合线性模型... / Fitting linear model...")
  fit <- limma::lmFit(eset, design)
  
  # 构建对比矩阵 / Build contrast matrix
  if (!is.null(contrast)) {
    contrast_matrix <- limma::makeContrasts(contrasts = contrast, levels = design)
    fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  } else {
    fit2 <- fit
  }
  
  # 贝叶斯平滑 / Bayesian smoothing
  fit2 <- limma::eBayes(fit2)
  
  # 提取结果 / Extract results
  message("提取差异基因... / Extracting differentially expressed genes...")
  results <- limma::topTable(fit2, number = Inf, adjust.method = adjust_method)
  
  # 添加差异表达状态 / Add differential expression status
  results$DEG_status <- "Not Significant"
  results$DEG_status[results$adj.P.Val < p_threshold & results$logFC > lfc_threshold] <- "Up"
  results$DEG_status[results$adj.P.Val < p_threshold & results$logFC < -lfc_threshold] <- "Down"
  results$DEG_status <- factor(results$DEG_status, levels = c("Down", "Not Significant", "Up"))
  
  # 统计结果 / Summary statistics
  n_up <- sum(results$DEG_status == "Up")
  n_down <- sum(results$DEG_status == "Down")
  
  message(paste0("上调基因 / Up-regulated: ", n_up))
  message(paste0("下调基因 / Down-regulated: ", n_down))
  
  # 生成可视化 / Generate visualizations
  plots <- list()
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    plots <- generate_deg_plots(eset, results, group_factor, output_dir, 
                                lfc_threshold, p_threshold)
  }
  
  # 返回结果 / Return results
  result <- list(
    deg_table = results,
    fit = fit2,
    design = design,
    contrast = contrast,
    summary = list(
      n_total = nrow(results),
      n_up = n_up,
      n_down = n_down,
      n_significant = n_up + n_down,
      lfc_threshold = lfc_threshold,
      p_threshold = p_threshold,
      adjust_method = adjust_method
    ),
    plots = plots
  )
  
  message("差异分析完成！/ Differential analysis complete!")
  
  return(result)
}


#' 构建设计矩阵
#' Build Design Matrix
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param group 分组因子 / Group factor
#' @param design_type 设计类型 / Design type
#' @param paired_var 配对变量 / Paired variable
#'
#' @return 列表，包含设计矩阵 / List containing design matrix
build_design_matrix <- function(eset, group, design_type, paired_var = NULL) {
  
  design <- switch(design_type,
    "two_group" = {
      model.matrix(~ 0 + group)
    },
    "multi_group" = {
      model.matrix(~ 0 + group)
    },
    "paired" = {
      if (is.null(paired_var)) {
        stop("配对设计需要提供paired_var / Paired design requires paired_var")
      }
      if (is.character(paired_var) && paired_var %in% varLabels(eset)) {
        block <- factor(pData(eset)[[paired_var]])
      } else {
        block <- factor(paired_var)
      }
      model.matrix(~ 0 + group + block)
    },
    "continuous" = {
      model.matrix(~ group)
    }
  )
  
  colnames(design) <- gsub("group", "", colnames(design))
  
  return(list(design = design))
}


#' 生成差异分析可视化图
#' Generate Differential Expression Plots
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param results 差异分析结果 / DEG results
#' @param group 分组因子 / Group factor
#' @param output_dir 输出目录 / Output directory
#' @param lfc_threshold logFC阈值 / logFC threshold
#' @param p_threshold P值阈值 / P-value threshold
#'
#' @return 列表，图文件路径 / List of plot file paths
generate_deg_plots <- function(eset, results, group, output_dir, lfc_threshold, p_threshold) {
  
  plots <- list()
  
  # 1. 火山图 / Volcano plot
  tryCatch({
    volcano_file <- file.path(output_dir, "volcano_plot.png")
    png(volcano_file, width = 1000, height = 800, res = 100)
    
    volcano_df <- data.frame(
      logFC = results$logFC,
      negLogP = -log10(results$adj.P.Val),
      Status = results$DEG_status,
      Gene = rownames(results)
    )
    
    # 标记top基因 / Label top genes
    top_genes <- head(volcano_df[order(volcano_df$negLogP, decreasing = TRUE), ], 10)
    
    p <- ggplot(volcano_df, aes(x = logFC, y = negLogP, color = Status)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("Down" = "#3366CC", 
                                    "Not Significant" = "grey60", 
                                    "Up" = "#CC3366"),
                        name = "Expression / 表达状态") +
      geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), 
                 linetype = "dashed", color = "grey40") +
      geom_hline(yintercept = -log10(p_threshold), 
                 linetype = "dashed", color = "grey40") +
      labs(title = "Volcano Plot / 火山图",
           x = "log2 Fold Change",
           y = "-log10(adj. P-value)") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # 添加基因标签 / Add gene labels
    if (nrow(top_genes) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = top_genes,
        aes(label = Gene),
        size = 3,
        max.overlaps = 20
      )
    }
    
    print(p)
    dev.off()
    plots$volcano <- volcano_file
  }, error = function(e) {
    message(paste0("火山图生成失败 / Volcano plot failed: ", conditionMessage(e)))
  })
  
  # 2. MA图 / MA plot
  tryCatch({
    ma_file <- file.path(output_dir, "ma_plot.png")
    png(ma_file, width = 1000, height = 800, res = 100)
    
    ma_df <- data.frame(
      A = results$AveExpr,
      M = results$logFC,
      Status = results$DEG_status
    )
    
    p <- ggplot(ma_df, aes(x = A, y = M, color = Status)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("Down" = "#3366CC", 
                                    "Not Significant" = "grey60", 
                                    "Up" = "#CC3366"),
                        name = "Expression / 表达状态") +
      geom_hline(yintercept = c(-lfc_threshold, 0, lfc_threshold), 
                 linetype = c("dashed", "solid", "dashed"), 
                 color = c("grey40", "black", "grey40")) +
      labs(title = "MA Plot / MA图",
           x = "Average Expression",
           y = "log2 Fold Change") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    print(p)
    dev.off()
    plots$ma <- ma_file
  }, error = function(e) {
    message(paste0("MA图生成失败 / MA plot failed: ", conditionMessage(e)))
  })
  
  # 3. 热图（top差异基因）/ Heatmap (top DEGs)
  tryCatch({
    heatmap_file <- file.path(output_dir, "deg_heatmap.png")
    
    # 获取top差异基因 / Get top DEGs
    sig_genes <- rownames(results)[results$DEG_status != "Not Significant"]
    
    if (length(sig_genes) > 0) {
      # 限制基因数量 / Limit number of genes
      n_genes <- min(50, length(sig_genes))
      top_sig <- head(sig_genes[order(abs(results[sig_genes, "logFC"]), decreasing = TRUE)], n_genes)
      
      expr_matrix <- exprs(eset)[top_sig, , drop = FALSE]
      
      # 标准化用于热图 / Scale for heatmap
      expr_scaled <- t(scale(t(expr_matrix)))
      
      # 注释列 / Annotation columns
      annotation_col <- data.frame(
        Group = group,
        row.names = colnames(expr_matrix)
      )
      
      png(heatmap_file, width = 1200, height = 1000, res = 100)
      
      pheatmap(expr_scaled,
               annotation_col = annotation_col,
               color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = nrow(expr_scaled) <= 50,
               show_colnames = ncol(expr_scaled) <= 30,
               main = paste0("Top ", n_genes, " DEGs Heatmap / 差异基因热图"),
               fontsize_row = 8,
               fontsize_col = 8)
      
      dev.off()
      plots$heatmap <- heatmap_file
    }
  }, error = function(e) {
    message(paste0("热图生成失败 / Heatmap failed: ", conditionMessage(e)))
  })
  
  # 4. P值分布图 / P-value distribution
  tryCatch({
    pval_file <- file.path(output_dir, "pvalue_distribution.png")
    png(pval_file, width = 800, height = 600, res = 100)
    
    p <- ggplot(data.frame(pvalue = results$P.Value), aes(x = pvalue)) +
      geom_histogram(bins = 50, fill = "steelblue", color = "white") +
      labs(title = "P-value Distribution / P值分布",
           x = "P-value",
           y = "Frequency / 频数") +
      theme_minimal()
    
    print(p)
    dev.off()
    plots$pvalue_dist <- pval_file
  }, error = function(e) {
    message(paste0("P值分布图生成失败 / P-value distribution plot failed"))
  })
  
  message(paste0("图表已保存到 / Plots saved to: ", output_dir))
  
  return(plots)
}


#' 筛选差异表达基因
#' Filter Differentially Expressed Genes
#'
#' @param deg_result 差异分析结果 / DEG analysis result
#' @param lfc_threshold logFC阈值 / logFC threshold
#' @param p_threshold 调整后P值阈值 / Adjusted P-value threshold
#' @param direction 方向 / Direction
#'   可选: "both", "up", "down"
#'
#' @return 数据框，筛选后的差异基因 / Filtered DEG data frame
#' @export
filter_degs <- function(deg_result, 
                        lfc_threshold = 1, 
                        p_threshold = 0.05,
                        direction = c("both", "up", "down")) {
  
  direction <- match.arg(direction)
  
  results <- deg_result$deg_table
  
  # 基础筛选 / Basic filtering
  filtered <- results[results$adj.P.Val < p_threshold, ]
  
  # 方向筛选 / Direction filtering
  filtered <- switch(direction,
    "both" = filtered[abs(filtered$logFC) > lfc_threshold, ],
    "up" = filtered[filtered$logFC > lfc_threshold, ],
    "down" = filtered[filtered$logFC < -lfc_threshold, ]
  )
  
  # 按logFC排序 / Sort by logFC
  filtered <- filtered[order(abs(filtered$logFC), decreasing = TRUE), ]
  
  message(paste0("筛选出 ", nrow(filtered), " 个差异基因 / Filtered ", 
                 nrow(filtered), " DEGs"))
  
  return(filtered)
}


#' 导出差异分析结果
#' Export Differential Analysis Results
#'
#' @param deg_result 差异分析结果 / DEG analysis result
#' @param output_file 输出文件路径 / Output file path
#' @param format 输出格式 / Output format
#'
#' @return 无返回值 / No return value
#' @export
export_deg_results <- function(deg_result, output_file, format = c("csv", "xlsx", "tsv")) {
  
  format <- match.arg(format)
  
  results <- deg_result$deg_table
  
  switch(format,
    "csv" = {
      write.csv(results, output_file, row.names = TRUE)
    },
    "xlsx" = {
      if (requireNamespace("openxlsx", quietly = TRUE)) {
        openxlsx::write.xlsx(results, output_file, rowNames = TRUE)
      } else {
        warning("openxlsx包未安装，使用CSV格式 / openxlsx not installed, using CSV format")
        write.csv(results, gsub("\\.xlsx$", ".csv", output_file), row.names = TRUE)
      }
    },
    "tsv" = {
      write.table(results, output_file, sep = "\t", row.names = TRUE, quote = FALSE)
    }
  )
  
  message(paste0("结果已导出到 / Results exported to: ", output_file))
}


#' 多组比较分析
#' Multi-group Comparison Analysis
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param group 分组变量 / Group variable
#' @param contrasts 对比列表 / List of contrasts
#' @param ... 其他参数传递给run_differential_analysis
#'
#' @return 列表，每个对比的结果 / List of results for each contrast
#' @export
multi_group_analysis <- function(eset, group, contrasts, ...) {
  
  results <- list()
  
  for (i in seq_along(contrasts)) {
    contrast <- contrasts[i]
    message(paste0("\n===== 分析对比 / Analyzing contrast: ", contrast, " ====="))
    
    results[[contrast]] <- tryCatch({
      run_differential_analysis(eset, group = group, contrast = contrast, ...)
    }, error = function(e) {
      warning(paste0("对比 ", contrast, " 分析失败 / Contrast analysis failed: ", 
                     conditionMessage(e)))
      NULL
    })
  }
  
  return(results)
}


#' 获取差异分析摘要
#' Get Differential Analysis Summary
#'
#' @param deg_result 差异分析结果 / DEG analysis result
#'
#' @return 字符串，摘要文本 / Character, summary text
#' @export
get_deg_summary <- function(deg_result) {
  
  summary <- deg_result$summary
  
  text <- paste0(
    "差异表达分析摘要 / Differential Expression Analysis Summary\n",
    "================================================\n",
    "总基因数 / Total genes: ", summary$n_total, "\n",
    "上调基因 / Up-regulated: ", summary$n_up, "\n",
    "下调基因 / Down-regulated: ", summary$n_down, "\n",
    "显著差异基因 / Significant DEGs: ", summary$n_significant, "\n",
    "logFC阈值 / logFC threshold: ", summary$lfc_threshold, "\n",
    "P值阈值 / P-value threshold: ", summary$p_threshold, "\n",
    "校正方法 / Adjustment method: ", summary$adjust_method, "\n"
  )
  
  return(text)
}
