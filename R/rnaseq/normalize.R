#' ============================================================================
#' GEO RNA-seq数据标准化模块
#' RNA-seq Data Normalization Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - DESeq2标准化（VST/rlog）/ DESeq2 normalization (VST/rlog)
#' - edgeR标准化（TMM/RLE）/ edgeR normalization (TMM/RLE)
#' - limma-voom转换 / limma-voom transformation
#' - CPM/TPM/FPKM计算 / CPM/TPM/FPKM calculation
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

#' RNA-seq数据标准化
#' Normalize RNA-seq Data
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param method 标准化方法 / Normalization method
#' @param phenotype 表型数据（可选）/ Phenotype data (optional)
#' @param design 设计公式 / Design formula
#'
#' @return 列表，包含标准化数据和相关对象 / List with normalized data and objects
#' @export
#'
#' @examples
#' \dontrun{
#' # DESeq2 VST标准化 / DESeq2 VST normalization
#' result <- normalize_rnaseq(counts, method = "vst")
#' 
#' # edgeR TMM标准化 / edgeR TMM normalization
#' result <- normalize_rnaseq(counts, method = "tmm")
#' }
normalize_rnaseq <- function(counts,
                              method = c("vst", "rlog", "tmm", "rle", "voom", "cpm"),
                              phenotype = NULL,
                              design = NULL) {
  
  method <- match.arg(method)
  
  message(paste0("执行 ", toupper(method), " 标准化... / Performing ", 
                 toupper(method), " normalization..."))
  
  result <- switch(method,
    "vst" = normalize_deseq2_vst(counts, phenotype, design),
    "rlog" = normalize_deseq2_rlog(counts, phenotype, design),
    "tmm" = normalize_edger_tmm(counts),
    "rle" = normalize_edger_rle(counts),
    "voom" = normalize_limma_voom(counts, phenotype, design),
    "cpm" = list(normalized = calculate_cpm_normalized(counts), method = "cpm")
  )
  
  result$method <- method
  
  message(paste0(toupper(method), " 标准化完成！/ ", toupper(method), " normalization complete!"))
  
  return(result)
}


#' DESeq2 VST标准化
#' DESeq2 VST Normalization
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param phenotype 表型数据 / Phenotype data
#' @param design 设计公式 / Design formula
#' @return 列表，标准化结果 / List with normalization results
normalize_deseq2_vst <- function(counts, phenotype, design) {
  
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("需要DESeq2包 / DESeq2 package required")
  }
  
  # 创建DESeqDataSet / Create DESeqDataSet
  if (is.null(phenotype)) {
    phenotype <- data.frame(
      sample = colnames(counts),
      row.names = colnames(counts)
    )
    design <- ~ 1
  }
  
  if (is.null(design)) {
    design <- ~ 1
  }
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = phenotype,
    design = design
  )
  
  # VST转换 / VST transformation
  vst_data <- DESeq2::vst(dds, blind = TRUE)
  
  return(list(
    normalized = SummarizedExperiment::assay(vst_data),
    dds = dds,
    vst_object = vst_data
  ))
}


#' DESeq2 rlog标准化
#' DESeq2 rlog Normalization
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param phenotype 表型数据 / Phenotype data
#' @param design 设计公式 / Design formula
#' @return 列表，标准化结果 / List with normalization results
normalize_deseq2_rlog <- function(counts, phenotype, design) {
  
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("需要DESeq2包 / DESeq2 package required")
  }
  
  # 创建DESeqDataSet / Create DESeqDataSet
  if (is.null(phenotype)) {
    phenotype <- data.frame(
      sample = colnames(counts),
      row.names = colnames(counts)
    )
    design <- ~ 1
  }
  
  if (is.null(design)) {
    design <- ~ 1
  }
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = phenotype,
    design = design
  )
  
  # rlog转换 / rlog transformation
  rlog_data <- DESeq2::rlog(dds, blind = TRUE)
  
  return(list(
    normalized = SummarizedExperiment::assay(rlog_data),
    dds = dds,
    rlog_object = rlog_data
  ))
}


#' edgeR TMM标准化
#' edgeR TMM Normalization
#'
#' @param counts Counts矩阵 / Counts matrix
#' @return 列表，标准化结果 / List with normalization results
normalize_edger_tmm <- function(counts) {
  
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("需要edgeR包 / edgeR package required")
  }
  
  # 创建DGEList / Create DGEList
  dge <- edgeR::DGEList(counts = counts)
  
  # TMM标准化 / TMM normalization
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  
  # 获取标准化的counts / Get normalized counts
  normalized <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
  
  return(list(
    normalized = normalized,
    dge = dge,
    norm_factors = dge$samples$norm.factors
  ))
}


#' edgeR RLE标准化
#' edgeR RLE Normalization
#'
#' @param counts Counts矩阵 / Counts matrix
#' @return 列表，标准化结果 / List with normalization results
normalize_edger_rle <- function(counts) {
  
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("需要edgeR包 / edgeR package required")
  }
  
  # 创建DGEList / Create DGEList
  dge <- edgeR::DGEList(counts = counts)
  
  # RLE标准化 / RLE normalization
  dge <- edgeR::calcNormFactors(dge, method = "RLE")
  
  # 获取标准化的counts / Get normalized counts
  normalized <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
  
  return(list(
    normalized = normalized,
    dge = dge,
    norm_factors = dge$samples$norm.factors
  ))
}


#' limma-voom转换
#' limma-voom Transformation
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param phenotype 表型数据 / Phenotype data
#' @param design 设计公式或矩阵 / Design formula or matrix
#' @return 列表，标准化结果 / List with normalization results
normalize_limma_voom <- function(counts, phenotype, design) {
  
  if (!requireNamespace("limma", quietly = TRUE) || 
      !requireNamespace("edgeR", quietly = TRUE)) {
    stop("需要limma和edgeR包 / limma and edgeR packages required")
  }
  
  # 创建DGEList并TMM标准化 / Create DGEList and TMM normalize
  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  
  # 创建设计矩阵 / Create design matrix
  if (is.null(design)) {
    design_matrix <- matrix(1, ncol(counts), 1)
  } else if (is.matrix(design)) {
    design_matrix <- design
  } else {
    design_matrix <- model.matrix(design, data = phenotype)
  }
  
  # voom转换 / voom transformation
  v <- limma::voom(dge, design_matrix, plot = FALSE)
  
  return(list(
    normalized = v$E,
    weights = v$weights,
    dge = dge,
    voom_object = v
  ))
}


#' 计算标准化CPM
#' Calculate Normalized CPM
#'
#' @param counts Counts矩阵 / Counts matrix
#' @return 矩阵，标准化CPM / Normalized CPM matrix
calculate_cpm_normalized <- function(counts) {
  
  lib_sizes <- colSums(counts)
  cpm <- t(t(counts) / lib_sizes * 1e6)
  log_cpm <- log2(cpm + 1)
  
  return(log_cpm)
}


#' 比较不同标准化方法
#' Compare Different Normalization Methods
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param methods 要比较的方法 / Methods to compare
#' @param output_dir 输出目录 / Output directory
#'
#' @return 列表，各方法结果 / List with results from each method
#' @export
compare_normalization_methods <- function(counts,
                                           methods = c("vst", "tmm", "voom", "cpm"),
                                           output_dir = NULL) {
  
  results <- list()
  
  for (method in methods) {
    message(paste0("\n处理方法 / Processing method: ", method))
    tryCatch({
      results[[method]] <- normalize_rnaseq(counts, method = method)
    }, error = function(e) {
      warning(paste0(method, " 失败 / failed: ", conditionMessage(e)))
      results[[method]] <<- NULL
    })
  }
  
  # 生成比较图 / Generate comparison plots
  if (!is.null(output_dir)) {
    generate_normalization_comparison_plots(results, output_dir)
  }
  
  return(results)
}


#' 生成标准化比较图
#' Generate Normalization Comparison Plots
#'
#' @param results 标准化结果列表 / List of normalization results
#' @param output_dir 输出目录 / Output directory
generate_normalization_comparison_plots <- function(results, output_dir) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 箱线图比较 / Boxplot comparison
  tryCatch({
    n_methods <- sum(!sapply(results, is.null))
    png(file.path(output_dir, "normalization_comparison.png"),
        width = 400 * n_methods, height = 400, res = 100)
    
    par(mfrow = c(1, n_methods))
    
    for (method in names(results)) {
      if (!is.null(results[[method]])) {
        boxplot(results[[method]]$normalized[, 1:min(20, ncol(results[[method]]$normalized))],
                main = paste0(toupper(method), " Normalization"),
                las = 2, col = "steelblue")
      }
    }
    
    dev.off()
    message(paste0("比较图已保存 / Comparison plot saved to: ", output_dir))
  }, error = function(e) {
    message("比较图生成失败 / Comparison plot failed")
  })
}


#' 获取可用的标准化方法
#' Get Available Normalization Methods
#'
#' @return 字符向量，可用方法 / Character vector of available methods
#' @export
get_rnaseq_normalization_methods <- function() {
  
  methods <- c("cpm")
  
  if (requireNamespace("DESeq2", quietly = TRUE)) {
    methods <- c(methods, "vst", "rlog")
  }
  
  if (requireNamespace("edgeR", quietly = TRUE)) {
    methods <- c(methods, "tmm", "rle")
  }
  
  if (requireNamespace("limma", quietly = TRUE) && 
      requireNamespace("edgeR", quietly = TRUE)) {
    methods <- c(methods, "voom")
  }
  
  return(methods)
}
