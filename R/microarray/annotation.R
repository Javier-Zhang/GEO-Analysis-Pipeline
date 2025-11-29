#' ============================================================================
#' GEO 微阵列数据注释模块
#' Microarray Data Annotation Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 调用多种注释方法 / Call multiple annotation methods
#' - 探针到基因映射 / Probe to gene mapping
#' - 多探针合并策略 / Multiple probe aggregation strategies
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

suppressPackageStartupMessages({
  library(Biobase)
})

#' 微阵列数据注释
#' Annotate Microarray Data
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param method 注释方法 / Annotation method
#'   可选: "bioconductor", "biomart", "gpl"
#' @param platform GPL平台ID / GPL platform ID
#' @param organism 物种 / Organism
#' @param id_type 基因ID类型 / Gene ID type
#' @param collapse_method 多探针合并方法 / Multiple probe collapse method
#'
#' @return 列表，包含注释后的ExpressionSet和注释表 / List with annotated ExpressionSet and annotation table
#' @export
#'
#' @examples
#' \dontrun{
#' # 使用Bioconductor注释 / Annotate using Bioconductor
#' result <- annotate_microarray(eset, method = "bioconductor", organism = "human")
#' }
annotate_microarray <- function(eset,
                                 method = c("bioconductor", "biomart", "gpl"),
                                 platform = NULL,
                                 organism = "human",
                                 id_type = "SYMBOL",
                                 collapse_method = c("maxIQR", "mean", "median")) {
  
  method <- match.arg(method)
  collapse_method <- match.arg(collapse_method)
  
  message(paste0("开始微阵列注释... / Starting microarray annotation..."))
  message(paste0("方法 / Method: ", method))
  
  # 获取探针ID / Get probe IDs
  probe_ids <- featureNames(eset)
  
  # 获取平台信息 / Get platform information
  if (is.null(platform)) {
    platform <- tryCatch({
      annotation(eset)
    }, error = function(e) NULL)
  }
  
  # 根据方法选择注释函数 / Choose annotation function based on method
  annotation_table <- switch(method,
    "bioconductor" = {
      source_path <- system.file("R/annotation/bioconductor.R", package = "GEOAnalysis")
      if (source_path == "") {
        source_path <- file.path(dirname(sys.frame(1)$ofile), "../annotation/bioconductor.R")
      }
      if (file.exists(source_path)) source(source_path)
      annotate_with_bioconductor(probe_ids, organism, id_type, platform)
    },
    "biomart" = {
      source_path <- file.path(dirname(sys.frame(1)$ofile), "../annotation/biomart.R")
      if (file.exists(source_path)) source(source_path)
      annotate_with_biomart(probe_ids, organism, id_type, platform)
    },
    "gpl" = {
      if (is.null(platform)) {
        stop("GPL注释需要提供平台ID / GPL annotation requires platform ID")
      }
      source_path <- file.path(dirname(sys.frame(1)$ofile), "../annotation/gpl_parser.R")
      if (file.exists(source_path)) source(source_path)
      annotate_with_gpl(probe_ids, platform)
    }
  )
  
  # 如果注释函数不可用，使用默认方法 / If annotation function unavailable, use default
  if (is.null(annotation_table) || nrow(annotation_table) == 0) {
    message("使用默认注释方法... / Using default annotation method...")
    annotation_table <- get_default_annotation(probe_ids, platform)
  }
  
  # 合并多探针 / Collapse multiple probes
  if (collapse_method != "none") {
    message(paste0("使用 ", collapse_method, " 方法合并多探针... / Collapsing probes using ", collapse_method, "..."))
    result <- collapse_probes(eset, annotation_table, collapse_method)
  } else {
    result <- list(
      eset = eset,
      annotation = annotation_table
    )
  }
  
  message("注释完成！/ Annotation complete!")
  
  return(result)
}


#' 获取默认注释
#' Get Default Annotation
#'
#' @param probe_ids 探针ID向量 / Vector of probe IDs
#' @param platform GPL平台ID / GPL platform ID
#'
#' @return 数据框，注释表 / Data frame with annotation table
get_default_annotation <- function(probe_ids, platform = NULL) {
  
  annotation_table <- data.frame(
    PROBEID = probe_ids,
    SYMBOL = NA_character_,
    ENTREZID = NA_character_,
    GENENAME = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # 尝试从GEO获取注释 / Try to get annotation from GEO
  if (!is.null(platform) && requireNamespace("GEOquery", quietly = TRUE)) {
    tryCatch({
      gpl <- GEOquery::getGEO(platform)
      gpl_table <- GEOquery::Table(gpl)
      
      # 查找基因符号列 / Find gene symbol column
      symbol_cols <- grep("symbol|gene.symbol|gene_symbol", 
                         colnames(gpl_table), ignore.case = TRUE, value = TRUE)
      id_col <- grep("^ID$|probe.id|probeid", 
                    colnames(gpl_table), ignore.case = TRUE, value = TRUE)[1]
      
      if (length(symbol_cols) > 0 && !is.na(id_col)) {
        symbol_col <- symbol_cols[1]
        gpl_subset <- gpl_table[, c(id_col, symbol_col)]
        colnames(gpl_subset) <- c("PROBEID", "SYMBOL")
        
        annotation_table <- merge(annotation_table[, "PROBEID", drop = FALSE], 
                                  gpl_subset, 
                                  by = "PROBEID", 
                                  all.x = TRUE)
      }
    }, error = function(e) {
      message(paste0("无法从GEO获取注释 / Cannot get annotation from GEO: ", conditionMessage(e)))
    })
  }
  
  return(annotation_table)
}


#' 合并多探针到基因
#' Collapse Multiple Probes to Gene
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param annotation_table 注释表 / Annotation table
#' @param method 合并方法 / Collapse method
#'
#' @return 列表，包含合并后的ExpressionSet和注释 / List with collapsed ExpressionSet and annotation
#' @export
collapse_probes <- function(eset, annotation_table, method = c("maxIQR", "mean", "median")) {
  
  method <- match.arg(method)
  
  expr_matrix <- exprs(eset)
  
  # 确保注释表包含必要列 / Ensure annotation table has required columns
  if (!"PROBEID" %in% colnames(annotation_table)) {
    annotation_table$PROBEID <- rownames(annotation_table)
  }
  
  if (!"SYMBOL" %in% colnames(annotation_table)) {
    message("注释表缺少SYMBOL列，无法合并探针 / Annotation table missing SYMBOL column, cannot collapse probes")
    return(list(eset = eset, annotation = annotation_table))
  }
  
  # 移除无效注释 / Remove invalid annotations
  valid_probes <- !is.na(annotation_table$SYMBOL) & 
                  annotation_table$SYMBOL != "" &
                  annotation_table$PROBEID %in% rownames(expr_matrix)
  
  annotation_valid <- annotation_table[valid_probes, ]
  
  if (nrow(annotation_valid) == 0) {
    warning("没有有效的探针-基因映射 / No valid probe-gene mappings")
    return(list(eset = eset, annotation = annotation_table))
  }
  
  # 获取匹配的表达矩阵 / Get matching expression matrix
  expr_matched <- expr_matrix[annotation_valid$PROBEID, , drop = FALSE]
  
  # 根据方法合并 / Collapse based on method
  unique_genes <- unique(annotation_valid$SYMBOL)
  collapsed_matrix <- matrix(NA, 
                             nrow = length(unique_genes), 
                             ncol = ncol(expr_matrix))
  rownames(collapsed_matrix) <- unique_genes
  colnames(collapsed_matrix) <- colnames(expr_matrix)
  
  for (gene in unique_genes) {
    probe_idx <- which(annotation_valid$SYMBOL == gene)
    probe_ids <- annotation_valid$PROBEID[probe_idx]
    
    if (length(probe_ids) == 1) {
      collapsed_matrix[gene, ] <- expr_matched[probe_ids, ]
    } else {
      probe_expr <- expr_matched[probe_ids, , drop = FALSE]
      
      collapsed_matrix[gene, ] <- switch(method,
        "maxIQR" = {
          # 选择IQR最大的探针 / Select probe with max IQR
          iqrs <- apply(probe_expr, 1, IQR, na.rm = TRUE)
          best_probe <- names(iqrs)[which.max(iqrs)]
          probe_expr[best_probe, ]
        },
        "mean" = {
          colMeans(probe_expr, na.rm = TRUE)
        },
        "median" = {
          apply(probe_expr, 2, median, na.rm = TRUE)
        }
      )
    }
  }
  
  # 创建新的ExpressionSet / Create new ExpressionSet
  eset_collapsed <- ExpressionSet(
    assayData = collapsed_matrix,
    phenoData = phenoData(eset)
  )
  
  # 创建基因级注释 / Create gene-level annotation
  gene_annotation <- annotation_valid[!duplicated(annotation_valid$SYMBOL), ]
  rownames(gene_annotation) <- gene_annotation$SYMBOL
  gene_annotation <- gene_annotation[rownames(collapsed_matrix), ]
  
  message(paste0("探针合并完成：", nrow(expr_matrix), " 探针 -> ", 
                 nrow(collapsed_matrix), " 基因 / Probe collapse complete: ",
                 nrow(expr_matrix), " probes -> ", nrow(collapsed_matrix), " genes"))
  
  return(list(
    eset = eset_collapsed,
    annotation = gene_annotation,
    probe_gene_map = annotation_valid
  ))
}


#' 添加注释信息到ExpressionSet
#' Add Annotation Information to ExpressionSet
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param annotation_table 注释表 / Annotation table
#'
#' @return ExpressionSet，带注释 / ExpressionSet with annotation
#' @export
add_feature_annotation <- function(eset, annotation_table) {
  
  # 匹配探针 / Match probes
  probe_ids <- featureNames(eset)
  
  if ("PROBEID" %in% colnames(annotation_table)) {
    rownames(annotation_table) <- annotation_table$PROBEID
  }
  
  # 重排注释表 / Reorder annotation table
  annotation_matched <- annotation_table[probe_ids, , drop = FALSE]
  
  # 创建AnnotatedDataFrame / Create AnnotatedDataFrame
  fdata <- new("AnnotatedDataFrame", data = annotation_matched)
  
  featureData(eset) <- fdata
  
  return(eset)
}


#' 获取探针到基因映射表
#' Get Probe to Gene Mapping Table
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param annotation_table 注释表 / Annotation table
#'
#' @return 数据框，映射关系 / Data frame with mappings
#' @export
get_probe_gene_map <- function(eset, annotation_table = NULL) {
  
  if (!is.null(annotation_table)) {
    return(annotation_table[, c("PROBEID", "SYMBOL", "ENTREZID"), drop = FALSE])
  }
  
  # 尝试从featureData获取 / Try to get from featureData
  if (nrow(fData(eset)) > 0) {
    fdata <- fData(eset)
    
    # 查找相关列 / Find relevant columns
    symbol_col <- grep("symbol|gene", colnames(fdata), ignore.case = TRUE, value = TRUE)[1]
    entrez_col <- grep("entrez|gene.id", colnames(fdata), ignore.case = TRUE, value = TRUE)[1]
    
    probe_gene_map <- data.frame(
      PROBEID = featureNames(eset),
      SYMBOL = if (!is.na(symbol_col)) fdata[[symbol_col]] else NA,
      ENTREZID = if (!is.na(entrez_col)) fdata[[entrez_col]] else NA,
      stringsAsFactors = FALSE
    )
    
    return(probe_gene_map)
  }
  
  return(data.frame(
    PROBEID = featureNames(eset),
    SYMBOL = NA,
    ENTREZID = NA,
    stringsAsFactors = FALSE
  ))
}


#' 过滤低表达探针
#' Filter Low Expression Probes
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param threshold 表达阈值 / Expression threshold
#' @param min_samples 最小样本数 / Minimum number of samples
#'
#' @return ExpressionSet，过滤后 / Filtered ExpressionSet
#' @export
filter_low_expression <- function(eset, threshold = 1, min_samples = 2) {
  
  expr_matrix <- exprs(eset)
  
  # 计算每个探针在多少样本中高于阈值 / Count samples above threshold for each probe
  above_threshold <- rowSums(expr_matrix > threshold, na.rm = TRUE)
  
  keep_probes <- above_threshold >= min_samples
  
  eset_filtered <- eset[keep_probes, ]
  
  message(paste0("过滤低表达探针：", sum(!keep_probes), " 个被移除，",
                 sum(keep_probes), " 个保留 / Filtered low expression probes: ",
                 sum(!keep_probes), " removed, ", sum(keep_probes), " kept"))
  
  return(eset_filtered)
}
