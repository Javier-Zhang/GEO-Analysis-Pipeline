#' ============================================================================
#' GEO RNA-seq注释模块
#' RNA-seq Annotation Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 基因ID转换 / Gene ID conversion
#' - 注释信息添加 / Add annotation information
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

#' RNA-seq数据注释
#' Annotate RNA-seq Data
#'
#' @param counts Counts矩阵或基因ID向量 / Counts matrix or gene ID vector
#' @param from_type 输入ID类型 / Input ID type
#' @param to_type 输出ID类型 / Output ID type
#' @param organism 物种 / Organism
#' @param method 注释方法 / Annotation method
#'
#' @return 数据框，注释结果 / Data frame with annotation results
#' @export
#'
#' @examples
#' \dontrun{
#' # 从Ensembl转换到Symbol / Convert from Ensembl to Symbol
#' annotation <- annotate_rnaseq(gene_ids, from_type = "ENSEMBL", to_type = "SYMBOL")
#' }
annotate_rnaseq <- function(counts,
                             from_type = c("ENSEMBL", "ENTREZID", "SYMBOL", "REFSEQ"),
                             to_type = c("SYMBOL", "ENTREZID", "ENSEMBL", "GENENAME"),
                             organism = c("human", "mouse"),
                             method = c("bioconductor", "biomart")) {
  
  from_type <- match.arg(from_type)
  to_type <- match.arg(to_type)
  organism <- match.arg(organism)
  method <- match.arg(method)
  
  message(paste0("注释RNA-seq数据... / Annotating RNA-seq data..."))
  message(paste0("转换: ", from_type, " -> ", to_type))
  
  # 获取基因ID / Get gene IDs
  if (is.matrix(counts) || is.data.frame(counts)) {
    gene_ids <- rownames(counts)
  } else {
    gene_ids <- counts
  }
  
  # 清理Ensembl ID版本号 / Clean Ensembl ID version numbers
  if (from_type == "ENSEMBL") {
    gene_ids_clean <- gsub("\\.\\d+$", "", gene_ids)
  } else {
    gene_ids_clean <- gene_ids
  }
  
  # 根据方法注释 / Annotate based on method
  annotation <- switch(method,
    "bioconductor" = annotate_with_org_db(gene_ids_clean, from_type, to_type, organism),
    "biomart" = annotate_with_biomart_rnaseq(gene_ids_clean, from_type, to_type, organism)
  )
  
  # 保留原始ID / Keep original IDs
  annotation$original_id <- gene_ids
  
  message(paste0("注释完成！成功率 / Annotation complete! Success rate: ",
                 round(sum(!is.na(annotation[[to_type]])) / length(gene_ids) * 100, 1), "%"))
  
  return(annotation)
}


#' 使用org.db注释
#' Annotate Using org.db
#'
#' @param gene_ids 基因ID向量 / Gene ID vector
#' @param from_type 输入类型 / Input type
#' @param to_type 输出类型 / Output type
#' @param organism 物种 / Organism
#' @return 数据框，注释结果 / Data frame with annotation results
annotate_with_org_db <- function(gene_ids, from_type, to_type, organism) {
  
  # 选择注释包 / Select annotation package
  org_db <- switch(organism,
    "human" = {
      if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        stop("需要org.Hs.eg.db包 / org.Hs.eg.db package required")
      }
      org.Hs.eg.db::org.Hs.eg.db
    },
    "mouse" = {
      if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
        stop("需要org.Mm.eg.db包 / org.Mm.eg.db package required")
      }
      org.Mm.eg.db::org.Mm.eg.db
    }
  )
  
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("需要AnnotationDbi包 / AnnotationDbi package required")
  }
  
  # ID映射 / ID mapping
  tryCatch({
    mapping <- AnnotationDbi::select(
      org_db,
      keys = gene_ids,
      columns = c(from_type, to_type, "GENENAME"),
      keytype = from_type
    )
    
    # 处理多对一映射 / Handle many-to-one mappings
    mapping <- mapping[!duplicated(mapping[[from_type]]), ]
    
    return(mapping)
    
  }, error = function(e) {
    warning(paste0("注释失败 / Annotation failed: ", conditionMessage(e)))
    return(data.frame(
      gene_id = gene_ids,
      to_type = NA_character_,
      stringsAsFactors = FALSE
    ))
  })
}


#' 使用biomaRt注释RNA-seq
#' Annotate RNA-seq Using biomaRt
#'
#' @param gene_ids 基因ID向量 / Gene ID vector
#' @param from_type 输入类型 / Input type
#' @param to_type 输出类型 / Output type
#' @param organism 物种 / Organism
#' @return 数据框，注释结果 / Data frame with annotation results
annotate_with_biomart_rnaseq <- function(gene_ids, from_type, to_type, organism) {
  
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    warning("biomaRt包未安装，使用备选方法 / biomaRt not installed, using alternative")
    return(annotate_with_org_db(gene_ids, from_type, to_type, organism))
  }
  
  # 选择数据集 / Select dataset
  dataset <- switch(organism,
    "human" = "hsapiens_gene_ensembl",
    "mouse" = "mmusculus_gene_ensembl"
  )
  
  # 属性映射 / Attribute mapping
  attr_map <- list(
    "ENSEMBL" = "ensembl_gene_id",
    "ENTREZID" = "entrezgene_id",
    "SYMBOL" = "external_gene_name",
    "GENENAME" = "description",
    "REFSEQ" = "refseq_mrna"
  )
  
  tryCatch({
    mart <- biomaRt::useMart("ensembl", dataset = dataset)
    
    annotation <- biomaRt::getBM(
      attributes = c(attr_map[[from_type]], attr_map[[to_type]], "description"),
      filters = attr_map[[from_type]],
      values = gene_ids,
      mart = mart
    )
    
    colnames(annotation) <- c(from_type, to_type, "GENENAME")
    
    return(annotation)
    
  }, error = function(e) {
    warning(paste0("biomaRt查询失败 / biomaRt query failed: ", conditionMessage(e)))
    return(annotate_with_org_db(gene_ids, from_type, to_type, organism))
  })
}


#' 添加注释到counts矩阵
#' Add Annotation to Counts Matrix
#'
#' @param counts Counts矩阵 / Counts matrix
#' @param annotation 注释数据框 / Annotation data frame
#' @param collapse 是否合并同基因行 / Whether to collapse rows with same gene
#' @param collapse_method 合并方法 / Collapse method
#'
#' @return 列表，包含注释后的counts / List with annotated counts
#' @export
add_annotation_to_counts <- function(counts,
                                      annotation,
                                      collapse = TRUE,
                                      collapse_method = c("mean", "sum", "max")) {
  
  collapse_method <- match.arg(collapse_method)
  
  # 匹配注释 / Match annotation
  if (!"original_id" %in% colnames(annotation)) {
    annotation$original_id <- annotation[[1]]
  }
  
  idx <- match(rownames(counts), annotation$original_id)
  matched_annotation <- annotation[idx, ]
  rownames(matched_annotation) <- rownames(counts)
  
  # 如果需要合并 / If collapse needed
  if (collapse && "SYMBOL" %in% colnames(matched_annotation)) {
    # 只保留有效注释的基因 / Keep only genes with valid annotation
    valid <- !is.na(matched_annotation$SYMBOL) & matched_annotation$SYMBOL != ""
    counts_valid <- counts[valid, ]
    symbols <- matched_annotation$SYMBOL[valid]
    
    # 合并 / Collapse
    unique_symbols <- unique(symbols)
    collapsed <- matrix(0, nrow = length(unique_symbols), ncol = ncol(counts_valid))
    rownames(collapsed) <- unique_symbols
    colnames(collapsed) <- colnames(counts_valid)
    
    for (sym in unique_symbols) {
      rows <- which(symbols == sym)
      if (length(rows) == 1) {
        collapsed[sym, ] <- counts_valid[rows, ]
      } else {
        collapsed[sym, ] <- switch(collapse_method,
          "mean" = colMeans(counts_valid[rows, , drop = FALSE]),
          "sum" = colSums(counts_valid[rows, , drop = FALSE]),
          "max" = apply(counts_valid[rows, , drop = FALSE], 2, max)
        )
      }
    }
    
    message(paste0("合并后基因数 / Genes after collapse: ", nrow(collapsed)))
    
    return(list(
      counts = collapsed,
      annotation = matched_annotation[valid, ][!duplicated(symbols), ]
    ))
  }
  
  return(list(
    counts = counts,
    annotation = matched_annotation
  ))
}


#' 检测基因ID类型
#' Detect Gene ID Type
#'
#' @param gene_ids 基因ID向量 / Gene ID vector
#' @return 字符串，ID类型 / Character, ID type
#' @export
detect_gene_id_type <- function(gene_ids) {
  
  sample_ids <- head(gene_ids, 100)
  
  # Ensembl格式 / Ensembl format
  if (mean(grepl("^ENS[A-Z]*G\\d+", sample_ids)) > 0.5) {
    return("ENSEMBL")
  }
  
  # Entrez（纯数字）/ Entrez (pure numbers)
  if (mean(grepl("^\\d+$", sample_ids)) > 0.5) {
    return("ENTREZID")
  }
  
  # RefSeq格式 / RefSeq format
  if (mean(grepl("^N[MR]_\\d+", sample_ids)) > 0.5) {
    return("REFSEQ")
  }
  
  # 默认为Symbol / Default to Symbol
  return("SYMBOL")
}


#' 获取基因长度
#' Get Gene Lengths
#'
#' @param gene_ids 基因ID向量 / Gene ID vector
#' @param organism 物种 / Organism
#'
#' @return 数值向量，基因长度 / Numeric vector of gene lengths
#' @export
get_gene_lengths <- function(gene_ids, organism = "human") {
  
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    warning("biomaRt未安装，无法获取基因长度 / biomaRt not installed")
    return(NULL)
  }
  
  dataset <- switch(organism,
    "human" = "hsapiens_gene_ensembl",
    "mouse" = "mmusculus_gene_ensembl"
  )
  
  tryCatch({
    mart <- biomaRt::useMart("ensembl", dataset = dataset)
    
    id_type <- detect_gene_id_type(gene_ids)
    filter <- switch(id_type,
      "ENSEMBL" = "ensembl_gene_id",
      "SYMBOL" = "external_gene_name",
      "ENTREZID" = "entrezgene_id"
    )
    
    # 清理ID / Clean IDs
    clean_ids <- gsub("\\.\\d+$", "", gene_ids)
    
    results <- biomaRt::getBM(
      attributes = c(filter, "transcript_length"),
      filters = filter,
      values = clean_ids,
      mart = mart
    )
    
    # 取每个基因的最长转录本 / Take longest transcript for each gene
    lengths <- tapply(results$transcript_length, results[[filter]], max)
    
    gene_lengths <- lengths[clean_ids]
    names(gene_lengths) <- gene_ids
    
    return(gene_lengths)
    
  }, error = function(e) {
    warning(paste0("获取基因长度失败 / Failed to get gene lengths: ", conditionMessage(e)))
    return(NULL)
  })
}
