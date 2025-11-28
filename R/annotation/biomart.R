#' ============================================================================
#' biomaRt在线注释模块
#' biomaRt Online Annotation Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 连接Ensembl数据库 / Connect to Ensembl database
#' - 批量查询注释 / Batch annotation query
#' - 支持多物种 / Support multiple species
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

#' 使用biomaRt进行注释
#' Annotate Using biomaRt
#'
#' @param ids ID向量 / ID vector
#' @param organism 物种 / Organism
#' @param id_type ID类型 / ID type
#' @param attributes 要获取的属性 / Attributes to retrieve
#'
#' @return 数据框，注释结果 / Data frame with annotation results
#' @export
#'
#' @examples
#' \dontrun{
#' # 使用Ensembl ID注释 / Annotate using Ensembl IDs
#' annotation <- annotate_with_biomart(ensembl_ids, organism = "human", id_type = "ensembl")
#' }
annotate_with_biomart <- function(ids,
                                   organism = c("human", "mouse", "rat"),
                                   id_type = c("ensembl", "symbol", "entrez", "refseq"),
                                   attributes = NULL) {
  
  organism <- match.arg(organism)
  id_type <- match.arg(id_type)
  
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    warning("biomaRt包未安装 / biomaRt package not installed")
    return(create_empty_biomart_annotation(ids))
  }
  
  message(paste0("使用biomaRt注释 (", organism, ")... / Annotating with biomaRt..."))
  
  # 获取数据集 / Get dataset
  dataset <- get_biomart_dataset(organism)
  
  # 获取mart连接 / Get mart connection
  mart <- connect_to_biomart(dataset)
  
  if (is.null(mart)) {
    return(create_empty_biomart_annotation(ids))
  }
  
  # 获取filter和attributes / Get filter and attributes
  filter <- get_biomart_filter(id_type)
  
  if (is.null(attributes)) {
    attributes <- c(
      "ensembl_gene_id",
      "external_gene_name",
      "entrezgene_id",
      "description",
      "chromosome_name",
      "start_position",
      "end_position",
      "gene_biotype"
    )
  }
  
  # 清理ID / Clean IDs
  clean_ids <- clean_gene_ids(ids, id_type)
  
  # 批量查询 / Batch query
  annotation <- batch_biomart_query(mart, clean_ids, filter, attributes)
  
  # 添加原始ID / Add original IDs
  if (!is.null(annotation) && nrow(annotation) > 0) {
    annotation$original_id <- ids[match(annotation[[filter]], clean_ids)]
  }
  
  return(annotation)
}


#' 获取biomaRt数据集名称
#' Get biomaRt Dataset Name
#'
#' @param organism 物种 / Organism
#' @return 字符串，数据集名 / Character, dataset name
get_biomart_dataset <- function(organism) {
  
  switch(organism,
    "human" = "hsapiens_gene_ensembl",
    "mouse" = "mmusculus_gene_ensembl",
    "rat" = "rnorvegicus_gene_ensembl"
  )
}


#' 连接biomaRt
#' Connect to biomaRt
#'
#' @param dataset 数据集名 / Dataset name
#' @return mart对象 / mart object
connect_to_biomart <- function(dataset) {
  
  # 尝试多个镜像 / Try multiple mirrors
  mirrors <- c("www", "useast", "uswest", "asia")
  
  for (mirror in mirrors) {
    tryCatch({
      message(paste0("尝试连接 / Trying mirror: ", mirror))
      
      mart <- biomaRt::useMart(
        biomart = "ensembl",
        dataset = dataset,
        host = paste0(mirror, ".ensembl.org")
      )
      
      message("连接成功 / Connection successful")
      return(mart)
      
    }, error = function(e) {
      message(paste0("镜像 ", mirror, " 失败 / Mirror failed"))
    })
  }
  
  warning("无法连接biomaRt / Cannot connect to biomaRt")
  return(NULL)
}


#' 获取biomaRt filter名称
#' Get biomaRt Filter Name
#'
#' @param id_type ID类型 / ID type
#' @return 字符串，filter名 / Character, filter name
get_biomart_filter <- function(id_type) {
  
  switch(id_type,
    "ensembl" = "ensembl_gene_id",
    "symbol" = "external_gene_name",
    "entrez" = "entrezgene_id",
    "refseq" = "refseq_mrna"
  )
}


#' 清理基因ID
#' Clean Gene IDs
#'
#' @param ids ID向量 / ID vector
#' @param id_type ID类型 / ID type
#' @return 字符向量，清理后的ID / Character vector of cleaned IDs
clean_gene_ids <- function(ids, id_type) {
  
  clean_ids <- ids
  
  # 移除Ensembl版本号 / Remove Ensembl version numbers
  if (id_type == "ensembl") {
    clean_ids <- gsub("\\.\\d+$", "", clean_ids)
  }
  
  # 移除空白 / Remove whitespace
  clean_ids <- trimws(clean_ids)
  
  return(clean_ids)
}


#' 批量biomaRt查询
#' Batch biomaRt Query
#'
#' @param mart mart对象 / mart object
#' @param ids ID向量 / ID vector
#' @param filter filter名 / Filter name
#' @param attributes 属性列表 / Attribute list
#' @param batch_size 批次大小 / Batch size
#' @return 数据框，查询结果 / Data frame with query results
batch_biomart_query <- function(mart, ids, filter, attributes, batch_size = 500) {
  
  n_ids <- length(ids)
  n_batches <- ceiling(n_ids / batch_size)
  
  results <- list()
  
  for (i in seq_len(n_batches)) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, n_ids)
    batch_ids <- ids[start_idx:end_idx]
    
    message(paste0("查询批次 / Querying batch ", i, "/", n_batches))
    
    tryCatch({
      batch_result <- biomaRt::getBM(
        attributes = attributes,
        filters = filter,
        values = batch_ids,
        mart = mart
      )
      
      results[[i]] <- batch_result
      
    }, error = function(e) {
      warning(paste0("批次 ", i, " 查询失败 / Batch ", i, " query failed"))
    })
    
    # 避免请求过快 / Avoid too fast requests
    Sys.sleep(0.5)
  }
  
  # 合并结果 / Combine results
  if (length(results) > 0) {
    combined <- do.call(rbind, results)
    combined <- combined[!duplicated(combined), ]
    return(combined)
  }
  
  return(NULL)
}


#' 创建空biomaRt注释
#' Create Empty biomaRt Annotation
#'
#' @param ids ID向量 / ID vector
#' @return 数据框，空注释 / Data frame with empty annotation
create_empty_biomart_annotation <- function(ids) {
  
  data.frame(
    ensembl_gene_id = NA_character_,
    external_gene_name = NA_character_,
    entrezgene_id = NA_integer_,
    description = NA_character_,
    original_id = ids,
    stringsAsFactors = FALSE
  )
}


#' 获取基因坐标信息
#' Get Gene Coordinate Information
#'
#' @param gene_ids 基因ID / Gene IDs
#' @param organism 物种 / Organism
#' @param id_type ID类型 / ID type
#'
#' @return 数据框，坐标信息 / Data frame with coordinate information
#' @export
get_gene_coordinates <- function(gene_ids,
                                  organism = "human",
                                  id_type = "symbol") {
  
  attributes <- c(
    get_biomart_filter(id_type),
    "chromosome_name",
    "start_position",
    "end_position",
    "strand"
  )
  
  annotation <- annotate_with_biomart(
    ids = gene_ids,
    organism = organism,
    id_type = id_type,
    attributes = attributes
  )
  
  return(annotation)
}


#' 获取GO注释
#' Get GO Annotation
#'
#' @param gene_ids 基因ID / Gene IDs
#' @param organism 物种 / Organism
#' @param id_type ID类型 / ID type
#' @param go_type GO类型 / GO type
#'
#' @return 数据框，GO注释 / Data frame with GO annotation
#' @export
get_go_annotation <- function(gene_ids,
                               organism = "human",
                               id_type = "symbol",
                               go_type = c("biological_process", "molecular_function", "cellular_component")) {
  
  go_type <- match.arg(go_type)
  
  attributes <- c(
    get_biomart_filter(id_type),
    "go_id",
    "name_1006",
    "namespace_1003"
  )
  
  annotation <- annotate_with_biomart(
    ids = gene_ids,
    organism = organism,
    id_type = id_type,
    attributes = attributes
  )
  
  # 筛选GO类型 / Filter GO type
  if (!is.null(annotation) && nrow(annotation) > 0 && "namespace_1003" %in% colnames(annotation)) {
    annotation <- annotation[annotation$namespace_1003 == go_type, ]
  }
  
  return(annotation)
}


#' 获取同源基因
#' Get Homologous Genes
#'
#' @param gene_ids 基因ID / Gene IDs
#' @param from_organism 源物种 / Source organism
#' @param to_organism 目标物种 / Target organism
#' @param id_type ID类型 / ID type
#'
#' @return 数据框，同源基因映射 / Data frame with homolog mapping
#' @export
get_homologs <- function(gene_ids,
                          from_organism = "human",
                          to_organism = "mouse",
                          id_type = "symbol") {
  
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    return(NULL)
  }
  
  # 获取源数据集 / Get source dataset
  from_dataset <- get_biomart_dataset(from_organism)
  from_mart <- connect_to_biomart(from_dataset)
  
  if (is.null(from_mart)) {
    return(NULL)
  }
  
  # 获取同源基因属性名 / Get homolog attribute name
  to_prefix <- switch(to_organism,
    "human" = "hsapiens",
    "mouse" = "mmusculus",
    "rat" = "rnorvegicus"
  )
  
  filter <- get_biomart_filter(id_type)
  attributes <- c(
    filter,
    paste0(to_prefix, "_homolog_ensembl_gene"),
    paste0(to_prefix, "_homolog_associated_gene_name"),
    paste0(to_prefix, "_homolog_orthology_type")
  )
  
  tryCatch({
    homologs <- biomaRt::getBM(
      attributes = attributes,
      filters = filter,
      values = gene_ids,
      mart = from_mart
    )
    
    return(homologs)
    
  }, error = function(e) {
    warning(paste0("获取同源基因失败 / Failed to get homologs: ", conditionMessage(e)))
    return(NULL)
  })
}


#' 列出可用的biomaRt数据集
#' List Available biomaRt Datasets
#'
#' @param pattern 筛选模式 / Filter pattern
#' @return 数据框，可用数据集 / Data frame with available datasets
#' @export
list_biomart_datasets <- function(pattern = NULL) {
  
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    return(NULL)
  }
  
  tryCatch({
    ensembl <- biomaRt::useMart("ensembl")
    datasets <- biomaRt::listDatasets(ensembl)
    
    if (!is.null(pattern)) {
      datasets <- datasets[grep(pattern, datasets$dataset, ignore.case = TRUE), ]
    }
    
    return(datasets)
    
  }, error = function(e) {
    warning("无法获取数据集列表 / Cannot get dataset list")
    return(NULL)
  })
}
