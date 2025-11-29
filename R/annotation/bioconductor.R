#' ============================================================================
#' Bioconductor注释模块
#' Bioconductor Annotation Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 使用org.Hs.eg.db, org.Mm.eg.db等进行注释
#' - Use org.Hs.eg.db, org.Mm.eg.db etc. for annotation
#' - 自动检测物种 / Auto-detect organism
#' - ID转换（Symbol/Entrez/Ensembl）/ ID conversion
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

#' 使用Bioconductor org.db包注释
#' Annotate Using Bioconductor org.db Packages
#'
#' @param probe_ids 探针/基因ID向量 / Vector of probe/gene IDs
#' @param organism 物种 / Organism
#' @param to_type 目标ID类型 / Target ID type
#' @param platform GPL平台（可选）/ GPL platform (optional)
#'
#' @return 数据框，注释结果 / Data frame with annotation results
#' @export
#'
#' @examples
#' \dontrun{
#' # 人类基因注释 / Human gene annotation
#' annotation <- annotate_with_bioconductor(gene_ids, organism = "human")
#' }
annotate_with_bioconductor <- function(probe_ids,
                                        organism = c("human", "mouse", "rat"),
                                        to_type = "SYMBOL",
                                        platform = NULL) {
  
  organism <- match.arg(organism)
  
  message(paste0("使用Bioconductor注释 (", organism, ")... / ",
                 "Annotating with Bioconductor (", organism, ")..."))
  
  # 获取org.db包 / Get org.db package
  org_db <- get_org_db(organism)
  
  if (is.null(org_db)) {
    warning("无法加载org.db包 / Cannot load org.db package")
    return(create_empty_annotation(probe_ids))
  }
  
  # 检测输入ID类型 / Detect input ID type
  from_type <- detect_id_type(probe_ids)
  message(paste0("检测到输入ID类型 / Detected input ID type: ", from_type))
  
  # 如果是Affymetrix探针，需要特殊处理 / If Affymetrix probes, need special handling
  if (!is.null(platform) && grepl("^GPL", platform)) {
    annotation <- annotate_affymetrix_probes(probe_ids, platform, organism)
    if (!is.null(annotation) && nrow(annotation) > 0) {
      return(annotation)
    }
  }
  
  # 使用org.db注释 / Annotate using org.db
  annotation <- perform_org_db_annotation(probe_ids, org_db, from_type, to_type)
  
  return(annotation)
}


#' 获取org.db数据库
#' Get org.db Database
#'
#' @param organism 物种 / Organism
#' @return org.db对象 / org.db object
get_org_db <- function(organism) {
  
  tryCatch({
    org_db <- switch(organism,
      "human" = {
        if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
          message("尝试安装org.Hs.eg.db / Trying to install org.Hs.eg.db...")
          return(NULL)
        }
        org.Hs.eg.db::org.Hs.eg.db
      },
      "mouse" = {
        if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
          return(NULL)
        }
        org.Mm.eg.db::org.Mm.eg.db
      },
      "rat" = {
        if (!requireNamespace("org.Rn.eg.db", quietly = TRUE)) {
          return(NULL)
        }
        org.Rn.eg.db::org.Rn.eg.db
      }
    )
    return(org_db)
  }, error = function(e) {
    warning(paste0("加载org.db失败 / Failed to load org.db: ", conditionMessage(e)))
    return(NULL)
  })
}


#' 检测ID类型
#' Detect ID Type
#'
#' @param ids ID向量 / ID vector
#' @return 字符串，ID类型 / Character, ID type
detect_id_type <- function(ids) {
  
  sample_ids <- head(ids, 100)
  
  # Ensembl格式 / Ensembl format
  if (mean(grepl("^ENS[A-Z]*[GT]\\d+", sample_ids)) > 0.3) {
    return("ENSEMBL")
  }
  
  # Entrez（纯数字）/ Entrez (pure numbers)
  if (mean(grepl("^\\d+$", sample_ids)) > 0.5) {
    return("ENTREZID")
  }
  
  # Affymetrix探针 / Affymetrix probes
  if (mean(grepl("_at$|_s_at$|_x_at$", sample_ids)) > 0.3) {
    return("PROBEID")
  }
  
  # RefSeq / RefSeq
  if (mean(grepl("^N[MR]_", sample_ids)) > 0.3) {
    return("REFSEQ")
  }
  
  # 默认为Symbol / Default to Symbol
  return("SYMBOL")
}


#' 执行org.db注释
#' Perform org.db Annotation
#'
#' @param ids ID向量 / ID vector
#' @param org_db org.db对象 / org.db object
#' @param from_type 输入类型 / Input type
#' @param to_type 输出类型 / Output type
#' @return 数据框，注释结果 / Data frame with annotation results
perform_org_db_annotation <- function(ids, org_db, from_type, to_type) {
  
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    return(create_empty_annotation(ids))
  }
  
  # 定义要获取的列 / Define columns to retrieve
  columns <- c("SYMBOL", "ENTREZID", "GENENAME")
  if (to_type %in% c("ENSEMBL", "REFSEQ")) {
    columns <- c(columns, to_type)
  }
  
  # 过滤可用列 / Filter available columns
  available_cols <- AnnotationDbi::columns(org_db)
  columns <- columns[columns %in% available_cols]
  
  # 如果from_type不在可用的keytype中，尝试其他方式 / If from_type not available, try alternatives
  available_keys <- AnnotationDbi::keytypes(org_db)
  
  if (!from_type %in% available_keys) {
    if ("ALIAS" %in% available_keys) {
      from_type <- "ALIAS"
    } else {
      from_type <- "SYMBOL"
    }
  }
  
  tryCatch({
    annotation <- AnnotationDbi::select(
      org_db,
      keys = ids,
      columns = columns,
      keytype = from_type
    )
    
    # 处理重复 / Handle duplicates
    annotation <- annotation[!duplicated(annotation[[from_type]]), ]
    
    # 确保有PROBEID列 / Ensure PROBEID column exists
    if (!"PROBEID" %in% colnames(annotation)) {
      annotation$PROBEID <- annotation[[from_type]]
    }
    
    return(annotation)
    
  }, error = function(e) {
    warning(paste0("org.db注释失败 / org.db annotation failed: ", conditionMessage(e)))
    return(create_empty_annotation(ids))
  })
}


#' 注释Affymetrix探针
#' Annotate Affymetrix Probes
#'
#' @param probe_ids 探针ID / Probe IDs
#' @param platform GPL平台 / GPL platform
#' @param organism 物种 / Organism
#' @return 数据框，注释结果 / Data frame with annotation results
annotate_affymetrix_probes <- function(probe_ids, platform, organism) {
  
  # 尝试加载平台特定的注释包 / Try to load platform-specific annotation package
  affy_pkg <- get_affy_annotation_package(platform)
  
  if (!is.null(affy_pkg) && requireNamespace(affy_pkg, quietly = TRUE)) {
    tryCatch({
      db <- get(affy_pkg)
      
      annotation <- AnnotationDbi::select(
        db,
        keys = probe_ids,
        columns = c("PROBEID", "SYMBOL", "ENTREZID", "GENENAME"),
        keytype = "PROBEID"
      )
      
      return(annotation)
      
    }, error = function(e) {
      message(paste0("平台注释包查询失败 / Platform annotation failed: ", conditionMessage(e)))
    })
  }
  
  return(NULL)
}


#' 获取Affymetrix注释包名
#' Get Affymetrix Annotation Package Name
#'
#' @param platform GPL平台ID / GPL platform ID
#' @return 字符串，包名 / Character, package name
get_affy_annotation_package <- function(platform) {
  
  # 常见平台到注释包的映射 / Common platform to package mapping
  pkg_map <- list(
    "GPL570" = "hgu133plus2.db",
    "GPL571" = "hgu133a2.db",
    "GPL96" = "hgu133a.db",
    "GPL97" = "hgu133b.db",
    "GPL1261" = "mouse4302.db",
    "GPL6244" = "hugene10sttranscriptcluster.db",
    "GPL6246" = "mogene10sttranscriptcluster.db"
  )
  
  return(pkg_map[[platform]])
}


#' 创建空注释表
#' Create Empty Annotation Table
#'
#' @param ids ID向量 / ID vector
#' @return 数据框，空注释表 / Data frame with empty annotation
create_empty_annotation <- function(ids) {
  
  data.frame(
    PROBEID = ids,
    SYMBOL = NA_character_,
    ENTREZID = NA_character_,
    GENENAME = NA_character_,
    stringsAsFactors = FALSE
  )
}


#' ID批量转换
#' Batch ID Conversion
#'
#' @param ids ID向量 / ID vector
#' @param from_type 输入类型 / Input type
#' @param to_type 输出类型 / Output type
#' @param organism 物种 / Organism
#'
#' @return 数据框，转换结果 / Data frame with conversion results
#' @export
convert_gene_ids <- function(ids,
                              from_type = c("SYMBOL", "ENTREZID", "ENSEMBL", "REFSEQ"),
                              to_type = c("SYMBOL", "ENTREZID", "ENSEMBL", "REFSEQ"),
                              organism = "human") {
  
  from_type <- match.arg(from_type)
  to_type <- match.arg(to_type)
  
  org_db <- get_org_db(organism)
  
  if (is.null(org_db)) {
    return(data.frame(from = ids, to = NA, stringsAsFactors = FALSE))
  }
  
  tryCatch({
    result <- AnnotationDbi::select(
      org_db,
      keys = ids,
      columns = to_type,
      keytype = from_type
    )
    
    colnames(result) <- c("from", "to")
    result <- result[!duplicated(result$from), ]
    
    return(result)
    
  }, error = function(e) {
    warning(paste0("ID转换失败 / ID conversion failed: ", conditionMessage(e)))
    return(data.frame(from = ids, to = NA, stringsAsFactors = FALSE))
  })
}


#' 自动检测物种
#' Auto-detect Organism
#'
#' @param gene_ids 基因ID向量 / Gene ID vector
#' @return 字符串，物种 / Character, organism
#' @export
detect_organism <- function(gene_ids) {
  
  # 尝试在人类数据库中查找 / Try to find in human database
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    human_symbols <- tryCatch({
      keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL")
    }, error = function(e) character(0))
    
    if (sum(gene_ids %in% human_symbols) / length(gene_ids) > 0.3) {
      message("检测到人类基因 / Human genes detected")
      return("human")
    }
  }
  
  # 尝试在小鼠数据库中查找 / Try to find in mouse database
  if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    mouse_symbols <- tryCatch({
      keys(org.Mm.eg.db::org.Mm.eg.db, keytype = "SYMBOL")
    }, error = function(e) character(0))
    
    if (sum(gene_ids %in% mouse_symbols) / length(gene_ids) > 0.3) {
      message("检测到小鼠基因 / Mouse genes detected")
      return("mouse")
    }
  }
  
  # 默认返回人类 / Default to human
  message("无法确定物种，默认使用人类 / Cannot determine organism, defaulting to human")
  return("human")
}


#' 获取基因注释信息
#' Get Gene Annotation Information
#'
#' @param gene_ids 基因ID（Symbol）/ Gene IDs (Symbol)
#' @param organism 物种 / Organism
#' @param info_types 要获取的信息类型 / Information types to retrieve
#'
#' @return 数据框，基因信息 / Data frame with gene information
#' @export
get_gene_info <- function(gene_ids,
                          organism = "human",
                          info_types = c("SYMBOL", "GENENAME", "ENTREZID", "CHR", "GENETYPE")) {
  
  org_db <- get_org_db(organism)
  
  if (is.null(org_db)) {
    return(NULL)
  }
  
  available_cols <- AnnotationDbi::columns(org_db)
  info_types <- info_types[info_types %in% available_cols]
  
  if (length(info_types) == 0) {
    info_types <- c("SYMBOL", "GENENAME")
  }
  
  tryCatch({
    annotation <- AnnotationDbi::select(
      org_db,
      keys = gene_ids,
      columns = info_types,
      keytype = "SYMBOL"
    )
    
    return(annotation)
    
  }, error = function(e) {
    warning(paste0("获取基因信息失败 / Failed to get gene info: ", conditionMessage(e)))
    return(NULL)
  })
}
