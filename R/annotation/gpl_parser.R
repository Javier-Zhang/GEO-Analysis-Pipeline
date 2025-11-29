#' ============================================================================
#' GPL平台文件解析模块
#' GPL Platform File Parser Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 解析GPL平台文件 / Parse GPL platform files
#' - 提取探针-基因对应关系 / Extract probe-gene mappings
#' - 处理多种GPL格式 / Handle multiple GPL formats
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

#' 使用GPL平台文件进行注释
#' Annotate Using GPL Platform File
#'
#' @param probe_ids 探针ID向量 / Vector of probe IDs
#' @param platform GPL平台ID / GPL platform ID
#' @param dest_dir 下载目录 / Download directory
#'
#' @return 数据框，注释结果 / Data frame with annotation results
#' @export
#'
#' @examples
#' \dontrun{
#' # 使用GPL570注释 / Annotate using GPL570
#' annotation <- annotate_with_gpl(probe_ids, "GPL570")
#' }
annotate_with_gpl <- function(probe_ids, platform, dest_dir = tempdir()) {
  
  message(paste0("使用GPL平台文件注释: ", platform, " / Annotating with GPL: ", platform))
  
  # 获取GPL信息 / Get GPL information
  gpl_table <- get_gpl_table(platform, dest_dir)
  
  if (is.null(gpl_table) || nrow(gpl_table) == 0) {
    warning("无法获取GPL信息 / Cannot get GPL information")
    return(create_empty_gpl_annotation(probe_ids))
  }
  
  # 解析GPL表格 / Parse GPL table
  annotation <- parse_gpl_table(gpl_table, probe_ids)
  
  return(annotation)
}


#' 获取GPL表格
#' Get GPL Table
#'
#' @param platform GPL平台ID / GPL platform ID
#' @param dest_dir 下载目录 / Download directory
#' @return 数据框，GPL表格 / Data frame with GPL table
get_gpl_table <- function(platform, dest_dir) {
  
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    warning("需要GEOquery包 / GEOquery package required")
    return(NULL)
  }
  
  tryCatch({
    gpl <- GEOquery::getGEO(platform, destdir = dest_dir)
    gpl_table <- GEOquery::Table(gpl)
    
    message(paste0("GPL表格包含 ", nrow(gpl_table), " 行, ", 
                   ncol(gpl_table), " 列 / GPL table has ", 
                   nrow(gpl_table), " rows, ", ncol(gpl_table), " columns"))
    
    return(gpl_table)
    
  }, error = function(e) {
    warning(paste0("获取GPL失败 / Failed to get GPL: ", conditionMessage(e)))
    return(NULL)
  })
}


#' 解析GPL表格
#' Parse GPL Table
#'
#' @param gpl_table GPL表格 / GPL table
#' @param probe_ids 探针ID / Probe IDs
#' @return 数据框，解析后的注释 / Data frame with parsed annotation
parse_gpl_table <- function(gpl_table, probe_ids) {
  
  # 标准化列名 / Standardize column names
  colnames(gpl_table) <- toupper(colnames(gpl_table))
  
  # 查找ID列 / Find ID column
  id_col <- find_id_column(gpl_table)
  
  # 查找Gene Symbol列 / Find Gene Symbol column
  symbol_col <- find_symbol_column(gpl_table)
  
  # 查找Entrez ID列 / Find Entrez ID column
  entrez_col <- find_entrez_column(gpl_table)
  
  # 查找Gene Name列 / Find Gene Name column
  name_col <- find_gene_name_column(gpl_table)
  
  message(paste0("ID列 / ID column: ", id_col))
  message(paste0("Symbol列 / Symbol column: ", symbol_col))
  
  # 创建注释表 / Create annotation table
  annotation <- data.frame(
    PROBEID = probe_ids,
    stringsAsFactors = FALSE
  )
  
  # 从GPL表格匹配 / Match from GPL table
  if (!is.null(id_col)) {
    idx <- match(probe_ids, gpl_table[[id_col]])
    
    if (!is.null(symbol_col)) {
      annotation$SYMBOL <- gpl_table[[symbol_col]][idx]
      # 清理Symbol / Clean Symbol
      annotation$SYMBOL <- clean_gene_symbol(annotation$SYMBOL)
    }
    
    if (!is.null(entrez_col)) {
      annotation$ENTREZID <- gpl_table[[entrez_col]][idx]
    }
    
    if (!is.null(name_col)) {
      annotation$GENENAME <- gpl_table[[name_col]][idx]
    }
  }
  
  return(annotation)
}


#' 查找ID列
#' Find ID Column
#'
#' @param gpl_table GPL表格 / GPL table
#' @return 字符串，列名 / Character, column name
find_id_column <- function(gpl_table) {
  
  # 可能的ID列名 / Possible ID column names
  id_patterns <- c("^ID$", "^PROBE.ID$", "^PROBEID$", "^PROBE_ID$", 
                   "^SPOT.ID$", "^SPOT_ID$", "^ID_REF$")
  
  for (pattern in id_patterns) {
    matches <- grep(pattern, colnames(gpl_table), ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) {
      return(matches[1])
    }
  }
  
  # 如果找不到，返回第一列 / If not found, return first column
  return(colnames(gpl_table)[1])
}


#' 查找Symbol列
#' Find Symbol Column
#'
#' @param gpl_table GPL表格 / GPL table
#' @return 字符串，列名 / Character, column name
find_symbol_column <- function(gpl_table) {
  
  symbol_patterns <- c("^GENE.SYMBOL$", "^GENE_SYMBOL$", "^SYMBOL$", 
                       "^GENE.NAME$", "^GENENAME$", "^GENE$",
                       "SYMBOL", "GENE.ASSIGNMENT", "ORF")
  
  for (pattern in symbol_patterns) {
    matches <- grep(pattern, colnames(gpl_table), ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) {
      # 检查列是否有有效数据 / Check if column has valid data
      col <- matches[1]
      if (sum(!is.na(gpl_table[[col]]) & gpl_table[[col]] != "") > nrow(gpl_table) * 0.1) {
        return(col)
      }
    }
  }
  
  return(NULL)
}


#' 查找Entrez ID列
#' Find Entrez ID Column
#'
#' @param gpl_table GPL表格 / GPL table
#' @return 字符串，列名 / Character, column name
find_entrez_column <- function(gpl_table) {
  
  entrez_patterns <- c("^ENTREZ", "^GENE.ID$", "^GENEID$", "^GENE_ID$",
                       "ENTREZ.GENE.ID", "ENTREZ_GENE_ID", "LOCUSLINK")
  
  for (pattern in entrez_patterns) {
    matches <- grep(pattern, colnames(gpl_table), ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) {
      return(matches[1])
    }
  }
  
  return(NULL)
}


#' 查找Gene Name列
#' Find Gene Name Column
#'
#' @param gpl_table GPL表格 / GPL table
#' @return 字符串，列名 / Character, column name
find_gene_name_column <- function(gpl_table) {
  
  name_patterns <- c("^GENE.TITLE$", "^GENE_TITLE$", "^DESCRIPTION$",
                     "^GENE.DESCRIPTION$", "^GENE_DESCRIPTION$",
                     "DEFINITION", "GENE.NAME")
  
  for (pattern in name_patterns) {
    matches <- grep(pattern, colnames(gpl_table), ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) {
      return(matches[1])
    }
  }
  
  return(NULL)
}


#' 清理基因符号
#' Clean Gene Symbol
#'
#' @param symbols 基因符号向量 / Gene symbol vector
#' @return 字符向量，清理后的符号 / Character vector of cleaned symbols
clean_gene_symbol <- function(symbols) {
  
  if (is.null(symbols)) return(NULL)
  
  # 转为字符 / Convert to character
  symbols <- as.character(symbols)
  
  # 处理分隔符分隔的多个基因 / Handle multiple genes separated by delimiters
  # 只保留第一个基因 / Keep only first gene
  symbols <- sapply(strsplit(symbols, " /// |///| // |//|;|,"), function(x) {
    if (length(x) > 0) trimws(x[1]) else NA
  })
  
  # 移除空白和NA字符串 / Remove whitespace and NA strings
  symbols[symbols == ""] <- NA
  symbols[symbols == "NA"] <- NA
  symbols[symbols == "---"] <- NA
  symbols[symbols == "N/A"] <- NA
  
  return(symbols)
}


#' 创建空GPL注释
#' Create Empty GPL Annotation
#'
#' @param probe_ids 探针ID / Probe IDs
#' @return 数据框，空注释 / Data frame with empty annotation
create_empty_gpl_annotation <- function(probe_ids) {
  
  data.frame(
    PROBEID = probe_ids,
    SYMBOL = NA_character_,
    ENTREZID = NA_character_,
    GENENAME = NA_character_,
    stringsAsFactors = FALSE
  )
}


#' 获取GPL元信息
#' Get GPL Meta Information
#'
#' @param platform GPL平台ID / GPL platform ID
#' @param dest_dir 下载目录 / Download directory
#'
#' @return 列表，GPL元信息 / List with GPL meta information
#' @export
get_gpl_info <- function(platform, dest_dir = tempdir()) {
  
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    return(NULL)
  }
  
  tryCatch({
    gpl <- GEOquery::getGEO(platform, destdir = dest_dir)
    meta <- GEOquery::Meta(gpl)
    
    info <- list(
      platform_id = platform,
      title = meta$title,
      organism = meta$organism,
      technology = meta$technology,
      manufacturer = meta$manufacturer,
      distribution = meta$distribution,
      n_probes = nrow(GEOquery::Table(gpl))
    )
    
    return(info)
    
  }, error = function(e) {
    warning(paste0("获取GPL信息失败 / Failed to get GPL info: ", conditionMessage(e)))
    return(NULL)
  })
}


#' 列出常见GPL平台
#' List Common GPL Platforms
#'
#' @return 数据框，常见平台信息 / Data frame with common platform info
#' @export
list_common_platforms <- function() {
  
  platforms <- data.frame(
    GPL = c(
      "GPL570", "GPL571", "GPL96", "GPL97", "GPL1261",
      "GPL6244", "GPL6246", "GPL6883", "GPL6884", "GPL10558",
      "GPL4133", "GPL6480", "GPL13497"
    ),
    Name = c(
      "Affymetrix HG-U133 Plus 2.0", "Affymetrix HG-U133A 2.0",
      "Affymetrix HG-U133A", "Affymetrix HG-U133B", "Affymetrix Mouse 430 2.0",
      "Affymetrix HuGene 1.0 ST", "Affymetrix MoGene 1.0 ST",
      "Illumina HumanHT-12", "Illumina HumanWG-6", "Illumina HumanHT-12 V4",
      "Agilent Human 4x44K", "Agilent Human 8x60K", "Agilent Mouse 8x60K"
    ),
    Technology = c(
      rep("Affymetrix", 7),
      rep("Illumina", 3),
      rep("Agilent", 3)
    ),
    Organism = c(
      rep("Human", 4), "Mouse", rep("Human", 2),
      rep("Human", 3),
      rep("Human", 2), "Mouse"
    ),
    stringsAsFactors = FALSE
  )
  
  return(platforms)
}


#' 检测平台类型
#' Detect Platform Type
#'
#' @param platform GPL平台ID / GPL platform ID
#' @return 字符串，平台类型 / Character, platform type
#' @export
detect_platform_type_from_gpl <- function(platform) {
  
  common <- list_common_platforms()
  
  if (platform %in% common$GPL) {
    return(common$Technology[common$GPL == platform])
  }
  
  # 根据GPL信息推断 / Infer from GPL info
  info <- get_gpl_info(platform)
  
  if (!is.null(info)) {
    if (grepl("affymetrix", tolower(info$manufacturer))) return("Affymetrix")
    if (grepl("illumina", tolower(info$manufacturer))) return("Illumina")
    if (grepl("agilent", tolower(info$manufacturer))) return("Agilent")
  }
  
  return("Unknown")
}
