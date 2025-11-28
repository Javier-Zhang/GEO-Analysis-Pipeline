#' ============================================================================
#' 表型文件标准化模块
#' Phenotype File Standardization Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 自动解析GSM characteristics / Auto-parse GSM characteristics
#' - 标准化临床信息格式 / Standardize clinical information format
#' - 生成统一表型数据框 / Generate unified phenotype data frame
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

#' 解析和标准化表型数据
#' Parse and Standardize Phenotype Data
#'
#' @param eset ExpressionSet对象或表型数据框 / ExpressionSet or phenotype data frame
#'
#' @return 数据框，标准化的表型数据 / Data frame with standardized phenotype
#' @export
#'
#' @examples
#' \dontrun{
#' # 从ExpressionSet提取并标准化表型 / Extract and standardize from ExpressionSet
#' phenotype <- standardize_phenotype(eset)
#' }
standardize_phenotype <- function(eset) {
  
  message("解析表型数据... / Parsing phenotype data...")
  
  # 获取表型数据 / Get phenotype data
  if (inherits(eset, "ExpressionSet")) {
    pdata <- pData(eset)
  } else if (is.data.frame(eset)) {
    pdata <- eset
  } else {
    stop("输入必须是ExpressionSet或数据框 / Input must be ExpressionSet or data frame")
  }
  
  # 解析characteristics列 / Parse characteristics columns
  parsed_chars <- parse_characteristics(pdata)
  
  # 合并解析结果 / Merge parsed results
  standardized <- merge_phenotype(pdata, parsed_chars)
  
  # 清理列名 / Clean column names
  standardized <- clean_phenotype_columns(standardized)
  
  message(paste0("解析完成，共 ", ncol(standardized), " 列 / ",
                 "Parsing complete, ", ncol(standardized), " columns"))
  
  return(standardized)
}


#' 解析characteristics列
#' Parse Characteristics Columns
#'
#' @param pdata 表型数据框 / Phenotype data frame
#' @return 数据框，解析后的characteristics / Data frame with parsed characteristics
parse_characteristics <- function(pdata) {
  
  # 查找characteristics列 / Find characteristics columns
  char_cols <- grep("characteristics", colnames(pdata), ignore.case = TRUE, value = TRUE)
  
  if (length(char_cols) == 0) {
    message("未找到characteristics列 / No characteristics columns found")
    return(NULL)
  }
  
  parsed <- data.frame(row.names = rownames(pdata))
  
  for (col in char_cols) {
    col_data <- as.character(pdata[[col]])
    
    # 尝试解析"key: value"格式 / Try to parse "key: value" format
    for (i in seq_along(col_data)) {
      value <- col_data[i]
      
      if (grepl(":", value)) {
        parts <- strsplit(value, ":\\s*")[[1]]
        if (length(parts) >= 2) {
          key <- clean_column_name(parts[1])
          val <- trimws(paste(parts[-1], collapse = ":"))
          
          if (!key %in% colnames(parsed)) {
            parsed[[key]] <- NA_character_
          }
          parsed[i, key] <- val
        }
      }
    }
  }
  
  return(parsed)
}


#' 合并表型数据
#' Merge Phenotype Data
#'
#' @param pdata 原始表型数据 / Original phenotype data
#' @param parsed_chars 解析的characteristics / Parsed characteristics
#' @return 数据框，合并后的表型 / Data frame with merged phenotype
merge_phenotype <- function(pdata, parsed_chars) {
  
  # 保留关键列 / Keep key columns
  key_patterns <- c("title", "geo_accession", "source", "organism", 
                    "platform", "description", "data_processing")
  
  keep_cols <- character()
  for (pattern in key_patterns) {
    matches <- grep(pattern, colnames(pdata), ignore.case = TRUE, value = TRUE)
    keep_cols <- c(keep_cols, matches)
  }
  
  keep_cols <- unique(keep_cols)
  
  # 创建基础数据框 / Create base data frame
  if (length(keep_cols) > 0) {
    result <- pdata[, keep_cols, drop = FALSE]
  } else {
    result <- data.frame(row.names = rownames(pdata))
  }
  
  # 添加解析的characteristics / Add parsed characteristics
  if (!is.null(parsed_chars) && ncol(parsed_chars) > 0) {
    for (col in colnames(parsed_chars)) {
      if (!col %in% colnames(result)) {
        result[[col]] <- parsed_chars[[col]]
      }
    }
  }
  
  return(result)
}


#' 清理列名
#' Clean Column Names
#'
#' @param name 原始名称 / Original name
#' @return 字符串，清理后的名称 / Character, cleaned name
clean_column_name <- function(name) {
  
  # 转小写 / Convert to lowercase
  name <- tolower(name)
  
  # 替换特殊字符 / Replace special characters
  name <- gsub("[^a-z0-9]", "_", name)
  
  # 移除多余下划线 / Remove extra underscores
  name <- gsub("_+", "_", name)
  name <- gsub("^_|_$", "", name)
  
  return(name)
}


#' 清理表型数据列
#' Clean Phenotype Data Columns
#'
#' @param pdata 表型数据 / Phenotype data
#' @return 数据框，清理后的表型 / Data frame with cleaned phenotype
clean_phenotype_columns <- function(pdata) {
  
  # 清理列名 / Clean column names
  colnames(pdata) <- sapply(colnames(pdata), clean_column_name)
  
  # 移除完全NA的列 / Remove all-NA columns
  na_cols <- sapply(pdata, function(x) all(is.na(x)))
  pdata <- pdata[, !na_cols, drop = FALSE]
  
  # 移除只有一个唯一值的列（除了关键列）/ Remove columns with single unique value
  # key_cols <- c("geo_accession", "title")
  # for (col in setdiff(colnames(pdata), key_cols)) {
  #   if (length(unique(na.omit(pdata[[col]]))) <= 1) {
  #     pdata[[col]] <- NULL
  #   }
  # }
  
  return(pdata)
}


#' 提取分组信息
#' Extract Group Information
#'
#' @param pdata 表型数据框 / Phenotype data frame
#' @param group_pattern 分组列名模式 / Group column name pattern
#'
#' @return 因子，分组变量 / Factor, group variable
#' @export
extract_group <- function(pdata, group_pattern = NULL) {
  
  # 常见的分组列名 / Common group column names
  group_patterns <- c("disease", "condition", "group", "treatment", 
                      "tissue", "cell_type", "status", "type",
                      "sample_type", "diagnosis")
  
  if (!is.null(group_pattern)) {
    group_patterns <- c(group_pattern, group_patterns)
  }
  
  # 查找分组列 / Find group column
  for (pattern in group_patterns) {
    matches <- grep(pattern, colnames(pdata), ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) {
      group_col <- matches[1]
      group <- factor(pdata[[group_col]])
      
      message(paste0("使用分组列 / Using group column: ", group_col))
      message(paste0("分组水平 / Group levels: ", paste(levels(group), collapse = ", ")))
      
      return(group)
    }
  }
  
  warning("未找到分组列 / No group column found")
  return(NULL)
}


#' 创建设计信息
#' Create Design Information
#'
#' @param pdata 表型数据框 / Phenotype data frame
#' @param group_col 分组列名 / Group column name
#' @param batch_col 批次列名（可选）/ Batch column name (optional)
#'
#' @return 列表，设计信息 / List with design information
#' @export
create_design_info <- function(pdata, group_col, batch_col = NULL) {
  
  design_info <- list(
    n_samples = nrow(pdata),
    group_col = group_col,
    batch_col = batch_col
  )
  
  # 分组信息 / Group information
  if (group_col %in% colnames(pdata)) {
    group <- factor(pdata[[group_col]])
    design_info$group <- group
    design_info$group_levels <- levels(group)
    design_info$group_sizes <- table(group)
    design_info$n_groups <- length(levels(group))
  }
  
  # 批次信息 / Batch information
  if (!is.null(batch_col) && batch_col %in% colnames(pdata)) {
    batch <- factor(pdata[[batch_col]])
    design_info$batch <- batch
    design_info$batch_levels <- levels(batch)
    design_info$n_batches <- length(levels(batch))
  }
  
  return(design_info)
}


#' 生成样本摘要表
#' Generate Sample Summary Table
#'
#' @param pdata 表型数据框 / Phenotype data frame
#' @param group_cols 要摘要的列 / Columns to summarize
#'
#' @return 数据框，摘要表 / Data frame with summary
#' @export
generate_sample_summary <- function(pdata, group_cols = NULL) {
  
  if (is.null(group_cols)) {
    # 自动检测分类变量 / Auto-detect categorical variables
    group_cols <- sapply(pdata, function(x) {
      is.factor(x) || is.character(x) && length(unique(x)) < nrow(pdata) / 2
    })
    group_cols <- names(group_cols)[group_cols]
  }
  
  summary_list <- list()
  
  for (col in group_cols) {
    if (col %in% colnames(pdata)) {
      summary_list[[col]] <- table(pdata[[col]], useNA = "ifany")
    }
  }
  
  return(summary_list)
}


#' 导出表型数据
#' Export Phenotype Data
#'
#' @param pdata 表型数据框 / Phenotype data frame
#' @param output_file 输出文件 / Output file
#' @param format 输出格式 / Output format
#'
#' @return 无返回值 / No return value
#' @export
export_phenotype <- function(pdata, output_file, format = c("csv", "tsv", "xlsx")) {
  
  format <- match.arg(format)
  
  switch(format,
    "csv" = {
      write.csv(pdata, output_file, row.names = TRUE)
    },
    "tsv" = {
      write.table(pdata, output_file, sep = "\t", row.names = TRUE, quote = FALSE)
    },
    "xlsx" = {
      if (requireNamespace("openxlsx", quietly = TRUE)) {
        openxlsx::write.xlsx(pdata, output_file, rowNames = TRUE)
      } else {
        write.csv(pdata, gsub("\\.xlsx$", ".csv", output_file), row.names = TRUE)
      }
    }
  )
  
  message(paste0("表型数据已导出 / Phenotype exported: ", output_file))
}


#' 检测和处理缺失值
#' Detect and Handle Missing Values
#'
#' @param pdata 表型数据框 / Phenotype data frame
#'
#' @return 列表，缺失值信息 / List with missing value information
#' @export
check_missing_values <- function(pdata) {
  
  missing_info <- list(
    total_missing = sum(is.na(pdata)),
    columns_with_missing = character(),
    missing_by_column = list()
  )
  
  for (col in colnames(pdata)) {
    n_missing <- sum(is.na(pdata[[col]]))
    if (n_missing > 0) {
      missing_info$columns_with_missing <- c(missing_info$columns_with_missing, col)
      missing_info$missing_by_column[[col]] <- list(
        n_missing = n_missing,
        pct_missing = round(n_missing / nrow(pdata) * 100, 1)
      )
    }
  }
  
  return(missing_info)
}
