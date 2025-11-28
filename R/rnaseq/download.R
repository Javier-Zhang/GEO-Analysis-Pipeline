#' ============================================================================
#' GEO RNA-seq数据下载模块
#' RNA-seq Data Download Module for GEO
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 下载GEO RNA-seq数据集 / Download GEO RNA-seq datasets
#' - 获取counts matrix / Get counts matrix
#' - 支持supplementary files / Support supplementary files
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(utils)
})

#' 下载GEO RNA-seq数据
#' Download GEO RNA-seq Data
#'
#' @param gse_id 字符串，GEO数据集ID / Character, GEO dataset ID
#' @param dest_dir 字符串，下载目标目录 / Character, destination directory
#' @param download_supp 逻辑值，是否下载supplementary文件 / Logical, whether to download supplementary files
#' @param timeout 数值，下载超时时间（秒）/ Numeric, download timeout in seconds
#'
#' @return 列表，包含下载的数据对象和文件路径 / List containing downloaded data and file paths
#' @export
#'
#' @examples
#' \dontrun{
#' # 下载RNA-seq数据集 / Download RNA-seq dataset
#' result <- download_rnaseq_data("GSE63514", dest_dir = "data")
#' }
download_rnaseq_data <- function(gse_id,
                                  dest_dir = "data",
                                  download_supp = TRUE,
                                  timeout = 600) {
  
  # 参数验证 / Parameter validation
  if (!grepl("^GSE[0-9]+$", gse_id)) {
    stop("无效的GSE ID格式 / Invalid GSE ID format")
  }
  
  message(paste0("开始下载RNA-seq数据: ", gse_id, " / Downloading RNA-seq data..."))
  
  # 创建目录结构 / Create directory structure
  dirs <- create_rnaseq_directory_structure(dest_dir, gse_id)
  
  # 设置超时 / Set timeout
  old_timeout <- getOption("timeout")
  options(timeout = timeout)
  on.exit(options(timeout = old_timeout))
  
  result <- list(
    gse_id = gse_id,
    dirs = dirs,
    gse_object = NULL,
    counts_matrix = NULL,
    phenotype_data = NULL,
    supp_files = NULL,
    status = "success"
  )
  
  tryCatch({
    # 1. 下载Series Matrix / Download Series Matrix
    message("[1/3] 下载Series Matrix... / Downloading Series Matrix...")
    result$gse_object <- GEOquery::getGEO(
      gse_id,
      destdir = dirs$matrix_dir,
      GSEMatrix = TRUE,
      AnnotGPL = TRUE
    )
    
    # 2. 提取表型数据 / Extract phenotype data
    message("[2/3] 提取表型数据... / Extracting phenotype data...")
    if (is.list(result$gse_object)) {
      result$phenotype_data <- pData(result$gse_object[[1]])
    } else {
      result$phenotype_data <- pData(result$gse_object)
    }
    
    # 3. 下载supplementary files / Download supplementary files
    if (download_supp) {
      message("[3/3] 下载supplementary文件... / Downloading supplementary files...")
      result$supp_files <- download_rnaseq_supp_files(gse_id, dirs$supp_dir)
      
      # 尝试解析counts matrix / Try to parse counts matrix
      if (length(result$supp_files) > 0) {
        result$counts_matrix <- parse_counts_from_supp(result$supp_files)
      }
    }
    
    message("RNA-seq数据下载完成！/ RNA-seq data download complete!")
    
  }, error = function(e) {
    result$status <<- "error"
    result$error_message <<- conditionMessage(e)
    warning(paste0("下载出错 / Download error: ", conditionMessage(e)))
  })
  
  # 保存元数据 / Save metadata
  save_rnaseq_metadata(result, dirs$base_dir)
  
  return(result)
}


#' 创建RNA-seq目录结构
#' Create RNA-seq Directory Structure
#'
#' @param dest_dir 基础目录 / Base directory
#' @param gse_id GSE ID
#' @return 列表，目录路径 / List of directory paths
create_rnaseq_directory_structure <- function(dest_dir, gse_id) {
  
  dirs <- list(
    base_dir = file.path(dest_dir, gse_id),
    matrix_dir = file.path(dest_dir, gse_id, "matrix"),
    supp_dir = file.path(dest_dir, gse_id, "supplementary"),
    counts_dir = file.path(dest_dir, gse_id, "counts"),
    processed_dir = file.path(dest_dir, gse_id, "processed"),
    qc_dir = file.path(dest_dir, gse_id, "qc"),
    results_dir = file.path(dest_dir, gse_id, "results")
  )
  
  for (dir_path in dirs) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
      message(paste0("创建目录 / Created: ", dir_path))
    }
  }
  
  return(dirs)
}


#' 下载RNA-seq supplementary文件
#' Download RNA-seq Supplementary Files
#'
#' @param gse_id GSE ID
#' @param dest_dir 目标目录 / Destination directory
#' @return 字符向量，文件路径 / Character vector of file paths
download_rnaseq_supp_files <- function(gse_id, dest_dir) {
  
  files <- character()
  
  tryCatch({
    supp <- GEOquery::getGEOSuppFiles(
      gse_id,
      baseDir = dest_dir,
      makeDirectory = FALSE,
      fetch_files = TRUE
    )
    
    if (!is.null(supp) && nrow(supp) > 0) {
      files <- rownames(supp)
      
      # 解压文件 / Extract files
      for (f in files) {
        if (grepl("\\.gz$", f) && file.exists(f)) {
          message(paste0("解压 / Extracting: ", basename(f)))
          tryCatch({
            R.utils::gunzip(f, remove = FALSE, overwrite = TRUE)
          }, error = function(e) NULL)
        }
      }
    }
  }, error = function(e) {
    message(paste0("Supplementary下载失败 / Supp download failed: ", conditionMessage(e)))
  })
  
  return(files)
}


#' 从supplementary文件解析counts矩阵
#' Parse Counts Matrix from Supplementary Files
#'
#' @param supp_files supplementary文件路径 / Supplementary file paths
#' @return 矩阵，counts数据 / Counts matrix
parse_counts_from_supp <- function(supp_files) {
  
  counts <- NULL
  
  # 查找可能的counts文件 / Find potential counts files
  counts_patterns <- c("counts", "raw", "htseq", "featurecounts", "rsem", "kallisto")
  
  for (f in supp_files) {
    f_lower <- tolower(basename(f))
    
    # 检查是否可能是counts文件 / Check if potentially a counts file
    is_counts <- any(sapply(counts_patterns, function(p) grepl(p, f_lower)))
    is_matrix <- grepl("\\.txt$|\\.tsv$|\\.csv$|\\.tab$", f_lower) ||
                 grepl("\\.txt\\.gz$|\\.tsv\\.gz$|\\.csv\\.gz$", tolower(f))
    
    if ((is_counts || is_matrix) && file.exists(f)) {
      message(paste0("尝试解析counts文件 / Trying to parse: ", basename(f)))
      
      tryCatch({
        # 尝试读取文件 / Try to read file
        if (grepl("\\.csv", f_lower)) {
          df <- read.csv(f, row.names = 1, check.names = FALSE)
        } else {
          df <- read.delim(f, row.names = 1, check.names = FALSE)
        }
        
        # 检查是否为counts数据（整数值）/ Check if counts data (integer values)
        if (is.numeric(as.matrix(df))) {
          counts <- as.matrix(df)
          message(paste0("成功解析counts矩阵 / Successfully parsed counts matrix: ",
                        nrow(counts), " genes x ", ncol(counts), " samples"))
          break
        }
      }, error = function(e) NULL)
    }
  }
  
  return(counts)
}


#' 保存RNA-seq元数据
#' Save RNA-seq Metadata
#'
#' @param result 下载结果 / Download result
#' @param dest_dir 目标目录 / Destination directory
save_rnaseq_metadata <- function(result, dest_dir) {
  
  metadata <- list(
    gse_id = result$gse_id,
    download_time = as.character(Sys.time()),
    status = result$status,
    has_counts = !is.null(result$counts_matrix),
    n_samples = if (!is.null(result$phenotype_data)) nrow(result$phenotype_data) else NA,
    supp_files = result$supp_files,
    r_version = R.version.string
  )
  
  metadata_file <- file.path(dest_dir, "rnaseq_metadata.rds")
  saveRDS(metadata, metadata_file)
  message(paste0("元数据已保存 / Metadata saved: ", metadata_file))
}


#' 从GEO获取counts数据
#' Get Counts Data from GEO
#'
#' @param gse_id GSE ID
#' @param dest_dir 目标目录 / Destination directory
#'
#' @return 矩阵，counts数据 / Counts matrix
#' @export
get_geo_counts <- function(gse_id, dest_dir = "data") {
  
  # 检查是否已下载 / Check if already downloaded
  counts_file <- file.path(dest_dir, gse_id, "counts", "counts_matrix.rds")
  
  if (file.exists(counts_file)) {
    message("加载已有counts数据 / Loading existing counts data...")
    return(readRDS(counts_file))
  }
  
  # 下载数据 / Download data
  result <- download_rnaseq_data(gse_id, dest_dir)
  
  if (!is.null(result$counts_matrix)) {
    # 保存counts / Save counts
    saveRDS(result$counts_matrix, counts_file)
    return(result$counts_matrix)
  }
  
  warning("无法获取counts数据 / Cannot get counts data")
  return(NULL)
}


#' 检查RNA-seq数据下载状态
#' Check RNA-seq Data Download Status
#'
#' @param gse_id GSE ID
#' @param dest_dir 数据目录 / Data directory
#'
#' @return 列表，状态信息 / List with status information
#' @export
check_rnaseq_download <- function(gse_id, dest_dir = "data") {
  
  base_dir <- file.path(dest_dir, gse_id)
  
  status <- list(
    exists = dir.exists(base_dir),
    has_matrix = FALSE,
    has_counts = FALSE,
    has_supp = FALSE,
    metadata = NULL
  )
  
  if (status$exists) {
    status$has_matrix <- length(list.files(file.path(base_dir, "matrix"))) > 0
    status$has_counts <- file.exists(file.path(base_dir, "counts", "counts_matrix.rds"))
    status$has_supp <- length(list.files(file.path(base_dir, "supplementary"))) > 0
    
    meta_file <- file.path(base_dir, "rnaseq_metadata.rds")
    if (file.exists(meta_file)) {
      status$metadata <- readRDS(meta_file)
    }
  }
  
  return(status)
}


#' 获取推荐的RNA-seq数据集
#' Get Recommended RNA-seq Datasets
#'
#' @return 数据框，推荐数据集 / Data frame with recommended datasets
#' @export
get_rnaseq_example_datasets <- function() {
  
  datasets <- data.frame(
    GSE_ID = c(
      "GSE63514",
      "GSE89657",
      "GSE107399"
    ),
    Description = c(
      "Cervical cancer RNA-seq / 宫颈癌RNA-seq",
      "Cancer transcriptome study / 癌症转录组研究",
      "Gene expression profiling / 基因表达谱分析"
    ),
    Platform = c("Illumina", "Illumina", "Illumina"),
    stringsAsFactors = FALSE
  )
  
  return(datasets)
}
