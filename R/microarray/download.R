#' ============================================================================
#' GEO 微阵列数据下载模块
#' Microarray Data Download Module for GEO
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 自动下载GEO数据集（原始数据和Series Matrix）
#' - Auto-download GEO datasets (raw data and Series Matrix)
#' - 支持多种平台（Affymetrix, Illumina, Agilent等）
#' - Support multiple platforms (Affymetrix, Illumina, Agilent, etc.)
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

# 加载必要的包 / Load required packages
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(utils)
})

#' 下载GEO数据集
#' Download GEO Dataset
#'
#' @param gse_id 字符串，GEO数据集ID（如"GSE6791"）/ Character, GEO dataset ID
#' @param dest_dir 字符串，下载目标目录 / Character, destination directory
#' @param download_raw 逻辑值，是否下载原始数据 / Logical, whether to download raw data
#' @param download_matrix 逻辑值，是否下载Series Matrix / Logical, whether to download Series Matrix
#' @param download_supp 逻辑值，是否下载supplementary文件 / Logical, whether to download supplementary files
#' @param timeout 数值，下载超时时间（秒）/ Numeric, download timeout in seconds
#' 
#' @return 列表，包含下载的数据对象和文件路径 / List containing downloaded data objects and file paths
#' @export
#'
#' @examples
#' \dontrun{
#' # 下载宫颈癌相关数据集 / Download cervical cancer dataset
#' result <- download_geo_data("GSE6791", dest_dir = "data")
#' }
download_geo_data <- function(gse_id, 
                               dest_dir = "data",
                               download_raw = TRUE,
                               download_matrix = TRUE,
                               download_supp = TRUE,
                               timeout = 600) {
  
  # 参数验证 / Parameter validation
  if (!grepl("^GSE[0-9]+$", gse_id)) {
    stop("无效的GSE ID格式。请使用如 'GSE6791' 的格式。/ Invalid GSE ID format. Please use format like 'GSE6791'.")
  }
  
  # 创建目录结构 / Create directory structure
  dirs <- create_directory_structure(dest_dir, gse_id)
  
  # 设置下载超时 / Set download timeout
  old_timeout <- getOption("timeout")
  options(timeout = timeout)
  on.exit(options(timeout = old_timeout))
  
  # 初始化结果列表 / Initialize result list
  result <- list(
    gse_id = gse_id,
    dirs = dirs,
    gse_object = NULL,
    platform_info = NULL,
    raw_files = NULL,
    matrix_files = NULL,
    supp_files = NULL,
    download_log = list(),
    status = "success"
  )
  
  tryCatch({
    # 1. 下载Series Matrix / Download Series Matrix
    if (download_matrix) {
      message(paste0("\n[1/3] 正在下载 Series Matrix... / Downloading Series Matrix for ", gse_id, "..."))
      result$gse_object <- download_series_matrix(gse_id, dirs$matrix_dir)
      result$download_log$matrix <- "success"
      message("Series Matrix 下载完成！/ Series Matrix download complete!")
    }
    
    # 2. 获取平台信息 / Get platform information
    if (!is.null(result$gse_object)) {
      result$platform_info <- get_platform_info(result$gse_object)
      message(paste0("检测到平台 / Detected platform: ", result$platform_info$platform_id))
    }
    
    # 3. 下载原始数据 / Download raw data
    if (download_raw) {
      message(paste0("\n[2/3] 正在下载原始数据... / Downloading raw data for ", gse_id, "..."))
      result$raw_files <- download_raw_data(gse_id, dirs$raw_dir)
      result$download_log$raw <- "success"
      message("原始数据下载完成！/ Raw data download complete!")
    }
    
    # 4. 下载supplementary文件 / Download supplementary files
    if (download_supp) {
      message(paste0("\n[3/3] 正在下载补充文件... / Downloading supplementary files for ", gse_id, "..."))
      result$supp_files <- download_supplementary_files(gse_id, dirs$supp_dir)
      result$download_log$supp <- "success"
      message("补充文件下载完成！/ Supplementary files download complete!")
    }
    
  }, error = function(e) {
    result$status <<- "error"
    result$error_message <<- conditionMessage(e)
    warning(paste0("下载过程中出错 / Error during download: ", conditionMessage(e)))
  })
  
  # 生成下载报告 / Generate download report
  result$summary <- generate_download_summary(result)
  
  # 保存元数据 / Save metadata
  save_download_metadata(result, dirs$base_dir)
  
  return(result)
}


#' 创建目录结构
#' Create Directory Structure
#'
#' @param dest_dir 基础目录 / Base directory
#' @param gse_id GSE ID
#' @return 列表，包含各子目录路径 / List containing subdirectory paths
create_directory_structure <- function(dest_dir, gse_id) {
  
  dirs <- list(
    base_dir = file.path(dest_dir, gse_id),
    raw_dir = file.path(dest_dir, gse_id, "raw"),
    matrix_dir = file.path(dest_dir, gse_id, "matrix"),
    supp_dir = file.path(dest_dir, gse_id, "supplementary"),
    processed_dir = file.path(dest_dir, gse_id, "processed"),
    qc_dir = file.path(dest_dir, gse_id, "qc"),
    results_dir = file.path(dest_dir, gse_id, "results")
  )
  
  # 创建所有目录 / Create all directories
  for (dir_path in dirs) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
      message(paste0("创建目录 / Created directory: ", dir_path))
    }
  }
  
  return(dirs)
}


#' 下载Series Matrix
#' Download Series Matrix
#'
#' @param gse_id GSE ID
#' @param dest_dir 目标目录 / Destination directory
#'
#' @return ExpressionSet或列表 / ExpressionSet or list
download_series_matrix <- function(gse_id, dest_dir) {
  
  # 使用GEOquery下载 / Download using GEOquery
  gse <- GEOquery::getGEO(
    gse_id, 
    destdir = dest_dir,
    GSEMatrix = TRUE,
    AnnotGPL = TRUE,
    getGPL = TRUE
  )
  
  # 如果返回列表（多平台），保持原样；否则包装成列表
  # If returns list (multiple platforms), keep as is; otherwise wrap in list
  if (!is.list(gse)) {
    gse <- list(gse)
    names(gse) <- gse_id
  }
  
  return(gse)
}


#' 下载原始数据
#' Download Raw Data
#'
#' @param gse_id GSE ID
#' @param dest_dir 目标目录 / Destination directory
#'
#' @return 字符向量，下载的文件路径 / Character vector of downloaded file paths
download_raw_data <- function(gse_id, dest_dir) {
  
  downloaded_files <- character()
  
  tryCatch({
    # 获取supplementary文件URL / Get supplementary file URLs
    supp_files <- GEOquery::getGEOSuppFiles(
      gse_id,
      baseDir = dest_dir,
      makeDirectory = FALSE,
      fetch_files = TRUE
    )
    
    if (!is.null(supp_files) && nrow(supp_files) > 0) {
      downloaded_files <- rownames(supp_files)
      
      # 自动解压tar和gz文件 / Auto-extract tar and gz files
      downloaded_files <- auto_extract_files(downloaded_files, dest_dir)
    }
    
  }, error = function(e) {
    warning(paste0("原始数据下载失败 / Raw data download failed: ", conditionMessage(e)))
  })
  
  return(downloaded_files)
}


#' 下载补充文件
#' Download Supplementary Files
#'
#' @param gse_id GSE ID
#' @param dest_dir 目标目录 / Destination directory
#'
#' @return 字符向量，下载的文件路径 / Character vector of downloaded file paths
download_supplementary_files <- function(gse_id, dest_dir) {
  
  downloaded_files <- character()
  
  tryCatch({
    # 尝试获取补充文件 / Try to get supplementary files
    supp_files <- GEOquery::getGEOSuppFiles(
      gse_id,
      baseDir = dest_dir,
      makeDirectory = FALSE,
      fetch_files = TRUE
    )
    
    if (!is.null(supp_files) && nrow(supp_files) > 0) {
      downloaded_files <- rownames(supp_files)
    }
    
  }, error = function(e) {
    message(paste0("无补充文件或下载失败 / No supplementary files or download failed: ", conditionMessage(e)))
  })
  
  return(downloaded_files)
}


#' 自动解压文件
#' Auto Extract Files
#'
#' @param file_paths 文件路径向量 / Vector of file paths
#' @param dest_dir 解压目标目录 / Extraction destination directory
#'
#' @return 字符向量，所有文件路径（包括解压后的）/ Character vector of all file paths
auto_extract_files <- function(file_paths, dest_dir) {
  
  all_files <- file_paths
  
  for (file_path in file_paths) {
    if (file.exists(file_path)) {
      # 解压tar.gz文件 / Extract tar.gz files
      if (grepl("\.tar\.gz$|\.tgz$", file_path, ignore.case = TRUE)) {
        message(paste0("正在解压 / Extracting: ", basename(file_path)))
        tryCatch({
          utils::untar(file_path, exdir = dest_dir)
          extracted <- list.files(dest_dir, full.names = TRUE, recursive = TRUE)
          all_files <- c(all_files, extracted)
        }, error = function(e) {
          warning(paste0("解压失败 / Extraction failed: ", conditionMessage(e)))
        })
      }
      # 解压gz文件 / Extract gz files
      else if (grepl("\.gz$", file_path, ignore.case = TRUE) && \
               !grepl("\.tar\.gz$", file_path, ignore.case = TRUE)) {
        message(paste0("正在解压 / Extracting: ", basename(file_path)))
        tryCatch({
          out_file <- gsub("\.gz$", "", file_path)
          R.utils::gunzip(file_path, destname = out_file, remove = FALSE, overwrite = TRUE)
          all_files <- c(all_files, out_file)
        }, error = function(e) {
          # 如果R.utils不可用，尝试base R方法
          tryCatch({
            con <- gzfile(file_path, "rb")
            out_file <- gsub("\.gz$", "", file_path)
            writeLines(readLines(con), out_file)
            close(con)
            all_files <- c(all_files, out_file)
          }, error = function(e2) {
            warning(paste0("解压失败 / Extraction failed: ", conditionMessage(e2)))
          })
        })
      }
    }
  }
  
  return(unique(all_files))
}


#' 获取平台信息
#' Get Platform Information
#'
#' @param gse_object GSE对象 / GSE object
#'
#' @return 列表，包含平台详细信息 / List containing platform details
get_platform_info <- function(gse_object) {
  
  # 处理多平台情况 / Handle multiple platforms
  if (is.list(gse_object)) {
    gse <- gse_object[[1]]
  } else {
    gse <- gse_object
  }
  
  platform_info <- list(
    platform_id = annotation(gse),
    platform_type = detect_platform_type(annotation(gse)),
    n_samples = ncol(exprs(gse)),
    n_features = nrow(exprs(gse)),
    sample_ids = sampleNames(gse)
  )
  
  # 获取更多平台详情 / Get more platform details
  tryCatch({
    gpl <- GEOquery::getGEO(platform_info$platform_id)
    platform_info$platform_title <- Meta(gpl)$title
    platform_info$platform_organism <- Meta(gpl)$organism
    platform_info$platform_technology <- Meta(gpl)$technology
  }, error = function(e) {
    platform_info$platform_title <- "Unknown"
  })
  
  return(platform_info)
}


#' 检测平台类型
#' Detect Platform Type
#'
#' @param platform_id GPL ID
#' @return 字符串，平台类型 / Character, platform type
detect_platform_type <- function(platform_id) {
  
  # 常见平台对照表 / Common platform mapping
  affymetrix_platforms <- c(
    "GPL570", "GPL571", "GPL96", "GPL97", "GPL571", "GPL1261",
    "GPL6244", "GPL6246", "GPL17586", "GPL17692", "GPL16686"
  )
  
  illumina_platforms <- c(
    "GPL6883", "GPL6884", "GPL6947", "GPL10558", "GPL6102",
    "GPL6104", "GPL6887", "GPL13376", "GPL14951"
  )
  
  agilent_platforms <- c(
    "GPL4133", "GPL6480", "GPL13497", "GPL14550", "GPL17077",
    "GPL6848", "GPL7264"
  )
  
  if (platform_id %in% affymetrix_platforms) {
    return("Affymetrix")
  } else if (platform_id %in% illumina_platforms) {
    return("Illumina")
  } else if (platform_id %in% agilent_platforms) {
    return("Agilent")
  } else if (grepl("^GPL", platform_id)) {
    # 尝试从GPL信息推断 / Try to infer from GPL info
    return("Unknown")
  } else {
    return("Unknown")
  }
}


#' 生成下载摘要
#' Generate Download Summary
#'
#' @param result 下载结果对象 / Download result object
#'
#' @return 列表，下载摘要 / List, download summary
generate_download_summary <- function(result) {
  
  summary <- list(
    gse_id = result$gse_id,
    download_time = Sys.time(),
    status = result$status
  )
  
  if (!is.null(result$platform_info)) {
    summary$platform <- result$platform_info$platform_id
    summary$platform_type <- result$platform_info$platform_type
    summary$n_samples <- result$platform_info$n_samples
    summary$n_features <- result$platform_info$n_features
  }
  
  summary$raw_files_count <- length(result$raw_files)
  summary$supp_files_count <- length(result$supp_files)
  
  return(summary)
}


#' 保存下载元数据
#' Save Download Metadata
#'
#' @param result 下载结果对象 / Download result object
#' @param dest_dir 目标目录 / Destination directory
save_download_metadata <- function(result, dest_dir) {
  
  metadata <- list(
    gse_id = result$gse_id,
    download_time = as.character(Sys.time()),
    status = result$status,
    platform_info = result$platform_info,
    raw_files = result$raw_files,
    supp_files = result$supp_files,
    summary = result$summary,
    r_version = R.version.string,
    geo_query_version = as.character(packageVersion("GEOquery"))
  )
  
  # 保存为JSON格式 / Save as JSON format
  metadata_file <- file.path(dest_dir, "download_metadata.json")
  tryCatch({
    jsonlite::write_json(metadata, metadata_file, pretty = TRUE, auto_unbox = TRUE)
    message(paste0("元数据已保存 / Metadata saved to: ", metadata_file))
  }, error = function(e) {
    # 如果jsonlite不可用，保存为RDS / If jsonlite unavailable, save as RDS
    rds_file <- file.path(dest_dir, "download_metadata.rds")
    saveRDS(metadata, rds_file)
    message(paste0("元数据已保存 / Metadata saved to: ", rds_file))
  })
}


#' 检查下载状态
#' Check Download Status
#'
#' @param gse_id GSE ID
#' @param dest_dir 数据目录 / Data directory
#'
#' @return 列表，下载状态信息 / List, download status information
check_download_status <- function(gse_id, dest_dir = "data") {
  
  base_dir <- file.path(dest_dir, gse_id)
  
  status <- list(
    exists = dir.exists(base_dir),
    has_matrix = FALSE,
    has_raw = FALSE,
    has_supp = FALSE,
    metadata = NULL
  )
  
  if (status$exists) {
    # 检查各子目录 / Check subdirectories
    status$has_matrix <- length(list.files(file.path(base_dir, "matrix"))) > 0
    status$has_raw <- length(list.files(file.path(base_dir, "raw"))) > 0
    status$has_supp <- length(list.files(file.path(base_dir, "supplementary"))) > 0
    
    # 读取元数据 / Read metadata
    metadata_file <- file.path(base_dir, "download_metadata.json")
    rds_file <- file.path(base_dir, "download_metadata.rds")
    
    if (file.exists(metadata_file)) {
      status$metadata <- jsonlite::read_json(metadata_file)
    } else if (file.exists(rds_file)) {
      status$metadata <- readRDS(rds_file)
    }
  }
  
  return(status)
}


#' 批量下载多个GEO数据集
#' Batch Download Multiple GEO Datasets
#'
#' @param gse_ids 字符向量，GSE ID列表 / Character vector of GSE IDs
#' @param dest_dir 目标目录 / Destination directory
#' @param ... 其他参数传递给download_geo_data / Other parameters passed to download_geo_data
#'
#' @return 列表，每个数据集的下载结果 / List of download results for each dataset
#' @export
batch_download_geo <- function(gse_ids, dest_dir = "data", ...) {
  
  results <- list()
  
  for (i in seq_along(gse_ids)) {
    gse_id <- gse_ids[i]
    message(paste0("\n========== [", i, "/", length(gse_ids), "] 下载 / Downloading: ", gse_id, " =========="))
    
    results[[gse_id]] <- tryCatch({
      download_geo_data(gse_id, dest_dir = dest_dir, ...)
    }, error = function(e) {
      list(
        gse_id = gse_id,
        status = "error",
        error_message = conditionMessage(e)
      )
    })
  }
  
  # 生成批量下载报告 / Generate batch download report
  success_count <- sum(sapply(results, function(x) x$status == "success"))
  message(paste0("\n========== 批量下载完成 / Batch download complete =========="))
  message(paste0("成功 / Success: ", success_count, "/", length(gse_ids)))
  
  return(results)
}


#' 宫颈癌相关GEO数据集推荐列表
#' Recommended Cervical Cancer GEO Datasets
#'
#' @return 数据框，包含推荐的数据集信息 / Data frame with recommended dataset information
#' @export
get_cervical_cancer_datasets <- function() {
  
  datasets <- data.frame(
    GSE_ID = c(
      "GSE6791",    # 经典宫颈癌数据集
      "GSE7803",    # 宫颈癌 vs 正常
      "GSE9750",    # 宫颈癌基因表达
      "GSE26511",   # 宫颈癌亚型
      "GSE39001",   # HPV相关宫颈癌
      "GSE63514",   # 宫颈癌进展
      "GSE67522"    # 宫颈鳞癌
    ),
    Platform = c(
      "GPL570",
      "GPL570", 
      "GPL96",
      "GPL570",
      "GPL570",
      "GPL570",
      "GPL570"
    ),
    Samples = c(42, 28, 66, 39, 108, 128, 82),
    Description = c(
      "Cervical cancer vs Normal (经典数据集)",
      "Cervical cancer expression profiles",
      "Cervical cancer gene expression",
      "Cervical cancer subtypes",
      "HPV-associated cervical cancer",
      "Cervical cancer progression",
      "Cervical squamous cell carcinoma"
    ),
    stringsAsFactors = FALSE
  )
  
  return(datasets)
}