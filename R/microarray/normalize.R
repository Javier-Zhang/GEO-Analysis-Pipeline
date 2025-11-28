#' ============================================================================
#' GEO 微阵列数据标准化模块
#' Microarray Data Normalization Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 支持多种标准化方法（RMA, GCRMA, MAS5, VSN, Quantile）
#' - Support multiple normalization methods (RMA, GCRMA, MAS5, VSN, Quantile)
#' - 背景校正选项 / Background correction options
#' - ComBat批次效应校正 / ComBat batch effect correction
#' - 标准化前后对比可视化 / Before/after normalization visualization
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

suppressPackageStartupMessages({
  library(Biobase)
  library(affy)
  library(limma)
  library(ggplot2)
})

#' 执行微阵列数据标准化
#' Perform Microarray Data Normalization
#'
#' @param eset ExpressionSet对象或AffyBatch对象 / ExpressionSet or AffyBatch object
#' @param method 标准化方法 / Normalization method
#'   可选: "rma", "gcrma", "mas5", "quantile", "vsn"
#' @param background_correct 是否进行背景校正 / Whether to perform background correction
#' @param batch 批次变量向量（用于ComBat）/ Batch variable vector (for ComBat)
#' @param output_dir 输出目录 / Output directory
#'
#' @return ExpressionSet，标准化后的数据 / Normalized ExpressionSet
#' @export
#'
#' @examples
#' \dontrun{
#' # 使用RMA标准化 / Normalize with RMA
#' eset_norm <- normalize_microarray(eset, method = "rma")
#' 
#' # 使用ComBat批次校正 / With ComBat batch correction
#' eset_norm <- normalize_microarray(eset, method = "quantile", batch = pData(eset)$batch)
#' }
normalize_microarray <- function(eset, 
                                  method = c("rma", "gcrma", "mas5", "quantile", "vsn"),
                                  background_correct = TRUE,
                                  batch = NULL,
                                  output_dir = NULL) {
  
  method <- match.arg(method)
  
  message(paste0("开始 ", toupper(method), " 标准化... / Starting ", toupper(method), " normalization..."))
  
  # 根据输入类型选择处理方式 / Choose processing based on input type
  if (inherits(eset, "AffyBatch")) {
    eset_norm <- normalize_affybatch(eset, method, background_correct)
  } else if (inherits(eset, "ExpressionSet")) {
    eset_norm <- normalize_expressionset(eset, method, background_correct)
  } else {
    stop("输入必须是ExpressionSet或AffyBatch对象 / Input must be ExpressionSet or AffyBatch object")
  }
  
  # ComBat批次效应校正 / ComBat batch effect correction
  if (!is.null(batch)) {
    message("执行ComBat批次效应校正... / Performing ComBat batch correction...")
    eset_norm <- combat_batch_correction(eset_norm, batch)
  }
  
  # 生成对比图 / Generate comparison plots
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    generate_normalization_plots(eset, eset_norm, method, output_dir)
  }
  
  message(paste0(toupper(method), " 标准化完成！/ ", toupper(method), " normalization complete!"))
  
  return(eset_norm)
}


#' AffyBatch标准化
#' Normalize AffyBatch Object
#'
#' @param affy_batch AffyBatch对象 / AffyBatch object
#' @param method 标准化方法 / Normalization method
#' @param background_correct 背景校正 / Background correction
#'
#' @return ExpressionSet对象 / ExpressionSet object
normalize_affybatch <- function(affy_batch, method, background_correct) {
  
  tryCatch({
    eset_norm <- switch(method,
      "rma" = {
        affy::rma(affy_batch, background = background_correct)
      },
      "gcrma" = {
        if (requireNamespace("gcrma", quietly = TRUE)) {
          gcrma::gcrma(affy_batch)
        } else {
          warning("gcrma包未安装，使用RMA代替 / gcrma package not installed, using RMA instead")
          affy::rma(affy_batch, background = background_correct)
        }
      },
      "mas5" = {
        affy::mas5(affy_batch)
      },
      "quantile" = {
        # 先RMA背景校正，再quantile标准化
        affy::rma(affy_batch, background = background_correct, normalize = TRUE)
      },
      "vsn" = {
        if (requireNamespace("vsn", quietly = TRUE)) {
          vsn::vsnrma(affy_batch)
        } else {
          warning("vsn包未安装，使用RMA代替 / vsn package not installed, using RMA instead")
          affy::rma(affy_batch, background = background_correct)
        }
      }
    )
    return(eset_norm)
  }, error = function(e) {
    stop(paste0("标准化失败 / Normalization failed: ", conditionMessage(e)))
  })
}


#' ExpressionSet标准化
#' Normalize ExpressionSet Object
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param method 标准化方法 / Normalization method
#' @param background_correct 背景校正 / Background correction
#'
#' @return ExpressionSet对象 / ExpressionSet object
normalize_expressionset <- function(eset, method, background_correct) {
  
  expr_matrix <- exprs(eset)
  
  tryCatch({
    norm_matrix <- switch(method,
      "rma" = {
        # 对于已处理的数据，使用log2和quantile
        if (max(expr_matrix, na.rm = TRUE) > 100) {
          expr_matrix <- log2(expr_matrix + 1)
        }
        limma::normalizeBetweenArrays(expr_matrix, method = "quantile")
      },
      "gcrma" = {
        # 对于ExpressionSet，等同于quantile
        if (max(expr_matrix, na.rm = TRUE) > 100) {
          expr_matrix <- log2(expr_matrix + 1)
        }
        limma::normalizeBetweenArrays(expr_matrix, method = "quantile")
      },
      "mas5" = {
        # MAS5缩放
        limma::normalizeBetweenArrays(expr_matrix, method = "scale")
      },
      "quantile" = {
        limma::normalizeBetweenArrays(expr_matrix, method = "quantile")
      },
      "vsn" = {
        if (requireNamespace("vsn", quietly = TRUE)) {
          fit <- vsn::vsn2(expr_matrix)
          vsn::predict(fit, newdata = expr_matrix)
        } else {
          warning("vsn包未安装，使用quantile代替 / vsn package not installed, using quantile instead")
          limma::normalizeBetweenArrays(expr_matrix, method = "quantile")
        }
      }
    )
    
    exprs(eset) <- norm_matrix
    return(eset)
    
  }, error = function(e) {
    stop(paste0("标准化失败 / Normalization failed: ", conditionMessage(e)))
  })
}


#' ComBat批次效应校正
#' ComBat Batch Effect Correction
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param batch 批次变量 / Batch variable
#' @param mod 协变量模型矩阵（可选）/ Model matrix for covariates (optional)
#' @param par.prior 使用参数先验 / Use parametric priors
#'
#' @return ExpressionSet，批次校正后 / Batch-corrected ExpressionSet
#' @export
combat_batch_correction <- function(eset, batch, mod = NULL, par.prior = TRUE) {
  
  if (!requireNamespace("sva", quietly = TRUE)) {
    warning("sva包未安装，跳过批次校正 / sva package not installed, skipping batch correction")
    return(eset)
  }
  
  if (length(unique(batch)) < 2) {
    warning("批次数量少于2，无需校正 / Less than 2 batches, no correction needed")
    return(eset)
  }
  
  expr_matrix <- exprs(eset)
  
  tryCatch({
    # 运行ComBat / Run ComBat
    corrected <- sva::ComBat(
      dat = expr_matrix,
      batch = batch,
      mod = mod,
      par.prior = par.prior
    )
    
    exprs(eset) <- corrected
    message("ComBat批次校正完成 / ComBat batch correction complete")
    
    return(eset)
    
  }, error = function(e) {
    warning(paste0("ComBat校正失败，返回原始数据 / ComBat correction failed, returning original data: ", 
                   conditionMessage(e)))
    return(eset)
  })
}


#' 生成标准化前后对比图
#' Generate Normalization Comparison Plots
#'
#' @param eset_before 标准化前的ExpressionSet / ExpressionSet before normalization
#' @param eset_after 标准化后的ExpressionSet / ExpressionSet after normalization
#' @param method 使用的标准化方法 / Normalization method used
#' @param output_dir 输出目录 / Output directory
#'
#' @return 列表，图文件路径 / List of plot file paths
generate_normalization_plots <- function(eset_before, eset_after, method, output_dir) {
  
  plots <- list()
  
  expr_before <- exprs(eset_before)
  expr_after <- exprs(eset_after)
  
  # 1. 箱线图对比 / Boxplot comparison
  tryCatch({
    boxplot_file <- file.path(output_dir, paste0("normalization_boxplot_", method, ".png"))
    png(boxplot_file, width = 1400, height = 600, res = 100)
    
    par(mfrow = c(1, 2), mar = c(10, 4, 4, 2))
    
    # 限制样本数量用于可视化 / Limit samples for visualization
    n_samples <- min(30, ncol(expr_before))
    
    boxplot(expr_before[, 1:n_samples], las = 2, col = "lightblue",
            main = "Before Normalization / 标准化前",
            ylab = "Expression / 表达值")
    
    boxplot(expr_after[, 1:n_samples], las = 2, col = "lightgreen",
            main = paste0("After ", toupper(method), " / ", toupper(method), "标准化后"),
            ylab = "Expression / 表达值")
    
    dev.off()
    plots$boxplot <- boxplot_file
  }, error = function(e) {
    message(paste0("箱线图生成失败 / Boxplot generation failed: ", conditionMessage(e)))
  })
  
  # 2. 密度图对比 / Density plot comparison
  tryCatch({
    density_file <- file.path(output_dir, paste0("normalization_density_", method, ".png"))
    png(density_file, width = 1400, height = 600, res = 100)
    
    par(mfrow = c(1, 2))
    
    # 标准化前 / Before normalization
    plot(density(expr_before[, 1], na.rm = TRUE), 
         main = "Before Normalization / 标准化前",
         xlab = "Expression / 表达值", col = 1, ylim = c(0, 1))
    for (i in 2:min(20, ncol(expr_before))) {
      lines(density(expr_before[, i], na.rm = TRUE), col = i)
    }
    
    # 标准化后 / After normalization
    plot(density(expr_after[, 1], na.rm = TRUE),
         main = paste0("After ", toupper(method), " / ", toupper(method), "标准化后"),
         xlab = "Expression / 表达值", col = 1, ylim = c(0, 1))
    for (i in 2:min(20, ncol(expr_after))) {
      lines(density(expr_after[, i], na.rm = TRUE), col = i)
    }
    
    dev.off()
    plots$density <- density_file
  }, error = function(e) {
    message(paste0("密度图生成失败 / Density plot generation failed: ", conditionMessage(e)))
  })
  
  # 3. MA图 / MA plot
  tryCatch({
    ma_file <- file.path(output_dir, paste0("normalization_ma_", method, ".png"))
    png(ma_file, width = 1400, height = 600, res = 100)
    
    par(mfrow = c(1, 2))
    
    # 计算中位数参考 / Calculate median reference
    ref_before <- apply(expr_before, 1, median, na.rm = TRUE)
    ref_after <- apply(expr_after, 1, median, na.rm = TRUE)
    
    # 标准化前MA图 / MA plot before
    if (ncol(expr_before) > 1) {
      A_before <- (expr_before[, 1] + ref_before) / 2
      M_before <- expr_before[, 1] - ref_before
      plot(A_before, M_before, pch = ".", col = "grey50",
           main = "Before Normalization / 标准化前",
           xlab = "A (Average)", ylab = "M (Difference)")
      abline(h = 0, col = "red", lwd = 2)
    }
    
    # 标准化后MA图 / MA plot after
    if (ncol(expr_after) > 1) {
      A_after <- (expr_after[, 1] + ref_after) / 2
      M_after <- expr_after[, 1] - ref_after
      plot(A_after, M_after, pch = ".", col = "grey50",
           main = paste0("After ", toupper(method), " / ", toupper(method), "标准化后"),
           xlab = "A (Average)", ylab = "M (Difference)")
      abline(h = 0, col = "red", lwd = 2)
    }
    
    dev.off()
    plots$ma <- ma_file
  }, error = function(e) {
    message(paste0("MA图生成失败 / MA plot generation failed: ", conditionMessage(e)))
  })
  
  message(paste0("标准化对比图已保存到 / Normalization plots saved to: ", output_dir))
  
  return(plots)
}


#' 对数转换表达数据
#' Log Transform Expression Data
#'
#' @param eset ExpressionSet对象 / ExpressionSet object
#' @param base 对数底数 / Log base
#' @param offset 偏移量（避免log(0)）/ Offset to avoid log(0)
#'
#' @return ExpressionSet，对数转换后 / Log-transformed ExpressionSet
#' @export
log_transform <- function(eset, base = 2, offset = 1) {
  
  expr_matrix <- exprs(eset)
  
  # 检查是否已经对数转换 / Check if already log-transformed
  if (max(expr_matrix, na.rm = TRUE) < 30) {
    message("数据可能已经对数转换 / Data appears to be already log-transformed")
    return(eset)
  }
  
  # 对数转换 / Log transform
  expr_log <- log(expr_matrix + offset, base = base)
  
  exprs(eset) <- expr_log
  
  message(paste0("对数转换完成（base=", base, "）/ Log transformation complete (base=", base, ")"))
  
  return(eset)
}


#' 检测数据是否已对数转换
#' Check if Data is Log-Transformed
#'
#' @param eset ExpressionSet对象或表达矩阵 / ExpressionSet or expression matrix
#' @return 逻辑值，TRUE表示已对数转换 / Logical, TRUE if log-transformed
#' @export
is_log_transformed <- function(eset) {
  
  if (inherits(eset, "ExpressionSet")) {
    expr_matrix <- exprs(eset)
  } else {
    expr_matrix <- eset
  }
  
  # 基于最大值判断 / Judge based on max value
  max_val <- max(expr_matrix, na.rm = TRUE)
  
  if (max_val < 30) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


#' 获取可用的标准化方法
#' Get Available Normalization Methods
#'
#' @return 字符向量，可用的方法 / Character vector of available methods
#' @export
get_normalization_methods <- function() {
  
  methods <- c("rma", "quantile", "mas5")
  
  if (requireNamespace("gcrma", quietly = TRUE)) {
    methods <- c(methods, "gcrma")
  }
  
  if (requireNamespace("vsn", quietly = TRUE)) {
    methods <- c(methods, "vsn")
  }
  
  return(methods)
}
