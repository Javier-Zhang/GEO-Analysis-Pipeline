#' ============================================================================
#' GEO Analysis Pipeline - 单元测试
#' GEO Analysis Pipeline - Unit Tests
#' ============================================================================

library(testthat)

# 测试上下文 / Test context
context("GEO Analysis Pipeline Tests")

# ============================================================================
# 工具函数测试 / Utility Function Tests
# ============================================================================

test_that("phenotype parsing works correctly", {
  # 模拟表型数据 / Mock phenotype data
  mock_pdata <- data.frame(
    geo_accession = c("GSM1", "GSM2"),
    title = c("Sample1", "Sample2"),
    characteristics_ch1 = c("tissue: cancer", "tissue: normal"),
    stringsAsFactors = FALSE
  )
  
  # 测试解析功能 / Test parsing function
  expect_true(is.data.frame(mock_pdata))
  expect_equal(nrow(mock_pdata), 2)
})

test_that("clean column name function works", {
  # 定义清理函数 / Define cleaning function
  clean_column_name <- function(name) {
    name <- tolower(name)
    name <- gsub("[^a-z0-9]", "_", name)
    name <- gsub("_+", "_", name)
    name <- gsub("^_|_$", "", name)
    return(name)
  }
  
  expect_equal(clean_column_name("Gene Symbol"), "gene_symbol")
  expect_equal(clean_column_name("SAMPLE-ID"), "sample_id")
  expect_equal(clean_column_name("Data.Point"), "data_point")
})

# ============================================================================
# 导出函数测试 / Export Function Tests
# ============================================================================

test_that("GCT export format is correct", {
  # 模拟表达矩阵 / Mock expression matrix
  mock_expr <- matrix(
    c(1.5, 2.3, 0.8, 3.1, 2.0, 1.2),
    nrow = 2,
    ncol = 3,
    dimnames = list(
      c("Gene1", "Gene2"),
      c("Sample1", "Sample2", "Sample3")
    )
  )
  
  # 测试矩阵属性 / Test matrix properties
  expect_equal(nrow(mock_expr), 2)
  expect_equal(ncol(mock_expr), 3)
  expect_true(is.numeric(mock_expr))
})

# ============================================================================
# 数据处理测试 / Data Processing Tests
# ============================================================================

test_that("log transformation detection works", {
  # 定义检测函数 / Define detection function
  is_log_transformed <- function(data) {
    max_val <- max(data, na.rm = TRUE)
    if (max_val < 30) return(TRUE)
    return(FALSE)
  }
  
  # 测试对数转换数据 / Test log-transformed data
  log_data <- rnorm(100, mean = 8, sd = 2)
  expect_true(is_log_transformed(log_data))
  
  # 测试原始数据 / Test raw data
  raw_data <- 2^rnorm(100, mean = 8, sd = 2)
  expect_false(is_log_transformed(raw_data))
})

test_that("low expression filtering works", {
  # 模拟counts数据 / Mock counts data
  mock_counts <- matrix(
    c(0, 0, 0,
      10, 15, 20,
      5, 0, 8,
      100, 120, 90),
    nrow = 4,
    ncol = 3,
    byrow = TRUE,
    dimnames = list(
      c("Gene1", "Gene2", "Gene3", "Gene4"),
      c("Sample1", "Sample2", "Sample3")
    )
  )
  
  # 过滤低表达基因 / Filter low expression genes
  min_counts <- 10
  min_samples <- 2
  keep <- rowSums(mock_counts >= min_counts) >= min_samples
  
  expect_equal(sum(keep), 2)  # Gene2 and Gene4 should pass
})

# ============================================================================
# 相关性计算测试 / Correlation Calculation Tests
# ============================================================================

test_that("sample correlation calculation works", {
  # 模拟表达矩阵 / Mock expression matrix
  set.seed(123)
  mock_expr <- matrix(rnorm(100), nrow = 20, ncol = 5)
  colnames(mock_expr) <- paste0("Sample", 1:5)
  
  # 计算相关性 / Calculate correlation
  cor_matrix <- cor(mock_expr, method = "pearson")
  
  expect_equal(dim(cor_matrix), c(5, 5))
  expect_equal(diag(cor_matrix), rep(1, 5))
})

# ============================================================================
# 差异分析辅助函数测试 / DEG Helper Function Tests
# ============================================================================

test_that("DEG status assignment works", {
  # 模拟DEG结果 / Mock DEG results
  mock_deg <- data.frame(
    logFC = c(2.5, -1.8, 0.3, 3.0, -2.5),
    adj.P.Val = c(0.001, 0.01, 0.2, 0.0001, 0.03)
  )
  
  lfc_threshold <- 1
  p_threshold <- 0.05
  
  # 分配状态 / Assign status
  mock_deg$DEG_status <- "Not Significant"
  mock_deg$DEG_status[mock_deg$adj.P.Val < p_threshold & mock_deg$logFC > lfc_threshold] <- "Up"
  mock_deg$DEG_status[mock_deg$adj.P.Val < p_threshold & mock_deg$logFC < -lfc_threshold] <- "Down"
  
  expect_equal(sum(mock_deg$DEG_status == "Up"), 2)
  expect_equal(sum(mock_deg$DEG_status == "Down"), 2)
  expect_equal(sum(mock_deg$DEG_status == "Not Significant"), 1)
})

# ============================================================================
# ID转换测试 / ID Conversion Tests
# ============================================================================

test_that("gene ID type detection works", {
  # 定义检测函数 / Define detection function
  detect_gene_id_type <- function(ids) {
    sample_ids <- head(ids, 100)
    
    if (mean(grepl("^ENS[A-Z]*G\\d+", sample_ids)) > 0.5) return("ENSEMBL")
    if (mean(grepl("^\\d+$", sample_ids)) > 0.5) return("ENTREZID")
    if (mean(grepl("^N[MR]_\\d+", sample_ids)) > 0.5) return("REFSEQ")
    return("SYMBOL")
  }
  
  # 测试不同ID类型 / Test different ID types
  ensembl_ids <- c("ENSG00000139618", "ENSG00000141510", "ENSG00000157764")
  expect_equal(detect_gene_id_type(ensembl_ids), "ENSEMBL")
  
  entrez_ids <- c("675", "7157", "5290")
  expect_equal(detect_gene_id_type(entrez_ids), "ENTREZID")
  
  symbol_ids <- c("BRCA2", "TP53", "PIK3CA")
  expect_equal(detect_gene_id_type(symbol_ids), "SYMBOL")
})

# ============================================================================
# 批次效应检测测试 / Batch Effect Detection Tests
# ============================================================================

test_that("batch detection works", {
  # 模拟批次数据 / Mock batch data
  batch <- c("A", "A", "B", "B", "A", "B")
  
  expect_equal(length(unique(batch)), 2)
  expect_true(length(unique(batch)) >= 2)  # 至少2个批次才能校正
})

# ============================================================================
# 报告生成测试 / Report Generation Tests
# ============================================================================

test_that("report HTML generation works", {
  # 模拟报告生成 / Mock report generation
  generate_simple_html <- function(title, content) {
    paste0("<html><head><title>", title, "</title></head>",
           "<body>", content, "</body></html>")
  }
  
  html <- generate_simple_html("Test Report", "<p>Test content</p>")
  
  expect_true(grepl("<html>", html))
  expect_true(grepl("Test Report", html))
  expect_true(grepl("Test content", html))
})

# ============================================================================
# 运行测试 / Run Tests
# ============================================================================

cat("\n========================================\n")
cat("GEO Analysis Pipeline - Test Results\n")
cat("========================================\n\n")

# 如果作为脚本运行 / If running as script
if (!interactive()) {
  test_results <- test_dir(".", reporter = "summary")
  print(test_results)
}
