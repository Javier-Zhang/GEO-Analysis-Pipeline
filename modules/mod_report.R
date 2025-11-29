#' ============================================================================
#' Shiny报告模块
#' Shiny Report Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 报告预览 / Report preview
#' - 下载按钮 / Download button
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

#' 报告模块UI
#' Report Module UI
#'
#' @param id 模块ID / Module ID
#' @return tagList，UI元素 / tagList, UI elements
#' @export
mod_report_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(12,
        box(
          title = "报告设置 / Report Settings",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          fluidRow(
            column(4,
              textInput(
                ns("report_title"),
                label = "报告标题 / Report Title",
                value = "GEO Data Analysis Report / GEO数据分析报告"
              )
            ),
            column(4,
              textInput(
                ns("author"),
                label = "作者 / Author",
                value = ""
              )
            ),
            column(4,
              selectInput(
                ns("report_format"),
                label = "输出格式 / Output Format",
                choices = c(
                  "HTML" = "html",
                  "PDF" = "pdf"
                ),
                selected = "html"
              )
            )
          ),
          
          hr(),
          
          fluidRow(
            column(12,
              h4("包含内容 / Include Content:"),
              fluidRow(
                column(3,
                  checkboxInput(ns("include_qc"), "质控结果 / QC Results", TRUE)
                ),
                column(3,
                  checkboxInput(ns("include_norm"), "标准化结果 / Normalization", TRUE)
                ),
                column(3,
                  checkboxInput(ns("include_annot"), "注释结果 / Annotation", TRUE)
                ),
                column(3,
                  checkboxInput(ns("include_deg"), "差异分析 / DEG Analysis", TRUE)
                )
              )
            )
          ),
          
          hr(),
          
          fluidRow(
            column(6,
              actionButton(
                ns("generate_report"),
                "生成报告 / Generate Report",
                icon = icon("file-alt"),
                class = "btn-primary",
                width = "100%"
              )
            ),
            column(6,
              downloadButton(
                ns("download_report"),
                "下载报告 / Download Report",
                class = "btn-success",
                style = "width: 100%"
              )
            )
          )
        )
      )
    ),
    
    # 报告预览 / Report Preview
    fluidRow(
      column(12,
        box(
          title = "报告预览 / Report Preview",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          uiOutput(ns("report_preview"))
        )
      )
    ),
    
    # 分析摘要 / Analysis Summary
    fluidRow(
      column(12,
        box(
          title = "分析摘要 / Analysis Summary",
          status = "success",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          verbatimTextOutput(ns("analysis_summary"))
        )
      )
    ),
    
    # 导出选项 / Export Options
    fluidRow(
      column(12,
        box(
          title = "数据导出 / Data Export",
          status = "warning",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          fluidRow(
            column(3,
              downloadButton(ns("export_expression"), "表达矩阵 / Expression")
            ),
            column(3,
              downloadButton(ns("export_phenotype"), "表型数据 / Phenotype")
            ),
            column(3,
              downloadButton(ns("export_deg"), "差异基因 / DEG Table")
            ),
            column(3,
              downloadButton(ns("export_gct"), "GCT (GSEA)")
            )
          )
        )
      )
    )
  )
}


#' 报告模块Server
#' Report Module Server
#'
#' @param id 模块ID / Module ID
#' @param r 响应式值对象 / Reactive values object
#' @return 响应式值 / Reactive values
#' @export
mod_report_server <- function(id, r = reactiveValues()) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # 报告内容 / Report content
    report_html <- reactiveVal(NULL)
    
    # 生成报告 / Generate report
    observeEvent(input$generate_report, {
      
      withProgress(message = "生成报告 / Generating report...", {
        
        tryCatch({
          incProgress(0.3, detail = "收集数据 / Collecting data...")
          
          # 生成报告HTML / Generate report HTML
          html <- generate_report_html(
            title = input$report_title,
            author = input$author,
            gse_id = r$gse_id,
            data_type = r$data_type,
            qc_results = if (input$include_qc) r$qc_results else NULL,
            deg_results = if (input$include_deg) r$deg_result else NULL
          )
          
          report_html(html)
          
          incProgress(0.9, detail = "完成 / Complete...")
          
          showNotification("报告生成完成！/ Report generated!", type = "message")
          
        }, error = function(e) {
          showNotification(paste0("报告生成失败 / Report generation failed: ", conditionMessage(e)),
                          type = "error")
        })
      })
    })
    
    # 报告预览 / Report preview
    output$report_preview <- renderUI({
      html <- report_html()
      
      if (is.null(html)) {
        return(tags$div(
          class = "alert alert-info",
          "请先生成报告 / Please generate report first"
        ))
      }
      
      HTML(html)
    })
    
    # 分析摘要 / Analysis summary
    output$analysis_summary <- renderText({
      summary_text <- paste0(
        "GEO数据分析摘要 / GEO Data Analysis Summary\n",
        "=" %>% rep(50) %>% paste(collapse = ""), "\n\n",
        "数据集 / Dataset: ", ifelse(is.null(r$gse_id), "N/A", r$gse_id), "\n",
        "数据类型 / Data Type: ", ifelse(is.null(r$data_type), "N/A", r$data_type), "\n",
        "分析时间 / Analysis Time: ", as.character(Sys.time()), "\n\n"
      )
      
      # QC摘要 / QC summary
      if (!is.null(r$qc_results)) {
        summary_text <- paste0(
          summary_text,
          "质控结果 / QC Results:\n",
          "-" %>% rep(30) %>% paste(collapse = ""), "\n",
          "  样本数 / Samples: ", r$qc_results$summary$n_samples, "\n",
          "  异常样本 / Outliers: ", r$qc_results$summary$n_outliers, "\n",
          "  平均相关性 / Mean Correlation: ", round(r$qc_results$summary$mean_correlation, 3), "\n\n"
        )
      }
      
      # 差异分析摘要 / DEG summary
      if (!is.null(r$deg_result)) {
        summary_text <- paste0(
          summary_text,
          "差异分析结果 / DEG Results:\n",
          "-" %>% rep(30) %>% paste(collapse = ""), "\n",
          "  总基因数 / Total Genes: ", r$deg_result$summary$n_total, "\n",
          "  上调基因 / Up-regulated: ", r$deg_result$summary$n_up, "\n",
          "  下调基因 / Down-regulated: ", r$deg_result$summary$n_down, "\n",
          "  显著基因 / Significant: ", r$deg_result$summary$n_significant, "\n"
        )
      }
      
      return(summary_text)
    })
    
    # 下载报告 / Download report
    output$download_report <- downloadHandler(
      filename = function() {
        paste0(r$gse_id %||% "analysis", "_report.", input$report_format)
      },
      content = function(file) {
        html <- report_html()
        
        if (is.null(html)) {
          html <- "<html><body><p>No report generated</p></body></html>"
        }
        
        writeLines(html, file)
      }
    )
    
    # 导出表达矩阵 / Export expression
    output$export_expression <- downloadHandler(
      filename = function() paste0(r$gse_id, "_expression.csv"),
      content = function(file) {
        if (!is.null(r$eset_annotated %||% r$eset_normalized %||% r$eset)) {
          eset <- r$eset_annotated %||% r$eset_normalized %||% r$eset
          write.csv(exprs(eset), file)
        } else if (!is.null(r$counts_normalized %||% r$counts)) {
          write.csv(r$counts_normalized %||% r$counts, file)
        }
      }
    )
    
    # 导出表型数据 / Export phenotype
    output$export_phenotype <- downloadHandler(
      filename = function() paste0(r$gse_id, "_phenotype.csv"),
      content = function(file) {
        if (!is.null(r$eset)) {
          write.csv(pData(r$eset), file)
        } else if (!is.null(r$download_result$phenotype_data)) {
          write.csv(r$download_result$phenotype_data, file)
        }
      }
    )
    
    # 导出差异基因 / Export DEG
    output$export_deg <- downloadHandler(
      filename = function() paste0(r$gse_id, "_deg.csv"),
      content = function(file) {
        if (!is.null(r$deg_result)) {
          write.csv(r$deg_result$deg_table, file)
        }
      }
    )
    
    # 导出GCT / Export GCT
    output$export_gct <- downloadHandler(
      filename = function() paste0(r$gse_id, ".gct"),
      content = function(file) {
        if (!is.null(r$eset_annotated %||% r$eset_normalized %||% r$eset)) {
          eset <- r$eset_annotated %||% r$eset_normalized %||% r$eset
          source("R/utils/export.R", local = TRUE)
          export_gct(exprs(eset), file)
        } else if (!is.null(r$counts_normalized %||% r$counts)) {
          source("R/utils/export.R", local = TRUE)
          export_gct(r$counts_normalized %||% r$counts, file)
        }
      }
    )
    
    return(r)
  })
}


#' 生成报告HTML
#' Generate Report HTML
#'
#' @param title 标题 / Title
#' @param author 作者 / Author
#' @param gse_id GSE ID
#' @param data_type 数据类型 / Data type
#' @param qc_results QC结果 / QC results
#' @param deg_results 差异分析结果 / DEG results
#' @return 字符串，HTML内容 / Character, HTML content
generate_report_html <- function(title, author, gse_id, data_type, 
                                  qc_results = NULL, deg_results = NULL) {
  
  html <- paste0('
<div class="report-container">
  <h2>', title, '</h2>
  <p><strong>Dataset / 数据集:</strong> ', gse_id, '</p>
  <p><strong>Data Type / 数据类型:</strong> ', data_type, '</p>
  <p><strong>Author / 作者:</strong> ', author, '</p>
  <p><strong>Date / 日期:</strong> ', as.character(Sys.Date()), '</p>
  <hr>
')
  
  # QC部分 / QC section
  if (!is.null(qc_results)) {
    html <- paste0(html, '
  <h3>质控结果 / QC Results</h3>
  <table class="table table-striped">
    <tr><th>指标 / Metric</th><th>值 / Value</th></tr>
    <tr><td>样本数 / Samples</td><td>', qc_results$summary$n_samples, '</td></tr>
    <tr><td>异常样本 / Outliers</td><td>', qc_results$summary$n_outliers, '</td></tr>
    <tr><td>平均相关性 / Mean Correlation</td><td>', round(qc_results$summary$mean_correlation, 3), '</td></tr>
    <tr><td>PC1方差 / PC1 Variance</td><td>', round(qc_results$summary$pca_var_pc1, 1), '%</td></tr>
  </table>
')
  }
  
  # DEG部分 / DEG section
  if (!is.null(deg_results)) {
    html <- paste0(html, '
  <h3>差异分析结果 / Differential Expression Results</h3>
  <table class="table table-striped">
    <tr><th>指标 / Metric</th><th>值 / Value</th></tr>
    <tr><td>总基因数 / Total Genes</td><td>', deg_results$summary$n_total, '</td></tr>
    <tr><td>上调基因 / Up-regulated</td><td>', deg_results$summary$n_up, '</td></tr>
    <tr><td>下调基因 / Down-regulated</td><td>', deg_results$summary$n_down, '</td></tr>
    <tr><td>显著基因 / Significant</td><td>', deg_results$summary$n_significant, '</td></tr>
    <tr><td>logFC阈值 / logFC Threshold</td><td>', deg_results$summary$lfc_threshold, '</td></tr>
    <tr><td>P值阈值 / P-value Threshold</td><td>', deg_results$summary$p_threshold, '</td></tr>
  </table>
')
  }
  
  html <- paste0(html, '
  <hr>
  <p><em>Generated by GEO-Analysis-Pipeline</em></p>
</div>
')
  
  return(html)
}

# 辅助函数 / Helper function
`%||%` <- function(a, b) if (!is.null(a)) a else b
