#' ============================================================================
#' Shiny下载模块
#' Shiny Download Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - GSE ID输入 / GSE ID input
#' - 平台选择 / Platform selection
#' - 下载进度显示 / Download progress display
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

#' 下载模块UI
#' Download Module UI
#'
#' @param id 模块ID / Module ID
#' @return tagList，UI元素 / tagList, UI elements
#' @export
mod_download_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(12,
        box(
          title = "数据下载 / Data Download",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          fluidRow(
            column(4,
              textInput(
                ns("gse_id"),
                label = "GSE ID",
                placeholder = "例如 / e.g., GSE6791",
                width = "100%"
              )
            ),
            column(4,
              selectInput(
                ns("data_type"),
                label = "数据类型 / Data Type",
                choices = c(
                  "微阵列 / Microarray" = "microarray",
                  "RNA-seq" = "rnaseq"
                ),
                selected = "microarray",
                width = "100%"
              )
            ),
            column(4,
              br(),
              actionButton(
                ns("download_btn"),
                "下载数据 / Download",
                icon = icon("download"),
                class = "btn-primary",
                width = "100%"
              )
            )
          ),
          
          hr(),
          
          fluidRow(
            column(6,
              checkboxInput(
                ns("download_raw"),
                "下载原始数据 / Download raw data",
                value = TRUE
              )
            ),
            column(6,
              checkboxInput(
                ns("download_supp"),
                "下载supplementary文件 / Download supplementary files",
                value = TRUE
              )
            )
          )
        )
      )
    ),
    
    fluidRow(
      column(12,
        box(
          title = "下载进度 / Download Progress",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          
          verbatimTextOutput(ns("download_log")),
          
          hr(),
          
          uiOutput(ns("download_status"))
        )
      )
    ),
    
    fluidRow(
      column(12,
        box(
          title = "数据信息 / Data Information",
          status = "success",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          
          uiOutput(ns("data_info"))
        )
      )
    ),
    
    # 示例数据集面板 / Example datasets panel
    fluidRow(
      column(12,
        box(
          title = "示例数据集 / Example Datasets (宫颈癌 / Cervical Cancer)",
          status = "warning",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          
          tableOutput(ns("example_datasets")),
          
          actionButton(
            ns("use_example"),
            "使用示例数据 / Use Example",
            icon = icon("flask"),
            class = "btn-warning"
          )
        )
      )
    )
  )
}


#' 下载模块Server
#' Download Module Server
#'
#' @param id 模块ID / Module ID
#' @param r 响应式值对象 / Reactive values object
#' @return 响应式值 / Reactive values
#' @export
mod_download_server <- function(id, r = reactiveValues()) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # 下载日志 / Download log
    download_log <- reactiveVal("")
    
    # 示例数据集 / Example datasets
    output$example_datasets <- renderTable({
      data.frame(
        GSE_ID = c("GSE6791", "GSE7803", "GSE9750", "GSE63514", "GSE39001"),
        Platform = c("GPL570", "GPL570", "GPL96", "GPL570", "GPL570"),
        Samples = c(42, 28, 66, 128, 108),
        Description = c(
          "宫颈癌 vs 正常 / Cervical cancer vs Normal",
          "宫颈癌表达谱 / Cervical cancer expression",
          "宫颈癌基因表达 / Cervical cancer genes",
          "宫颈癌进展 / Cervical cancer progression",
          "HPV相关宫颈癌 / HPV-associated cervical cancer"
        ),
        stringsAsFactors = FALSE
      )
    })
    
    # 使用示例数据 / Use example data
    observeEvent(input$use_example, {
      updateTextInput(session, "gse_id", value = "GSE6791")
      showNotification("已选择示例数据集 GSE6791 / Selected example dataset GSE6791",
                       type = "message")
    })
    
    # 下载数据 / Download data
    observeEvent(input$download_btn, {
      req(input$gse_id)
      
      # 验证GSE ID / Validate GSE ID
      if (!grepl("^GSE\\d+$", input$gse_id)) {
        showNotification("无效的GSE ID格式 / Invalid GSE ID format",
                        type = "error")
        return()
      }
      
      # 更新日志 / Update log
      log_text <- paste0(
        "开始下载 / Starting download: ", input$gse_id, "\n",
        "数据类型 / Data type: ", input$data_type, "\n",
        "时间 / Time: ", Sys.time(), "\n",
        "-----------------------------------\n"
      )
      download_log(log_text)
      
      # 显示进度 / Show progress
      withProgress(message = paste0("下载中 / Downloading: ", input$gse_id), {
        
        tryCatch({
          # 下载数据 / Download data
          incProgress(0.2, detail = "连接GEO / Connecting to GEO...")
          
          if (input$data_type == "microarray") {
            # 微阵列数据下载 / Microarray data download
            # 注意：在Shiny应用中，这些函数应通过app.R中的source()预加载
            # Note: In Shiny app, these functions should be pre-loaded via source() in app.R
            if (!exists("download_geo_data")) {
              source("R/microarray/download.R", local = FALSE)
            }
            
            result <- download_geo_data(
              gse_id = input$gse_id,
              dest_dir = "data",
              download_raw = input$download_raw,
              download_supp = input$download_supp
            )
            
          } else {
            # RNA-seq数据下载 / RNA-seq data download
            if (!exists("download_rnaseq_data")) {
              source("R/rnaseq/download.R", local = FALSE)
            }
            
            result <- download_rnaseq_data(
              gse_id = input$gse_id,
              dest_dir = "data",
              download_supp = input$download_supp
            )
          }
          
          incProgress(0.8, detail = "处理数据 / Processing data...")
          
          # 保存结果 / Save result
          r$download_result <- result
          r$gse_id <- input$gse_id
          r$data_type <- input$data_type
          
          # 更新日志 / Update log
          log_text <- paste0(
            download_log(),
            "下载完成！/ Download complete!\n",
            "状态 / Status: ", result$status, "\n"
          )
          download_log(log_text)
          
          showNotification("下载完成！/ Download complete!", type = "message")
          
        }, error = function(e) {
          log_text <- paste0(
            download_log(),
            "错误 / Error: ", conditionMessage(e), "\n"
          )
          download_log(log_text)
          
          showNotification(paste0("下载失败 / Download failed: ", conditionMessage(e)),
                          type = "error")
        })
      })
    })
    
    # 输出下载日志 / Output download log
    output$download_log <- renderText({
      download_log()
    })
    
    # 下载状态 / Download status
    output$download_status <- renderUI({
      if (!is.null(r$download_result)) {
        if (r$download_result$status == "success") {
          tags$div(
            class = "alert alert-success",
            icon("check-circle"),
            "下载成功！/ Download successful!"
          )
        } else {
          tags$div(
            class = "alert alert-danger",
            icon("exclamation-circle"),
            "下载出错 / Download error"
          )
        }
      }
    })
    
    # 数据信息 / Data information
    output$data_info <- renderUI({
      req(r$download_result)
      
      result <- r$download_result
      
      if (!is.null(result$platform_info)) {
        info <- result$platform_info
        
        tagList(
          tags$p(tags$strong("平台 / Platform: "), info$platform_id),
          tags$p(tags$strong("平台类型 / Platform type: "), info$platform_type),
          tags$p(tags$strong("样本数 / Samples: "), info$n_samples),
          tags$p(tags$strong("特征数 / Features: "), info$n_features)
        )
      } else {
        tags$p("暂无数据信息 / No data information available")
      }
    })
    
    return(r)
  })
}
