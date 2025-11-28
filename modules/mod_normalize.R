#' ============================================================================
#' Shiny标准化模块
#' Shiny Normalization Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 方法选择 / Method selection
#' - 参数调整 / Parameter adjustment
#' - 前后对比可视化 / Before/after comparison visualization
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

#' 标准化模块UI
#' Normalization Module UI
#'
#' @param id 模块ID / Module ID
#' @return tagList，UI元素 / tagList, UI elements
#' @export
mod_normalize_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(12,
        box(
          title = "标准化设置 / Normalization Settings",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          fluidRow(
            column(4,
              selectInput(
                ns("norm_method"),
                label = "标准化方法 / Method",
                choices = NULL,  # 动态更新 / Dynamic update
                selected = NULL
              )
            ),
            column(4,
              checkboxInput(
                ns("background_correct"),
                "背景校正 / Background Correction",
                value = TRUE
              )
            ),
            column(4,
              checkboxInput(
                ns("batch_correct"),
                "批次校正（ComBat）/ Batch Correction",
                value = FALSE
              )
            )
          ),
          
          conditionalPanel(
            condition = paste0("input['", ns("batch_correct"), "'] == true"),
            selectInput(
              ns("batch_var"),
              label = "批次变量 / Batch Variable",
              choices = NULL
            )
          ),
          
          hr(),
          
          actionButton(
            ns("run_norm"),
            "运行标准化 / Run Normalization",
            icon = icon("magic"),
            class = "btn-primary"
          )
        )
      )
    ),
    
    # 标准化前后对比 / Before/After Comparison
    fluidRow(
      column(6,
        box(
          title = "标准化前 / Before Normalization",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          plotlyOutput(ns("before_boxplot"), height = "350px"),
          hr(),
          plotlyOutput(ns("before_density"), height = "300px")
        )
      ),
      column(6,
        box(
          title = "标准化后 / After Normalization",
          status = "success",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          plotlyOutput(ns("after_boxplot"), height = "350px"),
          hr(),
          plotlyOutput(ns("after_density"), height = "300px")
        )
      )
    ),
    
    # 标准化统计 / Normalization Statistics
    fluidRow(
      column(12,
        box(
          title = "标准化统计 / Normalization Statistics",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          verbatimTextOutput(ns("norm_stats"))
        )
      )
    )
  )
}


#' 标准化模块Server
#' Normalization Module Server
#'
#' @param id 模块ID / Module ID
#' @param r 响应式值对象 / Reactive values object
#' @return 响应式值 / Reactive values
#' @export
mod_normalize_server <- function(id, r = reactiveValues()) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # 更新方法选择 / Update method choices
    observe({
      if (!is.null(r$data_type)) {
        if (r$data_type == "microarray") {
          methods <- c(
            "RMA" = "rma",
            "GCRMA" = "gcrma",
            "MAS5" = "mas5",
            "Quantile" = "quantile",
            "VSN" = "vsn"
          )
        } else {
          methods <- c(
            "VST (DESeq2)" = "vst",
            "rlog (DESeq2)" = "rlog",
            "TMM (edgeR)" = "tmm",
            "RLE (edgeR)" = "rle",
            "voom (limma)" = "voom",
            "CPM" = "cpm"
          )
        }
        
        updateSelectInput(session, "norm_method", choices = methods)
      }
    })
    
    # 更新批次变量选择 / Update batch variable choices
    observe({
      if (!is.null(r$eset)) {
        vars <- varLabels(r$eset)
        updateSelectInput(session, "batch_var", choices = vars)
      }
    })
    
    # 运行标准化 / Run normalization
    observeEvent(input$run_norm, {
      req(input$norm_method)
      
      withProgress(message = "运行标准化 / Running normalization...", {
        
        tryCatch({
          incProgress(0.2, detail = "准备数据 / Preparing data...")
          
          if (r$data_type == "microarray") {
            req(r$eset)
            source("R/microarray/normalize.R", local = TRUE)
            
            batch <- NULL
            if (input$batch_correct && !is.null(input$batch_var)) {
              batch <- pData(r$eset)[[input$batch_var]]
            }
            
            incProgress(0.5, detail = paste0("运行 ", toupper(input$norm_method), "..."))
            
            r$eset_normalized <- normalize_microarray(
              eset = r$eset,
              method = input$norm_method,
              background_correct = input$background_correct,
              batch = batch,
              output_dir = file.path("data", r$gse_id, "normalized")
            )
            
          } else {
            req(r$counts)
            source("R/rnaseq/normalize.R", local = TRUE)
            
            incProgress(0.5, detail = paste0("运行 ", toupper(input$norm_method), "..."))
            
            norm_result <- normalize_rnaseq(
              counts = r$counts,
              method = input$norm_method
            )
            
            r$counts_normalized <- norm_result$normalized
            r$norm_result <- norm_result
          }
          
          incProgress(0.9, detail = "完成 / Complete...")
          
          showNotification("标准化完成！/ Normalization complete!", type = "message")
          
        }, error = function(e) {
          showNotification(paste0("标准化失败 / Normalization failed: ", conditionMessage(e)),
                          type = "error")
        })
      })
    })
    
    # 标准化前箱线图 / Before boxplot
    output$before_boxplot <- renderPlotly({
      if (!is.null(r$eset)) {
        expr <- exprs(r$eset)
      } else if (!is.null(r$counts)) {
        expr <- log2(r$counts + 1)
      } else {
        return(NULL)
      }
      
      n_samples <- min(20, ncol(expr))
      expr_subset <- expr[, 1:n_samples, drop = FALSE]
      
      plot_ly(type = "box") %>%
        add_trace(y = ~as.vector(expr_subset),
                  x = ~rep(colnames(expr_subset), each = nrow(expr_subset)),
                  marker = list(color = "steelblue")) %>%
        layout(title = "Before / 标准化前",
               xaxis = list(title = "", tickangle = 45),
               yaxis = list(title = "Expression / 表达值"),
               showlegend = FALSE)
    })
    
    # 标准化后箱线图 / After boxplot
    output$after_boxplot <- renderPlotly({
      if (!is.null(r$eset_normalized)) {
        expr <- exprs(r$eset_normalized)
      } else if (!is.null(r$counts_normalized)) {
        expr <- r$counts_normalized
      } else {
        return(plotly_empty() %>% layout(title = "请先运行标准化 / Run normalization first"))
      }
      
      n_samples <- min(20, ncol(expr))
      expr_subset <- expr[, 1:n_samples, drop = FALSE]
      
      plot_ly(type = "box") %>%
        add_trace(y = ~as.vector(expr_subset),
                  x = ~rep(colnames(expr_subset), each = nrow(expr_subset)),
                  marker = list(color = "green")) %>%
        layout(title = "After / 标准化后",
               xaxis = list(title = "", tickangle = 45),
               yaxis = list(title = "Expression / 表达值"),
               showlegend = FALSE)
    })
    
    # 标准化前密度图 / Before density
    output$before_density <- renderPlotly({
      if (!is.null(r$eset)) {
        expr <- exprs(r$eset)
      } else if (!is.null(r$counts)) {
        expr <- log2(r$counts + 1)
      } else {
        return(NULL)
      }
      
      n_samples <- min(10, ncol(expr))
      
      p <- plot_ly()
      for (i in 1:n_samples) {
        dens <- density(expr[, i], na.rm = TRUE)
        p <- p %>% add_trace(x = dens$x, y = dens$y, type = "scatter", mode = "lines",
                             name = colnames(expr)[i])
      }
      
      p %>% layout(title = "Density / 密度分布",
                   xaxis = list(title = "Expression / 表达值"),
                   yaxis = list(title = "Density / 密度"),
                   showlegend = FALSE)
    })
    
    # 标准化后密度图 / After density
    output$after_density <- renderPlotly({
      if (!is.null(r$eset_normalized)) {
        expr <- exprs(r$eset_normalized)
      } else if (!is.null(r$counts_normalized)) {
        expr <- r$counts_normalized
      } else {
        return(plotly_empty())
      }
      
      n_samples <- min(10, ncol(expr))
      
      p <- plot_ly()
      for (i in 1:n_samples) {
        dens <- density(expr[, i], na.rm = TRUE)
        p <- p %>% add_trace(x = dens$x, y = dens$y, type = "scatter", mode = "lines",
                             name = colnames(expr)[i])
      }
      
      p %>% layout(title = "Density / 密度分布",
                   xaxis = list(title = "Expression / 表达值"),
                   yaxis = list(title = "Density / 密度"),
                   showlegend = FALSE)
    })
    
    # 标准化统计 / Normalization statistics
    output$norm_stats <- renderText({
      if (!is.null(r$eset_normalized)) {
        expr <- exprs(r$eset_normalized)
        paste0(
          "标准化方法 / Method: ", input$norm_method, "\n",
          "样本数 / Samples: ", ncol(expr), "\n",
          "特征数 / Features: ", nrow(expr), "\n",
          "表达值范围 / Expression range: ", round(min(expr, na.rm = TRUE), 2), 
          " - ", round(max(expr, na.rm = TRUE), 2), "\n"
        )
      } else if (!is.null(r$counts_normalized)) {
        expr <- r$counts_normalized
        paste0(
          "标准化方法 / Method: ", input$norm_method, "\n",
          "样本数 / Samples: ", ncol(expr), "\n",
          "基因数 / Genes: ", nrow(expr), "\n",
          "表达值范围 / Expression range: ", round(min(expr, na.rm = TRUE), 2), 
          " - ", round(max(expr, na.rm = TRUE), 2), "\n"
        )
      } else {
        "请先运行标准化 / Run normalization first"
      }
    })
    
    return(r)
  })
}
