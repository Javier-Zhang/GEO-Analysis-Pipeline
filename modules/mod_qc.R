#' ============================================================================
#' Shiny质控模块
#' Shiny QC Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 交互式质控图表（plotly）/ Interactive QC charts (plotly)
#' - 异常样本标记 / Outlier sample marking
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

#' 质控模块UI
#' QC Module UI
#'
#' @param id 模块ID / Module ID
#' @return tagList，UI元素 / tagList, UI elements
#' @export
mod_qc_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(12,
        box(
          title = "质控设置 / QC Settings",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          fluidRow(
            column(4,
              sliderInput(
                ns("outlier_threshold"),
                label = "异常值阈值 / Outlier Threshold",
                min = 0.5,
                max = 0.95,
                value = 0.8,
                step = 0.05
              )
            ),
            column(4,
              selectInput(
                ns("correlation_method"),
                label = "相关性方法 / Correlation Method",
                choices = c("pearson", "spearman"),
                selected = "pearson"
              )
            ),
            column(4,
              br(),
              actionButton(
                ns("run_qc"),
                "运行质控 / Run QC",
                icon = icon("chart-line"),
                class = "btn-primary",
                width = "100%"
              )
            )
          )
        )
      )
    ),
    
    # QC摘要 / QC Summary
    fluidRow(
      column(3,
        valueBoxOutput(ns("n_samples"), width = 12)
      ),
      column(3,
        valueBoxOutput(ns("n_outliers"), width = 12)
      ),
      column(3,
        valueBoxOutput(ns("mean_correlation"), width = 12)
      ),
      column(3,
        valueBoxOutput(ns("pc1_variance"), width = 12)
      )
    ),
    
    # QC图表 / QC Plots
    fluidRow(
      column(6,
        box(
          title = "表达分布 / Expression Distribution",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          plotlyOutput(ns("boxplot"), height = "400px")
        )
      ),
      column(6,
        box(
          title = "样本相关性热图 / Correlation Heatmap",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          plotOutput(ns("correlation_heatmap"), height = "400px")
        )
      )
    ),
    
    fluidRow(
      column(6,
        box(
          title = "PCA图 / PCA Plot",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          plotlyOutput(ns("pca_plot"), height = "400px")
        )
      ),
      column(6,
        box(
          title = "平均相关性 / Mean Correlation",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          plotlyOutput(ns("mean_cor_plot"), height = "400px")
        )
      )
    ),
    
    # 异常样本表格 / Outlier Samples Table
    fluidRow(
      column(12,
        box(
          title = "异常样本 / Outlier Samples",
          status = "warning",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          DT::dataTableOutput(ns("outlier_table")),
          
          hr(),
          
          actionButton(
            ns("remove_outliers"),
            "移除异常样本 / Remove Outliers",
            icon = icon("trash"),
            class = "btn-warning"
          )
        )
      )
    )
  )
}


#' 质控模块Server
#' QC Module Server
#'
#' @param id 模块ID / Module ID
#' @param r 响应式值对象 / Reactive values object
#' @return 响应式值 / Reactive values
#' @export
mod_qc_server <- function(id, r = reactiveValues()) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # 运行QC / Run QC
    observeEvent(input$run_qc, {
      req(r$download_result)
      
      withProgress(message = "运行质控分析 / Running QC analysis...", {
        
        tryCatch({
          incProgress(0.2, detail = "准备数据 / Preparing data...")
          
          # 获取表达数据 / Get expression data
          if (r$data_type == "microarray") {
            gse <- r$download_result$gse_object
            if (is.list(gse)) {
              eset <- gse[[1]]
            } else {
              eset <- gse
            }
            
            # 运行微阵列QC / Run microarray QC
            source("R/microarray/qc.R", local = TRUE)
            
            incProgress(0.5, detail = "计算质控指标 / Calculating QC metrics...")
            
            r$qc_results <- run_microarray_qc(
              eset = eset,
              output_dir = file.path("data", r$gse_id, "qc")
            )
            r$eset <- eset
            
          } else {
            # RNA-seq QC
            counts <- r$download_result$counts_matrix
            
            if (!is.null(counts)) {
              source("R/rnaseq/qc.R", local = TRUE)
              
              incProgress(0.5, detail = "计算质控指标 / Calculating QC metrics...")
              
              r$qc_results <- run_rnaseq_qc(
                counts = counts,
                output_dir = file.path("data", r$gse_id, "qc")
              )
              r$counts <- counts
            } else {
              showNotification("无counts数据 / No counts data available",
                              type = "error")
              return()
            }
          }
          
          incProgress(0.9, detail = "完成 / Complete...")
          
          showNotification("质控完成！/ QC complete!", type = "message")
          
        }, error = function(e) {
          showNotification(paste0("质控失败 / QC failed: ", conditionMessage(e)),
                          type = "error")
        })
      })
    })
    
    # 值框：样本数 / Value box: Samples
    output$n_samples <- renderValueBox({
      n <- if (!is.null(r$qc_results)) r$qc_results$summary$n_samples else "-"
      valueBox(n, "样本数 / Samples", icon = icon("users"), color = "blue")
    })
    
    # 值框：异常样本 / Value box: Outliers
    output$n_outliers <- renderValueBox({
      n <- if (!is.null(r$qc_results)) r$qc_results$summary$n_outliers else "-"
      color <- if (is.numeric(n) && n > 0) "red" else "green"
      valueBox(n, "异常样本 / Outliers", icon = icon("exclamation-triangle"), color = color)
    })
    
    # 值框：平均相关性 / Value box: Mean correlation
    output$mean_correlation <- renderValueBox({
      val <- if (!is.null(r$qc_results)) round(r$qc_results$summary$mean_correlation, 3) else "-"
      valueBox(val, "平均相关性 / Mean Cor.", icon = icon("link"), color = "purple")
    })
    
    # 值框：PC1方差 / Value box: PC1 variance
    output$pc1_variance <- renderValueBox({
      val <- if (!is.null(r$qc_results)) paste0(round(r$qc_results$summary$pca_var_pc1, 1), "%") else "-"
      valueBox(val, "PC1方差 / PC1 Var.", icon = icon("chart-pie"), color = "orange")
    })
    
    # 箱线图 / Boxplot
    output$boxplot <- renderPlotly({
      req(r$qc_results)
      
      # 简化数据用于展示 / Simplify data for display
      if (!is.null(r$eset)) {
        expr <- exprs(r$eset)
      } else if (!is.null(r$counts)) {
        expr <- log2(r$counts + 1)
      } else {
        return(NULL)
      }
      
      # 限制样本数 / Limit samples
      n_samples <- min(30, ncol(expr))
      expr_subset <- expr[, 1:n_samples, drop = FALSE]
      
      plot_ly(type = "box") %>%
        add_trace(y = ~as.vector(expr_subset), 
                  x = ~rep(colnames(expr_subset), each = nrow(expr_subset)),
                  color = ~rep(colnames(expr_subset), each = nrow(expr_subset))) %>%
        layout(title = "Expression Distribution / 表达分布",
               xaxis = list(title = "Sample / 样本"),
               yaxis = list(title = "Expression / 表达值"),
               showlegend = FALSE)
    })
    
    # 相关性热图 / Correlation heatmap
    output$correlation_heatmap <- renderPlot({
      req(r$qc_results)
      
      pheatmap::pheatmap(
        r$qc_results$correlation_matrix,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        main = "Sample Correlation / 样本相关性",
        show_rownames = ncol(r$qc_results$correlation_matrix) <= 30,
        show_colnames = ncol(r$qc_results$correlation_matrix) <= 30
      )
    })
    
    # PCA图 / PCA plot
    output$pca_plot <- renderPlotly({
      req(r$qc_results)
      
      pca_df <- data.frame(
        PC1 = r$qc_results$pca_results$scores[, 1],
        PC2 = r$qc_results$pca_results$scores[, 2],
        Sample = rownames(r$qc_results$pca_results$scores)
      )
      
      var_exp <- r$qc_results$pca_results$var_explained
      
      # 标记异常样本 / Mark outlier samples
      pca_df$IsOutlier <- pca_df$Sample %in% r$qc_results$outliers$all_outliers
      
      plot_ly(pca_df, x = ~PC1, y = ~PC2, text = ~Sample, 
              color = ~IsOutlier, 
              colors = c("steelblue", "red"),
              type = "scatter", mode = "markers") %>%
        layout(
          title = "PCA Plot / PCA图",
          xaxis = list(title = paste0("PC1 (", round(var_exp[1], 1), "%)")),
          yaxis = list(title = paste0("PC2 (", round(var_exp[2], 1), "%)"))
        )
    })
    
    # 平均相关性图 / Mean correlation plot
    output$mean_cor_plot <- renderPlotly({
      req(r$qc_results)
      
      mean_cors <- r$qc_results$outliers$mean_correlations
      
      cor_df <- data.frame(
        Sample = names(mean_cors),
        Correlation = mean_cors
      )
      cor_df$IsOutlier <- cor_df$Sample %in% r$qc_results$outliers$correlation_outliers
      cor_df <- cor_df[order(cor_df$Correlation), ]
      cor_df$Sample <- factor(cor_df$Sample, levels = cor_df$Sample)
      
      plot_ly(cor_df, x = ~Sample, y = ~Correlation, 
              color = ~IsOutlier,
              colors = c("steelblue", "red"),
              type = "bar") %>%
        add_trace(y = rep(input$outlier_threshold, nrow(cor_df)),
                  type = "scatter", mode = "lines",
                  line = list(color = "red", dash = "dash"),
                  name = "Threshold / 阈值") %>%
        layout(
          title = "Mean Sample Correlation / 平均样本相关性",
          xaxis = list(title = "Sample / 样本", tickangle = 45),
          yaxis = list(title = "Mean Correlation / 平均相关性"),
          showlegend = TRUE
        )
    })
    
    # 异常样本表格 / Outlier table
    output$outlier_table <- DT::renderDataTable({
      req(r$qc_results)
      
      outlier_df <- r$qc_results$outliers$outlier_summary
      
      DT::datatable(
        outlier_df,
        options = list(
          pageLength = 10,
          scrollX = TRUE
        ),
        rownames = FALSE
      )
    })
    
    # 移除异常样本 / Remove outliers
    observeEvent(input$remove_outliers, {
      req(r$qc_results)
      
      outliers <- r$qc_results$outliers$all_outliers
      
      if (length(outliers) == 0) {
        showNotification("没有异常样本需要移除 / No outliers to remove",
                        type = "message")
        return()
      }
      
      tryCatch({
        if (!is.null(r$eset)) {
          source("R/microarray/qc.R", local = TRUE)
          r$eset <- remove_outlier_samples(r$eset, outliers)
          showNotification(paste0("已移除 ", length(outliers), " 个异常样本 / Removed ", 
                                  length(outliers), " outliers"),
                          type = "message")
        } else if (!is.null(r$counts)) {
          source("R/rnaseq/qc.R", local = TRUE)
          r$counts <- filter_low_quality_samples(r$counts, outliers)
          showNotification(paste0("已移除 ", length(outliers), " 个异常样本 / Removed ", 
                                  length(outliers), " outliers"),
                          type = "message")
        }
      }, error = function(e) {
        showNotification(paste0("移除失败 / Removal failed: ", conditionMessage(e)),
                        type = "error")
      })
    })
    
    return(r)
  })
}
