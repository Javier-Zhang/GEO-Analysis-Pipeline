#' ============================================================================
#' Shiny差异分析模块
#' Shiny Differential Analysis Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 分组设置 / Group settings
#' - 阈值调整 / Threshold adjustment
#' - 交互式火山图 / Interactive volcano plot
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

#' 差异分析模块UI
#' Differential Analysis Module UI
#'
#' @param id 模块ID / Module ID
#' @return tagList，UI元素 / tagList, UI elements
#' @export
mod_differential_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(12,
        box(
          title = "差异分析设置 / Differential Analysis Settings",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          fluidRow(
            column(3,
              selectInput(
                ns("group_var"),
                label = "分组变量 / Group Variable",
                choices = NULL
              )
            ),
            column(3,
              selectInput(
                ns("control_group"),
                label = "对照组 / Control Group",
                choices = NULL
              )
            ),
            column(3,
              selectInput(
                ns("treatment_group"),
                label = "实验组 / Treatment Group",
                choices = NULL
              )
            ),
            column(3,
              selectInput(
                ns("deg_method"),
                label = "分析方法 / Method",
                choices = NULL
              )
            )
          ),
          
          hr(),
          
          fluidRow(
            column(4,
              sliderInput(
                ns("lfc_threshold"),
                label = "log2FC阈值 / log2FC Threshold",
                min = 0,
                max = 3,
                value = 1,
                step = 0.1
              )
            ),
            column(4,
              sliderInput(
                ns("pvalue_threshold"),
                label = "adj.P-value阈值 / adj.P-value Threshold",
                min = 0.001,
                max = 0.1,
                value = 0.05,
                step = 0.001
              )
            ),
            column(4,
              br(),
              actionButton(
                ns("run_deg"),
                "运行差异分析 / Run Analysis",
                icon = icon("flask"),
                class = "btn-primary",
                width = "100%"
              )
            )
          )
        )
      )
    ),
    
    # 结果统计 / Result Statistics
    fluidRow(
      column(3,
        valueBoxOutput(ns("n_total"), width = 12)
      ),
      column(3,
        valueBoxOutput(ns("n_up"), width = 12)
      ),
      column(3,
        valueBoxOutput(ns("n_down"), width = 12)
      ),
      column(3,
        valueBoxOutput(ns("n_significant"), width = 12)
      )
    ),
    
    # 可视化 / Visualization
    fluidRow(
      column(6,
        box(
          title = "火山图 / Volcano Plot",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          plotlyOutput(ns("volcano_plot"), height = "500px")
        )
      ),
      column(6,
        box(
          title = "MA图 / MA Plot",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          plotlyOutput(ns("ma_plot"), height = "500px")
        )
      )
    ),
    
    fluidRow(
      column(12,
        box(
          title = "差异基因热图 / DEG Heatmap",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          plotOutput(ns("heatmap"), height = "500px")
        )
      )
    ),
    
    # 结果表格 / Result Table
    fluidRow(
      column(12,
        box(
          title = "差异基因表 / DEG Table",
          status = "success",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          fluidRow(
            column(4,
              downloadButton(ns("download_all"), "下载全部 / Download All")
            ),
            column(4,
              downloadButton(ns("download_up"), "下载上调 / Download Up")
            ),
            column(4,
              downloadButton(ns("download_down"), "下载下调 / Download Down")
            )
          ),
          
          hr(),
          
          DT::dataTableOutput(ns("deg_table"))
        )
      )
    )
  )
}


#' 差异分析模块Server
#' Differential Analysis Module Server
#'
#' @param id 模块ID / Module ID
#' @param r 响应式值对象 / Reactive values object
#' @return 响应式值 / Reactive values
#' @export
mod_differential_server <- function(id, r = reactiveValues()) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # 更新分组变量选择 / Update group variable choices
    observe({
      if (!is.null(r$eset_annotated %||% r$eset_normalized %||% r$eset)) {
        eset <- r$eset_annotated %||% r$eset_normalized %||% r$eset
        vars <- varLabels(eset)
        updateSelectInput(session, "group_var", choices = vars)
      }
    })
    
    # 更新分组水平选择 / Update group level choices
    observeEvent(input$group_var, {
      req(input$group_var)
      
      eset <- r$eset_annotated %||% r$eset_normalized %||% r$eset
      req(eset)
      
      if (input$group_var %in% varLabels(eset)) {
        levels <- unique(as.character(pData(eset)[[input$group_var]]))
        updateSelectInput(session, "control_group", choices = levels)
        updateSelectInput(session, "treatment_group", choices = levels)
      }
    })
    
    # 更新分析方法选择 / Update method choices
    observe({
      if (!is.null(r$data_type)) {
        if (r$data_type == "microarray") {
          methods <- c("limma" = "limma")
        } else {
          methods <- c(
            "DESeq2" = "deseq2",
            "edgeR" = "edger",
            "limma-voom" = "limma_voom"
          )
        }
        updateSelectInput(session, "deg_method", choices = methods)
      }
    })
    
    # 运行差异分析 / Run differential analysis
    observeEvent(input$run_deg, {
      req(input$group_var, input$control_group, input$treatment_group)
      
      if (input$control_group == input$treatment_group) {
        showNotification("对照组和实验组不能相同 / Control and treatment groups must differ",
                        type = "error")
        return()
      }
      
      withProgress(message = "运行差异分析 / Running differential analysis...", {
        
        tryCatch({
          incProgress(0.2, detail = "准备数据 / Preparing data...")
          
          if (r$data_type == "microarray") {
            eset <- r$eset_annotated %||% r$eset_normalized %||% r$eset
            req(eset)
            
            source("R/microarray/differential.R", local = TRUE)
            
            contrast <- paste0(input$treatment_group, "-", input$control_group)
            
            incProgress(0.5, detail = "运行limma分析 / Running limma analysis...")
            
            r$deg_result <- run_differential_analysis(
              eset = eset,
              group = input$group_var,
              contrast = contrast,
              lfc_threshold = input$lfc_threshold,
              p_threshold = input$pvalue_threshold,
              output_dir = file.path("data", r$gse_id, "deg")
            )
            
          } else {
            counts <- r$counts_normalized %||% r$counts
            phenotype <- r$download_result$phenotype_data
            req(counts, phenotype)
            
            source("R/rnaseq/differential.R", local = TRUE)
            
            contrast <- c(input$group_var, input$treatment_group, input$control_group)
            
            incProgress(0.5, detail = paste0("运行 ", input$deg_method, " 分析..."))
            
            r$deg_result <- run_rnaseq_differential(
              counts = counts,
              phenotype = phenotype,
              group_col = input$group_var,
              contrast = contrast,
              method = input$deg_method,
              lfc_threshold = input$lfc_threshold,
              p_threshold = input$pvalue_threshold,
              output_dir = file.path("data", r$gse_id, "deg")
            )
          }
          
          incProgress(0.9, detail = "完成 / Complete...")
          
          showNotification("差异分析完成！/ Differential analysis complete!", type = "message")
          
        }, error = function(e) {
          showNotification(paste0("差异分析失败 / Analysis failed: ", conditionMessage(e)),
                          type = "error")
        })
      })
    })
    
    # 值框 / Value boxes
    output$n_total <- renderValueBox({
      n <- if (!is.null(r$deg_result)) r$deg_result$summary$n_total else "-"
      valueBox(n, "总基因 / Total", icon = icon("dna"), color = "blue")
    })
    
    output$n_up <- renderValueBox({
      n <- if (!is.null(r$deg_result)) r$deg_result$summary$n_up else "-"
      valueBox(n, "上调 / Up", icon = icon("arrow-up"), color = "red")
    })
    
    output$n_down <- renderValueBox({
      n <- if (!is.null(r$deg_result)) r$deg_result$summary$n_down else "-"
      valueBox(n, "下调 / Down", icon = icon("arrow-down"), color = "green")
    })
    
    output$n_significant <- renderValueBox({
      n <- if (!is.null(r$deg_result)) r$deg_result$summary$n_significant else "-"
      valueBox(n, "显著 / Significant", icon = icon("star"), color = "purple")
    })
    
    # 火山图 / Volcano plot
    output$volcano_plot <- renderPlotly({
      req(r$deg_result)
      
      deg <- r$deg_result$deg_table
      
      # 统一列名 / Standardize column names
      logfc_col <- if ("logFC" %in% colnames(deg)) "logFC" else "log2FoldChange"
      pval_col <- if ("adj.P.Val" %in% colnames(deg)) "adj.P.Val" else "padj"
      
      volcano_df <- data.frame(
        logFC = deg[[logfc_col]],
        negLogP = -log10(deg[[pval_col]]),
        Status = deg$DEG_status,
        Gene = rownames(deg)
      )
      
      volcano_df$negLogP[is.infinite(volcano_df$negLogP)] <- 
        max(volcano_df$negLogP[is.finite(volcano_df$negLogP)], na.rm = TRUE) + 1
      
      plot_ly(volcano_df, x = ~logFC, y = ~negLogP, 
              color = ~Status,
              colors = c("Down" = "#3366CC", "Not Significant" = "grey60", "Up" = "#CC3366"),
              text = ~Gene,
              type = "scatter", mode = "markers",
              marker = list(size = 5, opacity = 0.6)) %>%
        add_trace(x = c(-input$lfc_threshold, -input$lfc_threshold),
                  y = c(0, max(volcano_df$negLogP, na.rm = TRUE)),
                  type = "scatter", mode = "lines",
                  line = list(dash = "dash", color = "grey"),
                  showlegend = FALSE) %>%
        add_trace(x = c(input$lfc_threshold, input$lfc_threshold),
                  y = c(0, max(volcano_df$negLogP, na.rm = TRUE)),
                  type = "scatter", mode = "lines",
                  line = list(dash = "dash", color = "grey"),
                  showlegend = FALSE) %>%
        add_trace(x = c(min(volcano_df$logFC, na.rm = TRUE), max(volcano_df$logFC, na.rm = TRUE)),
                  y = rep(-log10(input$pvalue_threshold), 2),
                  type = "scatter", mode = "lines",
                  line = list(dash = "dash", color = "grey"),
                  showlegend = FALSE) %>%
        layout(
          title = "Volcano Plot / 火山图",
          xaxis = list(title = "log2 Fold Change"),
          yaxis = list(title = "-log10(adj. P-value)")
        )
    })
    
    # MA图 / MA plot
    output$ma_plot <- renderPlotly({
      req(r$deg_result)
      
      deg <- r$deg_result$deg_table
      
      logfc_col <- if ("logFC" %in% colnames(deg)) "logFC" else "log2FoldChange"
      expr_col <- if ("AveExpr" %in% colnames(deg)) "AveExpr" else "baseMean"
      
      ma_df <- data.frame(
        A = if (expr_col == "baseMean") log2(deg[[expr_col]] + 1) else deg[[expr_col]],
        M = deg[[logfc_col]],
        Status = deg$DEG_status
      )
      
      plot_ly(ma_df, x = ~A, y = ~M,
              color = ~Status,
              colors = c("Down" = "#3366CC", "Not Significant" = "grey60", "Up" = "#CC3366"),
              type = "scatter", mode = "markers",
              marker = list(size = 5, opacity = 0.6)) %>%
        layout(
          title = "MA Plot / MA图",
          xaxis = list(title = "Average Expression"),
          yaxis = list(title = "log2 Fold Change")
        )
    })
    
    # 热图 / Heatmap
    output$heatmap <- renderPlot({
      req(r$deg_result)
      
      deg <- r$deg_result$deg_table
      sig_genes <- rownames(deg)[deg$DEG_status != "Not Significant"]
      
      if (length(sig_genes) == 0) {
        plot.new()
        text(0.5, 0.5, "No significant genes / 无显著基因")
        return()
      }
      
      # 限制基因数 / Limit gene count
      n_genes <- min(50, length(sig_genes))
      logfc_col <- if ("logFC" %in% colnames(deg)) "logFC" else "log2FoldChange"
      top_genes <- head(sig_genes[order(abs(deg[sig_genes, logfc_col]), decreasing = TRUE)], n_genes)
      
      # 获取表达数据 / Get expression data
      if (!is.null(r$eset_annotated %||% r$eset_normalized %||% r$eset)) {
        eset <- r$eset_annotated %||% r$eset_normalized %||% r$eset
        expr <- exprs(eset)
      } else if (!is.null(r$counts_normalized %||% r$counts)) {
        expr <- r$counts_normalized %||% r$counts
      } else {
        return()
      }
      
      # 匹配基因 / Match genes
      common_genes <- intersect(top_genes, rownames(expr))
      
      if (length(common_genes) > 0) {
        expr_subset <- expr[common_genes, , drop = FALSE]
        expr_scaled <- t(scale(t(expr_subset)))
        
        pheatmap::pheatmap(
          expr_scaled,
          color = colorRampPalette(c("blue", "white", "red"))(100),
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          show_rownames = length(common_genes) <= 50,
          main = paste0("Top ", length(common_genes), " DEGs / 差异基因热图")
        )
      }
    })
    
    # DEG表格 / DEG table
    output$deg_table <- DT::renderDataTable({
      req(r$deg_result)
      
      deg <- r$deg_result$deg_table
      deg$Gene <- rownames(deg)
      deg <- deg[, c("Gene", setdiff(colnames(deg), "Gene"))]
      
      DT::datatable(
        deg,
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          order = list(list(2, "desc"))
        ),
        rownames = FALSE
      ) %>%
        DT::formatRound(columns = c(2:5), digits = 3)
    })
    
    # 下载处理 / Download handlers
    output$download_all <- downloadHandler(
      filename = function() paste0(r$gse_id, "_DEG_all.csv"),
      content = function(file) {
        write.csv(r$deg_result$deg_table, file, row.names = TRUE)
      }
    )
    
    output$download_up <- downloadHandler(
      filename = function() paste0(r$gse_id, "_DEG_up.csv"),
      content = function(file) {
        deg <- r$deg_result$deg_table
        write.csv(deg[deg$DEG_status == "Up", ], file, row.names = TRUE)
      }
    )
    
    output$download_down <- downloadHandler(
      filename = function() paste0(r$gse_id, "_DEG_down.csv"),
      content = function(file) {
        deg <- r$deg_result$deg_table
        write.csv(deg[deg$DEG_status == "Down", ], file, row.names = TRUE)
      }
    )
    
    return(r)
  })
}

# 辅助函数 / Helper function
`%||%` <- function(a, b) if (!is.null(a)) a else b
