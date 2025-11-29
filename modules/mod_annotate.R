#' ============================================================================
#' Shiny注释模块
#' Shiny Annotation Module
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 注释方法选择 / Annotation method selection
#' - 结果预览表格 / Result preview table
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

#' 注释模块UI
#' Annotation Module UI
#'
#' @param id 模块ID / Module ID
#' @return tagList，UI元素 / tagList, UI elements
#' @export
mod_annotate_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(12,
        box(
          title = "注释设置 / Annotation Settings",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          fluidRow(
            column(4,
              selectInput(
                ns("annot_method"),
                label = "注释方法 / Method",
                choices = c(
                  "Bioconductor (org.db)" = "bioconductor",
                  "biomaRt (Ensembl)" = "biomart",
                  "GPL平台文件 / GPL Platform" = "gpl"
                ),
                selected = "bioconductor"
              )
            ),
            column(4,
              selectInput(
                ns("organism"),
                label = "物种 / Organism",
                choices = c(
                  "人类 / Human" = "human",
                  "小鼠 / Mouse" = "mouse",
                  "大鼠 / Rat" = "rat"
                ),
                selected = "human"
              )
            ),
            column(4,
              selectInput(
                ns("id_type"),
                label = "目标ID类型 / Target ID Type",
                choices = c(
                  "Gene Symbol" = "SYMBOL",
                  "Entrez ID" = "ENTREZID",
                  "Ensembl ID" = "ENSEMBL"
                ),
                selected = "SYMBOL"
              )
            )
          ),
          
          hr(),
          
          fluidRow(
            column(6,
              selectInput(
                ns("collapse_method"),
                label = "多探针合并方法 / Probe Collapse Method",
                choices = c(
                  "最大IQR / Max IQR" = "maxIQR",
                  "平均值 / Mean" = "mean",
                  "中位数 / Median" = "median"
                ),
                selected = "maxIQR"
              )
            ),
            column(6,
              br(),
              actionButton(
                ns("run_annot"),
                "运行注释 / Run Annotation",
                icon = icon("tags"),
                class = "btn-primary",
                width = "100%"
              )
            )
          )
        )
      )
    ),
    
    # 注释统计 / Annotation Statistics
    fluidRow(
      column(4,
        valueBoxOutput(ns("total_probes"), width = 12)
      ),
      column(4,
        valueBoxOutput(ns("annotated_probes"), width = 12)
      ),
      column(4,
        valueBoxOutput(ns("unique_genes"), width = 12)
      )
    ),
    
    # 注释预览 / Annotation Preview
    fluidRow(
      column(12,
        box(
          title = "注释预览 / Annotation Preview",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          DT::dataTableOutput(ns("annotation_table"))
        )
      )
    ),
    
    # 注释质量 / Annotation Quality
    fluidRow(
      column(6,
        box(
          title = "注释率统计 / Annotation Rate",
          status = "success",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          plotlyOutput(ns("annotation_rate_plot"), height = "300px")
        )
      ),
      column(6,
        box(
          title = "探针-基因映射 / Probe-Gene Mapping",
          status = "success",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          plotlyOutput(ns("mapping_stats_plot"), height = "300px")
        )
      )
    )
  )
}


#' 注释模块Server
#' Annotation Module Server
#'
#' @param id 模块ID / Module ID
#' @param r 响应式值对象 / Reactive values object
#' @return 响应式值 / Reactive values
#' @export
mod_annotate_server <- function(id, r = reactiveValues()) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # 运行注释 / Run annotation
    observeEvent(input$run_annot, {
      
      withProgress(message = "运行注释 / Running annotation...", {
        
        tryCatch({
          incProgress(0.2, detail = "准备数据 / Preparing data...")
          
          if (r$data_type == "microarray") {
            req(r$eset_normalized %||% r$eset)
            
            eset <- r$eset_normalized %||% r$eset
            
            source("R/microarray/annotation.R", local = TRUE)
            
            # 获取平台信息 / Get platform info
            platform <- tryCatch(annotation(eset), error = function(e) NULL)
            
            incProgress(0.5, detail = "查询注释 / Querying annotation...")
            
            annot_result <- annotate_microarray(
              eset = eset,
              method = input$annot_method,
              platform = platform,
              organism = input$organism,
              id_type = input$id_type,
              collapse_method = input$collapse_method
            )
            
            r$annotation_result <- annot_result
            r$eset_annotated <- annot_result$eset
            
          } else {
            req(r$counts_normalized %||% r$counts)
            
            counts <- r$counts_normalized %||% r$counts
            
            source("R/rnaseq/annotation.R", local = TRUE)
            
            incProgress(0.5, detail = "查询注释 / Querying annotation...")
            
            # 检测ID类型 / Detect ID type
            from_type <- detect_gene_id_type(rownames(counts))
            
            annot_result <- annotate_rnaseq(
              counts = counts,
              from_type = from_type,
              to_type = input$id_type,
              organism = input$organism,
              method = input$annot_method
            )
            
            r$annotation_result <- annot_result
          }
          
          incProgress(0.9, detail = "完成 / Complete...")
          
          showNotification("注释完成！/ Annotation complete!", type = "message")
          
        }, error = function(e) {
          showNotification(paste0("注释失败 / Annotation failed: ", conditionMessage(e)),
                          type = "error")
        })
      })
    })
    
    # 值框：总探针数 / Value box: Total probes
    output$total_probes <- renderValueBox({
      n <- "-"
      if (!is.null(r$annotation_result)) {
        if (is.list(r$annotation_result) && "annotation" %in% names(r$annotation_result)) {
          n <- nrow(r$annotation_result$annotation)
        } else {
          n <- nrow(r$annotation_result)
        }
      }
      valueBox(n, "总探针/基因 / Total", icon = icon("dna"), color = "blue")
    })
    
    # 值框：已注释数 / Value box: Annotated
    output$annotated_probes <- renderValueBox({
      n <- "-"
      if (!is.null(r$annotation_result)) {
        annot <- if (is.list(r$annotation_result) && "annotation" %in% names(r$annotation_result)) {
          r$annotation_result$annotation
        } else {
          r$annotation_result
        }
        
        symbol_col <- "SYMBOL"
        if (symbol_col %in% colnames(annot)) {
          n <- sum(!is.na(annot[[symbol_col]]) & annot[[symbol_col]] != "")
        }
      }
      valueBox(n, "已注释 / Annotated", icon = icon("check"), color = "green")
    })
    
    # 值框：唯一基因数 / Value box: Unique genes
    output$unique_genes <- renderValueBox({
      n <- "-"
      if (!is.null(r$annotation_result)) {
        annot <- if (is.list(r$annotation_result) && "annotation" %in% names(r$annotation_result)) {
          r$annotation_result$annotation
        } else {
          r$annotation_result
        }
        
        symbol_col <- "SYMBOL"
        if (symbol_col %in% colnames(annot)) {
          symbols <- annot[[symbol_col]]
          n <- length(unique(symbols[!is.na(symbols) & symbols != ""]))
        }
      }
      valueBox(n, "唯一基因 / Unique Genes", icon = icon("fingerprint"), color = "purple")
    })
    
    # 注释表格 / Annotation table
    output$annotation_table <- DT::renderDataTable({
      req(r$annotation_result)
      
      annot <- if (is.list(r$annotation_result) && "annotation" %in% names(r$annotation_result)) {
        r$annotation_result$annotation
      } else {
        r$annotation_result
      }
      
      # 显示前1000行 / Show first 1000 rows
      annot_display <- head(annot, 1000)
      
      DT::datatable(
        annot_display,
        options = list(
          pageLength = 15,
          scrollX = TRUE
        ),
        rownames = FALSE
      )
    })
    
    # 注释率图 / Annotation rate plot
    output$annotation_rate_plot <- renderPlotly({
      req(r$annotation_result)
      
      annot <- if (is.list(r$annotation_result) && "annotation" %in% names(r$annotation_result)) {
        r$annotation_result$annotation
      } else {
        r$annotation_result
      }
      
      total <- nrow(annot)
      annotated <- sum(!is.na(annot$SYMBOL) & annot$SYMBOL != "", na.rm = TRUE)
      not_annotated <- total - annotated
      
      plot_ly(
        labels = c("Annotated / 已注释", "Not Annotated / 未注释"),
        values = c(annotated, not_annotated),
        type = "pie",
        marker = list(colors = c("#27ae60", "#e74c3c"))
      ) %>%
        layout(title = "Annotation Rate / 注释率")
    })
    
    # 映射统计图 / Mapping stats plot
    output$mapping_stats_plot <- renderPlotly({
      req(r$annotation_result)
      
      annot <- if (is.list(r$annotation_result) && "annotation" %in% names(r$annotation_result)) {
        r$annotation_result$annotation
      } else {
        r$annotation_result
      }
      
      if (!"SYMBOL" %in% colnames(annot)) {
        return(plotly_empty())
      }
      
      # 计算每个基因对应的探针数 / Count probes per gene
      symbols <- annot$SYMBOL[!is.na(annot$SYMBOL) & annot$SYMBOL != ""]
      probe_counts <- table(symbols)
      
      count_dist <- table(probe_counts)
      
      plot_ly(
        x = as.numeric(names(count_dist)),
        y = as.numeric(count_dist),
        type = "bar",
        marker = list(color = "steelblue")
      ) %>%
        layout(
          title = "Probes per Gene / 每基因探针数",
          xaxis = list(title = "Probes / 探针数"),
          yaxis = list(title = "Genes / 基因数")
        )
    })
    
    return(r)
  })
}

# 辅助函数：空值合并 / Helper: null coalescing
`%||%` <- function(a, b) if (!is.null(a)) a else b
