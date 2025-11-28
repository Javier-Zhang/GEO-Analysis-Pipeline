#' ============================================================================
#' GEO Analysis Pipeline - Shiny主应用
#' GEO Analysis Pipeline - Main Shiny Application
#' ============================================================================
#' 
#' 功能说明 / Description:
#' - 完整的GEO数据分析流程 / Complete GEO data analysis pipeline
#' - 支持微阵列和RNA-seq数据 / Support microarray and RNA-seq data
#' - 交互式可视化和报告生成 / Interactive visualization and report generation
#' 
#' 示例数据集（宫颈癌）/ Example Datasets (Cervical Cancer):
#' - GSE6791: 宫颈癌 vs 正常 (GPL570)
#' - GSE7803: 宫颈癌表达谱
#' - GSE9750: 宫颈癌基因表达
#' - GSE63514: 宫颈癌进展研究
#' - GSE39001: HPV相关宫颈癌
#' 
#' @author GEO-Analysis-Pipeline
#' ============================================================================

# 加载必要的包 / Load required packages
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyWidgets)
  library(DT)
  library(plotly)
  library(ggplot2)
  library(Biobase)
})

# 加载Shiny模块 / Load Shiny modules
source("modules/mod_download.R")
source("modules/mod_qc.R")
source("modules/mod_normalize.R")
source("modules/mod_annotate.R")
source("modules/mod_differential.R")
source("modules/mod_report.R")

# UI定义 / UI Definition
ui <- dashboardPage(
  
  # 头部 / Header
  dashboardHeader(
    title = tags$span(
      tags$img(src = "https://www.ncbi.nlm.nih.gov/geo/img/geo_main.gif", height = "30px"),
      "GEO Analysis Pipeline"
    ),
    titleWidth = 350
  ),
  
  # 侧边栏 / Sidebar
  dashboardSidebar(
    width = 280,
    
    # 数据类型选择 / Data type selection
    div(
      style = "padding: 15px;",
      radioGroupButtons(
        inputId = "data_type",
        label = "数据类型 / Data Type",
        choices = c(
          "微阵列 / Microarray" = "microarray",
          "RNA-seq" = "rnaseq"
        ),
        selected = "microarray",
        justified = TRUE,
        size = "sm"
      )
    ),
    
    hr(),
    
    # 导航菜单 / Navigation menu
    sidebarMenu(
      id = "main_menu",
      
      menuItem(
        "数据下载 / Download",
        tabName = "download",
        icon = icon("download")
      ),
      
      menuItem(
        "质量控制 / QC",
        tabName = "qc",
        icon = icon("check-circle")
      ),
      
      menuItem(
        "标准化 / Normalize",
        tabName = "normalize",
        icon = icon("balance-scale")
      ),
      
      menuItem(
        "注释 / Annotate",
        tabName = "annotate",
        icon = icon("tags")
      ),
      
      menuItem(
        "差异分析 / DEG",
        tabName = "differential",
        icon = icon("chart-line")
      ),
      
      menuItem(
        "报告 / Report",
        tabName = "report",
        icon = icon("file-alt")
      ),
      
      hr(),
      
      menuItem(
        "帮助 / Help",
        tabName = "help",
        icon = icon("question-circle")
      )
    ),
    
    # 底部信息 / Footer info
    div(
      style = "position: absolute; bottom: 10px; left: 15px; right: 15px;",
      tags$small(
        style = "color: #95a5a6;",
        "GEO-Analysis-Pipeline v1.0",
        br(),
        "Cervical Cancer Research / 宫颈癌研究"
      )
    )
  ),
  
  # 主体内容 / Body
  dashboardBody(
    
    # 引入自定义CSS / Include custom CSS
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
      tags$style(HTML("
        .content-wrapper {
          background-color: #f5f5f5;
        }
        .main-header .logo {
          font-size: 16px;
        }
      "))
    ),
    
    # 标签页内容 / Tab content
    tabItems(
      
      # 下载页面 / Download page
      tabItem(
        tabName = "download",
        mod_download_ui("download")
      ),
      
      # QC页面 / QC page
      tabItem(
        tabName = "qc",
        mod_qc_ui("qc")
      ),
      
      # 标准化页面 / Normalization page
      tabItem(
        tabName = "normalize",
        mod_normalize_ui("normalize")
      ),
      
      # 注释页面 / Annotation page
      tabItem(
        tabName = "annotate",
        mod_annotate_ui("annotate")
      ),
      
      # 差异分析页面 / Differential analysis page
      tabItem(
        tabName = "differential",
        mod_differential_ui("differential")
      ),
      
      # 报告页面 / Report page
      tabItem(
        tabName = "report",
        mod_report_ui("report")
      ),
      
      # 帮助页面 / Help page
      tabItem(
        tabName = "help",
        fluidRow(
          column(12,
            box(
              title = "使用说明 / User Guide",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              
              h4("GEO数据分析流程 / GEO Data Analysis Pipeline"),
              
              tags$ol(
                tags$li(
                  tags$strong("数据下载 / Download: "),
                  "输入GSE ID，选择数据类型，下载GEO数据"
                ),
                tags$li(
                  tags$strong("质量控制 / QC: "),
                  "进行样本质控，检测异常样本，生成QC图表"
                ),
                tags$li(
                  tags$strong("标准化 / Normalization: "),
                  "选择合适的标准化方法处理数据"
                ),
                tags$li(
                  tags$strong("注释 / Annotation: "),
                  "将探针/基因ID转换为标准基因符号"
                ),
                tags$li(
                  tags$strong("差异分析 / DEG: "),
                  "设置分组，进行差异表达分析"
                ),
                tags$li(
                  tags$strong("报告 / Report: "),
                  "生成并下载分析报告"
                )
              ),
              
              hr(),
              
              h4("示例数据集（宫颈癌）/ Example Datasets (Cervical Cancer)"),
              
              tableOutput("example_table"),
              
              hr(),
              
              h4("参考资料 / References"),
              
              tags$ul(
                tags$li(
                  tags$a(
                    href = "https://www.ncbi.nlm.nih.gov/geo/",
                    target = "_blank",
                    "NCBI GEO Database"
                  )
                ),
                tags$li(
                  tags$a(
                    href = "https://bioconductor.org/",
                    target = "_blank",
                    "Bioconductor"
                  )
                ),
                tags$li(
                  tags$a(
                    href = "https://www.bioconductor.org/packages/release/bioc/html/limma.html",
                    target = "_blank",
                    "limma Package"
                  )
                ),
                tags$li(
                  tags$a(
                    href = "https://bioconductor.org/packages/release/bioc/html/DESeq2.html",
                    target = "_blank",
                    "DESeq2 Package"
                  )
                )
              )
            )
          )
        )
      )
    )
  )
)


# Server定义 / Server Definition
server <- function(input, output, session) {
  
  # 共享响应式值 / Shared reactive values
  r <- reactiveValues(
    gse_id = NULL,
    data_type = NULL,
    download_result = NULL,
    eset = NULL,
    eset_normalized = NULL,
    eset_annotated = NULL,
    counts = NULL,
    counts_normalized = NULL,
    qc_results = NULL,
    annotation_result = NULL,
    deg_result = NULL
  )
  
  # 监听数据类型变化 / Listen to data type changes
  observeEvent(input$data_type, {
    r$data_type <- input$data_type
  })
  
  # 调用模块服务器 / Call module servers
  mod_download_server("download", r)
  mod_qc_server("qc", r)
  mod_normalize_server("normalize", r)
  mod_annotate_server("annotate", r)
  mod_differential_server("differential", r)
  mod_report_server("report", r)
  
  # 示例数据表 / Example data table
  output$example_table <- renderTable({
    data.frame(
      GSE_ID = c("GSE6791", "GSE7803", "GSE9750", "GSE63514", "GSE39001"),
      Platform = c("GPL570", "GPL570", "GPL96", "GPL570", "GPL570"),
      Samples = c(42, 28, 66, 128, 108),
      Description = c(
        "Cervical cancer vs Normal / 宫颈癌vs正常",
        "Cervical cancer expression / 宫颈癌表达谱",
        "Cervical cancer genes / 宫颈癌基因表达",
        "Cervical cancer progression / 宫颈癌进展",
        "HPV-associated / HPV相关宫颈癌"
      ),
      stringsAsFactors = FALSE
    )
  })
}


# 运行应用 / Run application
shinyApp(ui = ui, server = server)
