#' app demo for SingleR application
#' @import shiny
#' @import shinytoastr
#' @import celldex
#' @importFrom SingleR SingleR
#' @importFrom scRNAseq MuraroPancreasData
#' @import scran
#' @import scater
#' @import scuttle
#' @import BiocParallel
#' @import Seurat
#' @import patchwork
#' @import readr
#' @note We work just with the MuraroPancreasData, leaving
#' more general upload capability as a project.
#' @export
scapp <- function() {
  # derived from ls("package:celldex")
  optfuns <- c("None", "BlueprintEncodeData", "DatabaseImmuneCellExpressionData", 
               "HumanPrimaryCellAtlasData", "ImmGenData", "MonacoImmuneData", 
               "MouseRNAseqData", "NovershternHematopoieticData")
  
  ui <- fluidPage(
    shinytoastr::useToastr(),
    sidebarLayout(
      sidebarPanel(
        helpText("App for labeling single cells with selected references"),
        radioButtons("ref", "refs", optfuns, selected = "None"),
        numericInput("ncomp", "npcs", min = 2, max = 5, value = 2),
        fileInput("file", "Upload Data File", accept = c(".rds")),
        helpText("Upload an .rds file containing your data"), width = 2
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("main", plotOutput("view")),
          tabPanel("interact", plotly::plotlyOutput("viewly")),
          tabPanel("author", plotOutput("auth")),
          tabPanel("ref comp", verbatimTextOutput("called"))
        )
      )
    )
  )
  
  server <- function(input, output) {
    # Set the maximum file upload size to 1GB
    options(shiny.maxRequestSize = 1024^3)
    
    library(scater)
    library(SingleCellExperiment)
    library(readr)
    library(Seurat)
    library(patchwork)
    
    # Reactive value to store the processed SingleCellExperiment data
    data <- reactiveVal(NULL)
    
    observeEvent(input$file, {
      req(input$file)
      file <- input$file$datapath
      
      # Read the data
      tryCatch({
        # Assuming the file is in .rds format
        seurat_data <- readRDS(file)
        
        # Convert to SingleCellExperiment
        sce_data <- as.SingleCellExperiment(seurat_data)
        
        # Normalize and perform PCA
        sce_data <- scuttle::logNormCounts(sce_data)
        sce_data <- scater::runPCA(sce_data)
        
        # Store in reactive value
        data(sce_data)
        
        shinytoastr::toastr_success("Data uploaded and processed successfully!")
      }, error = function(e) {
        shinytoastr::toastr_error(paste("Failed to process data:", e$message))
      })
    })
    
    run_SingleR <- reactive({
      given <- data()
      req(given)
      
      if (input$ref != "None") {
        ref2use <- get(input$ref)()
        myb <- BiocParallel::MulticoreParam(4)
        shinytoastr::toastr_info("Starting SingleR")
        sing <- SingleR::SingleR(given, ref2use, ref2use$label.main, BPPARAM = myb)
        shinytoastr::toastr_info("Done")
        given$celltype <- sing$labels
      }
      given
    })
    
    output$view <- renderPlot({
      given <- run_SingleR()
      colour_by <- ifelse("celltype" %in% colnames(colData(given)), "celltype", "label")
      scater::plotPCA(given, colour_by = colour_by, ncomponents = input$ncomp, theme_size = 14)
    })
    
    output$viewly <- plotly::renderPlotly({
      given <- run_SingleR()
      dfr <- SingleCellExperiment::reducedDim(given)
      colour_by <- ifelse("celltype" %in% colnames(colData(given)), "celltype", "label")
      mydf <- data.frame(PC1 = dfr[,1], PC2 = dfr[,2], type = colData(given)[[colour_by]])
      gg <- ggplot2::ggplot(mydf, aes(x = PC1, y = PC2, text = type, colour = type)) +
        ggplot2::geom_point()
      plotly::ggplotly(gg)
    })
    
    output$auth <- renderPlot({
      scater::plotPCA(data(), colour_by = "label", ncomponents = input$ncomp, theme_size = 14)
    })
    
    output$called <- renderPrint({
      if (input$ref != "None") {
        ref2use <- get(input$ref)()
        sort(table(ref2use$label.main))
      } else {
        "No reference selected"
      }
    })
  }
  
  runApp(list(ui = ui, server = server))
}
