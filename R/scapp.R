#' app demo for SingleR application
#' @import shiny
#' @import shinytoastr
#' @import shinythemes
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
#' @import PCAtools
#' @note We work just with the MuraroPancreasData, leaving
#' more general upload capability as a project.
#' @export

scapp = function() {
  # derived from ls("package:celldex")
  optfuns = c("None", "BlueprintEncodeData", "DatabaseImmuneCellExpressionData", 
              "HumanPrimaryCellAtlasData", "ImmGenData", "MonacoImmuneData", 
              "MouseRNAseqData", "NovershternHematopoieticData")
  optplots = c("PCA", "elbow", "TSNE", "UMAP")
  
  ui = fluidPage(
    shinytoastr::useToastr(),
    shinythemes::themeSelector(), # theme selection is possible to try out different ones
    theme = shinythemes::shinytheme("sandstone"),  # choice of a single theme once we have decided 
    titlePanel("Diamond24 scapp"),
    sidebarLayout(
      sidebarPanel(
        helpText("App for labeling single cells with selected references"),
        helpText(sprintf("Version %s", packageVersion("diamonds24"))),
        radioButtons("ref", "References", optfuns, selected = "None"),
        numericInput("ncomp", "Number of PCs in PCA plot", min = 2, max = 5, value = 2),
        numericInput("ncomp2", "Number of PCs in UMAP", min = 2, max = 10, value = 2),
        radioButtons("plot", "Plot type to download", optplots, selected = "PCA"),
        downloadButton('downloadPlot', 'Download Plot'),
        fileInput("file", "Upload Data File", accept = c(".rds")),
        helpText("Upload an .rds file containing your data"), width = 4
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("main PCA plot", plotOutput("view")),
          tabPanel("interactive PCA", sidebarLayout(
            mainPanel(plotly::plotlyOutput("viewly")),
            sidebarPanel(checkboxGroupInput('checkly', label = NULL, choices = c()))
          )),
          tabPanel("author", plotOutput("auth")),
          tabPanel("ref comp", verbatimTextOutput("called")),
          tabPanel("ref comp2", verbatimTextOutput("called2")),
          tabPanel("elbow", plotOutput("elbow")),
          tabPanel("tsne", plotOutput("tsne")),
          tabPanel("umap", plotOutput("umap")),
          tabPanel("diagnostics", plotly::plotlyOutput("diagnostics"))
        )
      )
    )
  )
  
  server = function(input, output) {
    # Set the maximum file upload size to 1GB
    options(shiny.maxRequestSize = 1024^3)
    # given = scRNAseq::MuraroPancreasData() 
    
    library(scater)
    library(SingleCellExperiment)
    library(readr)
    library(Seurat)
    library(patchwork)
    
    # Reactive value to store the processed SingleCellExperiment data
    data = reactiveVal(NULL)
    
    observeEvent(input$file, {
      req(input$file)
      file = input$file$datapath
      
      # Read the data
      tryCatch({
        # Assuming the file is in .rds format
        seurat_data = readRDS(file)
        
        # Convert to SingleCellExperiment
        sce_data = as.SingleCellExperiment(seurat_data)
        # sce_data = seurat_data
        
        # Normalize and perform PCA
        sce_data = scuttle::logNormCounts(sce_data)
        sce_data = scater::runPCA(sce_data)
        sce_data = scater::runTSNE(sce_data)
        
        # Store in reactive value
        data(sce_data)
        
        shinytoastr::toastr_success("Data uploaded and processed successfully!")
      }, error = function(e) {
        shinytoastr::toastr_error(paste("Failed to process data:", e$message))
      })
    })
    
    run_SingleR = reactive({
      given = data()
      req(given)
      
      if (input$ref != "None") {
        ref2use = get(input$ref)()
        myb = BiocParallel::MulticoreParam(4)
        shinytoastr::toastr_info("Starting SingleR")
        sing = SingleR::SingleR(given, ref2use, ref2use$label.main, BPPARAM = myb)
        shinytoastr::toastr_info("Done")
        given$celltype = sing$labels
        updateCheckboxGroupInput(
          inputId = 'checkly',
          choices = unique(given$celltype),
          selected = unique(given$celltype)
        )
      }
      given
    })
    
    output$view = renderPlot({
      given = run_SingleR()
      colour_by = ifelse("celltype" %in% colnames(colData(given)), "celltype", "label")
      scater::plotPCA(given, colour_by = colour_by, ncomponents = input$ncomp, theme_size = 14)
    })
    
    output$viewly = plotly::renderPlotly({
      given = run_SingleR()
      given = given[, which(given$celltype %in% input$checkly)]
      dfr = SingleCellExperiment::reducedDim(given)
      mydf = data.frame(PC1 = dfr[, 1], PC2 = dfr[, 2], type = given$celltype)
      gg = ggplot2::ggplot(mydf, aes(x = PC1, y = PC2, text = type, colour = type)) +
        ggplot2::geom_point() +
        ggplot2::theme(legend.position = 'none')
      plotly::ggplotly(gg)
    })
    
    output$auth = renderPlot({
      given = run_SingleR()
      scater::plotPCA(given, colour_by = "label", ncomponents = input$ncomp, theme_size = 14)
    })
    
    output$called = renderPrint({
      if (input$ref != "None") {
        ref2use = get(input$ref)()
        sort(table(ref2use$label.main))
      } else {
        "No reference selected"
      }
    })
    
    output$called2 = renderPrint({
      if (input$ref != "None") {
        ref2use = get(input$ref)()
        sort(table(ref2use$label.fine))
      } else {
        "No reference selected"
      }
    })
    
    output$elbow = renderPlot({
      given = run_SingleR()
      percent.var = attr(SingleCellExperiment::reducedDim(given), "percentVar")
      chosen.elbow = PCAtools::findElbowPoint(percent.var)
      plot(percent.var, xlab = "PC", ylab = "% variance explained", type = "b")  
      abline(v = chosen.elbow, col = "red", lty = 2)
    })
    
    output$tsne = renderPlot({
      given = run_SingleR()
      scater::plotTSNE(given, colour_by = "celltype")
    })
    
    output$umap = renderPlot({
      given = run_SingleR()
      given = scater::runUMAP(given, ncomponents = input$ncomp2)
      scater::plotUMAP(given, colour_by = "celltype")
    })
    
    output$diagnostics = renderPlot({
      given = run_SingleR()
      given = scuttle::addPerCellQC(given, subsets = list(Mt = rowData(given)$featureType == "mito"))
      qc = scuttle::quickPerCellQC(colData(given), sub.fields = "subsets_Mt_percent")
      given$discard = qc$discard
      scater::plotColData(given, x = "sum", y = "subsets_Mt_percent", colour_by = "discard")
    })
    
    output$downloadPlot = downloadHandler(
      filename = function() {
        paste("plot-", input$plot, ".png", sep = "")
      },
      content = function(file) {
        png(file)
        given = run_SingleR()
        if (input$plot == "PCA") {
          scater::plotPCA(given, colour_by = "label", ncomponents = input$ncomp, theme_size = 14)
        } else if (input$plot == "elbow") {
          percent.var = attr(SingleCellExperiment::reducedDim(given), "percentVar")
          chosen.elbow = PCAtools::findElbowPoint(percent.var)
          plot(percent.var, xlab = "PC", ylab = "% variance explained", type = "b")
          abline(v = chosen.elbow, col = "red", lty = 2)
        } else if (input$plot == "TSNE") {
          scater::plotTSNE(given, colour_by = "celltype")
        } else if (input$plot == "UMAP") {
          scater::plotUMAP(given, colour_by = "celltype")
        }
        dev.off()
      }
    )
  }
  
  runApp(list(ui = ui, server = server))
}
