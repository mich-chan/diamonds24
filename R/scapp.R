
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
#' @note We work just with the MuraroPancreasData, leaving
#' more general upload capability as a project.
#' @export
scapp = function() {
# derived from ls("package:celldex")
 optfuns = c("BlueprintEncodeData", "DatabaseImmuneCellExpressionData", 
"HumanPrimaryCellAtlasData", "ImmGenData", "MonacoImmuneData", 
"MouseRNAseqData", "NovershternHematopoieticData")
 optplots = c("PCA", "elbow","TSNE", "UMAP")
 ui = fluidPage(
  shinytoastr::useToastr(),
  shinythemes::themeSelector(), # theme selection is possible to try out different ones
  theme = shinythemes::shinytheme("sandstone"),  # choice of a single theme once we have decided 
  titlePanel("Diamond24 scapp") ,
  sidebarLayout(
   sidebarPanel(
    helpText("app for labeling single cells with selected references"),
    helpText(sprintf("version %s", packageVersion("diamonds24"))),
    radioButtons("ref", "refs", optfuns, selected="HumanPrimaryCellAtlasData"),
# consider option for label.main, label.fine
    numericInput("ncomp", "number of PCs in PCA plot", min=2, max=5, value=2),
    numericInput("ncomp2", "number of PCs in UMAP", min=2, max=10, value=2),
    radioButtons("plot", "plot type to download", optplots, selected="PCA"),
    downloadButton('downloadPlot', 'Download Plot'),
    helpText("provide an upload function here"), width=4
    ),
    mainPanel(
     tabsetPanel(
      tabPanel("main PCA plot", plotOutput("view")),
      tabPanel("interactive PCA", plotly::plotlyOutput("viewly")),
      tabPanel("author", plotOutput("auth")),
      tabPanel("ref comp", verbatimTextOutput("called")),
      tabPanel("ref comp2", verbatimTextOutput("called2")),
      tabPanel("elbow", plotOutput("elbow")),
      tabPanel("tsne", plotOutput("tsne")),
      tabPanel("umap", plotOutput("umap")),
      tabPanel("diagnostics", plotlyOutput("diagnostics"))
               
     )
    )
   )
  )
 server = function(input, output) {
  # needs help here -- 1) use upload method, 2) verify gene symbols present, get from rowData if not
  # setup, runPCA once, defines "given"
   given = scRNAseq::MuraroPancreasData() 
   rownames(given) = rowData(given)$symbol
   dups = which(duplicated(rownames(given)))
   if (length(dups)>0) given = given[-dups,]
   given = scuttle::logNormCounts(given)
   given = scater::runPCA(given)
   given <- scater::runTSNE(given)
  #
  run_SingleR = reactive({
   ref2use = get(input$ref)()
   myb = BiocParallel::MulticoreParam(4)
   shinytoastr::toastr_info("starting SingleR")
   sing = SingleR::SingleR(given, ref2use, ref2use$label.main, BPPARAM=myb)
   shinytoastr::toastr_info("done")
   given$celltype = sing$labels
   #scater::runPCA(given)
   given
   })
  output$view = renderPlot({
   given = run_SingleR()
   scater::plotPCA(given, colour_by = "celltype", 
        ncomponents=input$ncomp, theme_size=14)
   })
  output$viewly = plotly::renderPlotly({
   given = run_SingleR()
   dfr = SingleCellExperiment::reducedDim(given)
   mydf = data.frame(PC1 = dfr[,1], PC2=dfr[,2], type=given$celltype)
   gg = ggplot2::ggplot(mydf, aes(x=PC1, y=PC2, text=type,
      colour=type)) +
     ggplot2::geom_point() +
     ggplot2::theme_minimal()
   plotly::ggplotly(gg)
   })
  output$auth = renderPlot({
  scater::plotPCA(given, colour_by = "label",
        ncomponents=input$ncomp, theme_size=14)
   })
  output$called = renderPrint({
   ref2use = get(input$ref)()
   sort(table(ref2use$label.main))
   })
  output$called2 = renderPrint({
   ref2use = get(input$ref)()
   sort(table(ref2use$label.fine))
   })
  output$elbow = renderPlot({
    percent.var <- attr(SingleCellExperiment::reducedDim(given),"percentVar")
    chosen.elbow <- findElbowPoint(percent.var)
    plot(percent.var,  xlab="PC", ylab="% variance explained", type="b")  
    abline(v=chosen.elbow, col="red", lty=2)
    #scater::plotElbow(given)
    })
  output$tsne = renderPlot({
   given = run_SingleR()
   scater::plotTSNE(given, colour_by = "celltype")
   })
  output$umap = renderPlot({
   given = run_SingleR()
   given <- scater::runUMAP(given, ncomponents=input$ncomp2)
   scater::plotUMAP(given, colour_by = "celltype")
   
   })
  output$diagnostics = renderPlot({
   given <- scuttle::addPerCellQC(given, subsets = list(Mt=rowData(given)$featureType=="mito"))
   qc <- scuttle::quickPerCellQC(colData(given), sub.fields="subsets_Mt_percent")
  given$discard <- qc$discard
  scater::plotColData(given, x="sum", y="subsets_Mt_percent",colour_by = "discard")
   })
  output$downloadPlot = downloadHandler(
   filename = function() {
    paste("plot-", input$plot, ".png", sep="")
   },
   content = function(file) {
    png(file)
    if (input$plot == "PCA") {
      scater::plotPCA(given, colour_by = "label",
                      ncomponents=input$ncomp, theme_size=14)
    } else if (input$plot == "elbow") {
      plot(percent.var,  xlab="PC", ylab="% variance explained", type="b")  
            abline(v=chosen.elbow, col="red", lty=2)
    } else if (input$plot == "TSNE") {
      scater::plotTSNE(given, colour_by = "celltype")
    } else if (input$plot == "UMAP") {
      scater::plotUMAP(given, colour_by = "celltype")
    }
    dev.off()
   }
   
   )
  
 }
 runApp(list(ui=ui, server=server))
}
