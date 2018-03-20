library("shiny")
library("RCy3")
library("rcytoscapejs")
library("DT")
library("ggplot2")

options(DT.options = list(pageLength = 5))


#load data
setwd("/Users/fernandobiase/Documents/AUBURN/research/eet_endo_interaction/shiny_app/eet_endo")

#LOAD NECESSARY FILES
d18_EET<-read.delim("2017_06_17_fpkm_d18_EET_filtered.txt.bz2" ,header = TRUE, sep = "\t", stringsAsFactors =FALSE)
d18_CAR<-read.delim("2017_06_17_fpkm_d18_ENDO_CAR_filtered.txt.bz2" ,header = TRUE, sep = "\t", stringsAsFactors =FALSE)
d18_ICAR<-read.delim("2017_06_17_fpkm_d18_ENDO_ICAR_filtered.txt.bz2" ,header = TRUE, sep = "\t", stringsAsFactors =FALSE)
cor_d18_EET_CAR_a<-read.delim("2017-12-28_pearson_correlation_COR_PVALUE_EET_CAR_0.95_shyni.txt.bz2" ,header = TRUE, sep = "\t",stringsAsFactors=FALSE,comment.char="#",
                            colClasses=c( "character","character","numeric","character","character"))
cor_d18_EET_ICAR_a<-read.delim("2017-12-28_pearson_correlation_COR_PVALUE_EET_ICAR_0.95_shyni.txt.bz2" ,header = TRUE, sep = "\t", stringsAsFactors =FALSE ,comment.char="#",
                             colClasses=c( "character","character","numeric","character","character"))
day_18_EET_genes<-sort(unique(c(as.character(cor_d18_EET_CAR_a$external_gene_name.eet), as.character(cor_d18_EET_ICAR_a$external_gene_name.eet))))
day_18_CAR_EET_genes<-sort(unique(as.character(cor_d18_EET_CAR_a$external_gene_name.end)))
day_18_ICAR_EET_genes<-sort(unique(as.character(cor_d18_EET_ICAR_a$external_gene_name.end)))


map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

# Define UI ----
ui <- fluidPage(
  
  sidebarLayout(
    sidebarPanel(
      
      h5("Extraembryonic tissue (EET)"),
      h5("Caruncle tissue (EET)"),
      h5("Intercaruncle tissue (EET)"),
      
      h3("EET"),
      
       selectInput("varEET", 
                   label = "Choose a gene to display a network",
                   choices = day_18_EET_genes,
                   selected = "IFN-tau-c1",
                   selectize=FALSE),
       
    h3("CAR"),
    
    selectInput("varCAR", 
                label = "Choose a gene to display a network",
                choices = day_18_CAR_EET_genes,
                selected = "GRK3",
                selectize=FALSE),
   
    
    h3("ICAR"),
    
    selectInput("varICAR", 
                label = "Choose a gene to display a network",
                choices = day_18_ICAR_EET_genes,
                selected = "AHSP",
                selectize=FALSE)
    
   
               
    ),
               
  mainPanel(
            
    fluidRow(
            column(6,
            verbatimTextOutput("selected_EET"),
        
            rcytoscapejsOutput("plotEET_CAR", height="300px",width ="300px"),
            
            DTOutput("DataTable_EET_CAR", height = "30%")
            ),
            
            column(6,
            verbatimTextOutput("selected_EET2"),
                   
            rcytoscapejsOutput("plotEET_ICAR", height="300px",width ="300px"),    
                   
            DTOutput("DataTable_EET_ICAR", height = "30%")
            )
            ),
    br(),
    br(),
    fluidRow(
      column(6,
            verbatimTextOutput("selected_CAR"),
            
            rcytoscapejsOutput("plotCAR", height="300px",width ="300px"),
            
            DTOutput("DataTable_CAR_EET", height = "30%")
            ),
  
           column(6,
                  
            verbatimTextOutput("selected_ICAR"),
            
            rcytoscapejsOutput("plotICAR", height="300px",width ="300px"),
            
            DTOutput("DataTable_ICAR_EET", height = "30%")
           )
           )
         )
    )
)

# Define server logic ----
server <- function(input, output) {
  
  
  output$selected_EET <- renderText({ 
    paste("Network of ", input$varEET, " expressed in EET with genes expressed in CAR", sep="")
  })
  
  output$selected_EET2 <- renderText({ 
    paste("Network of ", input$varEET," expressed in EET with genes expressed in ICAR", sep="")
  })
  
  output$selected_CAR <- renderText({ 
    paste("Network of ", input$varCAR," expressed in CAR with genes expressed in EET ", sep="")
  })
  
  output$selected_ICAR <- renderText({ 
    paste("Network of ", input$varICAR," expressed in ICAR with genes expressed in EET", sep="")
  })
  
  output$plotEET_CAR <- renderRcytoscapejs({
  network<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$external_gene_name.eet == input$varEET, c(4,5,3)]
  mypal <- colorRampPalette( c( 'blue', 'red' ) )( 100 )
  id <- unique(c(as.character(network$external_gene_name.eet), as.character(network$external_gene_name.end)))
  name <- id
  name <- id
  nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
  nodeData$color<-ifelse(nodeData$name %in% input$varEET, "#CC79A7","#009E73"  )
  source <- network$external_gene_name.eet
  target <-  network$external_gene_name.end
  color<-map2color(network$pearson_cor,mypal,c(-1,1))
  edgeData <- data.frame(source, target, color, stringsAsFactors=FALSE)  
  cyNetwork <- createCytoscapeJsNetwork(nodeData, edgeData, nodeColor = "#888888",  edgeSourceShape = "none", edgeTargetShape = "none")
  rcytoscapejs(cyNetwork$nodes, cyNetwork$edges, showPanzoom=TRUE)
     })
  
  
  output$DataTable_EET_CAR <- renderDT({
    DataTable_EET_CAR<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$external_gene_name.eet == input$varEET, c(4,5,3)]
    colnames(DataTable_EET_CAR)<-c("EET gene symbol", "CAR gene symbol", "correlation")
    return(DataTable_EET_CAR)
  })
  
  
  
  output$plotEET_ICAR <- renderRcytoscapejs({
    network<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$external_gene_name.eet == input$varEET, c(4,5,3)]
    mypal <- colorRampPalette( c( 'blue', 'red' ) )( 100 )
    id <- unique(c(as.character(network$external_gene_name.eet), as.character(network$external_gene_name.end)))
    name <- id
    name <- id
    nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
    nodeData$color<-ifelse(nodeData$name %in% input$varEET, "#CC79A7","#0072B2"  )
    source <- network$external_gene_name.eet
    target <-  network$external_gene_name.end
    color<-map2color(network$pearson_cor,mypal,c(-1,1))
    edgeData <- data.frame(source, target, color, stringsAsFactors=FALSE)  
    cyNetwork <- createCytoscapeJsNetwork(nodeData, edgeData, nodeColor = "#888888",  edgeSourceShape = "none", edgeTargetShape = "none")
    rcytoscapejs(cyNetwork$nodes, cyNetwork$edges, showPanzoom=TRUE)
  })

  output$DataTable_EET_ICAR <- renderDT({
    DataTable_EET_ICAR<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$external_gene_name.eet == input$varEET, c(4,5,3)]
    colnames(DataTable_EET_ICAR)<-c("EET gene symbol", "ICAR gene symbol", "correlation")
    return(DataTable_EET_ICAR)
  })
  
  
  output$plotCAR <- renderRcytoscapejs({
    network<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$external_gene_name.end == input$varCAR, c(4,5,3)]
    mypal <- colorRampPalette( c( 'blue', 'red' ) )( 100 )
    id <- unique(c(as.character(network$external_gene_name.eet), as.character(network$external_gene_name.end)))
    name <- id
    name <- id
    nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
    nodeData$color<-ifelse(nodeData$name %in% input$varCAR, "#009E73", "#CC79A7"  )
    source <- network$external_gene_name.eet
    target <-  network$external_gene_name.end
    color<-map2color(network$pearson_cor,mypal,c(-1,1))
    edgeData <- data.frame(source, target, color, stringsAsFactors=FALSE)  
    cyNetwork <- createCytoscapeJsNetwork(nodeData, edgeData, nodeColor = "#888888",  edgeSourceShape = "none", edgeTargetShape = "none")
    rcytoscapejs(cyNetwork$nodes, cyNetwork$edges, showPanzoom=TRUE)
  })
  
  output$DataTable_CAR_EET <- renderDT({
    DataTable_CAR_EET<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$external_gene_name.end == input$varCAR, c(4,5,3)]
    colnames(DataTable_CAR_EET)<-c("EET gene symbol", "CAR gene symbol", "correlation")
    return(DataTable_CAR_EET)
  })
  
  output$plotICAR <- renderRcytoscapejs({
    network<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$external_gene_name.end == input$varICAR, c(4,5,3)]
    mypal <- colorRampPalette( c( 'blue', 'red' ) )( 100 )
    id <- unique(c(as.character(network$external_gene_name.eet), as.character(network$external_gene_name.end)))
    name <- id
    name <- id
    nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
    nodeData$color<-ifelse(nodeData$name %in% input$varICAR, "#0072B2", "#CC79A7"  )
    source <- network$external_gene_name.eet
    target <-  network$external_gene_name.end
    color<-map2color(network$pearson_cor,mypal,c(-1,1))
    edgeData <- data.frame(source, target, color, stringsAsFactors=FALSE)  
    cyNetwork <- createCytoscapeJsNetwork(nodeData, edgeData, nodeColor = "#888888",  edgeSourceShape = "none", edgeTargetShape = "none")
    rcytoscapejs(cyNetwork$nodes, cyNetwork$edges, showPanzoom=TRUE)
  })
  
  
  output$DataTable_ICAR_EET <- renderDT({
    DataTable_ICAR_EET<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$external_gene_name.end == input$varICAR, c(4,5,3)]
    colnames(DataTable_ICAR_EET)<-c("EET gene symbol", "ICAR gene symbol", "correlation")
    return(DataTable_ICAR_EET)
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)
