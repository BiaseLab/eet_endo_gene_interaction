library("shiny")
library("RCy3")
library("rcytoscapejs")
library("DT")
library("ggplot2")
library("ggpubr")

options(DT.options = list(pageLength = 5, rownames= FALSE))


#load data

#LOAD NECESSARY FILES
d18_EET<-read.delim("2017_06_17_fpkm_eet_day_18_AI_annotated.txt.bz2" ,header = TRUE, sep = "\t", stringsAsFactors =FALSE, row.names=1,
                    colClasses=c("character", "numeric","numeric","numeric","numeric","numeric","character"))
d18_CAR<-read.delim("2017_06_17_fpkm_endo_day_18_C_AI_annotated.txt.bz2" ,header = TRUE, sep = "\t", stringsAsFactors =FALSE, row.names=1,
                    colClasses=c("character", "numeric","numeric","numeric","numeric","numeric","character"))
d18_ICAR<-read.delim("2017_06_17_fpkm_endo_day_18_IC_AI_annotated.txt.bz2" ,header = TRUE, sep = "\t", stringsAsFactors =FALSE, row.names=1,
                     colClasses=c("character", "numeric","numeric","numeric","numeric","numeric","character"))
cor_d18_EET_CAR_a<-read.delim("2017-12-28_pearson_correlation_COR_PVALUE_EET_CAR_0.95_shiny.txt.bz2" ,header = TRUE, sep = "\t",stringsAsFactors=FALSE,comment.char="#",
                            colClasses=c( "character","character","numeric","character","character"))
cor_d18_EET_ICAR_a<-read.delim("2017-12-28_pearson_correlation_COR_PVALUE_EET_ICAR_0.95_shiny.txt.bz2" ,header = TRUE, sep = "\t", stringsAsFactors =FALSE ,comment.char="#",
                             colClasses=c( "character","character","numeric","character","character"))
day_18_EET_genes<-sort(unique(c(as.character(cor_d18_EET_CAR_a$external_gene_name.eet), as.character(cor_d18_EET_ICAR_a$external_gene_name.eet))))
day_18_CAR_EET_genes<-sort(unique(as.character(cor_d18_EET_CAR_a$external_gene_name.end)))
day_18_ICAR_EET_genes<-sort(unique(as.character(cor_d18_EET_ICAR_a$external_gene_name.end)))


map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

# Define UI ----
ui <- fluidPage( titlePanel("Inter tissue gene co-expression - Extra-embryonic tissue and endometrium"),
                      
                      sidebarLayout(
                        sidebarPanel(
                          
                          h5("Extraembryonic tissue (EET)"),
                          h5("Caruncle tissue (CAR)"),
                          h5("Intercaruncle tissue (ICAR)"),
                          
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
                          
                          ,
                          
                          br(),
                          
                          sliderInput("cor_threshold_pos",
                                      "Set positive correlation threshold",
                                      value = 0.97,
                                      step = 0.01,
                                      round = -2,
                                      min = 0.95,
                                      max = 1)
                          
                          ,
                          
                          br(),
                          sliderInput("cor_threshold_neg",
                                      "Set negative correlation threshold",
                                      value = -0.97,
                                      step = 0.01,
                                      round = -2,
                                      min = -1,
                                      max = -0.95)
                          
                          ,
                          
                          br(),
                          
                          radioButtons("choiceScatterPlot", "Create Scaterplots?",
                                       c("EET  - CAR" =  "S_P_EET_CAR",
                                         "EET  - ICAR"=  "S_P_EET_ICAR",
                                         "CAR  - EET" =  "S_P_CAR_EET",
                                         "ICAR - EET "=  "S_P_ICAR_EET"),
                                       selected = "S_P_EET_CAR")
                        , width=3),
                        
                        mainPanel(
                          
                          tabsetPanel(type = "tabs",
                                      
                                      tabPanel("Network", 
                                               
                                               fluidRow(
                                                 column(6,
                                                        verbatimTextOutput("selected_EET"),
                                                        
                                                        rcytoscapejsOutput("plotEET_CAR", height="300px",width ="300px")
                                                        
                                                       ),
                                                 
                                                 column(6,
                                                        verbatimTextOutput("selected_EET2"),
                                                        
                                                        rcytoscapejsOutput("plotEET_ICAR", height="300px",width ="300px")
                                                        )
                                                        ),
                                               br(),
                                               br(),
                                               fluidRow(
                                                 column(6,
                                                        verbatimTextOutput("selected_CAR"),
                                                        
                                                        rcytoscapejsOutput("plotCAR", height="300px",width ="300px")
                                                        ),
                                                 
                                                 column(6,
                                                        
                                                        verbatimTextOutput("selected_ICAR"),
                                                        
                                                        rcytoscapejsOutput("plotICAR", height="300px",width ="300px")
                                                       )
                                                       )
                                                       ),
                                      
                                      tabPanel("Tables",
                                               fluidRow(
                                                 column(6,
                                                        verbatimTextOutput("table_selected_EET"),
                                                        
                                                        DTOutput("DataTable_EET_CAR", height = "30%")
                                                 ),
                                                 
                                                 column(6,
                                                        verbatimTextOutput("table_selected_EET2"),
                                                        
                                                        DTOutput("DataTable_EET_ICAR", height = "30%")
                                                 )
                                               ),
                                               br(),
                                               br(),
                                               fluidRow(
                                                 column(6,
                                                        verbatimTextOutput("table_selected_CAR"),
                                                        
                                                        DTOutput("DataTable_CAR_EET", height = "30%")
                                                 ),
                                                 
                                                 column(6,
                                                        
                                                        verbatimTextOutput("table_selected_ICAR"),
                                                        
                                                        DTOutput("DataTable_ICAR_EET", height = "30%")
                                                 )
                                                 )
                                                 ),
                                      tabPanel("Scatterplots", 
                                               verbatimTextOutput("text_tab_scatterplot"),
                                               conditionalPanel(
                                                 condition = "input.choiceScatterPlot == 'S_P_EET_CAR'", plotOutput("scatterplot1")),
                                               conditionalPanel(
                                                 condition = "input.choiceScatterPlot == 'S_P_EET_ICAR'", plotOutput("scatterplot2")),
                                               conditionalPanel(
                                                 condition = "input.choiceScatterPlot == 'S_P_CAR_EET'", plotOutput("scatterplot3")),
                                               conditionalPanel(
                                                 condition = "input.choiceScatterPlot == 'S_P_ICAR_EET'", plotOutput("scatterplot4"))
                                               )
                                      
                          )
                        , width=9)
                        
  ),
  
  hr(),
  div(class = "footer", includeHTML("footer.html"))
)

# Define server logic ----
server <- function(input, output) {
  
  
  output$selected_EET <- renderText({ 
    paste("Network of ", input$varEET, " expressed in EET correlated with genes expressed in CAR", sep="")
  })
  
  output$selected_EET2 <- renderText({ 
    paste("Network of ", input$varEET," expressed in EET correlated with genes expressed in ICAR", sep="")
  })
  
  output$selected_CAR <- renderText({ 
    paste("Network of ", input$varCAR," expressed in CAR correlated with genes expressed in EET ", sep="")
  })
  
  output$selected_ICAR <- renderText({ 
    paste("Network of ", input$varICAR," expressed in ICAR correlated with genes expressed in EET", sep="")
  })
  
  output$selected_ICAR <- renderText({ 
    paste("Network of ", input$varICAR," expressed in ICAR correlated with genes expressed in EET", sep="")
  })
  
  
  output$table_selected_EET <- renderText({ 
    paste("Table of ", input$varEET, " expressed in EET correlated with genes expressed in  correlated with genes", sep="")
  })
  
  output$table_selected_EET2 <- renderText({ 
    paste("Table of ", input$varEET," expressed in EET correlated with genes expressed in ICAR", sep="")
  })
  
  output$table_selected_CAR <- renderText({ 
    paste("Table of ", input$varCAR," expressed in CAR correlated with genes expressed in EET ", sep="")
  })
  
  output$text_tab_scatterplot <- renderText({ 
    paste("Consider displaying 10 or less plots ;) . Each dot represents data from one pregnancy. On y-axis you have the genes expressed in EET and on x-axis you have the genes expressed in the endometrium")
  })
  
  output$plotEET_CAR <- renderRcytoscapejs({
    
    cor_d18_EET_CAR_a<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$external_gene_name.eet == input$varEET, c(4,5,3)]
    cor_d18_EET_CAR_b<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$pearson_cor >= input$cor_threshold_pos, ]
    cor_d18_EET_CAR_c<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$pearson_cor <= input$cor_threshold_neg, ]
    
  network<-rbind(cor_d18_EET_CAR_b,cor_d18_EET_CAR_c)
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
    
    cor_d18_EET_CAR_a<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$external_gene_name.eet == input$varEET, c(4,5,3)]
    cor_d18_EET_CAR_b<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$pearson_cor >= input$cor_threshold_pos, ]
    cor_d18_EET_CAR_c<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$pearson_cor <= input$cor_threshold_neg, ]
    
    DataTable_EET_CAR<-rbind(cor_d18_EET_CAR_b,cor_d18_EET_CAR_c)
    colnames(DataTable_EET_CAR)<-c("EET gene symbol", "CAR gene symbol", "correlation")
    return(DataTable_EET_CAR)
  })
  
  
  
  output$plotEET_ICAR <- renderRcytoscapejs({
    
    cor_d18_EET_ICAR_a<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$external_gene_name.eet == input$varEET, c(4,5,3)]
    cor_d18_EET_ICAR_b<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$pearson_cor >= input$cor_threshold_pos, ]
    cor_d18_EET_ICAR_c<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$pearson_cor <= input$cor_threshold_neg, ]
    
    network<-rbind(cor_d18_EET_ICAR_b,cor_d18_EET_ICAR_c)
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
    cor_d18_EET_ICAR_a<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$external_gene_name.eet == input$varEET, c(4,5,3)]
    cor_d18_EET_ICAR_b<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$pearson_cor >= input$cor_threshold_pos, ]
    cor_d18_EET_ICAR_c<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$pearson_cor <= input$cor_threshold_neg, ]
    
    DataTable_EET_ICAR<-rbind(cor_d18_EET_ICAR_b,cor_d18_EET_ICAR_c)
    colnames(DataTable_EET_ICAR)<-c("EET gene symbol", "ICAR gene symbol", "correlation")
    return(DataTable_EET_ICAR)
  })

  
  output$plotCAR <- renderRcytoscapejs({
    
    cor_d18_EET_CAR_a<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$external_gene_name.end == input$varCAR, c(4,5,3)]
    cor_d18_EET_CAR_b<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$pearson_cor >= input$cor_threshold_pos, ]
    cor_d18_EET_CAR_c<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$pearson_cor <= input$cor_threshold_neg, ]
    
    network<-rbind(cor_d18_EET_CAR_b,cor_d18_EET_CAR_c)
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
    
    cor_d18_EET_CAR_a<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$external_gene_name.end == input$varCAR, c(4,5,3)]
    cor_d18_EET_CAR_b<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$pearson_cor >= input$cor_threshold_pos, ]
    cor_d18_EET_CAR_c<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$pearson_cor <= input$cor_threshold_neg, ]
    
  
    DataTable_CAR_EET<-rbind(cor_d18_EET_CAR_b,cor_d18_EET_CAR_c)
    colnames(DataTable_CAR_EET)<-c("EET gene symbol", "CAR gene symbol", "correlation")
    return(DataTable_CAR_EET)
  })
  
  output$plotICAR <- renderRcytoscapejs({
    
    cor_d18_EET_ICAR_a<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$external_gene_name.end == input$varICAR, c(4,5,3)]
    cor_d18_EET_ICAR_b<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$pearson_cor >= input$cor_threshold_pos, ]
    cor_d18_EET_ICAR_c<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$pearson_cor <= input$cor_threshold_neg, ]
    
    network<-rbind(cor_d18_EET_ICAR_b,cor_d18_EET_ICAR_c)
    
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
    
    cor_d18_EET_ICAR_a<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$external_gene_name.end == input$varICAR, c(4,5,3)]
    cor_d18_EET_ICAR_b<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$pearson_cor >= input$cor_threshold_pos, ]
    cor_d18_EET_ICAR_c<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$pearson_cor <= input$cor_threshold_neg, ]
    
    DataTable_ICAR_EET<-rbind(cor_d18_EET_ICAR_b,cor_d18_EET_ICAR_c)
    colnames(DataTable_ICAR_EET)<-c("EET gene symbol", "ICAR gene symbol", "correlation")
    return(DataTable_ICAR_EET)
  })
  
  output$scatterplot1 <- renderPlot({
      
      cor_d18_EET_CAR_a<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$external_gene_name.eet == input$varEET, c(4,5,3)]
      cor_d18_EET_CAR_b<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$pearson_cor >= input$cor_threshold_pos, ]
      cor_d18_EET_CAR_c<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$pearson_cor <= input$cor_threshold_neg, ]
      cor_table<-rbind(cor_d18_EET_CAR_b,cor_d18_EET_CAR_c)
    
    data_chart_1<-data.frame( stringsAsFactors=FALSE)
    data_chart_2<-data.frame( stringsAsFactors=FALSE)
    
    for (i in seq(dim(cor_table)[1])){
      
      gene_eet<-cor_table[i,1]
      gene_endo<-cor_table[i,2]
      
      fpkm_gene_eet<-d18_EET[d18_EET$external_gene_name==gene_eet,]
      fpkm_gene_endo<-d18_CAR[d18_CAR$external_gene_name==gene_endo,]
      
      if( (dim(fpkm_gene_eet)[1]==1) & (dim(fpkm_gene_endo)[1]==1) ){
      
      data_chart_1<-data.frame(t(fpkm_gene_eet[1:5]))
      data_chart_1<-cbind(data_chart_1,data.frame(t(fpkm_gene_endo[1:5])))
      data_chart_1$chart<-i
      data_chart_1$gene1<-gene_eet
      data_chart_1$gene2<-gene_endo
      
      colnames(data_chart_1)<-c("gene_EET", "gene_ENDO", "chart", "gene1","gene2")
      
      data_chart_2<-rbind(data_chart_2,data_chart_1)
      }
    }

    plots <- list()
    k<-1
    for (j in unique(data_chart_2$chart)){
      data_chart_3<-data_chart_2[data_chart_2$chart %in% j , ]
      plot<-ggplot(data=data_chart_3, aes(x=gene_ENDO,y=gene_EET ))+
        geom_point(size=0.4)+
        scale_y_continuous(name=data_chart_3$gene1[1])+
        scale_x_continuous(name=data_chart_3$gene2[1])+
        theme(aspect.ratio = 1,
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background  = element_blank(),
              axis.text.x = element_text( colour = 'black' ,size = 7),
              axis.text.y = element_text( colour = 'black',size = 7),
              axis.title= element_text( colour = 'black' ,size = 7, face="italic"),
              axis.ticks = element_line(size=0.1),
              panel.spacing = unit(0, "mm"),
              legend.position="none",
              axis.line=element_line(size = 0.1, colour = "black")
        )
      
      plots[[k]] <- plot
      k<-k+1
    }
    
    ggarrange( plotlist=plots)
    
  })
 
  
  output$scatterplot2 <- renderPlot({
    
    cor_d18_EET_ICAR_a<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$external_gene_name.eet == input$varEET, c(4,5,3)]
    cor_d18_EET_ICAR_b<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$pearson_cor >= input$cor_threshold_pos, ]
    cor_d18_EET_ICAR_c<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$pearson_cor <= input$cor_threshold_neg, ]
    cor_table<-rbind(cor_d18_EET_ICAR_b,cor_d18_EET_ICAR_c)
    
    data_chart_1<-data.frame( stringsAsFactors=FALSE)
    data_chart_2<-data.frame( stringsAsFactors=FALSE)
    
    for (i in seq(dim(cor_table)[1])){
      
      gene_eet<-cor_table[i,1]
      gene_endo<-cor_table[i,2]
      
      fpkm_gene_eet<-d18_EET[d18_EET$external_gene_name==gene_eet,]
      fpkm_gene_endo<-d18_ICAR[d18_ICAR$external_gene_name==gene_endo,]
      
      if( (dim(fpkm_gene_eet)[1]==1) & (dim(fpkm_gene_endo)[1]==1) ){
        
        data_chart_1<-data.frame(t(fpkm_gene_eet[1:5]))
        data_chart_1<-cbind(data_chart_1,data.frame(t(fpkm_gene_endo[1:5])))
        data_chart_1$chart<-i
        data_chart_1$gene1<-gene_eet
        data_chart_1$gene2<-gene_endo
        
        colnames(data_chart_1)<-c("gene_EET", "gene_ENDO", "chart", "gene1","gene2")
        
        data_chart_2<-rbind(data_chart_2,data_chart_1)
      }
    }
    
    plots <- list()
    k<-1
    for (j in unique(data_chart_2$chart)){
      data_chart_3<-data_chart_2[data_chart_2$chart %in% j , ]
      plot<-ggplot(data=data_chart_3, aes(x=gene_ENDO,y=gene_EET ))+
        geom_point(size=0.4)+
        scale_y_continuous(name=data_chart_3$gene1[1])+
        scale_x_continuous(name=data_chart_3$gene2[1])+
        theme(aspect.ratio = 1,
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background  = element_blank(),
              axis.text.x = element_text( colour = 'black' ,size = 7),
              axis.text.y = element_text( colour = 'black',size = 7),
              axis.title= element_text( colour = 'black' ,size = 7, face="italic"),
              axis.ticks = element_line(size=0.1),
              panel.spacing = unit(0, "mm"),
              legend.position="none",
              axis.line=element_line(size = 0.1, colour = "black")
        )
      
      plots[[k]] <- plot
      k<-k+1
    }
    
    ggarrange( plotlist=plots)
    
  })
  
  
  
  output$scatterplot3 <- renderPlot({
    
    cor_d18_EET_CAR_a<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$external_gene_name.end == input$varCAR, c(4,5,3)]
    cor_d18_EET_CAR_b<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$pearson_cor >= input$cor_threshold_pos, ]
    cor_d18_EET_CAR_c<-cor_d18_EET_CAR_a[cor_d18_EET_CAR_a$pearson_cor <= input$cor_threshold_neg, ]
    cor_table<-rbind(cor_d18_EET_CAR_b,cor_d18_EET_CAR_c)
    
    data_chart_1<-data.frame( stringsAsFactors=FALSE)
    data_chart_2<-data.frame( stringsAsFactors=FALSE)
    
    for (i in seq(dim(cor_table)[1])){
      
      gene_eet<-cor_table[i,1]
      gene_endo<-cor_table[i,2]
      
      fpkm_gene_eet<-d18_EET[d18_EET$external_gene_name==gene_eet,]
      fpkm_gene_endo<-d18_CAR[d18_CAR$external_gene_name==gene_endo,]
      
      if( (dim(fpkm_gene_eet)[1]==1) & (dim(fpkm_gene_endo)[1]==1) ){
        
        data_chart_1<-data.frame(t(fpkm_gene_eet[1:5]))
        data_chart_1<-cbind(data_chart_1,data.frame(t(fpkm_gene_endo[1:5])))
        data_chart_1$chart<-i
        data_chart_1$gene1<-gene_eet
        data_chart_1$gene2<-gene_endo
        
        colnames(data_chart_1)<-c("gene_EET", "gene_ENDO", "chart", "gene1","gene2")
        
        data_chart_2<-rbind(data_chart_2,data_chart_1)
      }
    }
    
    plots <- list()
    k<-1
    for (j in unique(data_chart_2$chart)){
      data_chart_3<-data_chart_2[data_chart_2$chart %in% j , ]
      plot<-ggplot(data=data_chart_3, aes(x=gene_ENDO,y=gene_EET ))+
        geom_point(size=0.4)+
        scale_y_continuous(name=data_chart_3$gene1[1])+
        scale_x_continuous(name=data_chart_3$gene2[1])+
        theme(aspect.ratio = 1,
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background  = element_blank(),
              axis.text.x = element_text( colour = 'black' ,size = 7),
              axis.text.y = element_text( colour = 'black',size = 7),
              axis.title= element_text( colour = 'black' ,size = 7, face="italic"),
              axis.ticks = element_line(size=0.1),
              panel.spacing = unit(0, "mm"),
              legend.position="none",
              axis.line=element_line(size = 0.1, colour = "black")
        )
      
      plots[[k]] <- plot
      k<-k+1
    }
    
    ggarrange( plotlist=plots)
    
  })
  
  
  output$scatterplot4 <- renderPlot({
    
    cor_d18_EET_ICAR_a<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$external_gene_name.end == input$varICAR, c(4,5,3)]
    cor_d18_EET_ICAR_b<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$pearson_cor >= input$cor_threshold_pos, ]
    cor_d18_EET_ICAR_c<-cor_d18_EET_ICAR_a[cor_d18_EET_ICAR_a$pearson_cor <= input$cor_threshold_neg, ]
    cor_table<-rbind(cor_d18_EET_ICAR_b,cor_d18_EET_ICAR_c)
    
    data_chart_1<-data.frame( stringsAsFactors=FALSE)
    data_chart_2<-data.frame( stringsAsFactors=FALSE)
    
    for (i in seq(dim(cor_table)[1])){
      
      gene_eet<-cor_table[i,1]
      gene_endo<-cor_table[i,2]
      
      fpkm_gene_eet<-d18_EET[d18_EET$external_gene_name==gene_eet,]
      fpkm_gene_endo<-d18_ICAR[d18_ICAR$external_gene_name==gene_endo,]
      
      if( (dim(fpkm_gene_eet)[1]==1) & (dim(fpkm_gene_endo)[1]==1) ){
        
        data_chart_1<-data.frame(t(fpkm_gene_eet[1:5]))
        data_chart_1<-cbind(data_chart_1,data.frame(t(fpkm_gene_endo[1:5])))
        data_chart_1$chart<-i
        data_chart_1$gene1<-gene_eet
        data_chart_1$gene2<-gene_endo
        
        colnames(data_chart_1)<-c("gene_EET", "gene_ENDO", "chart", "gene1","gene2")
        
        data_chart_2<-rbind(data_chart_2,data_chart_1)
      }
    }
    
    plots <- list()
    k<-1
    for (j in unique(data_chart_2$chart)){
      data_chart_3<-data_chart_2[data_chart_2$chart %in% j , ]
      plot<-ggplot(data=data_chart_3, aes(x=gene_ENDO,y=gene_EET ))+
        geom_point(size=0.4)+
        scale_y_continuous(name=data_chart_3$gene1[1])+
        scale_x_continuous(name=data_chart_3$gene2[1])+
        theme(aspect.ratio = 1,
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background  = element_blank(),
              axis.text.x = element_text( colour = 'black' ,size = 7),
              axis.text.y = element_text( colour = 'black',size = 7),
              axis.title= element_text( colour = 'black' ,size = 7, face="italic"),
              axis.ticks = element_line(size=0.1),
              panel.spacing = unit(0, "mm"),
              legend.position="none",
              axis.line=element_line(size = 0.1, colour = "black")
        )
      
      plots[[k]] <- plot
      k<-k+1
    }
    
    ggarrange( plotlist=plots)
    
  })  
  
}

# Run the app ----
shinyApp(ui = ui, server = server)
