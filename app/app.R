
library(shiny)
library(shinythemes)
library(GEOquery)  ## go to https://github.com/seandavi/GEOquery for installation details
library(R.utils)
library(reshape2)
library(ggplot2)
library(limma) ## to look at linear model in matrix annotation
library(dplyr)
library(ggbiplot)
library(factoextra)
library(tidyverse)
library(DT)
library(grDevices)
library(BiocManager)
options(repos = BiocManager::repositories())


temp = read.csv('GSE131179_expression_matrix.csv')


ui <- fluidPage(
    theme = shinytheme("flatly"),
    tags$title("Differential gene expression analysis"),
    
    div(
        tags$header(p("Differential gene expression analysis", style="font-size:40px"),
                    p("SID:480139690", style="font-size:30px")),
        align = "center", style="color:#ffffff; background-color: #4d728d") 
    ,
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            
            sliderInput("range",
                        label = "Size of Data",
                        min = 0, max = 200, value = c(0,200), step = 10)
            
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(
                tabPanel("data table",
                         DT::dataTableOutput("table")
                ),
                tabPanel("heat map output",
                         plotOutput("heatmapImage")
                )
            ),
            
        )
    )
)

server <- function(input, output) {
    
    
    
    output$table <- DT::renderDataTable({
        temp = temp[input$range[1]:input$range[2],]
        DT::datatable(temp, options = list(lengthMenu = c(5, 10, 20), pageLength = 10))
        
    })
    
    
    output$heatmapImage <- renderPlot({
        temp = temp[input$range[1]:input$range[2],]
        groupname <- c(unlist(rep('ACR',16)),unlist(rep('Normal',18)))
        design <- model.matrix(~ groupname + 0)  ## design matrix
        gse = log2(temp[,2:length(temp)]+1)
        fit <- lmFit(gse, design) 
        cont.matrix <- makeContrasts(groupnameACR-groupnameNormal, levels=design)  ## pairwise comparison
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2) 
        df<- topTable(fit2, number=nrow(fit2), genelist=rownames(gse))
        sig50 <- c(order(fit2$F,decreasing=TRUE)[1:50])
        label = temp[sig50,1]
        heatmap(as.matrix(gse[sig50,]),labRow =label)
    })
    
    
}


# Run the application 
shinyApp(ui = ui, server = server)
