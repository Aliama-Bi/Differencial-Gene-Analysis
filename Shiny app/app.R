
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
                        min = 0, max = 5000, value = c(0,5000), step = 500)
            
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(
                tabPanel("data table",
                         DT::dataTableOutput("table")
                ),
                tabPanel("MA Plot output",
                         plotOutput("MA_PlotImage")
                ),
                tabPanel("volcano output",
                         plotOutput("volcanoImage") 
                ),
                tabPanel("heat map output",
                         plotOutput("heatmapImage")
                ),
                tabPanel("RF boxplot",
                         plotOutput("rfboxplot")
                )
            ),
            
        )
    )
)

server <- function(input, output) {
    
    
    
    output$table <- DT::renderDataTable({
        temp = temp[input$range[1]:input$range[2],]
        DT::datatable(temp, options = list(lengthMenu = c(25, 25, 5), pageLength = 10))
        
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
    
    
    output$MA_PlotImage <- renderPlot({
        temp = temp[input$range[1]:input$range[2],]
        groupname <- c(unlist(rep('ACR',16)),unlist(rep('Normal',18)))
        design <- model.matrix(~ groupname + 0)  ## design matrix
        gse = log2(temp[,2:length(temp)]+1)
        fit <- lmFit(gse, design) 
        cont.matrix <- makeContrasts(groupnameACR-groupnameNormal, levels=design)  ## pairwise comparison
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2) 
        df<- topTable(fit2, number=nrow(fit2), genelist=rownames(gse))
        
        ggplot(df, aes(x = AveExpr, y = logFC))+
            geom_point(aes(colour=-log10(P.Value)), alpha=1/3, size=1) +
            scale_colour_gradient(low="blue",high="red")+
            ylab("log2 fold change") + xlab("Average expression")
            
            
    })
    
    output$rfboxplot <- renderPlot({
        temp = temp[input$range[1]:input$range[2],]
        X = t(temp)[2:35,]
        y = c(unlist(rep('ACR',16)),unlist(rep('Normal',18)))
        cvK = 5  # number of CV folds
        cv_50acc5_knn = cv_50acc5_svm = cv_50acc5_rf = c()
        cv_acc_knn = cv_acc_svm = cv_acc_rf = c()
        n_sim = 5 ## number of repeats
        
        for (i in 1:n_sim) {
            
            cvSets = cvTools::cvFolds(nrow(X), cvK)  # permute all the data, into 5 folds
            cv_acc_knn = cv_acc_svm = cv_acc_rf = c()
            
            for (j in 1:cvK) {
                test_id = cvSets$subsets[cvSets$which == j]
                X_test = X[test_id, ]
                X_train = X[-test_id, ]
                y_test = y[test_id]
                y_train = y[-test_id]
                
                ## RandomForest
                rf_res <- randomForest::randomForest(x = X_train, y = as.factor(y_train))
                fit <- predict(rf_res, X_test)
                cv_acc_rf[j] = table(fit, y_test) %>% diag %>% sum %>% `/`(length(y_test))
            }
            cv_50acc5_rf <- append(cv_50acc5_rf, mean(cv_acc_rf))
        } ## end for
        RandomForest = as.data.frame(cv_50acc5_rf)
        boxplot(RandomForest, main = "Mean Accuracy of Random Forest")
    })
    
    output$volcanoImage <- renderPlot({
        temp = temp[input$range[1]:input$range[2],]
        groupname <- c(unlist(rep('ACR',16)),unlist(rep('Normal',18)))
        design <- model.matrix(~ groupname + 0)  ## design matrix
        gse = log2(temp[,2:length(temp)]+1)
        fit <- lmFit(gse, design) 
        cont.matrix <- makeContrasts(groupnameACR-groupnameNormal, levels=design)  ## pairwise comparison
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2) 
        df<- topTable(fit2, number=nrow(fit2), genelist=rownames(gse))
        ggplot(df, aes(logFC,-log10(P.Value)))+
            geom_point(aes(colour=-log10(P.Value)), alpha=1/3, size=1) +
            scale_colour_gradient(low="blue",high="red")+
            xlab("log2 fold change") + ylab("-log10 p-value")
    })
    
}


# Run the application 
shinyApp(ui = ui, server = server)
