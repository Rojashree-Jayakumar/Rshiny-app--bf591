## Author: Rojashree Jayakumar
## rshreej@bu.edu
## BU BF591
## Final Project
## R shiny app to view summary, counts normalization, differential expression and gene set enrichment analysis results (tables and plots)
## Tested using the data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150450  

library(shiny)
library(ggplot2)
library(colourpicker)
library(DT)
library('tidyverse')
library(ggbeeswarm)
library("ggplot2")
library("gridExtra") 
library('tidyverse')
library(ggbeeswarm)
library("ggplot2")
library("gridExtra")  
library(DESeq2)
options(shiny.maxRequestSize=30*1024^2) 

# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("Final Project"),
    p("To use this application, download the CSV metadata.csv from the data directory of this app's repository."),
    tabsetPanel(
        tabPanel("Samples",
                 sidebarPanel(
                     fileInput("file","Load metadata", accept = ".csv"),
                     actionButton("submit", label = "Submit", width = "400px",
                                  style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                                  style = "material-flat")
                 ),
                 mainPanel(
                     tabsetPanel(
                         tabPanel("Summary",tableOutput("summary")),
                         tabPanel("Table", DT::dataTableOutput("table")),
                         tabPanel("Plot",
                                  sidebarPanel(
                                      radioButtons("button", "Choose the column for the x-axis",
                                                   choices = c("lifestage", "timepoint", "treatment", "sex")),
                                      actionButton("plot_submit", label = "Plot", width = "230px",
                                                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                                                   style = "material-flat")
                                  ),
                                  mainPanel(plotOutput("histogram"))
                         )
                     )
                 )
        ),
        tabPanel("Counts",
                 sidebarPanel(
                     fileInput("file2", "Load normalised counts",accept = ".csv"),
                     sliderInput(inputId = "slider1", min = 0, max = 100,
                                 label = "Select the threshold for non-zero genes", value = 0, step = 10),
                     sliderInput(inputId = "slider2", min = 1, max = 100,
                                 label = "Select the threshold for variance", value = 0, step = 10)
                     #actionButton("submit2", label = "Submit", width = "400px",
                     #            style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                     #           style = "material-flat")
                 ),
                 mainPanel(
                     tabsetPanel(
                         tabPanel("Filter",
                                  sidebarPanel(
                                      actionButton("summary_submit", label = "Display summary statistics", width = "230px",
                                                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                                                   style = "material-flat")
                                  ),
                                  mainPanel(
                                      tableOutput("filter_counts")
                                  )
                         ),
                         tabPanel("Diagnostic_scatter_plots", label = "Diagnostic scatter plots",
                                  sidebarPanel(
                                      actionButton("plot2", label = "Plot", width = "230px",
                                                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                                                   style = "material-flat")
                                  ),
                                  mainPanel(
                                      plotOutput("plot2"),
                                      plotOutput("plot3")
                                  )
                         )
                     )
                 )
        )
    )
)
    

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    #tab 1 data
    
    load_data <- reactive({
        out<-read.csv(input$file$datapath, stringsAsFactors = TRUE)
        return(out)
    })
    
    #tab 1.1 summary
    
    summary_table<-function (data) {
        numeric_data<-data[,c(3,4,7)]
        mean<-t(apply(numeric_data,2,mean))
        sd<-t(apply(numeric_data,2,sd))
        summary<-colnames(numeric_data)
        type<-rep("numeric", 3)
        mean<-(paste0(mean," +/- ",sd))
        summary<-as_tibble(cbind(summary,type,mean))
        colnames(summary)<-c("Metadata", "type", "Mean (sd) or Distinct Values (frequency)")
        #2.processing factor data
        #subsetting only factor columns
        factor_data<-data[,-c(1,3,4,6,7,13,14,23)]
        levels(factor_data$sex) <- c("female","male","unknown")
        factor_data$sex[factor_data$sex == ''] <- 'unknown'
        y<-data.frame(summary(factor_data))
        m<-as_tibble(t(y%>%drop_na()))
        m<- m[-c(1,2),]
        m$lifestage <- paste(m[[12]], ",", m[[13]])
        m$sex<- paste(m[[28]], ",", m[[29]] , ",", m[[27]])
        m$treatment<-paste(m[[25]], ",", m[[26]])
        m$timepoint<- paste(m[[19]], ",", m[[20]], ",", m[[21]], ",",m[[22]], ",",m[[23]], ",",m[[24]])
        m<-m[,-c(12,13,19:29)]
        colnames<-c("Assay.Type","BioProject","Center.Name","Consent","DATASTORE.filetype","DATASTORE.provider", "DATASTORE.region","Instrument", "LibraryLayout", "LibrarySelection", "LibrarySource", "Organism", "Platform", "ReleaseDate", "source_name", "SRA.Study", "lifestage", "sex", "treatment", "timepoint" )
        type<-rep("factor", 20)
        summary_n<-as_tibble(t(rbind(colnames,type,m)))
        colnames(summary_n)<-c("Metadata", "type", "Mean (sd) or Distinct Values (frequency)")
        #3.processing character data
        character_data<-data[,c(1,6,13,14,23)]
        summary_c<-as_tibble(cbind(colnames(character_data),rep("character",5),rep(NA,5)))
        colnames(summary_c)<-c("Metadata", "type", "Mean (sd) or Distinct Values (frequency)")
        summary_final<-rbind(summary_n,summary,summary_c)
        
        return (summary_final)
    }
    
    # tab 1.3  boxplots
    
    plot_bases<-function(data,x_name){
        levels(data$sex) <- c("female","male", "unknown")
        data$sex[data$sex == ''] <- 'unknown'
        out<-ggplot(data) +
            geom_boxplot(aes(x=!!sym(x_name),y= Bases, fill= !!sym(x_name)))+ 
            theme_grey()+
            ggtitle("Boxplot of bases")
        return(out)
    }
    
    #tab 2.1
    
    #tab 2 data load
    
    load_norm_counts <- reactive({
        out<-read_csv(input$file2$datapath)
        return(out)
    })
    
    #tab 2 filter
    #tab 2.1
    non_zero_func<- function (data){
        return (sum(data!=0))
    }
    
    filter<-function (norm_counts, slider1, slider2) {
        non_zero <- apply(norm_counts[,-1], 1, non_zero_func)
        variance <- apply(norm_counts[,-1], 1, var)
        non_zero_filtered<-norm_counts[non_zero>slider1 & rank(variance)>(17491*slider2)/100,]
        out <- tibble(
            `number of genes` = 17491,
            `number of samples` = 81,
            `number of genes passing filter` = nrow(non_zero_filtered),
            `number of genes not passing filter` = 17491 - nrow(non_zero_filtered),
            `percentage of genes passing filter` = (nrow(non_zero_filtered)/17491)*100,
            `percentage of genes not passing filter` = ((17491 - nrow(non_zero_filtered))/17491)*100,
        )
        return(out)
    }
    
    #tab 2.2
    plot_variance_vs_median <- function(data, slider) {
        median <- apply(data[,-1], 1, median)
        variances <- apply(data[,-1], 1, var)
        tbl <- tibble(
            Median = median,
            Variance = variances
        )
        tbl$Rank <- rank(tbl$Median)
        tbl$filter<-case_when(rank(tbl$Variance)>(17491*slider)/100 ~ 'genes passing the filter', rank(tbl$Variance)<(17491*slider)/100 ~ 'genes not passing filter')
        tbl<-tbl[tbl$Rank>870,] #filtering extremely low ranks
        out <- ggplot2::ggplot(tbl, aes(x=Rank, y=Variance, color = filter)) +
            geom_point(alpha=0.5) +
            geom_smooth(method='gam', formula = y ~ s(x, bs = "cs")) +
            xlab("Rank(median)") +
            ylab("Variance") +
            theme_bw()+
            ggtitle('Variance vs Median, Log10(Normalized counts)') +
            scale_y_log10()
        return (out)
    }
    
    null_count<- function (data){
        return (sum(data==0))
    }
    
    plot_median_vs_null <- function(data, slider) {
        median <- apply(data[,-1], 1, median)
        null <- apply(data[,-1], 1, null_count)
        non_zero <- apply(data[,-1], 1, non_zero_func)
        tbl <- tibble(
            Median = median,
            null = null,
            Non_zero = non_zero
        )
        tbl$Rank <- rank(tbl$Median)
        tbl<-tbl[tbl$Rank>870,]
        tbl$filter<-case_when(tbl$Non_zero > slider ~ 'genes passing the filter', tbl$Non_zero < slider ~ 'genes not passing filter')
        #tbl<-tbl[tbl$Rank>870,] #filtering extremely low ranks
        out <- ggplot2::ggplot(tbl, aes(x=Rank, y=null, color = filter)) +
            geom_point(alpha=0.5) +
            xlab("Rank(median)") +
            ylab("number of zeros") +
            theme_bw()+
            ggtitle('Number of zero counts vs Median')
        return (out)
    }
    
    #output tab 1
    observeEvent(input$submit,{
        data<-load_data()
        summary<-summary_table(data)
        output$summary<-renderTable(summary)
        output$table<-DT::renderDataTable({data})})

    observeEvent(input$plot_submit,{
        output$histogram<-renderPlot(plot_bases(load_data(),input$button))
    })
    
    #output tab 2.1
    observeEvent(input$summary_submit,{
        data2<-read_csv(input$file2$datapath)
        summary<-filter(data2, input$slider1, input$slider2)
        output$filter_counts<-renderTable(summary)
    })
    
    observeEvent(input$plot2,{
        data2<-read_csv(input$file2$datapath)
        varvsmed<-plot_variance_vs_median(data2,input$slider2)
        countvsmed<-plot_median_vs_null(data2,input$slider1)
        output$plot2<-renderPlot(varvsmed)
        output$plot3<-renderPlot(countvsmed)
    })
    

}
# Run the application 
shinyApp(ui = ui, server = server)
