## Author: Rojashree Jayakumar
## rshreej@bu.edu
## BU BF591
## Final Project
## R shiny app to view summary, counts normalization, differential expression and gene set enrichment analysis results (tables and plots)

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
library(gplots)
library('GSEABase')
library('fgsea')
library('biomaRt')
options(shiny.maxRequestSize=30*1024^2) 


ui <- fluidPage(
    
    titlePanel("Final Project"),
    
    h3("To use this application, upload the CSV files into the app from the Project directory of this app's git repository."),
    
    tabsetPanel(
        
        #Tab 1
        tabPanel("Samples",
                 h5("This tab displays the metadata characterstics. Only csv files will be 
                   accepted. Use the submit button to submit the metadata. The 'metadata.csv' in the 
                   git repository should be used."),
                 sidebarPanel(
                     fileInput(
                         "file","Load metadata", accept = ".csv"
                     ),
                     actionButton(
                         "submit", label = "Submit", width = "400px",
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                         style = "material-flat"
                     )
                 ),
                 
                 mainPanel(
                     
                     tabsetPanel(
                         
                         tabPanel(
                             "Summary",
                             h5("Table summarising the type and values of the metadata fields"),
                             tableOutput("summary")
                         ),
                         
                         tabPanel(
                             "Table",
                             h5("Metadata table with sortable columns and search bar to identify samples"),
                             DT::dataTableOutput("table")
                         ),
                         
                         tabPanel(
                             "Plot",
                             h5("Violin plots of the bases with respect to different traits"),
                             p("Choose a trait and click on plot"),
                             sidebarPanel(
                                 radioButtons(
                                     "button", "Choose the metadata trait",
                                     choices = c("lifestage", "timepoint", "treatment", "sex")
                                 ),
                                 actionButton(
                                     "plot_submit", label = "Plot", width = "230px",
                                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                                     style = "material-flat"
                                 )
                             ),
                             
                             mainPanel(
                                 plotOutput("histogram")
                             )
                         )
                     )
                 )
        ),
        
        #Tab 2
        
        tabPanel(
            "Counts",
            h5("This tab displays the characterstics of normalised counts of Drosophila from the RNA-seq experiment. Only csv files will be accepted. 
               The 'Normcounts.csv' in the git repository should be used."),
            p("Normalization was performed using DESeq2"),    
            
            sidebarPanel(
                p("The sliders are used to filter the uploaded data"),
                fileInput("file2", "Load normalised counts",accept = ".csv"),
                sliderInput(
                    inputId = "slider1", min = 0, max = 100,
                    label = "Select the threshold for non-zero genes", value = 50, step = 10
                ),
                sliderInput(
                    inputId = "slider2", min = 0, max = 100,
                    label = "Select the threshold for variance", value = 50, step = 10
                )
            ),
            
            mainPanel(
                
                tabsetPanel(
                    
                    tabPanel(
                        "Filter",
                        
                        p("Use the 'Display summary statistics' button to view the output"),
                        
                        sidebarPanel(
                            actionButton(
                                "summary_submit", label = "Display summary statistics", width = "230px",
                                style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                                style = "material-flat"
                            )
                        ),
                        
                        mainPanel(
                            h5("Table displaying the filtered data statistics of the normalised counts"),
                            tableOutput("filter_counts")
                        )
                        
                    ),
                    
                    tabPanel(
                        "Diagnostic_scatter_plots", label = "Diagnostic scatter plots",
                        p("Click on 'Plot' button to view scatter plots based on filtering"),
                        sidebarPanel(
                            actionButton(
                                "plot2", label = "Plot", width = "230px",
                                style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                                style = "material-flat"
                            )
                        ),
                        
                        mainPanel(
                            h5("Diagnostic scatter plots"),
                            plotOutput("plot2"),
                            plotOutput("plot3")
                        )
                        
                    ),
                    
                    tabPanel(
                        "Heatmap", label = "Heatmap",
                        p("Click on 'Heatmap' button to view heatmap based on filtering"),
                        
                        sidebarPanel(
                            actionButton("Heatmap_plot", label = "Heatmap", width = "230px",
                                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                                         style = "material-flat")
                        ),
                        
                        mainPanel(
                            h5("Heatmap of the filtered normalised counts"),
                            plotOutput("Heatmap_plot1")
                        )
                    ),
                    
                    tabPanel(
                        "PCA", label = "PCA",
                        p("Click on 'Plot' button to view heatmap based on filtering"),
                        sidebarPanel(
                            actionButton("PCA_button", label = "Plot", width = "230px",
                                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                                         style = "material-flat"),
                            radioButtons(
                                "button1", "Choose the column for the x-axis",
                                choices = c("PC1", "PC2", "PC3", "PC4","PC5","PC6"), selected = "log2FoldChange"
                            ),
                            radioButtons(
                                "button2", "Choose the column for the y-axis",
                                choices = c("PC1", "PC2", "PC3", "PC4","PC5","PC6"), selected = "padj"
                            )
                        ),
                        
                        mainPanel(
                            h5("Principal component analysis of the filtered normalised counts"),
                            plotOutput("PCA1")
                        )
                    )
                )
            )
        ),
        
        #Tab 3
        
        tabPanel(
            "Differential Expression",
            
            h5("This tab displays the characterstics of differential expression values of Drosophila from the RNA-seq experiment. 
            Only csv files will be accepted. Use the differential_expression.csv in the git repository."),
            p("Differential expression was performed using DESeq2"),
           
            sidebarPanel(
                fileInput(
                    "diff_expr","Load differential expression results", accept = ".csv"
                ),
                radioButtons(
                    "buttons1", "Choose the column for the x-axis",
                    choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj"), selected = "log2FoldChange"
                ),
                radioButtons(
                    "buttons2", "Choose the column for the y-axis",
                    choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj"), selected = "padj"
                ),
                colourInput(
                    "color",label = "Base point color", value = "#DD3916", showColour = c("both"), palette = c("square"), allowedCols = NULL, allowTransparent = FALSE, returnName = FALSE, closeOnClick = FALSE
                ),
                colourInput(
                    "colour",label = "Highlight point color", value = "#16E0C8", showColour = c("both"), palette = c("square"), allowedCols = NULL, allowTransparent = FALSE, returnName = FALSE, closeOnClick = FALSE
                ),
                sliderInput(
                    inputId = "padjusted", min = -300, max = 0,
                    label = "Select the magnitude of the p adjusted coloring:", value = -150, step = 30
                ),
                actionButton(
                    "diff_submit", label = "Plot", width = "400px",
                    style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                    style = "material-flat",
                    icon = icon("sliders")
                )
            ),
            
            mainPanel(
                
                tabsetPanel(
                    tabPanel(
                        "Results",
                        
                        h5("Differential expression results which can be sorted and genes can be searched."),
                        
                        DT::dataTableOutput("results")
                    ),
                    
                    tabPanel(
                        "Table",
                        
                        h5("Filtered differential expression results based on the slider value at the bottom."),
                        
                        tableOutput("table2")
                    ),
                    
                    tabPanel(
                        "Volcano", 
                        h5("Volcano plot colored by genes passing the filter by highlight point color"),
                        plotOutput("volcano")
                    ),
                )
            )
        ),
        
        #tab 4
        
        tabPanel(
            "FGSEA",
            h5("This tab displays the characterstics of gene set enrichment results of Drosophila from the RNA-seq experiment. Only csv files will be accepted. 
               Use the fgsea.csv in the git repository."),
            p("fgsea was performed on the raw counts and the gmt file for Drosophilla from http://ge-lab.org/gskb/ was used"),
            
            sidebarPanel(
                fileInput("fgsea", "Load FGSEA results",accept = ".csv"),
            ),
            
            mainPanel(
                
                tabsetPanel(
                    
                    tabPanel(
                        "Top NES",
                        h5("Top pathways as per the selected threshold"),
                        sidebarPanel(
                            p("Click on submit to view the top pathways"),
                            sliderInput(
                                inputId = "slider_padj", min = 0, max = 10,
                                label = "Select the threshold", value = 5, step = 1
                            ),
                            actionButton(
                                "submit_fgsea", label = "Submit", width = "230px",
                                style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                                style = "material-flat"
                            )
                        ),
                        
                        mainPanel(
                            h5("Bar chart displaying the top pathways"),
                            plotOutput("pathways")
                        )
                    ),
                    
                    tabPanel(
                        "Table",
                        
                        sidebarPanel(
                            p("Click on submit to view the filtered table based on positive/negative NES and p-adjusted threshold."),
                            sliderInput(
                                inputId = "slider_filter", min = 0.00000001, max = 1,
                                label = "Select the p-adjusted threshold", value = 0.5, step = 0.000001),
                            radioButtons(
                                "fgsea_button", "Choose the NES pathways",
                                choices = c("positive", "negative")),
                            actionButton(
                                "submit_fgseatable", label = "Submit", width = "230px",
                                style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                                style = "material-flat"
                            )
                        ),
                        
                        mainPanel(
                            h5("Sortable data table displaying the results" ),
                        
                            DT::dataTableOutput("fgsea_table"),
                            
                            downloadLink("downloadData", "Download")
                        )
                    ),
                    
                    tabPanel(
                        "Scatter Plot",
                        p("Click on plot to view the scatterplot"),
                        
                        sidebarPanel(
                            sliderInput(inputId = "slider10", min = 0.00000001, max = 1,
                                        label = "Select the p-adjusted threshold", value = 0.5, step = 0.000001),
                            actionButton("submit_scatter", label = "Plot", width = "230px",
                                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                                         style = "material-flat"
                            )
                        ),
                        
                        mainPanel(
                            h5("Scatter plot of NES and -log10 adjusted p-value"),    
                            
                            plotOutput("scatter_plot")
                        )
                    )
                )
            )
        )
    )
)


server <- function(input, output) {
    
    #tab 1 data
    
    load_data <- reactive({
        out<-read.csv(input$file$datapath, stringsAsFactors = TRUE)
        return(out)
    })
    
    #tab 1.1 summary
    
    #function to create a summary table showing the type of data and also statistics for numeric data
    #factors for factor types
    #@param: metadata.csv and @return tibble\
    
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
    
    # tab 1.3  violin plot
    
    #larvae have no gender, so empty strings are give value as larvae. 
    #@param x_name: UI radiobutton input and data is metadata.csv
    #@return violin plots
    
    plot_bases<-function(data,x_name){
        levels(data$sex) <- c("female","male", "larvae")
        data$sex[data$sex == ''] <- 'larvae'
        out<-ggplot(data) +
            geom_violin(aes(x=!!sym(x_name),y= Bases, fill= !!sym(x_name)))+ 
            theme_grey()+
            ggtitle("VIOLIN PLOT OF BASES")
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
    
    #get the sum of non zero values of each gene from input of normalised counts
    
    non_zero_func<- function (data){
        return (sum(data!=0))
    }
    
    #filtering the normalized counts. The percentile of variance is used
    
    filter<-function (norm_counts, slider1, slider2) {
        non_zero <- apply(norm_counts[,-1], 1, non_zero_func)
        variance <- apply(norm_counts[,-1], 1, var)
        non_zero_filtered<-norm_counts[non_zero>slider1 & rank(variance)>(17491*slider2)/100,]
        out <- tibble(
            `Number of genes` = 17491,
            `Number of samples` = 81,
            `Number of genes passing filter` = nrow(non_zero_filtered),
            `Number of genes not passing filter` = 17491 - nrow(non_zero_filtered),
            `Percentage of genes passing filter` = (nrow(non_zero_filtered)/17491)*100,
            `Percentage of genes not passing filter` = ((17491 - nrow(non_zero_filtered))/17491)*100,
        )
        return(out)
    }
    
    #tab 2.2
    
    #scatter plot of variance vs median of filtered normalized counts, colored by genes passing 
    #the filter
    #@param data is normalised counts and slider is the variance slider input from UI
    #@results is scatter plot
    
    plot_variance_vs_median <- function(data, slider) {
        median <- apply(data[,-1], 1, median)
        variances <- apply(data[,-1], 1, var)
        tbl <- tibble(
            Median = median,
            Variance = variances
        )
        tbl$Rank <- rank(tbl$Median) #ranking median so the plot is not distorted
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
    
    #function which return sum of null values from input of normalised data frame
    
    null_count<- function (data){
        return (sum(data==0))
    }
    
    #scatter plot of null values vs median of filtered normalized counts, colored by genes passing 
    #the filter
    #@param data is normalised counts and slider is the non-zero slider input from UI
    #@results is scatter plot
    
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
        tbl<-tbl[tbl$Rank>870,] #filtering extremely low ranks
        out <- ggplot2::ggplot(tbl, aes(x=Rank, y=null, color = filter)) +
            geom_point(alpha=0.5) +
            xlab("Rank(median)") +
            ylab("number of zeros") +
            theme_bw()+
            ggtitle('Number of zero counts vs Median')
        return (out)
    }
    
    #tab 2.3
    
    #heatmap function of the genes passing the filters
    #@params norm_counts is the data frame of normalised counts
    #@params slider1 and slider2 are variance and non-zero count filters from UI slider inputs
    #@return heatmap
    
    heatmap_func<-function(norm_counts, slider1, slider2){
        non_zero <- apply(norm_counts[,-1], 1, non_zero_func)
        variance <- apply(norm_counts[,-1], 1, var)
        non_zero_filtered<-norm_counts[non_zero>slider1 & rank(variance)>(17491*slider2)/100,]
        x<-non_zero_filtered[,-1]
        row.names(x)<-non_zero_filtered$gene
        plotting<-as.matrix(log(x))
        plotting[plotting == -Inf] <- 0
        return (heatmap.2(plotting, xlab='Patient Samples', ylab='Genes', 
                          main='Gene Expression Across Samples',trace='none', density.info = 'none',
                          key.xlab='Expression Level', scale='row', margins=c(5,5), cexRow = 0.5, cexCol = 0.2, key = 'TRUE'))
    }
    
    #tab 2.4
    
    #PCA 
    #@params norm_counts is the data frame of normalised counts
    #@params button1 and button2 are Prinicipal components from UI slider inputs
    #@return PCA plot
    
    PCA <- function (norm_counts,button1,button2){
        norm_counts<-as.data.frame(norm_counts)
        row.names(norm_counts)<-norm_counts$gene
        norm_counts<-t(norm_counts[-1])
        pca <- prcomp(norm_counts, center=TRUE, scale=TRUE)
        ve_percent<-(pca$sdev**2/sum(pca$sdev**2))*100
        ve1 <- unlist(strsplit(button1, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=TRUE)) #regexp to strip the input$button, so the numeric can be used for indexing
        ve1<-as.numeric(ve1[2])
        ve2 <- unlist(strsplit(button2, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=TRUE))
        ve2<-as.numeric(ve2[2])
        xlabel <- paste0(button1, ' % variance is ', round(ve_percent[ve1]))
        ylabel <- paste0(button2, ' % variance is ', round(ve_percent[ve2]))
        return(as_tibble(pca$x) %>%
                   ggplot(aes(x=!!sym(button1),y=!!sym(button2))) + 
                   xlab(xlabel)+
                   ylab(ylabel)+
                   geom_point())
    }
    
    #output tab 1
    #displays table
    
    observeEvent(input$submit,{
        data<-load_data()
        summary<-summary_table(data)
        output$summary<-renderTable(summary,striped = TRUE, bordered = TRUE,  
                                    hover = TRUE)
        output$table<-DT::renderDataTable({data})
    })
    
    #displays plot
    
    observeEvent(input$plot_submit,{
        output$histogram<-renderPlot(plot_bases(load_data(),input$button))
    })
    
    #output tab 2.1
    #displays table
    
    observeEvent(input$summary_submit,{
        data2<-read_csv(input$file2$datapath)
        summary<-filter(data2, input$slider1, input$slider2)
        output$filter_counts<-renderTable(summary,striped = TRUE, bordered = TRUE,  
                                          hover = TRUE)
    })
    
    #output tab 2.2
    
    #displays plot
    
    observeEvent(input$plot2,{
        data2<-read_csv(input$file2$datapath)
        varvsmed<-plot_variance_vs_median(data2,input$slider2)
        countvsmed<-plot_median_vs_null(data2,input$slider1)
        output$plot2<-renderPlot(varvsmed)
        output$plot3<-renderPlot(countvsmed)
    })
    
    #output tab 2.3
    
    #displays plot heatmap
    
    observeEvent(input$Heatmap_plot,{
        data3<-read_csv(input$file2$datapath)
        output$Heatmap_plot1<-renderPlot(heatmap_func(data3,input$slider1,input$slider2))
    })
    
    #output tab 2.4
    
    #displays plot PCA
    
    observeEvent(input$PCA_button,{
        data4<-read_csv(input$file2$datapath)
        output$PCA1<-renderPlot(PCA(data4,input$button1,input$button2))
    })
    
    #tab 3
    
    #reactive to load data of differential expression
    
    load_data2 <- reactive({
        if(is.null(input$diff_expr$datapath))     
            return(data.frame(matrix(ncol = 0, nrow = 0)))
        out<-read.csv(input$diff_expr$datapath, header = TRUE)
        names(out)[1] <- "Gene"
        return(out)
    })
    
    #volcano plot
    #@params: dataf: differential expression tibble
    #@param x_name and y_name radio buttons from UI
    #@param color1 and color2 color inputs from UI
    #@return plot
    
    volcano_plot <-
        function(dataf, x_name, y_name, slider, color1, color2) {
            if(is.null(dataf))     
                return(NULL)
            dataf<-na.omit(dataf)
            out<- ggplot(dataf,aes(x= !!sym(x_name),y= -log10(!!sym(y_name)),color= !!sym(y_name)<(10^slider)))+ 
                theme_bw()+
                theme(legend.position="bottom")+
                ggtitle('Volcano plot')+
                scale_color_manual(name= paste0("padj < 1 x 10^",toString(slider)), values=c(color1, color2))+
                geom_point()
            return(out)
        }
    
    #filtering the table by padj
    #@param dataf - differential expression reusults
    #@param slider - UI slider input
    #@return table
    
    draw_table <- function(dataf, slider) {
        dataf<-na.omit(dataf)
        out<-dataf[dataf$padj<(1*10^slider),]
        out$padj<-formatC(out$padj, digits = 7)
        out$pvalue<-formatC(out$pvalue, digits = 7)
        return(out)
    }
    
    #reactive funtion 
    
    load_table2 <- reactive({
        out<-read.csv(input$diff_expr$datapath, stringsAsFactors = TRUE)
        return(out)
    })
    
    #data table output
    
    observeEvent(input$diff_submit,{
        csv<-load_data2()
        data6<-load_table2()
        volcano <- volcano_plot(csv, input$buttons1, input$buttons2, input$padjusted, input$color, input$colour)
        table_2 <- draw_table(csv,input$padjusted)
        output$volcano <- renderPlot(volcano)
        output$table2 <-  renderTable(table_2, striped = TRUE, bordered = TRUE,  
                                      hover = TRUE)
        output$results <-DT::renderDataTable(as.matrix(data6))
    })
    
    #tab 4
    
    #barchart for the top NES pathways
    #@param fgsea_results - fgsea results
    #@param num_paths - int value from slider in UI
    #@retun plot
    
    top_pathways <- function(fgsea_results, num_paths){
        
        top_pos <- fgsea_results %>% slice_max(padj, n=num_paths) %>% pull(pathway)
        top_neg <- fgsea_results %>% slice_min(padj, n=num_paths) %>% pull(pathway)
        
        subset <-fgsea_results[fgsea_results$pathway %in% c(top_pos, top_neg),] 
        subset<- subset%>%
            mutate(pathway = factor(pathway)) %>%
            mutate(plot_name = str_replace_all(pathway, '_', ' '))
        
        plot <- subset %>% 
            mutate(plot_name = forcats::fct_reorder(factor(plot_name), NES)) %>%
            ggplot() +
            geom_bar(aes(x=plot_name, y=NES, fill = NES > 0), stat='identity', show.legend = FALSE) +
            scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
            theme_minimal(base_size = 8) +
            ggtitle('fgsea results for Drosophila_melanogaster') +
            ylab('Normalized Enrichment Score (NES)') +
            xlab('') +
            scale_x_discrete(labels = function(x) str_wrap(x, width = 80)) +
            coord_flip()
        return(plot)
    }
    
    #filtering based on padj and NES
    #@param data - fgsea results
    #@param parameter - radiobuttons in UI, positive or negative 
    #@param padjusted slider input of a numeric from UI
    #@return table
    
    table_fgsea<- function (data,parameter,padjusted){
        
        if (parameter == "positive") {
            return (dplyr::filter(data,padj>padjusted & NES>0))
        }
        return (dplyr::filter(data,padj>padjusted & NES<0)) 
    }
    
    #scatter plot for NES vs -lpg10(padj)
    #@param fgsea is fgsea results
    #@param threshold from UI slider input 
    
    scatter<-function(fgsea,threshold){
        ggplot(fgsea,mapping=aes(x=NES, y=-log10(padj), color = padj<threshold)) +
            geom_point()+
            ggtitle("Scatter plot of NES on x-axis and -log10 adjusted p-value on y-axis")
    }
    
    #output tab 4.1
    
    #bar chart of top pathways
    
    observeEvent(input$submit_fgsea,{
        data7<-read_csv(input$fgsea$datapath)
        hm<-top_pathways(data7,input$slider_padj)
        output$pathways<-renderPlot(hm)
    })
    
    #data table with sortable columns and gene search of FGSEA results that are filtered
    #download handler
    
    observeEvent(input$submit_fgseatable,{
        data8<-read_csv(input$fgsea$datapath)
        fgsea_table<-table_fgsea(data8,input$fgsea_button,input$slider_filter)
        output$fgsea_table<-DT::renderDataTable({fgsea_table})
        output$downloadData <- downloadHandler(
        filename = paste("FGSEA", ".csv", sep = ""),
        content = function(file) {
            write.csv(fgsea_table, file)
        }
    )
    })


    #scatter plot of NES vs -log10 padj
    
    observeEvent(input$submit_scatter,{
        data9<-read_csv(input$fgsea$datapath)
        scatterplot<-scatter(data9,input$slider10)
        output$scatter_plot<-renderPlot(scatterplot)
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)