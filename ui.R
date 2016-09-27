library(shiny)
library(ggplot2)
library(zoo)
#options(warn=0)

source("functions.R")

# check if this applications runs inside the docker container
if(file.exists("/home/eden/eden.sh")){
  # we are inside the docker container
  packrat::on()
  tar.path <<- "/home/eden/data/tar"
  annotation.path <<- "/home/eden/data/annotation"
  tmp.path <<- "/srv/shiny-server/eden-visualizer/tmp/"
  folder.path <<- "data"
  dir.create(folder.path)
} else {
  # we are online hosted
  packrat::on()
  folder.path <<- "data"
  annotation.path <<- "annotation"
  tar.path <<- "examples"
  tmp.path <<- "tmp"
  dir.create("tmp")
  dir.create(folder.path)
}

# extract tar file on startup
unTar(tar.path, folder.path)

dataset <<- readData(paste(folder.path, list.dirs(path = folder.path, full.names = FALSE, recursive = FALSE)[1], sep="/")) 



headerPanel_2 <- function(title, h, windowTitle=title) {    
  tagList(
    tags$head(tags$title(windowTitle)),
    h(title)
  )
}

shinyUI(
  fluidPage(
    headerPanel_2(
      HTML('eden selection viewer v. 0.0.2'
      ), h3, "test"
    ),
    fluidRow(column(4,
                    wellPanel(
                      
                      # ----------- table panel
                      conditionalPanel( 
                        condition="input.tsp=='map'",
                        helpText("If you want to use the results of a previous eden run, select the analysis here"), 
                     
                        selectInput("dataset", "Import a previous eden run", choices= list.dirs(path = "data", full.names = FALSE, recursive = FALSE), selected=list.dirs(path = "csv", full.names = FALSE, recursive = FALSE)[1], multiple=F, width="100%"),
                         helpText("The selection of the dataset and p-value threshold will filter the data for all Panels. Try it out and select also the sample obese_m and obese_f"),
                        
                        uiOutput('filters_UI'),
                        
                        
                        
                        #selectInput("samples", "Choose one or more datasets:", choices=levels(factor(dataset$sample)), selected=c(levels(factor(dataset$sample))[1]), multiple=T, width="100%"),
                        
                        textInput("pval", label = "Filter by p-value", value = "1"),
                        conditionalPanel(
                          condition="input.tsp=='map'",
                          actionButton(  'resetSelection',label = "Click to reset row selection"),
                          helpText("You can select rows on the table and a p-value will be calculated to check if the selected families have a significant higher ratio than the background."),
                          downloadButton("dlTable", "Download Table", class="btn-block btn-primary")
                        )
                      ),
                      
                      # ----------- histogram panel
                      conditionalPanel( 
                        condition="input.tsp=='ts'",
                        sliderInput('binSize', 'Number of bins', min=10, max=500,
                                    value=min(10, 500), step=10, round=0),
                     
                        checkboxInput('facet', 'Facet by sample'),
                        helpText("Specific gene families will be shown here if they are selected selected on the Data Table panel")
                      ),
                      conditionalPanel(
                        condition="input.tsp=='ts'",
                        downloadButton("dlCurPlot", "Download Graphic")
                      ), 
       
                      conditionalPanel(
                        condition="input.tsp=='map'",helpText("You have to close this application before you can stop eden"),
                        tags$button(
                          id = 'close',
                          type = "button",
                          class = "btn action-button",
                          onclick = "setTimeout(function(){window.close();},500);",  # close browser
                          "Close eden Visualizer"
                        )
                      ), 
                      
                      # ----------- boxplot panel
                      conditionalPanel( 
                        condition="input.tsp=='box'",
                        selectInput("oderchoice", label = "Order by", 
                                    choices = list("Dataset" = "default", "Mean ratio" = "mean", "p-value" = "pvalue"), 
                                    selected = "default"),
                        checkboxInput('highlightbox', 'Highlight mean of selected elements')
                      ),
                      conditionalPanel(
                        condition="input.tsp=='box'",
                        downloadButton("dlCurBoxPlot", "Download Graphic")
                      ),
                      
                      # ----------- alignmentplot panel
                      conditionalPanel( 
                        condition="input.tsp=='ap'",
                        checkboxInput('gap', 'color by gap proportion'),
                        checkboxInput('epitopes', 'highlight epitopes')
                      )
                      
                    )),
             
             column(8,
                    tabsetPanel(
                      
                      tabPanel("Data Table",h4(""),div(DT::dataTableOutput("table"), style = "font-size:80%"), verbatimTextOutput("summary"), value="map"),
                      tabPanel("Alignment Plot", plotOutput("plot3", width="100%", height="auto"), value="ap"),
                      tabPanel("Histogram",h4(""), plotOutput("plot1", width="100%", height="auto"), value="ts"),
                      tabPanel("Boxplot",h4(""), plotOutput("plot4", width="100%", height="auto"), value="box"),
                id="tsp")
                    
                    
             )
             
    )
  ))