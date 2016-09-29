################################################################################
# shiny application to visualize EDEN output .tar files 
# by Philipp C. Münch (pmu15@helmholtz-hzi.de)
################################################################################

# import functions
library(shiny)
library(ggplot2)
library(zoo)
library(gridExtra)
source("functions.R")
source("version.R")

# no warnings
options(warn=0)

# set path variable
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

# initialize dataset
dataset <<- readData(paste(folder.path,
                           list.dirs(path = folder.path,
                                     full.names = FALSE,
                                     recursive = FALSE)[1], 
                           sep="/")) 

# define header
headerPanel_2 <- function(title, h, windowTitle=title){    
  tagList(
    tags$head(tags$title(windowTitle)),
    h(title)
  )
}

shinyUI(fluidPage(headerPanel_2(
      HTML(paste('eden selection viewer', eden.version)), h3, "test"),
      fluidRow(column(3,
                      wellPanel(
                        
# ----------- table panel
conditionalPanel(
  condition="input.tsp=='overview' || input.tsp=='annotation' || input.tsp=='alignment' ",
  helpText("Please select the analysis you want to import."),
  selectInput("dataset", "Select run:", 
              choices= list.dirs(path = "data", 
                                 full.names = FALSE, recursive = FALSE), 
              selected=list.dirs(path = "csv", full.names = FALSE, 
                                 recursive = FALSE)[1], multiple=F, width="100%"),
  uiOutput('filters_UI'),
  selectInput('dofiltering', 
              label='Choose a way to filter matrix', 
              choices=c("pvalue", "ratio"), 
              selected="no filtering"),
  uiOutput("dependentselection"),
  actionButton('resetSelection',label = "Click to reset row selection")),



# ----------- histogram panel
conditionalPanel(
  condition="input.tsp=='ts'",
  sliderInput('binSize', 'Number of bins', min=10, max=500,
              value=min(10, 500), step=10, round=0),
  checkboxInput('facet', 'Facet by sample'),
  helpText("Specific gene families will be shown here if they are selected selected on the Data Table panel")),

conditionalPanel(
  condition="input.tsp=='ts'",
  downloadButton("dlCurPlot", "Download Graphic")), 
       
conditionalPanel(
  condition="input.tsp=='map'",
  tags$button(
    id = 'close',
    type = "button",
    class = "btn action-button",
    onclick = "setTimeout(function(){window.close();},500);", 
    "Close Application")), 
                      
# ----------- boxplot panel
conditionalPanel(
  condition="input.tsp=='box'",
  selectInput("oderchoice", label = "Order by", 
              choices = list("Dataset" = "default", "Mean ratio" = "mean", "p-value" = "pvalue"),
              selected = "default"),
  checkboxInput('highlightbox', 'Highlight mean of selected elements')),

conditionalPanel(
  condition="input.tsp=='box'",
  downloadButton("dlCurBoxPlot", "Download Graphic")),
                      

# ----------- alignmentplot panel
conditionalPanel(
  condition="input.tsp=='annot'",
  checkboxInput('navalues', 'remove NA', value=TRUE),
  checkboxInput('showmean', 'plot mean value', value=TRUE),
  checkboxInput('bysamplefacet', 'facet by sample'),
  checkboxInput('bysamplecolor', 'color by sample'),
  checkboxInput('showmeanselected', 'plot mean of selected families'),
  selectInput("sortannotation", label = "Order by", 
              choices = list("ratio" = "ratio", "p-value" = "pvalue"), 
              selected = "ratio"),
  downloadButton("dlCurAnnotationplot", "Download Graphic"))),
             

 wellPanel(
  conditionalPanel(condition="input.tsp=='overview'",
                   downloadButton("dlTable", "Download Table", class="btn-block btn-primary")
                   ),
  conditionalPanel(condition="input.tsp=='annotation'",
                   downloadButton("dlAnnotationPlot", "Download Figure", class="btn-block btn-primary")
  ),
  # ----------- alignmentplot panel
  conditionalPanel(
    condition="input.tsp=='alignment'",
   
    checkboxInput('points', 'show points', value=TRUE),
    uiOutput('colorpoints'),
    downloadButton("dlCurSequenceplot", "Download Graphic"))
  )),


column(9,tabsetPanel(

tabPanel("Overview",htmlOutput("overview_hint"),
         div(DT::dataTableOutput("table"),style = "font-size:80%"), 
         htmlOutput("overview_table"), 
         htmlOutput("overview_fisher"),
         verbatimTextOutput("summary"),  value="overview"),

tabPanel("Annotation",htmlOutput("annotation_hint", inline = FALSE), 
         plotOutput("annotationplotglobal", width="100%", height="auto"),
         htmlOutput("annotation_figure"),  value="annotation"),

tabPanel("Alignment Plot",
         htmlOutput("alignment_hint"), 
         plotOutput("alignmentplot", width="100%", height="auto"), 
         htmlOutput("alignment_figure"),
         value="alignment"),

tabPanel("Categories",  
         div(DT::dataTableOutput("table_annotaion"), style = "font-size:80%"), 
         plotOutput("annotationplot", width="100%", height="auto"), value="annot"),

tabPanel("Histogram",h4(""), 
         plotOutput("plot1", width="100%", height="auto"), value="ts"),
                      
tabPanel("Boxplot",h4(""),  
         div(DT::dataTableOutput("table_sample"), style = "font-size:80%"), 
         plotOutput("plot4", width="100%", height="auto"), value="box"),id="tsp")
 
))
                    

))