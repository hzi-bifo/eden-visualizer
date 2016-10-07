################################################################################
# shiny application to visualize EDEN output .tar files 
# by Philipp C. Münch (pmu15@helmholtz-hzi.de)
################################################################################

# import functions
library(shiny)
library(ggplot2)
library(zoo)
library(gridExtra)
library(shinyjs)
source("functions.R")
source("version.R")

# no warnings
options(warn=-1)

# set path variable
if(file.exists("/home/eden/eden.sh")){
  # we are inside the docker container
 # packrat::on()
  tar.path <<- "/home/eden/data/tar"
  annotation.path <<- "/home/eden/data/annotation"
  fasta.path <<- "/home/eden/data/fasta"
  tmp.path <<- "/srv/shiny-server/eden-visualizer/tmp/"
  faa.path <<- "/home/eden/data/faa"
  ffn.path <<- "/home/eden/data/ffn"
  folder.path <<- "data"
  log.path <<- "/home/eden/log.txt"
  dir.create(folder.path)
} else {
  # we are online hosted
  folder.path <<- "data"
  annotation.path <<- "annotation"
  tar.path <<- "examples"
  tmp.path <<- "tmp"
  fasta.path <<- "fasta"
  log.path <<- "log.txt"
  faa.path <<- "faa"
  ffn.path <<- "ffn"
  dir.create("tmp")
  dir.create(folder.path)
}

# extract tar file on startup
try(unTar(tar.path, folder.path))

# initialize dataset
#dataset <<- readData(paste(folder.path,
#                           list.dirs(path = folder.path,
#                                     full.names = FALSE,
#                                     recursive = FALSE)[1], 
#                           sep="/")) 

# define header
headerPanel_2 <- function(title, h, windowTitle=title){    
  tagList(
    tags$head(tags$title(windowTitle)),
    h(title)
  )
}

shinyUI(fluidPage(headerPanel_2(
      HTML(paste('eden ', eden.version)), h3, "eden"),
      fluidRow( useShinyjs(), column(4,
                      wellPanel(


# ----------- table panel
conditionalPanel(condition="input.tsp=='start'",
                 helpText("Choose if you want to start a new eden run or analyse a completed run."),
                 selectInput("runtype", label = "Select mode", 
                             choices = list("Start new run" = "newstart", "Select previous run" = "previousstart"),
                             selected = "newstart")
                 ),

conditionalPanel(condition="input.tsp=='start'",
                 uiOutput("startdown_UI")),

conditionalPanel(
  condition="input.tsp=='overview' ||
input.tsp=='annotation' || 
  input.tsp=='alignment' || 
  input.tsp=='histogram' || 
  input.tsp=='box' ||
  input.tsp=='categories' ",
  helpText("Please select the analysis you want to import."),

  
  uiOutput('filters_UI'),
  selectInput('dofiltering', 
              label='Choose a way to filter matrix', 
              choices=c("pvalue", "ratio"), 
              selected="no filtering"),
  uiOutput("dependentselection"),
  actionButton('resetSelection',label = "Click to reset row selection")),

conditionalPanel(
  condition="input.tsp=='map'",
  tags$button(
    id = 'close',
    type = "button",
    class = "btn action-button",
    onclick = "setTimeout(function(){window.close();},500);", 
    "Close Application"))
),

wellPanel(
  conditionalPanel(condition="input.tsp=='start'",
                   uiOutput("start_UI")
                  
                   ),
  
  conditionalPanel(condition="input.tsp=='overview'",
                   downloadButton("dlTable", "Download Table", class="btn-block btn-primary")),
  
  conditionalPanel(condition="input.tsp=='annotation'",
                   downloadButton("dlAnnotationPlot", "Download Figure", class="btn-block btn-primary")),
  
  conditionalPanel(
    condition="input.tsp=='alignment'",
    checkboxInput('points', 'show points', value=TRUE),
    uiOutput('colorpoints'),
    downloadButton("dlCurSequenceplot", "Download Graphic",  class="btn-block btn-primary")),
  
  conditionalPanel(
    condition="input.tsp=='histogram'",
    sliderInput('binSize', 'Number of bins', min=10, max=500,
                value=min(10, 500), step=10, round=0),
    checkboxInput('facet', 'Facet by sample'),
    helpText("Specific gene families will be shown here if they are selected selected on the Data Table panel"),
    downloadButton("dlCurPlot", "Download Graphic", class="btn-block btn-primary")),
  
  conditionalPanel(
    condition="input.tsp=='categories'",
    checkboxInput('navalues', 'remove NA', value=TRUE),
    checkboxInput('showmean', 'plot mean value', value=TRUE),
    checkboxInput('bysamplefacet', 'facet by sample'),
    checkboxInput('bysamplecolor', 'color by sample'),
    checkboxInput('showmeanselected', 'plot mean of selected families'),
    selectInput("sortannotation", label = "Order by", 
                choices = list("ratio" = "ratio", "p-value" = "pvalue"), 
                selected = "ratio"),
    downloadButton("dlCurAnnotationplot", "Download Graphic", class="btn-block btn-primary")),
  
  conditionalPanel(
    condition="input.tsp=='box'",
    selectInput("oderchoice", label = "Order by", 
                choices = list("Dataset" = "default", "Mean ratio" = "mean", "p-value" = "pvalue"),
                selected = "default"),
    checkboxInput('highlightbox', 'Highlight mean of selected elements'),
    downloadButton("dlCurBoxPlot", "Download Graphic"))
  )),


column(8,tabsetPanel(
tabPanel("eden log",
         htmlOutput("log_hint"),
         tableOutput("log"),
#         tableOutput("filetable_faa")
 #        tableOutput("filetable_ffn"), 
#         htmlOutput("nTextupload"), 
#         htmlOutput("nTextcheck"), 
#         htmlOutput("nTexteden"),
         value="start"), 
  
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
         plotOutput("annotationplot", width="100%", height="auto"), 
         div(DT::dataTableOutput("table_annotaion"), style = "font-size:80%"), 
         value="categories"),

tabPanel("Histogram",h4(""), 
         plotOutput("plot1", width="100%", height="auto"), value="histogram"),
                      
tabPanel("Boxplot",h4(""),
         htmlOutput("boxplot_hint"),
         plotOutput("plot4", width="100%", height="auto"),
         htmlOutput("boxplot_figure"),
         div(DT::dataTableOutput("table_sample"), style = "font-size:80%"),
         htmlOutput("boxplot_table"),
         value="box"),id="tsp")
 
))
))