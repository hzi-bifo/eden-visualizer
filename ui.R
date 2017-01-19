
# import functions
library(shiny)
library(ggplot2)
library(zoo)
library(gridExtra)
library(shinyjs)
library(reshape2)
source("functions.R")
source("version.R")

# no warnings
options(warn = -1)
 
no_tar <<- TRUE
new_files <<- TRUE

# set path variable
if (file.exists("/home/eden/eden.sh")) {
  # we are inside the docker container
  # packrat::on()
  csv.path <<- "/home/eden/data/csv" # folder where processed .tar files are located in csv format (one file for a sample)
  tar.path <<- "/home/eden/data/tar" # folder where .tar files are located (eden output)
  raw.path <<- "/home/eden/data/raw" # folder where unpacked .tar files are located
  annotation.path <<- "/home/eden/tigr_data" # folder werhere TIGR_ROLE_NAMES and TIGRFAMS_ROLE_LINK are located
  dir.create(csv.path)
  dir.create(raw.path)
  dir.create(tar.path)
} else {
  # we are online hosted
  csv.path <<- "csv"
  tar.path <<- "tar"
  raw.path <<- "raw"
  dir.create(csv.path)
  dir.create(raw.path)
  annotation.path <<- "annotation"
}

# check if .csv files are located
# TODO(pmuench): if no_csv == TRUE dont show tabs
dirs <- list.dirs(csv.path,recursive=FALSE)
files.list <- list.files(dirs, ".*\\.csv", recursive=TRUE, full.names=TRUE)
if (length(files.list) > 0){
  no_csv <<- FALSE
} else {
  no_csv <<- TRUE
}

# check if .tar files are located
files.list <- list.files(tar.path, ".*\\.tar", recursive=TRUE, full.names=TRUE)
if (length(files.list) > 0){
  no_tar <<- FALSE
} else {
  no_tar <<- TRUE
}

# check if there are new tar files
tar.list <- basename(list.files(tar.path, ".*\\.tar", recursive=TRUE, full.names=TRUE))
csv.list <- basename(list.dirs(csv.path,recursive=FALSE))
if (any(is.na(match(tar.list, csv.list)))) { # true if the match is not complete (one tar file name dont match to csv file name)
  new_files <<- TRUE
} else {
  new_files <<- FALSE
}

# define header
headerPanel_2 <- function(title, h, windowTitle = title) {
  tagList(tags$head(tags$title(windowTitle)),
          tags$h1(a(href="www.someURLlogoLinksto.com")))
}

shinyUI(fluidPage(theme = "edentheme.css",
  
                  
       
                             
  
  headerPanel_2(HTML(paste(
    'eden visualizer ', eden.version
  )), h3, "eden visualizer"),
  
  
  #navbarPage("EDEN",
 
          
  
  
  
  fluidRow(
    useShinyjs(),
    column(
      4,
      wellPanel(
        conditionalPanel(condition = "input.tsp=='start' || input.tsp=='log'",
                         uiOutput("startdown_UI")),
        uiOutput("main_ui"),
        conditionalPanel(
          condition = "input.tsp=='map'",
          tags$button(
            id = 'close',
            type = "button",
            class = "btn action-button",
            onclick = "setTimeout(function(){window.close();},500);",
            "Close Application"
          )
        )
      ),
      
      wellPanel(
        conditionalPanel(condition = "input.tsp=='start' || input.tsp=='log'",
                         uiOutput("start_UI")),
        
        conditionalPanel(
          condition = "input.tsp=='start'",
          htmlOutput("reloadstatus"), 
          actionButton('reloadButton', label = "Reload/Import files")
        ),
        
        conditionalPanel(
          condition = "input.tsp=='overview'",
          
          ##s
          downloadButton("dlTable", "Download filtered table")
        ),
        
        conditionalPanel(
          condition = "input.tsp=='annotation'",
          downloadButton("dlAnnotationPlot", "Download barplot")
        ),
        
        conditionalPanel(
          condition = "input.tsp=='alignment'",
          checkboxInput('points', 'show points', value =
                          TRUE),
          uiOutput('colorpoints'),
          downloadButton("dlCurSequenceplot", "Download sequenceplot")
        ),
        
        conditionalPanel(
          condition = "input.tsp=='histogram'",
          sliderInput(
            'binSize',
            'Number of bins',
            min = 10,
            max = 500,
            value = min(10, 500),
            step = 10,
            round = 0
          ),
          checkboxInput('facet', 'Facet by sample'),
          downloadButton("dlCurPlot", "Download histogram")
        ),
     #   conditionalPanel(condition = "input.tsp=='start'", uiOutput("reload_ui")), 
        
        conditionalPanel(
          condition = "input.tsp=='categories'",
          checkboxInput('navalues', 'remove NA', value =
                          TRUE),
          checkboxInput('showmean', 'plot mean value', value =
                          TRUE),
          checkboxInput('bysamplefacet', 'facet by sample'),
          checkboxInput('bysamplecolor', 'color by sample'),
          checkboxInput('showmeanselected', 'plot mean of selected families'),
          selectInput(
            "sortannotation",
            label = "Order by",
            choices = list("ratio" = "ratio", "p-value (not implemented)" = "pvalue"),
            selected = "ratio"
          ),
          downloadButton("dlCurAnnotationplot", "Download boxplot")
        ),
        conditionalPanel(
          condition = "input.tsp=='box'",
          selectInput(
            "oderchoice",
            label = "Order by",
            choices = list("Dataset name" = "default", "Mean ratio" = "mean"),
            selected = "default"
          ),
          checkboxInput('highlightbox', 'Highlight mean of selected elements'),
          downloadButton("dlCurBoxPlot", "Download boxplot")
        )
      )
    ),
    
    column(
      8,
      tabsetPanel(
        tabPanel(
          "Start",
          textOutput("reloadmsg"),
          htmlOutput("welcome"),
          # htmlOutput("tar_check"), 
          # htmlOutput("newtar"), 
          # htmlOutput("csv_check"), 
          # htmlOutput("selected_dataset"),
          # htmlOutput("selected_samples"),
           value = "start"
        ),
        
        tabPanel(
          "Overview",
          htmlOutput("overview_hint"),
          div(DT::dataTableOutput("table"), style = "font-size:80%"),
          htmlOutput("overview_table"),
          htmlOutput("summary2"),
          value = "overview"
        ),
        
        tabPanel(
          "Annotation",
          htmlOutput("annotation_hint", inline = FALSE),
          plotOutput("annotationplotglobal", width = "100%", height =
                       "auto"),
          htmlOutput("annotation_figure"),
          value = "annotation"
        ),
        
        tabPanel(
          "Alignment Plot",
          htmlOutput("alignment_hint"),
          plotOutput("alignmentplot", width = "100%", height =
                       "auto"),
          htmlOutput("alignment_figure"),
          value = "alignment"
        ),
        
        tabPanel(
          "Categories",
          plotOutput("annotationplot", width = "100%", height =
                       "auto"),
          div(DT::dataTableOutput("table_annotaion"), style = "font-size:80%"),
          value = "categories"
        ),
        
        tabPanel(
          "Histogram",
          h4(""),
          plotOutput("plot1", width = "100%", height = "auto"),
          value = "histogram"
        ),
        
        tabPanel(
          "Boxplot",
          h4(""),
          htmlOutput("boxplot_hint"),
          plotOutput("plot4", width = "100%", height = "auto"),
          
          div(DT::dataTableOutput("table_sample"), style = "font-size:80%"),
          
          value = "box"
        ),
        id = "tsp"
      )
    )
) # fluid row

))
