################################################################################
# shiny application to visualize EDEN output .tar files 
# by Philipp C. Münch (pmu15@helmholtz-hzi.de)
################################################################################

# set max file size for file upload
options(shiny.maxRequestSize = 50*1024^2)

shinyServer(function(input, output, session) {

  # check if dataset is there and only then show the tab panels
  observe({
    toggle(condition = (input$runtype!="newstart"), selector = "#tsp li a[data-value=annotation]")
    toggle(condition = (input$runtype!="newstart"), selector = "#tsp li a[data-value=overview]")
    toggle(condition = (input$runtype!="newstart"), selector = "#tsp li a[data-value=alignment]")
    toggle(condition = (input$runtype!="newstart"), selector = "#tsp li a[data-value=categories]")
    toggle(condition = (input$runtype!="newstart"), selector = "#tsp li a[data-value=histogram]")
    toggle(condition = (input$runtype!="newstart"), selector = "#tsp li a[data-value=box]")
  })
  
  #################
  # Reactive
  #################
  dataset <- reactive({
      if(input$dataset[1]!=""){
        readData(paste(folder.path, input$dataset, sep="/"))
      } else {
      #  dataset <- NULL
      }
     
  })
  
  #################
  # render table functions
  #################
  # check if annoation is provided and only then show the annotation panels
#  observe({
#    toggle(condition = input$cond, selector = "#tsp li a[data-value=annotation]")
#  })
 
  output$table_filtered <- DT::renderDataTable(
    DT::datatable(dataset, options = list(paging = 25))
  )
  
  output$table_annotaion <- DT::renderDataTable(DT::datatable({
    require(pander)
    data <- readData(paste(folder.path, input$dataset, sep="/"))
    data <- data[which(data$fdr <= input$pval),]
    if(length(input$samples) > 1){
      subset <- NULL
      data_pool <- NULL
      for(i in 1:length(input$samples)){
        subset <- data[which(data$sample == input$samples[i]),]
        data_pool <- rbind(subset,data_pool)
      }
      data <- data_pool
    } else {
      data <- data[which(data$sample == input$samples),]
    }
    
      df <- NULL
      df <- data.frame(term=unique(data$term), pval=rep(-1,length(unique(data$term))), elements=rep(0,length(unique(data$term))))
      df <- df[which(!is.na(df$term)),]
      i <- 1
      for(term in df$term){
        data.term <- data[which(data$term == term),]
        data.nonterm <- data[which(data$term != term),]
        test.mat <-
          matrix(c(sum(data.nonterm$sum_pN), sum(data.term$sum_pN), sum(data.term$sum_pS), sum(data.nonterm$sum_pS)),
                 nrow = 2,
                 dimnames =
                   list(c("background", "selected"),
                        c("dN", "dS")))
        df[i,]$pval <- fisher.test(test.mat, alternative = "greater")$p.value
        df[i,]$elements <- df[i,]$elements  + nrow(data.term)
        i <- i + 1 
      }
      df$fdr <- p.adjust(df$pval, method="fdr")
      df
    }))
  
  output$table_sample <- DT::renderDataTable(DT::datatable({
    require(pander)
    data <- readData(paste(folder.path, input$dataset, sep="/"))
    data <- data[which(data$fdr <= input$pval),]
    if(length(input$samples) > 1){
      subset <- NULL
      data_pool <- NULL
      for(i in 1:length(input$samples)){
        subset <- data[which(data$sample == input$samples[i]),]
        data_pool <- rbind(subset,data_pool)
      }
      data <- data_pool
    } else {
      data <- data[which(data$sample == input$samples),]
    }
    
    df <- NULL
    df <- data.frame(sample=unique(data$sample), pval=rep(-1,length(unique(data$sample))))
    df <- df[which(!is.na(df$sample)),]
    i <- 1
    for(sample in df$sample){
      data.sample <- data[which(data$sample == sample),]
      data.nonsample <- data[which(data$sample != sample),]
      test.mat <-
        matrix(c(sum(data.sample$sum_pN), sum(data.sample$sum_pS), sum(data.nonsample$sum_pN), sum(data.nonsample$sum_pS)),
               nrow = 2,
               dimnames =
                 list(c("dN", "dS"),
                      c("selected", "background")))
      df[i,]$pval <- fisher.test(test.mat, alternative = "greater")$p.value
      i <- i + 1 
    }
    df$fdr <- p.adjust(df$pval, method="fdr")
    df
}))
  
  output$table <- DT::renderDataTable(
    DT::datatable(dataset, options = list(pageLength = 25))
  )
  
  # Filter data based on selections
  output$table <- DT::renderDataTable(DT::datatable({
    require(pander)
    data <- readData(paste(folder.path, input$dataset, sep="/"))
    data <- data[which(data$fdr <= input$pval),]
    if(length(input$samples) > 1){
      subset <- NULL
      data_pool <- NULL
      for(i in 1:length(input$samples)){
        subset <- data[which(data$sample == input$samples[i]),]
        data_pool <- rbind(subset,data_pool)
      }
      data <- data_pool
    } else {
      data <- data[which(data$sample == input$samples),]
    }
    data$stars <- add.significance.stars(data$fdr)
    data$sum_pN <- NULL
    data$sum_pS <- NULL
    data$role <- NULL
    data$pvalue <- NULL
    data$ratio <- round(data$ratio, digits=3)
    # names(data)[1] <- c("family")
    #  names(data)[2] <- c("dN/dS ratio")
    #  names(data)[6] <- c("annotation")
    data$fdr <- round(data$fdr, digits=5)
    num.name <<- nrow(data)
    num.meanratio <<- round(mean(data$ratio, na.rm=T), digits=2)
    num.sd <<- round(sd(data$ratio, na.rm=T), digits=3)
    downloadObj <<- data # generate downloadable table
    input$resetSelection 
    data
  }))
  
  #################
  # Render UI
  #################
    
  # render again if input$player_name changes
  output$filters_UI <- renderUI({
    dataset <- readData(paste(folder.path, input$dataset, sep="/"))
#    selectInput("samples", "Choose one or more samples:", 
#    choices=levels(factor(dataset$sample)), 
#    selected=c(levels(factor(dataset()$sample))[1:length(unique(data$sample))]), 
#    multiple=T, width="100%")
    selectInput("samples", "Choose one or more samples:", 
                choices=levels(factor(dataset$sample)), 
                selected=c(levels(factor(dataset()$sample))[1]), 
                multiple=T, width="100%")
  })
  
  # render again if input$dofiltering changes
  output$dependentselection <- renderUI({
    if (input$dofiltering == "pvalue"){
      sliderInput("pval", label = "p-value threshold", min = .001, 
                  max = 1, value = 1)
    } else {
      sliderInput("ratio", label = "select ratio range to display", 
                  min = round(min(data$ratio), digits = 2), 
                  max = round(max(data$ratio), digits = 2), 
                  value = c(round(min(data$ratio), digits = 2), 
                            round(max(data$ratio), digits = 2)))
    }
  })
  
  # render again if input$player_name changes
  output$start_UI <- renderUI({
  if(input$runtype == "newstart"){
    conditionalPanel(condition="input.tsp=='start'",
    helpText("To start a eden run please upload .fasta files"),
    fileInput('files_faa', 'Choose .faa file to upload',
              accept = c(
                '.faa'
              ), multiple=TRUE),
    fileInput('files_ffn', 'Choose .ffn files to upload',
              accept = c(
                '.ffn'
              ), multiple=TRUE),
   actionButton('uploadButton',label = "Add files"),
   actionButton('checkButton',label = "Check files"),
   actionButton('goButton',label = "Start analysis")
  )}
  })
  
  output$startdown_UI <- renderUI({
    if(input$runtype != "newstart"){
      conditionalPanel(condition="input.tsp=='start'",#
                       helpText("You can also select an previous eden run from the dropdown list"),
                       selectInput("dataset", "Select run:", 
                                   choices= list.dirs(path = "data", 
                                                      full.names = FALSE, recursive = FALSE), 
                                   selected=list.dirs(path = "csv", full.names = FALSE, 
                                                      recursive = FALSE)[1], multiple=F, width="100%"))
      }
  })
  
  
  
  # render again if input$dofiltering changes
  output$colorpoints <- renderUI({
    if (input$points){
      checkboxInput('gap', 'color by gap proportion', value=TRUE)
    } 
  })
  
 
  #################
  # ggplot functions
  #################
  
  # generates a histogram
  doPlotHistogram <- function(dataset, input){
    p <- ggplot(dataset, aes(ratio, fill = sample)) +
      geom_histogram(bins = input$binSize) + theme_bw()
    if (input$facet)
      p <- p + facet_grid(sample ~ .)
    
    if (length(input$table_rows_selected)) {
      # get the ratio of selected rows
      mark.ratio <- dataset[input$table_rows_selected,]$ratio
      mark.name <- dataset[input$table_rows_selected,]$name
      p <- p + geom_vline(xintercept = mark.ratio)
    }
    return(p)
  }
  
  # this functions calls the create_msa_plot() function multiple times based
  # on selected protein families
  doAlignmentPlot <- function(data, input){
    require(ggplot2)
    require(grid)
    require(gridExtra)
    fam_ids <- data$name
    dnds <- paste(tmp.path,"/",input$dataset[1],"/",input$samples,"/dnds/", fam_ids,".txt.DnDsRatio.txt", sep="")
    gap <- paste(tmp.path,"/", input$dataset[1],"/",input$samples,"/gap/", fam_ids,".gap.txt", sep="")
    if(input$points){
      if(input$gap){
        p <- list()
        for(i in 1:length(dnds)){
          p[[i]] <- create_msa_plot(dnds_path = dnds[i], 
                                    gap_path = gap[i],  gapcolor=T)
        }
      } else {
        p <- list()
        for(i in 1:length(dnds)){
          p[[i]] <- create_msa_plot(dnds_path = dnds[i], 
                                    gap_path = gap[i], gapcolor=F)
        }
      }
    } else {
      p <- list()
      for(i in 1:length(dnds)){
        p[[i]] <- create_msa_plot(dnds_path = dnds[i], 
                                  gap_path = gap[i], gapcolor=F, points=F)
      }
    }
    
    do.call(grid.arrange,p)
    return(p)
  }
  
  # show TIGRFAM annotation for selected samples
  doPlotAnnotationGlobal <- function(data, input){
    if(substring(data$name, 1,4)[1] == "TIGR"){
      require(gridExtra)
      num <- as.data.frame(table(data$term))
      num <- num[which(num$Freq > 0),]
      p <- ggplot(num, aes(reorder(Var1, Freq), Freq)) + coord_flip()
      p <- p + geom_bar(stat = "identity", fill="grey80", width=0.8) + theme_classic()
      #p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
      p <- p + labs(x = "", y="number of protein families annotated") 
      if(input$pval < 1){
        p <- p + ggtitle(paste("Number of protein families (p-value less than: ", input$pval, ")", sep=""))
      } else {
        p <- p + ggtitle("Number of protein families in dataset")
      }
      return(p)
    }
  }
  
  doPlotSample <- function(data, input){
    require(ggplot2)
    num <- as.data.frame(table(data$sample))
    num$selected <- FALSE
    for(i in 1:length(input$samples)){
      num[which(num$Var1 == input$samples[i]),]$selected <- TRUE
    }
    p <- ggplot(num, aes(Var1, Freq, fill=selected)) 
    p <- p + geom_bar(stat = "identity", width=0.8) + theme_classic()
    p <- p + scale_fill_manual(breaks = c(TRUE, FALSE), values = c("grey80", "black")) + guides(fill=FALSE)
    p <- p + labs(x = "Samples", y="# protein families") + ggtitle("Selected samples")
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    return(p)
  }
  
  
  doAnnotationplot <- function(data, input){
    require(ggplot2)
    if(input$navalues){
      # remove NA values
      data <- data[which(!is.na(data$term)),]
    }
    if(input$sortannotation == "ratio"){
      if(input$bysamplecolor){
        p <- ggplot(data, aes(x=reorder(term, ratio), y=ratio, fill=sample))
        p <- p + geom_boxplot(width=0.3) + theme_classic() + coord_flip()
      } else {
        p <- ggplot(data, aes(x=reorder(term, ratio), y=ratio))
        p <- p + geom_boxplot(width=0.3, fill="grey80") + theme_classic() + coord_flip()
      }
    } else {
      p <- ggplot(data, aes(x=reorder(term, -fdr), y=ratio))
    }
    p <- p + ylab("dN/dS ratio") + xlab("functional group")
    if(input$showmean){
      p <- p + geom_hline(yintercept = mean(data$ratio, na.rm=T))
    }
    if(input$showmeanselected){
      p <- p + geom_hline(yintercept = mean(data[input$table_rows_selected,]$ratio, na.rm=T), color="red")
    }
    if(input$bysamplefacet){
      p <- p + facet_wrap(~ sample)
    }
    
    return(p)
  }
  
  doPlotBox <- function(data, input){
    require(ggplot2)
  
    if(input$oderchoice == "mean"){
      p <- ggplot(data, aes(x=reorder(sample, ratio),y=ratio))
    } else {
      p <- ggplot(data, aes(x=sample,y=ratio))
    }
    p <- p + geom_boxplot(fill="grey80", width=0.8) + theme_classic() + coord_flip()
    if (input$highlightbox){
      mark.ratio <- data[input$table_rows_selected,]$ratio
      p <- p + geom_hline(yintercept = mean(mark.ratio, na.rm=T))
    }
     
    return(p)
  }
  
  #################
  # render graphics functions
  #################
  
  output$plot1 <- renderPlot({
    data <-readData(paste(folder.path, input$dataset, sep="/"))
    data <- data[which(data$fdr <= input$pval),]
    data <- data[which(data$sample == input$samples),]
    p <- doPlotHistogram(data, input)
    print(p)
  }, height=700)
  
  # Density plot (unused)
  output$plot2 <- renderPlot({
    p <- ggplot(dataset(), aes(ratio, fill = sample)) +
      geom_density(adjust=input$densityBw)
    print(p)
  }, height=700)
  
  observe({
    if (input$close > 0) stopApp() # stop shiny
  })
  
  # boxplot
  output$sampleplot <- renderPlot({
    data <- readData(paste(folder.path, input$dataset, sep="/"))
    p <- doPlotSample(data, input)
    downloadableSamplePlot <<- p
    print(p)
  }, height=200)
  
  # Alignment plot
  output$alignmentplot <- renderPlot({
    data <-readData(paste(folder.path, input$dataset, sep="/"))
    data <- data[which(data$fdr <= input$pval),]
    data <- data[which(data$sample == input$samples),]
    data <- data[input$table_rows_selected,]
    if (ncol(data)>0){
      # get the dnds an gap paths
      p <- doAlignmentPlot(data, input)
      downloadableAlignmentPlot <<- p
    } else {
      df <- data.frame()
      p <- ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)
      downloadableAlignmentPlot <<- p
    }
    # else write an error msg that the user have to select some rows
  }, height=700)
  
  # boxplot
  output$plot4 <- renderPlot({
    data <- readData(paste(folder.path, input$dataset, sep="/"))
    data <- data[which(data$fdr <= input$pval),]
    if(length(input$samples) > 1){
      subset <- NULL
      data_pool <- NULL
      for(i in 1:length(input$samples)){
        subset <- data[which(data$sample == input$samples[i]),]
        data_pool <- rbind(subset,data_pool)
      }
      data <- data_pool
    } else {
      data <- data[which(data$sample == input$samples),]
    }
    
    p <- doPlotBox(data, input)
    downloadableBoxplot <<- p
    print(p)
  }, height=300)
  
  # boxplot
  output$annotationplot <- renderPlot({
    data <- readData(paste(folder.path, input$dataset, sep="/"))
    data <- data[which(data$fdr <= input$pval),]
  #  data <- data[which(data$sample == input$samples),]
    
    if(length(input$samples) > 1){
      subset <- NULL
      data_pool <- NULL
      for(i in 1:length(input$samples)){
        subset <- data[which(data$sample == input$samples[i]),]
        data_pool <- rbind(subset,data_pool)
      }
      data <- data_pool
    } else {
      data <- data[which(data$sample == input$samples),]
    }
    p <- doAnnotationplot(data, input)
    downloadableAnnotaionplot <<- p
    print(p)
  }, height=400)
  
  # annotation plot
  output$annotationplotglobal <- renderPlot({
    data <- readData(paste(folder.path, input$dataset, sep="/"))
    data <- data[which(data$fdr <= input$pval),]
    if(length(input$samples) > 1){
      subset <- NULL
      data_pool <- NULL
      for(i in 1:length(input$samples)){
        subset <- data[which(data$sample == input$samples[i]),]
        data_pool <- rbind(subset,data_pool)
      }
      data <- data_pool
    } else {
      data <- data[which(data$sample == input$samples),]
    }
    p <- doPlotAnnotationGlobal(data, input)
    downloadableAnnotaionplot <<- p
    print(p)
  }, height=400)
  

  #################
  # RENDER TEXT
  #################
  
  # print the selected indices
  output$selected = renderPrint({
    s = input$table_rows_selected
    if (length(s)) {
      cat('These rows are selected:\n\n')
      cat(s, sep = ', ')
    }
  })
  
  
  # print summary and selected rows
  output$summary <- renderPrint({  
    data <-readData(paste(folder.path, input$dataset, sep="/"))
    data <- data[which(data$fdr <= input$pval),]
    data <- data[which(data$sample == input$samples),]
    summary(data$ratio) 
    cat(paste(length(input$samples), "dataset(s) are selected\n"))
    cat(paste("mean ratio of selected datasets:", round(mean(data$ratio),digits = 3), "+-",round(sd(data$ratio),digits = 3) ,"SD"))
    cat("\n")
    s = input$table_rows_selected
    if (length(s)) {
      cat("\n")
      cat(paste(length(s),' protein families are selected:\n'))
      # subset the data based on selection
      data.selection <<- data[input$table_rows_selected, ]
      cat(paste("mean ratio:", round(mean(data.selection$ratio),digits = 3), "+-",round(sd(data.selection$ratio),digits = 3) ,"SD"))
      # perform fisher test 
      cat("\n")
      print(data.selection$name)
   
      test.mat <-
        matrix(c(sum(data.selection$sum_pN), sum(data.selection$sum_pS), sum(data$sum_pN), sum(data$sum_pS)),
               nrow = 2,
               dimnames =
                 list(c("dN", "dS"),
                      c("selected", "background")))
      cat("\n")
      pval <- fisher.test(test.mat, alternative = "greater")$p.value
      cat("one sided Fisher's test:")
      cat("\n")
      ifelse(pval < 0.05, print(paste("difference is significant (p-value=",round(pval, digits = 5), ")", sep="")), print(paste("difference is not significant (p-value=", round(pval, digits=3),")", sep="")))
      cat("\n")
    } else {
      cat("\nTo test if specific protein families are significant under selection please select them on the table above.")
    }
    
  })
  
  
  # print summary and selected rows
  output$alignment <- renderPrint({  
    data <-readData(paste(folder.path, input$dataset, sep="/"))
    data <- data[which(data$fdr <= input$pval),]
    data <- data[which(data$sample == input$samples),]
    
    if(length(input$samples)>1){
      cat(paste("Error: More than one dataset selected! A plot can only be created for one dataset. Please go back to the Data Table tab and deselect the dataset\n"))
    }
    s = input$table_rows_selected
    if (length(s)) {
      cat("\n")
      cat(paste(length(s),' protein families were selected:\n'))
    } else {
      cat("Error: No gene families selected! Please go to the Data Table tab and select one or more rows in the data table.")
    }
    
    cat(length(input$samples))
  })
  
  
  #################
  # RENDER html
  #################
  
  # tab 1
  # overview  tab
  output$overview_hint <- renderText({paste("</br><font color=\"#008080\"><b>", "Hint: You can include or remove samples from your analysis. For this use the 'Choose one or more samples' form on the widget on the left side.</b></font></br></br>")})
  
  output$overview_hint2 <- renderText({paste("</br><font color=\"#008080\"><b>", "Hint: You can perform a one-sided Fisher's exact test by selecting protein families by clicking on the table.</b></font></br>") })
  
  output$overview_table <- renderText({paste("Table:",num.name, " protein families found in", length(input$samples)," samples with a mean dN/dS ratio of", num.meanratio, " +- ", num.sd, "(SD). Categories based on HMM match with E-value 0.01. Only protein families with a FDR adjusted p-value of less than ", input$pval," are shown. p-value(s) as: one star for value below 0.05, two for 0.01 and three for 0.001. Table generated with EDEN v.1.</br></br>")})
  
  output$overview_fisher <- renderText({paste("Basic statistics and one-sided Fisher's test")})

  # tab 2
  # annotation tab
  output$annotation_hint <- renderText({paste("<font color=\"#008080\"><b>", "Hint: You can specify a filter to show protein families that have a significant or high dn/ds ratio. For this use the slider 'p-value threshold' on the widget on the left side.</b></font>")})
  
  output$annotation_figure <- renderText({paste("Figure: Number of protein families found in", length(input$samples)," samples for each category. Only protein families with a FDR adjusted p-value of less than ", input$pval," are shown. Figure generated with EDEN v.1.")})
  
  # tab 3
  # overview  tab
  output$alignment_hint <- renderText({paste("</br><font color=\"#008080\"><b>", "Visualisation of the dN/dS ratio for selected protein families. If multiple samples are selected the first sample will be plottet.</b></font></br></br>")})
  # overview  tab
  output$alignment_figure <- renderText({paste("Figure: Sequence clusters of residues under positive selection in selected protein families. Dots indicate dN/dS ratio for a given position in the protein sequence, and their color corresponds to the proportion of gaps in the multiple sequence alignment (MSA). Gray-shaded areas indicate significant clusters of residues under positive selection.</br></br>")})
  
  
  # boxplot tab
  output$boxplot_hint <- renderText({paste("</br><font color=\"#008080\"><b>", "Hint: Select multiple samples on the widget on the right side to compare the dN/dS ratio of the samples.</b></font></br></br>")})
  output$boxplot_figure <- renderText({paste("Figure: some figure text here</br></br>")})
  output$boxplot_table <- renderText({paste("Table: some table text here</br></br>")})
  
  
  
  
  #################
  # download handlers
  #################
  
  
  # download main table on the first tab
  output$dlTable <- downloadHandler(
    filename="table.csv",
    content = function(file) {
      write.csv(downloadObj, file)
    }
  )
  
  # download histogram
  output$dlCurPlot <- downloadHandler(
    filename = 'histogram.pdf',
    content = function(file){
      pdf(file = file, width=11, height=8.5)
      p <- doPlotHistogram(dataset(), input)
      print(p)
      dev.off()
    }
  )
  
  # download sequenceplot
  output$dlCurSequenceplot <- downloadHandler(
    filename = 'sequenceplot.pdf',
    content = function(file){
      pdf(file = file, width=11, height=8.5)
      print(downloadableAlignmentPlot)
      dev.off()
    }
  )
  
  # download boxplot
  output$dlCurBoxPlot <- downloadHandler(
    filename = 'boxplot.pdf',
    content = function(file){
      pdf(file = file, width=11, height=8.5)
      print(downloadableBoxplot)
      dev.off()
    }
  )
  
  # download sequenceplot
  output$dlAnnotationPlot <- downloadHandler(
    filename = 'annotationplot.pdf',
    content = function(file){
      pdf(file = file, width=11, height=8.5)
      print(downloadableAnnotaionplot)
      dev.off()
    }
  )
  
  # move uploaded files (tmp) to destination folder dest folders must be chmod 777
  ntextupload <- eventReactive(input$uploadButton, {
    # process faa
    dir.create(faa.path)
    Sys.chmod(faa.path, mode = "0777", use_umask = TRUE)
    infiles_faa <- as.data.frame(input$files_faa)
    infiles_faa$dest <- paste(faa.path, infiles_faa$name, sep="/")
    for (i in 1:nrow(infiles_faa)){
      cmd <- paste("mv ",infiles_faa$datapath[i]," ", infiles_faa$dest[i], sep="")
      err <- system(cmd,  intern = TRUE)
    }
    out <- paste(err)
    # process ffn
    dir.create(ffn.path)
    Sys.chmod(ffn.path, mode = "0777", use_umask = TRUE)
    infiles_ffn <- as.data.frame(input$files_ffn)
    infiles_ffn$dest <- paste(ffn.path, infiles_ffn$name, sep="/")
    for (i in 1:nrow(infiles_ffn)){
      cmd <- paste("mv ",infiles_ffn$datapath[i]," ", infiles_ffn$dest[i], sep="")
      err <- system(cmd,  intern = TRUE)
    }
    out <- paste(err)
  })
  
  ntextcheck <- eventReactive(input$checkButton, {
    std <- system2("/home/eden/run_check.sh", stdout=TRUE,stderr=TRUE)
  })
  
  ntexteden <- eventReactive(input$goButton, {
    #out <- system("/home/eden/eden.sh --docker --cpu_number 4 --gap_threshold 0.8 --test --name shiny", intern=TRUE)
  
    })
  
  
  
  # Function to get new log entries
  get_new_log <- function(){
    data <- read.table(log.path)
    return(data)
  }
  
  # Initialize log
  my_log<<- get_new_log()
  
  # Function to update my_data
  update_log <- function(){
    my_log <<- my_log
  }
  
  output$log = renderTable({
    invalidateLater(millis=5000, session)
    update_log()
  })
  
  output$nTextupload <- renderText({
    ntextupload()
  })
  
  
  output$nTextcheck <- renderText({
    ntextcheck()
  })
  
  output$nTexteden <- renderText({
    ntexteden()
  })
  

  output$filetable_faa<- reactiveTable(function() {
    if (is.null(input$files_faa)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    input$files_faa
  })

  output$filetable_ffn<- reactiveTable(function() {
    if (is.null(input$files_ffn)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    input$files_ffn
  })
  
  ##### log
  ##### 

  
})