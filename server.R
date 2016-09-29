shinyServer(function(input, output) {
  source("download.functions.R")

  #################
  # Reactive
  #################
  dataset <- reactive({
      # choose selected one 
      readData(paste(folder.path, input$dataset, sep="/"))
  })
  
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
    dnds <- paste(tmp.path,"/",input$data,"/",input$samples,"/dnds/", fam_ids,".txt.DnDsRatio.txt", sep="")
    gap <- paste(tmp.path,"/", input$data,"/",input$samples,"/gap/", fam_ids,".gap.txt", sep="")
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
    df <- NULL
    df <- data.frame(term=unique(data$term), pval=rep(-1,length(unique(data$term))))
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
      df[i,]$pval <- fisher.test(test.mat, alternative = "less")$p.value
      i <- i + 1 
    }
    df$fdr <- p.adjust(df$pval, method="fdr")
    annotation <<- df
    
    if(input$navalues){
      # remove NA values
      data <- data[which(!is.na(data$term)),]
    }
    if(input$sortannotation == "ratio"){
      if(input$bysamplecolor){
        p <- ggplot(data, aes(x=reorder(term, ratio), y=ratio, fill=sample))
      } else {
        p <- ggplot(data, aes(x=reorder(term, ratio), y=ratio))
      }
    } else {
      p <- ggplot(data, aes(x=reorder(term, -fdr), y=ratio))
    }
    p <- p + geom_boxplot(width=0.3) + theme_bw() + coord_flip()
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
    df <- NULL
    df <- data.frame(sample=unique(data$sample), pval=rep(-1,length(unique(data$sample))))
    df <- df[which(!is.na(df$sample)),]
    i <- 1
    for(sample in df$sample){
      data.sample <- data[which(data$sample == sample),]
      data.nonsample <- data[which(data$sample != sample),]
      test.mat <-
        matrix(c(sum(data.nonsample$sum_pN), sum(data.sample$sum_pN), sum(data.sample$sum_pS), sum(data.nonsample$sum_pS)),
               nrow = 2,
               dimnames =
                 list(c("background", "selected"),
                      c("dN", "dS")))
      df[i,]$pval <- fisher.test(test.mat, alternative = "less")$p.value
      i <- i + 1 
    }
    df$fdr <- p.adjust(df$pval, method="fdr")
    sample_pval <<- df
    
    if(input$oderchoice == "mean"){
      p <- ggplot(data, aes(x=reorder(sample, ratio),y=ratio)) +
        geom_boxplot() + theme_bw() + coord_flip()
    } else {
      p <- ggplot(data, aes(x=sample,y=ratio)) +
        geom_boxplot() + theme_bw() + coord_flip()
    }
    
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
  output$plot3 <- renderPlot({
    data <-readData(paste(folder.path, input$dataset, sep="/"))
    data <- data[which(data$fdr <= input$pval),]
    data <- data[which(data$sample == input$samples),]
    data <- data[input$table_rows_selected,]
    if (ncol(data)>0){
      # get the dnds an gap paths
      p <- doAlignmentPlot(data, input)
      downloadableAlignmentPlot <<- p
    }
    # else write an error msg that the user have to select some rows
  }, height=700)
  
  # boxplot
  output$plot4 <- renderPlot({
    data <- readData(paste(folder.path, input$dataset, sep="/"))
    data <- data[which(data$fdr <= input$pval),]
    data <- data[which(data$sample == input$samples),]
    p <- doPlotBox(data, input)
    downloadableBoxplot <<- p
    print(p)
  }, height=700)
  
  # boxplot
  output$annotationplot <- renderPlot({
    data <- readData(paste(folder.path, input$dataset, sep="/"))
    data <- data[which(data$fdr <= input$pval),]
    data <- data[which(data$sample == input$samples),]
    p <- doAnnotationplot(data, input)
    downloadableAnnotaionplot <<- p
    print(p)
  }, height=700)
  
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
    downloadableAnnotationGlobal <<- p
    print(p)
  }, height=400)
  
  #################
  # render table functions
  #################
  
  output$table_filtered <- DT::renderDataTable(
    DT::datatable(dataset, options = list(paging = 25))
  )
  
  output$table_annotaion <- DT::renderDataTable(
    DT::datatable(annotation, options = list(paging = 25))
  )
  
  output$table_sample <- DT::renderDataTable(
    DT::datatable(sample_pval, options = list(paging = 25))
  )
  
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
    #  colnames(data) <- c("family name", "dnds ratio", "significance lvl","adjusted pvalue", "dataset", "family description")
    data$fdr <- round(data$fdr, digits=5)
    num.name <<- nrow(data)
    num.meanratio <<- round(mean(data$ratio, na.rm=T), digits=2)
    num.sd <<- round(sd(data$ratio, na.rm=T), digits=2)
    downloadObj <<- data # generate downloadable table
    input$resetSelection 
    data
  }))
  
  
  
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
    cat(paste("Summary of selected", length(input$samples), "dataset(s):\n"))
    cat(paste("mean ratio:", round(mean(data$ratio),digits = 3), "+-",round(sd(data$ratio),digits = 3) ,"SD"))
    cat("\n")
    s = input$table_rows_selected
    if (length(s)) {
      cat("\n")
      cat(paste(length(s),' rows were selected:\n'))
      # subset the data based on selection
      data.selection <- data[input$table_rows_selected, ]
      cat(paste("mean ratio:", round(mean(data.selection$ratio),digits = 3), "+-",round(sd(data.selection$ratio),digits = 3) ,"SD"))
      # perform fisher test 
      cat("\n")
      print(data.selection$name)
      cat("\n")
      test.mat <-
        matrix(c(sum(data$sum_pN), sum(data.selection$sum_pN), sum(data.selection$sum_pS), sum(data$sum_pS)),
               nrow = 2,
               dimnames =
                 list(c("background", "selected"),
                      c("dN", "dS")))
      cat("\n\n")
      pval <- fisher.test(test.mat, alternative = "less")$p.value
      cat("one sided Fisher's test:")
      cat("\n")
      ifelse(pval < 0.05, print(paste("difference is significant",pval)), print(paste("difference is not significant (pvalue=", pval,")", sep="")))
      cat("\n")
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
  
  output$overview_table <- renderText({paste("Table:",num.name, " protein families found in", length(input$samples)," samples with a mean dN/dS ratio of", num.meanratio, " +- ", num.sd, "(SD). Categories based on HMM match with E-value 0.01. Only protein families with a FDR adjusted p-value of less than ", input$pval," are shown. Table generated with EDEN v.1.</br>")})
  
  # tab 2
  # annotation tab
  output$annotation_hint <- renderText({paste("<font color=\"#008080\"><b>", "Hint: You can specify a filter to show protein families that have a significant or high dn/ds ratio. For this use the slider 'p-value threshold' on the widget on the left side.</b></font>")})
  
  output$annotation_figure <- renderText({paste("Figure: Number of protein families found in", length(input$samples)," samples for each category. Only protein families with a FDR adjusted p-value of less than ", input$pval," are shown. Figure generated with EDEN v.1.")})
  
  #################
  # download functions
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
  output$dlCurAnnotationplot <- downloadHandler(
    filename = 'annotationplot.pdf',
    content = function(file){
      pdf(file = file, width=11, height=8.5)
      print(downloadableAnnotaionplot)
      dev.off()
    }
  )
  
})