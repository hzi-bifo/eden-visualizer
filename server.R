library(shiny)
library(ggplot2)

shinyServer(function(input, output) {
  #################
  # Load data from results folder
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
  #  input$dataset
  dataset <- readData(paste(folder.path, input$dataset, sep="/"))
    selectInput("samples", "Choose one or more datasets:", choices=levels(factor(dataset$sample)), selected=c(levels(factor(dataset$sample))[1]), multiple=T, width="100%")
  })
  
  #################
  # RENDER PLOTS
  #################
  
  # Histogram
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
      #      p <- p + geom_text(aes(x=mark.ratio, label=mark.name, y=rep(1, mark.num)), colour="black", angle=90)
    }
    return(p)
  }
  
  doAlignmentPlot <- function(data, input){
    require(ggplot2)
    require(grid)
    require(gridExtra)
    fam_ids <- data$name
    dnds <- paste(tmp.path,"/",input$dataset,"/",input$samples,"/dnds/", fam_ids,".txt.DnDsRatio.txt", sep="")
    gap <- paste(tmp.path,"/", input$dataset,"/",input$samples,"/gap/", fam_ids,".gap.txt", sep="")
    
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
    p <- try(doPlotBox(data, input))
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
  
  #################
  # RENDER TEXT
  #################
  
  # print the selected indices
  output$selected = renderPrint({
    s = input$table_rows_selected
    if (length(s)) {
      cat('These rows were selected:\n\n')
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
      ifelse(pval < 0.05, print(paste("difference is significant",pval)), print(paste("difference is not significant", pval)))
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
  # RENDER TABLES
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
      print("multiple samples selected")
      subset <- NULL
      data_pool <- NULL
      for(i in 1:length(input$samples)){
        print(paste("check to", input$samples[i], sep=" "))
        subset <- data[which(data$sample == input$samples[i]),]
        data_pool <- rbind(subset,data_pool)
        print(paste("subset length is", nrow(subset), sep=" "))
      }
      data <- data_pool
    } else {
      data <- data[which(data$sample == input$samples),]
      print("length is one")
    }
 #   data$stars <- add.significance.stars(data$fdr)
    data$sum_pN <- NULL
    data$sum_pS <- NULL
    data$role <- NULL
    data$pvalue <- NULL
    data$ratio <- round(data$ratio, digits=3)
  #  colnames(data) <- c("family name", "dnds ratio", "significance lvl","adjusted pvalue", "dataset", "family description")
    data$fdr <- round(data$fdr, digits=5)
  #  downloadObj <<- data # generate downloadable table
  #  input$resetSelection 
    data
  }))
  
  #################
  # DOWNLOAD
  #################
  
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
  
  output$dlTable <- downloadHandler(
    filename="table.csv",
    content = function(file) {
      write.csv(downloadObj, file)
    }
  )
  
  
  
})