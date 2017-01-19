
shinyServer(function(input, output, session) {
  

  dataset <- reactive({
    if (input$dataset[1] != "") {
      readCsv(paste(csv.path, input$dataset, sep = "/"))
    }
  })

  output$table_filtered <- DT::renderDataTable(DT::datatable(dataset, options = list(paging = 25)))
  output$table_annotaion <- DT::renderDataTable(DT::datatable({
    require(pander)
    data <-   readCsv(paste(csv.path, input$dataset, sep = "/"))
    if (length(input$samples) > 1) {
      subset <- NULL
      data_pool <- NULL
      for (i in 1:length(input$samples)) {
        subset <- data[which(data$sample == input$samples[i]), ]
        data_pool <- rbind(subset, data_pool)
      }
      data <- data_pool
    } else {
      data <- data[which(data$sample == input$samples), ]
    }
    
    df <- NULL
    df <-
      data.frame(
        term = unique(data$term),
        pval = rep(-1, length(unique(data$term))),
        elements = rep(0, length(unique(data$term)))
      )
    df <- df[which(!is.na(df$term)), ]
    i <- 1
    for (term in df$term) {
      data.term <- data[which(data$term == term), ]
      data.nonterm <- data[which(data$term != term), ]
      test.mat <-
        matrix(c(
          sum(data.term$sum_pN),
          sum(data.term$sum_pS),
          sum(data.nonterm$sum_pN),
          sum(data.nonterm$sum_pS)
        ),
        nrow = 2,
        dimnames =
          list(c("background", "selected"),
               c("dN", "dS")))
      df[i, ]$pval <-
        fisher.test(test.mat, alternative = "greater")$p.value
      df[i, ]$elements <- df[i, ]$elements  + nrow(data.term)
      i <- i + 1
    }
    df$fdr <- p.adjust(df$pval, method = "fdr")
    df$fdr <- round(df$fdr, digits = 6)
    df$star <- add.significance.stars(df$fdr)
    df$pval <- NULL
    df
  }))
  
  output$table_sample <- DT::renderDataTable(DT::datatable({
    require(pander)
    data <-   readCsv(paste(csv.path, input$dataset, sep = "/"))
    data <- data[which(data$fdr <= input$pval), ]
    if (length(input$samples) > 1) {
      subset <- NULL
      data_pool <- NULL
      for (i in 1:length(input$samples)) {
        subset <- data[which(data$sample == input$samples[i]), ]
        data_pool <- rbind(subset, data_pool)
      }
      data <- data_pool
    } else {
      data <- data[which(data$sample == input$samples), ]
    }
    
    df <- NULL
    df <-
      data.frame(sample = unique(data$sample),
                 pval = rep(-1, length(unique(data$sample))))
    df <- df[which(!is.na(df$sample)), ]
    i <- 1
    for (sample in df$sample) {
      data.sample <- data[which(data$sample == sample), ]
      data.nonsample <- data[which(data$sample != sample), ]
      test.mat <-
        matrix(c(
          sum(data.sample$sum_pN),
          sum(data.sample$sum_pS),
          sum(data.nonsample$sum_pN),
          sum(data.nonsample$sum_pS)
        ),
        nrow = 2,
        dimnames =
          list(c("dN", "dS"),
               c("selected", "background")))
      df[i, ]$pval <-
        fisher.test(test.mat, alternative = "greater")$p.value
      i <- i + 1
    }
    
    df$fdr <- p.adjust(df$pval, method = "fdr")
    df$pval <- NULL
    df$fdr <- round(df$fdr, digits = 6)
    df$star <- add.significance.stars(df$fdr)
    df
  }))
  
  output$table <- DT::renderDataTable(DT::datatable(dataset, options = list(pageLength = 25)))
  
  # Filter data based on selections
  output$table <- DT::renderDataTable(DT::datatable({
    require(pander)
    data <-  readCsv(paste(csv.path, input$dataset, sep = "/"))
    data <- data[which(data$fdr <= input$pval), ]
    if (length(input$samples) > 1) {
      subset <- NULL
      data_pool <- NULL
      for (i in 1:length(input$samples)) {
        subset <- data[which(data$sample == input$samples[i]), ]
        data_pool <- rbind(subset, data_pool)
      }
      data <- data_pool
    } else {
      data <- data[which(data$sample == input$samples), ]
    }
    data$stars <- add.significance.stars(data$fdr)
    data$sum_pN <- NULL
    data$sum_pS <- NULL
    data$role <- NULL
    data$pvalue <- NULL
    data$ratio <- round(data$ratio, digits = 3)
    data$fdr <- round(data$fdr, digits = 5)
    num.name <<- nrow(data)
    num.meanratio <<- round(mean(data$ratio, na.rm = T), digits = 2)
    num.sd <<- round(sd(data$ratio, na.rm = T), digits = 3)
    downloadObj <<- data # generate downloadable table
    input$resetSelection
    data
  }))
  
  #################
  # Render UI
  #################
  

  
  # render again if input$player_name changes
  output$filters_UI <- renderUI({
   # if (!no_csv){
    dataset <- readCsv(paste(csv.path, input$dataset, sep = "/"))
    selectInput(
      "samples",
      "Choose one or more samples:",
      choices = levels(factor(dataset$sample)),
      selected = c(levels(factor(
        dataset()$sample
      ))[1]),
      multiple = T,
      width = "100%"
    ) 
    #}
    })
  
  # render again if input$dofiltering changes
  output$dependentselection <- renderUI({
    if (input$dofiltering == "pvalue") {
      sliderInput(
        "pval",
        label = "adjusted p-value threshold",
        min = .001,
        max = 1,
        value = 1
      )
    } else {
      sliderInput(
        "ratio",
        label = "select ratio range to display",
        min = round(min(data$ratio), digits = 2),
        max = round(max(data$ratio), digits = 2),
        value = c(round(min(data$ratio), digits = 2),
                  round(max(data$ratio), digits = 2))
      )
    }
  })
  
  ### start main ui
  output$main_ui <- renderUI({
    conditionalPanel(
      condition = "input.tsp=='overview' ||
      input.tsp=='annotation' ||
      input.tsp=='alignment' ||
      input.tsp=='histogram' ||
      input.tsp=='box' ||
      input.tsp=='start' ||
      input.tsp=='categories' ",
      helpText("Select which analysis run you want to show"),
      selectInput(
        "dataset",
        "Select run:",
        choices = list.dirs(
          path = csv.path,
          full.names = FALSE,
          recursive = FALSE
        ),
        selected = list.dirs(
          path = "csv",
          full.names = FALSE,
          recursive = FALSE
        )[1],
        multiple = F,
        width = "100%"
      ),
      
      uiOutput('filters_UI'),
      
      selectInput(
        'dofiltering',
        label = 'Choose a way to filter matrix',
        choices = c("pvalue", "ratio (not implemented)"),
        selected = "no filtering"
      ),
      uiOutput("dependentselection"),
      actionButton('resetSelection', label = "Reset row selection",   class =
                     "btn-block btn-primary")
    )
  })
  
  ### end main ui
  
  
  # render again if input$player_name changes
  output$start_UI_samples <- renderUI({
    if (input$analysistype == "comparative") {
      conditionalPanel(
        condition = "input.tsp=='start' || input.tsp=='log'",
        helpText(
          "For a comparative analysis please provide information which samples should be pooled together"
        ),
        fileInput(
          'file_sample',
          'upload a sample description file',
          accept = c('.txt'),
          multiple = FALSE
        )#,
        #   checkboxInput("eden_use_mgm", "find ORFs with MetaGeneMark", FALSE)
      )
    }
  })
  
  
  # render again if input$dofiltering changes
  output$colorpoints <- renderUI({
    if (input$points) {
      checkboxInput('gap', 'color by gap proportion', value = TRUE)
    }
  })
  
  
  #################
  # ggplot functions
  #################
  
  # generates a histogram
  doPlotHistogram <- function(dataset, input) {
    p <- ggplot(dataset, aes(ratio, fill = sample)) +
      geom_histogram(bins = input$binSize) + theme_classic()
    p <-
      p + labs(x = "dN/dS ratio", y = "protein families") + ggtitle("Histogram")
    if (input$facet)
      p <- p + facet_grid(sample ~ .)
    
    if (length(input$table_rows_selected)) {
      # get the ratio of selected rows
      mark.ratio <- dataset[input$table_rows_selected, ]$ratio
      mark.name <- dataset[input$table_rows_selected, ]$name
      p <- p + geom_vline(xintercept = mark.ratio)
    }
    return(p)
  }
  
  # this functions calls the create_msa_plot() function multiple times based
  # on selected protein families
  doAlignmentPlot <- function(data, input) {
    require(ggplot2)
    require(grid)
    require(gridExtra)
    fam_ids <- data$name
    dnds <-
      paste(
        raw.path,
        "/",
        input$dataset[1],
        "/",
        input$samples,
        "/dnds/",
        fam_ids,
        ".txt.DnDsRatio.txt",
        sep = ""
      )
    gap <-
      paste(
        raw.path,
        "/",
        input$dataset[1],
        "/",
        input$samples,
        "/gap/",
        fam_ids,
        ".gap.txt",
        sep = ""
      )
    if (input$points) {
      if (input$gap) {
        # get list of ggplot obj with gap color
        p <- list()
        for (i in 1:length(dnds)) {
          p[[i]] <- create_msa_plot(
            dnds_path = dnds[i],
            gap_path = gap[i],
            gapcolor = T
          )
        }
      } else {
        # get list of ggplot obj without gap color
        p <- list()
        for (i in 1:length(dnds)) {
          p[[i]] <- create_msa_plot(
            dnds_path = dnds[i],
            gap_path = gap[i],
            gapcolor = F
          )
        }
      }
    } else {
      p <- list()
      for (i in 1:length(dnds)) {
        p[[i]] <- create_msa_plot(
          dnds_path = dnds[i],
          gap_path = gap[i],
          gapcolor = F,
          points = F
        )
      }
    }
    
    do.call(grid.arrange, p)
    return(p)
  }
  
  # show TIGRFAM annotation for selected samples
  doPlotAnnotationGlobal <- function(data, input) {
    if (substring(data$name, 1, 4)[1] == "TIGR") {
      require(gridExtra)
      num <- as.data.frame(table(data$term))
      num <- num[which(num$Freq > 0), ]
      p <-
        ggplot(num, aes(reorder(Var1, Freq), Freq)) + coord_flip()
      p <-
        p + geom_bar(stat = "identity",
                     fill = "grey80",
                     width = 0.8) + theme_classic()
      #p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
      p <-
        p + labs(x = "", y = "number of protein families annotated")
      if (input$pval < 1) {
        p <-
          p + ggtitle(
            paste(
              "Number of protein families (p-value less than: ",
              input$pval,
              ")",
              sep = ""
            )
          )
      } else {
        p <- p + ggtitle("Number of protein families in dataset")
      }
      return(p)
    }
  }
  
  doPlotSample <- function(data, input) {
    require(ggplot2)
    num <- as.data.frame(table(data$sample))
    num$selected <- FALSE
    for (i in 1:length(input$samples)) {
      num[which(num$Var1 == input$samples[i]), ]$selected <- TRUE
    }
    p <- ggplot(num, aes(Var1, Freq, fill = selected))
    p <-
      p + geom_bar(stat = "identity", width = 0.8) + theme_classic()
    p <-
      p + scale_fill_manual(breaks = c(TRUE, FALSE),
                            values = c("grey80", "black")) + guides(fill = FALSE)
    p <-
      p + labs(x = "Samples", y = "# protein families") + ggtitle("Selected samples")
    p <-
      p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    return(p)
  }
  
  doAnnotationplot <- function(data, input) {
    require(ggplot2)
    if (input$navalues) {
      # remove NA values
      data <- data[which(!is.na(data$term)), ]
    }
    if (input$sortannotation == "ratio") {
      if (input$bysamplecolor) {
        p <- ggplot(data, aes(
          x = reorder(term, ratio),
          y = ratio,
          fill = sample
        ))
        p <-
          p + geom_boxplot(width = 0.3) + theme_classic() + coord_flip()
      } else {
        p <- ggplot(data, aes(x = reorder(term, ratio), y = ratio))
        p <-
          p + geom_boxplot(width = 0.3, fill = "grey80") + theme_classic() + coord_flip()
      }
    } else {
      p <- ggplot(data, aes(x = reorder(term,-fdr), y = ratio))
    }
    p <- p + ylab("dN/dS ratio") + xlab("functional group")
    if (input$showmean) {
      p <- p + geom_hline(yintercept = mean(data$ratio, na.rm = T))
    }
    if (input$showmeanselected) {
      p <-
        p + geom_hline(yintercept = mean(data[input$table_rows_selected, ]$ratio, na.rm =
                                           T),
                       color = "red")
    }
    if (input$bysamplefacet) {
      p <- p + facet_wrap( ~ sample)
    }
    
    return(p)
  }

  doPlotBox <- function(data, input) {
    require(ggplot2)
    if (input$oderchoice == "mean") {
      p <- ggplot(data, aes(x = reorder(sample, ratio), y = ratio))
    }
    # if(input$oderchoice == "pvalue"){
    #    p <- ggplot(data, aes(x=reorder(sample, -fdr),y=ratio))
    #  }
    if (input$oderchoice == "default") {
      p <- ggplot(data, aes(x = sample, y = ratio))
    }
    p <- p + ylab("dN/dS ratio") + xlab("sample")
    p <-
      p + geom_boxplot(fill = "grey80", width = 0.8) + theme_classic() + coord_flip()
    if (input$highlightbox) {
      mark.ratio <- data[input$table_rows_selected, ]$ratio
      p <- p + geom_hline(yintercept = mean(mark.ratio, na.rm = T))
    }
    
    return(p)
  }

  
  output$plot1 <- renderPlot({
    data <-  readCsv(paste(csv.path, input$dataset, sep = "/"))
    data <- data[which(data$fdr <= input$pval), ]
    data <- data[which(data$sample == input$samples), ]
    p <- doPlotHistogram(data, input)
    histogram_obj <<- p
    print(p)
  }, height = 700)
  
  # Density plot (unused)
  output$plot2 <- renderPlot({
    p <- ggplot(dataset(), aes(ratio, fill = sample)) +
      geom_density(adjust = input$densityBw)
    print(p)
  }, height = 700)
  
  observe({
    if (input$close > 0)
      stopApp() # stop shiny
  })
  
  # boxplot
  output$sampleplot <- renderPlot({
    data <-  readCsv(paste(csv.path, input$dataset, sep = "/"))
    p <- doPlotSample(data, input)
    downloadableSamplePlot <<- p
    print(p)
  }, height = 200)
  
  output$alignmentplot <- renderPlot({
    data <-  readCsv(paste(csv.path, input$dataset, sep = "/"))
    data <- data[which(data$fdr <= input$pval), ]
    data <- data[which(data$sample == input$samples), ]
    data <- data[input$table_rows_selected, ]
    if (length(input$table_rows_selected) > 0) {
      # get the dnds an gap paths
      p <- doAlignmentPlot(data, input)
      downloadableAlignmentPlot <<- p
    } else {
      df <- data.frame()
      p <- ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)
      downloadableAlignmentPlot <<- p
    }
    # else write an error msg that the user have to select some rowss
  }, height = 700)
  ####
  #### BOXPLOT TAB
  ####
  
  # boxplot
  output$plot4 <- renderPlot({
    data <-  readCsv(paste(csv.path, input$dataset, sep = "/"))
    data <- data[which(data$fdr <= input$pval), ]
    if (length(input$samples) > 1) {
      subset <- NULL
      data_pool <- NULL
      for (i in 1:length(input$samples)) {
        subset <- data[which(data$sample == input$samples[i]), ]
        data_pool <- rbind(subset, data_pool)
      }
      data <- data_pool
    } else {
      data <- data[which(data$sample == input$samples), ]
    }
    
    p <- doPlotBox(data, input)
    downloadableBoxplot <<- p
    print(p)
  }, height = 300)
  
  # boxplot
  output$annotationplot <- renderPlot({
    data <-  readCsv(paste(csv.path, input$dataset, sep = "/"))
    data <- data[which(data$fdr <= input$pval), ]
    #  data <- data[which(data$sample == input$samples),]
    
    if (length(input$samples) > 1) {
      subset <- NULL
      data_pool <- NULL
      for (i in 1:length(input$samples)) {
        subset <- data[which(data$sample == input$samples[i]), ]
        data_pool <- rbind(subset, data_pool)
      }
      data <- data_pool
    } else {
      data <- data[which(data$sample == input$samples), ]
    }
    p <- doAnnotationplot(data, input)
    catplot_obj <<- p
    print(p)
  }, height = 400)
  
  # annotation plot
  output$annotationplotglobal <- renderPlot({
    data <-  readCsv(paste(csv.path, input$dataset, sep = "/"))
    data <- data[which(data$fdr <= input$pval), ]
    if (length(input$samples) > 1) {
      subset <- NULL
      data_pool <- NULL
      for (i in 1:length(input$samples)) {
        subset <- data[which(data$sample == input$samples[i]), ]
        data_pool <- rbind(subset, data_pool)
      }
      data <- data_pool
    } else {
      data <- data[which(data$sample == input$samples), ]
    }
    p <- doPlotAnnotationGlobal(data, input)
    downloadableAnnotaionplot <<- p
    print(p)
  }, height = 400)
  
  
  # print the selected indices
  output$selected = renderPrint({
    s = input$table_rows_selected
    if (length(s)) {
      cat('These rows are selected:\n\n')
      cat(s, sep = ', ')
    }
  })
  
  
  output$summary2 <-
    renderText({
      data <-  readCsv(paste(csv.path, input$dataset, sep = "/"))
      data <- data[which(data$fdr <= input$pval), ]
      data <- data[which(data$sample == input$samples), ]

      s = input$table_rows_selected
      if (length(s)) {
        data.selection <<- data[input$table_rows_selected,]
        
        test.mat <-
          matrix(c(
            sum(data.selection$sum_pN),
            sum(data.selection$sum_pS),
            sum(data$sum_pN),
            sum(data$sum_pS)
          ),
          nrow = 2,
          dimnames =
            list(c("dN", "dS"),
                 c("selected", "background")))
  
        pval <- fisher.test(test.mat, alternative = "greater")$p.value
        ## selected
        paste(
          "</br><div class='panel panel-default'>
        <div class='panel-heading'>Fisher's test</div>
          <div class='panel-body'>",
          "mean ratio of selected datasets: <span class='badge'>",
          
          round(mean(data$ratio, na.rm=TRUE), digits = 3),
          "+-",
          round(sd(data$ratio, na.rm=TRUE), digits = 3) ,
          "(SD) </span></br> compared to  mean ratio of selected families: <span class='badge'>",
        round(mean(data.selection$ratio, na.rm=TRUE), digits = 3),
        "+-",
        round(sd(data.selection$ratio, na.rm=TRUE), digits = 3) ,
        "(SD) </span></br> p-value (one-sided):  <span class='badge'>", round(pval, digits = 5),"</span>"
      )
      } else {

        paste(
          "</br><div class='panel panel-default'>
        <div class='panel-heading'>Fisher's test</div>
          <div class='panel-body'>",
          "please select rows from the table above to perform a fisher test")
              }
      

    })
  
  
  
  
  
  
  # print summary and selected rows
  output$alignment <- renderPrint({
    data <-  readCsv(paste(csv.path, input$dataset, sep = "/"))
    data <- data[which(data$fdr <= input$pval), ]
    data <- data[which(data$sample == input$samples), ]
    
    if (length(input$samples) > 1) {
      cat(
        paste(
          "Error: More than one dataset selected! A plot can only be created for one dataset. Please go back to the Data Table tab and deselect the dataset\n"
        )
      )
    }
    s = input$table_rows_selected
    if (length(s)) {
      cat("\n")
      cat(paste(length(s), ' protein families were selected:\n'))
    } else {
      cat(
        "Error: No gene families selected! Please go to the Data Table tab and select one or more rows in the data table."
      )
    }
    
    cat(length(input$samples))
  })
  
  

  
  output$reloadmsg <- renderPrint({
    if (input$reloadButton) {
      withProgress(message = 'Extract files, please wait', value = 0, {
        extractTar(tar.path, raw.path, csv.path, progress=TRUE)
      })
     # if(!no_csv){
        updateSelectInput(session,
                          "dataset",  label = "Select run", choices = list.dirs(
                            path = csv.path,
                            full.names = FALSE,
                            recursive = FALSE
                          )
       )
        updateSelectInput(session, "samples")
        new_files <<- FALSE
      #}
    }
  })

  #################
  # RENDER html
  #################
  
  output$welcome <-
    renderText({
     
        paste("</br><font color=\"#4db898\"><b>Welcome to eden!</b></br></font>",
              "Eden is a fast implementation of the widely used method for the detection of protein families that are under positive selection based on the ratio of amino acid replacement versus silent substitution rates (dN/dS) that can applied an large metagenomic samples.",
              "</br></br>", "<font color=\"#4db898\"><b>About the examples</b></br></font>", 
              "On the left panel you can select example datasets we have computed for you. In the bodysites example we used over 60 metagenomic samples from healthy individuals from the Human Microbiome Project. In a second example named bmi.tar we analyzed over 50 metagenomic samples from the gut of lean, overwight and obese individuals.",
              "</br></br>  <span class='label label-default'>v. 0.1.0</span>"
          
              )
    
    })
  
  output$reloadstatus <-
    renderText({
      if (new_files){
        paste("<font color=\"#4db898\"><b>New samples detected! You have to import these samples to make them visible",
              "</b></font>")
        
      } else 
        " "
    })
 
  # print if new files found 
  output$newtar <-
    renderText({
      paste("<font color=\"#4db898\"><b>New tar::",new_files,
            "</b></font>")
    })
  
  
  
  # print if csv file is found or not
  output$selected_dataset <-
    renderText({
        paste("<font color=\"#4db898\"><b>Dataset selected:",input$dataset,
              "</b></font>")
    })
  
  output$selected_samples <-
    renderText({
      paste("<font color=\"#4db898\"><b>Samples selected: ",input$samples,
            "</b></font>")
    })
  
  
  # print if csv file is found or not
  output$csv_check <-
    renderText({
      if(no_csv){
      paste("<font color=\"#4db898\"><b>",
            "No csv files found!</b></font>")
      } else {
        paste("<font color=\"#4db898\"><b>",
              "Csv files found!</b></font>")
      }
    })
  
  # print if .tar files are found or not
  output$tar_check <-
    renderText({
      if(no_tar){
        paste("<font color=\"#4db898\"><b>",
              "No tar files found!</b></font>")
      } else {
        paste("<font color=\"#4db898\"><b>",
              "Tar files found!</b></font>")
      }
    })
  
  output$start_hint_online <-
    renderText({
      paste("</br><font color=\"#4db898\"><b>",
            "Welcome text here</b></font></br></br>")
    })
  
  # tab 1
  # overview  tab
  output$overview_hint <-
    renderText({
      paste(
        
        "<div class='alert alert-dismissible alert-warning'>
          <button type='button' class='close' data-dismiss='alert'>&times;</button>
          <h4>Hint!</h4>
          <p>You can include or remove samples from your analysis. For this use the <strong>Choose one or more samples</strong> form on the widget on the left side.</p>
        </div>"
        
  
      )
    })
  
  output$overview_hint2 <-
    renderText({
      paste(
        
        "<div class='alert alert-dismissible alert-warning'>
          <button type='button' class='close' data-dismiss='alert'>&times;</button>
          <p>You can perform a one-sided Fisher's exact test by selecting protein families by clicking on the table</p>
        </div>"
        
      )
    })
  
  output$overview_table <-
    renderText({
      paste(
        "</br><div class='panel panel-default'>
        <div class='panel-heading'>Table description</div>
          <div class='panel-body'>",
        "<span class='badge'>",
        num.name,
        "</span> protein families found in <span class='badge'>",
        length(input$samples),
        "</span> samples with a mean dN/dS ratio of <span class='badge'>",
        num.meanratio,
        " +- ",
        num.sd,
        "(SD) </span>. Categories based on HMM match with E-value 0.01. Only protein families with a FDR adjusted p-value of less than ",
        input$pval,
        " are shown. p-value(s) as: one star for value below 0.05, two for 0.01 and three for 0.001. Table generated with eden <span class='label label-default'>v. 0.1.0</span>",
      "</div></div>"
      )
    })
  
  
  # tab 2
  # annotation tab
  output$annotation_hint <-
    renderText({
      paste(
        
        "<div class='alert alert-dismissible alert-warning'>
  <h4>Just want to show gene families that are significant?</h4>
          <button type='button' class='close' data-dismiss='alert'>&times;</button>
          <p>You can specify a filter to show protein families that have a significant or high dn/ds ratio. For this use the slider <strong>p-value threshold</strong> on the widget on the left side.</p>
        </div>"
        
      )
    })
  
  output$annotation_figure <-
    renderText({
      paste(
        
        "</br><div class='panel panel-default'>
        <div class='panel-heading'>Figure description</div>
          <div class='panel-body'>",
        
        "Number of protein families found in <span class='badge'>",
        length(input$samples),
        "</span> samples for each category. Only protein families with a FDR adjusted p-value of less than <span class='badge'>",
        input$pval,
        "</span> are shown. Figure generated with eden <span class='label label-default'>v. 0.1.0</span></div>"
      )
    })
  
  # tab 3
  # overview  tab
  output$alignment_hint <-
    renderText({
      if (length(input$table_rows_selected) < 1) {
        paste( "<div class='alert alert-dismissible alert-danger'>
  <h4>You need to select protein families first</h4>
          <button type='button' class='close' data-dismiss='alert'>&times;</button>
          <p>Just go back to the <strong>Overview</strong> tab and select rows you want to show here</p>
        </div>")
      }
    
    })
  # overview  tab
  output$alignment_figure <-
    renderText({
      if(length(input$table_rows_selected)>0){
        paste(
          "</br><div class='panel panel-default'>
        <div class='panel-heading'>Figure description</div>
          <div class='panel-body'>",
          
          "Sequence clusters of residues under positive selection in selected protein families. Dots indicate dN/dS ratio for a given position in the protein sequence, and their color corresponds to the proportion of gaps in the multiple sequence alignment (MSA). Gray-shaded areas indicate significant clusters of residues under positive selection.</div>"
        )

        
      }
    })
  
  
  # boxplot tab
  output$boxplot_hint <-
    renderText({
      paste(
        "<div class='alert alert-dismissible alert-warning'>
          <button type='button' class='close' data-dismiss='alert'>&times;</button>
          <h4>Hint!</h4>
          <p>You can include or remove samples from your analysis. For this use the <strong>Choose one or more samples</strong> form on the widget on the left side.</p>
        </div>"
      )
    })


  #################
  # download handlers
  #################
  
  # download main table on the first tab
  output$dlTable <- downloadHandler(
    filename = "table.csv",
    content = function(file) {
      write.csv(downloadObj, file)
    }
  )
  
  # download histogram
  output$dlCurPlot <- downloadHandler(
    filename = 'histogram.pdf',
    content = function(file) {
      pdf(file = file,
          width = 11,
          height = 8.5)
     
      print(histogram_obj)
      dev.off()
    }
  )
  
  # download sequenceplot
  output$dlCurSequenceplot <- downloadHandler(
    filename = 'sequenceplot.pdf',
    content = function(file) {
      pdf(file = file,
          width = 11,
          height = 8.5)
      print(downloadableAlignmentPlot)
      dev.off()
    }
  )
  
  # download boxplot
  output$dlCurBoxPlot <- downloadHandler(
    filename = 'boxplot.pdf',
    content = function(file) {
      pdf(file = file,
          width = 11,
          height = 8.5)
      print(downloadableBoxplot)
      dev.off()
    }
  )
  
  
  # download categorie plot
  output$dlCurAnnotationplot <- downloadHandler(
    filename = 'boxplot_categories.pdf',
    content = function(file) {
      pdf(file = file,
          width = 11,
          height = 8.5)
      print(catplot_obj)
      dev.off()
    }
  )
  
  # download sequenceplot
  output$dlAnnotationPlot <- downloadHandler(
    filename = 'annotationplot.pdf',
    content = function(file) {
      pdf(file = file,
          width = 11,
          height = 8.5)
      print(downloadableAnnotaionplot)
      dev.off()
    }
  )
})
