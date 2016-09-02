# iterate over input data and generate mat
readData <- function(path.summary){
  require(reshape2)
  file.names <- paste(path.summary, dir(path.summary, pattern =".csv"), sep="/")
  sample.names <- colsplit(string=dir(path.summary, pattern =".csv"), pattern=".csv", names=c("name", "ending"))$name
  mat <- NULL
  for(i in 1:length(file.names)){
    file <- read.table(file.names[i],header=TRUE, sep=";", stringsAsFactors=FALSE)
    file$sample <- NULL
    file$sample <- sample.names[i]
    mat <- rbind(mat, file)
  }
  return(mat)
}

# iterate over rows and run one sided fisher test with same sample as background
fisher_test <- function(dnds_summary){
  # input: name;sum_pN;sum_pS;ratio
  # output: name;sum_pN;sum_pS;ratio;pvalue;fdr;stars
  require(pander)
  dnds_summary$pvalue <- NULL
  dnds_summary$fdr <- NULL
  dnds_summary$stars <- "ns"
  for (i in 1:nrow(dnds_summary)) {
    TestMartix <-
      matrix(c(
        dnds_summary$sum_pN[i],
        dnds_summary$sum_pS[i], 
        sum(dnds_summary$sum_pN)-dnds_summary$sum_pN[i],
        sum(dnds_summary$sum_pS)-dnds_summary$sum_pS[i]),
        nrow = 2,
        dimnames = list(type = c("dN", "dS"),
                        sample = c("case", "control")))
    dnds_summary$pvalue[i] <- fisher.test(TestMartix, alternative = "greater")$p.value
  }
  # fdr correction
  dnds_summary$fdr <- p.adjust(dnds_summary$pvalue, method = "fdr")
  dnds_summary$stars <- add.significance.stars(dnds_summary$fdr)
  return(dnds_summary)
}

# This function creates start and end coordinates of clusters froma list of pvalues
cluster_matrix <- function(pvalues,
                           significance_level=0.05,
                           w_size=10){
  iteration <- 1 
  i <- 1
  clusters <- NULL
  
  while (!is.na (which(pvalues[i:length(pvalues)] <= significance_level)[1])) { # check if there are some significant elements
    if (i==1){
      i <- which(pvalues[i:length(pvalues)] <= significance_level)[1]  # first sig elem
      clusters$start[iteration] <- i
      i <- which(pvalues[i:length(pvalues)] > significance_level)[1] + clusters$start[iteration] - 1 # last sig. elem.
      clusters$end[iteration] <- i 
      iteration <- iteration + 1
    } else {
      i <- which(pvalues[i:length(pvalues)] <= significance_level)[1] + clusters$end[iteration-1] 
      clusters$start[iteration] <- i
      i <- which(pvalues[i:length(pvalues)] > significance_level)[1] + clusters$start[iteration] - 1
      clusters$end[iteration] <- i 
      iteration <- iteration + 1
    }
  }
  if(is.null(clusters)){ # add some zero values to prevent plotting issues
    clusters$start[1] <- 0
    clusters$end[1] <- 0
    clusters <- as.data.frame(clusters)
  } else {
    clusters <- as.data.frame(clusters)
    clusters$start <- clusters$start - w_size/2
    clusters$end <- clusters$end + w_size/2
  }
  return(clusters)
}

# extracts a legend from a plot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# loads dnds results and generate data frame
load_dnds <- function(in_folder = "dnds/",
                      in_colnames = c("pos", "Sd","Nd", "Sd_tip", "Nd_tip", "S", "N", "comparisons"),
                      gap_folder = "gap/",
                      dnds_ending = ".txt.DnDsRatio.txt",
                      gap_ending = ".gap.txt",
                      threshold = 1,
                      sum_threshold = 10,
                      quite= FALSE){
  require(reshape2)
  
  # check arguments
  if(!is.character(in_folder))
    stop("Error: input path is not a character")
  if(!is.character(gap_folder))
    stop("Error: gap path is not a character")
  if(!is.numeric(threshold))
    stop("Error: Threshold is not numeric")
  if(threshold>1 || threshold<0)
    stop("Error: Threshold must be between 0 and 1")
  
  # process arguments
  file_list <- list.files(in_folder)
  
  # create results matrix
  mat = data.frame(matrix(vector(), length(file_list), 4, 
                          dimnames=list(c(), c("name", "sum_pN", "sum_pS","ratio"))), stringsAsFactors=F) # name, sum(dN), sum(dS), sum(dN/dS) 
  
  # iterate over every file
  i <- 1
  for(file in file_list)
  {
    # exctract family identifier
    fam_name <- colsplit(string=file, pattern=dnds_ending, names=c("name", "ending"))$name
    if(!quite){
      #cat(paste("Load",fam_name,"\n"))
    }
    file_data <- read.table(paste(in_folder,file,sep=""), header=T)
    if (threshold < 1){
      # exclude dnds values based on gap proportion
      gap_data <- read.table(paste(gap_folder,fam_name,gap_ending,sep=""), header=F)
      file_data <- file_data[which(gap_data<threshold),]
    }
    
    # minimal pathway method with correction for multiple substitutions
    # take mean value for Sd and S based on number of comparisons
    file_data$pS <- (file_data$Sd + file_data$Sd_tip) / file_data$S
    file_data$pN <- (file_data$Nd + file_data$Nd_tip) / file_data$N
    
    # sum without NA and without Inf
    sum_pN <- sum(file_data[which(is.finite(file_data$pN)),]$pN, na.rm=T)
    sum_pS <- sum(file_data[which(is.finite(file_data$pS)),]$pS, na.rm=T)     
    
    if(is.finite(sum_pS)){ 
      if (sum_pS > sum_threshold){
        ratio <- sum_pN / sum_pS
      } else {
        ratio <- NA
      }
    }
    
    # add row to output matrix
    mat$name[i] <- as.character(fam_name)
    mat$sum_pN[i] <- sum_pN
    mat$sum_pS[i] <- sum_pS
    mat$ratio[i] <- ratio
    i <- i + 1   
  }
  return(mat)
}

# create msa plot
create_msa_plot <- function(dnds_path = "example.dnds.txt", 
                            gap_path = "example.gap.txt", 
                            only_legend = FALSE,
                            window_size = 10,
                            gap_threshold = 0.6){
  
  
  #### get dataset information ####
  #data_files <- basename(list.files(data_folder))
  # sample_name <-sub("^([^.]*).*", "\\1", data_files[1] )
  
  #### import data #### 
  gap_data <- try(read.table(gap_path,header=F))
  #if (inherits(gap_data, 'try-error')){  cat("No gap information provided. Please ad .gapdata \n")} 
  dnds_data <- try(read.table(dnds_path,header=T))
  #if (inherits(dnds_data, 'try-error')){  cat("No dnds information provided. Please ad .dndsdata \n")} 
  #js_data <- try(read.table(paste(data_folder, data_files[grep(".jsdata", data_files)] ,sep="")))
  #if (inherits(js_data, 'try-error')){  cat("No js information provided. Please ad .jsdata \n")} 
  #repeat_data <- try(read.table(paste(data_folder, data_files[grep(".repeatdata", data_files)],sep="")))
  #if (inherits(repeat_data, 'try-error')){  cat("No repeat information provided. Please ad .repeatdata \n")} 
  #tm_data <- try(read.table(paste(data_folder, data_files[grep(".tmdata", data_files)],sep=""), sep=";"))
  #if (inherits(tm_data, 'try-error')){  cat("No tm information provided. Please ad .tmdata \n")} 
  
  #### check data ####
  # check for length
  #print(head(gap_data))
  #print(head(dnds_data))
  
  
  # if (nrow(gap_data) != nrow(dnds_data)){ stop("Number of rows in dnds data is not equal number of rows of gap data")}  
  
  #### prepare data frame ####
  # process the dnds data
  dN <- dnds_data$Nd_tip + dnds_data$Nd 
  dS <- dnds_data$Sd_tip + dnds_data$Sd
  dNdS <- rep(0, length(dN))
  dNdS[which(dS > 0)] <- dN[which(dS > 0)]/dS[which(dS > 0)]
  
  # process the dnds data with sliding window
  pvalue_window <- sliding_window(dN, dS, gap_data =  gap_data, 
                                  w_size = window_size,
                                  g_threshold = gap_threshold)
  
  # generate a matrix with start and end positions for plotting of significant windows
  epitopes <- cluster_matrix(pvalues = pvalue_window, w_size = window_size, 
                             significance_level = 0.05)
  
  # create data frame with all informations
  df <- data.frame(position=dnds_data$pos,
                   selection=dNdS,
                   selection_smooth=smooth.spline(dnds_data$pos, dNdS, spar=0.35)$y,
                   gap=gap_data$V1,
                   gap_smooth=smooth.spline(dnds_data$pos, gap_data$V1, spar=0.35)$y,
                   stringsAsFactors=FALSE) 
  
  #if we have js data, uncomment this:
  #  df <- data.frame(position=dnds_data$pos,
  #                   selection=dNdS,
  #                   selection_smooth=smooth.spline(dnds_data$pos, dNdS, spar=0.35)$y,
  #                   gap=gap_data$V1,
  #                   gap_smooth=smooth.spline(dnds_data$pos, gap_data$V1, spar=0.35)$y,
  #                   js=js_data$V1, 
  #                   js_smooth=smooth.spline(dnds_data$pos, js_data$V1, spar=0.35)$y,
  #                   stringsAsFactors=FALSE) 
  #
  
  # process repeat informations if supplied 
  # generating a boolean matrix with TRUE if there is a repeat element on this position
  # if (!inherits(repeat_data, 'try-error')){ # if annotation is provided
  #   df$is_repeat <- rep(FALSE, length(dNdS))
  #   if (!is.null(repeat_data)){ # if there are repetitive elements
  #     for (segment in 1:nrow(repeat_data)){
  #       try(df$is_repeat[repeat_data$V1[segment]:repeat_data$V2[segment]] <- TRUE,TRUE)
  #     }
  #   }
  # } 
  
  #  # process membrane informations if supplied 
  #  if (!inherits(tm_data, 'try-error')){ # if annotation is provided
  #    df$is_inside_membrane <- rep(FALSE, length(dNdS))
  #    df$is_outside_membrane <- rep(FALSE, length(dNdS))
  #    df$is_tm_helix <- rep(FALSE, length(dNdS))
  #    if (!is.null(tm_data)){ # if there are repetitive elements
  #      for (segment in 1:nrow(tm_data)){
  #        #TODO: add tm infomation to data frame 
  #      }
  #    }
  #  } 
  
  #### genereate plot #### 
  # create sequence plot figure with dnds and gap informations
  fig_a <- ggplot(df, aes(position, selection))
  # highlight significant aereas
  fig_a <- fig_a + geom_rect(data=epitopes, aes(NULL, NULL, xmin = start, xmax = end), 
                             ymin = -Inf, ymax = Inf,  fill="grey80")
  
  fig_a <- fig_a + geom_ribbon(aes(x=position, ymax=selection_smooth,ymin=1), fill="grey60")
  fig_a <- fig_a + geom_point(aes(colour=gap), size = 1.7, alpha=3/4)
  fig_a <- fig_a + scale_colour_gradient(low = "green", high="red")
  fig_a <- fig_a + theme_bw() + xlab(" ") + ylab("dN/dS ratio") 
  fig_a <- fig_a + geom_line(aes(y=selection_smooth))
  #fig_a <- fig_a + labs(title=paste(sample_name,"; g_threshold=", gap_threshold,"; window_size=", window_size,sep=""))
  # fig_a <- fig_a + geom_hline(aes(yintercept=1))
  
  # plot repeat informations
  #  if (!inherits(repeat_data, 'try-error')){
  #    fig_a <- fig_a + geom_segment(data=repeat_data,aes(x = V1, y = -1.2, xend = V2, yend = -1.2), color = "grey")
  #  }
  
  # create sequence plot figure with JS informations
  # fig_b <- ggplot(df, aes(position, js_smooth)) 
  # fig_b <- fig_b + geom_line(aes(y=js_smooth))
  # fig_b <- fig_b + theme_bw() + xlab("Sequence Position (nt)") + ylab("JS score") 
  # fig_b <- fig_b + scale_y_continuous(breaks=c(0,1),limits=c(0,1)) 
  
  # plot membrane informations
  # if (!inherits(tm_data, 'try-error')){
  #   fig_b <- fig_b + geom_rect(data=tm_data,aes(NULL, NULL, xmin = V2, xmax = V3, fill=V1), ymin = -Inf, ymax = Inf,alpha=0.2)
  # }
  
  #  if (only_legend){
  #    fig_b_legend <- g_legend(fig_b)  
  #    fig_a_legend <- g_legend(fig_a)
  #   fig <- grid.arrange(fig_b_legend, fig_a_legend, nrow=2, heights=1:3)
  
  #  } else {
  #fig_b <- fig_b + theme(legend.position="none")
  #   fig_a <- fig_a + theme(legend.position="none")    
  #   fig <- grid.arrange(fig_a, fig_b, nrow=2, heights=c(3/4,1/4), ncol=1)    
  #  } 
  return(fig_a)
}
# create density plot
plot_density <- function(set){
  require("ggplot2")
  den_set <- as.data.frame(set)
  den_set$ratio <- as.numeric(as.matrix(set$ratio))
  d <- ggplot(den_set, aes(ratio))
  d <- d + geom_histogram(alpha=0.8,binwidth=.05) + theme_bw()
  d <- d + ggtitle("Histogram")
  d <- d + theme(plot.title = element_text(size=15, vjust=1, lineheight=0.6))
  d <- d + xlab("dN/dS ratio") +  ylab("number of protein families")
  return(d)  
}

# This function creates a p-value matrix indicating significant positions based
sliding_window <- function(dN, dS, gap_data, w_size=10, g_threshold=0.2){
  dNdS_window <- rep(0,length(dN))
  if (w_size>0){
    dN_window <- rollsum(dN, w_size, fill = list(NA, NULL, NA))
    dS_window <- rollsum(dS, w_size, fill = list(NA, NULL, NA))
    dNdS_window[which(dS_window > 0)] <- dN_window[which(dS_window > 0)] / 
      dS_window[which(dS_window > 0)]
    dNdS_gap <- rollmean(gap_data$V1, w_size, fill=list(NA,NULL,NA))
    dNdS_gap[which(is.na(dNdS_gap))] <- 1
  }
  
  # iterate over window and find significant windows with FDR correction
  dN_all <- sum(dN, na.rm=T)
  dS_all <- sum(dS, na.rm=T)
  window <- rep(1, length(dN))
  for (i in 1:length(dNdS_window)){
    if (dNdS_gap[i] < g_threshold){ 
      if (!is.na(dN_window[i]) && !is.na(dS_window[i])){
        test_matrix <-  matrix(c(dN_window[i], dN_all,
                                 dS_window[i], dS_all),nrow=2)
        pval <- fisher.test(test_matrix,alternative=c("greater"))[["p.value"]]
        window[i] <- p.adjust(pval, method="fdr",n = length(window))
      } else {
        window[i] <- 1
      }
    } else { # g_threshold
      window[i] <- 1
    }
  }
  return(window)
}
