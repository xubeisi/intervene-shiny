options(shiny.maxRequestSize=1000*1024^2)

#devtools::install_github("jrowen/rhandsontable")

source("pairwise_intersect.R")
library(RColorBrewer)
library(htmlwidgets)
library(gplots)
library(dendextend)
library(excelR)

myisna <- function(x){
  is.na(x) | x == '.' | x == ''
}

#' List of named vectors to UpSetR converter
#' 
#' @description A function to convert a list of named vectors to a data frame compatible with UpSetR.
#' @param input A list of named vectors to be converted to a data frame compatible with UpSetR
#' @note See "Basic Usage" vignette for an example on how to use this function in UpSetR.
#' @export 
fromList <- function(input){
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x){x <- as.vector(match(elements, x))}))
  data[myisna(data)] <- as.integer(0); data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) !=0), ]
  names(data) <- names(input)
  row.names(data) <- elements
  return(data)
}

Upset_comb <- function(data){
  zz <- data.frame(group=apply(data==1,1,function(a) paste0(colnames(data)[a], collapse = "&")))
  zz$gene <- row.names(zz)
  zz <- dcast(zz, gene~group)
  zz <- lapply(zz,function(x) zz[!is.na(x),1])
  zz <- data.frame(sapply(zz, "length<-", max(lengths(zz))),stringsAsFactors = FALSE, check.names = FALSE)
  row.names(zz) <- zz$gene
  zz$gene <- NULL
  zz
}

Counter <- function(data, empty_intersects = FALSE){
  temp_data <- list()
  Freqs <- data.frame()
  num_sets <- ncol(data)
  end_col <- as.numeric(num_sets)
  name_of_sets <- names(data)
  #gets indices of columns containing sets used
  for( i in 1:num_sets){
    temp_data[i] <- match(name_of_sets[i], colnames(data))
  }
  Freqs <- data.frame(count(data[ ,as.integer(temp_data)]))
  colnames(Freqs)[1:num_sets] <- name_of_sets
  #Adds on empty intersections if option is selected
  if(is.null(empty_intersects) == F){
    empty <- rep(list(c(0,1)), times = num_sets)
    empty <- data.frame(expand.grid(empty))
    colnames(empty) <- name_of_sets
    empty$freq <- 0
    all <- rbind(Freqs, empty)
    Freqs <- data.frame(all[!duplicated(all[1:num_sets]), ], check.names = F)
  }
  #Remove universal empty set
  Freqs <- Freqs[!(rowSums(Freqs[ ,1:num_sets]) == 0), ]
  #Aggregation by degree
  for(i in 1:nrow(Freqs)){
    Freqs$degree[i] <- rowSums(Freqs[ i ,1:num_sets])
  }
  Freqs <- Freqs[order(Freqs$freq,decreasing = TRUE),]
  for( i in 1:nrow(Freqs)){
    Freqs$x[i] <- i
  }

  nintersections = min(1000,nrow(Freqs))

  Freqs <- Freqs[1:nintersections, ]
  Freqs <- na.omit(Freqs)
  return(Freqs)
}

#sever code
shinyServer(function(input, output, session) {
  
  output$hot_venn = renderExcel({
    test <- venn_data()
    ncol <- length(test)
    test <- data.frame(sapply(test, "length<-", max(lengths(test))),stringsAsFactors = FALSE)
    #old_able_rename test <- rbind(names(test),test)
    #old_able_rename names(test) <- paste("V",1:ncol,sep="")
    excelTable(data = test,wordWrap=TRUE,search=TRUE,loadingSpin=TRUE,showToolbar=TRUE,autoWidth=TRUE,defaultColWidth=10)
  })

  output$hot_upset = renderExcel({
    test <- upset_data()
    test <- lapply(test, function(x) as.character(row.names(test)[x>0])) 
    ncol <- length(test)
    test <- data.frame(sapply(test, "length<-", max(lengths(test))),stringsAsFactors = FALSE)
    #old_able_rename test <- rbind(names(test),test)
    #old_able_rename names(test) <- paste("V",1:ncol,sep="")
    excelTable(data=test,wordWrap=TRUE,search=TRUE,loadingSpin=TRUE,showToolbar=TRUE,autoWidth=TRUE,defaultColWidth=10)
  })
  
  # output$hot_venn = renderRHandsontable({
  #   rhandsontable(test)
  # })  
  
  #====================================================#
  ## Venn module ####
  #====================================================#
  venn_type <- reactive({
    return(input$venn_type)
  })
  
  doWeights <- reactive({
    return(input$doWeights)
  })
  
  doEuler <- reactive({
    return(input$doEuler)
  })
  
  venn_size <- reactive({
    return(input$venn_size)
  })
  
  venn_lwd <- reactive({
    return(as.numeric(input$venn_lwd))
  })
  
  venn_labelsize <- reactive({
    return(as.numeric(input$venn_labelsize))
  })
  
  venn_color_type <- reactive({
    return(input$venn_color_type)
  })
  
  
  venn_cex <- reactive({
    return(as.numeric(input$venn_cex))
  })
  
  venn_lty <- reactive({
    return(as.numeric(input$venn_lty))
  })
  
  set1_color <- reactive({
    return(input$set1_color)
  })
  set2_color <- reactive({
    return(input$set2_color)
  })
  set3_color <- reactive({
    return(input$set3_color)
  })
  set4_color <- reactive({
    return(input$set4_color)
  })
  set5_color <- reactive({
    return(input$set5_color)
  })
  set6_color <- reactive({
    return(input$set6_color)
  })
  
  venn_data_excel <- eventReactive(input$update_tbl_venn,{
    input_type <- input$venn_input_type

    if (!is.null(input$hot_venn)) {
      #oldbeisi data = hot_to_r(input$hot_venn)
      data = excel_to_R(input$hot_venn)
      #old_able_rename names(data) <- data[1,,drop=F]
      #old_able_rename data <- data[-1,]
      if (input_type == 'list'){
        data <- lapply(data, function(x) x[!myisna(x)])
      } else if (input_type == 'binary'){
        data <- lapply(data, function(x) as.character(data[,1][x>0]))
        data[[1]] <- NULL
      }
      data
    } else {
      list()
    }
  },ignoreNULL = FALSE
  )
  
  venn_data <- reactive({
    inFile <- input$file_venn
    string <- input$venn_comb
    string <- gsub("\n", "", string)
    input_type <- input$venn_input_type
    if(string != ""){
      string <- as.list(unlist(strsplit(string, ",")))
      names <- lapply(string, function(x){x <- unlist(strsplit(x, "=")); x <- x[1]})
      names <- unlist(lapply(names, function(x){x <- gsub(" ", "", x)}))
      values <- as.numeric(unlist(lapply(string, function(x){x <- unlist(strsplit(x,"=")); x <- x[2]})))
      names(values) <- names
      venneuler <- fromExpression(values)
      data <- lapply(venneuler, function(x) as.character(row.names(venneuler)[x>0])) # df to list
    } else if(is.null(inFile) == F){
      if (input_type == 'list'){
        data <- read_delim(input$file_venn$datapath, input$sep_venn , escape_double = FALSE, trim_ws = TRUE, col_names = input$header_venn)
        data <- lapply(data, function(x) x[!myisna(x)])
      } else if (input_type == 'binary'){
        data <- read.csv(input$file_venn$datapath, header = input$header_venn,
                         sep = input$sep_venn, quote = input$quote)
        data <- lapply(data, function(x) as.character(data[,1][x>0]))
        data[[1]] <- NULL
      }
    } else if (length(venn_data_excel())){
      data <- venn_data_excel()
    }else{
      data <- read_delim('data/Whyte_et_al_2013_SEs_genes.csv', ",", escape_double = FALSE, trim_ws = TRUE, col_names = TRUE)
      data <- lapply(data, function(x) x[!myisna(x)])
    }
    if (input$sep_venn_row != "")
    {
      data <- lapply(data, function(x) unique(unlist(strsplit(x,input$sep_venn_row))))
    }
    
    return(data)
  })
  
  set_names <- reactive({
    names <- names(venn_data())
    return(names)
  })
  
  output$venn_sets <- renderUI({
    venn_sets <- selectInput('venn_sets', label = "Select sets",
                             choices = as.character(set_names()),
                             multiple = T, selectize = T, selected = as.character(set_names()[1:3]))
    return(venn_sets)
  })
  
  venn_selected_names <- reactive({
    venn_selected_names <- as.character(c(input$venn_sets))
  })
  
  venn_data_filtered <- reactive({
    
    data <- venn_data()
    if(is.null(input$venn_sets)){
      data <- data[names(data)[1:3]]
      return(data)
    }else{
      data <- data[c(venn_selected_names())]
      return(data)
    }
    return(data)
  })
  
  venn_combinations <- reactive({
    #string <- input$venn_comb
    string <- ""
    data <- venn_data_filtered()
    if (string !=""){
      return(data)
    } else {
      return(Venn(data))
    }
  })
  
  get_venn_gp <- reactive({
    venn_gp <- VennThemes(compute.Venn(venn_combinations()))
    venn_gp$SetText <- lapply(venn_gp$SetText,function(x) {x$fontsize<-venn_labelsize(); return(x)})
    venn_gp$FaceText <- lapply(venn_gp$FaceText,function(x) {x$cex<-venn_cex(); return(x)})
    venn_gp$Set <- lapply(venn_gp$Set,function(x) {x$lwd<-venn_lwd(); return(x)})
    venn_gp$Set <- lapply(venn_gp$Set,function(x) {x$lty<-venn_lty(); return(x)})
    
    if (venn_color_type () == 'custom'){
      venn_gp$Set$Set1$col <- set1_color()
      venn_gp$Set$Set2$col <- set2_color()
      venn_gp$Set$Set3$col <- set3_color()
      venn_gp$Set$Set4$col <- set4_color()
      venn_gp$Set$Set5$col <- set5_color()
      venn_gp$Set$Set6$col <- set6_color()
      
      venn_gp$SetText$Set1$col <- set1_color()
      venn_gp$SetText$Set2$col <- set2_color()
      venn_gp$SetText$Set3$col <- set3_color()
      venn_gp$SetText$Set4$col <- set4_color()
      venn_gp$SetText$Set5$col <- set5_color()
      venn_gp$SetText$Set6$col <- set6_color()
      
    }
    
    return(venn_gp)
  })
  
  data_size <- reactive({
    return(length(venn_data_filtered()))
  })
  
  get_venn_type <- reactive({
    if (venn_type() == 'Classical'){
      if(data_size() < 4)
        return("circles")
      else
        return("ellipses")
    }else if (venn_type() == 'ChowRuskey' && data_size() < 3){
      return("circles")
    }
    else{
      return(venn_type())
    }
  })
  
  output$vennPlot <- renderPlot({
    plot(compute.Venn(venn_combinations(), doWeights = doWeights(), doEuler = doEuler(), type = get_venn_type()),
         gp = get_venn_gp(),
         show = list(Universe = FALSE)
    )
  },
  width = venn_size,
  height = venn_size,
  outputArgs = list()
  )
  
  output$VennDown <- downloadHandler(
    filename = function(){
      paste("Venn_diagram", tolower(input$filetype_venn), sep =".")
    }, 
    content = function(file){
      width  <- venn_size()
      height <- venn_size()
      #width  <- session$clientData$output_plot_width
      #height <- ((session$clientData$output_plot_height)*1)
      #pixelratio <- session$clientData$pixelratio
      pixelratio <- 2
      
      if(input$filetype_venn == "PNG")
        png(file, width=width*pixelratio, height=height*pixelratio, units = "px", res=72*pixelratio)
      else if(input$filetype_venn == "SVG")
        svg(file, width=8, height=8)
      else if(input$filetype_venn == "TIFF")
        tiff(file, width=width*pixelratio, height=height*pixelratio, units = "px")
      else
        pdf(file, width = 8, height = 8)
      
      plot(venn_combinations(),
           doWeights = doWeights(),
           type = get_venn_type(),
           doEuler = doEuler(),
           show = list(Universe = FALSE)
           #venn_type()
      )
      dev.off()
    }
  )
  
  output$vennDownExcel <- downloadHandler(
    filename = function(){
      paste("Venn_diagram", tolower(input$filetype_venn_excel), "csv", sep =".")
    }, 
    content = function(file){
      data <- fromList(venn_data_filtered())
      if(input$filetype_venn_excel == "Freq"){
        data <- Counter(data)
      } else if(input$filetype_venn_excel == "Combinations"){
        data <- Upset_comb(data)
      }
      write.table(data.frame("Row"=row.names(data),data,check.names = FALSE),file,na = "",sep=",",row.names = FALSE)
    }
  )
  
  #====================================================#
  ## UpSet module ####
  #====================================================#
  #Some of the code for upset module is taken from
  #https://github.com/hms-dbmi/UpSetR-shiny
  
  output$plot_text <- renderUI({
    if(is.null(upset_data()) == T){
      h5("There is no data entered. Please upload your data to draw UpSet plot here!")
    }
    else{
      HTML(" ")
    }
  })
  
  upset_data_excel <- eventReactive(input$update_tbl_upset,{
    input_type <- input$upset_input_type
    if (!is.null(input$hot_upset)) {
      #oldbeisi data = hot_to_r(input$hot_upset)
      data = excel_to_R(input$hot_upset)
      #old_able_rename names(data) <- data[1,,drop=F]
      #old_able_rename data <- data[-1,]
      if (input_type == 'list'){
        data <- lapply(data, function(x) x[!myisna(x)])
      } else if (input_type == 'binary'){
        data <- lapply(data, function(x) as.character(data[,1][x>0]))
        data[[1]] <- NULL
      }
      fromList(data)
    } else {
      data.frame()
    }
  },ignoreNULL = FALSE
  )
  
  upset_data <- reactive({  
    inFile <- input$file_upset
    string <- input$comb_upset
    string <- gsub("\n", "", string)
    input_type <- input$upset_input_type

    if(string != ""){
      string <- as.list(unlist(strsplit(string, ",")))
      names <- lapply(string, function(x){x <- unlist(strsplit(x, "=")); x <- x[1]})
      names <- unlist(lapply(names, function(x){x <- gsub(" ", "", x)}))
      values <- as.numeric(unlist(lapply(string, function(x){x <- unlist(strsplit(x,"=")); x <- x[2]})))
      names(values) <- names
      data <- fromExpression(values)
    } else if(is.null(inFile) == F){
      if (input_type == 'list'){
        data <- read_delim(inFile$datapath, input$sep_upset , escape_double = FALSE, trim_ws = TRUE, col_names = input$header_upset)
        if (input$sep_row_upset != "")
        {
          data <- lapply(data, function(x) unique(unlist(strsplit(x,input$sep_row_upset))))
        }
        data <- fromList(lapply(as.list(data), function(x) x[!myisna(x)]))
      } else if (input_type == 'binary'){
        data <- read.csv(inFile$datapath, header = input$header_upset,
                         sep = input$sep_upset, quote = input$quote)
      }
    } else if (nrow(upset_data_excel())) {
      data <- upset_data_excel()
    }else{
      data<- fromExpression(c('H3K4me2&H3K4me3'=1632,'H3K4me2&H3K4me3&H3K27me3'=575,'H3K27me3'=2517,'H3K4me3&H3K27me3'=1553,'H3K4me3'=3296,'H3K4me2&H3K27me3'=1903,'H3K4me2'=6029,'H3K27ac&H3K4me2&H3K4me3&H3K27me3'=723,'H3K27ac&H3K4me2&H3K4me3'=1750,'H3K27ac&H3K4me2'=2134,'H3K27ac&H3K4me2&H3K27me3'=169,'H3K27ac&H3K4me3'=813,'H3K27ac&H3K4me3&H3K27me3'=29,'H3K27ac&H3K27me3'=760,'H3K27ac'=4216))
    }
    return(data)
  })
  
  FindStartEnd <- function(data){
    startend <- c()
    for(i in 1:ncol(data)){
      column <- data[, i]
      column <- (levels(factor(column)))
      if((column[1] == "0") && (column[2] == "1" && (length(column) == 2))){
        startend[1] <- i
        break
      }
      else{
        next
      }
    }
    for(i in ncol(data):1){
      column <- data[ ,i]
      column <- (levels(factor(column)))
      if((column[1] == "0") && (column[2] == "1") && (length(column) == 2)){
        startend[2] <- i
        break
      }
      else{
        next
      }
    }
    return(startend)
  }
  
  startEnd <- reactive({
    startEnd <- FindStartEnd(upset_data())
  })
  
  setSizes <- reactive({
    if(is.null(upset_data()) != T){
      sizes <- colSums(upset_data()[startEnd()[1]:startEnd()[2]])
      sizes <- sizes[order(sizes, decreasing = T)]
      names <- names(sizes); sizes <- as.numeric(sizes);
      maxchar <- max(nchar(names))
      total <- list()
      for(i in 1:length(names)){
        spaces <- as.integer((maxchar - nchar(names[i]))+1)
        spaces <- paste(rep(" ", each=spaces), collapse = "")
        total[[i]] <- paste(paste(names[i], ":", sep=""), spaces, sizes[i], "\n", sep="")
      }
      total <- unlist(total)
      total <- paste(total, collapse = " ")
      return(total)
    }
    else{
      return(NULL)
    }
  })
  
  output$setsizes <- renderText({
    if(is.null(setSizes()) != T){
      paste("---Set Sizes---\n", setSizes())
    }
    else{
      paste("---Set Sizes---\n", "\n No Data Entered")
    }
  })
  
  Specific_sets <- reactive({
    Specific_sets <- as.character(c(input$upset_sets))
  })
  
  output$sets <- renderUI({
    if(is.null(upset_data()) == T){
      sets <-  selectInput('upset_sets', label="Select at least two sets ",
                           choices = NULL,
                           multiple=TRUE, selectize=TRUE, selected = Specific_sets())
    }
    else{
      data <- upset_data()[startEnd()[1]:startEnd()[2]]
      topfive <- colSums(data)
      topfive <- as.character(head(names(topfive[order(topfive, decreasing = T)]), 5))
      sets <- selectInput('upset_sets', label="Select sets ",
                          choices = as.character(colnames(upset_data()[ , startEnd()[1]:startEnd()[2]])),
                          multiple=TRUE, selectize=TRUE, selected = topfive)
    }
    return(sets)
  })
  
  
  mat_prop <- reactive({
    mat_prop <- input$mbratio
  })
  upset_width <- reactive({
    return(input$upset_width)
  })
  upset_height <- reactive({
    return(input$upset_height)
  })
  
  bar_prop <- reactive({
    bar_prop <- (1 - input$mbratio)
  })
  
  orderdat <- reactive({
    orderdat <- as.character(input$order)
    if(orderdat == "degree"){
      orderdat <- c("degree")
    }
    else if(orderdat == "freq"){
      orderdat <- "freq"
    }
    return(orderdat)
  })
  
  show_numbers <- reactive({
    show_numbers <- input$show_numbers
    if(show_numbers){
      show_numbers <- "yes"
      return(show_numbers)
    }
    else{
      show_numbers <- FALSE
      return(show_numbers)
    }
    
  })
  
  main_bar_color <- reactive({
    mbcolor <- input$mbcolor
    return(mbcolor)
  })
  sets_bar_color <- reactive({
    sbcolor <- input$sbcolor
    return(sbcolor)
  })
  
  
  decrease <- reactive({
    decrease <- as.character(input$decreasing)
    if(decrease == "inc"){
      decrease <- FALSE
    }
    else if(decrease == "dec"){
      decrease <- TRUE
    }
    return(decrease)
  })
  
  number_angle <- reactive({
    angle <- input$angle
    return(angle)
  })
  
  line_size <- reactive({
    line_size <- input$linesize
    return(line_size)
  })
  
  emptyIntersects <- reactive({
    if(isTRUE(input$empty)){choice <- "on"
    return(choice)
    }
    else{
      return(NULL)
    }
  })
  
  scale.intersections <- reactive({
    return(input$scale.intersections)
  })
  
  scale.sets <- reactive({
    return(input$scale.sets)
  })
  
  keep.order <- reactive({
    return(input$keep.order)
  })
  
  # A plot of fixed size
  output$plot1 <- renderPlot({
    
    if(length(upset_data()) == 0){stop()}
    if(length(Specific_sets()) == 1){
      stop()
    }
    upset(data = upset_data(), 
          nintersects = input$nintersections,
          point.size = input$pointsize,
          line.size = line_size(),
          sets = Specific_sets(),
          order.by = orderdat(),
          main.bar.color= main_bar_color(),
          sets.bar.color= sets_bar_color(),
          decreasing = c(decrease()),
          show.numbers = show_numbers(),
          number.angles = number_angle(),
          scale.intersections = scale.intersections(),
          scale.sets = scale.sets(),
          keep.order = keep.order(),
          mb.ratio = c(as.double(bar_prop()), as.double(mat_prop())),
          empty.intersections = emptyIntersects(),
          text.scale = c(input$intersection_title_scale, input$intersection_ticks_scale,
                         input$set_title_scale, input$set_ticks_scale, input$names_scale,
                         input$intersection_size_numbers_scale))},
    #width  <- session$clientData$output_plot_width
    #height <- ((session$clientData$output_plot_height)*1.7)
    width = upset_width,
    height = upset_height
  )
  
  #outputOptions(output, "plot", suspendWhenHidden = FALSE)
  
  # observe({
  #   if(pushed$B != 0 && length(pushed$B) == 1){
  #     updateTabsetPanel(session, "main_panel", "upset_plot")
  #   }
  # })
  
  output$UpSetDown <- downloadHandler(
    
    filename = function(){
      paste("UpSet_plot", tolower(input$filetype), sep =".")
    }, 
    content = function(file){
      width <- upset_width()
      height <- upset_height()
      pixelratio <- 2
      
      if(input$filetype == "PNG")
        png(file, width=width*pixelratio, height=height*pixelratio,
            res=72*pixelratio, units = "px")
      else if(input$filetype == "SVG")
        svg(file, width = width/100, height = height/100)
      else if(input$filetype == "TIFF")
        tiff(file, width=width*pixelratio, height=height*pixelratio, units = "px")
      else
        pdf(file, width = width/100, height = height/100, onefile=FALSE)
      
      print(upset(data = upset_data(), 
                  nintersects = input$nintersections,
                  point.size = input$pointsize,
                  line.size = line_size(),
                  
                  sets = Specific_sets(),
                  order.by = orderdat(),
                  main.bar.color= main_bar_color(),
                  sets.bar.color= sets_bar_color(),
                  decreasing = c(decrease()),
                  number.angles = number_angle(),
                  show.numbers = show_numbers(),
                  scale.intersections = scale.intersections(),
                  scale.sets = scale.sets(),
                  keep.order = keep.order(),
                  mb.ratio = c(as.double(bar_prop()), as.double(mat_prop())),
                  empty.intersections = emptyIntersects(),
                  text.scale = c(input$intersection_title_scale, input$intersection_ticks_scale,
                                 input$set_title_scale, input$set_ticks_scale, input$names_scale,
                                 input$intersection_size_numbers_scale))
      )
      
      dev.off()
    }
  )
  
  output$upsetDownExcel <- downloadHandler(
    filename = function(){
      paste("UpSet", tolower(input$filetype_upset_excel), "csv", sep =".")
    }, 
    content = function(file){
      data <- upset_data()
      #browser()
      if(input$filetype_upset_excel == "Freq"){
        data <- Counter(data)
      } else if(input$filetype_upset_excel == "Combinations"){
        data <- Upset_comb(data)
      }
      write.table(data.frame("Row"=row.names(data),data,check.names = FALSE),file,na = "",sep=",",row.names = FALSE)
    }
  )
  
  #====================================================#
  ## Pairwise module ####
  #====================================================#
  output$plot_text_p <- renderUI({
    if(is.null(pairwiseMatrix()) == T){
      h5("There is no data entered. Please upload your data to draw pairwise heatmap here!")
    }
    else{
      HTML(" ")
    }
  })
  
  corplot_method <- reactive({
    return(input$corp_method)
  })
  
  corplot_type <- reactive({
    return(input$corp_type)
  })
  
  corplot_order <- reactive({
    return(input$corp_order)
  })
  
  corplot_diag <- reactive({
    return(input$corp_diag)
  })
  
  corplot_tl.col <- reactive({
    return(input$tl_col)
  })
  
  corplot_title <- reactive({
    return(input$corp_title)
  })
  
  heatmap_size <- reactive({
    return(input$heatmap_size)
  })
  
  addrect <- reactive({
    return(input$addrect)
  })
  
  rect_col <- reactive({
    return(input$rect_col)
  })
  
  hclust_method <- reactive({
    return(input$hclust_method)
  })
  
  
  tl_pos <- reactive({
    if (corplot_type() =="lower" && input$tl_pos != 'n')
    {
      return('ld')
      
    }else if(corplot_type() =="upper" && input$tl_pos != 'n'){
      return('td')
    }else{
      return(input$tl_pos)
    }
  })
  
  cl_pos <- reactive({
    if (corplot_type() =="lower" && input$cl.pos != 'n')
    {
      return('b')
      
    }else if(corplot_type() =="upper" && input$cl.pos != 'n'){
      return('r')
    }else{
      return(input$cl.pos)
    }
  })
  
  addgrid_col <- reactive({
    return(input$addgrid_col)
  })
  
  tl_srt <- reactive({
    return(input$tl.srt)
  })
  
  tl_cex <- reactive({
    return(input$tl.cex)
  })
  
  cl_cex <- reactive({
    return(input$cl.cex)
  })
  
  lower_colour <- reactive({
    return(input$lower_colour)
  })
  middle_colour <- reactive({
    return(input$middle_colour)
  })
  higher_colour <- reactive({
    return(input$higher_colour)
  })
  
  heamap_colors <- reactive({
    if(input$color_type == 'custom'){
      colors = colorRampPalette(c(lower_colour(), middle_colour(), higher_colour()))(100)
    }else{
      colors = colorRampPalette(brewer.pal(9,input$color_type))(100)
    }
    return(colors)
  })
  
  
  
  dendrogram <- reactive({
    return(input$dendrogram)
  })
  symm <- reactive({
    return(input$symm)
  })
  key <- reactive({
    return(input$key)
  })
  keysize <- reactive({
    return(input$keysize)
  })
  key.title <- reactive({
    return(input$key.title)
  })
  key.xlab <- reactive({
    return(input$key.xlab)
  })
  key.ylab <- reactive({
    return(input$key.ylab)
  })
  
  distance <- reactive({
    if(input$distance == 'none'){
      distance= as.dist(pairwiseMatrix())    
    }else{
      distance= dist(pairwiseMatrix(), method =input$distance)    
    }
    return(distance)
  })
  
  is_correlation<- reactive({
    isCor <- input$corp_cor
    if(isCor == 'non'){
      return(FALSE)
    }else
    {
      return(TRUE)
    }
  })
  
  pairwiseMatrix <- reactive({
    inFile <- input$file_p
    isCor <- input$corp_cor
    input_type <- input$pairwise_input_type
    if (is.null(inFile)){
      
      myMatrix <- as.matrix(read.table("data/frac_pairwise_matrix.txt"))
      if(isCor != 'non'){
        myMatrix <- cor(myMatrix, method=isCor)
      }
      return (myMatrix)
    }else{
      if(input_type == 'matrix'){
        
        myMatrix <- as.matrix(read.table(inFile$datapath))
        
      }else{
        #myMatrix <- as.matrix(pairwise_intersect(read.csv(inFile$datapath, header = input$header_p, sep=input$sep_p)))
        myMatrix <- as.matrix(pairwise_intersect(lapply(as.list(read_delim(inFile$datapath, input$sep_p , escape_double = FALSE, trim_ws = TRUE, col_names = input$header_p)), function(x) x[!myisna(x)])))
        #myMatrix <- as.matrix(pairwise_intersect(lapply(as.list(read.csv(inFile$datapath, header = input$header_p, sep=input$sep_p)), function(x) x[!myisna(x)])))
      }
      
      if(isCor != 'non'){
        myMatrix <- cor(myMatrix, method=isCor)
      }
      return(myMatrix)
    }
  })
  
  min_limit <- reactive({
    isCor <- input$corp_cor
    if(isCor == 'non'){
      return(0)
    }else
    {
      return(-1)
    }
  })
  
  max_limit <- reactive({
    isCor <- input$corp_cor
    input_type <- input$pairwise_input_type
    
    if(isCor == 'non' && input_type == 'list'){
      return(as.integer(max(pairwiseMatrix())))
    }else
    {
      return(1)
    }
  })
  
  output$pairwiseTable = DT::renderDataTable(
    round(pairwiseMatrix(),2), options = list(
      lengthChange = TRUE
    )
  )
  
  heatmap2_plot <- reactive({
    hcluster <- hclust(distance(), method =hclust_method())
    dend1 <- as.dendrogram(hcluster)
    # get some colors
    dend1 <- color_branches(dend1, k = addrect())
    col_labels <- get_leaves_branches_col(dend1)
    # order of the data!    
    col_labels <- col_labels[order(order.dendrogram(dend1))]
    
    plt <- heatmap.2(pairwiseMatrix(),
                     scale = "none",
                     #dendrogram = "both",
                     col = heamap_colors(),
                     cexRow = tl_cex(),
                     cexCol = tl_cex(),
                     #srtRow = tl_srt(),
                     srtCol = tl_srt(),
                     Rowv = dend1, 
                     Colv = dend1, 
                     main = corplot_title(),
                     dendrogram = dendrogram(),
                     symm = symm(),
                     #revC = TRUE,
                     key = key(),
                     keysize = keysize(),
                     key.title = key.title(),
                     key.xlab =  key.xlab(),
                     key.ylab = key.ylab(),
                     key.par = list(cex=cl_cex()),
                     sepwidth = c(0.05, 0.05),  # width of the borders
                     mar=c(6,6),
                     sepcolor = addgrid_col(),
                     colsep =1:ncol(pairwiseMatrix()),
                     rowsep =1:nrow(pairwiseMatrix()),
                     #offsetRow = 0.1,
                     #offsetCol = 0.1,
                     trace="none",
                     #RowSideColors = col_labels, #colored strips        
                     colRow = col_labels,
                     ColSideColors = col_labels, #colored strips        
                     colCol = col_labels
    )
    return(plt)
  })
  
  output$heatmap2_plot_out <- renderPlot(
    heatmap2_plot(),
    
    width= heatmap_size,
    height= heatmap_size
  )
  
  output$Heatmap2PlotDown <- downloadHandler(
    filename = function(){
      paste("Pairwise_heatmap2", tolower(input$filetype_heatmap), sep =".")
    }, 
    content = function(file){
      width  <- heatmap_size()
      height <- heatmap_size()
      #width  <- session$clientData$output_plot_width
      #height <- ((session$clientData$output_plot_height)*2)
      pixelratio <- 2
      if(input$filetype_heatmap == "PNG")
        png(file, width=width*pixelratio, height=height*pixelratio,
            res=72*pixelratio, units = "px")
      else if(input$filetype_heatmap == "SVG")
        svg(file, width=12, height=12)
      else if(input$filetype_heatmap == "TIFF")
        tiff(file, width=width*pixelratio, height=height*pixelratio, units = "px")
      else
        pdf(file, width = 12, height = 12)
      
      hcluster <- hclust(distance(), method =hclust_method())
      dend1 <- as.dendrogram(hcluster)
      # get some colors
      dend1 <- color_branches(dend1, k = addrect())
      col_labels <- get_leaves_branches_col(dend1)
      # order of the data!    
      col_labels <- col_labels[order(order.dendrogram(dend1))]
      
      heatmap.2(pairwiseMatrix(),
                scale = "none",
                #dendrogram = "both",
                col = heamap_colors(),
                cexRow = tl_cex(),
                cexCol = tl_cex(),
                #srtRow = tl_srt(),
                srtCol = tl_srt(),
                Rowv = dend1, 
                Colv = dend1, 
                main = corplot_title(),
                dendrogram = dendrogram(),
                symm = symm(),
                #revC = TRUE,
                key = key(),
                keysize = keysize(),
                key.title = key.title(),
                key.xlab =  key.xlab(),
                key.ylab = key.ylab(),
                key.par = list(cex=cl_cex()),
                sepwidth = c(0.05, 0.05),  # width of the borders
                mar=c(6,6),
                sepcolor = addgrid_col(),
                colsep =1:ncol(pairwiseMatrix()),
                rowsep =1:nrow(pairwiseMatrix()),
                #offsetRow = 0.1,
                #offsetCol = 0.1,
                trace="none",
                #RowSideColors = col_labels, #colored strips        
                colRow = col_labels,
                ColSideColors = col_labels, #colored strips        
                colCol = col_labels
      )
      dev.off()
    }
  )
  
  d3HM_plot <- reactive({
    hcluster = hclust(distance(), method =hclust_method())
    dend1 <- as.dendrogram(hcluster)
    
    # get some colors
    #cols_branches <- c("darkred", "forestgreen", "orange", "blue")
    #cols_branches <- brewer.pal(addrect(), "Set1")
    
    # Set the colors of 4 branches
    #dend1 <- color_branches(dend1, k = addrect(), col = cols_branches[1:addrect()])
    dend1 <- color_branches(dend1, k = addrect())
    
    col_labels <- get_leaves_branches_col(dend1)
    # But due to the way heatmap.2 works - we need to fix it to be in the 
    # order of the data!    
    col_labels <- col_labels[order(order.dendrogram(dend1))]
    d3heatmap(pairwiseMatrix(),
              show_grid = TRUE,
              scale = "none",
              dendrogram = dendrogram(),
              anim_duration = 0,
              k_row = addrect(),
              k_col = addrect(),
              Rowv = dend1, 
              Colv = dend1,
              symm = symm(),
              #revC = FALSE,
              #hclustfun = function(x) hclust(x,method = hclust_method()),
              #distfun = function(x) dist(x,method = hclust_method()),
              colors = heamap_colors(),
              xaxis_font_size = "12px",
              yaxis_font_size = "12px"
              #xaxis_height = 8,
              #yaxis_width = 8
    )
  })
  
  output$d3HM <- renderD3heatmap(d3HM_plot())
  
  output$HeatmapHTMLDown <- downloadHandler(
    filename = function(){
      paste("Interactive_pairwise_heatmap", "html", sep =".")
    }, 
    content = function(file){
      saveWidget(d3HM_plot(), file)
    }
  )
  
  output$corrplotHM <- renderPlot({
    corrplot(pairwiseMatrix(),
             method = corplot_method(),
             title = corplot_title(),
             tl.col= corplot_tl.col(),
             cl.lim=c(min_limit(),max_limit()),
             is.corr = is_correlation(),
             diag = corplot_diag(),
             order = corplot_order(),
             hclust.method = hclust_method(),
             type = corplot_type(),
             addrect = addrect(),
             tl.pos = tl_pos(),
             cl.pos = cl_pos(),
             rect.col = rect_col(),
             addgrid.col= addgrid_col(),
             tl.cex = tl_cex(),
             cl.cex = cl_cex(),
             tl.srt = tl_srt(),
             mar=c(0,0,2,2),
             col = heamap_colors()
             #addCoef.col = "red",
             #col = col1(300)
             #col = list(color = brewer.pal(20, "RdBu"))
             
    )},
    width= heatmap_size,
    height= heatmap_size
  )
  
  output$HeatmapCSVDown <- downloadHandler(
    filename = function(){
      paste("Pairwise_matrix", "csv", sep =".")
    }, 
    content = function(file){
      write.csv(pairwiseMatrix(), file)
    }
  )
  
  
  output$HeatmapDown <- downloadHandler(
    filename = function(){
      paste("Pairwise_heatmap", tolower(input$filetype_heatmap), sep =".")
    }, 
    content = function(file){
      width  <- heatmap_size()
      height <- heatmap_size()
      #width  <- session$clientData$output_plot_width
      #height <- ((session$clientData$output_plot_height)*2)
      pixelratio <- 2
      if(input$filetype_heatmap == "PNG")
        png(file, width=width*pixelratio, height=height*pixelratio,
            res=72*pixelratio, units = "px")
      else if(input$filetype_heatmap == "SVG")
        svg(file, width=12, height=12)
      else if(input$filetype_heatmap == "TIFF")
        tiff(file, width=width*pixelratio, height=height*pixelratio, units = "px")
      else
        pdf(file, width = 12, height = 12)
      
      corrplot(pairwiseMatrix(),
               method = corplot_method(),
               title = corplot_title(),
               tl.col= corplot_tl.col(),
               cl.lim=c(min_limit(),max_limit()),
               is.corr = is_correlation(),
               diag = corplot_diag(),
               order = corplot_order(),
               hclust.method = hclust_method(),
               type = corplot_type(),
               addrect = addrect(),
               tl.pos = tl_pos(),
               cl.pos = cl_pos(),
               rect.col = rect_col(),
               addgrid.col= addgrid_col(),
               tl.cex = tl_cex(),
               cl.cex = cl_cex(),
               tl.srt = tl_srt(),
               mar=c(0,0,2,1),
               #addCoef.col = "red",
               col = heamap_colors()
               #col = col3(100)
      )
      dev.off()
    }
  )
  
}
)
