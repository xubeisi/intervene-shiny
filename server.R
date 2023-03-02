options(shiny.maxRequestSize=1000*1024^2)

options(java.parameters = "-Xss2048k")
libs <- "excelR,graph,RColorBrewer,htmlwidgets,gplots,dendextend,shiny,shinydashboard,DT,plotly,ggplot2,gridExtra,plyr,UpSetR,colourpicker,corrplot,BBmisc,readr,excelR,reshape2,RBGL,caTools,remotes,d3heatmap,Vennerable"
libs <- unlist(strsplit(libs,","))
req<-unlist( lapply(libs,function(p) suppressPackageStartupMessages(require(p,character.only=TRUE)) ) )
need<-libs[req==FALSE]

base::print(need)

if (0){
  if ("pak" %in% need){
    install.packages("pak", repos = sprintf("https://r-lib.github.io/p/pak/stable/%s/%s/%s", .Platform$pkgType, R.Version()$os, R.Version()$arch))
  }
  options(install.packages.check.source = "no")

  for (lib in need) { pak::pkg_install(lib,ask = FALSE) }
  
  if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
    
  }
  BiocManager::install("BiocGenerics")
  if ("d3heatmap" %in% need){
    remotes::install_github("talgalili/d3heatmap",dep = FALSE)
  }
  if ("Vennerable" %in% need){
    remotes::install_github("js229/Vennerable",dep = FALSE)
  }
  
  library(excelR)
  library(graph)
  library(RColorBrewer)
  library(htmlwidgets)
  library(gplots)
  library(dendextend)
  library(shiny)
  library(shinydashboard)
  library(DT)
  #library(d3heatmap)
  library(plotly)
  library(ggplot2)
  library(gridExtra)
  library(plyr)
  library(UpSetR)
  library(colourpicker)
  library(corrplot)
  library(BBmisc)
  library(readr)
  library(excelR)
  library(reshape2)
}

source("pairwise_intersect.R")

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
  data <- data.frame(Rowname=elements,data)
  return(data)
}

Upset_comb <- function(data){
  zz <- data.frame(Rowname=data[,1],group=apply(data[,-1]==1,1,function(a) paste0(colnames(data[,-1])[a], collapse = "&")))
  zz <- dcast(zz, Rowname~group, value.var="Rowname", fun.aggregate=length)
  if ("Var.2" %in% names(zz)) zz[['Var.2']] <- NULL
  zz <- lapply(zz,function(x) unique(zz[x>0,1]))
  zz <- data.frame(sapply(zz, "length<-", max(lengths(zz))),stringsAsFactors = FALSE, check.names = FALSE)
  zz
}

Counter <- function(data, empty_intersects = FALSE){
  data <- data[,-1]
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

read.gmt = function(file){ #https://rdrr.io/bioc/qusage/src/R/qusage.R
  if(!grepl("\\.gmt$",file)[1]){stop("Pathway information must be a .gmt file")}
  geneSetDB = readLines(file)                                ##read in the gmt file as a vector of lines
  geneSetDB = strsplit(geneSetDB,"\t")                       ##convert from vector of strings to a list
  names(geneSetDB) = sapply(geneSetDB,"[",1)                 ##move the names column as the names of the list
  geneSetDB = lapply(geneSetDB, "[",-1:-2)                   ##remove name and description columns
  geneSetDB = lapply(geneSetDB, function(x){x[which(x!="")]})##remove empty strings
  return(geneSetDB)
}

Univ_reader <- function(input_type,inFile,string,sep_,header_,sep_row_,thequote,dataread,dedup=TRUE,example="exp",inputfile_local=NULL){
  if(string != ""){
    string <- as.list(unlist(strsplit(string, ",")))
    names <- lapply(string, function(x){x <- unlist(strsplit(x, "=")); x <- x[1]})
    names <- unlist(lapply(names, function(x){x <- gsub(" ", "", x)}))
    values <- as.numeric(unlist(lapply(string, function(x){x <- unlist(strsplit(x,"=")); x <- x[2]})))
    names(values) <- names
    data <- fromExpression(values)
    if(!"Rowname" %in% names(data)){
      data <- data.frame(Rowname=row.names(data),data)      
    }
  } else if(is.null(inFile) == F || input_type %in% c('local')){
    if (input_type == 'list'){
      if (grepl("\\.gmt$",inFile$datapath)[1]){
        data <- read.gmt(inFile$datapath)
      } else {
        data <- read_delim(inFile$datapath, sep_ , escape_double = FALSE, trim_ws = TRUE, col_names = header_)
        if (sep_row_ != "No")
        {
          data <- lapply(data, function(x) unique(unlist(strsplit(x,sep_row_))))
        }
      }
      data <- fromList(lapply(as.list(data), function(x) x[!myisna(x)]))
    } else if (input_type %in% c('binary','local')){
      if(thequote == "No")
      {
        thequote <- ""
      }
      if (input_type %in% c('local') && !is.null(inputfile_local)){
        ftoread <- inputfile_local
      } else {
        ftoread <- inFile$datapath
      }
      data <- read.csv(ftoread, header = header_,
                       sep = sep_, quote = thequote, row.names = NULL)
    }
  } else if (nrow(dataread)) {
    #browser()
    data <- dataread
  }else{
    if(example == "exp"){
      data<- fromExpression(c('H3K4me2&H3K4me3'=1632,'H3K4me2&H3K4me3&H3K27me3'=575,'H3K27me3'=2517,'H3K4me3&H3K27me3'=1553,'H3K4me3'=3296,'H3K4me2&H3K27me3'=1903,'H3K4me2'=6029,'H3K27ac&H3K4me2&H3K4me3&H3K27me3'=723,'H3K27ac&H3K4me2&H3K4me3'=1750,'H3K27ac&H3K4me2'=2134,'H3K27ac&H3K4me2&H3K27me3'=169,'H3K27ac&H3K4me3'=813,'H3K27ac&H3K4me3&H3K27me3'=29,'H3K27ac&H3K27me3'=760,'H3K27ac'=4216))
      if(!"Rowname" %in% names(data)){
        data <- data.frame(Rowname=row.names(data),data)      
      }
    } else {
      data <- read_delim('data/Whyte_et_al_2013_SEs_genes.csv', ",", escape_double = FALSE, trim_ws = TRUE, col_names = TRUE)
      if (sep_row_ != "No")
      {
        data <- lapply(data, function(x) unique(unlist(strsplit(x,sep_row_))))
      }
      data <- fromList(lapply(as.list(data), function(x) x[!myisna(x)]))
    }
  }
  return(data)
}

#sever code
shinyServer(function(input, output, session) {
  
  output$hot_venn = renderExcel({
    test <- venn_data_filtered() #list
    ncol <- length(test)
    test <- data.frame(sapply(test, "length<-", max(lengths(test))),stringsAsFactors = FALSE)
    excelTable(data=test,colHeaders=names(test),wordWrap=TRUE,search=TRUE,loadingSpin=TRUE,showToolbar=TRUE,autoWidth=TRUE,defaultColWidth=10,pagination=500)
  })
  
  output$hot_upset = renderExcel({
    test <- upset_data_filtered() # bin
    test <- lapply(test[,-1], function(x) as.character(test[,1][x>0])) 
    ncol <- length(test)
    test <- data.frame(sapply(test, "length<-", max(lengths(test))),stringsAsFactors = FALSE)
    excelTable(data=test,colHeaders=names(test),wordWrap=TRUE,search=TRUE,loadingSpin=TRUE,showToolbar=TRUE,autoWidth=TRUE,defaultColWidth=10,pagination=500)
  })
  
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
  set7_color <- reactive({
    return(input$set7_color)
  })
  set8_color <- reactive({
    return(input$set8_color)
  })
  set9_color <- reactive({
    return(input$set9_color)
  })
  set10_color <- reactive({
    return(input$set10_color)
  })
  
  venn_data_excel <- eventReactive(input$update_tbl_venn,{
    if (!is.null(input$hot_venn)) {
      data = excel_to_R(input$hot_venn)
      data <- lapply(data, function(x) x[!myisna(x)])
      fromList(data)
    } else {
      data.frame()
    }
  },ignoreNULL = FALSE
  )
  
  venn_data <- reactive({
    inFile <- input$file_venn
    string <- input$venn_comb
    inputfile_local <- input$inputfile_local
    string <- gsub("\n", "", string)
    input_type <- input$venn_input_type
    
    data <- Univ_reader(input_type,inFile,string,input$sep_venn,input$header_venn,input$sep_row_venn,input$quote_venn,venn_data_excel(),example="gene",inputfile_local=inputfile_local)
    return(data)
  })
  
  set_names_venn <- reactive({
    data <- venn_data()
    names0 <- as.character(names(data))
    if (input$venn_input_type == 'binary') {
      notbin <- apply(data,2,function(x)any(!x %in% c(0,1)))
      nameother <- names0[notbin]
      names0 <- names0[!notbin]
    } else {
      names0 <- names0[names0 != "Rowname"]
      nameother <- c("Rowname")
    }
    return(list(names0,nameother))
  })
  
  output$rowname_venn <- renderUI({
    if (input$venn_input_type == 'binary') {
      sets_venn <- selectInput('rowname_venn', label = "Rowname for binary data",
                               choices = set_names_venn()[[2]],
                               multiple = F, selected = set_names_venn()[[2]][1])
    } else {
      sets_venn <- NULL
    }
    return(sets_venn)
  })
  
  output$sets_venn <- renderUI({
    goodname <- set_names_venn()[[1]]
    goodnamesel <- goodname[1:3]
    sets_venn <- selectInput('sets_venn', label = "Select sets",
                             choices = goodname,
                             multiple = T, selectize = T, selected = goodnamesel)
    return(sets_venn)
  })
  
  selected_names_venn <- reactive({
    selected_names_venn <- as.character(c(input$sets_venn))
  })
  
  venn_data_filtered <- reactive({
    sep_row_ <- input$sep_row_venn
    dedup <- input$dedup_venn
    data <- venn_data() # bin w/ rowname
    
    i_rowname_venn <- 1
    if (input$venn_input_type == 'binary' && is.data.frame(data)) {
      if(!is.null(input$rowname_venn)){
        i_rowname_venn <- which(names(data) == input$rowname_venn)
      }
      data <- data[!myisna(data[,i_rowname_venn]),]
    }
    
    if(is.null(input$sets_venn)){
      thenames <- set_names_venn()[[1]][1:3]
    }else{
      thenames <- c(selected_names_venn())
    }
    data <- data.frame(Rowname=data[,i_rowname_venn],data[,names(data) %in% thenames])
    
    if (sep_row_ != "No")
    {
      aa <- data.frame(
        do.call(rbind,
                apply(data, 1, function(x) {
                  do.call(expand.grid, list(strsplit(x, sep_row_),stringsAsFactors = FALSE))
                })
        ))
      data <- as.data.frame(mapply(FUN = as,aa,sapply(data,class),SIMPLIFY = FALSE))
    }
    
    if (dedup)
    {
      touse4mean <- data[,-1]
      itouse4mean <- apply(touse4mean,2,is.numeric)
      themean <- apply(touse4mean[,itouse4mean],1,sum)
      genes <- data[,1]
      theorder <- order(themean,decreasing=T)
      genes <- genes[theorder]
      data <- data[theorder,]
      data <- data[ave(as.character(genes),as.character(genes),FUN=seq_along) == 1,]
    }  
    
    i_rowname_venn <- 1
    data <- lapply(data[,-i_rowname_venn,drop=F], function(x) as.character(data[,i_rowname_venn][x>0]))
    data <- lapply(data, function(x) x[!myisna(x)])
    
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
    if (venn_color_type() == 'custom'){
      venn_gp <- VennThemes(compute.Venn(venn_combinations()))
    } else {
      venn_gp <- VennThemes(compute.Venn(venn_combinations()),colourAlgorithm=venn_color_type())
    }
    
    venn_gp$SetText <- lapply(venn_gp$SetText,function(x) {x$fontsize<-venn_labelsize(); return(x)})
    venn_gp$FaceText <- lapply(venn_gp$FaceText,function(x) {x$cex<-venn_cex(); return(x)})
    venn_gp$Set <- lapply(venn_gp$Set,function(x) {x$lwd<-venn_lwd(); return(x)})
    venn_gp$Set <- lapply(venn_gp$Set,function(x) {x$lty<-venn_lty(); return(x)})
    
    #browser()
    
    if (venn_color_type() == 'custom'){
      inputcolors <- c(set1_color(),set2_color(),set3_color(),set4_color(),set5_color(),set6_color(),set7_color(),set8_color(),set9_color(),set10_color())
      for (i in 1:length(venn_gp$Set)){
        thisset <- names(venn_gp$Set)[i]
        venn_gp$Set[[thisset]]$col <- inputcolors[i]
        venn_gp$SetText[[thisset]]$col <- inputcolors[i]
      }
      if(length(venn_gp$Face) < 9){
        for (i in 2:length(venn_gp$Face)){
          thisset <- names(venn_gp$Face)[i]
          venn_gp$Face[[thisset]]$col <- inputcolors[i+2]
          venn_gp$Face[[thisset]]$fill <- inputcolors[i+2]
        }
      }
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
      pixelratio <- 1
      pixeldiv <- 72
      
      if(input$filetype_venn == "PNG")
        png(file, width=width*pixelratio, height=height*pixelratio, units = "px", res=72*pixelratio)
      else if(input$filetype_venn == "SVG")
        svg(file, width=width/pixeldiv, height=height/pixeldiv)
      else if(input$filetype_venn == "TIFF")
        tiff(file, width=width*pixelratio, height=height*pixelratio, units = "px")
      else
        pdf(file, width=width/pixeldiv, height=height/pixeldiv)
      
      plot(compute.Venn(venn_combinations(), doWeights = doWeights(), doEuler = doEuler(), type = get_venn_type()),
           gp = get_venn_gp(),
           show = list(Universe = FALSE)
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
    if(is.null(upset_data_filtered()) == T){
      h5("There is no data entered. Please upload your data to draw UpSet plot here!")
    }
    else{
      HTML(" ")
    }
  })
  
  upset_data_excel <- eventReactive(input$update_tbl_upset,{
    if (!is.null(input$hot_upset)) {
      data = excel_to_R(input$hot_upset)
      data <- lapply(data, function(x) x[!myisna(x)])
      fromList(data)
    } else {
      data.frame()
    }
  },ignoreNULL = FALSE
  )
  
  upset_data <- reactive({  
    inFile <- input$file_upset
    string <- input$comb_upset
    inputfile_local <- input$inputfile_local
    string <- gsub("\n", "", string)
    input_type <- input$upset_input_type
    
    data <- Univ_reader(input_type,inFile,string,input$sep_upset,input$header_upset,input$sep_row_upset,input$quote_upset,upset_data_excel(),example="exp",inputfile_local=inputfile_local)
    return(data)
  })
  
  set_names_upset <- reactive({
    data <- upset_data()
    names0 <- as.character(names(data))
    if (input$upset_input_type == 'binary') {
      notbin <- apply(data,2,function(x)any(!x %in% c(0,1)))
      nameother <- names0[notbin]
      names0 <- names0[!notbin]
    } else {
      names0 <- names0[names0 != "Rowname"]
      nameother <- c("Rowname")
    }
    return(list(names0,nameother))
  })
  
  output$rowname_upset <- renderUI({
    if (input$upset_input_type == 'binary') {
      sets_upset <- selectInput('rowname_upset', label = "Rowname for binary data",
                                choices = set_names_upset()[[2]],
                                multiple = F, selected = set_names_upset()[[2]][1])
    } else {
      sets_upset <- NULL
    }
    return(sets_upset)
  })
  
  output$sets_upset <- renderUI({
    goodname <- set_names_upset()[[1]]
    goodnamesel <- goodname[1:3]
    sets_upset <- selectInput('sets_upset', label = "Select sets",
                              choices = goodname,
                              multiple = T, selectize = T, selected = goodnamesel)
    return(sets_upset)
  })
  
  setSizes <- reactive({
    if(is.null(upset_data()) != T){
      sizes <- colSums(upset_data()[,set_names_upset()[[1]]])
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
  
  selected_names_upset <- reactive({
    selected_names_upset <- as.character(c(input$upset_sets))
  })
  
  output$sets_upset <- renderUI({
    if(is.null(upset_data()) == T){
      sets <-  selectInput('upset_sets', label="Select at least two sets ",
                           choices = NULL,
                           multiple=TRUE, selectize=TRUE, selected = selected_names_upset())
    } else {
      data <- upset_data()[,set_names_upset()[[1]]]
      topfive <- colSums(data)
      topfive <- as.character(head(names(topfive[order(topfive, decreasing = T)]), 5))
      sets <- selectInput('upset_sets', label="Select sets ",
                          choices = set_names_upset()[[1]],
                          multiple=TRUE, selectize=TRUE, selected = topfive)
    }
    return(sets)
  })
  
  upset_data_filtered <- reactive({
    sep_row_ <- input$sep_row_upset
    dedup <- input$dedup_upset   
    data <- upset_data() # df
    rowname_upset <- input$rowname_upset
    
    i_rowname_upset <- 1

    if (input$upset_input_type == 'binary' && is.data.frame(data)) {
      
      if(!is.null(rowname_upset)){
        i_rowname_upset <- which(names(data) == rowname_upset)
      }
      data <- data[!myisna(data[,i_rowname_upset]),]
    }
    
    if(is.null(input$upset_sets)){
      thenames <- set_names_upset()[[1]][1:3]
    }else{
      thenames <- c(selected_names_upset())
    }
    data <- data.frame(Rowname=data[,i_rowname_upset],data[,names(data) %in% thenames])
    
    
    if (sep_row_ != "No")
    {
      aa <- data.frame(
        do.call(rbind,
                apply(data, 1, function(x) {
                  do.call(expand.grid, list(strsplit(x, sep_row_),stringsAsFactors = FALSE))
                })
        ))
      data <- as.data.frame(mapply(FUN = as,aa,sapply(data,class),SIMPLIFY = FALSE))
    }
    
    if (dedup)
    {
      touse4mean <- data[,-1]
      itouse4mean <- apply(touse4mean,1,is.numeric)
      themean <- apply(touse4mean[,itouse4mean],1,sum)
      genes <- data[,1]
      theorder <- order(themean,decreasing=T)
      genes <- genes[theorder]
      data <- data[theorder,]
      data <- data[ave(as.character(genes),as.character(genes),FUN=seq_along) == 1,]
    }  
    
    return(data)
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
    
    if(length(upset_data_filtered()) == 0){stop()}
    if(length(selected_names_upset()) == 1){
      stop()
    }
    upset(data = upset_data_filtered(), 
          nintersects = input$nintersections,
          point.size = input$pointsize,
          line.size = line_size(),
          sets = selected_names_upset(),
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
    width = upset_width,
    height = upset_height
  )
  
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
      
      base::print(upset(data = upset_data_filtered(), 
                        nintersects = input$nintersections,
                        point.size = input$pointsize,
                        line.size = line_size(),
                        sets = selected_names_upset(),
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
                                       input$intersection_size_numbers_scale)
      )
      )
      dev.off()
    }
  )
  
  output$upsetDownExcel <- downloadHandler(
    filename = function(){
      paste("UpSet", tolower(input$filetype_upset_excel), "csv", sep =".")
    }, 
    content = function(file){
      data <- upset_data_filtered()
      #browser()
      if(input$filetype_upset_excel == "Freq"){
        data <- Counter(data)
      } else if(input$filetype_upset_excel == "Combinations"){
        data <- Upset_comb(data)
      }
      write.table(data.frame("Row"=row.names(data),data,check.names = FALSE),file,na = "",sep=",",row.names = FALSE)
    }
  )
  
  localdata <- reactive({
    ff <- "local.lst"
    if (file.exists(ff)){
      gtmp <- read.table(ff,header=FALSE)$V1
    } else { gtmp <- c() }
    gtmp <- as.character(gtmp)
    if (length(gtmp) > 0){
      ggg_symbol_choices <- selectInput('inputfile_local', label = "Local Data",
                                        choices = gtmp, multiple = F, selected = gtmp[1])
    } else {
      ggg_symbol_choices <- NULL
    }
    return(ggg_symbol_choices)
  })
  
  output$ggg_local_data_venn <- renderUI({
    return(localdata())
  })
  output$ggg_local_data_upset <- renderUI({
    return(localdata())
  })
  output$ggg_local_data_pairwise <- renderUI({
    return(localdata())
  })
  
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
