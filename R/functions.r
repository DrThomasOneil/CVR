#setup
loadPack <- function() {
  library(cli)
  cat("\n--------------------------------------\n")
  cat(style_bold(col_magenta("\n***Installing General Packages***\n\n")))
  not <- c(); not2 <- c()
  
  #FUTURE TOM: ADD PACKAGES HERE!
  packages1 <- c("ggplot2", "rstudioapi", "rmarkdown", 'tidyr', "cli", "knitr", "dplyr", "Seurat","SeuratObject",
                 "DoubletFinder", "SeuratDisk", "flipPlots",'stringr', "crayon","Matrix", "cowplot", 'scater', "BiocParallel",
                 "ComplexHeatmap","xlsx", "ggpubr")#, "Test")
  
  for (i in 1:length(packages1)){
    if(requireNamespace(packages1[i], quietly = TRUE)==F) {
      cat(paste(style_bold(col_red(packages1[i])), "has not been installed\n"))
      not <- c(not,i)
    } else {
      suppressWarnings(suppressMessages(library(as.character(packages1[i]), character.only = TRUE)))
      cat(col_yellow(packages1[i]), "is loaded!\n")
    }
  }
  cat("\n--------------------------------------\n")
  
  if (length(not) > 0){
    cat(style_bold(bg_red("\n  **IMPORTANT**  ")),
        style_bold(col_yellow("\n\nYou need to install: \n")),
        
        paste(paste(c(packages1[not]), collapse=", ")),
        "\n\n--------------------------------------",
        
        "\n\n Use:\n - install.packages(),\n - BiocManager::install() or, \n - use Google to find installation instructions.\n\n", style_bold(col_green("Then run this function again!\n\n")))
  } else {
    cat("",col_green(style_bold("\n All packages are loaded!\n\n Happy Coding! :)\n\n")))
  }
}
loadPack()
theme_set(theme_classic())

# Functions ---------------------------------------------------------------

plotSankey<-function(seuratObj,idvar=c("varRes.0.3","emt_res.0.3")){
  require(flipPlots)
  message('try install_github("Displayr/flipPlots") if this doesnt work')
  require(dplyr)
  seuratObj@meta.data[,match(idvar,colnames(seuratObj@meta.data))] %>% arrange(.[,1]) %>% group_by_all() %>% summarise(COUNT = n()) ->> my.data
  #my.data<-as.factor(my.data[,1])
  SankeyDiagram(my.data[, -grep("COUNT",colnames(my.data))],link.color = "Source",weights = my.data$COUNT,,max.categories = 100)
}

proportions <- function(data, ident.1, ident.2, position) {
  x<- FetchData(data,c(ident.1,ident.2))
  colnames(x) <- c('ident.2', 'ident.1')
  x%>% group_by(ident.1) %>%
    mutate(prop=1/length(ident.2)) %>%
    ungroup() %>%
    group_by(ident.2,ident.1) %>%
    summarise(totprop=sum(prop)) %>%
    ggplot(aes(x=ident.2,fill=ident.1,y=totprop)) +
    geom_bar(position=position, stat='identity') + theme(axis.text.x =
                                                           element_text(angle = 45,hjust=1))+scale_y_continuous(name="Cluster
    Proportion")+ theme_classic()
}

process <- function(dat=dat, dimuse = 1:15, features=800, verbose=F, reduction.name=NULL){
  dat <- NormalizeData(dat, verbose=verbose)
  dat <- FindVariableFeatures(dat, nfeatures=features, verbose=verbose)
  dat <- ScaleData(dat, features=VariableFeatures(dat), verbose=verbose)
  dat <- RunPCA(dat, verbose=verbose, reduction.name = paste0("pca", reduction.name))
  #dat <- RunUMAP(dat, dims=dimuse, verbose=verbose, reduction = paste0("pca", reduction.name), reduction.name = paste0("umap", reduction.name))
  #dat <- RunTSNE(dat, dims=dimuse, check_duplicates=F, verbose=verbose, reduction = paste0("pca", reduction.name),reduction.name = paste0("tsne", reduction.name))
  return(dat)
}

source("https://raw.githubusercontent.com/tomoneil58/LabCode/main/HPA/HPA.R")

predictionHeat <- function(ref, query, refID = "ident", queryID = "ident", norm=F, crow=T,ccol=T, return.plot=F, return.seurat=F, col.name="predictedID", var.feat="") {
  if(length(var.feat)==1) {
    predictions <- TransferData(
      anchorset = FindTransferAnchors(reference = ref, query = query, features=VariableFeatures(ref)),
      refdata = FetchData(ref, refID)[,1]
    )
  } else {
    predictions <- TransferData(
      anchorset = FindTransferAnchors(reference = ref, query = query, features=var.feat),
      refdata = FetchData(ref, refID)[,1]
    )
  }
  
  predictions$orig =FetchData(query, queryID)[,1]
  df <- as.data.frame(matrix(data=NA,ncol=length(unique(predictions$orig)), nrow=length(unique(predictions$predicted.id))))
  colnames(df) = unique(predictions$orig); 
  rownames(df) = paste0("prediction.score.",gsub(" ", ".",unique(predictions$predicted.id)))
  df2 <- as.data.frame(matrix(data=NA,ncol=length(unique(predictions$orig)), nrow=length(unique(predictions$predicted.id))))
  colnames(df2) = unique(predictions$orig); 
  rownames(df2) = paste0("prediction.score.",gsub(" ", ".",unique(predictions$predicted.id)))
  
  for(col in 1:ncol(df)) {
    for(row in 1:nrow(df)){
      x = mean(predictions[predictions$orig==colnames(df)[col], rownames(df)[row]])
      x2 = mean(predictions[predictions$orig==colnames(df)[col], 
                            rownames(df)[row]]/
                  predictions[predictions$orig==colnames(df)[col], 
                              "prediction.score.max"])
      
      if(is.na(x) | is.infinite(x)){
        x=0
      }
      if(is.na(x2) | is.infinite(x2)){
        x2=0
      }
      df[row,col] <-x 
      df2[row,col] <-x2
    }
  }
  if (norm) {
    df = df2
  }
  if(!ccol) {
    ccol = F
  }
  if(!crow) {
    crow =F
  }
  print(ComplexHeatmap::Heatmap(na.omit(df), cluster_columns = ccol, cluster_rows = crow))
  #dont return both a metatable
  if(return.plot*return.seurat ==1) {
    cat(style_bold(col_red("\n***ERROR***\n\n")))    
    cat(style_bold(col_yellow("\n***ERROR***\n\n")))

  } else {
    if(return.plot){
      plot = ComplexHeatmap::Heatmap(na.omit(df), cluster_columns = ccol, cluster_rows = crow)
      return(plot)
    } 
    if(return.seurat){
      query <- AddMetaData(query, metadata = predictions$predicted.id, col.name=col.name)
      return(query)
    } 
    
  }


}

topm <- function(data, min.diff.pct = 0.01, n=40) {
  FindAllMarkers(data, only.pos=T, min.diff.pct = min.diff.pct) %>%
    filter(p_val_adj <0.0001) %>%
    group_by(cluster) %>%
    top_n(n=n, wt = avg_log2FC)
}


