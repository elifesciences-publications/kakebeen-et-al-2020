## common wrangling and single cell functions to source

# find if gene is in data set
FindGene <- function(object, gene){
  grep(pattern = gene, rownames(object), value = TRUE, ignore.case = TRUE)
}

# id how many cells express a gene
gene.ncells <- function(object, gene){
  temp.table <- as.data.frame(GetAssayData(object = object)[gene,] >0)
  colnames(temp.table) <- "t"
  temp.table[,1] <- as.numeric(temp.table[,1])
  temp.table <- temp.table %>%
    filter(temp.table[,1] > 0)
  return(nrow(temp.table))
}

## To get all genes from a treemap or go output table and put into form for plotting in single cell RNA-Seq

concat.GO <- function(input){
  ## make list of call genes in nephron family.
  temp <- strsplit(as.character(input), "," ) # makes a list of all genes fron column 5
  temp <- stringi::stri_paste_list(temp, sep = ",", collapse = ",") # concatenates list of lists into one list
  temp <- strsplit(temp, ",") # split list to list of individual character strings
  temp <- unique(as.vector(temp[[1]])) # save out 1 list as character vector and get unique values to eliminate repeats
  temp <- tolower(temp)
  #assign(out, temp, envir = .GlobalEnv)
}

## To get all genes from treemap or go output table and put in form for searching via grep
concat.grep <- function(input){
  ## make list of call genes in nephron family.
  temp <- strsplit(as.character(input), "," ) # makes a list of all genes fron column 5
  temp <- stringi::stri_paste_list(temp, sep = ",", collapse = ",") # concatenates list of lists into one list
  temp <- strsplit(temp, ",") # split list to list of individual character strings
  temp <- unique(as.vector(temp[[1]])) # save out 1 list as character vector and get unique values to eliminate repeats
  temp <- tolower(temp) # lower case
  temp <- list(temp) # relist after find unique values
  temp <- unlist(stringi::stri_paste_list(temp, sep = "|", collapse = ",")) # unlist and seperate with "|" for grepping purposes
  #assign(out, temp, envir = .GlobalEnv)
}

## To get all genes from treemap or go output table and put in form for searching via grep *** for semicolon
concat.grep_sc <- function(input){
  ## make list of call genes in nephron family.
  temp <- strsplit(as.character(input), ";" ) # makes a list of all genes fron column 5
  temp <- stringi::stri_paste_list(temp, sep = ",", collapse = ",") # concatenates list of lists into one list
  temp <- strsplit(temp, ",") # split list to list of individual character strings
  temp <- unique(as.vector(temp[[1]])) # save out 1 list as character vector and get unique values to eliminate repeats
  temp <- tolower(temp) # lower case
  temp <- list(temp) # relist after find unique values
  temp <- unlist(stringi::stri_paste_list(temp, sep = "|", collapse = ",")) # unlist and seperate with "|" for grepping purposes
  #assign(out, temp, envir = .GlobalEnv)
}

## split violin plots
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

## get seurat data and put it in plottable format
seuratTOplot <- function(dataframe, gene, outs) {
  for(i in 1:length(gene)){
    temp.table <- data.frame(dataframe[,gene[i]])
    temp.table$gene <- gene[i]
    temp.table$condition <- rownames(dataframe)
    colnames(temp.table) <- c("expression", "gene", "condition")
    assign(outs[i], temp.table, envir = .GlobalEnv)
  }
}

## blank plot
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
library(scales)

## Taking a gProfiler2 object to a useable table
GO2table <- function(GOResults){
  temp.table <- GOResults$result # get table from list object
  #colnames(all_GOterms)
  temp.table <- temp.table[,c(3, 9, 11, 16)] # select for p_value, term_id, term_name, intersection
  temp.table$logp <- -1*log10(temp.table[,1]) # make log10pvalue column
  temp.table <- temp.table[order(-temp.table[,5]),] # order table by log10
  temp.table[,3] <- factor(temp.table[,3], levels = c((temp.table[, 3])))
  return(temp.table)
}
