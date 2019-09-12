Single Cell Figures
================
Anneke Kakebeen
13 Aug 2019

-   [Set up for scRNA-Seq Figures](#set-up-for-scrna-seq-figures)
-   [Read in scRNASeq Objects](#read-in-scrnaseq-objects)
-   [Figure 4a: UMAP of all neural cells](#figure-4a-umap-of-all-neural-cells)
-   [Figure 4b: UMAP of all neural cells colored for condition (uninjured or 24hpa)](#figure-4b-umap-of-all-neural-cells-colored-for-condition-uninjured-or-24hpa)
-   [Figure 4C: Marker dotplot for neural cell clusters](#figure-4c-marker-dotplot-for-neural-cell-clusters)
-   [Figure 4D: cell cluster proportions](#figure-4d-cell-cluster-proportions)
    -   [Set up data table](#set-up-data-table)
    -   [Pie Charts](#pie-charts)
    -   [UMAPS](#umaps)
    -   [Pie charts of cells using three state condition](#pie-charts-of-cells-using-three-state-condition)
-   [Figure 4E: Cell cycle phase prediction](#figure-4e-cell-cycle-phase-prediction)
    -   [Table set up](#table-set-up)
    -   [Pie Charts](#pie-charts-1)
    -   [UMAP](#umap)
-   [Figure 5C: average expression of meis1 and pbx3 over sc clusters](#figure-5c-average-expression-of-meis1-and-pbx3-over-sc-clusters)
-   [Figure 5D: featuremap to show expression of pbx3 and meis1 across neural clusters](#figure-5d-featuremap-to-show-expression-of-pbx3-and-meis1-across-neural-clusters)

Set up for scRNA-Seq Figures
----------------------------

``` r
# load packages
library(Seurat)
library(LaCroixColoR)
library(ggplot2)
library(cowplot)
```

    ## 
    ## ********************************************************

    ## Note: As of version 1.0.0, cowplot does not change the

    ##   default ggplot2 theme anymore. To recover the previous

    ##   behavior, execute:
    ##   theme_set(theme_cowplot())

    ## ********************************************************

``` r
library(tidyr)
library(DT)
```

    ## 
    ## Attaching package: 'DT'

    ## The following object is masked from 'package:Seurat':
    ## 
    ##     JS

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggpubr)
```

    ## Loading required package: magrittr

    ## 
    ## Attaching package: 'magrittr'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     extract

    ## 
    ## Attaching package: 'ggpubr'

    ## The following object is masked from 'package:cowplot':
    ## 
    ##     get_legend

``` r
library(rmarkdown)
library(kableExtra)
```

    ## 
    ## Attaching package: 'kableExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     group_rows

``` r
# load functions
source("~/Desktop/pax6 paper/Final Markdowns/common_source_functions.R")

# Set colors
my.color.neural <- c("#2CB11B", "#E9E4A6", "#DFCEE0" , "#E9A17C", "#1BB6AF" , "#172869" , "#FF3200")
condition.color <- c( "#6F0909", "#C0C0C0")
cc_colors <- c((lacroix_palette("Lime"))[c(1,3)], "black")

# Set working directory
setwd("~/Desktop/pax6 paper/Final Markdowns/SingleCell/analysis/")
```

Read in scRNASeq Objects
------------------------

``` r
# Read in pax.neural seurat object
pax.neural <- readRDS("~/Desktop/pax6 paper/Final Markdowns/SingleCell/UMAP/pax.neural_FINAL.RDS")

# read in plotting data
pax.neural.plotting<- read.table("~/Desktop/pax6 paper/Final Markdowns/SingleCell/UMAP/pax.neural_plotting.txt", sep = "\t")

# factor plotting table for plotting
pax.neural.plotting$AK_Names <- factor(pax.neural.plotting$AK_Names, levels = c("NSC", "Transition", "Neural.A","Neural.B", "Neural.C", "Neural.D", "Prdm14" ))
pax.neural.plotting$condition <- factor(pax.neural.plotting$condition, levels = c("uninj", "dmso", "iwr"))
pax.neural.plotting$three_state <- factor(pax.neural.plotting$three_state, levels = c("NSC", "Differentiating", "Neuron"))
pax.neural.plotting$Laevis <- factor(pax.neural.plotting$Laevis, levels = c("Spinal Cord Progenitor", "Differentiating Neuron" , "Interneurons",  "Motor Neuron" , "Vulnerable Motor Neuron", "Motor Neuron (leptin+)" , "Dopaminergic Neurons" ))
pax.neural.plotting$Phase <- factor(pax.neural.plotting$Phase, levels = c("S", "G2M", "G1"))
```

Figure 4a: UMAP of all neural cells
-----------------------------------

``` r
# change idents of suerat object
Idents(pax.neural) <- pax.neural$Laevis
# plot wit hseurat umap
DimPlot(pax.neural, cols = my.color.neural) + NoAxes() + NoLegend()
```

![](README_files/figure-markdown_github/4a-1.png)

Figure 4b: UMAP of all neural cells colored for condition (uninjured or 24hpa)
------------------------------------------------------------------------------

``` r
ggplot(pax.neural.plotting, aes(x=UMAP_1, y=UMAP_2, color=pax.neural.plotting$mock_24))+
  geom_point(stat="identity", aes(alpha = 0.6)) + 
  scale_color_manual(values = condition.color) + 
  Seurat:::NoAxes() + 
  blank_theme +
  theme(legend.position = "none")
```

![](README_files/figure-markdown_github/4b-1.png)

Figure 4C: Marker dotplot for neural cell clusters
--------------------------------------------------

``` r
# read in marker table
DEClusters.neural <- read.table("~/Desktop/pax6 paper/Final Markdowns/SingleCell/analysis/outs/DEClusters.neural.txt", sep = "\t", header =TRUE, row.names = 1)
# reassign class of genes to factors
DEClusters.neural$gene <- as.character(DEClusters.neural$gene)
# get
markers <- unique(c("sox2", (filter(DEClusters.neural, cluster =="DEClusters.NSC")[1:3,7]),
             (filter(DEClusters.neural, cluster == "DEClusters.Transition")[1:3,7]),
             (filter(DEClusters.neural, cluster == "DEClusters.Neural.A")[1:3,7]),
             (filter(DEClusters.neural, cluster == "DEClusters.Neural.B")[1:3,7]),
             (filter(DEClusters.neural, cluster == "DEClusters.Neural.C")[1:3,7]),
             (filter(DEClusters.neural, cluster == "DEClusters.Neural.D")[1:3,7]),
             "lepr",
              (filter(DEClusters.neural, cluster == "DEClusters.Prdm14")[1:3,7]),
             "prdm14", "prdm1"
             ))


DefaultAssay(pax.neural) <- "RNA" # change to RNA assay
pax.neural$AK_Names <- factor(pax.neural$AK_Names, levels = c("NSC", "Transition", "Neural.A", "Neural.B", "Neural.C", "Neural.D", "Prdm14")) # factor names in AK_Names for ordering purposes
Idents(pax.neural) <- pax.neural$Laevis # change object identity to AK_Names
pax.neural$Laevis <- factor(pax.neural$Laevis, levels = c("Spinal Cord Progenitor", "Differentiating Neuron", "Interneurons", "Vulnerable Motor Neuron", "Motor Neuron", "Motor Neuron (leptin+)", "Dopaminergic Neurons"))
Idents(pax.neural) <- factor(Idents(pax.neural), levels = (levels(pax.neural$Laevis))) # Reverse the order of names fo the plot

# plot dotplot
DotPlot(pax.neural, features = (markers), cols = c("#172869", "#F7AA14"))  + RotatedAxis() +  NoLegend() + coord_flip()
```

![](README_files/figure-markdown_github/4C-1.png)

Figure 4D: cell cluster proportions
-----------------------------------

### Set up data table

``` r
Idents(object = pax.neural) <- "celltype_condition"

# Get number of cells in each cluster for each condition into data frame
nCells.neural <- as.data.frame(table(pax.neural@active.ident))
# Make condition and cluster columns
nCells.neural$condition <- gsub(nCells.neural$Var1, pattern = "_.*", replacement = "")
nCells.neural$cluster <- gsub(nCells.neural$Var1, pattern = ".*_", replacement = "")
nCells.neural <- nCells.neural[order(nCells.neural$condition),]

# get total numner of cells for each condition
totals <-as.data.frame(nCells.neural %>%
    group_by(condition) %>% 
  summarise(totals = sum(Freq)))
totals$totals <- as.numeric(totals$totals)

# Make ratio columns in new table expressing percent cluster of all neural cells
nCells.neural$percent <- c((nCells.neural[1:7,2]/totals[1,2])*100, (nCells.neural[8:14,2]/totals[2,2])*100, (nCells.neural[15:21,2]/totals[3,2])*100)


# factor table
nCells.neural$condition <- factor(nCells.neural$condition, levels = c("uninj", "dmso", "iwr"))
nCells.neural$cluster <- factor(nCells.neural$cluster, levels = c("NSC", "Transition", "Neural.A", "Neural.B", "Neural.C", "Neural.D", "Prdm14"))
```

### Pie Charts

``` r
# Filter table for unjured
plot.uninj <- filter(nCells.neural, condition == "uninj")
# plot with ggplots
p.uninj <- ggplot(plot.uninj, aes(x="", y=percent, fill=cluster))+
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = my.color.neural) +
  coord_flip() +
  coord_polar("y", start = 0) +
  ggtitle("Uninjured") +
  blank_theme +
  theme(legend.position = "bottom")

# Filter table for dmso
plot.dmso <- filter(nCells.neural, condition == "dmso")
# plot
p.dmso <- ggplot(plot.dmso, aes(x="", y=percent, fill=cluster))+
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = my.color.neural) +
  coord_flip() +
  coord_polar("y", start = 0) +
  ggtitle("Uninjured") +
  blank_theme +
  theme(legend.position = "bottom")

# plot both pie charts side by side
plot_grid(p.uninj, p.dmso)
```

![](README_files/figure-markdown_github/uninjured%20cell%20proportions-1.png)

### UMAPS

``` r
# umap plot of cell types faceted by condition
ggplot(filter(pax.neural.plotting, condition %in% c("uninj", "dmso")), aes(x=UMAP_1, y=UMAP_2, color=AK_Names))+
  geom_point(stat="identity") +
  scale_color_manual(values = my.color.neural) +
  facet_grid(.~condition) +
  blank_theme +
  theme(legend.position = "right")
```

![](README_files/figure-markdown_github/uninjured%20cell%20proportions%202-1.png)

``` r
  Seurat:::NoAxes()
```

    ## List of 8
    ##  $ axis.line.x : list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  $ axis.line.y : list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  $ axis.title.x: list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  $ axis.title.y: list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  $ axis.text.x : list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  $ axis.text.y : list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  $ axis.ticks.x: list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  $ axis.ticks.y: list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  - attr(*, "class")= chr [1:2] "theme" "gg"
    ##  - attr(*, "complete")= logi FALSE
    ##  - attr(*, "validate")= logi TRUE

### Pie charts of cells using three state condition

``` r
three_state <- as.data.frame(table(pax.neural.plotting$three_state,  pax.neural.plotting$condition))[1:6,] # save out table of freuqency of cell types by three state criteria (only for uninj and dmso)
colnames(three_state) <- c("CellType", "condition", "Freq") # rename columns 
three_state_totals <- as.data.frame(table(pax.neural.plotting$condition))[1:2,] # get sum of each condition types
three_state$total <- c(rep(three_state_totals[1,2], 3), rep(three_state_totals[2,2], 3)) # make column in table for sums
three_state$percent <- three_state$Freq/three_state$total*100 # make percentage column

# plot three_state in pie chart
ggplot(three_state, aes(x="", y=percent, fill=CellType))+
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = my.color.neural) +
  coord_flip() +
  coord_polar("y", start = 0) +
  blank_theme +
  theme(legend.position = "bottom") +
  facet_grid(.~condition)
```

![](README_files/figure-markdown_github/three%20state%20pie-1.png)

Figure 4E: Cell cycle phase prediction
--------------------------------------

### Table set up

``` r
# extract data
cell.cycle.table1 <- as.data.frame(pax.neural[[]]) # get cell cycle data from seurat
cell.cycle.table <- as.data.frame(table(cell.cycle.table1$condition, cell.cycle.table1$AK_Names, cell.cycle.table1$Phase)) # get frequency of cells that have a specific phase in a cluster in a condition
cell.cycle.table <- cell.cycle.table[order(cell.cycle.table$Var1),] # order by cluster
colnames(cell.cycle.table) <- c("condition", "cluster", "phase", "frequency_cluster") # assign column names

# get sums of phases in a condition
cellcycle_sum <- as.data.frame(table(cell.cycle.table1$condition, cell.cycle.table1$Phase)) # get frequency of number of cells in a given phase and condition
cellcycle_sum <- cellcycle_sum[order(cellcycle_sum$Var1),] # order this by cluster
colnames(cellcycle_sum) <- c("condition", "phase", "frequency_phase") # assign column names
plot.table <- merge.data.frame(cell.cycle.table, cellcycle_sum, by = c("condition", "phase")) # merge two tables by condition and phase

# get sum of phases overal
condition_sum <- as.data.frame(table(cell.cycle.table1$condition)) # get frequency of number of cells in a phase
condition_sum <- condition_sum[order(condition_sum$Var1),]
colnames(condition_sum) <- c("condition","frequency_condition")
plot.table <- merge.data.frame(plot.table, condition_sum, by = c("condition"))

# make columns with percentages
plot.table$percent_phase <- (plot.table$frequency_phase/plot.table$frequency_condition)*100 # % cells in a phse per condition
plot.table$percent_cluster <- (plot.table$frequency_cluster/plot.table$frequency_condition)*100  # % cells in a cluster per condition

# factor phase for order in plot
plot.table$phase <- factor(plot.table$phase, levels = c("S", "G2M", "G1"))
```

### Pie Charts

``` r
# uninjured
# pie chart of uninjured cell cycle phases
plot.uninj <- unique(filter(plot.table, condition == "uninj")[,c(1,2,7),])
plot.uninj1 <- ggplot(plot.uninj, aes(x="", y=percent_phase, fill=phase))+
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = c(cc_colors)) +
  coord_polar("y") +
  blank_theme+
  theme(legend.position = "none")

# dmso
# stacked bar plot of cell cycle phases
plot.dmso <- unique(filter(plot.table, condition == "dmso")[,c(1,2,7),])
plot.dmso1 <- ggplot(plot.dmso, aes(x="", y=percent_phase, fill=phase))+
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = c(cc_colors)) +
  coord_polar("y") +
  blank_theme+
  theme(legend.position = "none")
```

### UMAP

``` r
# umap plot of cell types faceted by condition
ggplot(filter(pax.neural.plotting, condition %in% c("uninj", "dmso")), aes(x=UMAP_1, y=UMAP_2, color=Phase))+
  geom_point(stat="identity") +
  scale_color_manual(values = cc_colors) +
  facet_grid(.~condition) +
  blank_theme +
  theme(legend.position = "right")
```

![](README_files/figure-markdown_github/4e.3%20charts-1.png)

``` r
  Seurat:::NoAxes()
```

    ## List of 8
    ##  $ axis.line.x : list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  $ axis.line.y : list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  $ axis.title.x: list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  $ axis.title.y: list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  $ axis.text.x : list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  $ axis.text.y : list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  $ axis.ticks.x: list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  $ axis.ticks.y: list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  - attr(*, "class")= chr [1:2] "theme" "gg"
    ##  - attr(*, "complete")= logi FALSE
    ##  - attr(*, "validate")= logi TRUE

Figure 5C: average expression of meis1 and pbx3 over sc clusters
----------------------------------------------------------------

``` r
# Predicted GRN
" 6hpa_Gene:meis1 -> 24hpa_TF:meis1 -> 72hpa_TF:runx|etv1|klf9"
```

    ## [1] " 6hpa_Gene:meis1 -> 24hpa_TF:meis1 -> 72hpa_TF:runx|etv1|klf9"

``` r
" 6hpa_Gene:etv1 -> 24hpa_TF:etv1 -> 72hpa_TF:pbx3"
```

    ## [1] " 6hpa_Gene:etv1 -> 24hpa_TF:etv1 -> 72hpa_TF:pbx3"

``` r
# Define ins and outs for extraction function
Idents(pax.neural) <- "laevis_condition"
average_clusters <- t(as.data.frame(AverageExpression(pax.neural, assays = "RNA", use.counts = TRUE, features = c("meis1", "pbx3"))))
outs <- c("meis1_out", "pbx3_out")
genes <- c("meis1", "pbx3")

# Extract average expression data for plotting
seuratTOplot(average_clusters, gene = genes, outs = outs) # run to plot function (found in sourced functions)
objects <- ls(pattern=".*_out") # get all object names made in previous function
objects.list <- lapply(objects, get) # get outs
plot <- do.call(bind_rows, objects.list) # rowbind both gene data sets
# Make metdata columns
plot$treatment <- gsub(plot$condition, pattern = "RNA.", replacement = "") # split condition column
plot <- separate(plot, col = treatment, sep = "_", into = c("treament", "celltype")) # make seperate columns for treatment and celltype
plot$treament <- gsub(plot$treament, pattern = "dmso", replacement = "24hpa") # rename dmso to 24hpa

# Factor data to order
plot$gene <- factor(plot$gene, levels = c("meis1", "pbx3")) # order gene row for plot
plot$treament <- factor(plot$treament, levels = c("uninj", "24hpa", "iwr")) # order conditions
plot$celltype <- factor(plot$celltype, levels = c("Spinal.Cord.Progenitor", "Differentiating.Neuron", "Interneurons", "Vulnerable.Motor.Neuron", "Motor.Neuron", "Motor.Neuron..leptin..", "Dopaminergic.Neurons"))  

## plot average expression of genes across each cluster at each timepoint
time.plot <- filter(plot, treament %in% c("uninj", "24hpa"))

ggplot(time.plot, aes(x=treament, y=expression, group=treament, fill=celltype))+
  geom_bar(stat = "identity", position = "dodge2") +
  scale_fill_manual(values = my.color.neural) +
  theme_bw() +
  facet_grid(gene~celltype) +
  theme(legend.position = "none", text = element_text(size=20)) +
  Seurat:::RotatedAxis()
```

![](README_files/figure-markdown_github/grn-1.png)

Figure 5D: featuremap to show expression of pbx3 and meis1 across neural clusters
---------------------------------------------------------------------------------

``` r
# Plot feature plot of pbx3 and meis 1 and split by condition 
FeaturePlot(pax.neural , features = c("pbx3", "meis1"), split.by = "mock_24", order = TRUE, combine = TRUE, cols = c("grey89", "purple")) + NoAxes()
```

![](README_files/figure-markdown_github/5d-1.png)
