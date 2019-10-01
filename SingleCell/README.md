Single Cell Figures
================
Anneke Kakebeen
13 Aug 2019

-   [Set up for scRNA-Seq Figures](#set-up-for-scrna-seq-figures)
-   [Read in scRNASeq Objects](#read-in-scrnaseq-objects)
-   [Figure 4a: UMAP of all neural cells](#figure-4a-umap-of-all-neural-cells)
-   [Figure 4b: UMAP of all neural cells colored for condition (uninjured or 24hpa)](#figure-4b-umap-of-all-neural-cells-colored-for-condition-uninjured-or-24hpa)
-   [Figure 4C: Marker dotplot for neural cell clusters](#figure-4c-marker-dotplot-for-neural-cell-clusters)
-   [Figure 4D-G: cell cluster proportions](#figure-4d-g-cell-cluster-proportions)
    -   [Set up data table](#set-up-data-table)
    -   [Pie Charts](#pie-charts)
    -   [UMAPS](#umaps)
-   [Figure 4E: Cell cycle phase prediction](#figure-4e-cell-cycle-phase-prediction)
    -   [Table set up](#table-set-up)
    -   [Pie Charts](#pie-charts-1)
    -   [UMAP](#umap)
-   [Figure 5E: average expression of meis1 and pbx3 over sc clusters](#figure-5e-average-expression-of-meis1-and-pbx3-over-sc-clusters)
-   [Figure 5F: sunburst of number of cells that express pbx3 and meis 1](#figure-5f-sunburst-of-number-of-cells-that-express-pbx3-and-meis-1)

Set up for scRNA-Seq Figures
----------------------------

``` r
# load packages
library(Seurat)
library(LaCroixColoR)
library(ggplot2)
library(cowplot)
library(tidyr)
library(DT)
library(dplyr)
library(ggpubr)
library(rmarkdown)
library(knitr)

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

Figure 4a: UMAP of all neural cells
-----------------------------------

![](README_files/figure-markdown_github/4a-1.png)

Figure 4b: UMAP of all neural cells colored for condition (uninjured or 24hpa)
------------------------------------------------------------------------------

![](README_files/figure-markdown_github/4b-1.png)

Figure 4C: Marker dotplot for neural cell clusters
--------------------------------------------------

![](README_files/figure-markdown_github/4C-1.png)

Figure 4D-G: cell cluster proportions
-------------------------------------

### Set up data table

### Pie Charts

![](README_files/figure-markdown_github/4f/g%20uninjured%20cell%20proportions-1.png)

    ## # A tibble: 14 x 3
    ## # Groups:   condition [2]
    ##    condition cluster                    sum
    ##    <fct>     <fct>                    <dbl>
    ##  1 uninj     Spinal Cord Progenitor  38.5  
    ##  2 uninj     Differentiating Neuron  17.9  
    ##  3 uninj     Interneurons             5.61 
    ##  4 uninj     Vulnerable Motor Neuron 25.3  
    ##  5 uninj     Dopaminergic Neurons     5.99 
    ##  6 uninj     Motor Neuron (leptin+)   0.255
    ##  7 uninj     Motor Neuron             6.51 
    ##  8 dmso      Spinal Cord Progenitor  25.3  
    ##  9 dmso      Differentiating Neuron  10.1  
    ## 10 dmso      Interneurons            16.9  
    ## 11 dmso      Vulnerable Motor Neuron 20.9  
    ## 12 dmso      Dopaminergic Neurons     2.03 
    ## 13 dmso      Motor Neuron (leptin+)  21.3  
    ## 14 dmso      Motor Neuron             3.38

### UMAPS

![](README_files/figure-markdown_github/4d/e%20uninjured%20cell%20proportions%202-1.png)

Figure 4E: Cell cycle phase prediction
--------------------------------------

### Table set up

### Pie Charts

![](README_files/figure-markdown_github/4j/k-1.png)

    ## # A tibble: 6 x 3
    ## # Groups:   condition [2]
    ##   condition phase   sum
    ##   <fct>     <fct> <dbl>
    ## 1 uninj     S     30.4 
    ## 2 uninj     G2M   17.0 
    ## 3 uninj     G1    52.7 
    ## 4 dmso      S     16.2 
    ## 5 dmso      G2M    7.77
    ## 6 dmso      G1    76.0

### UMAP

![](README_files/figure-markdown_github/4h/i%20charts-1.png)

Figure 5E: average expression of meis1 and pbx3 over sc clusters
----------------------------------------------------------------

    ## [1] " 6hpa_Gene:meis1 -> 24hpa_TF:meis1 -> 72hpa_TF:runx|etv1|klf9"

    ## [1] " 6hpa_Gene:etv1 -> 24hpa_TF:etv1 -> 72hpa_TF:pbx3"

![](README_files/figure-markdown_github/grn-1.png)

Figure 5F: sunburst of number of cells that express pbx3 and meis 1
-------------------------------------------------------------------

![](README_files/figure-markdown_github/5f-1.png)![](README_files/figure-markdown_github/5f-2.png)
