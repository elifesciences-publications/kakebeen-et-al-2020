ATAC Figures
================
Anneke Kakebeen
13 Aug 2019

-   [Setting up project: ATAC](#setting-up-project-atac)
-   [Read in data for figures](#read-in-data-for-figures)
-   [Figure 2B](#figure-2b)
-   [Figure 2C: tss plots](#figure-2c-tss-plots)
-   [Figure 2E: \# of peaks called DA between conditions](#figure-2e-of-peaks-called-da-between-conditions)
    -   [Number of peaks](#number-of-peaks)
-   [Figure 2F: GO analysis of aggregate pax v all regions](#figure-2f-go-analysis-of-aggregate-pax-v-all-regions)
    -   [PAX &gt; All](#pax-all)
    -   [ALL v pax GO](#all-v-pax-go)
-   [Figure 2G: heat map of terms in neurogenesis](#figure-2g-heat-map-of-terms-in-neurogenesis)
-   [Figure 3A: heatmaps for all DA regions at 6hpa, 24hpa, 72hpa](#figure-3a-heatmaps-for-all-da-regions-at-6hpa-24hpa-72hpa)
    -   [6hpa Heatmap](#hpa-heatmap)
    -   [24hpa Heatmap](#hpa-heatmap-1)
    -   [72hpa Heatmap](#hpa-heatmap-2)
-   [Figure 3B: venn diagrams for \# of regions](#figure-3b-venn-diagrams-for-of-regions)
    -   [Venn diagram data set up](#venn-diagram-data-set-up)
    -   [6hpa venn diagram](#hpa-venn-diagram)
    -   [24hpa venn diagram](#hpa-venn-diagram-1)
    -   [24hpa venn diagram](#hpa-venn-diagram-2)
-   [Figure 3C: Gene ontology analysis](#figure-3c-gene-ontology-analysis)
    -   [6hpa](#hpa)
    -   [24hpa](#hpa-1)
    -   [72hpa](#hpa-2)
-   [Figure 5B: heatmap of all meis1 and pbx3 peaks](#figure-5b-heatmap-of-all-meis1-and-pbx3-peaks)
    -   [meis1 heatmap](#meis1-heatmap)
    -   [pbx3 heatmap](#pbx3-heatmap)

Setting up project: ATAC
------------------------

``` r
# Load packages
library(knitr)
library(gridExtra)
library(dplyr)
library(gprofiler2)
library(ggplot2)
library(AnnotationDbi)
library(biomaRt)
library(gplots)
library(tidyverse)
library(RColorBrewer)
library(DT)
library(LaCroixColoR)
library(Seurat)
library(VennDiagram)
library(cowplot)
library(imager)
library(magick)
library(pdftools)
library(limma)

# source functions
source("~/Desktop/R_working/common_source_functions.R")

## color palette for heatmaps
color.palette <- colorRampPalette(lacroix_palette("Lemon", n = 50, type = "continuous"))(256) 

# setwd
setwd("~/Desktop/pax6 paper/Final Markdowns/ATAC/ATAC_DA/")
```

Read in data for figures
------------------------

Figure 2B
---------

![](README_files/figure-markdown_github/2b-1.png)

Figure 2C: tss plots
--------------------

Figure 2E: \# of peaks called DA between conditions
---------------------------------------------------

### Number of peaks

![](README_files/figure-markdown_github/peaks%20between%20conditions-1.png)![](README_files/figure-markdown_github/peaks%20between%20conditions-2.png)

    ## <ggproto object: Class FacetGrid, Facet, gg>
    ##     compute_layout: function
    ##     draw_back: function
    ##     draw_front: function
    ##     draw_labels: function
    ##     draw_panels: function
    ##     finish_data: function
    ##     init_scales: function
    ##     map_data: function
    ##     params: list
    ##     setup_data: function
    ##     setup_params: function
    ##     shrink: TRUE
    ##     train_scales: function
    ##     vars: function
    ##     super:  <ggproto object: Class FacetGrid, Facet, gg>

Figure 2F: GO analysis of aggregate pax v all regions
-----------------------------------------------------

### PAX &gt; All

    ## function (data, options = list(), class = "display", callback = JS("return table;"), 
    ##     rownames, colnames, container, caption = NULL, filter = c("none", 
    ##         "bottom", "top"), escape = TRUE, style = "default", width = NULL, 
    ##     height = NULL, elementId = NULL, fillContainer = getOption("DT.fillContainer", 
    ##         NULL), autoHideNavigation = getOption("DT.autoHideNavigation", 
    ##         NULL), selection = c("multiple", "single", "none"), extensions = list(), 
    ##     plugins = NULL, editable = FALSE) 
    ## {
    ##     oop = base::options(stringsAsFactors = FALSE)
    ##     on.exit(base::options(oop), add = TRUE)
    ##     options = modifyList(getOption("DT.options", list()), if (is.function(options)) 
    ##         options()
    ##     else options)
    ##     params = list()
    ##     attr(params, "TOJSON_ARGS") = getOption("DT.TOJSON_ARGS")
    ##     if (crosstalk::is.SharedData(data)) {
    ##         params$crosstalkOptions = list(key = data$key(), group = data$groupName())
    ##         data = data$data(withSelection = FALSE, withFilter = TRUE, 
    ##             withKey = FALSE)
    ##     }
    ##     rn = if (missing(rownames) || isTRUE(rownames)) 
    ##         base::rownames(data)
    ##     else {
    ##         if (is.character(rownames)) 
    ##             rownames
    ##     }
    ##     hideDataTable = FALSE
    ##     if (is.null(data) || identical(ncol(data), 0L)) {
    ##         data = matrix(ncol = 0, nrow = NROW(data))
    ##         hideDataTable = TRUE
    ##     }
    ##     else if (length(dim(data)) != 2) {
    ##         str(data)
    ##         stop("'data' must be 2-dimensional (e.g. data frame or matrix)")
    ##     }
    ##     if (is.data.frame(data)) {
    ##         data = as.data.frame(data)
    ##         numc = unname(which(vapply(data, is.numeric, logical(1))))
    ##     }
    ##     else {
    ##         if (!is.matrix(data)) 
    ##             stop("'data' must be either a matrix or a data frame, and cannot be ", 
    ##                 classes(data), " (you may need to coerce it to matrix or data frame)")
    ##         numc = if (is.numeric(data)) 
    ##             seq_len(ncol(data))
    ##         data = as.data.frame(data)
    ##     }
    ##     if (!is.null(rn)) {
    ##         data = cbind(` ` = rn, data)
    ##         numc = numc + 1
    ##     }
    ##     if (length(numc)) {
    ##         undefined_numc = setdiff(numc - 1, classNameDefinedColumns(options, 
    ##             ncol(data)))
    ##         if (length(undefined_numc)) 
    ##             options = appendColumnDefs(options, list(className = "dt-right", 
    ##                 targets = undefined_numc))
    ##     }
    ##     if (is.null(options[["order"]])) 
    ##         options$order = list()
    ##     if (is.null(options[["autoWidth"]])) 
    ##         options$autoWidth = FALSE
    ##     if (is.null(options[["orderClasses"]])) 
    ##         options$orderClasses = FALSE
    ##     cn = base::colnames(data)
    ##     if (missing(colnames)) {
    ##         colnames = cn
    ##     }
    ##     else if (!is.null(names(colnames))) {
    ##         i = convertIdx(colnames, cn)
    ##         cn[i] = names(colnames)
    ##         colnames = cn
    ##     }
    ##     if (ncol(data) - length(colnames) == 1) 
    ##         colnames = c(" ", colnames)
    ##     if (length(colnames) && colnames[1] == " ") 
    ##         options = appendColumnDefs(options, list(orderable = FALSE, 
    ##             targets = 0))
    ##     style = match.arg(tolower(style), DTStyles())
    ##     if (style == "bootstrap") 
    ##         class = DT2BSClass(class)
    ##     if (style != "default") 
    ##         params$style = style
    ##     if (isTRUE(fillContainer)) 
    ##         class = paste(class, "fill-container")
    ##     if (is.character(filter)) 
    ##         filter = list(position = match.arg(filter))
    ##     filter = modifyList(list(position = "none", clear = TRUE, 
    ##         plain = FALSE), filter)
    ##     filterHTML = as.character(filterRow(data, !is.null(rn) && 
    ##         colnames[1] == " ", filter))
    ##     if (filter$position == "top") 
    ##         options$orderCellsTop = TRUE
    ##     params$filter = filter$position
    ##     if (filter$position != "none") 
    ##         params$filterHTML = filterHTML
    ##     if (missing(container)) {
    ##         container = tags$table(tableHeader(colnames, escape), 
    ##             class = class)
    ##     }
    ##     else {
    ##         params$class = class
    ##     }
    ##     attr(options, "escapeIdx") = escapeToConfig(escape, colnames)
    ##     if (is.list(extensions)) {
    ##         extensions = names(extensions)
    ##     }
    ##     else if (!is.character(extensions)) {
    ##         stop("'extensions' must be either a character vector or a named list")
    ##     }
    ##     params$extensions = if (length(extensions)) 
    ##         as.list(extensions)
    ##     if ("Responsive" %in% extensions) 
    ##         options$responsive = TRUE
    ##     params$caption = captionString(caption)
    ##     if (isTRUE(editable)) 
    ##         editable = "cell"
    ##     if (is.character(editable)) 
    ##         editable = list(target = editable, disable = list(columns = NULL))
    ##     if (is.list(editable)) 
    ##         params$editable = editable
    ##     if (!identical(class(callback), class(JS("")))) 
    ##         stop("The 'callback' argument only accept a value returned from JS()")
    ##     if (length(options$pageLength) && length(options$lengthMenu) == 
    ##         0) {
    ##         if (!isFALSE(options$lengthChange)) 
    ##             options$lengthMenu = sort(unique(c(options$pageLength, 
    ##                 10, 25, 50, 100)))
    ##         if (identical(options$lengthMenu, c(10, 25, 50, 100))) 
    ##             options$lengthMenu = NULL
    ##     }
    ##     if (!is.null(fillContainer)) 
    ##         params$fillContainer = fillContainer
    ##     if (!is.null(autoHideNavigation)) 
    ##         params$autoHideNavigation = autoHideNavigation
    ##     params = structure(modifyList(params, list(data = data, container = as.character(container), 
    ##         options = options, callback = if (!missing(callback)) JS("function(table) {", 
    ##             callback, "}"))), colnames = cn, rownames = length(rn) > 
    ##         0)
    ##     if (inShiny() || length(params$crosstalkOptions)) {
    ##         if (is.character(selection)) {
    ##             selection = list(mode = match.arg(selection))
    ##         }
    ##         selection = modifyList(list(mode = "multiple", selected = NULL, 
    ##             target = "row"), selection)
    ##         if (grepl("^row", selection$target) && is.character(selection$selected) && 
    ##             length(rn)) {
    ##             selection$selected = match(selection$selected, rn)
    ##         }
    ##         params$selection = selection
    ##     }
    ##     deps = list(DTDependency(style))
    ##     deps = c(deps, unlist(lapply(extensions, extDependency, style, 
    ##         options), recursive = FALSE))
    ##     if (params$filter != "none") 
    ##         deps = c(deps, filterDependencies())
    ##     if (isTRUE(options$searchHighlight)) 
    ##         deps = c(deps, list(pluginDependency("searchHighlight")))
    ##     if (length(plugins)) 
    ##         deps = c(deps, lapply(plugins, pluginDependency))
    ##     deps = c(deps, crosstalk::crosstalkLibs())
    ##     if (isTRUE(fillContainer)) {
    ##         width = NULL
    ##         height = NULL
    ##     }
    ##     htmlwidgets::createWidget("datatables", if (hideDataTable) 
    ##         NULL
    ##     else params, package = "DT", width = width, height = height, 
    ##         elementId = elementId, sizingPolicy = htmlwidgets::sizingPolicy(knitr.figure = FALSE, 
    ##             knitr.defaultWidth = "100%", knitr.defaultHeight = "auto"), 
    ##         dependencies = deps, preRenderHook = function(instance) {
    ##             data = instance[["x"]][["data"]]
    ##             if (object.size(data) > 1500000 && getOption("DT.warn.size", 
    ##                 TRUE)) 
    ##                 warning("It seems your data is too big for client-side DataTables. You may ", 
    ##                   "consider server-side processing: https://rstudio.github.io/DT/server.html")
    ##             data = escapeData(data, escape, colnames)
    ##             data = unname(data)
    ##             instance$x$data = data
    ##             instance
    ##         })
    ## }
    ## <bytecode: 0x7fa4337a7fd0>
    ## <environment: namespace:DT>

![](README_files/figure-markdown_github/2f.1%20-1.png)

### ALL v pax GO

Goal : Show that regions that are more accessible in pax libraries vs all-tissue libraries have neural character Method: Identify differential regions of the chromatin between Pax and AT, use these regions for GO Conclusions: Regions prioritized in all tissue show more generic GO terms. ![](README_files/figure-markdown_github/ALL%20-1.png)

Figure 2G: heat map of terms in neurogenesis
--------------------------------------------

    ##  [1] "All_Tissue_0hpa"      "All_Tissue_24hpa"     "All_Tissue_6hpa"     
    ##  [4] "All_Tissue_72hpa"     "All_Tissue_uninjured" "Pax_0hpa"            
    ##  [7] "Pax_24hpa"            "Pax_6hpa"             "Pax_72hpa"           
    ## [10] "Pax_uninjured"

![](README_files/figure-markdown_github/PaxVAll-1.png)

Figure 3A: heatmaps for all DA regions at 6hpa, 24hpa, 72hpa
------------------------------------------------------------

### 6hpa Heatmap

![](README_files/figure-markdown_github/6hpa%20heat-1.png)

### 24hpa Heatmap

![](README_files/figure-markdown_github/24hpa%20heat-1.png)

### 72hpa Heatmap

![](README_files/figure-markdown_github/72hpa%20heat-1.png)

Figure 3B: venn diagrams for \# of regions
------------------------------------------

### Venn diagram data set up

### 6hpa venn diagram

![](README_files/figure-markdown_github/3b.2-1.png)

    ## (polygon[GRID.polygon.202], polygon[GRID.polygon.203], polygon[GRID.polygon.204], polygon[GRID.polygon.205], text[GRID.text.206], text[GRID.text.207], text[GRID.text.208], text[GRID.text.209], text[GRID.text.210])

### 24hpa venn diagram

![](README_files/figure-markdown_github/3b.3-1.png)

    ## (polygon[GRID.polygon.211], polygon[GRID.polygon.212], polygon[GRID.polygon.213], polygon[GRID.polygon.214], text[GRID.text.215], text[GRID.text.216], text[GRID.text.217], text[GRID.text.218], text[GRID.text.219])

### 24hpa venn diagram

![](README_files/figure-markdown_github/3b.4-1.png)

    ## (polygon[GRID.polygon.220], polygon[GRID.polygon.221], polygon[GRID.polygon.222], polygon[GRID.polygon.223], text[GRID.text.224], text[GRID.text.225], text[GRID.text.226], text[GRID.text.227], text[GRID.text.228])

Figure 3C: Gene ontology analysis
---------------------------------

### 6hpa

#### GO

#### ReviGO

    ##       term_id      p_value
    ## 1  GO:0030182 2.548818e-10
    ## 2  GO:0022008 3.519637e-10
    ## 3  GO:0048699 1.888581e-09
    ## 4  GO:0007399 9.589671e-08
    ## 5  GO:0048666 3.986622e-07
    ## 6  GO:0048522 6.467150e-07
    ## 7  GO:0045935 1.717019e-06
    ## 8  GO:0009653 1.739057e-06
    ## 9  GO:0031175 5.886463e-06
    ## 10 GO:0048518 6.027687e-06
    ## 11 GO:0031328 6.899671e-06
    ## 12 GO:0045893 8.105129e-06
    ## 13 GO:0048869 1.127867e-05
    ## 14 GO:0010557 1.176648e-05
    ## 15 GO:0009891 2.011116e-05
    ## 16 GO:0071840 2.164779e-05
    ## 17 GO:0016043 2.249058e-05
    ## 18 GO:1903508 3.579506e-05
    ## 19 GO:1902680 3.695455e-05
    ## 20 GO:0032990 4.476258e-05
    ## 21 GO:0048523 4.959781e-05
    ## 22 GO:0009888 6.920488e-05
    ## 23 GO:0051254 7.254999e-05
    ## 24 GO:0000902 7.309745e-05
    ## 25 GO:0120039 8.879869e-05
    ## 26 GO:0032989 9.994818e-05
    ## 27 GO:0048858 1.081515e-04
    ## 28 GO:0048519 1.103739e-04
    ## 29 GO:0048812 1.103990e-04
    ## 30 GO:0050793 1.116761e-04
    ## 31 GO:0010628 1.136861e-04
    ## 32 GO:0030154 1.413590e-04
    ## 33 GO:0032502 1.542352e-04
    ## 34 GO:0007275 1.609296e-04
    ## 35 GO:0051173 1.893469e-04
    ## 36 GO:0051094 1.971396e-04
    ## 37 GO:0051239 2.014401e-04
    ## 38 GO:0048731 2.175160e-04
    ## 39 GO:0048856 2.414275e-04
    ## 40 GO:0010604 2.772573e-04
    ## 41 GO:0120036 3.685297e-04
    ## 42 GO:0048468 4.544928e-04
    ## 43 GO:0023051 5.958513e-04
    ## 44 GO:0060322 6.410714e-04
    ## 45 GO:0045944 6.530781e-04
    ## 46 GO:0010646 6.918512e-04
    ## 47 GO:0019219 9.906617e-04
    ## 48 GO:0030030 1.098061e-03
    ## 49 GO:2000026 1.107067e-03
    ## 50 GO:0031325 1.419093e-03
    ## 51 GO:0009893 1.433490e-03
    ## 52 GO:0051240 1.643579e-03
    ## 53 GO:0080090 1.718099e-03
    ## 54 GO:0045595 1.862481e-03
    ## 55 GO:0048646 2.251078e-03
    ## 56 GO:0009790 2.385804e-03
    ## 57 GO:0043009 2.403900e-03
    ## 58 GO:0007420 2.594750e-03
    ## 59 GO:0009792 5.840858e-03
    ## 60 GO:0007417 6.104407e-03
    ## 61 GO:0032879 6.747124e-03
    ## 62 GO:0000904 7.190087e-03
    ## 63 GO:0051252 7.741704e-03
    ## 64 GO:0072088 8.048957e-03
    ## 65 GO:0050767 8.065112e-03
    ## 66 GO:0050789 8.254871e-03
    ## 67 GO:0051960 9.003299e-03
    ## 68 GO:0072028 1.081568e-02
    ## 69 GO:0048667 1.292215e-02
    ## 70 GO:0019438 1.615936e-02
    ## 71 GO:0045597 1.654268e-02
    ## 72 GO:0030900 1.724287e-02
    ## 73 GO:0065007 1.817771e-02
    ## 74 GO:0048568 1.904146e-02
    ## 75 GO:2000112 2.028119e-02
    ## 76 GO:0006355 2.112200e-02
    ## 77 GO:0045664 2.263485e-02
    ## 78 GO:0072009 2.266228e-02
    ## 79 GO:0060284 2.409965e-02
    ## 80 GO:0034654 2.428372e-02
    ## 81 GO:0022603 2.567728e-02
    ## 82 GO:0009887 2.605589e-02
    ## 83 GO:0031326 2.865658e-02
    ## 84 GO:0061564 2.995153e-02
    ## 85 GO:0009966 3.127162e-02
    ## 86 GO:0018130 3.144557e-02
    ## 87 GO:0009889 3.202411e-02
    ## 88 GO:1901362 3.302883e-02
    ## 89 GO:0019222 3.348661e-02
    ## 90 GO:0072078 3.352561e-02
    ## 91 GO:1903506 3.752556e-02
    ## 92 GO:0090304 3.926926e-02
    ## 93 GO:2001141 4.073441e-02
    ## 94 GO:0072006 4.195966e-02
    ## 95 GO:0045934 4.430092e-02

![](README_files/figure-markdown_github/3c.2%206hpa-1.png)![](README_files/figure-markdown_github/3c.2%206hpa-2.png)

### 24hpa

#### GO

#### ReviGO

    ##       term_id      p_value
    ## 1  GO:0045778 0.0007202793
    ## 2  GO:0045669 0.0008388526
    ## 3  GO:0071495 0.0047108171
    ## 4  GO:0110151 0.0053909212
    ## 5  GO:0070169 0.0053909212
    ## 6  GO:0110149 0.0070897482
    ## 7  GO:0070167 0.0070897482
    ## 8  GO:0048645 0.0083406265
    ## 9  GO:0045893 0.0084430792
    ## 10 GO:0071310 0.0095625655
    ## 11 GO:0045944 0.0106831618
    ## 12 GO:0045935 0.0114157045
    ## 13 GO:0070887 0.0124694472
    ## 14 GO:0000902 0.0130013493
    ## 15 GO:0048812 0.0152597294
    ## 16 GO:0032989 0.0210264369
    ## 17 GO:0120039 0.0242644462
    ## 18 GO:0030278 0.0262590730
    ## 19 GO:1903508 0.0271850996
    ## 20 GO:0001649 0.0274044569
    ## 21 GO:0048858 0.0276176631
    ## 22 GO:1902680 0.0277409428
    ## 23 GO:0045667 0.0283729861
    ## 24 GO:0048729 0.0294494743
    ## 25 GO:0048666 0.0364367807
    ## 26 GO:0010033 0.0406704796
    ## 27 GO:0007399 0.0408672272
    ## 28 GO:0051173 0.0460073312
    ## 29 GO:0003002 0.0476348551
    ## 30 GO:0001503 0.0485672526

![](README_files/figure-markdown_github/3c.2%2024hpa-1.png)![](README_files/figure-markdown_github/3c.2%2024hpa-2.png)

### 72hpa

#### GO

#### ReviGO

    ##       term_id      p_value
    ## 1  GO:0009653 3.398761e-07
    ## 2  GO:0009790 1.044657e-06
    ## 3  GO:0009887 1.457045e-06
    ## 4  GO:0043009 1.082296e-05
    ## 5  GO:0009792 2.501293e-05
    ## 6  GO:0045944 1.216410e-04
    ## 7  GO:0010557 2.193700e-04
    ## 8  GO:0048568 3.055113e-04
    ## 9  GO:0001822 4.518147e-04
    ## 10 GO:1903508 4.740584e-04
    ## 11 GO:1902680 4.844947e-04
    ## 12 GO:0048513 5.104617e-04
    ## 13 GO:0001655 5.391639e-04
    ## 14 GO:0031328 6.303260e-04
    ## 15 GO:0045935 7.194544e-04
    ## 16 GO:0045893 8.085233e-04
    ## 17 GO:0051239 8.821247e-04
    ## 18 GO:0072001 1.127397e-03
    ## 19 GO:0010628 1.140531e-03
    ## 20 GO:0050793 1.185661e-03
    ## 21 GO:0051094 1.195827e-03
    ## 22 GO:0009891 1.252811e-03
    ## 23 GO:0001501 1.283131e-03
    ## 24 GO:0051254 1.340510e-03
    ## 25 GO:0009888 1.490622e-03
    ## 26 GO:2000026 1.652732e-03
    ## 27 GO:0071407 1.715271e-03
    ## 28 GO:0040019 2.031548e-03
    ## 29 GO:0048706 2.723179e-03
    ## 30 GO:0003013 3.808656e-03
    ## 31 GO:0010604 4.066619e-03
    ## 32 GO:1901700 4.078998e-03
    ## 33 GO:0007423 4.241169e-03
    ## 34 GO:0048731 4.610722e-03
    ## 35 GO:0071495 5.740837e-03
    ## 36 GO:0008015 8.248609e-03
    ## 37 GO:0048705 1.364965e-02
    ## 38 GO:0007507 1.395906e-02
    ## 39 GO:0048598 1.537970e-02
    ## 40 GO:0051240 1.642779e-02
    ## 41 GO:0072359 1.687225e-02
    ## 42 GO:0031325 2.421851e-02
    ## 43 GO:0009719 2.463327e-02
    ## 44 GO:0061351 2.541712e-02
    ## 45 GO:2000177 2.645629e-02
    ## 46 GO:0009893 2.720525e-02
    ## 47 GO:0010033 2.985669e-02
    ## 48 GO:0007275 3.074801e-02
    ## 49 GO:0048704 3.212064e-02
    ## 50 GO:1901701 3.335811e-02
    ## 51 GO:0048593 3.420338e-02
    ## 52 GO:0051173 3.501576e-02
    ## 53 GO:0071310 3.654910e-02
    ## 54 GO:0090596 3.690413e-02
    ## 55 GO:0043010 4.235992e-02
    ## 56 GO:0048646 4.291354e-02
    ## 57 GO:0048522 4.422501e-02
    ## 58 GO:0045165 4.715255e-02
    ## 59 GO:0003002 4.882245e-02
    ## 60 GO:0032502 4.951174e-02

![](README_files/figure-markdown_github/3c.2%2072hpa-1.png)![](README_files/figure-markdown_github/3c.2%2072hpa-2.png)

Figure 5B: heatmap of all meis1 and pbx3 peaks
----------------------------------------------

### meis1 heatmap

    ##  [1] "All_Tissue_0hpa"      "All_Tissue_24hpa"     "All_Tissue_6hpa"     
    ##  [4] "All_Tissue_72hpa"     "All_Tissue_uninjured" "Pax_0hpa"            
    ##  [7] "Pax_24hpa"            "Pax_6hpa"             "Pax_72hpa"           
    ## [10] "Pax_uninjured"

![](README_files/figure-markdown_github/5b%20meis-1.png)

    ## [1] "Pax_uninjured" "Pax_0hpa"      "Pax_6hpa"      "Pax_24hpa"    
    ## [5] "Pax_72hpa"

![](README_files/figure-markdown_github/5b%20meis-2.png)

### pbx3 heatmap

    ##  [1] "All_Tissue_0hpa"      "All_Tissue_24hpa"     "All_Tissue_6hpa"     
    ##  [4] "All_Tissue_72hpa"     "All_Tissue_uninjured" "Pax_0hpa"            
    ##  [7] "Pax_24hpa"            "Pax_6hpa"             "Pax_72hpa"           
    ## [10] "Pax_uninjured"

![](README_files/figure-markdown_github/5b%20pbx3-1.png)

    ## [1] "Pax_uninjured" "Pax_0hpa"      "Pax_6hpa"      "Pax_24hpa"    
    ## [5] "Pax_72hpa"

![](README_files/figure-markdown_github/5b%20pbx3-2.png)
