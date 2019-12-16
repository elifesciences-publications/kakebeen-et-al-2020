FACS Analysis of Neural Progenitor Cells Over Regeneration
================
Anneke Kakebeen
11/15/2019

-   [Set up](#set-up)
-   [Flowjo Analysis](#flowjo-analysis)
    -   [Read in all FACS data](#read-in-all-facs-data)
    -   [ATAC-Seq samples and Most recent FACS analysis for higher N](#atac-seq-samples-and-most-recent-facs-analysis-for-higher-n)
        -   [Subset table](#subset-table)
        -   [Plot table](#plot-table)

Reviewer 3's minor comment: "How does the number of pax6+ cells change during regeneration? According to the authors' conclusions, the number of NPCs would decrease at 6 hpa and increase at 72 hpa. A supplemental figure showing the number of GFP+ cells sorted and a description of how many animals were used to collect that many cells would also support this claim."

In response, I elected to analyze and plot existing data from the ATAC-Seq collections and most recent collections related to a different project, but including necessary samples.

Set up
======

Flowjo Analysis
===============

All pax regeneration samples that were run on FACS have been re-analyzed on flowjo to idenitfy more accurate numbers sorted (does not relate to number of cells collected).

Read in all FACS data
---------------------

ATAC-Seq samples and Most recent FACS analysis for higher N
-----------------------------------------------------------

These plots take into account each biological replicate sorted for atac-seq libraries as well as biological replicates measured in a different experiment performed summer of 2019.

### Subset table

Based on our ATAC-Seq analysis and scRNA-Seq analysis, we hypothesized that the 72hpa timepoint has more proliferative neural progenitor cells than at 6hpa or 24hpa. To quantify the number of cells per regenerating tail, we identified the GFP+ cells in each flow cytometry analysis and divided this number by the number of input tails.

### Plot table

![](README_files/figure-markdown_github/flowjo%20output%20plot-1.png)

From this data, we conclude that there is an increase in the number of GFP+ cells per tail detectet at 72hpa compared to 6hpa and 24hpa.
