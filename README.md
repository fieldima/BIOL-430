# README

## Background
In this project, I will be performing a differential expression meta-analysis on Mice infected by WNV.
Specifically, I will be looking for genes that are uniquely expressed upon infection by infective WNV
but are common across different neuronal tissue types.

I will be including the following experiments from the Gene Expression Omnibus into my analysis: 
GSE78888, GSE77193, GSE68380, GSE77192, and GSE67473.

To do this, I will first perform DEA on all the experiments individually, comparing WNV infection to
control condition for both WT and the deletion strain E218A (which attenuates WNV) at a single timepoint
of 24 hrs post infection.

I will then perform meta analysis using MetaVolcanoR with included experiments for WT WNV infection and
then perform a separate meta analysis with E218A WNV infection. I will compare the hit lists to identify
genes specific to the infection response but not inflammation in neuronal tissue as well as one 
lymph node.

By performing a meta-analysis I hope to identify genes with strong evidence differential expression upon
WNV infection brain-wide. 

## Directory Structure

├── BIOL-430.Rproj
├── Figures
│   ├── Heatmap.png # Heatmap of sample-sample correlations
│   ├── PCA.png # Principal Component Analysis of samples
│   ├── VennDiagram.png # Overlap of DEGs between tissue types
│   └── VolcanoPlot.png # Volcano plots of WT and #E218A DEGs
├── Outputs
│   ├── E218A_GO_Terms.csv # Gene Ontology terms for just E218A strain WNV
│   ├── Shared_GO_Terms.csv # Gene Ontology terms shared between strains
│   └── WT_GO_Terms.csv # Gene Ontology terms for just WT WNV
├── README.md
└── Scripts
    ├── Analysis.R # Script containing all analyses
    └── Install_Packages.R # Script to install dependencies