#!/usr/local/bin/Rscript

#from github
if(!require(gwasrapidd)) remotes::install_github("ramiromagno/gwasrapidd")
#from bioconductor
if(!require(rtracklayer)) BiocManager::install("rtracklayer")
if(!require(GenomicRanges)) BiocManager::install("GenomicRanges")
if(!require(limma)) BiocManager::install("limma")
# from cran
if(!require(GSA)) install.packages("GSA")
if(!require(RCurl)) install.packages("RCurl")
if(!require(stringr)) install.packages("stringr")
if(!require(gprofiler2)) install.packages("gprofiler2")
if(!require(ggplot2)) install.packages("ggplot2")

sessionInfo()

