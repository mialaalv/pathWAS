#!/usr/local/bin/Rscript

if(!require(gwasrapidd)) remotes::install_github("ramiromagno/gwasrapidd")

if(!require(stringr)) install.packages("stringr")

if(!require(rtracklayer)) BiocManager::install("rtracklayer")
if(!require(GenomicRanges)) BiocManager::install("GenomicRanges")

sessionInfo()

