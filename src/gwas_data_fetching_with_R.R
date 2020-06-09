#!/usr/local/bin/Rscript


library(gwasrapidd)
library(stringr)
library(rtracklayer)

##### get associations with trait asthma and its child terms:
# efo <- get_traits(efo_trait = "asthma")@traits$efo_id
efo <- "EFO_0000270" # asthma, used as example

##### path variables

associations_file <- file.path("data", "GWAS_data", efo, paste0(efo, "_associations.rds"))
studies_file <- file.path("data", "GWAS_data", efo, paste0(efo, "_studies.rds"))
variants_file <- file.path("data", "GWAS_data", efo, paste0(efo, "_variants.rds"))
traits_file <- file.path("data", "GWAS_data", efo, paste0(efo, "_traits.rds"))

summary_file <- file.path("data", "GWAS_data", efo, paste0(efo, "_GWAS_summary.rds"))
for_magma_file <- file.path("data", "GWAS_data", efo, paste0(efo, "_GWAS_summary_grch37_forMAGMA.tsv"))

##### get data from GWAS Catalog

efos <- c(efo, get_child_efo(efo_id = efo)[[1]])

my_traits <- get_traits(efo_id = efos)

# summary of the traits considered:
my_traits@traits
# # A tibble: 5 x 3
#   efo_id      trait                      uri                                 
#   <chr>       <chr>                      <chr>                               
# 1 EFO_0000270 asthma                     http://www.ebi.ac.uk/efo/EFO_0000270
# 2 EFO_0010638 atopic asthma              http://www.ebi.ac.uk/efo/EFO_0010638
# 3 EFO_0009759 Chronic Obstructive Asthma http://www.ebi.ac.uk/efo/EFO_0009759
# 4 EFO_1002011 adult onset asthma         http://www.ebi.ac.uk/efo/EFO_1002011
# 5 EFO_0004591 childhood onset asthma     http://www.ebi.ac.uk/efo/EFO_0004591

my_associations <- get_associations(efo_id = efos) # 2363 associations

association_ids <- my_associations@associations$association_id

my_variants <- get_variants(association_id = association_ids) # n 1541
my_studies <- get_studies(association_id = association_ids) # n 102

studies_from_associations <- sapply(association_ids, function(x) {
  a <- get_studies(association_id = x)@studies$study_id
  # if(length(a) > 1) print(paste0("looko at ", x, ":", a))
  return(a)
}) # checked that there's only 1 study per association

variants_from_associations <- sapply(association_ids, function(x) {
  a <- get_variants(association_id = x)@variants$variant_id
  # if(length(a) > 1) print(paste0("look at ", x, ":", a))
  return(a)
}) # returns a list with variants per association

# in case there is more than one variant per association, the rows are duplicated

# I only take the initial number of individuals (not replicas), ancestry_id = 1
nobs <- dplyr::filter(my_studies[studies_from_associations,]@ancestries, ancestry_id == 1)$number_of_individuals

df <- data.frame(association_id = association_ids, study = studies_from_associations, NOBS = nobs)

n_vars <- sapply(variants_from_associations, length) # number of variants per accession
df <- df[rep(association_ids, n_vars),]              # repeat each row according to the number of variants

df$SNP <- unlist(variants_from_associations)        # add the variants to the df

df$P <- my_associations[df$association_id,]@associations$pvalue
df$CHR38 <- my_variants[df$SNP,]@variants$chromosome_name
df$BP38 <- my_variants[df$SNP,]@variants$chromosome_position

# some variants do not have specified positions, but they are specified in the name of the variant. it can be parsed to rescue them:
# for the df rows with na in chr
for(x in 1:nrow(df)){
  if(is.na(df[x,"CHR38"]) & stringr::str_starts(string = df[x, "SNP"], pattern = "ch")){
    capture <- stringr::str_match(df[x, "SNP"], pattern = "ch[r]?(.*):([0-9]*)")

    df$CHR38[x] <- capture[1,2]
    df$BP38[x] <- capture[1,3]
  }
}

# now there are only two variants with no chromosome position annotations:
df[which(is.na(df$CHR38)),"SNP"]
# [1] "rs152271219" "rs67431028" 

# they are dropped
df <- df[complete.cases(df),]


####### Lift genomic positions to GRCH38:

path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)

ranges38 <- GRanges(seqnames = paste0("chr", df$CHR38), ranges = df$BP)
ranges37 <- liftOver(x = ranges38, chain = ch)

df$CHR <- stringr::str_replace(string = as.character(ranges37@unlistData@seqnames),pattern = "chr", replacement = "")
df$BP <- as.character(ranges37@unlistData@ranges@start)



###### saving results and data for future reference #######

dir.create(file.path("data", "GWAS_data", efo), recursive = T)

saveRsink(file.path("logs", "gwas_data_fetching_with_R_log.txt"))
DS(my_associations, associations_file)
saveRDS(my_studies, studies_file)
saveRDS(my_variants, variants_file)
saveRDS(my_traits, traits_file)

saveRDS(df, summary_file)

# rearrange columns to meet requested order: SNP CHR BP P NOBS
write.table(x = df[, c("SNP", "CHR", "BP", "P", "NOBS")], file = for_magma_file, sep="\t", quote=F, row.names = F, col.names = T)



##### Session Info
sessionInfo()


# R version 4.0.0 (2020-04-24)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 19.10
# 
# Matrix products: default
# BLAS:   /opt/R/4.0.0/lib/R/lib/libRblas.so
# LAPACK: /opt/R/4.0.0/lib/R/lib/libRlapack.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=es_ES.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=es_ES.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=es_ES.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=es_ES.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] graph_1.66.0         rtracklayer_1.48.0   GenomicRanges_1.40.0 GenomeInfoDb_1.24.0  IRanges_2.22.2       S4Vectors_0.26.1    
# [7] BiocGenerics_0.34.0  gwascat_2.20.1       gwasrapidd_0.99.8   
# 
# loaded via a namespace (and not attached):
#   [1] Biobase_2.48.0                          httr_1.4.1                              tidyr_1.1.0                            
# [4] bit64_0.9-7                             jsonlite_1.6.1                          assertthat_0.2.1                       
# [7] askpass_1.1                             BiocManager_1.30.10                     triebeard_0.3.0                        
# [10] urltools_1.7.3                          BiocFileCache_1.12.0                    blob_1.2.1                             
# [13] GenomeInfoDbData_1.2.3                  Rsamtools_2.4.0                         remotes_2.1.1                          
# [16] progress_1.2.2                          pillar_1.4.4                            RSQLite_2.2.0                          
# [19] lattice_0.20-41                         glue_1.4.1                              digest_0.6.25                          
# [22] XVector_0.28.0                          colorspace_1.4-1                        Matrix_1.2-18                          
# [25] plyr_1.8.6                              XML_3.99-0.3                            pkgconfig_2.0.3                        
# [28] biomaRt_2.44.0                          zlibbioc_1.34.0                         purrr_0.3.4                            
# [31] scales_1.1.1                            BiocParallel_1.22.0                     tibble_3.0.1                           
# [34] openssl_1.4.1                           generics_0.0.2                          ggplot2_3.3.1                          
# [37] ellipsis_0.3.1                          SummarizedExperiment_1.18.1             GenomicFeatures_1.40.0                 
# [40] cli_2.0.2                               magrittr_1.5                            crayon_1.3.4                           
# [43] memoise_1.1.0                           fansi_0.4.1                             tools_4.0.0                            
# [46] prettyunits_1.1.1                       hms_0.5.3                               lifecycle_0.2.0                        
# [49] matrixStats_0.56.0                      stringr_1.4.0                           munsell_0.5.0                          
# [52] DelayedArray_0.14.0                     AnnotationDbi_1.50.0                    Biostrings_2.56.0                      
# [55] compiler_4.0.0                          rlang_0.4.6                             grid_4.0.0                             
# [58] RCurl_1.98-1.2                          rstudioapi_0.11                         rappdirs_0.3.1                         
# [61] bitops_1.0-6                            gtable_0.3.0                            DBI_1.1.0                              
# [64] curl_4.3                                R6_2.4.1                                GenomicAlignments_1.24.0               
# [67] lubridate_1.7.8                         dplyr_1.0.0                             bit_1.1-15.2                           
# [70] utf8_1.1.4                              TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 stringi_1.4.6                          
# [73] Rcpp_1.0.4.6                            vctrs_0.3.0                             dbplyr_1.4.4                           
# [76] tidyselect_1.1.0    