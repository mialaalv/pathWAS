#!/usr/local/bin/Rscript

library(gprofiler2)
library(ggplot2)
library(stringr)

efo <- "EFO_0000270"
res_file <- file.path("results", paste0(efo, ".genes.out_with_FDR_correction_filtered_cutoff.tsv"))
fig1_file <- file.path("results", paste0(efo, "_gprofiler_unfiltered_plot.jpeg"))
fig2_file <- file.path("results", paste0(efo, "_gprofiler_filtered_plot.jpeg"))

pre_out_file <- file.path("results", paste0("gprofiler2_", efo, ".tsv"))
out_file <- file.path("results", paste0(efo, "_enrichment_for_cytoscape.tsv"))

res <- read.delim(res_file)

genes <- as.character(res$GENE[order(res$ZSTAT, decreasing = T)]) # ordered list of genes by decreasing Z-score of association to disease trait


### in asthma example, none are filtered out

# perform functional enrichment of the gene list
gprofres <- gprofiler2::gost(query = genes,
                             organism = "hsapiens",
                             ordered_query = T, # if the ids are ordered by biological importance (decreasing importance == increasing p-value)
                             correction_method = "fdr",  # if commented, corrects with default method developed by g:Profile team
                             exclude_iea = T,   # excludes automatic annotations (in silico)
                             sources = c("GO:BP", "REAC"), # only searches in GO biological process and Reactome
                             evcodes = T,                  # give extra columns with list of intersections trait-gene
                             numeric_ns = "ENTREZGENE_ACC" # entry id is entrez
)

p <- gostplot(gostres = gprofres, capped = FALSE, interactive = FALSE)
ggsave(plot = p, filename = fig1_file, device = "jpeg")


# keep only the annotations to pathways that are not too big or too small
ind_size <- gprofres$result$term_size >= 5 & gprofres$result$term_size <= 350
print(paste0("Dropped ", sum(!ind_size), " items. Kept annotations to pathways with sizes 3-350 items."))
gprofres2 <- gprofres$result[which(ind_size),]

# keep only annotations with intersection >3
ind_inter <- sapply(gprofres2$intersection, function(x) length(strsplit(x, split = ",")[[1]])>=3 )
print(paste0("Dropped ", sum(!ind_inter), " items. Kept annotations to pathways with 3 or more genes in our input."))
gprofres2 <- gprofres2[which(ind_inter),]

gprofres3 <- gprofres2
gprofres3$parents <- paste(gprofres3$parents, collapse = ";")

# make table with appropriate format:
gprofres3 <- data.frame(source = gprofres3$source,term_name = gprofres3$term_name,	term_id = gprofres3$term_id, adjusted_p_value=gprofres3$p_value,	negative_log10_of_adjusted_p_value = -log(base = 10, x = gprofres3$p_value),	term_size = gprofres3$term_size,	query_size = gprofres3$query_size,	intersection_size = gprofres3$intersection_size,	effective_domain_size = gprofres3$effective_domain_size,	intersections = gprofres3$intersection)

write.table(gprofres3, pre_out_file, sep="\t", quote = F, row.names = F)

#repeat plot after filtering
mock <- list(gprofres2, gprofres$meta)
names(mock) <- c("result", "meta")
q <- gostplot(mock, capped = FALSE, interactive = FALSE)
ggsave(plot = p, filename = fig2_file, device = "jpeg")


##### make table appropriate for EnrichmentMap
cytoscape_enrichment <- data.frame(GO.ID = gprofres3$term_id, Description = gprofres3$term_name, p.Val = gprofres3$adjusted_p_value, FDR = gprofres3$adjusted_p_value, Phenotype = 1, Genes = stringr::str_replace_all(string = gprofres3$intersections, pattern = " ", replacement = ","))
write.table(cytoscape_enrichment, out_file, sep="\t", quote = F, row.names = F)



sessionInfo()
