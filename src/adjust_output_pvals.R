# get stats

efo <- "EFO_0000270"
res_file <- file.path("results", paste0(efo, ".genes.out"))
out_file <- file.path("results", paste0(efo, ".genes.out_with_FDR_correction.tsv"))

res <- read.delim(res_file, sep="")
res$P_adjust <- p.adjust(p = res$P, method = "fdr")

write.table(x = res, quote = F, file = out_file, sep="\t", row.names = F)

sessionInfo()

