# get stats

efo <- "EFO_0000270"
res_file <- file.path("results", paste0(efo, ".genes.out"))
out_file <- file.path("results", paste0(efo, ".genes.out_with_FDR_correction.tsv"))
out_filtered_file <- file.path("results", paste0(efo, ".genes.out_with_FDR_correction_filtered_cutoff.tsv"))
list_entrezs_file <- file.path("results", paste0(efo, "_entrezs_list.txt"))


res <- read.delim(res_file, sep="")
res$P_adjust <- p.adjust(p = res$P, method = "fdr")

write.table(x = res, quote = F, file = out_file, sep="\t", row.names = F)

png("results/analisis_genes_cutoff.png")
plot(res$P_adjust, res$ZSTAT, xlab = "p-valor ajustado", ylab = "Z stat", main = "Análisis de genes con MAGMA")
cut0 <- mean(res$ZSTAT) - 0.9*sd(res$ZSTAT)
abline(h =  cut0)
dev.off()

res_filtered <- res[res$P_adjust<0.05,] # filter by adjusted p-value, consider valid only those <0.05

cut <- mean(res_filtered$ZSTAT) - 0.9*sd(res_filtered$ZSTAT)
res_filtered <- res_filtered[res_filtered$ZSTAT > cut,] # keep only values above cutoff, which is mean - 3 times the standard deviation

write.table(x = res_filtered, quote = F, file = out_filtered_file, sep="\t", row.names = F)

print("Number of genes dropped in filtering: ", sum(res$ZSTAT < cut))
write.table(x=res_filtered$GENE, list_entrezs_file, quote=F, row.names = F, col.names = F)

sessionInfo()

