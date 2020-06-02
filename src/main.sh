bash ./src/downloads.sh
# after download of all needed files (including MAGMA software from https://ctg.cncr.nl/software/magma)
Rscript ./src/gwas_data_preprocessing.R
# MAGMA analysis
# SNPs annotation
./magma/magma_v1.07bb/magma --annotate --snp-loc data/prueba.txt --gene-loc magma/NCBI38/NCBI38.gene.loc --out results/prueba_1

# Gene analysis from SNP p-values
./magma/magma_v1.07bb/magma --bfile ./magma/g1000_eur/g1000_eur --pval ./data/prueba.txt ncol=NOBS --gene-annot results/prueba_1.genes.annot --out results/prueba_1

# Gene set analysis:
./magma/magma_v1.07bb/magma --gene-results results/prueba_1.genes.raw --set-annot [file with annotation sets, could they be individual?] self-contained  --out results/prueba_fin


# needs fixing the reference genome data
# needs an annotation set file (try to make it one for each individual gene? check if it makes sense)
