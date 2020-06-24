#!/bin/bash

mkdir -p logs # create folder to store logs from executions
mkdir -p results # create folder to store result files

### Set up
echo MAGMA installation and download of resources
# Install MAGMA software from https://ctg.cncr.nl/software/magma
bash ./src/magma_setup.sh &> logs/magma_setup.log
echo DONE

echo Installation of R dependencies
# Install R dependencies
Rscript ./src/R_dependencies.R &> logs/R_dependencies.log
echo DONE

echo Fetching GWAS Catalog data
# fetch data from GWAS Catalog usin gwasrapidd bioconductor package. Lifts GRCh38 SNP locations to GRCh37 build.
Rscript ./src/gwas_data_fetching_with_R.R &> logs/gwas_data_fetching_with_R.log
echo DONE


### MAGMA analysis
echo MAGMA Analysis
echo SNPs annotation
# SNPs annotation
./magma/magma --annotate --snp-loc data/GWAS_data/EFO_0000270/EFO_0000270_GWAS_summary_grch37_forMAGMA.tsv --gene-loc magma/res/NCBI37.3/NCBI37.3.gene.loc --out results/EFO_0000270
# make sure the SNP file and gene location file are in the same genome build
echo DONE

echo Gene analysis
# Gene analysis from SNP p-values
./magma/magma --bfile ./magma/res/g1000_eur/g1000_eur --pval data/GWAS_data/EFO_0000270/EFO_0000270_GWAS_summary_grch37_forMAGMA.tsv ncol=NOBS --gene-annot results/EFO_0000270.genes.annot --out results/EFO_0000270
# make sure the reference panel is in the same genome build as the SNP pval data
# NOBS is the initial sample size of the study the data were retrieved from
# *.genes.annot is the output from the annotation step.
echo DONE

echo P-val adjusting
# Adjusting pvalues for MAGMA output and filtering out the worst associations:
Rscript ./src/adjust_output_pvals.R &> logs/adjust_output_pvals.log
echo DONE

echo EnrichmentMap preparation
Rscript ./src/prep_enrichmentmap.R &> logs/prep_enrichmentmap.log
echo DONE


