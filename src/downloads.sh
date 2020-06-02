#MAGMA v1.07bb

mkdir magma

# Gene locations, build 38
wget https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI38.zip magma/NCBI38.zip

# Reference population data (1000 genomes) to correct for LD ???? ----> still need to convert to GRCH38 or find an alternative!
wget https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip magma/g1000_eur.zip


# EFO_0000270 associations GWAS Catalog
# EFO_0000270 studies GWAS Catalog
https://www.ebi.ac.uk/gwas/efotraits/EFO_0000270
# (Export data > csv)
