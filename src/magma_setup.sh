# MAGMA v1.07bb installation and fetching of necessary resources

mkdir -p magma/res

# Install MAGMA
wget --no-check-certificate https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.07bb.zip -P magma
echo YA
unzip magma/magma_v1.07bb.zip -d magma

# Gene locations, build GRCh37
wget https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip -P magma/res
unzip magma/res/NCBI37.3.zip -d magma/res/NCBI37.3

# Reference population data (1000 genomes, GRCh37 genome build) to correct for LD
wget https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip -P magma/res
unzip magma/res/g1000_eur.zip -d magma/res/g1000_eur



