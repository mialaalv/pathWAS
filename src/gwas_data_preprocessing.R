#!/usr/local/bin/Rscript

# load data downloaded from GWAS Catalog, referring to trait EFO_0000270 (asthma)
data <- read.csv("data/efotraits_EFO_0000270-associations-2020-06-1.csv")
study <- read.csv("data/efotraits_EFO_0000270-studies-2020-06-1.csv")

# parse string to get variant name
SNP <- strsplit(data$Variant.and.risk.allele, "-") %>% sapply(function(x) x[1])

# parse SNP localization to split in chromosome and position
CHR <- strsplit(data$Location, ":") %>% sapply(function(x) x[1])
BP <- strsplit(data$Location, ":") %>% sapply(function(x) {x[2]
  # if(length(x) != 2) print(x)
  })
# there are 10 associated SNPs for which there is no mapping available to the reference genome. They will be removed later


# transform string in scientific notation to numeric value
P <- as.numeric(sub(" x 10", "e", data$P.value))

# parse the study information file to get the sample size of each study
study$sample_size <- sapply(study$Discovery.sample.description,                   # parsing the text, first remove all commas
            function(x) stringr::str_remove_all(string = x, pattern = ",")) %>%
  strsplit(" ") %>%                                                               # split by whitespaces
  sapply(function(x) sum(as.numeric(x), na.rm = T))                               # keep only elements coercible to numbers and save the sum

# store study sample size for each associated SNP 
NOBS <- study$sample_size[match(data$Study.accession,study$Study.accession)]

# checking...
all(data$Study.accession %in% study$Study.accession)
# [1] FALSE
data$Study.accession[which(!data$Study.accession %in% study$Study.accession)]
# [1] "GCST002342" "GCST002342"

# There are two SNPs whose study does not appear on the study information file (GCST002342).
# their entry is added by hand consulting the study entry on GWAS Catalog (https://www.ebi.ac.uk/gwas/publications/24486069)
df <- data.frame(SNP, CHR, BP, P, NOBS)
df[which(data$Study.accession == "GCST002342"), "NOBS"] <- 2

# removal of unmapped (no localization) associations
df <- df[which(!is.na(df$BP)),]

write.table(df, "data/EFO_0000270_pval_file.txt", quote = F, col.names =  F, row.names = F)
