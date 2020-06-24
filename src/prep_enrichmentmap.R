# download UP-TO-DATE gmt gene sets file from baderLab's repo

library(limma)
library(GSA)
library(RCurl)

working_dir <- getwd()
gmt_url = "http://download.baderlab.org/EM_Genesets/current_release/Human/Entrezgene/"
#list all the files on the server
filenames = getURL(gmt_url)
tc = textConnection(filenames)
contents = readLines(tc)
close(tc)

#get the gmt that has all the pathways and does not include terms inferred from electronic annotations(IEA)
rx = gregexpr("(?<=<a href=\")(.*.GOBP_AllPathways_no_GO_iea.*.)(.gmt)(?=\">)",
              contents, perl = TRUE)
gmt_file = unlist(regmatches(contents, rx))

dest_gmt_file <- file.path(working_dir, "data", paste("Supplementary_Table3_", gmt_file, sep="") )

download.file(
  paste(gmt_url,gmt_file,sep=""),
  destfile=dest_gmt_file
)

##################
gmt_file <- file.path("data", "Supplementary_Table3_Human_GOBP_AllPathways_no_GO_iea_June_01_2020_entrezgene.gmt") # downloaded in previous script
my_gmt_file <- file.path("data", "mi_gmt.gmt")

# get gene set file
migmt <- read.delim(gmt_file, sep = "#", header = F)

# parse it to turn it into a data frame
mi_gmt <- apply(migmt, MARGIN = 1, FUN=function(x){
  a <- strsplit(x, "\t")[[1]]
  return(data.frame(name = a[1], desc = a[2], genes = paste(a[3:length(a)], collapse = " ")))
})

mi_gmt <- do.call("rbind", mi_gmt)

replacement <- sapply(as.character(mi_gmt$name), function(x){
  
  a <- strsplit(x, "%")[[1]]
  return(a[length(a)])
  
})

mi_gmt$name <- replacement

write.table(x = mi_gmt, file = my_gmt_file, quote=F, row.names = F, col.names = T)


sessionInfo()
