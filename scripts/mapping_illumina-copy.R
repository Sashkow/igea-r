# install.packages("BiocManager")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
library(stringr)
source('/home/sashkoah/a/r/igea-r/scripts/aln_illumina/degs_utils.R')


library(illuminaHumanv4.db)

setwd('/home/sashkoah/a/r/igea-r')

getwd()

rawspath = 'raws/illumina'
prepath = 'preprocessed/illumina'
mappedpath = 'mapped/illumina'
pdatapath = 'pdata/'
plotsqcpath = paste(getwd(), 'plots/qc', sep='/')

# source("https://bioconductor.org/biocLite.R")


# Load studies description
studies <- read.table("general/illumina_placenta_studies.tsv", header = TRUE, sep = ",")

# load IGEA phenodata
igea = read.table('igea_tsv/samples.tsv',header = TRUE, sep = '\t', fill = TRUE)


# 60438_preprocessed_illumina_1.tsv - illuminaHumanv4.db
# 60438_preprocessed_illumina_2.tsv - illuminaHumanv4.db

exprs_path = "/home/sashkoah/a/r/igea-r/preprocessed/illumina/E-GEOD-30186_preprocessed_illumina.tsv"
exprs_filename = basename(exprs_path)
exprs_filename = str_replace(exprs_filename,'\\.tsv','_mapped\\.tsv') 
exprs_filename

exprs = read.table(exprs_path, header = TRUE, sep = "\t", quote = '"')
# rownames(exprs) = exprs$Reporter.Identifier
# colnames(exprs)


# drops <- c("X","merged", "Reporter.Identifier")
# exprs = exprs[ , !(names(exprs) %in% drops)]
ncol(exprs)
nrow(exprs)

# exprs.save = exprs
# map probes to engrez gene identifiers in exprs

file.db = illuminaHumanv4.db


file.db


rownames(exprs)
exprs = getUniqueProbesets(exprs,"illuminaHumanv4")

# if no .db file provided
# file.txt = read.table('/home/sashkoah/a/r/igea-r/annotations/illumina/HumanWG-6_V3_0_R3_11282955_A_probe_id_entrez.txt', header = TRUE, sep = "\t")
# file.txt = read.table('/home/sashkoah/a/r/igea-r/annotations/illumina/A-MEXP-930.adf_Illumina_Human-6_v2_Expression BeadChip_probe_id_entrez.txt', header = TRUE, sep = "\t")
# 
# exprs = getUniqueProbesetsTxt(exprs)


probeset_ids = as.character(rownames(exprs))
length(probeset_ids)

exprs_filename
write.table(exprs, file.path(mappedpath, exprs_filename), sep="\t", quote=FALSE)
# write.table(illuminaHumanv4.db, file.path(mappedpath, "exprs_filename"), sep="\t", quote=FALSE)

