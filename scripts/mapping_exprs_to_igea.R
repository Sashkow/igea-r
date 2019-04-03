# fill Array Data File field in igea table based on colnames in respective exprs files

library(stringr)



setwd('/home/sashkoah/a/r/igea-r')

getwd()

rawspath = 'raws/illumina'
prepath = 'preprocessed/illumina'
mappedpath = 'mapped/illumina'
pdatapath = 'pdata/illumina'
plotsqcpath = paste(getwd(), 'plots/qc', sep='/')

# Load studies description
studies <- read.table("general/illumina_placenta_studies.tsv", header = TRUE, sep = ",")

# load IGEA phenodata
igea = read.table('igea_tsv/samples.tsv',header = TRUE, sep = '\t', fill = TRUE)

filename = "E-GEOD-60438_preprocessed_illumina_1_mapped.tsv"
pdata_filename = "E-GEOD-60438_preprocessed_illumina_1_pdata.tsv"

exprs = read.table(file.path(mappedpath,filename), header = TRUE, sep = '\t')
pdata = read.table(file.path(pdatapath,pdata_filename), header = TRUE, sep = '\t', fill = TRUE)



igea_part = igea[which(igea$Sample.Name %in% pdata$Comment..Sample_description.),]
igea_part[]


setdiff(as.character(pdata$Comment..Sample_description.),as.character(igea_part$Sample.Name))


