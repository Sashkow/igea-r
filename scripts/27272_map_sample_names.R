library(readxl)

setwd("/home/sashkoah/a/r/igea-r/preprocessed/illumina/GSE27272_from_raw")
exprs = read.table("exprs_208_24525.tsv", sep='\t', header = TRUE)
length(colnames(exprs))

mapping = read.table("/home/sashkoah/a/r/igea-r/preprocessed/illumina/GSE27272_from_raw/GSE27272_CB_PL.csv", sep = ',', header = TRUE)

mapping$Microarray



# mapping = read_excel("GSE27272_CB_PL.xlsx", sheet = 1)
f = function(X) gsub("_","", X)
colnames(exprs) = lapply(colnames(exprs), f)
length(colnames(exprs))
colnames(exprs)



make.names(mapping$Microarray)
nrow(mapping)
matching_cols = intersect(colnames(exprs), make.names(mapping$Microarray))

length(matching_cols)

mapping[which(make.names(mapping$Microarray) %in% matching_cols),]

class(colnames(exprs))
class(matching_cols)
shorter_exprs = exprs[,matching_cols]

mapping$Microarray = make.names(mapping$Microarray)
mapping = mapping[order(mapping$Microarray),]
mapping
shorter_mapping = mapping[which(mapping$Microarray %in% matching_cols),]
shorter_mapping
shorter_exprs = shorter_exprs[,shorter_mapping$Microarray]
colnames(shorter_exprs) == shorter_mapping$Microarray
colnames(shorter_exprs) = shorter_mapping$IDSample

nrow(shorter_mapping[which(grepl("PL",shorter_mapping$IDSample)==TRUE),])
pl_colnames = as.character(shorter_mapping[which(grepl("PL",shorter_mapping$IDSample)==TRUE),]$IDSample)

nrow(shorter_exprs[,pl_colnames])

write.table(shorter_exprs[,pl_colnames], "/home/sashkoah/a/r/igea-r/preprocessed/illumina/GSE27272_from_raw/GSE27272_54_24525.tsv", sep = '\t', quote = FALSE)

dim(shorter_exprs[,pl_colnames])


nrow(shorter_mapping)
ncol(shorter_exprs)
setdiff(colnames(exprs),make.names(mapping$Microarray))






# 

X4671881095A


samplenames