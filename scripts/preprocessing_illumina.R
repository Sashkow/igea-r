# Load packages

library(illuminaHumanv4.db)
library(lumi)

# BiocManager::install("lumi")


# BiocManager::install("limma")
# BiocManager::install("sva")
# BiocManager::install("stringr")
# BiocManager::install("ggplot2")
# BiocManager::install("ggfortify")
# BiocManager::install("cowplot")
# BiocManager::install("affy")
# BiocManager::install("ArrayExpress")
# BiocManager::install("illuminaHumanv4.db")
# BiocManager::install("illuminaHumanv2.db", version = "3.8")
# BiocManager::install("WGCNA")
# BiocManager::install("RAM")


library(limma)
library(sva)
library(stringr) 
library(ggplot2)
library(ggfortify)
library(cowplot)
library(affy)
library(ArrayExpress)
library(illuminaHumanv2.db)
setwd('/home/darya/Documents/diploma/igea-r')
getwd()

rawspath = 'raws/illumina'
prepath = 'preprocessed/illumina'
pdatapath = 'pdata/'
plotsqcpath = paste(getwd(), 'plots/qc', sep='/')

# source("https://bioconductor.org/biocLite.R")


# Load studies description
studies <- read.table("general/smoking_illumina_placenta_studies.tsv", header = TRUE, sep = "\t")
studies
studies$accession
# load IGEA phenodata
igea = read.table('igea_tsv/samples.tsv',header = TRUE, sep = '\t', fill = TRUE)

# # install cdf annotation files for all listed microarray platforms
# for (array in levels(studies$platformAbbr)){
#   install.brainarray(array)
# }

BiocManager::install("illuminaHumanWGDASLv3.db")
library(illuminaHumanWGDASLv3.db)
i = 1

current_path = paste(rawspath, '/', studies$accession[[i]], sep='')

if (! dir.exists(current_path)){
  dir.create(current_path, recursive = TRUE)
}

studies$accession[[i]]
aeData = getAE(
  studies$accession[[i]],
  path = current_path,
  sourcedir=current_path,
  local = FALSE,
  type = 'raw')

# aeData$sdrf
aeData$path
list.files(aeData$path)
# readPhenoData() will not work unles you rename column
# Derived Array Data File into Array Data File
# (.*?)\t "$1"\t
# 
# "\t([^"]*?)\n "\t"$1"\n

library(readr)
sdrf <- read_file(paste(current_path, aeData$sdrf, sep='/'))
sdrf = read.table(paste(current_path, aeData$sdrf, sep='/'), sep = "\t", quote = '"', header = TRUE)
sdrf$Array.Data.File
# qsdrf = str_replace(sdrf, 'Derived Array Data File', 'Array Data File')
# write_file(sdrf, paste(current_path, aeData$sdrf, sep='/'))


z <- ArrayExpress:::readPhenoData(aeData$sdrf, aeData$path)
z@varMetadata
# merge ArrayExpress phenodata with IGEA phenodata
# pd = merge(z@data, igea, all.x = TRUE, by.x = 'Source.Name', by.y = 'Sample.Name')
# pd$Array.Data.File
pd=z@data
rownames(pd) = pd$Array.Data.File


# write.table(pd, paste(current_path, 'merged_phenodata', sep = '/'))

# pd[is.na(pd$Experiment),]$Array.Data.File.x

nrow(pd)

# 
# # batches
# batch1 = pd[pd$Array.Design.REF == 'A-GEOD-10558',]
# batch2 = pd[pd$Array.Design.REF == 'A-MEXP-1173',]
# pd = batch1
# nrow(batch2)
# 
# map probes to engrez gene identifiers in exprs


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


# source(paste(getwd(),'scripts/install.brainarray.R',sep='/'))
# annotationfiles = install.brainarray(studies$platformAbbr[[i]])
# library(annotationfiles[3], character.only=TRUE)
# file.db = get(annotationfiles[3])

file.db

# hugene10sthsentrezg.db
# s= select(get(annotationdb), keys(get(annotationdb))[1:10], columns(get(annotationdb)))
# keys(get(annotationdb))[1:10]
rownames(exprs)

t = read.table("/home/sashkoah/a/r/article-microarrays/differential_expression_from_literature/GSE9984/NIHMS101231-supplement-Suppl_3.csv", header = TRUE, sep = "\t", quote = '"')

t = read.table("/home/sashkoah/a/r/article-microarrays/differential_expression_from_literature/GSE9984/NIHMS101231-supplement-Suppl_3.csv", header = TRUE, sep = "\t", quote = '"')


probeset_ids = as.character(rownames(exprs))
length(probeset_ids)

# select 1:1 mapping of exprs's probeid onto entrezid
# annotation <- AnnotationDbi::select(file.db, rownames(exprs), "ENTREZID")

annotation <- AnnotationDbi::select(file.db, probeset_ids, "ENTREZID")

nrow(annotation)
nrow(exprs)

# remove probes that map onto NA entrez id
annotation = annotation[!is.na(annotation$ENTREZID),]

# remove rows that are not annotated with entrezid
exprs = exprs[rownames(exprs) %in% annotation$PROBEID,]

# check entrezid is unique
length(unique(annotation$ENTREZID)) == length(annotation$ENTREZID)
FALSE %in% (rownames(exprs) == annotation$PROBEID)

assertthat::are_equal(nrow(annotation), nrow(exprs))

rownames(exprs) = annotation$ENTREZID



write.table(exprs, paste(mappedpath, '/', studies$accession[[i]], "_mapped_affymetrix_no_adipose.tsv", sep=""), sep="\t", quote=FALSE)

# read.table("general/affymetrix_placenta_studies.tsv", header = TRUE, sep = "\t", fill=TRUE)




















# source(paste(getwd(),'/scripts/merge_lumi_tables.R',sep=''))
# merged_samples = mergeLumiTables(current_path)
# merged_samples = mergeSomeLumiTables(paste(current_path, batch1, sep = '/'))
# # View(merged_samples)
# write.table(merged_samples, paste(prepath, '/', studies$accession[[i]], "_preprocessed_illumina.tsv", sep=""), sep="\t", quote=FALSE)
# 
# probesetsID_EntrezID<-select(get(paste(studies[i,]$platformAbbr, ".db", sep="")), probesetsID, "ENTREZID")


# assayData = as.matrix(merged_samples)
# View(assayData) 
# GSE60438_non_normalized_WG6v3.txt
# 
# GSE44711_non-normalized_data.txt
# GSE30186_non_normalized.txt
x.lumi <- lumiR(file.path(current_path, "E-MTAB-6418.sdrf.txt"),
                sep = "\t", QC=FALSE, columnNameGrepPattern = list(exprs='AVG_Signal', se.exprs='BEAD_STD', detection='Detection', beadNum='Avg_NBEADS'))

ncol(x.lumi) 
# x.lumi <- x.lumi[, -which(grepl("replicate", colnames(x.lumi)))]
# Normaliz(e
lumi.T <- lumiT(x.lumi, method = 'log2')
View(lumi.T)
lumi.N <- lumiN(lumi.T, method = "quantile")

# align pheno and expression data
rownames(pd) = pd$Comment..Sample_title.
colnames(lumi.N) <- pd$Source.Name
pData(lumi.N) = pd




write.table(exprs(lumi.N), paste(prepath, '/', studies$accession[[i]], "_preprocessed_illumina.tsv", sep=""), sep="\t", quote=FALSE)
write.table(pData(lumi.N), paste(prepath, '/', studies$accession[[i]], "_preprocessed_illumina_pdata.tsv", sep=""), sep="\t", quote=FALSE)

class(lumi.N)

View(exprs(lumi.N))

arrayQualityMetrics::arrayQualityMetrics(
  lumi.N,
  outdir = paste(plotsqcpath,studies$accession[[i]],"1", sep='/'),
  force = TRUE,
  intgroup = "Diagnosis"
)


# 
# 
# files = list.files(current_path, pattern = "*sample_table.txt")
# 
# x = read.ilmn(paste(current_path, files[1], sep='/'))
# 
# arrayQualityMetrics::arrayQualityMetrics(eset)
# 
# x.lumi <- lumiQ(eset)
# 
# 
# 
# x.lumi <- x.lumi[, -which(grepl("replicate", colnames(x.lumi)))]
# # Normalize
# lumi.T <- lumiT(x.lumi, method = 'log2')
# lumi.N<-lumiN(lumi.T, method = "quantile")
# colnames(lumi.N) <- pd@data$SampleAccessionNumber
# # Save results into file
# write.table(exprs(lumi.N), paste("../preprocessed/", studies[i,]$ID, "_preprocessed_illumina.tsv", sep=""),
#             sep="\t", quote=FALSE)
# 
# # Add phenoData to ExpressionSet object
# lumi.N@phenoData = pd
# # Perform PCA and plot
# pca = prcomp(t(exprs(lumi.N)))
# title <- ggdraw() + draw_label(studies[i,]$ID, fontface='bold')  
# pl <- autoplot(pca, data = pData(lumi.N), colour="CancerType")
# pl <- plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))
# 
# # Save plot for manual quality control
# save_plot(paste("../plots/qc/", studies[i,]$ID, "_PCA.pdf", sep=""),
#           pl)
# 
# #probesetsID <- rownames(exprs(lumi.N))
# #probesetsID_EntrezID<-select(get(paste(studies[i,]$platformAbbr, ".db", sep="")), probesetsID, "ENTREZID")
# 
