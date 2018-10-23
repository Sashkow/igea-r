# Load packages
library(illuminaHumanv4.db)
library(lumi)
library(limma)
library(sva)
library(stringr) 
library(ggplot2)
library(ggfortify)
library(cowplot)
library(affy)
library(ArrayExpress)
setwd('/home/sashkoah/a/r/article-microarrays')
getwd()

rawspath = 'raws/illumina'
prepath = 'preprocessed/illumina'
pdatapath = 'pdata/'
plotsqcpath = paste(getwd(), 'plots/qc', sep='/')

# source("https://bioconductor.org/biocLite.R")


# Load studies description
studies <- read.table("general/illumina_placenta_studies.tsv", header = TRUE, sep = "\t")

# load IGEA phenodata
igea = read.table('igea_tsv/samples.tsv',header = TRUE, sep = '\t', fill = TRUE)

# install cdf annotation files for all listed microarray platforms
for (array in levels(studies$platformAbbr)){
  install.brainarray(array)
}

biocLite("illuminaHumanWGDASLv3.db")

# i = 1 two
# i = 2 not mrna
i = 1

current_path = paste(rawspath, '/', studies$accession[[i]], sep='')
if (! dir.exists(current_path)){
  dir.create(current_path, recursive = TRUE)
}


aeData = getAE(
  studies$accession[[i]],
  path = current_path,
  sourcedir=current_path,
  local = TRUE)
  # type = 'processed')

aeData$sdrf
aeData$path

# readPhenoData() will not work unles you rename column
# Derived Array Data File into Array Data File
(.*?)\t "$1"\t

"\t([^"]*?)\n "\t"$1"\n

library(readr)
sdrf <- read_file(paste(current_path, aeData$sdrf, sep='/'))
sdrf = read.table(paste(current_path, aeData$sdrf, sep='/'), sep = "\t", quote = '"', header = TRUE)
sdrf$Array.Data.File
qsdrf = str_replace(sdrf, 'Derived Array Data File', 'Array Data File')
write_file(sdrf, paste(current_path, aeData$sdrf, sep='/'))


z <- ArrayExpress:::readPhenoData(aeData$sdrf, aeData$path)

# merge ArrayExpress phenodata with IGEA phenodata
pd = merge(z@data, igea, all.x = TRUE, by.x = 'Source.Name', by.y = 'Sample.Name')
pd$Array.Data.File
rownames(pd) = pd$Array.Data.File.x
pd$Experiment

# write.table(pd, paste(current_path, 'merged_phenodata', sep = '/'))

pd[is.na(pd$Experiment),]$Array.Data.File.x

nrow(pd)


# batches
batch1 = pd[pd$Array.Design.REF == 'A-GEOD-10558',]
batch2 = pd[pd$Array.Design.REF == 'A-MEXP-1173',]
pd = batch1
nrow(batch2)


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
x.lumi <- lumiR(paste(current_path, "GSE60438_non_normalized_HT12v4.txt", sep="/"),
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
