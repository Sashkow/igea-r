# source("http://bioconductor.org/biocLite.R")
# biocLite("pdInfoBuilder")
# biocLite()
# biocLite("arrayQualityMetrics")
# biocLite("oligo")
# biocLite("Biobase")
library(oligo)
library(Biobase)


library(affy)
library(ArrayExpress)
library(affycoretools)


setwd('/home/sashkoah/a/r/article-microarrays')

rawspath = 'raws/affymetrix'
prepath = 'preprocessed/affymetrix'
pdatapath = 'pdata/'
# plotsqcpath = paste(wd, 'plots/qc', sep='/')


# Load studies description
studies <- read.table("general/affymetrix_placenta_studies.tsv", header = TRUE, sep = "\t")
studies


# load IGEA phenodata
igea = read.table('igea_tsv/samples.tsv',header = TRUE, sep = '\t', fill = TRUE)

# # install cdf annotation files for all listed microarray platforms
# for (array in levels(studies$platformAbbr)){
#   install.brainarray(array)
# }

i = 7
source(paste(getwd(),'scripts/install.brainarray.R',sep='/'))
install.brainarray(studies$platformAbbr[[i]])



# i = 6 E-GEOD-36083


current_path = paste(rawspath, '/', studies$accession[[i]], sep='')
if (! dir.exists(current_path)){
  dir.create(current_path)
}



aeData = getAE(
  studies$accession[[i]],
  path = current_path,
  sourcedir=current_path,
  local = TRUE,
  type = 'raw')

z = Biobase::read.AnnotatedDataFrame(filename = paste(current_path, aeData$sdrf, sep='/'))

z@data$Source.Name = paste(z@data$Assay.Name, " 1", sep = '')  # for E-GEOD-73374


# merge ArrayExpress phenodata with IGEA phenodata
pd = merge(z@data, igea, all.x = TRUE, by.x = 'Source.Name', by.y = 'Sample.Name')


pd$Experiment
rownames(pd) = pd$Array.Data.File.y
rownames(pd)


# gseUntaredFolder = paste(gseString,  "_untared", sep = "")
celFilesPath = current_path

setwd(current_path)

filePaths = rownames(pd)


install.packages("http://mbni.org/customcdf/22.0.0/entrezg.download/pd.hugene20st.hs.entrezg_22.0.0.tar.gz", repos = NULL)
oligoData = oligo::read.celfiles(filenames = filePaths, pkgname = "pd.hugene20st.hs.entrezg")
oligoData = oligo::read.celfiles(filenames = filePaths)
oligoData@annotation


setwd('/home/sashkoah/a/r/article-microarrays')
source(paste(getwd(), 'scripts/install.brainarray.R',sep='/'))
annotationfiles = install.brainarray(studies$platformAbbr[[i]])

library(annotationfiles[1], character.only=TRUE)
library(annotationfiles[3], character.only=TRUE)


pData(oligoData) = pd
pd$Diagnosis

raw_exprs = exprs(oligoData)
rownames(raw_exprs)

View(raw_exprs)
normalizedData = oligo::rma(oligoData)

protocol.data = normalizedData@protocolData@data
write.table(protocol.data, paste(protocolpath, '/', studies$accession[[i]], "_protocol_affymetrix.tsv", sep=""), sep="\t", quote=FALSE)
studies$accession[[i]]
exprs = exprs(normalizedData)





# map probes to engrez gene identifiers in exprs
source(paste(getwd(),'scripts/install.brainarray.R',sep='/'))
annotationfiles = install.brainarray(studies$platformAbbr[[i]])
library(annotationfiles[3], character.only=TRUE)
file.db = get(annotationfiles[3])





# hugene10sthsentrezg.db
# s= select(get(annotationdb), keys(get(annotationdb))[1:10], columns(get(annotationdb)))
# keys(get(annotationdb))[1:10]

# read diffexp genes found in literature
diffexp_from_article = read.table("/home/sashkoah/a/r/article-microarrays/differential_expression_from_literature/GSE73374/pone.0141294.s006-1.csv", header = TRUE, quote = '"', sep = ",")

# biocLite("hugene20sttranscriptcluster.db")
library(hugene20sttranscriptcluster.db)
columns(hugene20sttranscriptcluster.db)
selection = select(hugene20sttranscriptcluster.db, as.character(diffexp_from_article$Probe.ID), c("ENTREZID","SYMBOL"))
write.table(selection, "clusterfuck.csv", quote=FALSE, sep=',')
getwd()

nrow(exprs)
rownames(exprs)
# select 1:1 mapping of exprs's probeid onto entrezid
annotation <- select(file.db, rownames(exprs), "ENTREZID")
keytypes(file.db)

keys(file.db, keytype = "ENTREZID")


# remove probes that map onto NA entrez id
annotation = annotation[!is.na(annotation$ENTREZID),]

# remove rows that are not annotated with entrezid
exprs = exprs[rownames(exprs) %in% annotation$PROBEID,]

# check entrezid is unique
length(unique(annotation$ENTREZID)) == length(annotation$ENTREZID)

assertthat::are_equal(nrow(annotation), nrow(exprs))
rownames(exprs) = annotation$ENTREZID

mappedpath = 'mapped/affymetrix'

write.table(exprs, paste(mappedpath, '/', studies$accession[[i]], "_mapped_affymetrix.tsv", sep=""), sep="\t", quote=FALSE)












# # map probes to engrez gene identifiers in exprs
# 
# #get genbank id to entrez id map
# file.db = get(annotationfiles[3])
# annotation = select(file.db, keys(file.db), c("ACCNUM","ENTREZID"))
# 
# View(annotation)
# 
# #get transcript cluster id to genbank id map fbrom adf.file
# path = paste(current_path, "A-GEOD-16686.adf.txt", sep="/")
# adf = read.table(path, header = TRUE, sep = '\t', skip = 15) # first 15 lines are not table
# 
# # remove rows with no genbank id
# adf = adf[!(exprs$Reporter.Database.Entry..genbank. == ""),]
# 
# merged = merge(adf, annotation,  by.x = "Reporter.Database.Entry..genbank.", by.y = "ACCNUM")
# nrow(merged)
# 
# exprs.entrez = exprs[rownames(exprs) %in% merged$Reporter.Name,]
# 
# # reorder rownames in exprs.entrez as in merged
# exprs(exprs.entrez) = exprs(exprs.entrez)[match(rownames(exprs.entrez), merged$Reporter.Name),]
# 
# match(rownames(exprs.entrez), merged$Reporter.Name)[1:10]
# exprs(exprs.entrez)[1:10,]
# merged$Reporter.Name[1:10]
# rownames(newexprs) = aandb$ENTREZID 
# 
# nrow(exprs.entrez)
# 
# columns(get(annotationdb))
# nrow(exprs)
# # select 1:1 mapping of exprs's probeid onto entrezid
# annotation <- select(get(annotationdb), rownames(exprs), "ENTREZID")
# 
# # remove probes that map onto NA entrez id
# annotation = annotation[!is.na(annotation$ENTREZID),]
# annotation[is.na(aandb),]
# 
# # remove rows that are not annotated with entrezid
# exprs = exprs[rownames(exprs) %in% annotation$PROBEID,]
# 
# # check entrezid is unique
# length(unique(annotation$ENTREZID)) == length(annotation$ENTREZID)
# 
# assertthat::are_equal(nrow(annotation), nrow(exprs))
# rownames(exprs) = annotation$ENTREZID
# 
# mappedpath = 'mapped/affymetrix'
# 
# write.table(newexprs, paste(mappedpath, '/', studies$accession[[i]], "_mapped_affymetrix.tsv", sep=""), sep="\t", quote=FALSE)
# 
# 
# 



arrayQualityMetrics::arrayQualityMetrics(expressionset = normalizedData,
                                         outdir = paste(plotsqcpath, studies$accession[[i]], sep='/'),
                                         force = TRUE,
                                         intgroup = "Diagnosis")

setwd('/home/sashkoah/a/r/article-microarrays')
write.table(exprs(normalizedData), paste(prepath, '/', studies$accession[[i]], "_preprocessed_affymetrix.tsv", sep=""), sep="\t", quote=FALSE)
write.table(pData(normalizedData), paste(prepath, '/', studies$accession[[i]], "_preprocessed_affymetrix_pdata.tsv", sep=""), sep="\t", quote=FALSE)


