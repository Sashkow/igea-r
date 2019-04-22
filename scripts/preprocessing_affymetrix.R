# Load HuGene and hgu133 BrainArray packages
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
BiocManager::install("erer")

# BiocManager::install("affycoretools")
  # 
  # BiocManager::install("sva")
  # BiocManager::install("stringr")
  # BiocManager::install("ggplot2")
  # BiocManager::install("ggfortify")
  # BiocManager::install("affy") 
  # BiocManager::install("ArrayExpress")
  # BiocManager::install("arrayQualityMetrics")
  # BiocManager::install("httr")
  # BiocManager::install("AnnotationDbi")

library(org.Hs.eg.db)
library(hugene10sthsentrezgcdf)
library(affycoretools)
library(sva)
library(stringr)
library(ggplot2)
library(ggfortify)
library(affy) 
library(ArrayExpress)
library(arrayQualityMetrics)
library(httr)
library(AnnotationDbi)

setwd('/home/sashkoah/a/r/article-microarrays')
getwd()

source(paste(getwd(),'scripts/plots.R',sep='/'))


rawspath = 'raws/affymetrix'
prepath = 'preprocessed/affymetrix'
pdatapath = 'pdata/'
plotsqcpath = paste(getwd(), 'plots/qc', sep='/')
mappedpath = 'mapped/affymetrix'
protocolpath = 'protocol/affymetrix'

# Load studies description
studies <- read.table("general/affymetrix_placenta_studies.tsv", header = TRUE, sep = "\t", fill=TRUE)
studies


# load IGEA phenodata
igea = read.table('igea_tsv/samples.tsv',header = TRUE, sep = '\t', fill = TRUE)

# install cdf annotation files for all listed microarray platforms

# for (array in levels(studies$platformAbbr)){
#   install.brainarray(array)
# }

i = 9
studies[9,]






# i = 6 E-GEOD-36083


current_path = paste(rawspath, '/', studies$accession[[i]], sep='')
if (! dir.exists(current_path)){
  dir.create(current_path)
}

current_path

aeData = getAE(
  studies$accession[[i]],
  path = current_path,
  sourcedir=current_path,
  local = TRUE,
  type = 'raw')

z <- ArrayExpress:::readPhenoData(aeData$sdrf, aeData$path)

z@data$Extract.Name

# merge ArrayExpress phenodata with IGEA phenodata
pd = merge(z@data, igea, all.x = TRUE, by.x = 'Source.Name', by.y = 'Sample.Name')
pd$Gestational.Age.Category

rownames(pd) = as.character(pd$Array.Data.File.y)
# write.table(pd,"pdata.tsv", sep="\t", quote=FALSE)
# pd = read.table("pdata.tsv", sep="\t", header=TRUE, row.names = 1)



# decidua_chorion_only = pd[which(pd$Biological.Specimen == 'Chorion' | pd$Biological.Specimen == 'Decidua'),]

# &
                                # pd$Characteristics..outcome. %in% c('TL term with labor','TNL term no labor')),]
# nrow(decidua_chorion_only)
# pd = decidua_chorion_only



nrow(pd)

affyData = ReadAffy(phenoData=pd,
                           sampleNames=pd$Sample.Name,
                           filenames=pd$Array.Data.File.y,
                           celfile.path=paste("raws/affymetrix/",
                                              studies$accession[[i]],
                                              sep="")
)



affyData@cdfName <- paste(studies$platformAbbr[[i]], 'hsentrezgcdf', sep="")
affyData@cdfName
nrow(exprs(affyData))

affyData.rma = affy::rma(affyData)
nrow(exprs(affyData.rma))

nrow(pd)
ncol(affyData.rma)

columns = c('Diagnosis', 'Gestational.Age.Category', 'accession', 'Biological.Specimen','Array.Data.File.x')
pd = pd[,columns]
write.table(pd,"pdata.tsv", sep="\t", quote=FALSE)
pd = read.table("pdata.tsv", sep="\t", header=TRUE)

# write.table(exprs(affyData.rma), paste(prepath, '/', studies$accession[[i]], "_preprocessed_affymetrix.tsv", sep=""), sep="\t", quote=FALSE)-.4


exprs = read.table(paste(prepath, '/', studies$accession[[i]], "_preprocessed_affymetrix.tsv", sep=""), sep="\t")



ncol(exprs)
nrow(pd)

# exprs@phenoData = AnnotatedDataFrame(pd)



# exclude decidua and chorion samples that are too close to maternal blood samples
exclude = c("GSM1900776_4384_40294_2503CHRN_HuGene1.0st.CEL", "GSM1900833_4384_40351_5517DSRN_HuGene1.0st.CEL",
  "GSM1900840_4384_40358_5520DSRN_HuGene1.0st.CEL", "GSM1900855_4384_40373_5524DSRN_HuGene1.0st.CEL", 
  "GSM1900870_4384_40388_5534DSRN_HuGene1.0st.CEL", "GSM1900875_4384_40393_5535DSRN_HuGene1.0st.CEL", 
  "GSM1900883_4384_40401_5537DSRN_HuGene1.0st.CEL", "GSM1900943_4384_40460_8502DSRN_HuGene1.0st.CEL"
)

pd_excluded = pd[which(!(pd$Array.Data.File.x %in% exclude)),]
nrow(pd_excluded)

pd_excluded = pd[pd$Biological.Specimen!="Adipose Tissue",]
nrow(pd_excluded)

pd_excluded = pd_excluded[which((pd_excluded$Biological.Specimen == "Chorion" |
                                  pd_excluded$Biological.Specimen == "Decidua" |
                                  pd_excluded$Biological.Specimen == "Placenta") & 
                                  pd_excluded$Diagnosis == "Healthy"),]
nrow(pd_excluded)


exprs_excluded = exprs
# exprs_excluded = exprs
ncol(exprs_excluded)
propercolnames = as.character(make.names(pd_excluded$Array.Data.File))
colnames(exprs_excluded) = as.character(colnames(exprs_excluded))
exprs_excluded = exprs_excluded[,propercolnames]
colnames(exprs_excluded) == propercolnames



pca = prcomp(t(as.matrix(exprs_excluded)))


exprs_excluded = exprs[,rownames(exprs) %in% as.character(pd_excluded$Array.Data.File)]
ncol(exprs_excluded)

pl <- pcaPlots(pca, pd_excluded, c("Biological.Specimen", "Diagnosis", "Gestational.Age.Category"), ncol=2)
pl


dir.create(paste(plotsqcpath,"manual", sep='/'), showWarnings = FALSE)
save_plot(paste(plotsqcpath, "manual", "85_excluded_chor_dec_pla_healthy.pdf", sep='/'),
          base_height=5, base_aspect_ratio = pl[[2]], pl[[1]], nrow=2)


arrayQualityMetrics::arrayQualityMetrics(
  exprs,
  outdir = paste(plotsqcpath,studies$accession[[i]],"all", sep='/'),
  force = TRUE,
  intgroup = 'Biological.Specimen'
)

write.table(exprs_excluded, paste(prepath, '/', studies$accession[[i]], "_preprocessed_no_adipose.tsv", sep=""), sep="\t", quote=FALSE)
write.table(pd_excluded , paste(prepath, '/', studies$accession[[i]], "_preprocessed_no_adipose_pdata.tsv", sep=""), sep="\t", quote=FALSE)


exprs = exprs_excluded
pd = pd_excluded

i = 5
studies[i,]$accession

# read preprocessed tsf with probesets
path = paste(prepath, "/", studies$accession[[i]], "_preprocessed_affymetrix.tsv", sep="")
exprs = read.table(path, header = TRUE, sep = '\t')
nrow(exprs)
ncol(exprs)
pd = igea[igea$Array.Data.File %in% colnames(exprs),]
nrow(pd)



eset = ExpressionSet(as.matrix(exprs))
eset@phenoData = AnnotatedDataFrame(pd)
arrayQualityMetrics::arrayQualityMetrics(
  eset,
  outdir = paste(plotsqcpath,studies$accession[[i]],"scan", sep='/'),
  force = TRUE,
  intgroup = 'Scan.Date'
)
sub_pdata$arraydatafile_exprscolumnnames



# library(biomaRt)
# mart <- useMart('ensembl', dataset="hsapiens_gene_ensembl")
# 
# filters <- listFilters(mart)
# attributes <- listAttributes(mart)



# map probes to engrez gene identifiers in exprs
source(paste(getwd(),'scripts/install.brainarray.R',sep='/'))
annotationfiles = install.brainarray(studies$platformAbbr[[i]])
library(annotationfiles[3], character.only=TRUE)
file.db = get(annotationfiles[3])

file.db

# hugene10sthsentrezg.db
# s= select(get(annotationdb), keys(get(annotationdb))[1:10], columns(get(annotationdb)))
# keys(get(annotationdb))[1:10]
rownames(exprs)

t = read.table("/home/sashkoah/a/r/article-microarrays/differential_expression_from_literature/GSE9984/NIHMS101231-supplement-Suppl_3.csv", header = TRUE, sep = "\t", quote = '"')

t = read.table("/home/sashkoah/a/r/article-microarrays/differential_expression_from_literature/GSE9984/NIHMS101231-supplement-Suppl_3.csv", header = TRUE, sep = "\t", quote = '"')

exprs = mixtures

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



## third stage



i=1
studies$accession[[i]]
exprs = read.table(paste(mappedpath, '/', studies$accession[[i]], "_mapped_affymetrix.tsv", sep=""), sep="\t", header=TRUE, row.names = 1 )

colnames(exprs)

pd = igea[make.names(igea$Array.Data.File) %in% colnames(exprs),]
nrow(pd)
pd <- data.frame(lapply(pd, as.character), stringsAsFactors=FALSE)
pd$Scan.Date





source(paste(getwd(),'scripts/plots.R',sep='/'))
pca = prcomp(t(exprs))
pl <- pcaPlots(pca, pd, c("Diagnosis", "Estimated.Fetus.Sex", "Scan.Date"), ncol=2)
pl[[2]]

dir.create(paste(plotsqcpath,studies$accession[[i]],"manual", sep='/'), showWarnings = FALSE)

save_plot(paste(plotsqcpath,studies$accession[[i]],"manual","plots.pdf", sep='/'),
          base_height=5, base_aspect_ratio = pl[[2]], pl[[1]], nrow=2)


eset = ExpressionSet(as.matrix(exprs))
eset@phenoData = AnnotatedDataFrame(pd)

arrayQualityMetrics::arrayQualityMetrics(
  eset,
  outdir = paste(plotsqcpath,studies$accession[[i]],"sex", sep='/'),
  force = TRUE,
  intgroup = 'Estimated.Fetus.Sex'
)


eset@phenoData@data$Scan.Time=lapply(eset@phenoData@data$Scan.Time, function(x) substr(x,1,2))
eset@phenoData@data$Scan.Time

# eset@phenoData@data$Scan.Date=lapply(eset@phenoData@data$Scan.Date, function(x) substr(x,1,7))
# eset@phenoData@data$Scan.Date

arrayQualityMetrics::arrayQualityMetrics(
  eset,
  outdir = paste(plotsqcpath,studies$accession[[i]],"scan", sep='/'),
  force = TRUE,
  intgroup = 'Scan.Date'
)




# eset = affyData.rma
# 
# # QC
# # Add ScanDate to pdata
# pData(eset)$ScanDate <- str_replace_all(eset@protocolData@data$ScanDate, "T", " ")
# pData(eset)$ScanDate <- sapply(strsplit(pData(eset)$ScanDate, split=' ', fixed=TRUE), function(x) (x[1]))
# pData(eset.br)$ScanDate <- pData(eset)$ScanDate
# 
# # Sometimes the scan date is in strange format, try this also
# #pData(eset.br)$ScanDate <- substr(pData(eset.br)$ScanDate, 1, 7)
# #pData(eset)$ScanDate <- substr(pData(eset)$ScanDate, 1, 7)
# 
# # Perform PCA
# pca = prcomp(t(exprs(eset)))
# pca.br = prcomp(t(exprs(eset.br)))
# 
# title <- ggdraw() + draw_label(paste("Affymetrix Probesets Definitions.", studies[i,]$ID), fontface='bold')
# 
# pl1 <- autoplot(pca, data = pData(eset), colour="Characteristics.condition.")
# pl2 <- autoplot(pca, data = pData(eset), colour="ScanDate")
# pl <- plot_grid(pl1, pl2, ncol=2, align="hv")
# pl <- plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))
# title <- ggdraw() + draw_label(paste("Brainarray Probesets Definitions.", studies[i,]$ID), fontface='bold')
# 
# pl1 <- autoplot(pca.br, data = pData(eset.br), colour="Characteristics.condition.")
# pl2 <- autoplot(pca.br, data = pData(eset.br), colour="ScanDate")
# pl3 <- plot_grid(pl1, pl2, ncol=2, align="hv")
# pl3 <- plot_grid(title, pl3, ncol=1, rel_heights=c(0.1, 1))
# pl <- plot_grid(pl, pl3, nrow=2, align="hv")
# 
# # Save plot for manual quality control
# save_plot(paste(plotsqcpath, studies[i,]$accession, "_PCA_nobatch.pdf", sep=""),
#           pl1, base_width=10, nrow=2)
# 
# 
# # Batch-effect removal. Perform only if needed.
# 
# # Eliminate samples, which are the single representation of particular batch. Reload data after that
# # with new phenoData file, cause ExpressionSet subsetting works in a weird way
# pdata = pData(eset.br)
# n_occur <- data.frame(table(pdata$ScanDate))
# uniques <- n_occur[n_occur$Freq == 1,]$Var1
# pdata <- pdata[-which(pdata$ScanDate %in% uniques),]
# pd@data <- pd@data[pd@data$SampleAccessionNumber %in% pdata$SampleAccessionNumber,]
# write.table(pd@data, paste("../pdata/pd_", studies[i,]$ID, ".tsv", sep=""), sep="\t", quote=FALSE)
# 
# # Remove batch effect
# batch = pData(eset.br)$ScanDate
# mod = model.matrix(~as.factor(CancerType), data=pData(eset.br))
# combat_edata = ComBat(dat=exprs(eset.br), batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
# exprs(eset.br) <- combat_edata
# 
# combat_edata = ComBat(dat=exprs(eset), batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
# exprs(eset) <- combat_edata
# 
# # Save affy and brain expression sets after batch-effect removal
# write.table(exprs(eset), paste("../preprocessed/", studies[i,]$ID, "_preprocessed_affymetrix.tsv", sep=""), sep="\t", quote=FALSE)
# write.table(exprs(eset.br), paste("../preprocessed/", studies[i,]$ID, "_preprocessed_brainarray.tsv", sep=""), sep="\t", quote=FALSE)
# 
# # In case not all the TNBC samples are clustered together
# # perform manual step for outliers determination within TNBC subtype based on PCA plots
# # However, pdata files already contain Outlier variable, therefor you can skip this and use
# # ready pdata files instead
# pd.tnbc <- pd[pd@data$CancerType=="TNBC"]
# pca.tnbc <- pca.br$x[which(rownames(pca.br$x) %in% pd.tnbc@data$SampleAccessionNumber),]
# outliers <- rownames(pca.tnbc[which(pca.tnbc[,1]<(37)),])
# pd@data$Outliers <- rep("NA", length(pd@data$SampleAccessionNumber))
# pd@data$Outliers[which(pd@data$SampleAccessionNumber %in% outliers)] <- "Yes"
# autoplot(pca.br, data = pData(pd), colour="Outliers")
# write.table(pd@data, paste("../pdata/pd_", studies[i,]$ID, ".tsv", sep=""), sep="\t", quote=FALSE)
# 
# ReadAffy2 <- function (..., filenames = character(0), widget = getOption("BioC")$affy$use.widgets,
#                        compress = getOption("BioC")$affy$compress.cel, celfile.path = NULL,
#                        sampleNames = NULL, phenoData = NULL, description = NULL,
#                        notes = "", rm.mask = FALSE, rm.outliers = FALSE, rm.extra = FALSE,
#                        verbose = FALSE, sd = FALSE, cdfname = NULL)
# {
#     l <- AllButCelsForReadAffy(..., filenames = filenames, widget = widget,
#                                celfile.path = celfile.path, sampleNames = sampleNames,
#                                phenoData = phenoData, description = description)
# 
#     ret <- read.affybatch(filenames = l$filenames, phenoData = l$phenoData,
#                           description = l$description, notes = notes, compress = compress,
#                           rm.mask = rm.mask, rm.outliers = rm.outliers, rm.extra = rm.extra,
#                           verbose = verbose, sd = sd, cdfname = cdfname)
#     sampleNames(ret) <- l$sampleNames
#     return(ret)
# }
