library(illuminaHumanv3.db)
library(cowplot)
library(ggfortify)
library(beadarray)
library(arrayQualityMetrics)
library(sva)
library(dplyr)
library(AnnotationDbi)
# source("plots_utils.R")
# source("degs_utils.R")

### General variables

# studies <- read.table("../pdata/studies.tsv", header = TRUE, sep = "\t")

path = "/home/sashkoah/a/r/igea-r/raws/illumina/GSE27272/part1"
path = "/home/sashkoah/a/r/igea-r/raws/illumina/GSE27272/attempt2"

# Illumina HumanRef-8 v3.0 expression beadchip

pdata = read.table(file.path(path,"pdata.tsv"), sep=",", head=TRUE, stringsAsFactors = FALSE)
setwd(path)
path
raw.data <- readIllumina(dir=path, sampleSheet=file.path(path,"pdata.csv"),
                         illuminaAnnotation="Humanv3")


# files  = list.files(path, pattern = ".idat", recursive = TRUE)
# files
# lines = character()
# for (file in files){
#   folder_code = strsplit(as.character(file),'/')[[1]][1]
#   file_code = strsplit(as.character(file),'/')[[1]][2]
#   file_code = strsplit(as.character(file_code),"[.]")[[1]][1]
#   letter_code = strsplit(as.character(file_code),"_")[[1]][2]
#   lines = c(lines, paste(file_code,file_code,folder_code,letter_code,1,sep=","))
#   # print(file_code)
# }
# 
# lines=paste(lines,sep='\n')
# write(lines,"thefile")



## Normalize and QC

## P95 and P05 scores ratio
sID <- as.factor(pdata$Sentrix_ID)
ht12metrics <- c()
for (id in levels(sID)) {
  ht12metrics <- rbind(ht12metrics, read.table(paste(path, "/", id, "/Metrics.txt", sep=""),
                                               sep = "\t", header = TRUE, as.is = TRUE))
}

pl <- illuPRatio(ht12metrics)
save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PRatio.pdf", sep=""),
          base_height=3, base_aspect_ratio = 2.5, pl, ncol=1)
save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PRatio.svg", sep=""),
          base_height=3, base_aspect_ratio = 2.5, pl, ncol=1)

## Array image plot
#imageplot(raw.data, array=16,high="white",low="darkgreen",zlim=c(4,10))
#combinedControlPlot(raw.data,array=1)

## Summarization and mormalization
myMean = function(x) mean(x,na.rm=TRUE)
mySd = function(x) sd(x,na.rm=TRUE)
greenChannel = new("illuminaChannel", greenChannelTransform, illuminaOutlierMethod, myMean, mySd,"G")

eset.sum <- beadarray::summarize(raw.data, channelList = list(greenChannel))
boxplot(exprs(eset.sum),outline=FALSE)

assayDataElementReplace(eset.sum, "exprs", as.matrix(log2(exprs(eset.sum))))
dims(eset.sum)
dim(exprs(eset.sum))


dimvalue <- dim(as.matrix(log2(exprs(eset.sum))))
dimobj <- dim(exprs(eset.sum))[seq_along(dimvalue)]
all.equal(unname(dimvalue), unname(dimobj))

eset.norm <- normaliseIllumina(eset.sum, method="neqc", transform="none")
boxplot(exprs(eset.norm),outline=FALSE)
exprs <- exprs(eset.norm)
dim(exprs)
# if (studies[i,]$ID=="london") {
#   third = which(pdata$Trimester!="Third")
# } else {
#   third = c(1:nrow(pdata))
# }

#save(list=c("exprs", "pdata"), file="out")
#load("out")
# eset = ExpressionSet(assayData=exprs[,third],
#                      phenoData = AnnotatedDataFrame(pdata[third,]))

# arrayQualityMetrics(expressionset = eset,
#                     outdir = paste("../plots/QC/", studies[i,]$ID,"/AQM_report_", studies[i,]$ID, sep=""),
#                     force = TRUE,
#                     do.logtransform = FALSE,
#                     intgroup = c("Condition"))
# 
# exprs[which(is.na(exprs))] <- min(exprs[!is.na(exprs)])
# #which(is.na(exprs), arr.ind = T)
# pca = prcomp(t(exprs[,third]))
# pl <- pcaPlots(pca, pdata[third,], c("Condition", "Sentrix_ID"))
# 
# ## Save plot for manual quality control
# save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PCA.pdf", sep=""),
#           base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])
# save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PCA.svg", sep=""),
#           base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])
# 
# ## Remove outliers
# exl = which(pdata$QC=="Outlier")
# pdata <- pdata[-exl,]
# exprs <- exprs[,-exl]
# if (studies[i,]$ID=="london") {
#   third = which(pdata$Trimester!="Third")
# } else {
#   third = c(1:nrow(pdata))
# }
# 
# eset = ExpressionSet(assayData=exprs[,third],
#                      phenoData = AnnotatedDataFrame(pdata[third,]))
# arrayQualityMetrics(expressionset = eset,
#                     outdir = paste("../plots/QC/", studies[i,]$ID,"/AQM_report_", studies[i,]$ID, "_nooutliers", sep=""),
#                     force = TRUE,
#                     do.logtransform = FALSE,
#                     intgroup = c("Condition"))
# 
# 
# ### Batch-effect removal
# 
# if (studies[i,]$ID=="london") {
#   batch = as.factor(pdata$Batch)
# } else {
#   batch = as.factor(pdata$Sentrix_ID)
# }
# mod = model.matrix(~as.factor(Condition), data=pdata)
# exprs.nobatch = ComBat(dat=exprs, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
# 
# ## Perform PCA and create plots
# pca = prcomp(t(exprs[,third]))
# pca.nobatch = prcomp(t(exprs.nobatch[,third]))
# pl <- pcaPlots(pca, pdata[third,], c("Condition", "Sentrix_ID"))
# pl.nobatch <- pcaPlots(pca.nobatch, pdata[third,], c("Condition", "Sentrix_ID"))
# 
# ## Save plot for manual quality control
# save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PCA_nooutliers.pdf", sep=""),
#           base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])
# save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PCA_nooutliers.svg", sep=""),
#           base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])
# 
# save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PCA_nooutliers_nobatch.pdf", sep=""),
#           base_height=3, base_aspect_ratio = pl.nobatch[[2]], pl.nobatch[[1]])
# save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PCA_nooutliers_nobatch.svg", sep=""),
#           base_height=3, base_aspect_ratio = pl.nobatch[[2]], pl.nobatch[[1]])
# 
# if (studies[i,]$ID=="oslo") {
#   pca.nobatch = prcomp(t(exprs.nobatch[,third]))
#   pl.nobatch <- pcaPlots(pca.nobatch, pdata[third,], c("Condition", "ConditionDetailed", "Onset", "acogPE"), ncol=2)
#   save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PCA_nooutliers_nobatch_cond.pdf", sep=""),
#             base_height=3, base_aspect_ratio = pl.nobatch[[2]], pl.nobatch[[1]], nrow=2)
#   save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PCA_nooutliers_nobatch_cond.svg", sep=""),
#             base_height=3, base_aspect_ratio = pl.nobatch[[2]], pl.nobatch[[1]], nrow=2)
# }
# 
# eset = ExpressionSet(assayData=exprs.nobatch[,third],
#                      phenoData = AnnotatedDataFrame(pdata[third,]))
# arrayQualityMetrics(expressionset = eset,
#                     outdir = paste("../plots/QC/", studies[i,]$ID,"/AQM_report_", studies[i,]$ID, "_nooutliers_nobatch", sep=""),
#                     force = TRUE,
#                     do.logtransform = FALSE,
#                     intgroup = c("Condition"))

### Save all the data (without outliers)

write.table(exprs,"/home/sashkoah/a/r/igea-r/preprocessed/illumina/GSE27272_from_raw/exprs_new.tsv", sep="\t", quote=FALSE)

ncol(exprs)




  ## Save data with unique probesets only
detach("package:ggfortify", unload=TRUE)
exprs.unique <- getUniqueProbesets(exprs.nobatch, studies[i,]$platformAbbr)
write.table(exprs.unique, paste("../exprs/", studies[i,]$ID, "_exprs.tsv", sep=""), sep="\t", quote=FALSE)




