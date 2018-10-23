sexpath = 'sex/affymetrix'

# Load studies description
studies <- read.table("general/affymetrix_placenta_studies.tsv", header = TRUE, sep = "\t", fill=TRUE)
i = 9

# i = 8 pr


exprs = read.table(paste(mappedpath, '/', studies$accession[[i]], "_mapped_affymetrix.tsv", sep=""), header = TRUE, sep = '\t')


library(biomaRt)
mart <- useMart('ensembl', dataset="hsapiens_gene_ensembl")
filters <- listFilters(mart)
attributes <- listAttributes(mart)

platfrom_abbreviation = as.character(studies[i,]$martPlatformAbbr)
platfrom_abbreviation_with = paste("with_", platfrom_abbreviation, sep ="")

gene.attributes <-
  getBM(mart=mart, values=TRUE,
        filters=c(platfrom_abbreviation_with),
        attributes= c(platfrom_abbreviation, "entrezgene",
                      "chromosome_name", "start_position",
                      "end_position", "strand"))

nrow(gene.attributes)

all_genes = gene.attributes[which(!is.na(gene.attributes$entrezgene) & gene.attributes$chromosome_name=="Y"),]$entrezgene

length(all_genes)

all_genes = unique(all_genes)

length(intersect(rownames(exprs), all_genes))


library(massiR)



y = intersect(rownames(exprs), all_genes)

View(exprs[all_genes,])

genes.df = data.frame(matrix(, nrow=length(all_genes), ncol=0))
rownames(genes.df) = all_genes

massi.test.dataset = exprs
massi.test.probes = genes.df


massi.y.out <- massi_y(exprs, genes.df)
length(intersect(massi.y.out$id,y))
length(y)

massi_y_plot(massi.y.out)


massi.y.out <- massi_y(massi.test.dataset, massi.test.probes)




massi.select.out <-
  massi_select(massi.test.dataset, massi.test.probes, threshold=4)

nrow(massi.select.out)


results <- massi_cluster(massi.select.out)

results[[2]]

massi_cluster_plot(massi.select.out, results)

write.table(results[[2]], paste(sexpath, '/', studies$accession[[i]], "_sex.tsv", sep=""), sep="\t", quote=FALSE)


# 
# ### R code from vignette source 'massiR_Vignette.Rnw'
# ### Encoding: UTF-8
# 
# ###################################################
# ### code chunk number 1: load the example data
# ###################################################
# biocLite("massiR")
# library(massiR)
# 
# data(massi.test.dataset)
# massi.test.dataset
# 
# names(y.probes)
# 
# ###################################################
# ### code chunk number 2: load the test probe list
# ###################################################
# data(massi.test.probes)
# y = intersect(rownames(massi.test.dataset), rownames(massi.test.probes))
# 
# 
# ###################################################
# ### code chunk number 3: extract Y from test dataset calculate CV
# ###################################################
# massi.y.out <- massi_y(massi.test.dataset, massi.test.probes)
# intersect(massi.y.out,y)
# 
# ###################################################
# ### code chunk number 4: plot the results from massi.y
# ###################################################
# massi_y_plot(massi.y.out)
# 
# 
# ###################################################
# ### code chunk number 5: fig1too
# ###################################################
# massi.y.out <- massi_y(massi.test.dataset, massi.test.probes)
# 
# 
# ###################################################
# ### code chunk number 6: fig1
# ###################################################
# 
# barplot(height=massi.y.out[[2]], names.arg=massi.y.out[[1]], xpd=T,
#         cex.names=0.5, las=2, ylab="Probe CV (%)")
# # Get the quantile values from the massi.y output
# quantiles <- massi.y.out[[3]]
# 
# # add lines for the 0%, 25%, 50%, and 75% quartiles
# abline(h=quantiles[1:4], col=c("black", "red", "blue", "green"), lwd=2)
# legend("topleft",cex=0.7, title="Threshold (Quantile)",
#        col=c("black", "red", "blue", "green"),
#        fill= c("black", "red", "blue", "green"),
#        legend=c("1 (0%)", "2 (25%)",
#                 "3 (50%)", "4 (75%)"))
# 
# 
# ###################################################
# ### code chunk number 7: run massi.select
# ###################################################
# massi.select.out <-
#   massi_select(massi.test.dataset, massi.test.probes, threshold=4)
# 
# nrow(massi.select.out)
# 
# 
# ###################################################
# ### code chunk number 8: head massi select out
# ###################################################
# head(massi.select.out)[,1:5]
# 
# 
# ###################################################
# ### code chunk number 9: massiR_Vignette.Rnw:109-110
# ######s#############################################
# results <- massi_cluster(massi.select.out)
# 
# results[[2]]
# 
# massi_cluster_plot(massi.select.out, results)
# 
# ###################################################
# ### code chunk number 10: massiR_Vignette.Rnw:114-115
# ###################################################
# sample.results <- data.frame(results[[2]])
# 
# 
# ###################################################
# ### code chunk number 11: massiR_Vignette.Rnw:118-119
# ###################################################
# head(sample.results)
# 
# 
# ###################################################
# ### code chunk number 12: massiR_Vignette.Rnw:129-130 (eval = FALSE)
# ###################################################
# ## massi_cluster_plot(massi.select.out, results)
# 
# 
# ###################################################
# ### code chunk number 13: fig2
# ###################################################
# ord <- order(rowSums(abs(massi.select.out)),decreasing=T)
# heatmap.2(x=as.matrix(massi.select.out[ord,]), keysize=2, cexRow=0.7,
#           key=T, trace="none", dendrogram="row", col=redgreen(75), scale="row")
# 
# 
# ###################################################
# ### code chunk number 14: fig3
# ###################################################
# massi.cluster.results <- data.frame(results[[2]])
# massi.cluster.results.sort <- massi.cluster.results[order(massi.cluster.results$sex),] # sort data by sex
# probe.means <- massi.cluster.results.sort$mean_y_probes_value # samples probe mean values
# probe.sd <- massi.cluster.results.sort$y_probes_sd # sample probe sd values
# sample.names <- massi.cluster.results.sort$ID # set x-axis names
# plot.top <- ceiling(max(probe.means+probe.sd*1.1)) # set y-axis upper limit
# plot.bottom <- floor(min(probe.means-probe.sd*1.1)) # set y-axis lower limit
# sample.sex <- massi.cluster.results.sort$sex # set the factor for bar color
# # create the plot
# barCenters <- barplot(probe.means, xpd=F, names.arg=results$ID, cex.names=0.7,
#                       ylab="Chr.Y mean probe value +/- SD",
#                       xlab="",
#                       col=c("red", "green")[as.factor(sample.sex)],
#                       las=2, ylim=c(plot.bottom,plot.top))
# segments(barCenters, probe.means-probe.sd, # add the sd bars
#          barCenters, probe.means+probe.sd, lwd=0.8)
# legend("topleft", fill=c("red", "green"), title="predicted sex", ## add legend to plot
#        legend=c("female", "male"), cex=0.5, )
# 
# 
# ###################################################
# ### code chunk number 15: fig4
# ###################################################
# ## generate PC plot of clusters
# k.medoids.results <- results[[1]]
# clusplot(t(massi.select.out), k.medoids.results$clustering, color=TRUE, shade=FALSE, main="",cex.txt=0.5,
#          labels=2, lines=0)
# 
# 
# ###################################################
# ### code chunk number 16: massi.dip
# ###################################################
# dip.result <- massi_dip(massi.select.out)
# 
# 
# ###################################################
# ### code chunk number 17: fig5
# ###################################################
# dip.result <- massi_dip(massi.select.out)
# plot(dip.result[[3]])
# 
# 
# ###################################################
# ### code chunk number 18: fig6
# ###################################################
# dip.result <- massi_dip(massi.select.out)
# hist(x=dip.result[[2]], breaks=20)
# 
# 
# ###################################################
# ### code chunk number 19: bias_dip
# ###################################################
# male.ids <-
#   subset(sample.results$ID,
#          subset=sample.results$sex=="male")
# 
# female.ids <-
#   subset(sample.results$ID,
#          subset=sample.results$sex=="female")
# 
# 
# ###################################################
# ### code chunk number 20: massiR_Vignette.Rnw:246-249
# ###################################################
# bias.subset.ids <- c(female.ids[1:18], male.ids[1:2])
# 
# bias.subset <- massi.select.out[bias.subset.ids]
# 
# 
# ###################################################
# ### code chunk number 21: massiR_Vignette.Rnw:252-253
# ###################################################
# bias.dip <- massi_dip(bias.subset)
# 
# 
# ###################################################
# ### code chunk number 22: massiR_Vignette.Rnw:264-265
# ###################################################
# data(massi.eset, massi.test.probes)
# 
# 
# ###################################################
# ### code chunk number 23: massiR_Vignette.Rnw:269-274
# ###################################################
# eset.select.out <-
#   massi_select(massi.eset, massi.test.probes)
# 
# eset.results <-
#   massi_cluster(eset.select.out)
# 
# 
# ###################################################
# ### code chunk number 24: massiR_Vignette.Rnw:278-297
# ###################################################
# # Get the sex for each sample from the massi_cluster results
# eset.sample.results <-
#   data.frame(eset.results[[2]])
# 
# sexData <-
#   data.frame(eset.sample.results[c("ID", "sex")])
# 
# # Extract the order of samples in the ExpressionSet and match with results
# eset.names <-
#   colnames(exprs(massi.eset))
# 
# # match the sample order in massiR results to the same as the ExpressionSet object
# sexData <- sexData[match(eset.names, sexData$ID),]
# 
# # create an annotatedDataFrame to add to ExpressionSet
# pData <- new("AnnotatedDataFrame", data = sexData)
# 
# # add the annotatedDataFrame to the Expressionset as phenoData
# phenoData(massi.eset) <- pData
# 
# 
# ###################################################
# ### code chunk number 25: massiR_Vignette.Rnw:302-304
# ###################################################
# # check the phenodata is now within the ExpressionSet
# phenoData(massi.eset)
# 
# 
# ###################################################
# ### code chunk number 26: massiR_Vignette.Rnw:307-310
# ###################################################
# # check that all phenodata id's match expressionSet column names.
# # This must return "TRUE"
# all(massi.eset$ID == colnames(exprs(massi.eset)))
# 
# 
# ###################################################
# ### code chunk number 27: load the included probe lists
# ###################################################
# data(y.probes)
# 
# 
# ###################################################
# ### code chunk number 28: y.probes names
# ###################################################
# names(y.probes)
# 
# 
# ###################################################
# ### code chunk number 29: massiR_Vignette.Rnw:327-328
# ###################################################
# illumina.v2.probes <- data.frame(y.probes["illumina_humanwg_6_v2"])
# 
# 
# ###################################################
# ### code chunk number 30: massiR_Vignette.Rnw:340-351 (eval = FALSE)
# ###################################################
# ## library(biomaRt)
# ## mart <- useMart('ensembl', dataset="hsapiens_gene_ensembl")
# ## filters <- listFilters(mart)
# ## attributes <- listAttributes(mart)
# ##
# ## gene.attributes <-
# ##  getBM(mart=mart, values=TRUE,
# ##        filters=c("with_illumina_humanwg_6_v2"),
# ##        attributes= c("illumina_humanwg_6_v2", "entrezgene",
# ##                       "chromosome_name", "start_position",
# ##                       "end_position", "strand"))
# 
# 
# ###################################################
# ### code chunk number 31: massiR_Vignette.Rnw:354-356 (eval = FALSE)
# ###################################################
# ## unique.probe <-
# ##   subset(gene.attributes, subset=!duplicated(gene.attributes[,1]))
# 
# 
# ###################################################
# ### code chunk number 32: massiR_Vignette.Rnw:359-361 (eval = FALSE)
# ###################################################
# ## y.unique <-
# ##   subset(unique.probe, subset=unique.probe$chromosome_name == "Y")
# 
# 
# ###################################################
# ### code chunk number 33: massiR_Vignette.Rnw:365-367 (eval = FALSE)
# ###################################################
# ## illumina.v2.probes <-
# ##   data.frame(row.names=y.unique$illumina_humanwg_6_v2)
# 
