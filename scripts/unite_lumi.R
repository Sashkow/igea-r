library(affy) 
library(ArrayExpress)
library(arrayQualityMetrics)
library(sva)


# library(stringr)
# library(ggplot2)
# library(ggfortify)

library(limma)
library(stats)

source(paste(getwd(),'scripts/plots.R',sep='/'))

setwd('/home/sashkoah/a/r/igea-r')
getwd()
source('scripts/addon.R')


rawspath = 'raws/affymetrix'
prepath = 'preprocessed/affymetrix'
pdatapath = 'pdata/'
plotsqcpath = paste(getwd(), 'plots/qc', sep='/')
mappedpath = 'mapped/affymetrix'
mergedpath = 'merged'


# Load studies description
# studies <- read.table("general/affymetrix_placenta_studies.tsv", header = TRUE, sep = "\t")
# write.table(studies, "studies.csv", sep=",", quote=FALSE)





# load IGEA phenodata
igea = read.table('igea_tsv/samples.tsv',header = TRUE, sep = '\t', fill = TRUE)



exprs_files = list.files(mappedpath, pattern = "_mapped_affymetrix")
exprs_files

# read.table(paste(mappedpath, exprs_files[1], sep = '/'), header = TRUE, sep = '\t')

paste(mappedpath, exprs_files[2], sep = '/')

mrgd = read.table(paste(mappedpath, exprs_files[1], sep = '/'), header = TRUE, sep = '\t')

# one = read.table(paste(mappedpath, exprs_files[1], sep = '/'), header = TRUE, sep = '\t')
# two = read.table(paste(mappedpath, exprs_files[2], sep = '/'), header = TRUE, sep = '\t')
# 
exprs_files[1]
# exprs_files[2]
# intersect(colnames(one), colnames(two))

current_exprs = read.table(paste(mappedpath, exprs_files[1], sep = '/'), header = TRUE, sep = '\t')

# write.table(current_exprs1, "exprs.tsv", sep="\t", quote=FALSE)
# again = read.table("exprs.tsv", header = TRUE, sep = '\t')


for (exprs_file in exprs_files[2:length(exprs_files)]){
  current_exprs = read.table(paste(mappedpath, exprs_file, sep = '/'), header = TRUE, sep = '\t')
  # print(nrow(current_exprs))
  mrgd = merge(mrgd, current_exprs, by = "row.names")
  rownames(mrgd) = mrgd$Row.names
  mrgd = mrgd[,!(colnames(mrgd) == "Row.names")]
  print(paste(exprs_file, nrow(mrgd), sep = ' '))
}


# for (col in 1:ncol(mrgd)){
#   colnames(mrgd)[col] = substr(colnames(mrgd)[col],1,9) 
# } 


# levels(igea$Sample.Name)
# split = strsplit(as.character(igea$Sample.Name[[i]])," ")
# 
# if (length(split)
# 
# for (i in 1:length(levels(igea$Sample.Name))){
#   levels(igea$Sample.Name)[i] = strsplit(as.character(igea$Sample.Name[[i]])," ")
# 
#   # colnames(mrgd)[col] = substr(colnames(mrgd)[col],1,9) 
# }

write.table(mrgd, file.path("temp", "mrgd.tsv"), sep="\t", quote=FALSE)
mrgd = read.table(file.path("temp", "mrgd.tsv"), sep="\t", header=TRUE)
pdata = igea[make.names(igea$Array.Data.File) %in% colnames(mrgd),]
write.table(pdata,file.path("temp","pdata.tsv"), sep="\t", quote=FALSE)
pdata = read.table(file.path("temp","pdata.tsv"), sep="\t", header=TRUE)

head(mrgd)


# end, next ask.R



pdata[,c("Array.Data.File", "Biological.Specimen", "accession")]

ncol(mrgd)
nrow(pdata)
nrow(pdata[which(pdata$Biological.Specimen == "Chorion"),])
nrow(pdata[which(pdata$Biological.Specimen == "Decidua"),])
nrow(pdata[which(pdata$Biological.Specimen == "Placenta"),])
nrow(pdata[which(pdata$Biological.Specimen == "Basal Plate"),])

nrow(pdata[pdata$Estimated.Fetus.Sex == "Male",])
nrow(pdata[pdata$Estimated.Fetus.Sex == "Female",])


View(pdata[,c("Sex","Fetus.Sex", "Estimated.Fetus.Sex")])



columns = c('Diagnosis', 'Gestational.Age.Category', 'accession', 'Biological.Specimen','Array.Data.File')

pdata = pdata[,columns]
pdata$accession
propercolnames = as.character(make.names(pdata$Array.Data.File))

colnames(mrgd) = as.character(colnames(mrgd))
mrgd = mrgd[,propercolnames]

# difab = setdiff(colnames(mrgd), propercolnames)
# difab
# difba = setdiff(propercolnames,colnames(mrgd))
# difba

colnames(mrgd) == propercolnames

batch = as.factor(pdata$accession)

# mod = model.matrix(~ as.factor(Gestational.Age.Category) + as.factor(Diagnosis), data=pdata)
mod = model.matrix(~ 1, data=pdata)

exprs = ComBat(dat=mrgd, batch=batch, mod=mod, par.prior=TRUE, prior.plots=TRUE)

source(paste(getwd(),'scripts/plots.R',sep='/'))
pca = prcomp(t(exprs))
pl <- pcaPlots(pca, pd, columns(pd), ncol=2)
pl

eset = ExpressionSet(as.matrix(exprs))
eset@phenoData = AnnotatedDataFrame(pd)

arrayQualityMetrics::arrayQualityMetrics(
  eset,
  outdir = paste(plotsqcpath,"merged_all", sep='/'),
  force = TRUE,
  intgroup = 'accession'
)







healthy_term = pdata[which(pdata$Diagnosis == 'Healthy' & pdata$Gestational.Age.Category == 'Term' & pdata$Biological.Specimen != 'Adipose Tissue'),]
nrow(healthy_term)

ncol(exprs.nobatch)

healthy_term$Array.Data.File

healthy_term.exprs.nobatch = exprs.nobatch[, colnames(exprs.nobatch) %in% make.names(healthy_term$Array.Data.File)]

ncol(healthy_term.exprs.nobatch)
exprs = healthy_term.exprs.nobatch
# pd = healthy_term[healthy_term$accession %in% c('E-GEOD-73685'),]
# exprs = healthy_term.exprs.nobatch[,colnames(healthy_term.exprs.nobatch) %in% make.names(pd$Array.Data.File)]
# ncol(exprs)
nrow(pd)

source(paste(getwd(),'scripts/plots.R',sep='/'))
pca = prcomp(t(exprs))
pl <- pcaPlots(pca, pd, c("Biological.Specimen", "Experiment"), ncol=2)
pl





# healthy_term.exprs.nobatch = healthy_term.exprs.nobatch[c('1442', '1444', '1081', '8788','3814', '1443', '3283', '9506', '1588'),]
nrow(healthy_term.exprs.nobatch)

pd = healthy_term
write.table(pd,"pdata.tsv", sep="\t", quote=FALSE)
pd = read.table("pdata.tsv", sep="\t", header=TRUE)
exprs = healthy_term.exprs.nobatch
nrow(pd)
ncol(exprs)
# exprs = exprsall
nrow(exprs)

#kmeans clustering


pca = prcomp(t(exprs))
sample = sample(1:nrow(exprs),1000)
class(sample)
sample


pca$center

re = kmeans(pca$x, centers = 3, iter.max = 1000, nstart = 1, trace=FALSE)

# re = kmeans(as.matrix(t(exprs)), centers = 3, iter.max = 1000, nstart = 1, trace=FALSE)
re$size

re$size
cluster = as.data.frame(re$cluster)
pd[,"cluster"] = cluster$`re$cluster`
pl <- pcaPlots(pca, pd, c("Biological.Specimen", "Experiment"), ncol=2)
pl
dir.create(paste(plotsqcpath,"manual", sep='/'), showWarnings = FALSE)
save_plot(paste(plotsqcpath, "manual", "healty_term_three_clusters.pdf", sep='/'),
          base_height=5, base_aspect_ratio = pl[[2]], pl[[1]], nrow=2)

# pdall = pd
#exprsall = exprs

nrow(pdall)

col(exprsall)

d12 = pdall[cluster!=3,]
exprs12 = exprsall[,colnames(exprsall) %in% pd12$Array.Data.File]

pd13 = pdall[cluster!=2,]
exprs13 = exprsall[,colnames(exprsall) %in% pd13$Array.Data.File]

pd = pdall
exprs = exprsall
nrow(pd)
ncol(exprs)

# differential expression analysis
fit_mod = model.matrix(~ 0 + as.factor(cluster), data=pd)
colnames(fit_mod)
colnames(fit_mod) = c('a','b','c')
# check1 = makeContrasts(b-a, levels=fit_mod)
# check2 = makeContrasts(a-b, levels=fit_mod)
# contrast.matrix = check2

contrast.matrix <- makeContrasts(b-a,c-b,c-a, levels=fit_mod)
# # View(fit_mod)
fit <- lmFit(exprs, fit_mod)  # fit each probeset to model
fit2 <- contrasts.fit(fit, contrast.matrix)
efit <- eBayes(fit2)        # empirical Bayes adjustment
# 






rame(efit$coefficients)
coefs
t = topTable(efit, number = nrow(exprs))

# tcheck1 = -t
# tcheck2 = t

# tcheck1pval = tcheck1[tcheck1$adj.P.Val<.05,]
# tcheck1pvalfc = tcheck1pval[abs(tcheck1pval$logFC)>2,]
# nrow(tcheck1pvalfc)
# 
# tcheck2pval = tcheck2[tcheck2$adj.P.Val<.05,]
# tcheck2pvalfc = tcheck2pval[abs(tcheck2pval$logFC)>2,]
# nrow(tcheck2pvalfc)
# 
# length(intersect(rownames(tcheck1pvalfc),rownames(tcheck2pvalfc))) == nrow(tcheck1pvalfc)
# length(intersect(rownames(tcheck1pvalfc),rownames(tcheck2pvalfc))) == nrow(tcheck2pvalfc)


# how b differs from a
# how c differs from b
# how c differs from a
d = t

nrow(d)
d = d[d$adj.P.Val<.05,]
nrow(d)
ba = d[abs(d$b...a)>2,]
nrow(ba)
cb = d[abs(d$c...b)>2,]
nrow(cb)

ca = d[abs(d$c...a)>2,]
nrow(ca)

# how a differs from b and c
# how b differs from a and c
# how c differs from a and b
a = intersect(rownames(ba),rownames(ca))
b = intersect(rownames(ba),rownames(cb))
c = intersect(rownames(ca),rownames(cb))
length(a)
length(b)
length(c)
write(a,"")
d = d[with(d, order(abs(b...a))),]
nrow(d)

immune_response_genes = c(10093, 10096, 1191, 1508, 1514, 1520, 1604, 2207, 23118, 29979, 3075, 3113, 3117, 3122, 3320, 3326, 5062, 51324, 5290, 5295, 5359, 55914, 5641, 5682, 5683, 5687, 5690, 5691, 5692, 5698, 5701, 5788, 7099, 71, 710, 713, 718, 7314, 7322, 7456, 929, 966)

cgenes = c(100132386, 10015, 100506144, 10066, 10093, 10096, 10140, 10151, 10159, 10203, 10204, 10209, 10370, 10376, 10409, 10437, 10456, 10457, 10487, 10541, 10549, 10550, 10577, 10578, 10581, 10627, 10628, 10632, 10640, 10730, 10914, 10923, 10944, 10945, 10955, 10957, 10959, 10964, 10972, 10983, 11010, 11067, 11161, 11261, 114882, 116254, 1191, 12, 1200, 127933, 131578, 1345, 1362, 1396, 140606, 140739, 1410, 1462, 1471, 1476, 1486, 1508, 1514, 1520, 1536, 158471, 158586, 1588, 1595, 1603, 1604, 1634, 1652, 1656, 170622, 1727, 1964, 1965, 1968, 201895, 203068, 2098, 2124, 2149, 216, 2207, 220988, 221477, 221830, 2287, 22915, 23011, 23118, 23204, 23385, 23443, 23471, 23480, 23484, 23521, 23710, 23760, 23788, 2512, 253461, 253782, 253943, 25793, 25800, 25890, 25937, 25963, 25972, 25987, 26020, 26224, 2683, 2697, 27089, 27166, 27243, 27248, 27249, 27258, 2771, 2776, 2799, 28526, 28755, 2878, 2896, 28962, 29887, 2992, 29979, 29992, 30001, 3002, 3006, 3007, 3008, 3017, 3021, 3043, 3075, 3113, 3117, 3122, 3128, 3149, 3163, 3183, 3290, 3320, 3326, 339324, 3423, 3429, 3434, 3437, 3454, 3459, 347, 3484, 3486, 3487, 3488, 3489, 3553, 3554, 3572, 360023, 3696, 3703, 3732, 374395, 3842, 3853, 387, 388650, 388962, 3920, 3956, 397, 4026, 4048, 4052, 4069, 4071, 4189, 4256, 4257, 4313, 4502, 4666, 4673, 4697, 4701, 4717, 4738, 475, 476, 4837, 4869, 4905, 4907, 492311, 4925, 4953, 4958, 498, 4982, 5047, 5062, 51009, 51014, 51094, 51100, 51108, 51142, 51167, 51246, 51249, 51258, 51324, 51371, 51604, 51614, 51635, 51639, 51643, 51669, 5167, 51699, 51714, 5175, 51762, 51768, 5204, 5229, 5270, 5290, 5295, 5355, 5358, 5359, 5412, 541471, 54205, 54210, 5440, 54499, 54543, 54545, 54664, 54680, 54741, 5476, 54809, 54816, 54918, 54947, 5501, 55017, 5504, 55076, 5516, 55167, 55173, 55186, 55196, 55207, 55233, 55249, 55252, 55272, 55297, 55303, 5544, 5547, 55505, 5552, 5571, 55742, 55752, 55788, 55819, 55860, 55914, 55959, 5617, 5641, 5660, 56650, 5682, 5683, 5687, 56889, 5690, 5691, 5692, 5698, 56990, 5701, 57124, 57222, 5788, 5813, 58472, 58494, 58505, 5868, 5928, 59350, 5947, 5962, 5997, 6035, 6038, 60493, 6122, 6135, 6138, 6156, 6169, 6185, 6188, 6208, 6217, 6218, 6224, 6229, 6248, 6275, 6281, 6282, 6302, 6355, 6362, 6366, 639, 6414, 6421, 6515, 6526, 6535, 653784, 6648, 6695, 6702, 6727, 6742, 6745, 6746, 6747, 6748, 6767, 6845, 6880, 6935, 6950, 6990, 7009, 7026, 7076, 7099, 71, 710, 713, 7132, 7179, 718, 7280, 728411, 7295, 7305, 7314, 7322, 7324, 7341, 7358, 7431, 7456, 746, 7531, 7533, 7534, 7644, 767579, 7812, 7846, 7850, 7873, 79135, 79572, 79694, 79770, 79832, 79901, 79956, 80315, 8036, 80789, 80829, 80856, 81533, 81617, 81688, 821, 826, 8349, 8358, 8364, 8366, 8394, 8404, 84230, 84790, 84925, 84987, 8519, 85236, 8554, 858, 8611, 8649, 8668, 8743, 8763, 8804, 8821, 892, 8933, 8992, 901, 90701, 91851, 92597, 929, 9295, 93487, 9349, 9375, 94239, 9442, 950, 9528, 9551, 9556, 9562, 9588, 960, 963, 966, 967, 968, 9687, 972, 975, 976, 9802, 9867, 9934, 9969)




nrow(exprsall)
ncol(exprsall)

exprs = exprsall[!(rownames(exprsall) %in% cgenes),]
nrow(exprs)
exprsall = exprs_noimune

t = t[order(abs(t$as.factor.cluster.2), decreasing = TRUE),]
t12 = t12[order(abs(t12$logFC), decreasing = TRUE),]
# t = t[order(-t$logFC),]

head(tall)
head(t12)
t12 = t
tall = t
tall
t12
# library(org.Hs.eg.db)
# sel = AnnotationDbi::select(org.Hs.eg.db, rownames(t), c("SYMBOL","GENENAME"))
# sel$SYMBOL
# 
# tabl = topTable(efit, number = 100, coef = "as.factor(Biological.Specimen)Chorion")      # table of differentially expressed probesets
# tabl = tabl[order(-tabl$logFC),]
# chorion_genes = rownames(tabl)
# 
# 
# tabl = topTable(efit, number = 100, coef = "as.factor(Biological.Specimen)Decidua")      # table of differentially expressed probesets
# tabl = tabl[order(-tabl$logFC),]
# decidua_genes=  rownames(tabl)
# 
# tabl = topTable(efit, number = 100, coef = "as.factor(Biological.Specimen)Placenta")      # table of differentially expressed probesets
# tabl = tabl[order(-tabl$logFC),]
# placenta_genes =  rownames(tabl)
# 
# 

# library(RAM)
# 
# 
# foo <- c('a','b','c','d')
# baa <- c('a','e','f','g')
# group.venn(list(foo=foo, baa=baa), label=TRUE, 
#            fill = c("orange", "blue"),
#            cat.pos = c(0, 0),
#            lab.cex=1.1)
# 
# group.venn(list(chorion_genes=chorion_genes,
#                 decidua_genes=decidua_genes,
#                 placenta_genes=placenta_genes),
#            label=TRUE, 
#            fill = c("orange", "blue"),
#            lab.cex=1.1)
# 
# 
# 
# 
# chord = intersect(chorion_genes, decidua_genes)
# 
# 
# 
# intersect(placenta_genes, chorion_genes)
# intersect(placenta_genes, decidua_genes)
# 
# genes = exclusive_or(chorion_genes, decidua_genes)
# genes2 = exclusive_or_3(chorion_genes, decidua_genes,placenta_genes)
# setdiff(genes, genes2)
# paste(genes, collapse = ", ")
# 
# # sort by decending of differential expression
# 
# 
# 
# 
# rownames(tabl)
# 
# exprs = exprs[genes,]
# 
# View(exprs['6736',])







eset = ExpressionSet(as.matrix(exprs))
eset@phenoData = AnnotatedDataFrame(pd)
pd$Continental.Population.Groups

eset.exprs.nobatch = ExpressionSet(as.matrix(healthy_term.exprs.nobatch))
eset.exprs.nobatch@phenoData = AnnotatedDataFrame(healthy_term)



arrayQualityMetrics::arrayQualityMetrics(
  eset,
  outdir = paste(plotsqcpath,"merged_healthy_term_all", sep='/'),
  force = TRUE,
  intgroup = 'cluster'
)


library(ggfortify)
autoplot(prcomp(exprs.nobatch))



