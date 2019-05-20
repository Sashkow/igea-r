# BiocManager::install("Rtsne")
library(Rtsne)

library(affy) 
library(ArrayExpress)
library(arrayQualityMetrics)
library(sva)
library(limma)
library(stringr)
library(RColorBrewer)
library(factoextra)

setwd('/home/sashkoah/a/r/igea-r')
source(paste(getwd(),'scripts/plots.R',sep='/'))
plotsqcpath = paste(getwd(), 'plots/qc/illumina', sep='/')


# read merged expression data and metadata created with unite.R
mrgd = NULL
pdata = NULL
mrgd = read.table(file.path("temp", "mrgd.tsv"), sep="\t", header=TRUE)
pdata = read.table(file.path("temp","pdata.tsv"), sep="\t", header=TRUE)


ncol(mrgd)
nrow(pdata)


# pdata = read.table('home/sashkoah/a/r/igea-r/pure_mixed/igea-cibersort-python/73685_60438/combined_re.csv', sep="\t", header=TRUE)

# pdata = read.table(file.path("temp","pdata.csv"), sep=",", header=TRUE)
# pdata$Array.Data.File = pdata$Expression.Data.ID
mrgd
# columns = c('Diagnosis',
#             'Gestational.Age.Category',
#             'accession',
#             'Biological.Specimen',
#             'Array.Data.File',
#             'Gestational.Age', 
#             "exprs_column_names",
#             "platform_batch",
#             "arraydatafile_exprscolumnnames", # column with expression data sample names
#             "Fetus.Sex",
#             "Gravidity",
#             "Parity",
#             "medical_sample_name",
#             "Caesarean.Section")
# 
# # align metadata in pdata and expression data in mrgd so that rows in pdata match cols in mrgd
# pdata = pdata[,columns]

pdata$Biological.Specimen
# pdata$Array.Data.File = pdata$arraydatafile_exprscolumnnames
# pdata$Array.Data.File = pdata$Input.Sample
propercolnames = as.character(make.names(pdata$Array.Data.File))

rownames(pdata) = propercolnames
propercolnames
colnames(mrgd) = as.character(colnames(mrgd))
colnames(mrgd)
propercolnames
mrgd = mrgd[,propercolnames]
# check correctness
colnames(mrgd) == propercolnames
colnames(mrgd) == rownames(pdata)


# write columnames of expression data in mgrd to Array.Data.File in pdata for backwards compatability
pdata$Array.Data.File = colnames(mrgd)
pdata$Biological.Specimen


# remove technical batch effect caused by different datasets 
batch = as.factor(pdata$accession)
batch = as.factor(pdata$secondaryaccession)
# batch = as.factor(pdata$Batch)
pdata$smoking.status.ch1
# mod = model.matrix(~ as.factor(Gestational.Age.Category) + as.factor(Diagnosis), data=pdata)
# mod = model.matrix(~as.factor(Gestational.Age.Category), data=pdata)
mod = model.matrix(~as.factor(smoking.status.ch1), data=pdata)
colnames(mrgd) == rownames(pdata)
mod
# mod = model.matrix(~1)

# nrow(mrgd)
# mrgd1 <- mrgd[rowSums(mrgd) > 1, ]
# nrow(mrgd1)
exprs = mrgd
exprs = ComBat(dat=as.matrix(mrgd), batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)

# set trim: First Trimester, Second Trimester, Third Trimester
pdata$trim = with(pdata, ifelse(Gestational.Age.Category == "Term" | Gestational.Age.Category == "Late Preterm" | Gestational.Age.Category == "Early Preterm", "Third Trimester", as.character(Gestational.Age.Category)))
sub_pdata = pdata
nrow(pdata)
ncol(exprs)
sub_exprs = exprs




# sub_pdata$is_outlier = FALSE
# sub_pdata[which(sub_pdata$arraydatafile_exprscolumnnames %in% outlier_elements),]$is_outlier = TRUE

sub_pdata = pdata[which(pdata$Diagnosis=="Healthy"),]
sub_pdata = pdata[pdata$secondaryaccession=="GSE27272",]

sub_pdata = sub_pdata[which(sub_pdata$accession=='E-GEOD-73685' | sub_pdata$Decidua>.75),]
# sub_pdata = sub_pdata[which(sub_pdata$accession == "E-GEOD-73685"),]
sub_pdata$accession
                              
# sub_pdata = sub_pdata[which(sub_pdata$`tissue:ch1`=="term placenta"),]
# sub_pdata = pdata[which(!(pdata$Array.Data.File %in% c("B1_Cell_2", "B1_Cell_5", "B2_Cell_29"))),]
nrow(sub_pdata)
# sub_pdata = pdata[which(pdata$Diagnosis=="Healthy" & 

#                           pdata$trim=='Third Trimester' & 
#                           pdata$Biological.Specimen!="Umbilical Cord Blood" &
#                           pdata$Biological.Specimen!="Maternal Blood" &
#                           pdata$Biological.Specimen!="Amnion" &
#                           pdata$Biological.Specimen!="Lower Segment" &
#                           pdata$Biological.Specimen!="Uterus Fundus"),]
# sub_pdata = pdata[which(pdata$Diagnosis=="Pre-Eclampsia"),]

#
exprs = mrgd
source('scripts/sub_pd_to_sub_exprs.R')

ncol(sub_exprs) 
nrow(sub_pdata)

# how many samples in each dataset
as.data.frame(table(sub_pdata$accession))[all_studies$accession,]
as.data.frame(table(sub_pdata$secondaryaccession))



# # sort of takes sub_pdata and sub_exprs as inputs
# source('scripts/diff_exp.R')




# sub_pdata = pdata[which((pdata$Gestational.Age.Category == 'Term' | pdata$Gestational.Age.Category == 'Early Preterm' | pdata$Gestational.Age.Category == 'Late Preterm') & (pdata$Biological.Specimen=='Placenta' | pdata$Biological.Specimen=='Chorion')), ]

# 
# sub_pdata = pdata[which((pdata$Gestational.Age.Category == 'First Trimester' | pdata$Gestational.Age.Category == 'Second Trimester')),]
# 
# sub_pdata = pdata[which(pdata$Diagnosis == "Healthy" & pdata$Biological.Specimen == "Placenta" & (pdata$Gestational.Age.Category == 'Term' | pdata$Gestational.Age.Category == 'Early Preterm' | pdata$Gestational.Age.Category == 'Late Preterm')),]

# sub_pdata  = pdata[pdata$accession=="E-GEOD-24129",]
nrow(sub_pdata)
sub_pdata$accession
nrow(sub_pdata)
nrow(sub_pdata[sub_pdata$Diagnosis=='Healthy',])
nrow(sub_pdata[sub_pdata$Diagnosis=='Pre-Eclampsia',])
nrow(sub_pdata[sub_pdata$Diagnosis=='Severe Pre-Eclampsia',])
nrow(sub_pdata[sub_pdata$Diagnosis=='Fetal Growth Retardation',])

# write.table(sub_pdata,"sub_pdata.tsv", sep="\t", quote=FALSE)
# sub_pdata = read.table("sub_pdata.tsv", sep="\t", header=TRUE)
# 
# sub_pdata$Array.Data.File
# propercolnames = as.character(make.names(sub_pdata$Array.Data.File))
# 
# # sub_exprs = exprs
# colnames(sub_exprs) = as.character(colnames(sub_exprs))
# sub_exprs = sub_exprs[,propercolnames]
# colnames(sub_exprs) == propercolnames
# 
getwd()
write.table(sub_exprs,"sub_exprs.tsv", sep="\t", quote=FALSE)
write.table(sub_pdata,"sub_pdata.tsv", sep="\t", quote=FALSE)

sub_exprs = read.table("sub_exprs.tsv", sep="\t", header=TRUE)
sub_pdata = read.table("sub_pdata.tsv", sep="\t", header=TRUE)

eset = ExpressionSet(as.matrix(sub_exprs))
eset@phenoData = AnnotatedDataFrame(sub_pdata)

source(paste(getwd(),'scripts/plots.R',sep='/'))

# colnames(sub_exprs) = lookUp(as.character(colnames(sub_exprs)), 'org.Hs.eg', 'SYMBOL') if 
# sub_exprs = sub_exprs[which(colnames(sub_exprs)!=NA),]
# length(colnames(sub_exprs)==NA)


minrows
#_____________
# length(levels(sub_pdata$Biological.Specimen))
# re = kmeans(pca$x[,1:100],
#             centers = length(levels(sub_pdata$Biological.Specimen))-5,
#             iter.max = 1000,
#             nstart = 50,
#             # trace=FALSE
#             )
# re$iter
# re$ifault
# # re = kmeans(as.matrix(t(exprs)), centers = 3, iter.max = 1000, nstart = 1, trace=FALSE)
# re$size
# re
# 
# re$cluster
# cluster = as.data.frame(re$cluster)
# cluster
# sub_pdata$Array.Data.File = as.character(make.names(sub_pdata$Array.Data.File))
# sub_pdata$Array.Data.File == rownames(cluster)
# 
# sub_pdata[,"cluster"] = cluster$`re$cluster`
# 
# nrow(pca$x)
# nrow(sub_pdata)
#_______________

# f <- function(sub_pdata) {
#   print(sub_pdata)
#   if (sub_pdata[3]=='E-GEOD-60438'){
#     return('Outlier')
#   } else {
#     return(sub_pdata[4])
#   }
# }
# 
# sub_pdata$temp = apply(sub_pdata, 1, f)
# sub_pdata$temp


# pl <- pcaPlots(pca, axis(1,2),
#                sub_pdata,
#                c('accession','trim', 'Biological.Specimen', 'Fetus.Sex', 'Gravidity','Parity',"Caesarean.Section"),
#                title = 'current',
#                ncol=2)
# pl
 # pl = NULL

# pl <- fviz_pca(pca,
#                       geom=c("point", "text"),
#                       # label=c("ind", "ind.sup", "quali", "var", "quanti.sup"),
#                       label=c("var"),
#                       habillage=sub_pdata$Biological.Specimen,
#                       col.var = "contrib",
#                       # fill.ind = factor(pd$Donor),
#                       #palette = indCol,
#                       col.ind = "black",
#                       pointshape=21, pointsize = 3,
#                       palette="Dark2",
#                       repel = TRUE,     # Avoid text overlapping
#                       select.var = list(contrib=20),
#                       invisible = "quali"
# )

# sub_pdata$Biological.Specimen.Outliers = as.character(sub_pdata$Biological.Specimen)
# pca_df = as.data.frame(pca$x)
# # blood outliers
# blood_outliers = rownames(pca_df[pca_df$PC1 < -50,])[1:8]
# blood_outliers
# sub_pdata[,]
# s
# ncol(pca_df)
sub

pca = prcomp(t(na.omit(mrgd)))
pca = prcomp(t(na.omit(exprs)))

t(na.omit(sub_exprs))

pca = prcomp(t(sub_exprs))

pca = prcomp(t(sub_exprs[X$SD>800,]))
pca = prcomp(t(log2(sub_exprs + 1)))

pca_df = as.data.frame(pca$x)
pl = plot(pca_df$PC1,pca_df$PC2)
pl
sub_pdata$evt5

pl = fviz_pca_ind(pca, axes = c(1,2),
                      label=c("none"),
                      habillage=sub_pdata$smoking.status.ch1,
                      # repel = TRUE,
                      title = "GSE18044, GSE27272 and GSE7434",
                      palette = "Set4"
                      # addEllipses = TRUE,
                      # col.ind = sub_pdata$evt5,
                      # gradient.cols = rainbow(3)
                  
          
) 
pl
nrow(pdata[which(pdata$accession=="E-GEOD-60438" & pdata$Diagnosis=="Healthy"),])
nrow(sub_pdata[sub_pdata$accession=="E-GEOD-60438",])


these = sub_pdata[sub_pdata$accession=='E-GEOD-60438',]
these
View(sub_pdata)
# pdata[which(pdata$accession=='E-GEOD-60438' & pdata$Biological.Specimen,]
# + geom_point(aes(colour=as.character(sub_pdata$Biological.Specimen)), size=2)
pl



summary(pca)

library(rgl)

levels(sub_pdata$Biological.Specimen)
length(levels(sub_pdata$Gestational.Age.Category))
sub_pdata$temp = as.factor(sub_pdata$temp)
sub_pdata$Cluster = as.factor(sub_pdata$Cluster)
levels(sub_pdata$Cluster)
plot3d(pca$x[,1:3], col=rainbow(length(levels(sub_pdata$Cluster)))[factor(as.integer(sub_pdata$Cluster))], size=10)
legend3d("topright", legend = levels(sub_pdata$Cluster), pch = 16, col = rainbow(length(levels(sub_pdata$Custer))), cex=1, inset=c(0.02))

text3d(pc$scores[,1:3],texts=rownames(iris))
text3d(pc$loadings[,1:3], texts=rownames(pc$loadings), col="red")
coords <- NULL
for (i in 1:nrow(pc$loadings)) {
  coords <- rbind(coords, rbind(c(0,0,0),pc$loadings[i,1:3]))
}
lines3d(coords, col="red", lwd=4)





plot3d(pca$x[,3:5], col=rainbow(2)[factor(as.integer(sub_pdata$accession))], size=10)
legend3d("topright", legend = levels(sub_pdata$accession), pch = 16, col = rainbow(2), cex=1, inset=c(0.02))

pl = fviz_contrib(pca, choice = "var", axes = 1, top = 20)
length(pl$data[order(-pl$data$contrib),][1:50,]$name)
genes = pl$data[order(-pl$data$contrib),][1:50,]$name




rotation_df = as.data.frame(pca$rotation)
"NA.917" %in% rotation_df
genes = rownames(rotation_df[which((rownames(rotation_df) %in% genes) & rotation_df$PC1>0),])

length(genes)
symbols = lookUp(as.character(genes), 'org.Hs.eg', 'SYMBOL')

# genes = lookUp(c("CCR7", "TMEM71", "PHOSPHO1", "DEFA4", "EPB42", "ALAS2", "SULT1B1", "TRAC", "GYPA", "HBM", "HBZ", "IL7R", "CXCR1", "CXCR2", "AQP9", "MIR223", "IFIT1B", "NFE2", "TRAT1", "LEF1", "AHSP", "HEMGN", "S100A12", "SELL", "LOC644462", "SLC4A1", "SLC14A1", "CA1", "ADGRE3", "VNN2", "MGAM", "SELENBP1", "CD3G", "MS4A1"), 'org.Hs.eg', 'ENTREZ')

print(symbols, row.names = FALSE)
getwd()
write.table(as.data.frame(genes), file="genes1", row.names = FALSE)




iris_matrix = t(as.matrix(sub_exprs[X$SD>800,]))

iris_matrix = t(as.matrix(log(sub_exprs+1)))

iris_matrix = t(as.matrix(na.omit(sub_exprs)))

pca = prcomp(t(na.omit(mrgd)))
# Set a seed if you want reproducible results
set.seed(42)
tsne_out <- Rtsne(iris_matrix,perplexity=40,theta=0.0, normalize = FALSE) # Run TSNE

# Show the objects in the 2D tsne representation
plot(tsne_out$Y, col=sub_pdata$smoking.status.ch1)

set.seed(007)
X <- sub_exprs
X = transform(X, SD=apply(X,1, sd, na.rm = TRUE))
nrow(sub_exprs[X$SD>500,])

# Placenta vs uterus
# 1588
# 5670
# 5671
# 5672
# 5673
# 
# blood vs tissues
# 51327
# 6521
# 
# chorion(?)
# 4060
# 1410
# 
# decidua, uterus
# 8404
# 


pl
class(list(contrib=20))

pl$data$contrib[pl$data$contrib>1]



brewer.pal(8, "Set2")
ggloadings <- function(pca, num_pca=1, num_loading=30, xlab="genes") {
  loadings <- pca$rotation[,num_pca]
  loadings <- abs(loadings)
  loadings <- sort(loadings, decreasing = TRUE)
  label_names <- names(loadings)
  loadings <- loadings[1:num_loading]
  label_names <- label_names[1:num_loading]
  dat <- data.frame(pc = label_names, loadings = loadings)
  dat$pc <- factor(dat$pc, levels = unique(dat$pc))
  
  p <- ggplot(dat, aes(x = pc, y = loadings))
  p <- p + geom_bar(stat = "identity")
  p <- p + xlab(xlab) + ylab("contribution scores")
  return(p)
}

pl <- ggloadings(pca, num_pca=1)
pl

# outlier_rows = rownames(pl$data[which(pl$data$PC2>.2 | pl$data$PC1>.2),])
# outlier_elements = pl$data[outlier_rows,]$arraydatafile_exprscolumnnames
# outlier_elements

dir.create(paste(plotsqcpath,"manual", sep='/'), showWarnings = FALSE)

save_plot(paste(plotsqcpath, "manual" , "pca_2.png", sep='/'),
          base_height=5, base_aspect_ratio = pl[[2]], pl[[1]], nrow=2,limitsize = FALSE)

save_plot(paste(plotsqcpath, "manual" , "pca_arrows", sep='/'), pl)


library(org.Hs.eg.db)
library(annotate)
re = lookUp(c('1588', '5670', '5671', '5672', '5673', '51327', '6521', '4060', '1410', '8404'), 'org.Hs.eg', 'SYMBOL') 

re = lookUp(c("162466", "212", "2993", "3577", "3579", "407008", "439996", "4778", "51327", "6283", "6402", "6521", "6563", "759", "84658", "8875"),'org.Hs.eg', 'SYMBOL')
as.data.frame(re)
re = lookUp(c('212','51327','4778', '6521', '759', '8875','6283','6402'), 'org.Hs.eg', 'SYMBOL') 

re = lookUp(c( '5670', '9502'), 'org.Hs.eg', 'SYMBOL') 
re = lookUp(c( '5670', '9502'), 'org.Hs.eg', 'SYMBOL') 
re = lookUp(c( '388136', '2335'), 'org.Hs.eg', 'SYMBOL') 
re = lookUp(c( '1410', '1634','10979','4060',), 'org.Hs.eg', 'SYMBOL') 
re = lookUp(c( '3488','4256','22943','8404'), 'org.Hs.eg', 'SYMBOL') 
print(re, row.names = FALSE)

getwd()
as.data.frame()
write.table(as.data.frame(re), file="genes1", row.names = FALSE)
library(erer)

arrayQualityMetrics::arrayQualityMetrics(
  eset,
  outdir = "3_trim_placenta",
  force = TRUE,
  intgroup = 'accession'
)


#
# sub_exprs = mrgd
# sub_pdata = pdata
eset = ExpressionSet(as.matrix(sub_exprs))
eset@phenoData = AnnotatedDataFrame(sub_pdata)

getwd()
arrayQualityMetrics::arrayQualityMetrics(
  eset,
  do.logtransform = FALSE,
  outdir = "smoking",
  force = TRUE,
  intgroup = 'secondaryaccession'
)

sub_pdata$Array.Data.File
# # 
# # 


# for (level in levels(sub_pdata$Biological.Specimen)){
#   sub_pdata = pdata[which(pdata$Diagnosis=='Healthy' & pdata$Biological.Specimen=='Chorion' & 
#                       (pdata$trim=='First Trimester' | pdata$trim=='Second Trimester')),]
#   
#   write.table(sub_pdata,"sub_pdata.tsv", sep="\t", quote=FALSE)
#   sub_pdata = read.table("sub_pdata.tsv", sep="\t", header=TRUE)
# 
#   propercolnames = as.character(make.names(sub_pdata$Array.Data.File))
#   
#   # sub_exprs = exprs
#   colnames(sub_exprs) = as.character(colnames(sub_exprs))
#   sub_exprs = sub_exprs[,propercolnames]
#   colnames(sub_exprs) == propercolnames
#   
#   # sort of takes sub_pdata and sub_exprs as inputs
#   source('scripts/diff_exp.R')
# }
