# BiocManager::install("rgl")
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


mrgd = NULL
pdata = NULL
mrgd = read.table(file.path("temp", "mrgd.tsv"), sep="\t", header=TRUE)
pdata = read.table(file.path("temp","pdata.tsv"), sep="\t", header=TRUE)

colnames(pdata)
pdata$arraydatafile_exprscolumnnames
columns = c('Diagnosis',
            'Gestational.Age.Category',
            'accession',
            'Biological.Specimen',
            'Array.Data.File',
            'Gestational.Age', 
            "exprs_column_names",
            "platform_batch",
            "arraydatafile_exprscolumnnames",
            "Fetus.Sex",
            "Gravidity",
            "Parity",
            "medical_sample_name",
            "Caesarean.Section")

pdata = pdata[,columns]
pdata$Array.Data.File = pdata$arraydatafile_exprscolumnnames
propercolnames = as.character(make.names(pdata$Array.Data.File))
rownames(pdata) = propercolnames
propercolnames
colnames(mrgd) = as.character(colnames(mrgd))
propercolnames
mrgd = mrgd[,propercolnames]

colnames(mrgd) == propercolnames
colnames(mrgd) == rownames(pdata)
colnames(mrgd)

pdata$Array.Data.File = colnames(mrgd)
pdata$Biological.Specimen

batch = as.factor(pdata$accession)

batch
as.factor(pdata$Diagnosis)
mod = model.matrix(~ as.factor(Gestational.Age.Category) + as.factor(Diagnosis), data=pdata)
# mod = model.matrix(~as.factor(Gestational.Age.Category), data=pdata)

# exprs=mrgd
exprs = ComBat(dat=as.matrix(mrgd), batch=batch, mod=mod, par.prior=TRUE, prior.plots=TRUE)

pdata$trim = with(pdata, ifelse(Gestational.Age.Category == "Term" | Gestational.Age.Category == "Late Preterm" | Gestational.Age.Category == "Early Preterm", "Third Trimester", as.character(Gestational.Age.Category)))

sub_pdata = pdata
nrow(pdata)
sub_exprs = exprs
# 
pdata$trim
pdata$Biological.Specimen



# sub_pdata$is_outlier = FALSE
# sub_pdata[which(sub_pdata$arraydatafile_exprscolumnnames %in% outlier_elements),]$is_outlier = TRUE
pdata$trim
pdata$Biological.Specimen
sub_pdata = pdata[which(pdata$Diagnosis=="Healthy"),]
# sub_pdata = pdata[which(pdata$Diagnosis=="Healthy" & 
#                           pdata$trim=='Third Trimester' & 
#                           pdata$Biological.Specimen!="Umbilical Cord Blood" &
#                           pdata$Biological.Specimen!="Maternal Blood" &
#                           pdata$Biological.Specimen!="Amnion" &
#                           pdata$Biological.Specimen!="Lower Segment" &
#                           pdata$Biological.Specimen!="Uterus Fundus"),]
# sub_pdata = pdata[which(pdata$Diagnosis=="Pre-Eclampsia"),]

# 
source('scripts/sub_pd_to_sub_exprs.R')

ncol(sub_exprs) 
nrow(sub_pdata)

# how many samples in each dataset
as.data.frame(table(sub_pdata$accession))[all_studies$accession,]



# # sort of takes sub_pdata and sub_exprs as inputs
source('scripts/diff_exp.R')


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
write.table(sub_exprs,"sub_exprs.tsv", sep="\t", quote=FALSE)

eset = ExpressionSet(as.matrix(sub_exprs))
eset@phenoData = AnnotatedDataFrame(sub_pdata)

source(paste(getwd(),'scripts/plots.R',sep='/'))

# colnames(sub_exprs) = lookUp(as.character(colnames(sub_exprs)), 'org.Hs.eg', 'SYMBOL') if 
# sub_exprs = sub_exprs[which(colnames(sub_exprs)!=NA),]
# length(colnames(sub_exprs)==NA)

pca = prcomp(t(sub_exprs))

re = kmeans(pca$x, centers = 2, iter.max = 1000, nstart = 1, trace=FALSE)

# re = kmeans(as.matrix(t(exprs)), centers = 3, iter.max = 1000, nstart = 1, trace=FALSE)
re$size

re$cluster
cluster = as.data.frame(re$cluster)
cluster
sub_pdata$Array.Data.File = as.character(make.names(sub_pdata$Array.Data.File))
sub_pdata$Array.Data.File == rownames(cluster)

sub_pdata[,"cluster"] = cluster$`re$cluster`

nrow(pca$x)
nrow(sub_pdata)


f <- function(sub_pdata) {
  print(sub_pdata)
  if (sub_pdata[3]=='E-GEOD-60438'){
    return('Outlier')
  } else {
    return(sub_pdata[4])
  }
}

sub_pdata$temp = apply(sub_pdata, 1, f)
sub_pdata$temp

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

pl = fviz_pca_ind(pca, axes = c(1,2),
                      label=c("var"),
                      # habillage=sub_pdata$Biological.Specimen,
                      repel = TRUE,
                      title = "GEOD-60438 and GEOD-73685",
                      # addEllipses = TRUE,
                   
) + geom_point(aes(colour=sub_pdata$temp), size=2)
pl

summary(pca)

library(rgl)

levels(sub_pdata$Biological.Specimen)
length(levels(sub_pdata$Gestational.Age.Category))
plot3d(pca$x[,1:3], col=rainbow(length(levels(sub_pdata$temp)))[factor(as.integer(sub_pdata$temp))], size=10)
legend3d("topright", legend = levels(sub_pdata$temp), pch = 16, col = rainbow(length(levels(sub_pdata$temp))), cex=1, inset=c(0.02))

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
eset = ExpressionSet(as.matrix(sub_exprs))
eset@phenoData = AnnotatedDataFrame(sub_pdata)


arrayQualityMetrics::arrayQualityMetrics(
  eset,
  outdir = "60438",
  force = TRUE,
  intgroup = 'Array.Data.File'
)

sub_pdata$Array.Data.File
# # 
# # 


