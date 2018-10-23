library(affy) 
library(ArrayExpress)
library(arrayQualityMetrics)
library(sva)
library(limma)
library(stringr)

source(paste(getwd(),'scripts/plots.R',sep='/'))
plotsqcpath = paste(getwd(), 'plots/qc', sep='/')



mrgd = read.table("mrgd.tsv", sep="\t", header=TRUE)
pdata = read.table("pdata.tsv", sep="\t", header=TRUE)

columns = c('Diagnosis', 'Gestational.Age.Category', 'accession', 'Biological.Specimen','Array.Data.File', 'Gestational.Age')
pdata = pdata[,columns]
propercolnames = as.character(make.names(pdata$Array.Data.File))
colnames(mrgd) = as.character(colnames(mrgd))
mrgd = mrgd[,propercolnames]
colnames(mrgd) == propercolnames

pdata$Biological.Specimen


batch = as.factor(pdata$accession)
mod = model.matrix(~ as.factor(Gestational.Age.Category) + as.factor(Diagnosis), data=pdata)
# mod = model.matrix(~as.factor(Biological.Specimen), data=pdata)



exprs = ComBat(dat=mrgd, batch=batch, mod=mod, par.prior=TRUE, prior.plots=TRUE)

pdata$trim = with(pdata, ifelse(Gestational.Age.Category == "Term" | Gestational.Age.Category == "Late Preterm" | Gestational.Age.Category == "Early Preterm", "Third Trimester", as.character(Gestational.Age.Category)))

sub_pdata = pdata
sub_exprs = exprs

sub_pdata$trim

pdata$Diagnosis



# sub_pdata = pdata[which(pdata$Diagnosis == 'Pre-Eclampsia' & pdata$Gestational.Age != '_'),]
# sub_pdata = sub_pdata[as.double(as.character(sub_pdata$Gestational.Age)) > 34,]

# mean(as.double(as.character(sub_pdata$Gestational.Age)))
# median(as.double(as.character(sub_pdata$Gestational.Age)))

# sub_pdata_healthy = pdata[which(pdata$Diagnosis == 'Healthy' & pdata$Gestational.Age != '_'),]
# sub_pdata_healthy = sub_pdata_healthy[as.double(as.character(sub_pdata_healthy$Gestational.Age)) %in% as.double(as.character(sub_pdata$Gestational.Age)),]
# sub_pdata_healthy$Gestational.Age
# 
# sub_pdata = rbind(sub_pdata, sub_pdata_healthy)
pdata$trim

sub_pdata = pdata[which(pdata$Diagnosis=='Healthy' & pdata$Biological.Specimen=='Chorion' & 
                          (pdata$trim=='First Trimester' | pdata$trim=='Second Trimester')),]

source('scripts/sub_pd_to_sub_exprs.R')
# sort of takes sub_pdata and sub_exprs as inputs
source('scripts/diff_exp.R')


for (level in levels(sub_pdata$Biological.Specimen)){
  sub_pdata = pdata[which(pdata$Diagnosis=='Healthy' & pdata$Biological.Specimen=='Chorion' & 
                      (pdata$trim=='First Trimester' | pdata$trim=='Second Trimester')),]
  
  write.table(sub_pdata,"sub_pdata.tsv", sep="\t", quote=FALSE)
  sub_pdata = read.table("sub_pdata.tsv", sep="\t", header=TRUE)

  propercolnames = as.character(make.names(sub_pdata$Array.Data.File))
  
  # sub_exprs = exprs
  colnames(sub_exprs) = as.character(colnames(sub_exprs))
  sub_exprs = sub_exprs[,propercolnames]
  colnames(sub_exprs) == propercolnames
  
  # sort of takes sub_pdata and sub_exprs as inputs
  source('scripts/diff_exp.R')
}

unique(sub_pdata$Gestational.Age.Category)



sub_pdata = pdata[which((pdata$Gestational.Age.Category == 'Term' | pdata$Gestational.Age.Category == 'Early Preterm' | pdata$Gestational.Age.Category == 'Late Preterm') & (pdata$Biological.Specimen=='Placenta' | pdata$Biological.Specimen=='Chorion')), ]
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

write.table(sub_pdata,"sub_pdata.tsv", sep="\t", quote=FALSE)
sub_pdata = read.table("sub_pdata.tsv", sep="\t", header=TRUE)

sub_pdata$Array.Data.File
propercolnames = as.character(make.names(sub_pdata$Array.Data.File))

# sub_exprs = exprs
colnames(sub_exprs) = as.character(colnames(sub_exprs))
sub_exprs = sub_exprs[,propercolnames]
colnames(sub_exprs) == propercolnames

write.table(sub_exprs,"sub_exprs.tsv", sep="\t", quote=FALSE)


eset = ExpressionSet(as.matrix(sub_exprs))
eset@phenoData = AnnotatedDataFrame(sub_pdata)

source(paste(getwd(),'scripts/plots.R',sep='/'))
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

sub_pdata$h
pl <- pcaPlots(pca, sub_pdata, c('Diagnosis', 'Gestational.Age.Category', 'accession', 'Biological.Specimen', 'cluster', 'trim', 'Highly.mixed', 'Highly.pure'), ncol=2)

dir.create(paste(plotsqcpath,"manual", sep='/'), showWarnings = FALSE)

save_plot(paste(plotsqcpath, "manual", "all_mixed_pure.png", sep='/'),
          base_height=5, base_aspect_ratio = pl[[2]], pl[[1]], nrow=2,limitsize = FALSE)


arrayQualityMetrics::arrayQualityMetrics(
  eset,
  outdir = "3_trim_placenta",
  force = TRUE,
  intgroup = 'accession'
)
# 
# eset = ExpressionSet(as.matrix(mrgd))
# eset@phenoData = AnnotatedDataFrame(pd)
# 
# 
# 
# 
# arrayQualityMetrics::arrayQualityMetrics(
#   eset,
#   outdir = "before_combat_all",
#   force = TRUE,
#   intgroup = 'accession'
# )

