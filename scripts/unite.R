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


microarray_platform = 'illumina'

# rawspath = 'raws/affymetrix'
# prepath = 'preprocessed/affymetrix'
# pdatapath = 'pdata/'
# plotsqcpath = paste(getwd(), 'plots/qc', sep='/')
# mappedpath = 'mapped/affymetrix'
# mergedpath = 'merged'


rawspath = 'raws/illumina'
prepath = 'preprocessed/illumina'
pdatapath = 'pdata/illumina'
plotsqcpath = paste(getwd(), 'plots/illumina/qc', sep='/')
mappedpath = 'mapped/illumina'
mergedpath = 'merged'


# Load studies description
# studies <- read.table("general/affymetrix_placenta_studies.tsv", header = TRUE, sep = "\t")
studies <- read.table("general/illumina_placenta_studies.tsv", header = TRUE, sep = ",")

all_studies = read.table("general/accession_platform", header = TRUE, sep = '\t')
all_studies




# load sample metadata from csv.file from IGEA database which is based on ArrayExpress data mostly
igea = read.table('igea_tsv/samples.csv',header = TRUE, sep = ',', fill = TRUE)

decidua_exps = unique(igea[igea$Biological.Specimen=="Decidua",]$accession)

# get list of gene expression matrix files
exprs_files = list.files(mappedpath)
exprs_files
mrgd=NULL
mrgd = read.table(paste(mappedpath, exprs_files[1], sep = '/'), header = TRUE, sep = '\t')

exprs_files[1]

# read first expression matrix
current_exprs = read.table(paste(mappedpath, exprs_files[1], sep = '/'), header = TRUE, sep = '\t')

# read each expression matrix and concatenate its columns, leaving only common rows
for (exprs_file in exprs_files[2:length(exprs_files)]){
  current_exprs = read.table(paste(mappedpath, exprs_file, sep = '/'), header = TRUE, sep = '\t')
  # print(nrow(current_exprs))
  mrgd = merge(mrgd, current_exprs, by = "row.names")
  rownames(mrgd) = mrgd$Row.names
  mrgd = mrgd[,!(colnames(mrgd) == "Row.names")]
  print(paste(exprs_file, nrow(mrgd), sep = ' '))
}

# for affymetrix only
# pdata = igea[make.names(igea$Array.Data.File) %in% colnames(mrgd),]

# arrayexpress does not store processed exprs sample name column in a single place
# so this column was created manually for illumina
make.names(igea$exprs_column_names)
# arrayexpress stores affymetrix exprs raw file name in columns in Array.Data.File
make.names(igea$Array.Data.File)

# merge exprs samples file/colum names for illumina and affymetrix into a single column
new_column = character()
for (i in 1:nrow(igea)) {
  if (igea$Array.Data.File[i] != "_"){
    value = as.character(igea$Array.Data.File[i])
  } else {
    value = as.character(igea$exprs_column_names[i])
  }
  new_column = c(new_column, value)
}
new_column
igea$arraydatafile_exprscolumnnames = new_column

# get sample metadata only for expression data in mrgd
pdata = igea[make.names(igea$arraydatafile_exprscolumnnames) %in% colnames(mrgd),]
# read-write metadata to elliminate unneeded factor levels
write.table(pdata,file.path("temp","pdata.tsv"), sep="\t", quote=FALSE)
pdata = read.table(file.path("temp","pdata.tsv"), sep="\t", header=TRUE)

# check how many samples we have in expression data and metadata
# must be equal
ncol(mrgd)
nrow(pdata)

pdata$Diagnosis=="Healthy"
pdata = pdata[pdata$Diagnosis == "Healthy",]

# outliers = c("X11761", "X11420", "X11670", "X13497", "X13521", "X14258")
# yes yes no, yes.., no, yes!
# pdata = pdata[which((pdata$arraydatafile_exprscolumnnames %in% outliers)| pdata$accession=="E-GEOD-73685"),]

# pdata = pdata[which(pdata$accession=="E-GEOD-73685" &
#                     pdata$Diagnosis=="Healthy" &
#                       (pdata$Biological.Specimen == "Decidua" |
#                        pdata$Biological.Specimen == "Maternal Blood" |
#                        pdata$Biological.Specimen == "Umbilical Cord Blood")
# ),]
# pdata = pdata[which(pdata$Diagnosis=="Healthy" & 
#                     pdata$Gestational.Age.Category!="First Trimester" & pdata$Gestational.Age.Category!="Second Trimester" &
#                     pdata$Biological.Specimen!="Umbilical Cord Blood" &
#                     pdata$Biological.Specimen!="Maternal Blood" &
#                     pdata$Biological.Specimen!="Amnion" &
#                     pdata$Biological.Specimen!="Lower Segment" &
#                     pdata$Biological.Specimen!="Uterus Fundus")
# ,]

# allign expression data in mrgd with metadata in pdata
mrgd= mrgd[,make.names(pdata$arraydatafile_exprscolumnnames)]
# check colnames in mrgd match col arraydatafile_exprscolumnnames in pdata
setdiff(colnames(mrgd), make.names(pdata$arraydatafile_exprscolumnnames))
# 
write.table(pdata,file.path("temp","pdata.tsv"), sep="\t", quote=FALSE)
write.table(mrgd, file.path("temp", "mrgd.tsv"), sep="\t", quote=FALSE)
# 
nrow(mrgd)



# create Label Classes for CIBERSORT algorigthm for tissue mixture deconvolution
pdata = read.table(file.path("temp","pdata.tsv"), sep="\t", header=TRUE)
pdata$arraydatafile_exprscolumnnames
pdata$Expression.Data.ID
# decidua 1 1 1 1 0 0 0 0 0 
# blood1  0 0 0 0 1 1 0 0 0
# blood2  0 0 0 0 0 0 1 1 1
labels_df = data.frame(row.names = levels(as.factor(pdata$Cluster)))

for (tissue in rownames(labels_df)){
  for (sample in pdata$Expression.Data.ID) {
    if (pdata[pdata$Expression.Data.ID==sample,]$Cluster == tissue){
      labels_df[tissue,sample] = 1
    } else {
      labels_df[tissue,sample] = 2
    }
  }
}


labels_df
colnames(labels_df) == colnames(mrgd)

decidua_samples_pdata = pdata[pdata$Biological.Specimen=="Decidua",]$arraydatafile_exprscolumnnames 
decidua_samples_labels_df = colnames(labels_df["Decidua",which(labels_df["Decidua",] == 1)])
decidua_samples_labels_df
decidua_samples_pdata == decidua_samples_labels_df

write.table(labels_df,file.path("temp","labels_reference_cell_types_ctb123_evt4_dendric5_excluded6_stb7.tsv"), sep="\t", quote=FALSE)




 
# write.table(mrgd,file.path("temp","mrgd.tsv"), sep="\t", quote=FALSE)
# 
# setdiff(colnames(mrgd), make.names(igea$arraydatafile_exprscolumnnames))
# 
# 
# # end, next ask.R

