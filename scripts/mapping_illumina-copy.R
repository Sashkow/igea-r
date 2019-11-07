# install.packages("BiocManager")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")


library(stringr)
source('/home/sashkoah/a/r/igea-r/scripts/aln_illumina/degs_utils.R')


library(illuminaHumanv4.db)

setwd('/home/sashkoah/a/r/igea-r')


getwd()

rawspath = 'raws/illumina'
prepath = 'preprocessed/illumina'
mappedpath = 'mapped/illumina'
datapath = 'pdata/'
plotsqcpath = paste(getwd(), 'plots/qc', sep='/')

# source("https://bioconductor.org/biocLite.R")


# Load studies description
studies <- read.table("general/illumina_placenta_studies.tsv", header = TRUE, sep = ",")

# load IGEA phenodata
igea = read.table('igea_tsv/smoking_preeclampsia.csv',header = TRUE, sep = ',', fill = TRUE)


# 60438_preprocessed_illumina_1.tsv - illuminaHumanv4.db
# 60438_preprocessed_illumina_2.tsv - illuminaHumanv4.db


exprs_path = "/home/sashkoah/a/r/igea-r/preprocessed/illumina/GSE27272_from_raw/GSE27272_54_24525.tsv"

exprs_filename = basename(exprs_path)
exprs_filename = str_replace(exprs_filename,'\\.tsv','_mapped\\.tsv') 
exprs_filename

exprs = read.table(exprs_path, header = TRUE, sep = "\t")
dim(exprs)
exprs


# rownames(exprs) = exprs$Reporter.Identifier
# colnames(exprs)
# exprs=sub_exprs


# drops <- c("X","merged", "Reporter.Identifier")
# exprs = exprs[ , !(names(exprs) %in% drops)]
ncol(exprs)
nrow(exprs)
# exprs.save = exprs
# map probes to engrez gene identifiers in exprs


# file.db = illuminaHumanv4.db

file.db = illuminaHumanv3.db


rownames(exprs)
# exprs = getUniqueProbesets(exprs,"illuminaHumanv3")

# if no .db file provided
# file.txt = read.table('/home/sashkoah/a/r/igea-r/annotations/illumina/HumanWG-6_V3_0_R3_11282955_A_probe_id_entrez.txt', header = TRUE, sep = "\t")
# file.txt = read.table('/home/sashkoah/a/r/igea-r/annotations/illumina/A-MEXP-930.adf_Illumina_Human-6_v2_Expression BeadChip_probe_id_entrez.txt', header = TRUE, sep = "\t")
file.txt = read.table('/home/sashkoah/a/r/igea-r/annotations/illumina/A-MEXP-1172.adf_Illumina_HumanRef_8_v3.0_Expression_BeadChip.txt', header = TRUE, sep = "\t")
# file.txt = read.table('/home/sashkoah/a/r/igea-r/annotations/illumina/GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt', header = TRUE, sep = "\t", fill = TRUE)

# change file.txt inside this function before using
# exprs = getUniqueProbesetsTxt(exprs)

require(WGCNA)
## Get probeset to entrezid mapping
probesetsID <- rownames(exprs)
# probesetsID_EntrezID<-select(get(paste(platform, ".db", sep="")), probesetsID, "ENTREZID")
# file.txt = read.table('/home/sashkoah/a/r/igea-r/annotations/illumina/HumanWG-6_V3_0_R3_11282955_A_probe_id_entrez.txt', header = TRUE, sep = "\t")
file.txt = read.table('/home/sashkoah/a/r/igea-r/annotations/illumina/A-MEXP-1172.adf_Illumina_HumanRef_8_v3.0_Expression_BeadChip.txt', header = TRUE, sep = "\t")
# file.txt = read.table('/home/sashkoah/a/r/igea-r/annotations/illumina/GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt', header = TRUE, sep = "\t", fill = TRUE)

# file.txt = read.table('/home/sashkoah/a/r/igea-r/annotations/illumina/A-MEXP-930.adf_Illumina_Human-6_v2_Expression BeadChip_probe_id_entrez.txt', header = TRUE, sep = "\t")

probesetsID_EntrezID<-file.txt

probesetsID_EntrezID

## Replace probesetsIDs with gene IDs in expression data matrix

# Exclude NA probesets
probesetsID_EntrezID$ENTREZID = probesetsID_EntrezID$Reporter.Database.Entry.entrez.
probesetsID_EntrezID$PROBEID = probesetsID_EntrezID$Reporter.Name

probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]

# Exclude probesets mapped to different genes simultaneously
n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
probesetsID_EntrezID

uniques <- n_occur[n_occur$Freq == 1,]$Var1
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]
nrow(probesetsID_EntrezID)
# Filter expression matrix based on left probesets
exprs1 <- exprs[which(rownames(exprs) %in% probesetsID_EntrezID$PROBEID),]
probesetsID_EntrezID = probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% rownames(exprs1)),] 
nrow(exprs1)
nrow(probesetsID_EntrezID)


exprs = exprs1
# Select one probeset among the probesets mapped to the same gene based on maximum average value across the samples
collapsed = collapseRows(exprs, probesetsID_EntrezID$ENTREZID, probesetsID_EntrezID$PROBEID, method="MaxMean")  
exprs <- collapsed$datETcollapsed

dim(exprs)
exprs

write.table(exprs, "/home/sashkoah/a/r/igea-r/mapped/new_smoking_symbol/GSE27272_exprs_from_raw_pl_54_mapped.tsv", sep="\t", quote=FALSE)


mapa = AnnotationDbi::select(org.Hs.eg.db, rownames(exprs), columns = "SYMBOL", keytype = "ENTREZID")
FALSE %in% (mapa$ENTREZID == rownames(exprs))
mapa = mapa[which(mapa$SYMBOL!="NA"),]
exprs = exprs[mapa$ENTREZID,]
FALSE %in% (mapa$ENTREZID == rownames(exprs))
rownames(exprs) = mapa$SYMBOL

write.table(exprs, "/home/sashkoah/a/r/igea-r/mapped/new_smoking_symbol/GSE27272_exprs_from_raw_pl_54_mapped_symbol.tsv", sep="\t", quote=FALSE)

e = read.table("/home/sashkoah/a/r/igea-r/mapped/new_smoking_symbol/GSE27272_.tsv", sep="\t")

exprs_7434 = read.table("/home/sashkoah/a/r/igea-r/mapped/smoking/E-GEOD-7434_mapped_affymetrix.tsv", sep = '\t', header = TRUE)
dim(exprs_7434)
mapa = AnnotationDbi::select(org.Hs.eg.db, rownames(exprs_7434), columns = "SYMBOL", keytype = "ENTREZID")

FALSE %in% (mapa$ENTREZID == rownames(exprs_7434))
na_genes = mapa[which(mapa$SYMBOL=="NA"),]$ENTREZ
mapa = mapa[which(mapa$SYMBOL!="NA"),]
nrow(mapa)
exprs_7434 = exprs_7434[mapa$ENTREZID,] 
FALSE %in% (mapa$ENTREZID == rownames(exprs_7434))
rownames(exprs_7434) = mapa$SYMBOL

TRUE %in% duplicated(rownames(exprs_7434))

write.table(exprs_7434, "/home/sashkoah/a/r/igea-r/mapped/new_smoking_symbol/E-GEOD-7434_mapped_affymetrix_symbol.tsv", sep="\t", quote=FALSE)
exprs_7434 = read.table("/home/sashkoah/a/r/igea-r/mapped/new_smoking_symbol/E-GEOD-7434_mapped_affymetrix_symbol.tsv", sep = '\t', header = TRUE)

#done
sub_exprs = exprs
nrow(sub_pdata)
ncol(sub_exprs)

i=3
studies[i,]$secondaryaccession

exprs_file_name = paste(studies[i,]$secondaryaccession,"mapped_exprs.tsv", sep = "_")
pdata_file_name = paste(studies[i,]$secondaryaccession,"mapped_pdata.tsv", sep = "_")
exprs_path = file.path(mappedpath,studies[i,]$secondaryaccession)
if (! dir.exists(exprs_path)){
  dir.create(exprs_path, recursive = TRUE)
}

write.table(sub_exprs, file.path(exprs_path,exprs_file_name), sep = "\t", quote = FALSE)


write.table(exprs, "/home/sashkoah/a/r/igea-r/mapped/smoking/E-GEOD-27272_mapped_exprs.tsv", sep="\t", quote=FALSE)
file.path(mappedpath, exprs_filename)
write.table(exprs, file.path(mappedpath, exprs_filename), sep="\t", quote=FALSE)
# write.table(illuminaHumanv4.db, file.path(mappedpath, "exprs_filename"), sep="\t", quote=FALSE)

