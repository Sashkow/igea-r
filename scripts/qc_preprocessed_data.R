setwd('/home/sashkoah/a/article-microarrays')
getwd()

rawspath = 'raws/affymetrix'
prepath = 'preprocessed/affymetrix/'
pdatapath = 'pdata/'
plotsqcpath = paste(getwd(), 'plots/qc', sep='/')



# Load studies description
studies <- read.table("general/affymetrix_placenta_studies.tsv", header = TRUE, sep = "\t")
i = 1
path = paste(prepath, studies[i,]$accession, '_preprocessed_affymetrix.tsv', sep='')

# load IGEA phenodata
igea = read.table('igea_tsv/samples.tsv',header = TRUE, sep = '\t', fill = TRUE)


# hacky way to add CEL file name column to igea by inferring file names from sample names
celfilenames = lapply(as.character(igea$Sample.Name), function(x) paste(strsplit(x, " ")[[1]][1], ".CEL", sep = ""))

igea$Cel.File = as.factor(unlist(celfilenames))


#read preprocessed data into ExpressionSetw
exprs = read.table(path, header = TRUE, sep = '\t')

z = AnnotatedDataFrame(igea[igea$Cel.File %in% colnames(exprs),])

rownames(z) = colnames(exprs)

ncol(exprs)
nrow(z)

exprs.matrix = data.matrix(exprs)
preData = ExpressionSet(assayData = exprs.matrix, phenoData = z)

arrayQualityMetrics::arrayQualityMetrics(
  preData,
  outdir = paste(plotsqcpath,studies$accession[[i]], sep='/'),
  force = TRUE,
  intgroup = 'Diagnosis'
)
