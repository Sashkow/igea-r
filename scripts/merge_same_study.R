library(sva)
setwd('/home/sashkoah/a/r/article-microarrays')

path1 = 'preprocessed/illumina/E-GEOD-60438_preprocessed_illumina_1.tsv'
path2 = 'preprocessed/illumina/E-GEOD-60438_preprocessed_illumina_2.tsv'

exprs1 = read.table(path1, header = TRUE, sep = '\t')
exprs2 = read.table(path2, header = TRUE, sep = '\t')

ab = intersect(rownames(exprs1), rownames(exprs2))
nrow(exprs1)
nrow(exprs2)
length(ab)
length(ab) == nrow(exprs1)



