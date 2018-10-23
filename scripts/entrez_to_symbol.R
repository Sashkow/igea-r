library(org.Hs.eg.db)
setwd('/home/sashkoah/a/r/article-microarrays/pure_mixed')
filename = 'exprs_healthy_all.tsv'
mixtures = read.csv(filename, sep='\t')
# reference = read.csv('GSE89497_reference.tsv', sep = '\t')
# symbols <- rownames(reference)
# map = select(org.Hs.eg.db, keys = ,"ENTREZID", "ALIAS")
mapa = select(org.Hs.eg.db, rownames(mixtures), columns = "SYMBOL", keytype = "ENTREZID")
rownames(mixtures) = mapa$SYMBOL
write.table(mixtures, filename, sep = '\t', quote = FALSE)
