library(org.Hs.eg.db)
setwd('/home/sashkoah/a/r/igea-')

filename = 'mixture_symbol_E-GEOD-60438.tsv'
mixtures = read.csv("/home/sashkoah/a/r/igea-r/sub_exprs.tsv", sep='\t')
mixtures
# reference = read.csv('GSE89497_reference.tsv', sep = '\t')
# symbols <- rownames(reference)
# map = select(org.Hs.eg.db, keys = ,"ENTREZID", "ALIAS")

annotation <- AnnotationDbi::select(org.Hs.eg.db, rownames(mixtures),c("SYMBOL"), keytype = "ENTREZID" )

rownames(mixtures) = annotation$SYMBOL
mixtures

sub_pdata$arraydatafile_exprscolumnnames =  make.names(sub_pdata$arraydatafile_exprscolumnnames) 

colnames(mixtures) = str_replace_all(colnames(mixtures),"\\.","_")


write.table(mixtures, "/home/sashkoah/a/r/igea-r/sub_exprs.tsv", sep = '\t', quote = FALSE)



