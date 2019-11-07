pdata = read.table('/home/sashkoah/a/r/igea-r/igea_tsv/pdata.csv', sep = ',',header = TRUE)
exprs = read.table('/home/sashkoah/a/r/igea-r/mapped/illumina/E-GEOD-27272_mapped_correctly_illumina.tsv', sep = '\t',header = TRUE)
colnames(exprs)
ncol(exprs)

pdata$XRagul.Sample.Name
pdata$XRagul.Sample.Name = as.character(pdata$XRagul.Sample.Name)
pdata = pdata[which(pdata$XRagul.Sample.Name %in% colnames(exprs)),]
pdata
nrow(pdata)
exprs = exprs[,pdata$XRagul.Sample.Name]
colnames(exprs) == pdata$XRagul.Sample.Name
colnames(exprs) = pdata$geo_accession

colnames(exprs)
write.table(exprs, "/home/sashkoah/a/r/igea-r/mapped/illumina/E-GEOD-27272_mapped_correctly.tsv", sep="\t", quote=FALSE)