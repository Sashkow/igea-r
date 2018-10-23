write.table(sub_pdata,"sub_pdata.tsv", sep="\t", quote=FALSE)
sub_pdata = read.table("sub_pdata.tsv", sep="\t", header=TRUE)

propercolnames = as.character(make.names(sub_pdata$Array.Data.File))

# sub_exprs = exprs
colnames(sub_exprs) = as.character(colnames(sub_exprs))
sub_exprs = sub_exprs[,propercolnames]
colnames(sub_exprs) == propercolnames
