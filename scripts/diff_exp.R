
# 
# pd = sub_pdata
# exprs = sub_exprs

nrow(sub_pdata)
ncol(sub_exprs)



colnames(sub_exprs) == pd$Array.Data.File


# unipd$accession

# i=7
# exprs = read.table(paste(mappedpath, '/', studies$accession[[i]], "_mapped_affymetrix.tsv", sep=""), sep="\t", header=TRUE, row.names = 1 )

write.table(sub_pdata,"sub_pdata.tsv", sep="\t", quote=FALSE)
sub_pdata = read.table("sub_pdata.tsv", sep="\t", header=TRUE)
propercolnames = as.character(make.names(sub_pdata$Array.Data.File))

colnames(sub_exprs) = as.character(colnames(sub_exprs))
sub_exprs = sub_exprs[,propercolnames]
colnames(sub_exprs) == propercolnames


write.table(sub_exprs, "sub_exprs.tsv", sep="\t", quote=FALSE)
write.table(sub_pdata, "sub_pd.tsv", sep="\t", quote=FALSE)




# age = as.data.frame(pd$Gestational.Age.Category)
# as.character(pd$Gestational.Age.Category)
# pd$trim = with(pd, ifelse(Gestational.Age.Category == "Term" | Gestational.Age.Category == "Late Preterm", "Third Trimester", as.character(Gestational.Age.Category)))
# pd$trim

# differential expression analysis

# colnames(fit_mod)

# check1 = makeContrasts(b-a, levels=fit_mod)
# check2 = makeContrasts(a-b, levels=fit_mod)
# contrast.matrix = check2

# contrast.matrix <- makeContrasts(b-a,c-b,c-a, levels=fit_mod)

fit_mod = model.matrix(~ 0 + as.factor(trim), data=sub_pdata)

colnames(fit_mod)
colnames(fit_mod) = c('a','b')
# colnames(fit_mod) = c('healthy', 'preeclampsia')
contrast.matrix <- makeContrasts(b-a, levels=fit_mod)
fit <- lmFit(sub_exprs, fit_mod)  # fit each probeset to model
fit2 <- contrasts.fit(fit, contrast.matrix)
efit <- eBayes(fit2)        # empirical Bayes adjustment
# 


t = topTable(efit, number = nrow(sub_exprs))
t = t[order(-abs(t$logFC)),]



d = t

nrow(d)
d$adj.P.Val
d = d[d$adj.P.Val<.4,]
nrow(d)

diff = d[abs(d$logFC)>1.5,]
nrow(diff)
diff = d
library(org.Hs.eg.db)
sel = AnnotationDbi::select(org.Hs.eg.db, rownames(diff), c("SYMBOL","GENENAME"))
nrow(sel)


specimen = levels(sub_pdata$Biological.Specimen)[1]

trims = str_replace_all(paste(levels(as.factor(sub_pdata$trim)), collapse = '__'), ' ', '_')
filename = paste(specimen,trims, sep='__')

write.table(sel, file.path('healthy_diff_lists',filename), sep = '\t', quote=FALSE)


# write.table(diff, "/home/sashkoah/a/r/article-microarrays/diff12.tsv", sep="\t", quote=FALSE)


# entrez_identifiers = as.character(rownames(diff))
# anno = select(org.Hs.eg.db, entrez_identifiers, columns = "SYMBOL", keytype = "ENTREZID")
# write.table(anno, "/home/sashkoah/a/r/article-microarrays/differential_expression_from_data/3trim_placenta/3trim_placenta_symbol.tsv", sep="\t", quote=FALSE)
# 
# 
# write.table(diff, "/home/sashkoah/a/r/article-microarrays/differential_expression_from_data/3trim_placenta/3trim_placenta_.tsv", sep="\t", quote=FALSE)
# ba = d[abs(d$b...a)>2,]
# nrow(ba)
# 
# cb = d[abs(d$c...b)>2,]
# nrow(cb)
# ca = d[abs(d$c...a)>2,]
# nrow(ca)
# 
# write.table(cb,"cb.tsv", sep="\t", quote=FALSE)

# 
# # how a differs from b and c
# # how b differs from a and c
# # how c differs from a and b
# a = intersect(rownames(ba),rownames(ca))
# b = intersect(rownames(ba),rownames(cb))
# c = intersect(rownames(ca),rownames(cb))
# length(a)
# length(b)
# length(c)
# 
# d = d[with(d, order(abs(b...a))),]
# nrow(d)
# 
# 
# 
# t = t[order(abs(t$as.factor.cluster.2), decreasing = TRUE),]
# t12 = t12[order(abs(t12$logFC), decreasing = TRUE),]
# # t = t[order(-t$logFC),]
# 
# head(tall)
# head(t12)
# t12 = t
# tall = t
# tall
# t12


# 
# 
# tabl = topTable(efit, number = 100, coef = "as.factor(Biological.Specimen)Chorion")      # table of differentially expressed probesets
# tabl = tabl[order(-tabl$logFC),]
# chorion_genes = rownames(tabl)
# 
# 
# tabl = topTable(efit, number = 100, coef = "as.factor(Biological.Specimen)Decidua")      # table of differentially expressed probesets
# tabl = tabl[order(-tabl$logFC),]
# decidua_genes=  rownames(tabl)
# 
# tabl = topTable(efit, number = 100, coef = "as.factor(Biological.Specimen)Placenta")      # table of differentially expressed probesets
# tabl = tabl[order(-tabl$logFC),]
# placenta_genes =  rownames(tabl)
# 
# 
# 
# library(RAM)
# 
# 
# foo <- c('a','b','c','d')
# baa <- c('a','e','f','g')
# group.venn(list(foo=foo, baa=baa), label=TRUE,
#            fill = c("orange", "blue"),
#            cat.pos = c(0, 0),
#            lab.cex=1.1)
# 
# group.venn(list(chorion_genes=chorion_genes,
#                 decidua_genes=decidua_genes,
#                 placenta_genes=placenta_genes),
#            label=TRUE,
#            fill = c("orange", "blue"),
#            lab.cex=1.1)
# 
# 
# 
# 
# chord = intersect(chorion_genes, decidua_genes)
# 
# 
# 
# intersect(placenta_genes, chorion_genes)
# intersect(placenta_genes, decidua_genes)
# 
# genes = exclusive_or(chorion_genes, decidua_genes)
# genes2 = exclusive_or_3(chorion_genes, decidua_genes,placenta_genes)
# setdiff(genes, genes2)
# paste(genes, collapse = ", ")
# 
# # sort by decending of differential expression
# 
# 
# 
# 
# rownames(tabl)
# 
# exprs = exprs[genes,]
# 
# View(exprs['6736',])




