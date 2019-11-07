library(data.table)
library(openxlsx)
# library(xlsx)
# pd = sub_pdata
# exprs = sub_exprs
# /home/sashkoah/a/r/igea-r/differential_expression_from_data/lumiaffy/trims_healthy
nrow(sub_pdata)
ncol(sub_exprs)
lima
sub_pdata$Biological.Specimen

colnames(sub_exprs) == sub_pdata$Array.Data.File

# xxx <- as.data.frame(sub_exprs)
# library(org.Hs.eg.db)
# sel = AnnotationDbi::select(org.Hs.eg.db, rownames(xxx), c("SYMBOL","GENENAME"))
# sel = sel[!is.na(sel$SYMBOL),]
# xxx <- xxx[which(rownames(xxx) %in% sel$ENTREZID),]
# rownames(xxx) <- sel$SYMBOL
# 
# xxx <- as.data.frame(t(xxx))
# xxx$trim <- sub_pdata$trim_term
# colnames(xxx)

# model <- glm(LEP~trim, data=xxx)
# base::summary(model)

#difexp analysis two groups without contrasts
# fit_mod = model.matrix(~ as.factor(sub_pdata$trim_term), data=sub_pdata)
# colnames(fit_mod) = c('i','ii')
# fit <- lmFit(sub_exprs,fit_mod)
# fit <- eBayes(fit)
# t = topTable(fit, adjust="BH",coef= "ii", number = nrow(sub_exprs))



as.factor(sub_pdata$trim_term)
fit_mod = model.matrix(~ 0 + as.factor(sub_pdata$trim_term), data=sub_pdata)
colnames(fit_mod)
colnames(fit_mod) = c('i','ii','iii')
# colnames(fit_mod) = c('i','ii')
fit <- lmFit(sub_exprs, fit_mod)  # fit each probeset to model
contrast.matrix <- makeContrasts(ii-i,iii-ii, levels=fit_mod)
# contrast.matrix <- makeContrasts(ii-i, levels=fit_mod)
fit2 <- contrasts.fit(fit, contrast.matrix)
efit <- eBayes(fit2)        # empirical Bayes adjustment


decideTests()
methods = c("separate", "global", "hierarchical","nestedF")

# method = methods[4]
# results <- decideTests(efit, method = method)
# results = as.data.frame(results)
# results4= rownames(results[which(abs(results$`ii - i`)==1 & abs(results$`iii - ii`)==1),])

library(RAM)
# results_intersect = Reduce(intersect, list(results1,results2, results3, results4))


intersect.Vector()

group.venn(list(separate=results1, global=results2,hierarchical=results3, nestedF=results4), label=FALSE,
           fill = c("yellow", "blue", "red", "black"))



summary(results)
vennDiagram(results)
dfresults = as.data.frame(results@.Data)

vennCounts(dfresults)
# efit_just_two = efit
# efit_all_three = efit

# efit = efit_all_three


t = topTable(efit, number = nrow(sub_exprs))

# t_results = merge(t,dfresults, by = "row.names")
# rownames(t_results) = t_results$Row.names
# t_results = t_results[,!(colnames(t_results) == "Row.names")]
# nrow(t_results)
# rownames(t_results)

d = t[which(rownames(t) %in% results_intersect),]

nrow(t)
rownames(t)


# t = topTable(efit, coef = 1, number = nrow(sub_exprs))

# t = t[order(-abs(t$logFC)),]
# t = t[order(abs(t$adj.P.Val )),]
# t = t[order(-abs(t$p...h)),]
t

# d = t_results
d
nrow(d)

d$adj.P.Val

# d = d[d$P.Value <0.05,]

# d = d[d$adj.P.Val<0.05,]


nrow(d)
# d1 = d
# d2 = d



# dnew = d
diff = d
# d$ii...i
# diff = d[which(abs(d$ii...i)>1),]
# diff = d[which(abs(d$ii...i)< -1 & abs(d$iii...ii)>1),]
# diff
# 
# diff = d[abs(d$logFC)>1,]

# diff = d[which(abs(d$p...h)>.4 & abs(d$s...h)>.3),]
# rownames(diff[1:173,])
nrow(diff)
rownames(diff)
diff

cyp_genes = diff[grep(rownames(diff)),]
rownames(diff[grep("CYP",rownames(diff)),])
cyp_genes
# 
# gene    2-1      3-2
# HBE1   -2.761074 -3.744563
# HSD3B2  1.859300 -1.689840 
# LEP    -4.121517 2.499421 
# NRN1    1.912029 1.522200 
# PF4     1.925063 -1.66373



library(org.Hs.eg.db)

rownames(diff)
keys(org.Hs.eg.db)
sel = AnnotationDbi::select(org.Hs.eg.db, rownames(diff), c("SYMBOL","GENENAME"))
sel = sel[!is.na(sel$SYMBOL),]

sel2 = AnnotationDbi::select(org.Hs.eg.db, rownames(sub_exprs), c("SYMBOL","GENENAME"))
sel2 = sel2[!is.na(sel2$SYMBOL),]

# if (nrow(sel2[is.na(sel2$SYMBOL),])==1){
#   sel2[is.na(sel2$SYMBOL),]$SYMBOL = "C17orf47"
# }

length(unique(rownames(sub_exprs)))==length(rownames(sub_exprs))


length(sel2$SYMBOL) == length(unique(sel2$SYMBOL))

NA %in% sel2$SYMBOL
NA %in% rownames(sub_exprs)
sub_exprs = sub_exprs[sel2$ENTREZID,]

nrow(sub_exprs)== nrow(sel2)
!(FALSE %in% (rownames(sub_exprs) == sel2$ENTREZID))
nrow(diff)
rownames(sub_exprs) = sel2$SYMBOL
# write.table(sub_exprs, "merged_for_article_3_exprs.tsv", sep= '\t')

merged_sel = merge(sel,diff, by.x = "ENTREZID", by.y = "row.names")
nrow(merged_sel)


for (i in 1:length(levels(as.factor(sub_pdata$trim_term))))
{
  
  trim = levels(as.factor(sub_pdata$trim_term))[i]
  trim
  sub_exprs[which(rownames(sub_exprs) %in% merged_sel$SYMBOL),]
  difexp_exprs = sub_exprs[which(rownames(sub_exprs) %in% merged_sel$SYMBOL),]
  
  setdiff(rownames(difexp_exprs), merged_sel$SYMBOL)
  merged_sel[is.na(merged_sel$SYMBOL),]
  nrow(difexp_exprs)
  nrow(merged_sel)
  !(FALSE %in% (merged_sel$SYMBOL == rownames(difexp_exprs)))
  trim_sample_names = sub_pdata[sub_pdata$trim_term==trim,]$arraydatafile_exprscolumnnames
  difexp_exprs_trim = difexp_exprs[,which(colnames(difexp_exprs) %in% trim_sample_names)]
  as.numeric(difexp_exprs_trim[which(rownames(difexp_exprs_trim)=="NAT2"),])
  
  merged_sel[,paste(trim,"Average",sep = " ")] = rowMeans(difexp_exprs_trim)
  merged_sel
}

nrow(merged_sel)

# name_ending = "_2_1_E-GEOD-9984_E-GEOD-37901.tsv"
name_ending = "_term_2_and_2_1_full_with_averages_multiple_testing_methods_intersect.tsv"
write.table(merged_sel, paste("/home/sashkoah/a/r/igea-r/article_3/placenta/difexp", name_ending,sep=''), sep= '\t', row.names = FALSE)
write.table(sub_exprs, paste("/home/sashkoah/a/r/igea-r/article_3/placenta/sub_exprs", name_ending,sep=''), sep= '\t')
sub_sub_pdata = sub_pdata[,c("secondaryaccession","trim_term","Gestational.Age.Appr","Deviation.Gestational.Age", "Combined.Fetus.Sex","arraydatafile_exprscolumnnames")]
write.table(sub_sub_pdata, paste("/home/sashkoah/a/r/igea-r/article_3/placenta/sub_pdata", name_ending,sep=''), sep= '\t', row.names = FALSE)

library(openxlsx)
name_ending = "_term_2_and_2_1_full_with_averages_multiple_testing_methods_intersect.xlsx"
# name_ending = "_2_1_E-GEOD-9984_E-GEOD-37901.xlsx"

write.xlsx(merged_sel, paste("/home/sashkoah/a/r/igea-r/article_3/placenta/difexp", name_ending,sep=''), sheetName = "Data", 
           col.names = TRUE, row.names = FALSE, append = FALSE)
write.xlsx(sub_exprs, paste("/home/sashkoah/a/r/igea-r/article_3/placenta/sub_exprs", name_ending,sep=''), sheetName = "Data", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(sub_sub_pdata, paste("/home/sashkoah/a/r/igea-r/article_3/placenta/sub_pdata", name_ending,sep=''), sheetName = "Data", 
           col.names = TRUE, row.names = FALSE, append = FALSE)



# filename = paste(paste(as.character(unique(sub_pdata$secondaryaccession)),sep="_",collapse="_"),"_sorted_by_preeclampsia_04.tsv",sep = "")
# filename = paste(paste(as.character(unique(sub_pdata$secondaryaccession)),sep="_",collapse="_"),"_sorted_by_smoking_full_54.tsv",sep = "")
# gofilename

write.table(merged_sel, file.path("/home/sashkoah/a/r/igea-r/article_5/difexp", gofilename), sep = '\t', quote=TRUE)

write.table(merged_sel, file.path("/home/sashkoah/a/r/igea-r/article_5/difexp_go_summary",gofilename), sep = '\t', quote=TRUE)

write.table(sub_exprs, file.path("/home/sashkoah/a/r/igea-r/article_5/exprs_pdata",paste("exprs",filename,sep="_")), sep = '\t', quote=TRUE)
write.table(sub_pdata, file.path("/home/sashkoah/a/r/igea-r/article_5/exprs_pdata",paste("pdata",filename,sep="_")), sep = '\t', quote=TRUE)
filename = NULL
gofilename = NULL

merged_sel




nrow(sub_pdata[sub_pdata$Smoking.Status == "smoker",])




one = read.table("/home/sashkoah/a/r/igea-r/smoking_preeclampsia_diff_exp_with_27272/pe_vs_healthy_more.tsv", sep = '\t')
two = read.table("/home/sashkoah/a/r/igea-r/smoking_preeclampsia_diff_exp_with_27272/smoker_vs_healthy_more.tsv", sep = '\t')
two = merged_sel
intersect_intersect = intersect(intersect_symbol,merged_sel$SYMBOL)
# intersect_intersect
intersect_symbol = intersect(one$SYMBOL, two$SYMBOL)

intersect_entrez = intersect(one$ENTREZID, two$ENTREZID)

one_intersect = one[one$SYMBOL %in% intersect_symbol,]
two_intersect = two[two$SYMBOL %in% intersect_symbol,]
merged = merge(one_intersect,two_intersect, by.x="SYMBOL",by.y="SYMBOL")
merged[,c("logFC.x","logFC.y", "adj.P.Val.x",)]

merged$SYMBOL
write.table(merged, "/home/sashkoah/a/r/igea-r/smoking_preeclampsia_diff_exp_with_27272/merged.tsv", sep = '\t', quote=TRUE)

# sel_full = read.table("/home/sashkoah/a/r/igea-r/smoking_diff_exp/GSE18044.tsv", sep = '\t', quote='"')



sel_full
nrow(sel)
nrow(sel_full)
52-44


setdiff(sel$SYMBOL, sel_full$SYMBOL)


specimen = levels(sub_pdata$Biological.Specimen)[1]

trims = str_replace_all(paste(levels(as.factor(sub_pdata$trim)), collapse = '__'), ' ', '_')
filename = paste(specimen,trims, sep='__')

filename
write.table(sel, file.path('healthy_diff_lists_2019',filename), sep = '\t', quote=TRUE)
nrow(sel)
# read.table(file.path('healthy_diff_lists_2019',filename), sep = '\t', quote='"')

# write.table(diff, "/home/sashkoah/a/r/article-microarrays/diff12.tsv", sep="\t", quote=FALSE)







# BiocManager::install("officer")
library(data.table)

# library(xlsx)
# pd = sub_pdata
# exprs = sub_exprs
# /home/sashkoah/a/r/igea-r/differential_expression_from_data/lumiaffy/trims_healthy
nrow(sub_pdata)
ncol(sub_exprs)
sub_pdata$Biological.Specimen

colnames(sub_exprs) == sub_pdata$Array.Data.File


# unipd$accession

# i=7
# exprs = read.table(paste(mappedpath, '/', studies$accession[[i]], "_mapped_affymetrix.tsv", sep=""), sep="\t", header=TRUE, row.names = 1 )

# write.table(sub_pdata,"sub_pdata.tsv", sep="\t", quote=FALSE)
# 
# sub_pdata = read.table("sub_pdata.tsv", sep="\t", header=TRUE)

propercolnames = as.character(make.names(sub_pdata$Array.Data.File))

colnames(sub_exprs) = as.character(colnames(sub_exprs))
sub_exprs = sub_exprs[,propercolnames]
colnames(sub_exprs) == propercolnames




write.table(sub_exprs, "sub_exprs.tsv", sep="\t", quote=FALSE)
write.table(sub_pdata, "sub_pd.tsv", sep="\t", quote=FALSE)




nrow(sub_pdata[sub_pdata$trim_and_tissue == "First Trimester Chorion",])
nrow(sub_pdata[sub_pdata$trim_and_tissue == "Second Trimester Chorion",])
sub_pdata
trim_and_tissue = as.character(unique(sub_pdata$trim_and_tissue))
trim_and_tissue = trim_and_tissue[c(1,2,4,5,6,7)]
trim_and_tissue

i = 0
j = 0
library(org.Hs.eg.db)
df = NULL
difexp_df = NULL

df = rbind(df, c("Group1", "Group2", "Diff.Exp.Genes", "Centroid.Euclidian.Distances"))
dist_matrix = dist(t(sub_exprs), method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
sub_pdata$trim_and_tissue

length(trim_and_tissue)
for (i in 1:length(trim_and_tissue)){
  for (j in i:length(trim_and_tissue)){
    print(c(i,j)) 
    if (i==j){
      # df = rbind(df, c(trim_and_tissue[i], trim_and_tissue[j], 0))
      next
    }
    # i=1
    # j=2
    sub_sub_pdata = sub_pdata[which(sub_pdata$trim_and_tissue==trim_and_tissue[i] | sub_pdata$trim_and_tissue==trim_and_tissue[j]),]
    
    sub_sub_pdata$Biological.Specimen==sub_sub_pdata[1,]$Biological.Specimen
    
    different_tissue = FALSE %in% (sub_sub_pdata$Biological.Specimen == sub_sub_pdata[1,]$Biological.Specimen)
    different_trim = FALSE %in% (sub_sub_pdata$trim == sub_sub_pdata[1,]$trim)
    different_tissue
    different_trim
    if (different_tissue & different_trim){
      next
    }
    
    sub_sub_exprs = sub_exprs[,rownames(sub_sub_pdata)]
    fit_mod = model.matrix(~ 0 + as.factor(sub_sub_pdata$trim_and_tissue), data=sub_sub_pdata)
    ab = colnames(fit_mod)
    ab
    trim_and_tissue[i]
    trim_and_tissue[j]
    if (grepl(trim_and_tissue[i], colnames(fit_mod)[1])){
      colnames(fit_mod) = c('a','b')
    } else if (grepl(trim_and_tissue[i], colnames(fit_mod)[2])){
      colnames(fit_mod) = c('b','a')
    } else {
      colnames(fit_mod) = c('error')
    }
    contrast.matrix <- makeContrasts(b-a, levels=fit_mod)
    fit <- lmFit(sub_sub_exprs, fit_mod)  # fit each probeset to model
    fit2 <- contrasts.fit(fit, contrast.matrix)
    efit <- eBayes(fit2)        # empirical Bayes adjustment
    t = topTable(efit, number = nrow(sub_sub_exprs))
    t = t[order(-abs(t$logFC)),]
    d = t
    
    d = d[d$adj.P.Val<0.05,]
    d
    diff = d[abs(d$logFC)>1.5,]
    columns(org.Hs.eg.db)
    
    sel = AnnotationDbi::select(org.Hs.eg.db, rownames(diff), c("SYMBOL","GENENAME","GO"))
    library(GO.db)
    
    sel_go = AnnotationDbi::select(GO.db, sel$GO, "TERM")
    sel$Term = sel_go$TERM
    
    columns(org.Hs.eg.db)
    
    sel
    merged_sel = merge(sel,diff, by.x = "ENTREZID", by.y = "row.names")
    
    idx1 = rownames(sub_pdata[sub_pdata$trim_and_tissue==trim_and_tissue[i],])
    idx2 = rownames(sub_pdata[sub_pdata$trim_and_tissue==trim_and_tissue[j],])
    current_dist = dist_between_centroids(dist_matrix, idx1, idx2)
    if (nrow(merged_sel) > 0){
      # source("scripts/article3/gene_description.R")
      library(mygene)
      gene <- getGenes(merged_sel$ENTREZID, fields="summary")
      FALSE %in% (gene$`_id` == merged_sel$ENTREZID)
      merged_sel$summary = gene$summary
      merged_sel$up_down = ifelse(merged_sel$logFC > 0, "up", "down")  
      # print(merged_sel$ENTREZID)
      merged_sel$Group1 = as.character(trim_and_tissue[i])
      merged_sel$Group2 = trim_and_tissue[j]
      merged_sel$Centroid.Euclidian.Distance = current_dist
      merged_sel$Dif.Exp.Amount = nrow(merged_sel)
      difexp_df = rbind(difexp_df, merged_sel)
    }
    
    path = "/home/sashkoah/a/r/igea-r/third_article_difexp_all"
    if (!dir.exists(path)){
      dir.create(path)
    }
    
    ab = as.character(ab)
    b = str_replace_all(strsplit(ab[1], ")", fixed=TRUE)[[1]][2], " ", "_")
    a = str_replace_all(strsplit(ab[2], ")", fixed=TRUE)[[1]][2], " ", "_")
    filename = paste(b, a, "difexp.tsv", sep = '__')
    write.table(merged_sel, file.path(path,filename), sep = '\t', quote=TRUE, row.names = FALSE)
    
    df = rbind(df, c(trim_and_tissue[i], trim_and_tissue[j], nrow(merged_sel), current_dist))
  }
}

nrow(difexp_df[which(difexp_df$Group1 %in% c("First Trimester Placenta", "Third Trimester Placenta") & difexp_df$Group2%in% c("First Trimester Placenta", "Third Trimester Placenta")),])
write.table(difexp_df,"difexp_genes.tsv", sep = '\t', row.names = FALSE)


df
df = as.data.frame(df)
length(levels(df$V1))

new_df <- data.frame(matrix(ncol = length(levels(df$V1)), nrow = length(levels(df$V1))))
colnames(new_df) = as.character(levels(df$V1))
rownames(new_df) = as.character(levels(df$V1))
new_df
df[2,]
for (row in 1:nrow(df)){
  rowname = as.character(df[row,]$V1)
  colname = as.character(df[row,]$V2)
  new_df[rowname,colname] = as.character(df[row,]$V3)
}

new_df

write.table(new_df,"difexp_amounts_new.tsv", sep = '\t', row.names = TRUE)








