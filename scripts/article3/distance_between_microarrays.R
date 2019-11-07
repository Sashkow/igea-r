# needs sub_pdata and sub_exprs


dist_matrix = dist(t(sub_exprs), method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
# 
# 
# labels(dist_matrix)

trim_and_tissue = as.character(unique(sub_pdata$trim_and_tissue))
i = 0
j = 0

df = NULL
for (i in 1:length(trim_and_tissue)){
  for (j in 1:length(trim_and_tissue)){
    idx1 = rownames(sub_pdata[sub_pdata$trim_and_tissue==trim_and_tissue[i],])
    idx2 = rownames(sub_pdata[sub_pdata$trim_and_tissue==trim_and_tissue[j],])
    current_dist = dist_between_centroids(dist_matrix, idx1, idx2)
    print(c(trim_and_tissue[i],trim_and_tissue[j],current_dist))
    df =rbind(df,c(trim_and_tissue[i],trim_and_tissue[j],current_dist))
    
  }
}
sub_pdata[,c("trim","Gestational.Age", "Average.Gestational.Age")]

df = as.data.frame(df)
df$V3 = as.numeric(as.character(df$V3))
df = df[order(-as.numeric(df$V3)),]
df = df[order(df$V1),]
df
# write.table(df,"euclidian_distances.tsv", sep = '\t', row.names = FALSE)

# idx1 = rownames(sub_pdata[which(sub_pdata$Biological.Specimen=="Chorion" & sub_pdata$trim=="Third Trimester"),])
# idx2= rownames(sub_pdata[which(sub_pdata$Biological.Specimen=="Chorion" & sub_pdata$trim=="First Trimester"),])
# 
# dist_between_centroids(dist_matrix, idx1, idx2)

