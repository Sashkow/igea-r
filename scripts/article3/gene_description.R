# needs merged_sel

# BiocManager::install("mygene")

library(mygene)
gene <- getGenes(merged_sel$ENTREZID, fields="summary")
FALSE %in% (gene$`_id` == merged_sel$ENTREZID)
merged_sel$summary = gene$summary
