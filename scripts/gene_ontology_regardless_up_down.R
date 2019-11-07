# BiocManager::install("openxlsx")

install.packages("devtools") 
devtools::install_github("wjawaid/enrichR")

library(readxl)
# library(xlsx)
library(enrichR)
library(officer)
library(openxlsx)


# library(dplyr)
# library(biomaRt)
# library(karyoploteR) 
# library(regioneR)
# library(EnhancedVolcano)
# library(TissueEnrich)
# library(tidyr)
# library(pheatmap)
# library(org.Hs.eg.db)
# library(cowplot)

difexppath = "/home/sashkoah/a/r/igea-r/article_3/placenta/term vs second trimester logfc 1"
gopath = "/home/sashkoah/a/r/igea-r/article_3/placenta/term vs second trimester logfc 1/go"

difexpfiles = list.files(difexppath, full.names = TRUE)
difexpfiles
dbs <- c("GO_Biological_Process_2018")

difexpfiles[1]
for (difexpfile in difexpfiles){
  difexpfile = difexpfiles[1]
  difexpfile_no_path = tail(strsplit(difexpfile,'/')[[1]], n=1)
  # difexp = read_excel(difexpfile, sheet = 1)
  difexp = read.table(difexpfile,sep="\t", header = TRUE)
  genes = unique(difexp$SYMBOL)

  # genes = rownames(difexp)
  genes
  gofilename = difexpfile_no_path
  print(gofilename)
  genes
  enriched_full = enrichr(as.character(genes), dbs)
  enriched = enriched_full$GO_Biological_Process_2018
  enriched = enriched[,c('Term', 'Overlap', 'Adjusted.P.value', 'Odds.Ratio','Genes')]
  enriched$Adjusted.P.value
  enriched = enriched[which(as.numeric(enriched$Adjusted.P.value)<=0.05),]
  write.xlsx(enriched, file.path(gopath,gofilename), sheetName = "BiologicalProcess", 
             col.names = TRUE, row.names = TRUE, append = FALSE)
  
}

genes = read.table(path, sep='\t', header = TRUE)

# in 2 vs 1 trim
upgenes = genes[genes$up_down=="up",]
upgenes = unique(upgenes$SYMBOL)
upgenes = as.character(upgenes)

# dbs <- listEnrichrDbs()
dbs
geneontology_files = list.files(difexppath)
GO_Biological_Process_2018


enriched <- enrichr(upgenes, dbs)
enriched = enriched$GO_Biological_Process_2018[order(enriched$GO_Biological_Process_2018$P.value),]
enriched
enriched[which(enriched$P.value < 0.05),]
