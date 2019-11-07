# BiocManager::install("topGO")
library(enrichR)
library(openxlsx)
library(topGO)

logfc = 1

condition_up = function(expression_values){
  return(expression_values > logfc)
}

condition_down = function(expression_values){
  return(expression_values < -logfc)
}

condition_const = function(expression_values){
  return(expression_values >= -logfc & expression_values <= logfc)
}

condition_any = function(expression_values){
  return(expression_values > logfc | expression_values <= logfc)
}

kind_to_condition = function(kind){
  if (kind == "Up"){
    return(condition_up)
  } else if (kind == "Down"){
    return(condition_down)
  } else if (kind == "Const"){
    return(condition_const)
  } else if (kind == "Any"){
    return(condition_any)
  } else{
    print("Unacceptable Condition")
    return(NULL)
  }
}

get_gene_condition = function(difexp, column, kind){
  gene_condition = kind_to_condition(kind)
  return(which(gene_condition(difexp[,column])))
}

get_genes = function (difexp, column, kind){
  return(difexp[get_gene_condition(difexp,column,kind),])
}


difexp = read.table("/home/sashkoah/a/r/igea-r/article_3/placenta/multiple_testing_methods intersect/difexp_term_2_and_2_1_full_with_averages_multiple_testing_methods_intersect.tsv", sep = '\t', header = TRUE)

#test1 
column = "ii...i"
kind = "Up"
!(FALSE %in% (get_genes(difexp,column,kind)$ii...i > logfc))

#test1 
column = "ii...i"
kind = "Down"
!(FALSE %in% (get_genes(difexp,column,kind)$ii...i < logfc))

#test1 
column = "ii...i"
kind = "Const"
!(FALSE %in% (get_genes(difexp,column,kind)$ii...i >= -logfc))
!(FALSE %in% (get_genes(difexp,column,kind)$ii...i <= logfc))

#test2
column = "ii...i"
kind = "Any"
!(FALSE %in% (get_genes(difexp,column,kind)[,column] == difexp[,column]))


kinds = c("Up", "Down", "Const")
kinds = c("Up", "Down", "Const")
columns = c("ii...i", "iii...ii")

wb = createWorkbook("difexp_fil")
addWorksheet(wb, "Data")
writeData(wb,"Data", difexp)

i=0
for (kind1 in kinds){
  for (kind2 in kinds){
    i=i+1
    sub_difexp = get_genes(difexp,columns[1],kind1)
    sub_difexp = get_genes(sub_difexp,columns[2],kind2)
    worksheet_name = paste(kind1, kind2, sep="_")
    addWorksheet(wb, worksheet_name)
    writeData(wb,worksheet_name,sub_difexp)
    print(c(i,kind1, kind2, nrow(sub_difexp), round(100*nrow(sub_difexp)/nrow(difexp),2) ))
  }
}

saveWorkbook(wb,"/home/sashkoah/a/r/igea-r/article_3/placenta/multiple_testing_methods intersect/difexp_term_2_and_2_1_full_with_averages_multiple_testing_methods_intersect_multisheet.xlsx")
wb = loadWorkbook("/home/sashkoah/a/r/igea-r/article_3/placenta/multiple_testing_methods intersect/difexp_term_2_and_2_1_full_with_averages_multiple_testing_methods_intersect.tsv")
gowb = createWorkbook("go_difexp")

gopath = "/home/sashkoah/a/r/igea-r/article_3/placenta/multiple_testing_methods intersect/go"

dbs <- c("GO_Biological_Process_2018")
wb$sheet_names

for (sheet_name in wb$sheet_names){
  difexp<-read.xlsx(wb, sheet=sheet_name, colNames=T)
  
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
  enriched = enriched[which(as.numeric(enriched$Adjusted.P.value)<=0.1),]
  
  addWorksheet(gowb, sheet_name)
  writeData(gowb,sheet_name,enriched)
  
  
}
saveWorkbook(gowb, file.path(gopath,"go.xlsx"))


for (sheet in wb$tables){
  print(nrow(sheet))
}













