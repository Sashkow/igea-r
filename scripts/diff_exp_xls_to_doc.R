# BiocManager::install("GOfuncR")
library(GOfuncR)
library(officer)
library(magrittr)
library(openxlsx)
library(readxl)
my_doc <- read_docx() 
styles_info(my_doc)

# src <- tempfile(fileext = ".png")
# png(filename = src, width = 5, height = 6, units = 'in', res = 300)
# barplot(1:10, col = 1:10)
difexpfolder = "/home/sashkoah/a/r/igea-r/third_article_difexp_xlsx/placenta"
# dir.create("/home/sashkoah/a/r/igea-r/third_article_difexp_docx/placenta", recursive = TRUE)
docgopath = "/home/sashkoah/a/r/igea-r/third_article_difexp_docx/placenta"
difexplist = list.files(difexpfolder)


for (difexpfile in difexplist){
  # difexpfile = difexplist[1]
  my_doc <- read_docx()
  docfilename = paste(tail(strsplit(difexpfile,"/")[[1]],1) ,"docx", sep=".")
  difexppath = file.path(difexpfolder, difexpfile)
  difexp = read_excel(difexppath, sheet = 1)
  my_doc = my_doc %>% body_add_par(difexpfile, style = "Normal")
  my_doc = my_doc %>% body_add_par(paste("Genes total: ",difexp[1,]$Dif.Exp.Amount, sep=""), style = "Normal")
  difexp_up = difexp[which(difexp$up_down == "up"),]
  difexp_up = difexp_up[order(difexp_up$SYMBOL),]
  genes = difexp_up$SYMBOL
  genes
  my_doc = my_doc %>% body_add_par("Upregulated:", style = "heading 1")
  for (gene in genes){
    # gene = genes[5]
    my_doc = my_doc %>% body_add_par(gene, style = "heading 2")
    difexp_gene = difexp_up[which(as.character(difexp_up$SYMBOL)==gene),][1,]
    difexp_gene
    gene_description = paste(
      difexp_gene$GENENAME,
      difexp_gene$ENTREZID,
      paste("logfc", round(as.numeric(difexp_gene$logFC),3), sep=" = "),
      paste("adj.P.Val", round(difexp_gene$adj.P.Val,20), sep = " = "),
      difexp_gene$summary,
      sep = "; "
    )
    my_doc = my_doc %>% body_add_par(gene_description, style = "Normal")
  }
  #down
  my_doc = my_doc %>% body_add_par("Downregulated:", style = "heading 1")
  difexp_up = difexp[which(difexp$up_down == "down"),]
  difexp_up = difexp_up[order(difexp_up$SYMBOL),]
  genes = difexp_up$SYMBOL
  
  for (gene in genes){
    # gene = genes[1]
    my_doc = my_doc %>% body_add_par(gene, style = "heading 2")
    difexp_gene = difexp_up[which(as.character(difexp_up$SYMBOL)==gene),][1,]
    gene_description = paste(
      difexp_gene$GENENAME,
      difexp_gene$ENTREZID,
      paste("logfc", round(as.numeric(difexp_gene$logFC),3), sep=" = "),
      paste("adj.P.Val", round(difexp_gene$adj.P.Val,20), sep = " = "),
      difexp_gene$summary,
      sep = "; "
    )
    my_doc = my_doc %>% body_add_par(gene_description, style = "Normal")
  }
  
  print(docfilename)
  # dir.create(docgopath, recursive = TRUE)
  print(my_doc, target = file.path(docgopath,docfilename))
}




difexp = read.table(path,sep = '\t', header = TRUE)
term_ids = as.character(levels(difexp$GO))
length(term_ids)
difexp$GO = as.character(difexp$GO)
for (term_id in term_ids){
  sub_difexp = difexp[difexp$GO==term_id,]
  print(c(term_id,nrow(sub_difexp), sub_difexp$SYMBOL))
}


term_ids[1]

difexp[difexp$GO == "GO:0000015",]

difexp[is.na(difexp$SYMBOL),]
sub_difexp[is.na(sub_difexp$SYMBOL),]
sub_difexp

dev.off()

my_doc <- my_doc %>% 
  body_add_par("Hello world!", style = "heading 1") %>% 
  body_add_par("this", style = "Normal") %>%
  body_add_par("Hello world!", style = "heading 1") %>% 
  body_add_par("this", style = "Normal") %>% # blank paragraph
  
my_doc

dir.create('assets/docx', recursive = TRUE)
print(my_doc, target = "assets/docx/first_example.docx")


# go to docx

difexpfolder= "/home/sashkoah/a/r/igea-r/third_article_difexp_xlsx"
gofolder  = "/home/sashkoah/a/r/igea-r/go_article_3/go_article_3_v_2/chorion"
docgopath = "/home/sashkoah/a/r/igea-r/go_article_3/docx/chorion"
dir.create(docgopath, recursive = TRUE)

gopathlist = list.files(gofolder)
gopathlist
# list.files("/home/sashkoah/a/r/igea-r/third_article_difexp_xlsx/")

difexplist = c("First_Trimester_Chorion__Second_Trimester_Chorion__difexp.tsv.xlsx",
               "First_Trimester_Chorion__Third_Trimester_Chorion__difexp.tsv.xlsx",        
               "Second_Trimester_Chorion__Third_Trimester_Chorion__difexp.tsv.xlsx",       
               "First_Trimester_Chorion__Second_Trimester_Chorion__difexp.tsv.xlsx",  
               "First_Trimester_Chorion__Third_Trimester_Chorion__difexp.tsv.xlsx",          
               "Second_Trimester_Chorion__Third_Trimester_Chorion__difexp.tsv.xlsx")

for (i in 1:length(gopathlist)){
  gopath = file.path(gofolder,gopathlist[i])
  difexppath = file.path(difexpfolder,difexplist[i])
  
  docfilename = paste(tail(strsplit(gopath,"/")[[1]],1) ,"docx", sep=".")
  docfilename
  my_doc <- read_docx()
  styles_info(my_doc)


  difexp = read_excel(difexppath, sheet = 1)
  go = read_excel(gopath, sheet = 1)
  if (nrow(go) == 0){
    next
  }
  # strsplit(go$Genes,";")
  
  my_doc = my_doc %>% body_add_par(docfilename, style = "Normal")
  my_doc = my_doc %>% body_add_par(paste("Genes total: ",difexp[1,]$Dif.Exp.Amount, sep=""), style = "Normal")

  for (term_id in 1:nrow(go)){
    # print(term)
    # term = go[1,]
    term = go[term_id,]
    term_description = paste(
      term$Term,
      paste("overlap",term$Overlap,sep=" = "),
      paste("adj.pvalue",round(term$Adjusted.P.value,5),sep=" = "),
      paste("odds ratio",round(term$Odds.Ratio,2),sep=" = "),
      sep = "; "
    )
    my_doc = my_doc %>% body_add_par(term_description, style = "heading 1")
    genes = strsplit(term$Genes, ";")[[1]]
    genes
    for (gene in genes){
      # gene = genes[1]
      my_doc = my_doc %>% body_add_par(gene, style = "heading 2")
      difexp_gene = difexp[which(difexp$SYMBOL==gene),][1,]
      gene_description = paste(
        difexp_gene$GENENAME,
        difexp_gene$ENTREZID,
        paste("logfc", round(as.numeric(difexp_gene$logFC),3), sep=" = "),
        paste("adj.P.Val", round(difexp_gene$adj.P.Val,10), sep = " = "),
        difexp_gene$summary,
        sep = "; "
      )
      my_doc = my_doc %>% body_add_par(gene_description, style = "Normal")
    }
  }
  print(docfilename)
  dir.create(docgopath, recursive = TRUE)
  print(my_doc, target = file.path(docgopath,docfilename))
}


# short version 


for (i in 1:length(gopathlist)){
  gopath = file.path(gofolder,gopathlist[i])
  difexppath = file.path(difexpfolder,difexplist[i])
  
  docfilename = paste(tail(strsplit(gopath,"/")[[1]],1) ,"docx", sep=".")
  docfilename
  my_doc <- read_docx()
  styles_info(my_doc)
  
  
  difexp = read_excel(difexppath, sheet = 1)
  go = read_excel(gopath, sheet = 1)
  if (nrow(go) == 0){
    next
  }
  # strsplit(go$Genes,";")
  
  my_doc = my_doc %>% body_add_par(docfilename, style = "Normal")
  my_doc = my_doc %>% body_add_par(paste("Genes total: ",difexp[1,]$Dif.Exp.Amount, sep=""), style = "Normal")
  
  for (term_id in 1:nrow(go)){
    # print(term)
    # term = go[1,]
    term = go[term_id,]
    term_description = paste(
      term$Term,
      paste("overlap",term$Overlap,sep=" = "),
      paste("adj.pvalue",round(term$Adjusted.P.value,10),sep=" = "),
      paste("odds ratio",round(term$Odds.Ratio,2),sep=" = "),
      sep = "; "
    )
    my_doc = my_doc %>% body_add_par(term_description, style = "heading 1")
    genes = strsplit(term$Genes, ";")[[1]]
    genes
    for (gene in genes){
      # gene = genes[1]
      # my_doc = my_doc %>% body_add_par(gene, style = "heading 2")
      difexp_gene = difexp[which(difexp$SYMBOL==gene),][1,]
      gene_description = paste(
        gene,
        difexp_gene$GENENAME,
        difexp_gene$ENTREZID,
        paste("logfc", round(as.numeric(difexp_gene$logFC),3), sep=" = "),
        paste("adj.P.Val", round(difexp_gene$adj.P.Val,10), sep = " = "),
        sep = "; "
      )
      my_doc = my_doc %>% body_add_par(gene_description, style = "heading 2")
    }
  }
  print(docfilename)
  
  print(my_doc, target = file.path(docgopath,docfilename))
}





