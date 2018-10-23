setwd('/home/sashkoah/a/r/article-microarrays')
#' Merges all files ending with "_sample_table.txt" into one tsv file 
mergeLumiTables <- function(path, files=NULL){
  files = list.files(path, pattern = "*sample_table.txt")
  # get column with identifiers from the first sample_table
  f = read.csv(paste(path, files[1], sep='/'),sep = '\t') 
  # merged = file.create(paste(path, "merged_sample_table.csv", sep='/'))
  merged = data.frame(matrix(NA,nrow=nrow(f[1]),ncol=0))
  
  for (file in files){
    # if samplple_table.txt in file name
    if (grepl("sample_table.txt", file)){
      sample = read.csv(paste(path,file, sep='/'), sep='\t')
      colnames(sample)[2] = file
      merged = cbind(merged, sample[2])
    }
    print(paste("processed",file))
  }
  rownames(merged) = f[[1]]
  write.csv(merged, paste(path, "merged_sample_table.csv", sep='/'))
  return(merged)
}


# answ = mergeLumiTables('raws/illumina/E-GEOD-35574')
# answ
# f = read.csv('raws/illumina/E-GEOD-35574/GSM871046_sample_table.txt',sep = '\t')
# f2 = read.csv('raws/illumina/E-GEOD-35574/GSM871047_sample_table.txt',sep = '\t')
# 
# cbind(f[1:2],f2[2])
# class(f)
# 
# View(f)