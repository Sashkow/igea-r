#healthy chorion
a = read.table("/home/sashkoah/a/r/igea-r/healthy_diff_lists_2019/Placenta__First_Trimester__Third_Trimester", sep = '\t', quote='"')

b = read.table("/home/sashkoah/a/r/igea-r/differential_expression_from_literature/GSE9984/1_3_trim_healthy_placenta.csv", sep=",", quote='"')

nrow(a)
nrow(b)
a = as.character(a$V1)
b = as.character(b$V1)

group.venn(list(Literature=b,  Data=a), label=FALSE,
           fill = c("orange", "blue"),
           cat.pos = c(0, 0),
           cat.cex=5,
           cex=5)

nrow(pdata[which(pdata$Biological.Specimen=="Chorion" & pdata$Diagnosis=="Healthy"),])

nrow(pdata[which(pdata$Biological.Specimen=="Decidua" & pdata$Diagnosis=="Healthy"  & pdata$accession=="E-GEOD-60438"),])
nrow(pdata[which(pdata$Biological.Specimen=="Decidua" & pdata$Diagnosis=="Healthy"  & pdata$accession=="E-GEOD-60438"),])

95-65-18
55-18-12

12


nrow(pdata)

#literature global
a = read.table("/home/sashkoah/a/r/article-microarrays/differential_expression_from_literature/GSE73374/GSE73374_literature.csv", header = TRUE, quote = '"', sep = ",")

#literature local

b = read.table("/home/sashkoah/a/r/article-microarrays/differential_expression_from_literature/GSE73374/pone.0141294.s006-1.csv", header = TRUE, quote = '"', sep = ",")

# b = read.table("/home/sashkoah/a/r/article-microarrays/differential_expression_from_data/3trim_placenta/3trim_placenta.tsv",header = TRUE, sep = "\t")

#merged
c = read.table("/home/sashkoah/a/r/article-microarrays/differential_expression_from_data/3trim_placenta/3trim_placenta_symbol.tsv", header = TRUE, sep = "\t")



a = read.table("/home/sashkoah/a/r/igea-r/73685_60438_chor_dec/diff_exp_full", header = TRUE, quote = '"', sep = "\t")
b = read.table("/home/sashkoah/a/r/igea-r/73685_60438_chor_dec/diff_exp_50", header = TRUE, quote = '"', sep = "\t")
c = read.table("/home/sashkoah/a/r/igea-r/73685_60438_chor_dec/diff_exp_75", header = TRUE, quote = '"', sep = "\t")
nrow(a)
nrow(b)
nrow(c)
# my gse
d = read.table()

global_literature = as.character(a$Gene.Symbol)
merged = as.character(c$SYMBOL)
local_literature = as.character(b[b$Gene.Name!="---",]$Gene.Name)
local_literature = unique(local_literature)
length(local_literature)



14+24+22
24+22+6
22+6+16





# print(Sys.setenv(R_MAX_NUM_DLLS = "10000"))  # `A+C` could also be used
# Sys.getenv("R_MAX_NUM_DLLS")

library(RAM)

group.venn(list(Literature=global_literature,GSE73374=local_literature,  merged_data=merged), label=FALSE,
           fill = c("orange", "blue")
           # cat.pos = c(0, 0,0),
           lab.cex=36,
           cex=5)

group.venn(list(full=a$SYMBOL, fifty=b$SYMBOL, seventy_five=c$SYMBOL), label=FALSE,
           fill = c("orange", "blue", "green"),
           # cat.pos = c(0, 0,0),
           lab.cex=36,
           cex=5)


# group.venn(list(chorion_genes=chorion_genes,
#                 decidua_genes=decidua_genes,
#                 placenta_genes=placenta_genes),
#            label=TRUE,
#            fill = c("orange", "blue"),
#            lab.cex=1.1)


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
