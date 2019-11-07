their = read.table("/home/sashkoah/a/r/igea-r/article_5/27272_comparison/their.csv", sep = ",", header = TRUE)
raw = read.table("/home/sashkoah/a/r/igea-r/article_5/27272_comparison/from_raw.tsv", sep = "\t",header = TRUE)
proccessed = read.table("/home/sashkoah/a/r/igea-r/article_5/27272_comparison/from_preprocessed.tsv", sep = "\t",header = TRUE)
raw
nrow(their)
nrow(raw)
nrow(proccessed)
proccessed$SYMBOL

length(intersect(as.character(their$Symbol),as.character(rownames(raw))))
length(intersect(as.character(their$Symbol),as.character(proccessed$SYMBOL)))
length(intersect(as.character(rownames(raw)),as.character(proccessed$SYMBOL)))


intersect(as.character(their$Symbol),as.character(rownames()))

library(RAM)

raw = raw[which(raw$P.Value < .05),]


foo <- as.character(their$Symbol)
baa <- as.character(rownames(raw))
dee = as.character(proccessed$SYMBOL)
dee = dee[!is.na(dee)]


group.venn(list(their_329=foo, raw_173=baa, processed_468=dee), label = FALSE,
           fill = c("orange", "blue", "red"),
           # cat.pos = c(0, 0,0),
           lab.cex=1.1)



s1 = sample(1:20000, 329, replace=TRUE)
s2 = sample(1:20000, 173, replace=TRUE)
s3 = sample(1:20000, 468, replace=TRUE)

group.venn(list(their_329=s1, raw_173=s2, processed_468=s3), label = FALSE,
           fill = c("orange", "blue", "red"),
           # cat.pos = c(0, 0,0),
           lab.cex=1.1)

          