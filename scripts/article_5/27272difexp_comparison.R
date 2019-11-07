library(RAM)

from_raw = read.table("/home/sashkoah/a/r/igea-r/article_5/27272_comparison/from_raw.tsv", sep='\t', header = TRUE)
from_processed = read.table("/home/sashkoah/a/r/igea-r/article_5/27272_comparison/from_preprocessed.tsv", sep='\t', header = TRUE)
from_their = read.table("/home/sashkoah/a/r/igea-r/article_5/27272_comparison/their.csv", sep=',', header = TRUE)
from_processed = from_processed[!is.na(from_processed$SYMBOL),]

raw = as.character(from_raw$SYMBOL)
pro = as.character(from_processed$SYMBOL)
their = as.character(from_their$Symbol)


length(raw)
length(pro)
length(their)
NA %in% raw
NA %in% pro
NA %in% their

group.venn(list(raw_155=raw, processed_468=pro, article_329=their), label=FALSE,
           fill = c("red", "blue", "green"),
           # cat.pos = c(0, 0,0),
           lab.cex=30,
           cex=4)