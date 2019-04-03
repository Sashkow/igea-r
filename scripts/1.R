library(sva)
setwd('/home/sashkoah/a/r/igea-r')

# path = 'preprocessed/affymetrix/E-GEOD-14722_preprocessed_affymetrix.tsv'
# pathA = 'preprocessed/affymetrix/ab/E-GEOD-14722_preprocessed_affymetrixA.tsv'
# pathB = 'preprocessed/affymetrix/ab/E-GEOD-14722_preprocessed_affymetrixB.tsv'

# path = 'preprocessed/affymetrix/E-GEOD-14722_preprocessed_affymetrix.tsv'
pathA = 'preprocessed/illumina/1_2/E-GEOD-60438_preprocessed_illumina_1.tsv'
pathB = 'preprocessed/illumina/1_2/E-GEOD-60438_preprocessed_illumina_2.tsv'




# exprs = read.table(path, header = TRUE, sep = '\t')
exprsA = read.table(pathA, header = TRUE, sep = '\t')
exprsB = read.table(pathB, header = TRUE, sep = '\t')

exprsA$group = 'A'
exprsB$group = 'B'

ab = intersect(rownames(exprsA), rownames(exprsB))
ab
nrow(exprsA)
nrow(exprsB)
length(ab)

a.remove = exprsA[!rownames(exprsA) %in% ab, ]
nrow(a.remove)

intersect(rownames(a.remove), rownames(exprsB))

colnames(a.remove) = colnames(exprsB)
a.remove
exprsB
aandb = rbind(a.remove, exprsB)
nrow(aandb)

pdata = data.frame(aandb$group)
rownames(pdata) = rownames(aandb)

aandb <- subset(aandb, select = -c(group))



batch = as.factor(pdata$aandb.group)

aandb.t = t(aandb)

View(aandb)


### Remove batch caused by two different studies
mod = model.matrix(~1, data=pdata)
exprs.nobatch = ComBat(dat=aandb.t, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
exprs.nobatch

write.table(t(exprs.nobatch), path, sep="\t")
t = read.table(path, sep="\t")




exprs = cbind(t(exprs.nobatch),pdata)
View(exprs)
exprs = within(exprs, rm(aandb.group))

library(ggplot2)
library(reshape)

ncol(exprs)
melted = melt(exprs)
melted$group = 'A'
melted[1:nrow(melted)/2,]$group = 'B'



ggplot(melted, aes(x=value)) +
  geom_density(aes(colour=group, group=variable))




indA <- sample.int(12322, size=1232)
indB <- sample.int(9142, size=914)
meltedA <- melt(exprsA)
meltedA$group <- "A"
meltedB <- melt(exprsB)
meltedB$group <- "B"
melted <- rbind(meltedA, meltedB)

meltedA <- melt(exprsA[indA,])
meltedA$group <- "A"
meltedB <- melt(exprsB[indB,])
meltedB$group <- "B"

melted <- rbind(meltedA, meltedB)

ggplot(melted, aes(x=value)) +
  geom_density(aes(group=variable))


# do with annotation libraries what we did with expression matrix
# map propeids to gene entrez ids in exprs

source(paste(getwd(),'scripts/install.brainarray.R',sep='/'))
install.brainarray('hgu133a')
install.brainarray('hgu133b')
require("hgu133ahsentrezg.db")
require("hgu133bhsentrezg.db")


a <- AnnotationDbi::select(hgu133ahsentrezg.db, keys(hgu133ahsentrezg.db), "ENTREZID")
b = AnnotationDbi::select(hgu133bhsentrezg.db, keys(hgu133bhsentrezg.db), "ENTREZID")
ab = intersect(a$PROBEID, b$PROBEID)
justb = b[!(b$PROBEID %in% ab),]
# assert true
nrow(justb) + length(ab) == nrow(b)

aandb = rbind(a, justb)

# remove probes that map onto NA entrez id
aandb = aandb[!is.na(aandb$ENTREZID),]
aandb[is.na(aandb),]


length(unique(aandb$ENTREZID)) == length(aandb$ENTREZID)

newexprs = exprs[rownames(exprs) %in% aandb$PROBEID,]

nrow(newexprs) == nrow(aandb)
# reorder rows in aandb as in newexprs
aandb = aandb[match(rownames(newexprs), aandb$PROBEID),]
rownames(newexprs) = aandb$ENTREZID 

mappedpath = 'mapped/affymetrix'

write.table(newexprs, paste(mappedpath, '/', studies$accession[[i]], "_mapped_affymetrix.tsv", sep=""), sep="\t", quote=FALSE)


