BiocManager::install("ALL")
library(topGO)
library(ALL)
data("ALL")
affyLib <- paste(annotation(ALL), "db", sep = ".")
# BiocManager::install(affyLib)
data(geneList)
library(package = affyLib, character.only = TRUE)

class(geneList)


sampleGOdata <- new("topGOdata",
                    description = "Simple session", ontology = "BP",
                    allGenes = geneList, geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.db, affyLib = affyLib)
sampleGOdata
