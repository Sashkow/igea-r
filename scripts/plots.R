# biocLite("cowplot")
require(ggplot2)
require(cowplot)
library(cowplot)
require(ggfortify)
library(grid)
library(reshape2)
library(RColorBrewer)

pcaplotspath = 'plots/pca'

getAspectRatio <- function(p){
  gb <- ggplot_build(p)
  g <- ggplot_gtable(gb)
  +
  nullw <- sapply(g$widths, attr, "unit")
  nullh <- sapply(g$heights, attr, "unit")
  
  # ar of plot
  if(any(nullw == "null"))
    ar <- unlist(g$widths[nullw == "null"]) / unlist(g$heights[nullh == "null"])
  
  # ar of plot + legend
  g$fullwidth <- convertWidth(sum(g$widths), "in", valueOnly=TRUE)
  g$fullheight <- convertHeight(sum(g$heights), "in", valueOnly=TRUE)
  ar <- g$fullwidth / g$fullheight
  
  return(ar)
}



# pca.data = pca
# pheno.data = sub_pdata
# meta.vars = c('Diagnosis', 'Gestational.Age.Category', 'accession', 'Biological.Specimen', 'Fetus.Sex', 'is_outlier')
# ncol = ncol

pcaPlots <- function(pca.data, pheno.data, meta.vars, title, ncol) {
  pheno.data[] <- lapply(pheno.data, as.character)
  plots <- c()
  ar <- -100
  for (i in meta.vars) {
    pl <- fviz_pca_ind(pca.data, 
                         label=c("none"),
                         habillage=pheno.data[,i],
                         # palette="Dark2",
                         repel = TRUE,
                         title = i# Avoid text overlappingi
                         # select.va = list(name=c("162466", "212", "2993", "3577", "3579", "407008", "439996", "4778", "51327", "6283", "6402", "6521", "6563", "759", "84658", "8875"))
                         # col.var = "contrib",
    )
    filename = paste(c(title,"_",i,".png"), sep = "")
    
    print("here")
    print(filename)
    print(file.path(pcaplotspath, filename))

    # save_plot(file.path(pcaplotspath,filename),
    #           pl
    # )
    # pl <- autoplot(pca.data, data = pheno.data, colour=i) +
    #   coord_fixed()
    newar <- getAspectRatio(pl)
    if (ar<newar) {
      ar <- newar
    }
    plots <- c(plots, list(pl))
  }
  
  if(missing(ncol)) {
    pl <- plot_grid(plotlist = plots, ncol=length(meta.vars), align="hv")
  } else {
    pl <- plot_grid(plotlist = plots, ncol=ncol, align="hv")
  }
  if (!missing(title)) {
    title <- ggdraw() + draw_label(title, fontface='bold')  
    pl <- plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))
  }
  
  
  
  pl <- pl + theme(plot.margin=margin(t=10, r=10, b=10, l=10))
  return(list(pl, ar))
}

