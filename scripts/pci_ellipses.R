source("https://bioconductor.org/biocLite.R")
biocLite("StatPerMeCo")
library(StatPerMeCo)

library(rlang)
library(ggplot2)
library("factoextra")
X = matrix(rnorm(4000),ncol=4)

biocLite("tidyverse")
library(tidyverse)