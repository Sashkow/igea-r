data = read.table("/home/sashkoah/a/r/igea-r/pure_mixed/igea-cibersort-python/cell_types/re.tsv", sep = '\t', header = TRUE)
data = data[which(data$Biological.Specimen=='Chorion' | data$Biological.Specimen=='Decidua'),]
data$Biological.Specimen
data$is_decidua = ifelse(data$Biological.Specimen=='Decidua',2,3)
data$dendric4
model <- lm(is_decidua ~ evt5, data=data)
summary(model)

data$stb7



