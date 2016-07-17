install.packages(c("samr", "matrixStats", "GSA", "shiny", "shinyFiles", "openxlsx"))
source("http://bioconductor.org/biocLite.R")
biocLite("impute")
library(shiny)
library(shinyFiles)

runGitHub("SAM", "MikeJSeo")


View(eset)

library(hgu133a.db)
library(AnnotationDbi)
probe.ids <- rownames(tT)
inp.sam.df <- as.matrix(eset)
inp.sam.df <- select(hgu133a.db, probe.ids, c('SYMBOL'))
sample.class = c(2,1,2,1,1,2,1,2,1,1,2,1,2,1,2,1,2)
colnames(inp.sam.df) = sample.class
sum(is.na(match(rownames(inp.sam.df), rownames(tT.1) ))) == 0
gene.ids = as.character(tT.1$Gene.Symbol[match(rownames(inp.sam.df), rownames(tT.1) )])
length(gene.ids) == length(rownames(tT))
write.csv(inp.sam.df, 'inp.sam.csv')
