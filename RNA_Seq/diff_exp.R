##Using deseq package
##Using exp_prof and other variables from Feature_Extract.R

library(BiocStyle)
library(rmarkdown)
library(geneplotter)
library(ggplot2)
library(plyr)
library(LSD)
library(gplots)
library(RColorBrewer)
library(stringr)
library(topGO)
library(genefilter)
library(biomaRt)
library(dplyr)
library(EDASeq)
library(fdrtool)
library(DESeq2)
sample_info = data.frame(type = group.samples[diff.ids], names = labels[diff.ids])
rownames(sample_info) = sample.ids[diff.ids]
dds <- DESeqDataSetFromMatrix(countData = exp_prof[,diff.ids], colData = sample_info, design = ~0+type)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
multidensity( counts(dds, normalized = T),
              xlab="mean counts", xlim=c(0, 1000))
multiecdf( counts(dds, normalized = T),
           xlab="mean counts", xlim=c(0, 1000))



dds <- DESeq(dds, betaPrior = F)
res <- results(dds)
#library("BiocParallel")
#register(MulticoreParam(4))

#cols <- brewer.pal(8, "Set1")
#boxplot(exp_prof)

rld <- rlogTransformation(dds)

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              summary(res)
View(res)
res <- res[order(res$padj), ]
summary(res)
diff.expr <- rownames(res[res$padj < 0.01 & abs(res$log2FoldChange) > 2, ])

##Count Plot
d <- plotCounts(dds, gene=which.min(res$log2FoldChange), intgroup="type",
                returnData=TRUE)
ggplot(d, aes(x=type, y=count)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(25,100,400))

plotMA(res, ylim = c(-2,2))

select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
heatmap(assay(dds), Rowv = NA, Colv = NA, col = cm.colors(256), scale="column", margins=c(5,10))

d <- dist(t(counts(dds,normalized = T)[diff.expr,]))
plot(hclust(d))

library(FactoMineR)
library(factoextra)
a = t(counts(dds, normalized = T)[diff.expr, ])
rownames(a) = sample_info$type
pc = PCA(a)
fviz_pca_ind(pc)
###Using edgeR
library(edgeR)
y <- DGEList(counts=exp_prof, group=group.samples)
keep <- rowSums(cpm(y)>1) >= 2; y <- y[keep, ]
y <- calcNormFactors(y)
design <- model.matrix(~0+group, data=y$samples); colnames(design) <- levels(y$samples$group)
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
edgeglm <- as.data.frame(topTags(lrt, n=length(rownames(y))))
