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

diff.genes = list()
##Match-Control Samples
match.control = list()
match.control[['sample.names']] = sample.ids[diff.ids]
match.control[['type']] = group.samples[diff.ids]
match.control[['pat.ids']] = pat.ids[diff.ids]
match.control[['pat.short.ids']] = c()
match.control[['pat.comb']] = mapply(function(x,y){
                  paste(as.character(x),y,sep = '')},
                  match.control$pat.short.ids, match.control$type)
exp_prof_diff = exp_prof[,diff.ids]
match.control$pat.comb[85] = paste(match.control$pat.comb[85], 'B', sep='_')
match.control$pat.comb[84] = paste(match.control$pat.comb[84], 'A', sep='_')

match.control$pat.comb[97] = paste(match.control$pat.comb[97], 'A', sep='_')
match.control$pat.comb[98] = paste(match.control$pat.comb[98], 'B', sep='_')

colnames(exp_prof_diff) = match.control$pat.comb
j = 1
pat = match.control[['pat.ids']][1]
for(i in seq_along(diff.ids))
{
    if(match.control$pat.ids[i] != pat)
    {
      j = j + 1
      pat = match.control$pat.ids[i]
    }
  match.control$pat.short.ids[i] = j
}

sample_info = data.frame(match.control)
rownames(sample_info) = colnames(exp_prof_diff)

dds <- DESeqDataSetFromMatrix(countData = exp_prof_diff, colData = sample_info, design = ~pat.short.ids + type)
dds_sim <- DESeqDataSetFromMatrix(countData = exp_prof_diff, colData = sample_info, design = ~type)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
multidensity( counts(dds, normalized = T),
              xlab="mean counts", xlim=c(0, 1000))
multiecdf( counts(dds, normalized = T),
           xlab="mean counts", xlim=c(0, 1000))


                                                                                                                                                                                                                                                                                                                                                                                                                                                                        library('IHW')
dds <- DESeq(dds)
dds_sim <- DESeq(dds_sim)
res <- results(dds)
res_sim <- results(dds_sim)

#library("BiocParallel")
#register(MulticoreParam(4))

#cols <- brewer.pal(8, "Set1")
#boxplot(exp_prof)

rld <- rlogTransformation(dds)
rld_sim <- rlogTransformation(dds_sim)
summary(res)
View(res)
res <- res[order(res$padj), ]
summary(res)
diff.genes[['deseq']] <- rownames(res[res$padj < 0.01 & abs(res$log2FoldChange) > 2.5, ])

##Count Plot
d <- plotCounts(dds, gene=which.min(res$log2FoldChange), intgroup="type",
                returnData=TRUE)
ggplot(d, aes(x=type, y=count)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(25,100,400))

plotMA(res, ylim = c(-2,2))
nt <- normTransform(dds)
vs <- vst(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
###heatmap(assay(dds), Rowv = NA, Colv = NA, col = cm.colors(256), scale="column", margins=c(5,10))
library(FactoMineR)
library(factoextra)
library('vsn')

##VST
d <- dist(t(assay(vs)[diff.genes$deseq,]))
plot(hclust(d), labels = sample_info$type)
plotPCA(vs, intgroup = 'type')

######RLD
#Hie Clus
d <- dist(t(assay(rld)[diff.genes$deseq,]))
plot(hclust(d), labels = sample_info$type)
plotPCA(rld, intgroupzz)
#PCA
df <-t(assay(rld)[diff.genes$deseq,])
rownames(df) = sample_info$type
pc <- PCA(df)
fviz_pca_ind(pc)

hcc <-HCPC(pc)


#MDS Plot
fit <- cmdscale(d,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric	MDS",	type="n")
text(x, y, labels = sample_info$type, cex=.7)
library(MASS)
fit <- isoMDS(d, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric	MDS",	type="n")
text(x, y, labels = rownames(sample_info), cex=.7)

###Normalised
#Hie Clus
d <- dist(t(assay(nt)[diff.genes$deseq,]))
plot(hclust(d), labels = sample_info$type)

#PCA
df = t(assay(nt)[diff.genes$deseq,])
rownames(df) = sample_info$type
pc <- PCA(df)
fviz_pca_ind(pc)
#MDS Plot
fit <- cmdscale(d,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric	MDS",	type="n")
text(x, y, labels = sample_info$type, cex=.7)

library(pheatmap)
pheatmap(assay(nt)[diff.genes$deseq,], rowna)

library(GOstats)
library(GO.db)
library("AnnotationForge")
library("org.Hs.eg.db")
diff.genes$exact = read.delim('edgreR.txt', header = F)
diff.genes$exact = diff.genes$exact$V1
diff.genes.mod = list()
diff.genes.mod[['deseq']] = sapply(diff.genes$deseq, function(x) 
{
  unlist(strsplit(x, split = '.', fixed = T))[1]
})
diff.genes.mod[['edgeR']] = genes.entrez$entrezgene[match(diff.genes$exact,as.character(genes.entrez$ensembl_gene_id))]
universal.entrez = unique(genes.entrez$entrezgene[match(g,genes.entrez$ensembl_gene_id)])
ent.ids = genes.entrez$entrezgene[match(diff.genes.mod$deseq,as.character(genes.entrez$ensembl_gene_id))]
write(as.character(ent.ids), 'deseq.txt')
write(as.character(diff.genes.mod$edgeR), 'edgeR.txt')
params <- new('GOHyperGParams', geneIds = ent.ids, universeGeneIds = universal.entrez, 
              ontology = 'BP', pvalueCutoff = 0.5, annotation = 'org.Hs.eg.db')
hgOver <- hyperGTest(params)

params1 <- new('GOHyperGParams', geneIds = diff.genes.mod$edgeR, universeGeneIds = universal.entrez, 
              ontology = 'BP', pvalueCutoff = 0.1, annotation = 'org.Hs.eg.db')
hgOver1 <- hyperGTest(params1)

df <- summary(hgOver)
df1 <- summary(hgOver1)
write(as.character(universal.entrez), 'all_genes.txt')

common = intersect(diff.genes.mod$edgeR, ent.ids)

only_cancer = data.frame(type = group.samples[tumor.indexes])

dds_full = DESeqDataSetFromMatrix(exp_prof[,tumor.indexes], only_cancer, design = ~type)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nt_full = normTransform(dds_full)
d = dist(t(assay(nt_full)[diff.genes$deseq, ]))
plot(hclust(d), labels = group.samples)
