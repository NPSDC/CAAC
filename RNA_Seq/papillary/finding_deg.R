library(DESeq2)
library(pheatmap)
library(gplots)
library(ConsensusClusterPlus)

dds <- DESeqDataSetFromMatrix(countData = exp_prof_diff, colData = sample.info, design = ~pat.nums + type)
colData(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
rld <- rlogTransformation(dds)
vs <- vst(dds)
nt <- normTransform(dds)

res <- results(dds)
summary(res)
res <- res[order(res$padj), ]

diff.genes = list() ##Here name stands for log 2 fold change cut off
exp_fpqm = exp_fpqm[match(rownames(assay(dds)), rownames(exp_fpqm)),] 
sum(rownames(exp_fpqm) != rownames(assay(dds))) == 0
sum(rownames(exp_fpqm) == rownames(assay(dds))) == length(rownames(assay(dds)))

diff.genes[[2]] = rownames(res[abs(res$log2FoldChange) > 2 & res$padj < 0.01 ,])
diff.genes[[3]] = rownames(res[abs(res$log2FoldChange) > 3 & res$padj < 0.01 ,])
diff.genes[[4]] = rownames(res[abs(res$log2FoldChange) > 4 & res$padj < 0.01 ,])
diff.genes[[5]] = rownames(res[abs(res$log2FoldChange) > 5 & res$padj < 0.01 ,])
names(diff.genes) = c(2,3,4,5)

top.val.genes = list()
top.val.genes[['rld']] = get.imp.genes(1,assay(rld),5000)
top.val.genes[['nt']] = get.imp.genes(1,assay(nt),5000)
top.val.genes[['vs']] = get.imp.genes(1,assay(vs),5000)
top.val.genes[['fpqm']] = get.imp.genes(1,exp_fpqm[,diff.ids],5000)
top.val.genes[['fpqm_log']] = get.imp.genes(1,log2(exp_fpqm[,diff.ids]+1),5000)

top.var.genes = list()
top.var.genes[['rld']] = get.imp.genes(3,assay(rld),5000)
top.var.genes[['nt']] = get.imp.genes(3,assay(nt),5000)
top.var.genes[['vs']] = get.imp.genes(3,assay(vs),5000)
top.var.genes[['fpqm']] = get.imp.genes(3,exp_fpqm[,diff.ids],5000)
top.var.genes[['fpqm_log']] = get.imp.genes(3,log2(exp_fpqm[,diff.ids]+1),5000)

top.mad.genes = list()
top.mad.genes[['rld']] = get.imp.genes(2,assay(rld),5000)
top.mad.genes[['nt']] = get.imp.genes(2,assay(nt),5000)
top.mad.genes[['vs']] = get.imp.genes(2,assay(vs),5000)
top.mad.genes[['fpqm']] = get.imp.genes(2,exp_fpqm[,diff.ids],5000)
top.mad.genes[['fpqm_log']] = get.imp.genes(2,log2(exp_fpqm[,diff.ids]+1),5000)




View(t(counts(dds)[diff.genes[[5]],]))
d = dist(t(counts(dds)[diff.genes[[5]],]))
exp_fpqm = exp_fpqm[rowSums(exp_fpqm) > 1, ]
d1 = dist(t(as.matrix(exp_fpqm[,diff.ids][diff.genes[[5]],])))
plot(hclust(d1), labels = sample.info$type)
plotPCA(nt, intgroup = 'type')
d = dist(t(assay(rld)[diff.genes[[2]],]))
pheatmap(d1, cluster_rows = F, cluster_cols = T, show_colnames  = T, labels_col = sample.info$type )
heatmap.2(exp_fpqm[,diff.ids][diff.genes[[5]],], labCol = sample.info$type)
heatmap.2(t(assay(nt)[diff.genes[[3]],]), labCol = sample.info$type)  


##Consenus Clustering
mads = apply(exp_fpqm[,diff.ids], 1, mad)
d_fp=exp_fpqm[rev(order(mads))[1:5000],]
d_fp = sweep(d_fp, 1, apply(d_fp, 1, median, na.rm=T))
title = tempdir()
con_res = ConsensusClusterPlus(d_fp,maxK=6,reps=50,pItem=0.8,pFeature=1,title = title,
                               clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")

library(factoextra)
library(FactoMineR)
dt = t(exp_fpqm[,diff.ids])
rownames(dt) = sample.info$type
pc = PCA(dt)
fviz_pca_ind(pc)

