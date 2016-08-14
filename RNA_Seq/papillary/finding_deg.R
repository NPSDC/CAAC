library(DESeq2)
library(pheatmap)
library(gplots)
library(ConsensusClusterPlus)

dds <- DESeqDataSetFromMatrix(countData = exp_prof_diff, colData = sample.info, design = ~pat.nums + type)
colData(dds)
match_control_tum = list()
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
dds_data <- read.csv('dds.csv')
exp_fpqm = read.csv('fpqm.csv')
exp_fpqm = exp_fpqm[match(rownames(assay(dds)), rownames(exp_fpqm)),]
exp_fpqm <- exp_fpqm[ rowSums(exp_fpqm) > 1, ]

match.control.tumor = list()
match.control.tumor[['dfs']] = list()

match.control.tumor[['dfs']][['rld']] <- rlogTransformation(dds)
match.control.tumor[['dfs']][['vs']] <- vst(dds)
match.control.tumor[['dfs']][['nt']] <- normTransform(dds)
match.control.tumor[['dfs']][['fpqm']] <- exp_fpqm[,diff.ids]
match.control.tumor[['dfs']][['fpqm_log']] <- log2(exp_fpqm[,diff.ids]+1)

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

match.control.tumor[['genes']] = list()
match.control.tumor[['genes']][['top.val.genes']] = create.list.imp.genes(match.control.tumor[['dfs']], 1, 5000)
match.control.tumor[['genes']][['top.mad.genes']] = create.list.imp.genes(match.control.tumor[['dfs']], 2, 5000)
match.control.tumor[['genes']][['top.var.genes']] = create.list.imp.genes(match.control.tumor[['dfs']], 3, 5000)


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
dt = t(assay(vs)[top.mad.genes$vs[1:5000],])
rownames(dt) = sample.info$type
pc = PCA(dt)
fviz_pca_ind(pc)

