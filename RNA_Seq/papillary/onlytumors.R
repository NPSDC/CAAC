###Tumor samples
library(DESeq2)
pat.stages = rep('NA', length(sample.ids))
indexes.stages = match(pats[tumor.indexes], names(pat_stages))
pat.stages[tumor.indexes] = pat_stages[indexes.stages] #from extract_stage.R
names(pat.stages) = sample.ids

df.stages = data.frame(sample.id=sample.ids[tumor.indexes], patient.id = pats[tumor.indexes], stage =
                         unlist(pat_stages[indexes.stages]), short.id = c(1:289))
rownames(df.stages) = c(1:289)


exp_prof_tumor = exp_prof[,tumor.indexes]
colnames(exp_prof_tumor) = c(1:289)
dds_tumor = DESeqDataSetFromMatrix(countData = exp_prof_tumor, colData = df.stages, design = ~stage)
colData(dds_tumor)
dds_tumor = dds_tumor[match(rownames(dds), rownames(dds_tumor)),]
dds_tumor <- dds_tumor[ rowSums(counts(dds_tumor)) > 1, ]
dds_tumor <- DESeq(dds_tumor)

only.tumor.dfs[['rld']] = rlogTransformation(dds_tumor)
remove(exp_prof_tumor)

exp_fpqm_tumor <- exp_fpqm[match(rownames(dds_tumor), rownames(exp_fpqm)), tumor.indexes]

only.tumor = list()
only.tumor[['dfs']] = list()
only.tumor[['dfs']][['vs']] = vst(dds_tumor)
only.tumor[['dfs'']][['nt']] = normTransform(dds_tumor)
only.tumor[['dfs']][['rld']] = rlogTransformation(dds_tumor)

only.tumor[['genes']] = list()
only.tumor[['genes']][['top.mad.genes']] = create.list.imp.genes()
only.tumor[['genes']][['top.var.genes']] = list()
only.tumor[['genes']][['top.val.genes']] = list()

only.tumor.genes$top.mad.genes[['']]
top.mad.genes[['fpqm_log_tumor']] = get.imp.genes(3,log2(exp_fpqm[,tumor.indexes]+1),5000)
plotPCA(log2(exp_fpqm[,tumor.indexes]+1), intgroup = 'stage', title = 'overall', sample.info = df.stages)
pca_fpqm_log_tumor[['mad']] <- get.list.plotPCA(c(100,500,1000,2000,3000,4000,5000), log2(exp_fpqm[,tumor.indexes]+1), 
                                                top.mad.genes$fpqm_log_tumor, 'mad', type = 2, sample.info = df.stages,
                                                intgroup = 'stage')
multiplot(pca_fpqm_log_tumor[['mad']], cols = 3)
vs_tumor <- vst(dds_tumor)
plotPCA(vs_tumor[diff.genes[[2]],], intgroup = 'stage', title = 'overall')
