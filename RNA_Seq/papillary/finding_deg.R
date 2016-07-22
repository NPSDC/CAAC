library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = exp_prof_diff, colData = sample.info, design = ~pat.nums + type)
colData(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
rld <- rlogTransformation(dds)
