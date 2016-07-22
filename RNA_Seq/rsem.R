rsem = read.delim('~/Downloads/data/gdac.broadinstitute.org_PRAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0/PRAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt', header = T)
rsem = a
rownames(rsem) = rsem[,1]
rsem = rsem[,-1]
rsem = rsem[,seq(1,1650,3)]
rsem = rsem[-1,]
samp = list()
samp[['pat_ids']] = colnames(rsem)
samp[['lab']] = sapply(samp[['pat_ids']],function(x) ##Normal or Tumor 
{
  unlist(strsplit(x,split = '.', fixed = T))[4]
})
samp[['type']] = sapply(samp[['lab']], function(x)
{
  if(x == '11A' | x == '11B')
    x = 'N'
  else
    x = "T"
})
samp[['pat']] = sapply(samp[['pat_ids']],function(x)
{
  m = unlist(strsplit(x,split = '.', fixed = T))
  paste(m[1],m[2],m[3], sep = '-')
})
samp[['normal.indexes']] = which(samp[['type']] == 'N')
samp[['tumor.indexes']] = which(samp[['type']] == 'T')
samp[['diff.ids']] = c()

for(i in samp[['normal.indexes']])
{
  for(j in samp[['tumor.indexes']])
  {
    if(samp[['pat']][i] == samp[['pat']][j])
      samp[['diff.ids']] = c(samp[['diff.ids']],i,j)
  }
}

library(DESeq2)
rsem_req = rsem[,samp$diff.ids]
col_data = lapply(seq(4), function(x)
  {
  samp[[x]][samp$diff.ids]
})
names(col_data) = names(samp)[1:4]
col_data = data.frame(col_data)
j = 1
pat = col_data$pat[1]
for(i in seq_along(col_data$pat_ids))
{
  if(col_data$pat[i] != pat)
  {
    j = j + 1
    pat = col_data$pat[i]
  }
  col_data$pat.short.ids[i] = j
}
col_data$comb = mapply(function(x,y){
  paste(as.character(x),y,sep = '')},
  col_data$pat.short.ids, col_data$type)

for(i in seq_along(colnames(rsem_req)))
{
  rsem_req[,i] = as.integer(rsem_req[,i])
}
dds_rs <- DESeqDataSetFromMatrix(countData = rsem_req, colData = col_data, design = ~type)
dds_rs <- dds_rs[ rowSums(counts(dds_rs)) > 1, ]
dds_rs <- DESeq(dds_rs)
res_rs <- results(dds_rs)
res_rs <- res_rs[order(res_rs$padj), ]
diff_rs = rownames(res_rs[res_rs$padj < 0.05 & abs(res_rs$log2FoldChange) > 1,])
nt_rs = normTransform(dds_rs)


d <- dist(t(assay(nt_rs)[diff_rs,]))
plot(hclust(d), labels = col_data$type)
df <- t(assay(nt_rs)[diff_rs,])
rownames(df) = col_data$type
pr <- PCA(df)
fviz_pca_ind(pr)

