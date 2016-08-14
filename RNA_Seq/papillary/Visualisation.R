##This is a file for generating visualing plots
##input would be a list of respective dfs/assays from respective files containing of rlogtrans, normTrans,
##VarStabTrans, fpqm, log2fpqm

##Using variables from finding_deg or onlytumors.R
##Only match_control tumor


##PCA##
genes.list = match.control.tumor[['genes']]
dfs.list = match.control.tumor[['dfs']]

pca_diff = list()
pca_diff = lapply(dfs.list, function(x)
  {
  get.list.plotPCA(c(2:5), x , diff.genes, 'diff', type = 1, intgroup = 'type', sample.info = sample.info)
  
})
#rld
plotPCA(rld,intgroup = 'type', title = 'overall')
pca_rld <- list()
pca_rld[['diff']] <- get.list.plotPCA(c(2:5), dfs.list[['rld']] , diff.genes, 'diff') ##Gets us the list of 
##plotPCAs applied on assay for different degs based on log fold cut off
multiplot(pca_rld[['diff']], cols = 2)

pca_rld[['var']] <- get.list.plotPCA(c(100,500,1000,2000,3000,4000,5000), dfs.list[['rld']], 
                                     genes.list[['top.var.genes']]$rld, 'var', type = 2) ##Gets us the list of 
##plotPCAs applied on assay for top x variable genes based on mad
multiplot(pca_rld[['var']], cols = 3)

pca_rld[['val']] <- get.list.plotPCA(c(100,500,1000,2000,3000,4000,5000), rld, 
                                     top.val.genes$rld, 'val', type = 2) ##Gets us the list of 
##plotPCAs applied on assay for top x genes based on their values
multiplot(pca_rld[['val']], cols = 3)

pca_rld[['mad']] <- get.list.plotPCA(c(100,500,1000,2000,3000,4000,5000), rld, 
                                     top.mad.genes$rld, 'mad', type = 2) ##Gets us the list of 
##plotPCAs applied on assay for top x genes based on their values
multiplot(pca_rld[['mad']], cols = 3)

##nt
plotPCA(nt, intgroup = 'type', title = 'overall')
pca_nt = list()
pca_nt[['diff']] <- get.list.plotPCA(c(2:5), nt, diff.genes, 'diff')
multiplot(pca_nt[['diff']], cols = 2)

pca_nt[['var']] <- get.list.plotPCA(c(100,500,1000,2000,3000,4000,5000), nt, 
                                     top.var.genes$nt, 'var', type = 2)  
multiplot(pca_nt[['var']], cols = 3)
pca_nt[['val']] <- get.list.plotPCA(c(100,500,1000,2000,3000,4000,5000), nt, 
                                     top.val.genes$nt, 'val', type = 2) 
multiplot(pca_nt[['val']], cols = 3)

pca_nt[['mad']] <- get.list.plotPCA(c(100,500,1000,2000,3000,4000,5000), nt, 
                                    top.mad.genes$nt, 'mad', type = 2) 
multiplot(pca_nt[['mad']], cols = 3)

##vs
plotPCA(vs, intgroup = 'type', title = 'overall')
pca_vs = list()
pca_vs[['diff']] <- get.list.plotPCA(c(2:5), vs, diff.genes, 'diff')
multiplot(pca_vs[['diff']], cols = 2)
pca_vs[['var']] <- get.list.plotPCA(c(100,500,1000,2000,3000,4000,5000), vs, 
                                    top.var.genes$vs, 'var', type = 2)  
multiplot(pca_vs[['var']], cols = 3)
pca_vs[['val']] <- get.list.plotPCA(c(100,500,1000,2000,3000,4000,5000), vs, 
                                    top.val.genes$vs, 'val', type = 2)  
multiplot(pca_vs[['val']], cols = 3)

pca_vs[['mad']] <- get.list.plotPCA(c(100,500,1000,2000,3000,4000,5000), vs, 
                                    top.val.genes$vs, 'mad', type = 2)  
multiplot(pca_vs[['mad']], cols = 3)


##fpqm
plotPCA(exp_fpqm[,diff.ids], intgroup = 'type', title = 'overall', sample.info = sample.info)
pca_fpqm = list()
pca_fpqm[['diff']] <- get.list.plotPCA(c(2:5), exp_fpqm[,diff.ids], diff.genes, 'diff', sample.info=sample.info)
multiplot(pca_fpqm[['diff']], cols = 2)
pca_fpqm[['var']] <- get.list.plotPCA(c(100,500,1000,2000,3000,4000,5000), exp_fpqm[,diff.ids], 
                                    top.var.genes$fpqm, 'var', type = 2, sample.info = sample.info)  
multiplot(pca_fpqm[['var']], cols = 3)
pca_fpqm[['val']] <- get.list.plotPCA(c(100,500,1000,2000,3000,4000,5000), exp_fpqm[,diff.ids], 
                                      top.val.genes$fpqm, 'val', type = 2, sample.info = sample.info)  
multiplot(pca_fpqm[['val']], cols = 3)
pca_fpqm[['mad']] <- get.list.plotPCA(c(100,500,1000,2000,3000,4000,5000), exp_fpqm[,diff.ids], 
                                      top.mad.genes$fpqm, 'mad', type = 2, sample.info = sample.info)  
multiplot(pca_fpqm[['mad']], cols = 3)

##fpqm_log
plotPCA(log2(exp_fpqm[,diff.ids]+1), intgroup = 'type', title = 'overall', sample.info = sample.info)
pca_fpqm_log = list()
pca_fpqm_log[['diff']] <- get.list.plotPCA(c(2:5), log2(exp_fpqm[,diff.ids]+1), diff.genes, 'diff', sample.info=sample.info)
multiplot(pca_fpqm_log[['diff']], cols = 2)
pca_fpqm_log[['var']] <- get.list.plotPCA(c(100,500,1000,2000,3000,4000,5000), log2(exp_fpqm[,diff.ids]+1), 
                                      top.var.genes$fpqm_log, 'var', type = 2, sample.info = sample.info)  
multiplot(pca_fpqm_log[['var']], cols = 3)
pca_fpqm_log[['val']] <- get.list.plotPCA(c(100,500,1000,2000,3000,4000,5000), log2(exp_fpqm[,diff.ids]+1), 
                                      top.val.genes$fpqm_log, 'val', type = 2, sample.info = sample.info)  
multiplot(pca_fpqm_log[['val']], cols = 3)
pca_fpqm_log[['mad']] <- get.list.plotPCA(c(100,500,1000,2000,3000,4000,5000), log2(exp_fpqm[,diff.ids]+1), 
                                      top.mad.genes$fpqm_log, 'mad', type = 2, sample.info = sample.info)  
multiplot(pca_fpqm_log[['mad']], cols = 3)


#####HeatMaps###########
d = dist(t(assay(rld)[diff.genes[[2]],]))
pheatmap(d, cluster_rows = F, cluster_cols = T,  labels_col = sample.info$type, clustering_method = 'ward.D' )
heatmap.2(assay(vs)[diff.genes[[3]],], labCol = sample.info$type)  
plot(hclust(d), labels = sample.info$type)
