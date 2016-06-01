##########Working with DEGS########################
source('preprocess.R')
common.degs = Reduce(intersect, diff.expressed.all.canc.corr) #genes which are degs in all cancers
degs.int.bps = find.bp.deg.int(diff.expressed.all.canc.corr, bp.genes.ensembl) #intersection of degs of each cancer with the genes of biological process
l.bps = length(bp.genes.ensembl) #total different BPs
l.tumors = length(tissues.names) #total different tumors
lengths.bps.degs = list() #lengths of degs of each cancer which intersect with genes of BPs
for(i in seq(l.bps))
  lengths.bps.degs[[i]] = find.lengths(degs.int.bps[[i]])
names(lengths.bps.degs) = names(degs.int.bps)


degs.up.int.bps = find.bp.deg.int(diff.expressed.up.all.canc.corr, bp.genes.ensembl) #upregulated degs in a cancer that intersect with genes of biological process
degs.down.int.bps = find.bp.deg.int(diff.expressed.down.all.canc.corr, bp.genes.ensembl) #downregulated degs in a cancer that intersect with genes of biological process

df.lengths.bps.degs = create.df.length(lengths.bps.degs, l.bps, l.tumors, short.names)#data frame with the lengths of total degs within each BP
library(ggplot2)
##Plotting the lengths of degs in each cancer under each BP
ggplot.needs(df.lengths.bps.degs, colnames(df.lengths.bps.degs)[1], colnames(df.lengths.bps.degs)[2] , 'blue')
ggplot.needs(df.lengths.bps.degs, colnames(df.lengths.bps.degs)[1], colnames(df.lengths.bps.degs)[3] , 'red')
ggplot.needs(df.lengths.bps.degs, colnames(df.lengths.bps.degs)[1], colnames(df.lengths.bps.degs)[4] , 'green')
ggplot.needs(df.lengths.bps.degs, colnames(df.lengths.bps.degs)[1], colnames(df.lengths.bps.degs)[5] , 'yellow')
ggplot.needs(df.lengths.bps.degs, colnames(df.lengths.bps.degs)[1], colnames(df.lengths.bps.degs)[6] , 'pink')
ggplot.needs(df.lengths.bps.degs, colnames(df.lengths.bps.degs)[1], colnames(df.lengths.bps.degs)[7] , 'orange')
ggplot.needs(df.lengths.bps.degs, colnames(df.lengths.bps.degs)[1], colnames(df.lengths.bps.degs)[8] , 'purple')

common.bps.degs <- sapply(degs.int.bps, function(x) ##Checking under each BP whehther is there any
  {                                 #gene that is a DEG in all cancer for that process
    Reduce(intersect, x)
})

unique.degs.bps <- find.unique.degs.bps(degs.int.bps) #the genes which are unique for a cancer in a particular process
unique.degs.up.bps <- find.intersection(unique.degs.bps, diff.expressed.up.all.canc.corr)#above unique genes which are upregulated
unique.degs.down.bps <- find.intersection(unique.degs.bps, diff.expressed.down.all.canc.corr)#above genes which are downregulated

lengths.uni.degs.bp = list()
for(i in seq(l.bps))
    lengths.uni.degs.bp[[i]] = find.lengths(unique.degs.bps[[i]])
names(lengths.uni.degs.bp) = names(unique.degs.bps)

##lengths of the degs in each BP and corresponding df to them
df.lengths.unique.bps.degs = create.df.length(lengths.uni.degs.bp, l.bps, l.tumors, short.names)

##plotting the above number of degs of each cancer for each BP
ggplot.needs(df.lengths.unique.bps.degs, colnames(df.lengths.unique.bps.degs)[1], colnames(df.lengths.unique.bps.degs)[2] , 'blue')
ggplot.needs(df.lengths.unique.bps.degs, colnames(df.lengths.unique.bps.degs)[1], colnames(df.lengths.unique.bps.degs)[3] , 'red')
ggplot.needs(df.lengths.unique.bps.degs, colnames(df.lengths.unique.bps.degs)[1], colnames(df.lengths.unique.bps.degs)[4] , 'green')
ggplot.needs(df.lengths.unique.bps.degs, colnames(df.lengths.unique.bps.degs)[1], colnames(df.lengths.unique.bps.degs)[5] , 'yellow')
ggplot.needs(df.lengths.unique.bps.degs, colnames(df.lengths.unique.bps.degs)[1], colnames(df.lengths.unique.bps.degs)[6] , 'pink')
ggplot.needs(df.lengths.unique.bps.degs, colnames(df.lengths.unique.bps.degs)[1], colnames(df.lengths.unique.bps.degs)[7] , 'orange')
ggplot.needs(df.lengths.unique.bps.degs, colnames(df.lengths.unique.bps.degs)[1], colnames(df.lengths.unique.bps.degs)[8] , 'purple')

###Working with gene ids
genes.ids = read.csv('gene_names.csv')
tot.mapped = genes.ids$Ensembl.Gene.ID ##all genes with their gene ids

##Experiment with lymphoma
lymphoma.up.degs.id = as.character(genes.ids$Associated.Gene.Name[match(intersect(diff.expressed.up.all.canc.corr$Lymphoma, tot.mapped), genes.ids$Ensembl.Gene.ID)])
lymphoma.down.degs.id = as.character(genes.ids$Associated.Gene.Name[match(intersect(diff.expressed.down.all.canc.corr$Lymphoma, tot.mapped), genes.ids$Ensembl.Gene.ID)])
write(lymphoma.up.degs.id, 'lymph_gene_up.txt')
write(lymphoma.down.degs.id, 'lymph_gene_down.txt')
setwd('all_tissues/degs_and_bps/')
write.lists(degs.int.bps, 2)

##mapping various ensembl gene ids found above to gene names
genes.ids.map.up = map.genes.lists(degs.up.int.bps, genes.ids) 
genes.ids.map.down = map.genes.lists(degs.down.int.bps, genes.ids)
unique.genes.ids.map.up = map.genes.lists(unique.degs.up.bps, genes.ids)
unique.genes.ids.map.down = map.genes.lists(unique.degs.down.bps, genes.ids)

###writing the lists
setwd('up')
write.lists.2nd(genes.ids.map.up)
setwd('../down/')
write.lists.2nd(genes.ids.map.down)

###Cell Cycle
p53_up = check.for.pattern('P53', genes.ids.map.up$cell_cycle)
p53_down = check.for.pattern('P53', genes.ids.map.down$cell_cycle)

pten_up = check.for.pattern('PTEN', genes.ids.map.up$cell_cycle)
pten_down = check.for.pattern('PTEN', genes.ids.map.down$cell_cycle)

###Writing
write(tot.mapped, 'gene_symbs.txt')
bps.genes.ids = map.genes.lists(bp.genes.ensembl, genes.ids)
write.lists(bps.genes.ids, 2)
