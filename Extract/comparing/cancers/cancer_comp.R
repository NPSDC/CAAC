source('extract.R')
preprocess <- function(name, data, col.index, level.index)
{
  if(is.na(match(name, unique(data[,col.index]))))
  {
    print("wrong name")
    return()
  }
  type = data[data[,col.index] == name, ]
  type$Level = set.level(type[,level.index], levels(type[,level.index]), 2)
  return(type)
}

initialize.empty <- function(unique.genes, levels, times.to.repeat, no.unique.genes)
  ###Initializes the data frame creating an empty data frame specific for cancer as well 
  ###function reduce.levels
  ###Inputs are self explanatory
{
  no.unique.genes = length(unique.genes)
  canc.specific = data.frame(Gene = rep(unique.genes, each = times.to.repeat),  
                             Level = factor(rep(levels, no.unique.genes)), 
                             Count = rep(0L, no.unique.genes*times.to.repeat),
                             Total = rep(0L, no.unique.genes*times.to.repeat))  
  return(canc.specific)
}

reduce.levels <- function(canc.name, levels, combine.conc)
  ###Reudces the 4 levels of cancer to 2 levels
  ###canc.name- name cancer
  ###levels - levels of conc to keep(character vector)
  ###combine.conc - list concentration levels to be combined in each level(list)
{
  canc.type <- preprocess(canc.name, canc, 2, 3)
  l = length(levels)
  unique.genes = unique(canc.type$Gene)
  no.unique.genes = length(unique(canc.type$Gene))
  canc.specific = initialize.empty(unique.genes, levels, l, no.unique.genes)
  if(length(levels) != length(combine.conc))
  {
    print("Check Params")
    return()
  }
  
  indexes.levels = list() ##contains the indexes of the various levels
  for(i in seq(length(combine.conc)))
  {
    indexes.levels[[i]] = list()
    for(j in seq(length(combine.conc[[i]])))
        indexes.levels[[i]][[j]] = which(canc.type$Level == combine.conc[[i]][j])
  }
  k = 1
  for(gene in unique.genes)
  {
    gene.index = which(gene == canc.type$Gene) #gene index of the gene in the cancer set to be reduced
    gene.and.level.index = c()
    for(i in seq(l))
    {
      for(j in seq(length(indexes.levels[[i]])))
      {
        gene.and.level.index = c(gene.and.level.index,
                                 gene.index[!is.na(match(gene.index, indexes.levels[[i]][[j]]))])
      }
      canc.specific$Gene[k + i - 1] = gene
      canc.specific$Count[k + i - 1] = sum(canc.type$Count[gene.and.level.index])
      gene.and.level.index = c()
    }
    canc.specific$Total[k] = canc.specific$Total[k + 1] = 
      canc.specific$Count[k] + canc.specific$Count[k+1]
    
    levels.vals = sapply(c(0.5, 0.75, 0.9, 1), function(x)
    {
      calculate.level(x,canc.specific$Count[k], canc.specific$Count[k+1], levels)
    })
    canc.specific$actual.level.50[c(k, k + 1)] <- levels.vals[1]
    canc.specific$actual.level.75[c(k, k + 1)] <- levels.vals[2]
    canc.specific$actual.level.90[c(k, k + 1)] <- levels.vals[3]
    canc.specific$actual.level.100[c(k, k + 1)] <- levels.vals[4]
    k = k + l
      
  }
  return(canc.specific)
}
renal.mod.canc = reduce.levels('renal cancer', c('Present', 'Not detected'),
                               list('Present', 'Not detected'))
glioma.mod.canc = reduce.levels('glioma', c('Present', 'Not detected'),
                                list('Present', 'Not detected'))
lymphoma.mod.canc = reduce.levels('lymphoma', c('Present', 'Not detected'), 
                                  list('Present', 'Not detected'))


lymph.normal = preprocess('lymph node', normal, 2, 4)
renal.normal = preprocess('kidney', normal, 2, 4)
glioma.normal = preprocess('glial cells', normal, 3, 4)


find.diff.expressed <- function(normal, cancer, level.indexes, genes.indexes)
{
  genes = intersect(unique(normal[,genes.indexes[1]]), unique(cancer[,genes.indexes[2]]))
  diff.expressed = c()
  for(i in genes)
  {
    canc.level = as.character(cancer[match(i, cancer[,genes.indexes[2]]), level.indexes[2]])
    normal.gene.indexes = which(normal[,genes.indexes[1]] == i)
    normal.level = normal[,level.indexes[1]][normal.gene.indexes][1]
    if((sum(normal[,level.indexes[1]][normal.gene.indexes] == normal.level) == length(normal.gene.indexes))
       && canc.level != normal.level) 
      diff.expressed = c(diff.expressed, i)
  }
  return(diff.expressed)
}
bile.diff = find.diff.expressed(bile.duct, mod.canc.liver.2, c(4,5), c(1,1))
hep.diff = find.diff.expressed(hepatocytes_all, mod.canc.liver.2, c(4,5), c(1,1))
diff.expressed.4 = mapply(function(x,y)
  {
  find.diff.expressed(x,y,c(4,5), c(1,1))
}, list(lymph.normal, renal.normal, glioma.normal), list(lymphoma.mod.canc, renal.mod.canc, glioma.mod.canc))
diff.expressed.4[['liver']] = intersect(bile.diff, hep.diff)
names(diff.expressed.4) = c('lymph', 'renal', 'glioma', 'liver')
diff.expressed.all = intersect(intersect(diff.expressed.4$lymph, diff.expressed.4$renal), 
                               intersect(diff.expressed.4$glioma, diff.expressed.4$liver))
diff.expressed.all.indexes = match(intersect(diff.expressed.all, total.genes), total.genes)
res.mca.diff.expressed.all = MCA(data.cancer.gene[,diff.expressed.all.indexes])
fviz_mca_ind(res.mca.diff.expressed.all)
hcpc.diff.expressed.all = HCPC(res.mca.diff.expressed.all)
plot(hcpc.diff.expressed.all, choice = 'tree')
dev.off()
diff.genes.indexes.tt = c(26,27)

renal.index = match(diff.expressed.4$renal, total.genes)
res.mca.renal.diff = MCA(data.cancer.gene[,renal.index])
fviz_mca_ind(res.mca.renal.diff)

liver.index = match(intersect(diff.expressed.4$liver,diff.expressed.4$renal), total.genes)
res.mca.liver.diff = MCA(data.cancer.gene[,liver.index])
fviz_mca_ind(res.mca.liver.diff)

glioma.index = match(diff.expressed.4$liver, total.genes)
res.mca.glioma.diff = MCA(data.cancer.gene[,glioma.index])
fviz_mca_ind(res.mca.glioma.diff)

lymphoma.index = match(intersect(diff.expressed.4$lymph, total.genes), total.genes)
res.mca.lymph.diff = MCA(data.cancer.gene[, lymphoma.index])
fviz_mca_ind(res.mca.lymph.diff)
genes.only.lymph.index = c() #contains indexes of lymphoma.index which have a particular expression
#only across lymphoma
for(i in seq(length(lymphoma.index)))
{
  val = as.character(data.cancer.gene[10,lymphoma.index[i]])
  #print(val)
  if(sum(data.cancer.gene[c(c(1:9),c(11:20)), lymphoma.index[i]] == val) == 0)
  {
   # print('yes')
    genes.only.lymph.index = c(genes.only.lymph.index, i)
  }
}
#opposite of gene.only.lymph.index but contains lymphoma.index's indexes
not.only.lymph = setdiff(lymphoma.index, lymphoma.index[genes.only.lymph.index]) 
res.MCA.not.diff.lymph = MCA(data.cancer.gene[,not.only])
fviz_mca_ind(res.MCA.not.diff.lymph)
hc = HCPC(res.MCA.not.diff.lymph)
plot(hc,choice='tree')
dev.off()
res.MCA.diff.lymph = MCA(data.cancer.gene[,lymphoma.index[genes.only.lymph.index]])
fviz_mca_ind(res.MCA.diff.lymph)

all.cancers <- cancer.all.gene.wise[,c("Gene","Giloma.Cancer",'Liver.Cancer', 'Lymphoma.Cancer', 'Renal.Cancer')]
data.all.cancers <- t(all.cancers)
colnames(data.all.cancers) <- all.cancers$Gene
data.all.cancers <- data.all.cancers[-1,]

library(FactoMineR)
library(factoextra)

mca.4.canc <- MCA(data.all.cancers)
hcpc.4.canc <- HCPC(mca.4.canc)
find.sim <- function(data.frame.main,gene.index, whom.match.index, to.matched.index)
{
   ans <- lapply(to.matched.index, function(x)
   {
     data.frame.main[,gene.index][data.frame.main[,whom.match.index] == data.frame.main[,x]]
   }
     )
   names(ans) = colnames(data.frame.main)[to.matched.index]
   return(ans)
}
liver.match = find.sim(all.cancers, 1, 3, c(2,4,5))
giloma.match = find.sim(all.cancers, 1, 2, c(3,4,5))
lymp.match = find.sim(all.cancers, 1, 4, c(2,3,5))
renal.match = find.sim(all.cancers, 1, 5, c(2,3,4))


bp.data.frame = list()
for(i in seq(length(data.frame.mca)))
{
  bp.data.frame[[i]] = t(data.frame.mca[[i]])
  bp.data.frame[[i]] = data.frame(bp.data.frame[[i]])
  bp.data.frame[[i]]$Gene = rownames(bp.data.frame[[i]])
}
names(bp.data.frame) = biological.processes

bp.liver.match = lapply(bp.data.frame, function(x)
  {
    find.sim(x, 21, 8, c(6, 10, 15))
})
names(bp.liver.match) = biological.processes

bp.liver.intersect = lapply(bp.liver.match, function(x) #contains all the common genes for all bps for 
  {                                                     #4 bps
  intersect(intersect(x[[1]], x[[2]]), x[[3]])
})
names(bp.liver.intersect) = biological.processes

lengths.bp.intersect = lapply(bp.liver.intersect, length)
names(lengths.bp.intersect) = biological.processes

#contains genes common across 4 types in all biological processes
genes.all.common = intersect(intersect(colnames(data.cancer.gene),liver.match$Giloma.Cancer), 
                             intersect(liver.match$Lymphoma.Cancer, liver.match$Renal.Cancer))
genes.union = union(union(union(union(bp.liver.intersect$apoptosis, bp.liver.intersect$cell_adhesion), 
                                union(bp.liver.intersect$cell_cycle, bp.liver.intersect$immune_response)),
                          union(bp.liver.intersect$metabolic,bp.liver.intersect$signal_transduction)),
                    bp.liver.intersect$transport)
data.frame.gene.union = data.cancer.gene[, match(genes.union, colnames(data.cancer.gene))]

res.mca.union = MCA(data.frame.gene.union)
fviz_mca_ind(res.mca.union)

data.frame.common.bps = lapply(bp.liver.intersect, function(x)  ##contains the data frame which contains
  #genes that are common in 4 cancersprocess wise 
  {
  data.cancer.gene[,match(x, colnames(data.cancer.gene))]
})
names(data.frame.common.bps) = biological.processes

res.mca.common.bps = lapply(data.frame.common.bps, MCA)
names(res.mca.common.bps) = biological.processes

fviz_mca_ind(res.mca.common.bps$transport)


data.frame.gene.except.union = data.cancer.gene[, match(setdiff(colnames(data.cancer.gene), genes.union)
                                                        , colnames(data.cancer.gene))]
res.mca.except.union = MCA(data.frame.gene.except.union)
fviz_mca_ind(res.mca.except.union)



#####Finding difference b/w renal and liver #####################
diff.expressed.absent = mapply(setdiff, diff.expressed, diff.expressed.present)
liver.diff = lapply(liver.match, function(x) #contains disimilar genes b/w liver and other cancers
  {
    setdiff(all.cancers$Gene, x)
})
names(liver.diff) = names(liver.match)

liver.diff.up = lapply(liver.diff, function(x)
  {
  genes.diff.present = all.cancers$Gene[which(all.cancers$Liver.Cancer == 'Present')]
  x[!is.na(match(x, genes.diff.present))]
})
names(liver.diff.up) = names(liver.match)

liver.diff.down = mapply(setdiff, liver.diff, liver.diff.up)
names(liver.diff.down) = names(liver.match)

liver.diff.up.actual = lapply(liver.diff.up, function(x)
  {
  intersect(x, diff.expressed.present$`50`)
})
names(liver.diff.up.actual) = names(liver.match)

liver.diff.down.actual = lapply(liver.diff.down, function(x)
  {
    setdiff(x, diff.expressed.absent$`50`)
})
names(liver.diff.down.actual) = names(liver.match)

write(liver.diff.down$Renal.Cancer, 'comparing/cancers/liver/liver_renal_down.txt')
write(liver.diff.up$Renal.Cancer, 'comparing/cancers/liver/liver_renal_up.txt')
write(liver.diff.up.actual$Renal.Cancer, 'comparing/cancers/liver/liver_renal_up_actual.txt')
write(liver.diff.down.actual$Renal.Cancer, 'comparing/cancers/liver/liver_renal_down_actual.txt')

###Finding difference b/w  processes
common.cell.metabolic = intersect(bp.genes.ensembl$cell_cycle, bp.genes.ensembl$metabolic)
not.common.cell.metabolic = setdiff(union(bp.genes.ensembl$cell_cycle, bp.genes.ensembl$metabolic),
                                    common.cell.metabolic)
metabolic.unique = setdiff(bp.genes.ensembl$metabolic, common.cell.metabolic)
cell.cycle.unique = setdiff(bp.genes.ensembl$cell_cycle, common.cell.metabolic)
data.frame.mca.common.cellcycle.metabolic = 
  data.cancer.gene[, match(intersect(common.cell.metabolic, colnames(data.cancer.gene)),
                           colnames(data.cancer.gene))]
data.frame.mca.not.common.cellcycle.metabolic = 
  data.cancer.gene[, match(intersect(not.common.cell.metabolic, colnames(data.cancer.gene)),
                           colnames(data.cancer.gene))]
data.frame.mca.metabolic.unique = data.cancer.gene[,match(intersect(metabolic.unique, colnames(data.cancer.gene)),
                                                          colnames(data.cancer.gene))][1:200]
data.frame.mca.cell_cycle.unique = data.cancer.gene[, match(intersect(cell.cycle.unique, colnames(data.cancer.gene)), 
                                                            colnames(data.cancer.gene))][1:400]

res.mca.not.common.cellcycle.metabolic = MCA(data.frame.mca.not.common.cellcycle.metabolic)
res.mca.common.cellcycle.metabolic = MCA(data.frame.mca.common.cellcycle.metabolic)                                             
res.mca.common.metabolic.unique = MCA(data.frame.mca.metabolic.unique)
res.mca.common.cell_cycle.unqiue= MCA(data.frame.mca.cell_cycle.unique)
  
fviz_mca_ind(res.mca.common.cellcycle.metabolic)
fviz_mca_ind(res.mca.not.common.cellcycle.metabolic)
fviz_mca_ind(res.mca.common.metabolic.unique)
fviz_mca_ind(res.mca.common.cell_cycle.unqiue)

common.signal.metabolic = intersect(bp.genes.ensembl$metabolic, bp.genes.ensembl$signal_transduction)
not.common.signal.metabolic = setdiff(union(bp.genes.ensembl$metabolic, 
                                            bp.genes.ensembl$signal_transduction), 
                                      common.signal.metabolic)
signal.unique = setdiff(bp.genes.ensembl$signal_transduction, common.signal.metabolic)
  
data.frame.mca.not.common.signal.metabolic = 
  data.cancer.gene[, match(intersect(not.common.signal.metabolic, colnames(data.cancer.gene)),
                           colnames(data.cancer.gene))]
data.frame.mca.common.signal.metabolic = 
  data.cancer.gene[, match(intersect(common.signal.metabolic, colnames(data.cancer.gene)),
                           colnames(data.cancer.gene))]
data.frame.mca.signal.unique = 
  data.cancer.gene[, match(intersect(signal.unique, colnames(data.cancer.gene)),
                           colnames(data.cancer.gene))][1:550]

res.mca.not.common.signal.metabolic = MCA(data.frame.mca.not.common.signal.metabolic)
res.mca.common.signal.metabolic = MCA(data.frame.mca.common.signal.metabolic)
res.mca.signal.unique = MCA(data.frame.mca.signal.unique)

fviz_mca_ind(res.mca.not.common.signal.metabolic)
fviz_mca_ind(res.mca.common.signal.metabolic)
fviz_mca_ind(res.mca.signal.unique)

data.frame.mca.common.all = data.cancer.gene[, match(intersect(setdiff(bp.genes.ensembl$cell_cycle,int3), colnames(data.cancer.gene)),
                                                              colnames(data.cancer.gene))]

###random sequences####
res.mca.common.all = MCA(data.frame.mca.common.all)
fviz_mca_ind(res.mca.common.all)
indexes.rand = list()
indexes.rand.100 = list()
indexes.rand.500 = list()
for(i in seq(10))
{
  i = 1
  ind = sort(sample(seq(length(colnames(data.cancer.gene))), 500))
  data.frame.mca.trial = data.cancer.gene[, ind]
  res.mca.trial = MCA(data.frame.mca.trial, graph = F)
  indexes.rand.500[[i]] = ind
  fviz_mca_ind(res.mca.trial)
  dev.copy(png, paste(toString(i),'.png'))
  dev.off()
  i = i + 1
}
#cell.cycle.kegg = read.delim('Human_cellcycle_kegg.txt')
write(cell.cycle.kegg$GeneID, 'cell_cycle_kegg_ensg.txt', ncolumns = 1)
cell.cycle.ens= read.csv('comparing/cell_cycle_ensembl.csv')
cell.cycle.ens = cell.cycle.ens$Ensembl.Gene.ID
cell.cycle.ens =  as.character(cell.cycle.ens)
intersect(cell.cycle.ens, colnames(data.cancer.gene))
data.frame.mca.cell.cycle = data.cancer.gene[,
  sort(match(intersect(cell.cycle.ens, colnames(data.cancer.gene)), colnames(data.cancer.gene)))]
res.mca.cell.cycle = MCA(data.frame.mca.cell.cycle)
fviz_mca_ind(res.mca.cell.cycle)
hcpc.cell.cycle = HCPC(res.mca.cell.cycle)
plot(hcpc.cell.cycle, choice = 'tree')

indexes.cell.cycle = match(intersect(cell.cycle.ens,colnames(data.cancer.gene)), colnames(data.cancer.gene))
indexes.all.cycle = sort(union(indexes.match.same, indexes.cell.cycle))
res.mca.all.cycle = MCA(data.cancer.gene[, indexes.all.cycle[c(1:258,c(260:409),c(411:506),c(508:639), 
                        c(641:667), c(669:683),c(684:696), c(698:4213), c(4215:4251), 
                        c(4253,5834))]], graph = F)
fviz_mca_ind(res.mca.all.cycle)
hcpc.all.cycle = HCPC(res.mca.all.cycle, graph = F)
plot(hcpc.all.cycle, choice = 'tree')
dev.off()
gene.sep = c(259,410,507,640, 667, 697, 4214, 4252)#684 separates but not that much

total.genes = colnames(data.cancer.gene)

c1metabolism.genes = read.delim('c1metabolism.txt', header = F)
c1metabolism.genes =  as.character(c1metabolism.genes$V1)
c1metabolism.genes.indexes = match(intersect(c1metabolism.genes, total.genes), total.genes)
res.mca.c1metabolism = MCA(data.cancer.gene[,c1metabolism.genes.indexes])
fviz_mca_ind(res.mca.c1metabolism)

tca.genes = read.delim('TCA.txt', header = F)
tca.genes = as.character(tca.genes$V1)
tca.genes.indexes = match(intersect(tca.genes, total.genes), total.genes)
res.mca.tca = MCA(data.cancer.gene[,tca.genes.indexes])
fviz_mca_ind(res.mca.tca)

fattyacidelongation.genes = c('ENSG00000118402', 	'ENSG00000119915', 'ENSG00000164181',
                              'ENSG00000066322','ENSG00000170522')
fattyacidelongation.genes.indexes = match(intersect(fattyacidelongation.genes, total.genes), total.genes)
res.mca.fattyacidelongation = MCA(data.cancer.gene[,fattyacidelongation.genes.indexes])
fviz_mca_ind(res.mca.fattyacidelongation)
