##########Working with DEGS########################
source('preprocess.R')
common.degs = Reduce(intersect, diff.expressed.all.canc.corr) #genes which are degs in all cancers
degs.int.bps = list()
l.bps = length(bp.genes.ensembl)
l.tumors = length(tissues.names)
for(i in seq(l.bps))
{
  degs.int.bps[[i]] = sapply(diff.expressed.all.canc.corr, function(y)
    {
    intersect(bp.genes.ensembl[[i]],y)
  })
  names(degs.int.bps[[i]]) = names(diff.expressed.all.canc.corr) 
}
names(degs.int.bps) = names(bp.genes.ensembl)
lengths.bps.degs = list()
for(i in seq(l.bps))
  lengths.bps.degs[[i]] = find.lengths(degs.int.bps[[i]])
names(lengths.bps.degs) = names(degs.int.bps)

df.lengths.bps.degs = create.df.length(lengths.bps.degs, l.bps, l.tumors, short.names)
library(ggplot2)
ggplot.needs(df.lengths.bps.degs, colnames(df.lengths.bps.degs)[1], colnames(df.lengths.bps.degs)[2] , 'blue')
ggplot.needs(df.lengths.bps.degs, colnames(df.lengths.bps.degs)[1], colnames(df.lengths.bps.degs)[3] , 'red')
ggplot.needs(df.lengths.bps.degs, colnames(df.lengths.bps.degs)[1], colnames(df.lengths.bps.degs)[4] , 'green')
ggplot.needs(df.lengths.bps.degs, colnames(df.lengths.bps.degs)[1], colnames(df.lengths.bps.degs)[5] , 'yellow')
ggplot.needs(df.lengths.bps.degs, colnames(df.lengths.bps.degs)[1], colnames(df.lengths.bps.degs)[6] , 'pink')
ggplot.needs(df.lengths.bps.degs, colnames(df.lengths.bps.degs)[1], colnames(df.lengths.bps.degs)[7] , 'orange')
ggplot.needs(df.lengths.bps.degs, colnames(df.lengths.bps.degs)[1], colnames(df.lengths.bps.degs)[8] , 'purple')

common.bps.degs <- sapply(degs.int.bps, function(x)
  {
    Reduce(intersect, x)
})

unique.degs.bps <- list()
for(i in seq(l.bps))
{
 unique.degs.bps[[i]] = lapply(seq(l.tumors), function(x)
 {
                                   Reduce(setdiff, c(degs.int.bps[[i]][x], 
                                          degs.int.bps[[i]][-x]))
 })
 names(unique.degs.bps[[i]]) = tissues.names
} 
names(unique.degs.bps) = names(bp.genes.ensembl)

lengths.uni.degs.bp = list()
for(i in seq(l.bps))
    lengths.uni.degs.bp[[i]] = find.lengths(unique.degs.bps[[i]])
names(lengths.uni.degs.bp) = names(unique.degs.bps)


df.lengths.unique.bps.degs = create.df.length(lengths.uni.degs.bp, l.bps, l.tumors, short.names)
ggplot.needs(df.lengths.unique.bps.degs, colnames(df.lengths.unique.bps.degs)[1], colnames(df.lengths.unique.bps.degs)[2] , 'blue')
ggplot.needs(df.lengths.unique.bps.degs, colnames(df.lengths.unique.bps.degs)[1], colnames(df.lengths.unique.bps.degs)[3] , 'red')
ggplot.needs(df.lengths.unique.bps.degs, colnames(df.lengths.unique.bps.degs)[1], colnames(df.lengths.unique.bps.degs)[4] , 'green')
ggplot.needs(df.lengths.unique.bps.degs, colnames(df.lengths.unique.bps.degs)[1], colnames(df.lengths.unique.bps.degs)[5] , 'yellow')
ggplot.needs(df.lengths.unique.bps.degs, colnames(df.lengths.unique.bps.degs)[1], colnames(df.lengths.unique.bps.degs)[6] , 'pink')
ggplot.needs(df.lengths.unique.bps.degs, colnames(df.lengths.unique.bps.degs)[1], colnames(df.lengths.unique.bps.degs)[7] , 'orange')
ggplot.needs(df.lengths.unique.bps.degs, colnames(df.lengths.unique.bps.degs)[1], colnames(df.lengths.unique.bps.degs)[8] , 'purple')

