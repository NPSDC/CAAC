####Uses uniquely_exp from tissuesfrom normal and even1 and bp.genes.ensembl from exact file dont remember,
####total.genes from mca_cancer
###Work with both uniquely exp and uniquely diff
#######Renal####################

#contains the intersection of uniquely expressed degs in renal 50% and bps


int.renal.50.unique.bps = sapply(bp.genes.ensembl, function(x)
  {
    intersect(x, total.genes[uniquely.exp.canc.corr$Renal.Cancer])
})
names(int.renal.50.unique.bps) = names(bp.genes.ensembl)
unique.renal.genes.names = mapping[
  match(intersect(total.genes[uniquely.exp.canc.corr$Renal.Cancer], mapping$yourlist.M20150920158KD5HLKJ)
        ,mapping$yourlist.M20150920158KD5HLKJ), c(1,6)]
unique.renal.up = names(which(data.cancer.corr[15,uniquely.exp.canc.corr$Renal.Cancer] == 'Present'))
unique.renal.down = setdiff(total.genes[uniquely.exp.canc.corr$Renal.Cancer], unique.renal.up)

unique.renal.diff.down = intersect(uniquely.diff.corr$Renal.Cancer,
                                  cancer.all.gene.wise.corr$Gene[which(cancer.all.gene.wise.corr$Renal.Cancer == 'Not detected')])
unique.renal.diff.up = setdiff(uniquely.diff.corr$Renal.Cancer, unique.renal.diff.down)

int.renal.50.unique.bps.up = sapply(int.renal.50.unique.bps, function(x)
  {
  intersect(x, unique.renal.up)
})
names(int.renal.50.unique.bps.up) = names(bp.genes.ensembl)

int.renal.50.unique.bps.down = sapply(int.renal.50.unique.bps, function(x)
{
  intersect(x, unique.renal.down)
})
names(int.renal.50.unique.bps.down) = names(bp.genes.ensembl)
int.renal.50.diff.down = sapply(bp.genes.ensembl, function(x)
  {
  intersect(x, unique.renal.diff.down)
})
int.renal.50.diff.up = sapply(bp.genes.ensembl, function(x)
{
  intersect(x, unique.renal.diff.up)
})

write(unique.renal.up, 'up.txt')
met.prot = as.character(unique.renal.genes.names$Gene.names[match(int.renal.50.unique.bps$metabolic, unique.renal.genes.names$yourlist.M20150920158KD5HLKJ)])
write.lists(int.renal.50.unique.bps.down, 2)
write.lists(int.renal.50.unique.bps.up, 2)
setwd('../down')
write.lists(int.renal.50.diff.down, 2)
setwd('../up/')
write.lists(int.renal.50.diff.up, 2)

###Lymphoma
exp.lymp.down = intersect(diff.expressed.down.all.canc.corr$Lymphoma, total.genes[uniquely.exp.canc.corr$Lymphoma])
exp.lymp.up = intersect(diff.expressed.up.all.canc.corr$Lymphoma, total.genes[uniquely.exp.canc.corr$Lymphoma])
diff.lymph.down = intersect(diff.expressed.down.all.canc.corr$Lymphoma, uniquely.diff.corr$Lymphoma)
diff.lymph.up = intersect(diff.expressed.up.all.canc.corr$Lymphoma, uniquely.diff.corr$Lymphoma)
int.lymp.50.exp.bps =  sapply(bp.genes.ensembl, function(x)
{
  intersect(x, total.genes[uniquely.exp.canc.corr$Lymphoma])
})
int.lymp.50.exp.bps.down = sapply(bp.genes.ensembl, function(x)
{
  intersect(x, exp.lymp.down)
})
int.lymp.50.exp.bps.up = sapply(bp.genes.ensembl, function(x)
{
  intersect(x, exp.lymp.up)
})
int.lymp.50.diff.bps.up = sapply(bp.genes.ensembl, function(x)
  {
  intersect(x, diff.lymph.up)
})
int.lymp.50.diff.bps.down = sapply(bp.genes.ensembl, function(x)
{
  intersect(x, diff.lymph.down)
})
write(exp.lymp.up, 'exp.txt')
setwd('../down/')
write(exp.lymp.down, 'exp.txt')
setwd('../../diff/up/')
write(diff.lymph.up, 'diff.txt')
setwd('../down/')
write(diff.lymph.down,'diff_down.txt')
setwd('../../exp/down')
write.lists(int.lymp.50.exp.bps.down, 2)
write.lists(int.lymp.50.exp.bps.up, 2)
setwd('../../diff/up/')
write.lists(int.lymp.50.diff.bps.up, 2)
setwd('../down')
write.lists(int.lymp.50.diff.bps.down, 2)

library(ggplot2)
int.all.50.diff.bps.down= lapply(bp.genes.ensembl, function(x)
  {
     lapply(uniquely.diff.corr.down, function(y)
     {
       intersect(x,y)
     })
})

int.all.50.diff.bps.up= lapply(bp.genes.ensembl, function(x)
{
  lapply(uniquely.diff.corr.up, function(y)
  {
    intersect(x,y)
  })
})

int.all.50.exp.bps.up= lapply(bp.genes.ensembl, function(x)
{
  lapply(uniquely.exp.canc.corr.up, function(y)
  {
    intersect(x,y)
  })
})

int.all.50.exp.bps.down= lapply(bp.genes.ensembl, function(x)
{
  lapply(uniquely.exp.canc.corr.down, function(y)
  {
    intersect(x,y)
  })
})
diff.df.bps = list()
for(i in seq(length(bp.genes.ensembl)))
{
  diff.df.bps[[i]] = data.frame(name = short.names, exp.up =  find.lengths(int.all.50.exp.bps.up[[i]]),
                         exp.down = find.lengths(int.all.50.exp.bps.down[[i]]), 
                         diff.up = find.lengths(int.all.50.diff.bps.up[[i]]),
                         diff.down = find.lengths(int.all.50.diff.bps.down[[i]]))
}
names(diff.df.bps) = names(bp.genes.ensembl)
##Apoptosis
ggplot.needs(diff.df.bps[[1]], colnames(diff.df.bps[[1]])[1], colnames(diff.df.bps[[1]])[2] , 'blue')
ggplot.needs(diff.df.bps[[1]], colnames(diff.df.bps[[1]])[1], colnames(diff.df.bps[[1]])[3] , 'red')
ggplot.needs(diff.df.bps[[1]], colnames(diff.df.bps[[1]])[1], colnames(diff.df.bps[[1]])[4] , 'green')
ggplot.needs(diff.df.bps[[1]], colnames(diff.df.bps[[1]])[1], colnames(diff.df.bps[[1]])[5] , 'yellow')

##Cell Adhesion
ggplot.needs(diff.df.bps[[2]], colnames(diff.df.bps[[2]])[1], colnames(diff.df.bps[[2]])[2] , 'blue')
ggplot.needs(diff.df.bps[[2]], colnames(diff.df.bps[[2]])[1], colnames(diff.df.bps[[2]])[3] , 'red')
ggplot.needs(diff.df.bps[[2]], colnames(diff.df.bps[[2]])[1], colnames(diff.df.bps[[2]])[4] , 'green')
ggplot.needs(diff.df.bps[[2]], colnames(diff.df.bps[[2]])[1], colnames(diff.df.bps[[2]])[5] , 'yellow')


##Cell Cycle
ggplot.needs(diff.df.bps[[3]], colnames(diff.df.bps[[3]])[1], colnames(diff.df.bps[[3]])[2] , 'blue')
ggplot.needs(diff.df.bps[[3]], colnames(diff.df.bps[[3]])[1], colnames(diff.df.bps[[3]])[3] , 'red')
ggplot.needs(diff.df.bps[[3]], colnames(diff.df.bps[[3]])[1], colnames(diff.df.bps[[3]])[4] , 'green')
ggplot.needs(diff.df.bps[[3]], colnames(diff.df.bps[[3]])[1], colnames(diff.df.bps[[3]])[5] , 'yellow')

##Immune
ggplot.needs(diff.df.bps[[4]], colnames(diff.df.bps[[4]])[1], colnames(diff.df.bps[[4]])[2] , 'blue')
ggplot.needs(diff.df.bps[[4]], colnames(diff.df.bps[[4]])[1], colnames(diff.df.bps[[4]])[3] , 'red')
ggplot.needs(diff.df.bps[[4]], colnames(diff.df.bps[[4]])[1], colnames(diff.df.bps[[4]])[4] , 'green')
ggplot.needs(diff.df.bps[[4]], colnames(diff.df.bps[[4]])[1], colnames(diff.df.bps[[4]])[5] , 'yellow')

#Metabolic
ggplot.needs(diff.df.bps[[5]], colnames(diff.df.bps[[5]])[1], colnames(diff.df.bps[[5]])[2] , 'blue')
ggplot.needs(diff.df.bps[[5]], colnames(diff.df.bps[[5]])[1], colnames(diff.df.bps[[5]])[3] , 'red')
ggplot.needs(diff.df.bps[[5]], colnames(diff.df.bps[[5]])[1], colnames(diff.df.bps[[5]])[4] , 'green')
ggplot.needs(diff.df.bps[[5]], colnames(diff.df.bps[[5]])[1], colnames(diff.df.bps[[5]])[5] , 'yellow')

##Signal
ggplot.needs(diff.df.bps[[6]], colnames(diff.df.bps[[6]])[1], colnames(diff.df.bps[[6]])[2] , 'blue')
ggplot.needs(diff.df.bps[[6]], colnames(diff.df.bps[[6]])[1], colnames(diff.df.bps[[6]])[3] , 'red')
ggplot.needs(diff.df.bps[[6]], colnames(diff.df.bps[[6]])[1], colnames(diff.df.bps[[6]])[4] , 'green')
ggplot.needs(diff.df.bps[[6]], colnames(diff.df.bps[[6]])[1], colnames(diff.df.bps[[6]])[5] , 'yellow')

##Transport
ggplot.needs(diff.df.bps[[7]], colnames(diff.df.bps[[7]])[1], colnames(diff.df.bps[[7]])[2] , 'blue')
ggplot.needs(diff.df.bps[[7]], colnames(diff.df.bps[[7]])[1], colnames(diff.df.bps[[7]])[3] , 'red')
ggplot.needs(diff.df.bps[[7]], colnames(diff.df.bps[[7]])[1], colnames(diff.df.bps[[7]])[4] , 'green')
ggplot.needs(diff.df.bps[[7]], colnames(diff.df.bps[[7]])[1], colnames(diff.df.bps[[7]])[5] , 'yellow')


