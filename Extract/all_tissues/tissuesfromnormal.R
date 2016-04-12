source('preprocess.R')

tissue.indexes = list(c(4), c(8), c(9,31), c(11,12), c(7), c(25,32), c(21), c(5,22), c(23), c(26),
                      c(27), c(30), c(19), c(42,43), c(44), c(45), c(47))
tissue.names = list(c('breast'), c('cervix, uterine'), c('colon','rectum'), 
                    c('endometrium 1', 'endometrium 2'), c('cerbral cortex'), 
                    c('oral mucosa', 'salivary gland'), c('liver'), c('bronchus', 'lung'), 
                    c('lymph node'), c('ovary'), c('pancreas'), c('prostate'), c('kidney'), 
                    c('stomach 1', 'stomach 2'), c('testis'), c('thyroid gland'), c('urinary bladder'))
tissue.names = unique(normal$Tissue)
colorectal.tissue = preprocess(tissue.indexes[[3]], normal, 2, 4, tissue.names)

tumor.indexes = c(c(1),c(3:10),c(12:15),c(17:20))
tumor.names = as.character(unique(canc$Tumor))
tissues.names = unlist(sapply(tumor.names[tumor.indexes], alter.name))

tissues.normal = lapply(tissue.indexes, function(x)
  {
    preprocess(x, normal, 2, 4, tissue.names)  
})
names(tissues.normal) = tissues.names

tissues.canc = lapply(tumor.indexes, function(x)
  {
  reduce.levels(x, c('Present', 'Not detected'), list('Present', 'Not detected'), tumor.names)
})
names(tissues.canc) = tissues.names

diff.expressed.all.canc = mapply(function(x,y)
  {
    find.diff.expressed(x,y,c(4,5),c(1,1))
}, tissues.normal, tissues.canc)
names(diff.expressed.all.canc) = tissues.names

diff.expressed.up.all.canc = mapply( function(x,y)
  {
  find.up.regulated(x,5,1,y)
},diff.expressed.all.canc,tissues.canc)
names(diff.expressed.up.all.canc) = tissues.names

diff.expressed.down.all.canc = find.down.regulated(diff.expressed.all.canc, diff.expressed.up.all.canc, tissues.names)

###writing files
write.lists(tissues.normal, 1)
write.lists(tissues.canc, 1)
write.lists(diff.expressed.down.all.canc, 2)
write.lists(diff.expressed.up.all.canc, 2)
write.lists(uniquely.diff, 2)
write.lists(uniquely.exp.canc, 2)

###uniquely differentiated
uniquely.diff <- mapply(function(x,y)
  {
  find.unique(x, diff.expressed.all.canc, y)
}, diff.expressed.all.canc, seq(1,17))
names(uniquely.diff) = tissues.names
lengths.uniquely.diff <- sapply(uniquely.diff, length)

###uniquely expressed in cancer
uniquely.exp.canc <- mapply(function(x,y)
{
  find.unique.cancer(data.cancer.gene, x, y)
}, tumor.indexes, diff.expressed.all.canc)
names(uniquely.exp.canc) = tissues.names
lengths.uniquely.exp.canc <- sapply(uniquely.exp.canc, length)

library(ggplot2)
short.names = sapply(tissues.names, function(x)
  strsplit(x, '\\.')[[1]][1])
short.names = sapply(short.names, function(x)
  substr(x,1,3))
diff.df.50 = create.data.frame(short.names, diff.expressed.all.canc, diff.expressed.down.all.canc, diff.expressed.up.all.canc)
diff.df.50$uniquely_diff = lengths.uniquely.diff
diff.df.50$uniquely_exp = lengths.uniquely.exp.canc

ggplot.needs(diff.df.50, colnames(diff.df.50)[1], colnames(diff.df.50)[2] , 'blue')
ggplot.needs(diff.df.50, colnames(diff.df.50)[1], colnames(diff.df.50)[3] , 'red')
ggplot.needs(diff.df.50, colnames(diff.df.50)[1], colnames(diff.df.50)[4] , 'green')
ggplot.needs(diff.df.50, colnames(diff.df.50)[1], colnames(diff.df.50)[5] , 'pink')
ggplot.needs(diff.df.50, colnames(diff.df.50)[1], colnames(diff.df.50)[6] , 'orange')



#####Corrected################3
tissues.canc.corr = list() 
for(i in seq(length(tissues.canc.corr)))
{
  tissues.canc.corr[[i]] = recorrect.cancer(tissues.canc.corr[[i]], tissues.normal[[i]], 3, 4, c(5, 4, 9))
}
names(tissues.canc.corr) <- tissues.names

###correcting cancer.all.gene.wise
cancer.all.gene.wise.corr = cancer.all.gene.wise
cancer.all.gene.wise.corr = recorrect.all(tissues.canc.corr, cancer.all.gene.wise.corr, 5, tumor.indexes+1)

diff.expressed.all.canc.corr = mapply(function(x,y)
{
  find.diff.expressed(x,y,c(4,5),c(1,1))
}, tissues.normal, tissues.canc.corr)
names(diff.expressed.all.canc.corr) = tissues.names

diff.expressed.up.all.canc.corr = mapply( function(x,y)
{
  find.up.regulated(x,5,1,y)
},diff.expressed.all.canc.corr,tissues.canc.corr)
names(diff.expressed.up.all.canc.corr) = tissues.names

diff.expressed.down.all.canc.corr = find.down.regulated(diff.expressed.all.canc.corr, diff.expressed.up.all.canc.corr, tissues.names)

###writing files
write.lists(tissues.canc.corr, 1)
write.lists(diff.expressed.down.all.canc.corr, 2)
write.lists(diff.expressed.up.all.canc.corr, 2)
write.lists(uniquely.diff.corr, 2)
write.lists(uniquely.exp.canc.corr, 3)

###uniquely differentiated
uniquely.diff.corr <- mapply(function(x,y)
{
  find.unique(x, diff.expressed.all.canc.corr, y)
}, diff.expressed.all.canc.corr, seq(1,17))
names(uniquely.diff.corr) = tissues.names
lengths.uniquely.diff.corr <- sapply(uniquely.diff.corr, length)
names(lengths.uniquely.diff.corr) <- tissues.names

###uniquely expressed in cancer
data.cancer.corr = t(cancer.all.gene.wise.corr)
#data.cancer.corr = data.frame(data.cancer.corr)
colnames(data.cancer.corr) = data.cancer.corr[1,]
data.cancer.corr = data.cancer.corr[-1,]

uniquely.exp.canc.corr <- mapply(function(x,y)
{
  find.unique.cancer(data.cancer.corr, x, y)
}, tumor.indexes, diff.expressed.all.canc.corr)
names(uniquely.exp.canc.corr) = tissues.names
lengths.uniquely.exp.canc.corr <- sapply(uniquely.exp.canc.corr, length)
names(lengths.uniquely.exp.canc.corr) = tissues.names

uniquely.exp.canc.corr.all <- sapply(tumor.indexes,function(x)
{
  find.unique.cancer(data.cancer.corr, x, total.genes )
})
names(uniquely.exp.canc.corr.all) = tissues.names
lengths.uniquely.exp.canc.corr.all <- sapply(uniquely.exp.canc.corr.all, length)
names(uniquely.exp.canc.corr.all) = tissues.names

uniquely.exp.canc.corr.up = mapply(function(x,y){
  intersect(total.genes[x],y) 
}, uniquely.exp.canc.corr, diff.expressed.up.all.canc.corr)
lengths.uniquely.exp.canc.corr.up <- sapply(uniquely.exp.canc.corr.up, length)
names(lengths.uniquely.exp.canc.corr.up) = tissues.names

uniquely.exp.canc.corr.down = mapply(function(x,y){
  intersect(total.genes[x],y) 
}, uniquely.exp.canc.corr, diff.expressed.down.all.canc.corr)
lengths.uniquely.exp.canc.corr.down <- sapply(uniquely.exp.canc.corr.down, length)
names(lengths.uniquely.exp.canc.corr.down) = tissues.names

uniquely.diff.corr.down = mapply(function(x,y){
  intersect(x,y) 
}, uniquely.diff.corr, diff.expressed.down.all.canc.corr)
lengths.uniquely.diff.corr.down <- sapply(uniquely.diff.corr.down, length)
names(lengths.uniquely.diff.corr.down) = tissues.names

uniquely.diff.corr.up = mapply(function(x,y){
  intersect(x,y) 
}, uniquely.diff.corr, diff.expressed.up.all.canc.corr)
lengths.uniquely.diff.corr.up <- sapply(uniquely.diff.corr.up, length)
names(lengths.uniquely.diff.corr.up) = tissues.names
###MCA ANALYSIS##################
library(FactoMineR)
library(factoextra)
res.mca.corr = MCA(data.cancer.corr)
fviz_mca_ind(res.mca.corr)

###Plotting graphs###
library(ggplot2)
short.names = sapply(tissues.names, function(x)
  strsplit(x, '\\.')[[1]][1])
short.names = sapply(short.names, function(x)
  substr(x,1,3))
diff.df.50.corr = create.data.frame(short.names, diff.expressed.all.canc.corr, diff.expressed.down.all.canc.corr, diff.expressed.up.all.canc.corr)
diff.df.50.corr$uniquely_diff = lengths.uniquely.diff.corr
diff.df.50.corr$uniquely_exp = lengths.uniquely.exp.canc.corr
diff.df.50.corr$uniquely_exp.all = lengths.uniquely.exp.canc.corr.all
diff.df.50.corr$uniquely_exp_down = lengths.uniquely.exp.canc.corr.down
diff.df.50.corr$uniquely_exp_up = lengths.uniquely.exp.canc.corr.up
diff.df.50.corr$uniquely_diff_down = lengths.uniquely.diff.corr.down
diff.df.50.corr$uniquely_diff_up = lengths.uniquely.diff.corr.up

ggplot.needs(diff.df.50.corr, colnames(diff.df.50.corr)[1], colnames(diff.df.50.corr)[2] , 'blue')
ggplot.needs(diff.df.50.corr, colnames(diff.df.50.corr)[1], colnames(diff.df.50.corr)[3] , 'red')
ggplot.needs(diff.df.50.corr, colnames(diff.df.50.corr)[1], colnames(diff.df.50.corr)[4] , 'green')
ggplot.needs(diff.df.50.corr, colnames(diff.df.50.corr)[1], colnames(diff.df.50.corr)[5] , 'pink')
ggplot.needs(diff.df.50.corr, colnames(diff.df.50.corr)[1], colnames(diff.df.50.corr)[6] , 'orange')
ggplot.needs(diff.df.50.corr, colnames(diff.df.50.corr)[1], colnames(diff.df.50.corr)[7] , 'yellow')
ggplot.needs(diff.df.50.corr, colnames(diff.df.50.corr)[1], colnames(diff.df.50.corr)[8] , 'red')
ggplot.needs(diff.df.50.corr, colnames(diff.df.50.corr)[1], colnames(diff.df.50.corr)[9] , 'green')
ggplot.needs(diff.df.50.corr, colnames(diff.df.50.corr)[1], colnames(diff.df.50.corr)[10] , 'red')
ggplot.needs(diff.df.50.corr, colnames(diff.df.50.corr)[1], colnames(diff.df.50.corr)[11] , 'green')
