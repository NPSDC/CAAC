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
