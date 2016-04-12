source('preprocess.R')

canc.present.1 = add.vals(canc, 1, 2, 7)

canc$only1.present = rep('Present', length(canc$Gene))
for(i in seq(1,length(canc$Gene), 4))
{
  if(canc[i+3, 4] == canc[i+3,5])
  {
    canc[i,7] = 'Not detected'
    canc[i+1,7] = 'Not detected'
    canc[i+2,7] = 'Not detected'
    canc[i+3,7] = 'Not detected'
  }
}

library(FactoMineR)
library(factoextra)

data.present.1 = t(canc.present.1)
colnames(data.present.1) = data.present.1[1, ]
data.present.1 = data.present.1[-1, ]
res.mca.all.1.present = MCA(data.present.1)
fviz_mca_ind(res.mca.all.1.present)

g.indexes = c()
for(i in seq(1,length(canc[,1]), 4))
  g.indexes = c(g.indexes, i, i+1)
#tissues.canc ##using directly here otherwise could be read from cancers in all_tissues

tissues.canc.present1 = lapply(tissues.canc, function(x)
  {
  add.1.present(x, 3, 4)
})
names(tissues.canc.present1) = tissues.names

diff.expressed.all.canc.present1 = mapply(function(x,y)
{
  find.diff.expressed(x,y,c(4,9),c(1,1))
}, tissues.normal, tissues.canc)
names(diff.expressed.all.canc.present1) = tissues.names

diff.expressed.up.all.canc.present1 = mapply( function(x,y)
{
  find.up.regulated(x,9,1,y)
},diff.expressed.all.canc.present1,tissues.canc)
names(diff.expressed.up.all.canc.present1) = tissues.names

unique.exp.present1 <- sapply(tumor.indexes, function(x)
  {
  find.unique.cancer(data.present.1, x, total.genes )
})
names(unique.exp.present1) = tissues.names
lengths.unqiue.exp.present1 = sapply(unique.exp, length)

diff.expressed.down.all.canc.present1 = find.down.regulated(diff.expressed.all.canc.present1, diff.expressed.up.all.canc.present1, tissues.names)
write.lists(diff.expressed.down.all.canc.present1, 2)
write.lists(diff.expressed.up.all.canc.present1, 2)
write.lists(uniquely.diff.present1, 2)
write.lists(uniquely.exp.canc.present1, 2)

###uniquely differentiated
uniquely.diff.present1 <- mapply(function(x,y)
{
  find.unique(x, diff.expressed.all.canc.present1, y)
}, diff.expressed.all.canc.present1, seq(1,17))
names(uniquely.diff.present1) = tissues.names
lengths.uniquely.diff.present1 <- sapply(uniquely.diff.present1, length)

###uniquely expressed in cancer
uniquely.exp.canc.present1 <- mapply(function(x,y)
{
  find.unique.cancer(data.present.1, x, y)
}, tumor.indexes, diff.expressed.all.canc.present1)
names(uniquely.exp.canc.present1) = tissues.names
lengths.uniquely.exp.canc.present1 <- sapply(uniquely.exp.canc.present1, length)

library(ggplot2)
short.names = sapply(tissues.names, function(x)
  strsplit(x, '\\.')[[1]][1])
short.names = sapply(short.names, function(x)
  substr(x,1,3))
diff.df.1 = create.data.frame(short.names, diff.expressed.all.canc.present1, diff.expressed.down.all.canc.present1, diff.expressed.up.all.canc.present1)
diff.df.1$uniquely_diff = lengths.uniquely.diff.present1
diff.df.1$uniquely_exp = lengths.uniquely.exp.canc.present1

ggplot.needs(diff.df.1, colnames(diff.df.1)[1], colnames(diff.df.1)[2] , 'blue')
ggplot.needs(diff.df.1, colnames(diff.df.1)[1], colnames(diff.df.1)[3] , 'red')
ggplot.needs(diff.df.1, colnames(diff.df.1)[1], colnames(diff.df.1)[4] , 'green')
ggplot.needs(diff.df.1, colnames(diff.df.1)[1], colnames(diff.df.1)[5] , 'pink')
ggplot.needs(diff.df.1, colnames(diff.df.1)[1], colnames(diff.df.1)[6] , 'orange')


#####Corrected###############
canc.present.1.corrected = canc.present.1
canc.present.1.corrected = recorrect.all(tissues.canc.corr, canc.present.1.corrected, 9, tumor.indexes+1)

diff.expressed.all.canc.present1.corr = mapply(function(x,y)
{
  find.diff.expressed(x,y,c(4,9),c(1,1))
}, tissues.normal, tissues.canc.corr)
names(diff.expressed.all.canc.present1.corr) = tissues.names

diff.expressed.up.all.canc.present1.corr = mapply( function(x,y)
{
  find.up.regulated(x,9,1,y)
},diff.expressed.all.canc.present1.corr,tissues.canc.corr)
names(diff.expressed.up.all.canc.present1.corr) = tissues.names


diff.expressed.down.all.canc.present1.corr = find.down.regulated(diff.expressed.all.canc.present1.corr, diff.expressed.up.all.canc.present1.corr, tissues.names)
write.lists(diff.expressed.down.all.canc.present1.corr, 2)
write.lists(diff.expressed.up.all.canc.present1.corr, 2)
write.lists(uniquely.diff.present1.corr, 2)
write.lists(uniquely.exp.canc.present1.corr, 3)

###uniquely differentiated
uniquely.diff.present1.corr <- mapply(function(x,y)
{
  find.unique(x, diff.expressed.all.canc.present1.corr, y)
}, diff.expressed.all.canc.present1.corr, seq(1,17))
names(uniquely.diff.present1.corr) = tissues.names
lengths.uniquely.diff.present1.corr <- sapply(uniquely.diff.present1.corr, length)
names(lengths.uniquely.diff.present1.corr) = tissues.names

###uniquely expressed in cancer
data.present.1.corr = t(canc.present.1.corrected)
colnames(data.present.1.corr) = data.present.1.corr[1,]
data.present.1.corr = data.present.1.corr[-1, ]
data.present.1.corr = data.frame(data.present.1.corr)

uniquely.exp.canc.present1.corr <- mapply(function(x,y)
{
  find.unique.cancer(data.present.1.corr, x, y)
}, tumor.indexes, diff.expressed.all.canc.present1.corr)
names(uniquely.exp.canc.present1.corr) = tissues.names
lengths.uniquely.exp.canc.present1.corr <- sapply(uniquely.exp.canc.present1.corr, length)
names(lengths.uniquely.exp.canc.present1.corr) = tissues.names

uniquely.exp.present1.corr.all <- sapply(tumor.indexes, function(x)
{
  find.unique.cancer(data.present.1.corr, x, total.genes )
})
names(uniquely.exp.present1.corr.all) = tissues.names
lengths.uniquely.exp.present1.corr.all = sapply(unique.exp.corr, length)
names(lengths.uniquely.exp.present1.corr.all) = tissues.names
###MCA ANALYSIS##################
library(FactoMineR)
library(factoextra)
library(HCPC)
res.mca.1.corr = MCA(data.present.1.corr)
fviz_mca_ind(res.mca.1.corr)
hcc.1.corr = HCPC(res.mca.1.corr)
plot(hcc.1.corr, choice = 'tree' )

library(ggplot2)
short.names = sapply(tissues.names, function(x)
  strsplit(x, '\\.')[[1]][1])
short.names = sapply(short.names, function(x)
  substr(x,1,3))
diff.df.1.corr = create.data.frame(short.names, diff.expressed.all.canc.present1.corr, diff.expressed.down.all.canc.present1.corr, diff.expressed.up.all.canc.present1.corr)
diff.df.1.corr$uniquely_diff = lengths.uniquely.diff.present1.corr
diff.df.1.corr$uniquely_exp = lengths.uniquely.exp.canc.present1.corr
diff.df.1.corr$uniquely_exp.all = lengths.uniquely.exp.present1.corr.all
ggplot.needs(diff.df.1.corr, colnames(diff.df.1.corr)[1], colnames(diff.df.1.corr)[2] , 'blue')
ggplot.needs(diff.df.1.corr, colnames(diff.df.1.corr)[1], colnames(diff.df.1.corr)[3] , 'red')
ggplot.needs(diff.df.1.corr, colnames(diff.df.1.corr)[1], colnames(diff.df.1.corr)[4] , 'green')
ggplot.needs(diff.df.1.corr, colnames(diff.df.1.corr)[1], colnames(diff.df.1.corr)[5] , 'pink')
ggplot.needs(diff.df.1.corr, colnames(diff.df.1.corr)[1], colnames(diff.df.1.corr)[6] , 'orange')
ggplot.needs(diff.df.1.corr, colnames(diff.df.1.corr)[1], colnames(diff.df.1.corr)[7] , 'yellow')
