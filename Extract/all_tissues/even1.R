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

tissues.canc = lapply(tissues.canc, function(x)
  {
  add.1.present(x, 3, 4)
})
names(tissues.canc) = tissues.names

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

diff.expressed.down.all.canc.present1 = find.down.regulated(diff.expressed.all.canc.present1, diff.expressed.up.all.canc.present1, tissues.names)
write.lists(diff.expressed.down.all.canc.present1, 2)
write.lists(diff.expressed.up.all.canc.present1, 2)

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

