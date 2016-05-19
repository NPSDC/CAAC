####Deals with code leading to the creation of data set for doing MCA
tissues.mod.indexes = c(4,5,7,8,9,11,12,19,21,22,23,25,26,27,30,31,32,42,43,44,45,47)

mod.indexes = c()
for(i in tissues.mod.indexes)
{
  indexes =  which(as.character(normal$Tissue) == levels(normal$Tissue)[i])
  mod.indexes = c(mod.indexes, indexes)
}
sort.mod.index = sort(mod.indexes)
mod.normal = normal[sort.mod.index,]
write.csv( mod.normal, 'all_tissues/mca/normal/mod_normal.csv')
df.normal = data.frame(Gene = total.genes, bre1 = rep(c('P'), length(total.genes)),
                       bre2 = rep(c('P'), length(total.genes)), bre3 = rep(c('P'), length(total.genes)),
                       cerv1 = rep(c('P'), length(total.genes)), cerv2 = rep(c('P'), length(total.genes)),
                       col1 = rep(c('P'), length(total.genes)), col2 = rep(c('P'), length(total.genes)),
                       col3 = rep(c('P'), length(total.genes)),col4 = rep(c('P'), length(total.genes)),
                       end1 = rep(c('P'), length(total.genes)), end2 = rep(c('P'), length(total.genes)),
                       end3 = rep(c('P'), length(total.genes)), end4 = rep(c('P'), length(total.genes)),
                       gli1 = rep(c('P'), length(total.genes)), gli2 = rep(c('P'), length(total.genes)),
                       gli3 = rep(c('P'), length(total.genes)), gli4 = rep(c('P'), length(total.genes)),
                       head1 = rep(c('P'), length(total.genes)), head2 = rep(c('P'), length(total.genes)), 
                       liver1 = rep(c('P'), length(total.genes)),liver2 = rep(c('P'), length(total.genes)),
                       lung1 = rep(c('P'), length(total.genes)), lung2 = rep(c('P'), length(total.genes)), 
                       lung3 = rep(c('P'), length(total.genes)), lymph1 = rep(c('P'), length(total.genes)),
                       lymph2 = rep(c('P'), length(total.genes)), ovary1 = rep(c('P'), length(total.genes)),
                       ovary2 = rep(c('P'), length(total.genes)), panc1 = rep(c('P'), length(total.genes)),
                       panc2 = rep(c('P'), length(total.genes)), pros = rep(c('P'), length(total.genes)),
                       renal1 = rep(c('P'), length(total.genes)), renal2 = rep(c('P'), length(total.genes)),
                       stom1 = rep(c('P'), length(total.genes)), stom2 = rep(c('P'), length(total.genes)),
                       test1 = rep(c('P'), length(total.genes)), test2 =rep(c('P'), length(total.genes)),
                       thyr = rep(c('P'), length(total.genes)),  uroth = rep(c('P'), length(total.genes)),
                       stringsAsFactors = F
                       )
for(i in seq(length(tissues.normal.copy)))
{
  tissues.normal.copy[[i]] = mod.cell.type(tissues.normal.copy[[i]], 2, 3)
  tissues.normal.copy[[i]][,3] = as.factor(tissues.normal.copy[[i]][,3])
}
mod.cell.type <- function(tissue, tissue.col, cell.type.col)
{
  tissue[,tissue.col] = as.character(tissue[,tissue.col])
  tissue[,cell.type.col] = as.character(tissue[,cell.type.col])
  for(i in seq(length(tissue[,cell.type.col])))
  {
    tissue[,cell.type.col][i] = paste(substr(tissue[,tissue.col][i],1,3), tissue[,cell.type.col][i], sep = '.')
  }
  return(tissue)
}
mod.cell.type.rect <- function(tissue, tissue.col, cell.type.col)
{
  tissue[,tissue.col] = as.character(tissue[,tissue.col])
  tissue[,cell.type.col] = as.character(tissue[,cell.type.col])
  for(i in seq(length(tissue[,cell.type.col])))
  {
    l = nchar(tissue[,tissue.col][i])
    tissue[,cell.type.col][i] = paste(substr(tissue[,tissue.col][i],1,3),substr(tissue[,tissue.col][i],l,l), tissue[,cell.type.col][i], sep = '.')
  }
  tissue[,cell.type.col] = as.factor(tissue[,cell.type.col])
  return(tissue)
}
tissues.normal.copy$Endometrial.Cancer = mod.cell.type.rect(tissues.normal.copy$Endometrial.Cancer,2,3)
tissues.normal.copy$Stomach.Cancer = mod.cell.type.rect(tissues.normal.copy$Stomach.Cancer,2,3)
###Rectifying endometrial

gen.cols.normal <- function(df, df.index, tissue, tissue.level, tissue.cell.type)
{
  print(length(levels(tissue[,tissue.cell.type])))
  for(j in seq(length(levels(tissue[,tissue.cell.type]))))
  {
    levels.df = c()
    cell.type.indexes = which(as.character(tissue[,tissue.cell.type])== 
                                as.character(levels(tissue[,tissue.cell.type]))[j])
    print(cell.type.indexes)
    for(i in seq(length(total.genes)))
    {
      gene.indexes = which(tissue[,1] == total.genes[i])
      int.index = intersect(gene.indexes, cell.type.indexes)
      #print(length(int.index))
      if(length(int.index) == 0)
        levels.df = c(levels.df, 'NA')
      else if(length(int.index) == 1)
        levels.df = c(levels.df, as.character(tissue[,tissue.level][int.index]))
      else
        print('fu')
      
    }
    print(levels.df)
    df[,df.index[j]] = levels.df
  }
  return(df)
}

df.normal = gen.cols.normal(df.normal, c(2,3,4), tissues.normal.copy$Breast.Cancer, 4, 3)
df.normal = gen.cols.normal(df.normal, c(5,6), tissues.normal.copy$Cervical.Cancer, 4, 3)
df.normal = gen.cols.normal(df.normal, c(7,8,9,10), tissues.normal.copy$Colorectal.Cancer, 4, 3)
df.normal = gen.cols.normal(df.normal, c(11,12,13,14), tissues.normal.copy$Endometrial.Cancer, 4, 3)
df.normal = gen.cols.normal(df.normal, c(15,16,17,18), tissues.normal.copy$Glioma, 4, 3)
df.normal = gen.cols.normal(df.normal, c(19,20), tissues.normal.copy$Head.And.Neck.Cancer, 4, 3)
df.normal = gen.cols.normal(df.normal, c(21,22), tissues.normal.copy$Liver.Cancer, 4, 3)
df.normal = gen.cols.normal(df.normal, c(23,24,25), tissues.normal.copy$Lung.Cancer, 4, 3)
df.normal = gen.cols.normal(df.normal, c(26,27), tissues.normal.copy$Lymphoma, 4, 3)
df.normal = gen.cols.normal(df.normal, c(28,29), tissues.normal.copy$Ovarian.Cancer, 4, 3)
df.normal = gen.cols.normal(df.normal, c(30,31), tissues.normal.copy$Pancreatic.Cancer, 4, 3)
df.normal = gen.cols.normal(df.normal, c(32), tissues.normal.copy$Prostate.Cancer, 4, 3)
df.normal = gen.cols.normal(df.normal, c(33,34), tissues.normal.copy$Renal.Cancer, 4, 3)
df.normal = gen.cols.normal(df.normal, c(35,36), tissues.normal.copy$Stomach.Cancer, 4, 3)
df.normal = gen.cols.normal(df.normal, c(37,38), tissues.normal.copy$Testis.Cancer, 4, 3)
df.normal = gen.cols.normal(df.normal, c(39), tissues.normal.copy$Thyroid.Cancer, 4, 3)
df.normal = gen.cols.normal(df.normal, c(40), tissues.normal.copy$Urothelial.Cancer, 4, 3)


###MCA####
library(FactoMineR)
library(factoextra)
genes.na = c()
for(i in seq(length(total.genes)))
{
  indexes = which(df.normal[i,] == 'NA')
  if(length(indexes) > 0)
    df.normal[i,indexes] = NA
}

df.mca.normal = t(df.normal)
colnames(df.mca.normal) = df.mca.normal[1,]
df.mca.normal = df.mca.normal[-1, ]
df.mca.normal = data.frame(df.mca.normal)
res.mca.normal = MCA(df.mca.normal, na.method = 'NA')
fviz_mca_ind(res.mca.normal, axes = c(1,2))
df.mca.normal.rem = df.mca.normal[c(-8,-27,-28),]
res.mca.normal.rem = MCA(df.mca.normal.rem, na.method = 'NA')
fviz_mca_ind(res.mca.normal.rem)

combined = rbind(df.mca.normal, data.cancer.corr)
res.mca.combined = MCA(combined, na.method = 'NA')
fviz_mca_ind(res.mca.combined, axes = c(2,3))

combined
res.mca.combined1 = MCA(combined1, na.method = 'NA')
fviz_mca_ind(res.mca.combined1, axes = c(1,2))
combined2 = combined1[-8,]
res.mca.combined2 = MCA(combined2, na.method = 'NA')
fviz_mca_ind(res.mca.combined2, axes = c(1,2))

diff.expressed.atleast.4 <- get.diff.expressed.atleast.4(total.genes, diff.expressed.all.canc.corr)
get.diff.expressed.atleast.4 <- function(total.genes, diff.exp)
{
  diff.genes = c()
  for(i in total.genes)
  {
    count = 0
    for(j in seq(length(diff.exp)))
    {
      if(i %in% diff.exp[[j]])
        count = count + 1
      if(count == 4)
      {
        diff.genes = c(diff.genes, i)
        break
      }
    }
  }
  return(diff.genes)
}
  
###MCA on the genes that are differenitally expressed across atleast 4 cancers###
res.mca.atleast.4 <- MCA(data.cancer.corr[,match(diff.expressed.atleast.4, colnames(data.cancer.corr))])
fviz_mca_ind(res.mca.atleast.4, axes = c(2,3))
clus.atleast.4 = HCPC(res.mca.atleast.4, nb.clust = -1)
res.mca.met <- MCA(data.cancer.corr[,match(intersect(total.genes, bp.genes.ensembl$metabolic), colnames(data.cancer.corr))])
fviz_mca_ind(res.mca.met)
clus <- HCPC(res.mca.met, nb.clust = -1)

res.mca.combined.atleast.4 = MCA(combined[,match(diff.expressed.atleast.4, colnames(combined))])
fviz_mca_ind(res.mca.combined.atleast.4)
