setwd("~/honours/Extract/")
canc = read.csv('../cancer.csv')
normal = read.csv('../normal_tissue.csv')
unique_genes = unique(canc$Gene)
mapping = read.delim('map-swiss.tab', header = TRUE, sep = "\t")
cancer.liver = canc[canc$Tumor == 'liver cancer',]
normal.liver = normal[normal$Tissue == 'liver', ]
hepatocytes = normal.liver[which(normal.liver$Cell.type == "hepatocytes"),]
conc_hep = levels(hepatocytes$Level)
hepatocytes$Level = sapply(hepatocytes$Level, function(x)
{
  if(x == conc_hep[3] | x == conc_hep[1])
    x = conc_hep[1]
  else
    x = conc_hep[2]
  #return x;
}
)

calculate.level <- function(threshold, count.high, count.low)
{
  total <- count.high + count.low
  highest <- max(count.high, count.low)
  if(highest/total >= threshold)
  {
    if(highest == count.high)
      return('High')
    else
      return('Low')
  }   
  else
    return(NA)
}

l = length(unique(cancer.liver$Gene)) #Length of unique genes
mod_canc_liver = data.frame(Gene = rep(unique(cancer.liver$Gene),2),
                            Level = factor(rep(c("High", "Low"), l)),
                            count = rep(0L, 2*l), total = rep(0L, 2*l),
                            actual_level = factor(rep(c("High", "Low"), l)))
i = 0
conc = levels(cancer.liver$Level)
unique.cancer.genes = as.character(unique(cancer.liver$Gene))
high.indexes = cancer.liver$Level == conc[1]
med.indexes = cancer.liver$Level == conc[3]  #Medium
low.indexes = cancer.liver$Level == conc[2]  
not.det.indexes = cancer.liver$Level == conc[4] #Not Detected
for (x in unique.cancer.genes)
  {
    genes.indexes = cancer.liver$Gene == x
    mod.canc.liver$Gene[c(2*i +1 , 2*i + 2)] <- x
    mod.canc.liver$Level[c(2*i +1 , 2*i + 2)] <- c("High", "Low")
    mod.canc.liver$count[c(2*i +1 , 2*i + 2)] <- c(cancer.liver$Count.patients[which(genes.indexes & high.indexes)]
                                            + cancer.liver$Count.patients[which(genes.indexes & med.indexes)],
                                            cancer.liver$Count.patients[which(genes.indexes & low.indexes)] +
                                              cancer.liver$Count.patients[which(genes.indexes & not.det.indexes)] )
    mod.canc.liver$total[c(2*i +1 , 2*i + 2)] <- cancer.liver$Total.patients[which(genes.indexes)[1] ]
    levels = sapply(c(0.5, 0.75, 0.9, 1), function(x)
    {
      calculate.level(x, mod_canc_liver$count[c(2*i + 1)], mod.canc.liver$count[c(2*i + 2)])
    })
    #mod_canc_liver$paste0('actual.level.', c(50, 75, 90, 100))[c(2*i +1 , 2*i + 2)] = levels
    mod.canc.liver$actual_level_50[c(2*i +1 , 2*i + 2)] <- levels[1]
    mod.canc.liver$actual_level_75[c(2*i +1 , 2*i + 2)] <- levels[2]
    mod.canc.liver$actual_level_90[c(2*i +1 , 2*i + 2)] <- levels[3]
    mod.canc.liver$actual_level_100[c(2*i +1 , 2*i + 2)] <- levels[4]
    i = i + 1
}
hep.missing <- as.character(unique.cancer.genes[which(is.na(match(unique.cancer.genes, hepatocytes$Gene) == T))])
hep.matching <- as.character(unique.cancer.genes[match(hepatocytes$Gene, unique.cancer.genes)])
indexes <- sapply(hep.matching ,function(x)
  {
   x <- match(x , mod_canc_liver$Gene)
  
})
l.hep = length(hepatocytes$Gene)
hepatocytes$Canc.Level.50 = mod_canc_liver$actual_level_50[indexes]
hepatocytes$Canc.Level.75 = mod_canc_liver$actual_level_75[indexes]
hepatocytes$Canc.Level.90 =  mod_canc_liver$actual_level_90[indexes]
hepatocytes$Canc.Level.100 = mod_canc_liver$actual_level_100[indexes]
 

find.mismatch <- function(hep.prop.level, no.of.unique.genes)
{
  ans = c()
  for (id in 1:no.of.unique.genes)
  {
    
    if(!is.na(hep.prop.level[id]))
    {
      if(hepatocytes$Level[id] != hep.prop.level[id])
      {
        ans = c(ans, as.character(hepatocytes$Gene[id]))
      }
    }
    
  }
  return(ans)
}


ans = sapply(list(hepatocytes$Canc.Level.50, hepatocytes$Canc.Level.75, hepatocytes$Canc.Level.90, 
       hepatocytes$Canc.Level.100), function(x)
       {
          find.mismatch(x, l.hep)
         })



