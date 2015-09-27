canc = read.csv('../cancer.csv')
normal = read.csv('../normal_tissue.csv')
unique_genes = unique(canc$Gene)
setwd("~/honours/Extract/")
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

calculate.max <- function(threshold, count.high, count.low)
{
  total <- count.high + count.low
  highest <- max(count.high, count.low)
  if(highest/total > threshold)
  {
    if(highest == count.high)
      return('high')
    else
      return('low')
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
med.indexes = cancer.liver$Level == conc[3]
low.indexes = cancer.liver$Level == conc[2]
not_det.indexes = cancer.liver$Level == conc[4]
for (x in unique_genes)
  {
    genes.indexes = cancer.liver$Gene == x
    mod_canc_liver$Gene[c(2*i +1 , 2*i + 2)] <- x
    mod_canc_liver$Level[c(2*i +1 , 2*i + 2)] <- c("High", "Low")
    mod_canc_liver$count[c(2*i +1 , 2*i + 2)] <- c(cancer.liver$Count.patients[which(genes.indexes & high.indexes)]
                                            + cancer.liver$Count.patients[which(genes.indexes & med.indexes)],
                                            cancer.liver$Count.patients[which(genes.indexes & low.indexes)] +
                                              cancer.liver$Count.patients[which(genes.indexes & not_det.indexes)] )
    mod_canc_liver$total[c(2*i +1 , 2*i + 2)] <- cancer.liver$Total.patients[which(genes.indexes)[1] ]
    mod_canc_liver$actual_level_50[c(2*i +1 , 2*i + 2)] <- mod_canc_liver$Level[which.max(mod_canc_liver$count[c(2*i +1 , 2*i + 2)])]
    mod_canc_liver$actual_level_75[c(2*i +1 , 2*i + 2)] <- mod_canc_liver$Level[which.max(mod_canc_liver$count[c(2*i +1 , 2*i + 2)])]
    mod_canc_liver$actual_level_90[c(2*i +1 , 2*i + 2)] <- mod_canc_liver$Level[which.max(mod_canc_liver$count[c(2*i +1 , 2*i + 2)])]
    mod_canc_liver$actual_level_100[c(2*i +1 , 2*i + 2)] <- mod_canc_liver$Level[which.max(mod_canc_liver$count[c(2*i +1 , 2*i + 2)])]
    i = i + 1
    
}
hep.missing <- as.character(unique.cancer.genes[which(is.na(match(unique.cancer.genes, hepatocytes$Gene) == T))])
hep.matching <- as.character(unique.cancer.genes[match(hepatocytes$Gene, unique.cancer.genes)])
indexes <- sapply(hep.matching ,function(x)
  {
   x <- which(x == mod_canc_liver$Gene)[1]
  
})

hepatocytes$Canc_Level = mod_canc_liver$actual_level[indexes ]
hep.not.matching = as.character(hepatocytes$Gene[which(hepatocytes$Level != hepatocytes$Canc_Level)])

