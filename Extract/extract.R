setwd("~/honours/Extract/")
canc = read.csv('../cancer.csv') #The original cancer data
normal = read.csv('../normal_tissue.csv') # The normal tissue data
unique_genes = unique(canc$Gene) #The unique genes in cancer
mapping = read.delim('map-swiss.tab', header = TRUE, sep = "\t") #Protein mappings from their ids
cancer.liver = canc[canc$Tumor == 'liver cancer',] #Contains only liver cancer data
normal.liver = normal[normal$Tissue == 'liver', ] #Contains only liver cells
hepatocytes = normal.liver[which(normal.liver$Cell.type == "hepatocytes"),] #Contains only hepatocytes
conc.hep = levels(hepatocytes$Level) # Various modes of conc such as high;low,med;not det
hepatocytes$Level = sapply(hepatocytes$Level, function(x)
{
  #Sets the hepatocytes to only high or low
  if(x == conc.hep[3] | x == conc.hep[1])
    x = conc.hep[1]
  else
    x = conc.hep[2]
  #return x;
}
)

calculate.level <- function(threshold, count.high, count.low)
{
  #Sets the level of a gene based on no of threshold
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
mod.canc.liver = data.frame(Gene = rep(unique(cancer.liver$Gene),2),  #Contains the modified levels of liver cancer
                            Level = factor(rep(c("High", "Low"), l)), #and their levels w.r.t proportions
                            count = rep(0L, 2*l), total = rep(0L, 2*l))
                            
i = 0
conc.cancer = levels(cancer.liver$Level)
unique.cancer.genes = as.character(unique(cancer.liver$Gene))
high.indexes = cancer.liver$Level == conc.cancer[1] 
med.indexes = cancer.liver$Level == conc.cancer[3]  #Medium
low.indexes = cancer.liver$Level == conc.cancer[2]  
not.det.indexes = cancer.liver$Level == conc.cancer[4] #Not Detected

for (x in unique.cancer.genes) #Sets the modified cancer liver with thresholds and counts
  {
    genes.indexes = cancer.liver$Gene == x
    mod.canc.liver$Gene[c(2*i +1 , 2*i + 2)] <- x
    mod.canc.liver$Level[c(2*i +1 , 2*i + 2)] <- c("High", "Low")
    mod.canc.liver$count[c(2*i +1 , 2*i + 2)] <- c(cancer.liver$Count.patients[which(genes.indexes & high.indexes)]
                                            + cancer.liver$Count.patients[which(genes.indexes & med.indexes)],
                                            cancer.liver$Count.patients[which(genes.indexes & low.indexes)] +
                                              cancer.liver$Count.patients[which(genes.indexes & not.det.indexes)] )
    mod.canc.liver$total[c(2*i +1 , 2*i + 2)] <- mod.canc.liver$count[2*i + 1] + mod.canc.liver$count[2*i + 2]
    levels = sapply(c(0.5, 0.75, 0.9, 1), function(x)
    {
      calculate.level(x, mod_canc_liver$count[c(2*i + 1)], mod.canc.liver$count[c(2*i + 2)])
    })
    mod.canc.liver$actual.level.50[c(2*i +1 , 2*i + 2)] <- levels[1]
    mod.canc.liver$actual.level.75[c(2*i +1 , 2*i + 2)] <- levels[2]
    mod.canc.liver$actual.level.90[c(2*i +1 , 2*i + 2)] <- levels[3]
    mod.canc.liver$actual.level.100[c(2*i +1 , 2*i + 2)] <- levels[4]
    i = i + 1
}

#Missing genes in hepatocyte
hep.missing <- as.character(unique.cancer.genes[which(is.na(match(unique.cancer.genes, hepatocytes$Gene) == T))])

#Genes common in hepatocyte and liver cancer
hep.present <- as.character(unique.cancer.genes[match(hepatocytes$Gene, unique.cancer.genes)])
indexes <- sapply(hep.present ,function(x) #Indexess of common genes in mod.canc.liver
  {
   x <- match(x , mod.canc.liver$Gene)
  
})
l.hep = length(hepatocytes$Gene) #Total number of genes in hepatocytes
hepatocytes$Canc.Level.50 = mod.canc.liver$actual.level.50[indexes] #Copying data(of proportions) 
hepatocytes$Canc.Level.75 = mod.canc.liver$actual.level.75[indexes] #from mod.canc.liver to hepatocytes
hepatocytes$Canc.Level.90 =  mod.canc.liver$actual.level.90[indexes]
hepatocytes$Canc.Level.100 = mod.canc.liver$actual.level.100[indexes]
 

find.mismatch <- function(hep.prop.level, no.of.unique.genes) #Genes differentially expressed w.r.t proportion
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
names(ans) <- c("50", "75", "90", "100")

write.csv(cancer.liver, "LiverCancer.csv")
write.csv(normal.liver, "LiverNormal.csv")
write.csv(mod.canc.liver, "LiverCancerUpdated.csv")
write.csv(hepatocytes, "LiverNormalWithCancerLevel.csv")
