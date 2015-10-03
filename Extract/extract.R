#####Functions###########
find.missing <- function(comp1, comp2)
{
  ##Inputs- 2 Gene id sets
  #Outputs - The Gene Ids which are present in comp1 but not comp2
  return (as.character(comp1[which(is.na(match(comp1, comp2) == T))]))
}

set.level <- function(cell.level, cell.conc)
{
  ##Input - Cell whose level has to be set; cell.conc is the types conc
  #Output - sets it to high or low
  cell.level = sapply(cell.level, function(x)
  {
    #Sets the Cell$Level to only high or low
    if(x == cell.conc[3] | x == cell.conc[1])
      x = cell.conc[1]
    else
      x = cell.conc[2]
    #return x;
  }
  )  
}

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

find.mismatch <- function(type.Level,type.Gene, prop.level, no.of.unique.genes) 
{
  #Finds Genes differentially expressed w.r.t proportion
  ans = c() 
  for (id in 1:no.of.unique.genes)
  {
    
    if(!is.na(prop.level[id]))
    {
      if(type.Level[id] != prop.level[id])
      {
        ans = c(ans, as.character(type.Gene[id]))
      }
    }
    
  }
  return(ans)
}

get.mismatch <- function(type.Level, type.Gene, list.cols, l)
{
  ans = sapply(list.cols, function(x)
  {
    find.mismatch(type.Level, type.Gene , x, l)
  })
  names(ans) <- c("50", "75", "90", "100")
  return(ans)
}
###Normal Script
setwd("~/honours/Extract/")
canc = read.csv('../cancer.csv') #The original cancer data
normal = read.csv('../normal_tissue.csv') # The normal tissue data
unique_genes = unique(canc$Gene) #The unique genes in cancer
mapping = read.delim('map-swiss.tab', header = TRUE, sep = "\t") #Protein mappings from their ids
cancer.liver = canc[canc$Tumor == 'liver cancer',] #Contains only liver cancer data
normal.liver = normal[normal$Tissue == 'liver', ] #Contains only liver cells
hepatocytes = normal.liver[which(normal.liver$Cell.type == "hepatocytes"),] #Contains only hepatocytes
bile.duct = normal.liver[which(normal.liver$Cell.type == "bile duct cells"),] #Contains only bile duct cells
conc.hep = levels(hepatocytes$Level) # Various modes of conc such as high;low,med;not det
conc.bile = levels(bile.duct$Level)
hepatocytes$Level = set.level(hepatocytes$Level, conc.hep)
bile.duct$Level = set.level(bile.duct$Level, conc.bile)


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

#Missing genes in hepatocyte, bile.duct
hep.missing <- find.missing(unique.cancer.genes, hepatocytes$Gene)
bile.duct.missing <- find.missing(unique.cancer.genes, bile.duct$Gene)

#Genes common in hepatocyte and liver cancer
hep.present <- as.character(unique.cancer.genes[match(hepatocytes$Gene, unique.cancer.genes)])
bile.present <- as.character(unique.cancer.genes[match(bile.duct$Gene, unique.cancer.genes)])
indexes.hep <- match(hep.present, mod.canc.liver$Gene)
indexes.bile <- match(bile.present, mod.canc.liver$Gene)
l.hep = length(hepatocytes$Gene) #Total number of genes in hepatocytes
l.bile = length(bile.duct$Gene)
hepatocytes$Canc.Level.50 = mod.canc.liver$actual.level.50[indexes.hep] #Copying data(of proportions) 
hepatocytes$Canc.Level.75 = mod.canc.liver$actual.level.75[indexes.hep] #from mod.canc.liver to hepatocytes
hepatocytes$Canc.Level.90 =  mod.canc.liver$actual.level.90[indexes.hep]
hepatocytes$Canc.Level.100 = mod.canc.liver$actual.level.100[indexes.hep]
bile.duct$Canc.Level.50 = mod.canc.liver$actual.level.50[indexes.bile] #Copying data(of proportions) 
bile.duct$Canc.Level.75 = mod.canc.liver$actual.level.75[indexes.bile] #from mod.canc.liver to hepatocytes
bile.duct$Canc.Level.90 =  mod.canc.liver$actual.level.90[indexes.bile]
bile.duct$Canc.Level.100 = mod.canc.liver$actual.level.100[indexes.bile]

ans.hep = get.mismatch(hepatocytes$Level, hepatocytes$Gene, list(hepatocytes$Canc.Level.50, 
                        hepatocytes$Canc.Level.75, hepatocytes$Canc.Level.90, hepatocytes$Canc.Level.100), l.hep)
ans.bile = get.mismatch(bile.duct$Level, bile.duct$Gene,list(bile.duct$Canc.Level.50, bile.duct$Canc.Level.75, bile.duct$Canc.Level.90, 
                             bile.duct$Canc.Level.100), l.bile)



write(ans.hep$`100`, 'Expression100.txt')
write.csv(hepatocytes[match(ans.hep$`100`, hepatocytes$Gene), ], 'hep100.csv')
write.csv(bile.duct[match(ans.bile$`100`, bile.duct$Gene), ], 'bile100.csv')
write.csv(cancer.liver, "LiverCancer.csv")
write.csv(normal.liver, "LiverNormal.csv")
write.csv(mod.canc.liver, "LiverCancerUpdated.csv")
write.csv(hepatocytes, "LiverNormalWithCancerLevel.csv")
