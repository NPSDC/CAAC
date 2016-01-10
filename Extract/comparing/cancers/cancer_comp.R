source('extract.R')
prepocess.cancer <- function(cancer.name)
{
  if(is.na(match(cancer.name, unique(canc$Tumor))))
  {
    print("wrong cancer name")
    return()
  }
  canc.type = canc[canc$Tumor == cancer.name, ]
  canc.type$Level = set.level(renal.cancer$Level, levels(canc.type$Level), 2)
  return(canc.type)
}

initialize.empty <- function(unique.genes, levels, times.to.repeat, no.unique.genes)
  ###Initializes the data frame creating an empty data frame specific for cancer as well 
  ###function reduce.levels
  ###Inputs are self explanatory
{
  no.unique.genes = length(unique.genes)
  canc.specific = data.frame(Gene = rep(unique.genes, each = times.to.repeat),  
                             Level = factor(rep(levels, no.unique.genes)), 
                             Count = rep(0L, no.unique.genes*times.to.repeat),
                             Total = rep(0L, no.unique.genes*times.to.repeat))  
  return(canc.specific)
}

reduce.levels <- function(canc.name, levels, combine.conc)
  ###Reudces the 4 levels of cancer to 2 levels
  ###canc.name- name cancer
  ###levels - levels of conc to keep(character vector)
  ###combine.conc - list concentration levels to be combined in each level(list)
{
  canc.type <- prepocess.cancer(canc.name)
  l = length(levels)
  unique.genes = unique(canc.type$Gene)
  no.unique.genes = length(unique(canc.type$Gene))
  canc.specific = initialize.empty(unique.genes, levels, l, no.unique.genes)
  if(length(levels) != length(combine.conc))
  {
    print("Check Params")
    return()
  }
  
  indexes.levels = list() ##contains the indexes of the various levels
  for(i in seq(length(combine.conc)))
  {
    indexes.levels[[i]] = list()
    for(j in seq(length(combine.conc[[i]])))
        indexes.levels[[i]][[j]] = which(canc.type$Level == combine.conc[[i]][j])
  }
  k = 1
  for(gene in unique.genes)
  {
    gene.index = which(gene == canc.type$Gene) #gene index of the gene in the cancer set to be reduced
    gene.and.level.index = c()
    for(i in seq(l))
    {
      for(j in seq(length(indexes.levels[[i]])))
      {
        gene.and.level.index = c(gene.and.level.index,
                                 gene.index[!is.na(match(gene.index, indexes.levels[[i]][[j]]))])
      }
      canc.specific$Gene[k + i - 1] = gene
      canc.specific$Count[k + i - 1] = sum(canc.type$Count[gene.and.level.index])
      gene.and.level.index = c()
    }
    canc.specific$Total[k] = canc.specific$Total[k + 1] = 
      canc.specific$Count[k] + canc.specific$Count[k+1]
    
    levels.vals = sapply(c(0.5, 0.75, 0.9, 1), function(x)
    {
      calculate.level(x,canc.specific$Count[k], canc.specific$Count[k+1], levels)
    })
    canc.specific$actual.level.50[c(k, k + 1)] <- levels.vals[1]
    canc.specific$actual.level.75[c(k, k + 1)] <- levels.vals[2]
    canc.specific$actual.level.90[c(k, k + 1)] <- levels.vals[3]
    canc.specific$actual.level.100[c(k, k + 1)] <- levels.vals[4]
    k = k + l
      
  }
  return(canc.specific)
}
renal.mod.canc = reduce.levels('renal cancer', c('Present', 'Not detected'), list('Present', 'Not detected'))
