####only functions/methods

set.level <- function(cell.level, cell.conc, type)
{
  ##Input - Cell whose level has to be set; cell.conc is the types conc, type-tells what to club
  #Output - sets it to high or low
  if(type == 1)
  {  
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
  else
  {
    cell.level = sapply(cell.level, function(x)
    {
      #Sets the Cell$Level to only Present or Not detected
      if(x == cell.conc[4] )
        x = cell.conc[4]
      else
        x = "Present"
      
    })
  }
  return(cell.level)
}

preprocess <- function(name.indexes, data, col.index, level.index, types)
{
  ##It creates a data frame from an exisitng data frame based on name of tissue/tumor
  #name.indexes - Indexes of tissue/tumour
  #data - data.frame from which new data.frame will be created
  #col.index - index of column with which we will compare name
  #level.index - index of level in data.frame
  #types - types of input
  row.indexes = c()
  for(i in seq(length(name.indexes)))
  {
    if(is.na(match(types[name.indexes[i]], types)))
    {
      print("wrong name")
      return()
    }
    row.indexes = c(which(data[,col.index] == types[name.indexes[i]]), row.indexes)
  }
  row.indexes = sort(row.indexes)
  type = data[row.indexes, ]
  type[,level.index] = set.level(type[,level.index], levels(type[,level.index]), 2)
  return(type)
}
reduce.levels <- function(canc.index, levels, combine.conc, types.cancer)
{
  ###Reudces the 4 levels of cancer to 2 levels
  ###canc.index- index of cancer
  ###levels - levels of conc to keep(character vector)
  ###combine.conc - list concentration levels to be combined in each level(list)
  ###types.cancer - different types of cancer
  canc.type <- preprocess(canc.index, canc, 2, 3, types.cancer)
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

find.diff.expressed <- function(normal, cancer, level.indexes, genes.indexes)
{
  genes = intersect(unique(normal[,genes.indexes[1]]), unique(cancer[,genes.indexes[2]]))
  diff.expressed = c()
  for(i in genes)
  {
    canc.level = as.character(cancer[match(i, cancer[,genes.indexes[2]]), level.indexes[2]])
    normal.gene.indexes = which(normal[,genes.indexes[1]] == i)
    normal.level = normal[,level.indexes[1]][normal.gene.indexes][1]
    if((sum(normal[,level.indexes[1]][normal.gene.indexes] == normal.level) == length(normal.gene.indexes))
       && canc.level != normal.level) 
      diff.expressed = c(diff.expressed, i)
  }
  return(diff.expressed)
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

alter.name <- function(name)
{
  name = simpleCap(name)
  tissue = unlist(strsplit(name, ' '))
  
  if(length(tissue) >= 2)
  {
    name = tissue[1]
    for(i in seq(2,length(tissue)))
      name = paste(name, tissue[i], sep = '.')
  }
  
  return(name)
}

find.up.regulated <- function(diff.genes, level.index, genes.index, canc)
{
  ##finds the up regulated genes
  diff.genes.indexes = match(diff.genes, canc[,genes.index])
  return(as.character(canc[diff.genes.indexes, genes.index][as.character(canc[diff.genes.indexes,level.index]) == 'Present']))
}
find.down.regulated <- function(diff.genes, up.genes, names)
{
  down.regulated = mapply(function(x,y)
    {
    setdiff(x,y)
  },diff.genes,up.genes)
  names(down.regulated) = names
  return(down.regulated)
}

add.vals <- function(data,gene.index, tissue.index, level.index) 
  #adds the 2 levels for each gene and each cancer
{
  unique.cancer.genes = unique(data[,gene.index])
  canc.add = data.frame(Gene = unique.cancer.genes, 
                        Breast.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Carcinoid = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Cervical.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Colorectal.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Endometrial.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Glioma = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Head.And.Neck.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Liver.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Lung.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Lymphoma = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Melanoma = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Ovarian.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Pancreatic.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Prostate.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Renal.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Skin.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Stomach.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Testis.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Thyroid.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        Urothelial.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                        stringsAsFactors = F)
  gene.old = data[,gene.index][1]
  j = 1
  for(i in seq(1,length(data[,gene.index]),4))
  {
    print(i)
    gene.new = data[,gene.index][i]
    
    if(gene.new != gene.old)
    {
      gene.old = gene.new
      j = j + 1
    }
    tissue = as.character(data[i,tissue.index])
    tissue = alter.name(tissue)
    #print(colnames(canc.add)[20])
      tissue.col = match(tissue, colnames(canc.add))
      #print(tissue)
      #print(tissue.col)
      canc.add[j,tissue.col] <- data[i, level.index]
      
     # print(canc.add[j,tissue.col])
    }
    return(canc.add)
}

add.1.present <- function(tumor.data, count.index, total.index)
{
  #This function adds the column only.present.1 column from canc to exisiting tissue specific cancer dataframe
  #tumor.col contains the tumor column in data
  #tumor.index index of the tumor in types of tumor
  #tumor.data exisitng tumor data frame corresponding to tumor.data
  #only1.index contains the index corresponding to column index in data
  #types contains the cancer types
  #gene.index contains gene indexes tumor.data
  #g.indexes indexes of data to be considered
  actual_val = c()
  for(i in seq(1, length(tumor.data[,1]), 2))
  {
    if(tumor.data[i+1,count.index] == tumor.data[i+1,total.index])
      actual_val = c(actual_val, 'Not detected', 'Not detected')
    else
      actual_val = c(actual_val, 'Present', 'Present')
  }
  tumor.data$only.present1 = actual_val
  #View(tumor.data)
  return(tumor.data)
}

write.lists <- function(lists, type)
{
  #lists <- a list of items to be written
  #type <- csv or txt
  l = length(lists)
  if(type == 1) #csv
  {
    for(i in seq(l))
      write.csv(lists[[i]],paste(names(lists)[i],'.csv', sep = ''))
      
  }
  if(type == 2) #text
  {
    for(i in seq(l))
      write(lists[[i]],paste(names(lists)[i],'.txt'))
  }
  if(type == 3)#text
  {
    for(i in seq(l))
      write(total.genes[lists[[i]]],paste(names(lists)[i],'.txt'))
  }
}

find.lengths <- function(lists)
{
  lengths <- lapply(lists, length)
  names(lengths) <- names(lists)
  #print(names(lengths))
  return(lengths)
}
create.data.frame <- function(names, all, down, up)
{
  #name <- names of cancers
  #rest differentially expressed genes
  all.lengths <- find.lengths(all)
  down.lengths <- find.lengths(down)
  up.lengths <- find.lengths(up)
  print(names)
  diff.df <- data.frame(cancer = names, all.reg = all.lengths, down.reg = down.lengths, up.reg = up.lengths)
  return(diff.df)
}

find.unique <- function(diff.exp.tissue, diff.exp.all, tissue.index)
{
  #finds the genes which show differential expression for only that cancer
  common = c()
  for(i in seq(length(diff.exp.all)))
  {
    if(i != tissue.index)
      common = union(common, intersect(diff.exp.tissue, diff.exp.all[[i]]))
  }
  only = setdiff(diff.exp.tissue, common)
  return(only)
}

find.unique.cancer <- function(data, tumor.index, diff.expressed.tumor.genes)
{
  total.genes = unique(colnames(data))
  
  diff.indexes = match(intersect(total.genes, diff.expressed.tumor.genes), total.genes)
  #print(length(diff.indexes))
  unique.indexes = c()
  for(i in seq(length(diff.indexes)))
  {
    val = as.character(data[tumor.index, diff.indexes[i]])
    if(sum(data[setdiff(seq(20), tumor.index), diff.indexes[i]] == val) == 0)
      unique.indexes = c(unique.indexes, i)
  }
  return(diff.indexes[unique.indexes])
}

ggplot.needs <- function(df, x, y, color)
{
  ggplot(data = df, aes_string(x,y))  +
    geom_bar(stat = 'identity',position  ="dodge", width=0.5, fill = color) +
    geom_text(aes_string(label = y), vjust = -0.3) +  theme_minimal() 
  
}

recorrect.cancer <- function(tissue.canc, tissue.normal, count.col, total.col, level.indexes)
{
    ##Corrects the cancer for both 50% and only 1 by looking at the normal data
    ##For both I determine the level of normal and then try to set opposite in cancer enabling max diff
    ##level.indexes contains indexes for level of 50%(1) and 1(3) in cancer and (2) for normal
    #print(j)
    half.indexes = which(tissue.canc[,count.col] == tissue.canc[,total.col]/2)
    unique.half.genes = unique(tissue.canc$Gene[half.indexes])
    unique.genes = unique(tissue.canc$Gene)
    #print(half.indexes)
    for(i in unique.half.genes)
    {
      #print(k)
      indexes.normal = which(tissue.normal$Gene == i) 
      indexes.canc = which(tissue.canc$Gene == i)
      l.present.indexes = sum(which(tissue.normal[indexes.normal,level.indexes[2]] == 'Present'))
      l.absent.indexes = sum(which(tissue.normal[indexes.normal,level.indexes[2]] == 'Not detected'))
      #print(l.present.indexes)
      #print(l.absent.indexes)
      if(l.present.indexes > l.absent.indexes)
      {
        #print('yes')
        tissue.canc[indexes.canc,level.indexes[1]] = 'Not detected'  
      }
    }
    for(i in unique.genes)
    {
      indexes.normal = which(tissue.normal$Gene == i) 
      indexes.canc = which(tissue.canc$Gene == i)
      l.present.indexes = sum(which(tissue.normal[indexes.normal,level.indexes[2]] == 'Present'))
      l.absent.indexes = sum(which(tissue.normal[indexes.normal,level.indexes[2]] == 'Not detected'))
      level = 'Present'
      if(l.absent.indexes > l.present.indexes)
      {
        #print('yes')
        level = 'Not detected'
      }
      #print(level)
      if(level == 'Not detected' & tissue.canc[indexes.canc[1], count.col] > 0 )
        tissue.canc[indexes.canc, level.indexes[3]] = 'Present'
      if(level == 'Present' & tissue.canc[indexes.canc[2], count.col] > 0 )
        tissue.canc[indexes.canc, level.indexes[3]] = 'Not detected'
    }
  return(tissue.canc)
}

recorrect.all <- function(tissues.canc, all.canc, col.tissue, col.all)
{
  for(i in seq(length(col.all)))
  {
    tissue.genes = unique(tissues.canc[[i]]$Gene)
    indexes = match(tissue.genes, all.canc$Gene)
    level.indexes = seq(1,2*length(tissue.genes), 2)
    #print(level.indexes)
    #print(indexes)
    #print(length(tissues.canc[level.indexes, col.tissue]))
    all.canc[indexes,col.all[i]] = tissues.canc[[i]][level.indexes, col.tissue]
  }
  return(all.canc)
}

create.df.length <- function(lengths.list, l.bps, l.tumors, short.names)
{
  lengths.df <- data.frame(setNames(replicate(l.bps+1, seq(l.tumors), simplify = F),
                                    c('name', names(lengths.list))))
  for(i in seq(l.bps))
    lengths.df[,i+1] = unlist(lengths.list[[i]])
  lengths.df$name = short.names
  return(lengths.df)
}