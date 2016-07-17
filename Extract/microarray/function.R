filter.genes <- function(df, basis.col, gene.col)  
{
  rem.index = c()
  indexes.viewd = c()
  df[,gene.col] = as.character(df[,gene.col])
  for(i in seq(length(df[,gene.col])))
  {
    if(df[,gene.col][i] == '')
    {
      rem.index = c(rem.index, i)
      next
    }
    if(i %in% indexes.viewd)
      next
    
    indexes = which(df[,gene.col][i] == df[,gene.col])
    indexes.viewd = c(indexes.viewd, indexes)
    if(sum(df[,basis.col][indexes] >= 0) != length(indexes) && sum(df[,basis.col][indexes] <= 0) != length(indexes))
    {
      print(i)
      next
    }
    #print(indexes)
    #if(df[,basis.col])
    # print(indexes)
    max.index = indexes[which.max(abs(df[,basis.col][indexes]))]
    #print(max.index)
    #print(req.index)
    req.index = setdiff(indexes, max.index)
    rem.index = c(rem.index, setdiff(req.index, rem.index))
  }
  rem.index = sort(rem.index)
  #print(rem.index)
  return(temp[-rem.index,])
}

add.gene.name <- function(df, gene.col, sym.col, map.df)
{
  prev.id = df[,gene.col][1]
  prev.name = 'TSPAN6'
  df[,gene.col] = as.character(df[,gene.col])
  map.df
  for(i in seq(length(df[,gene.col])))
  {
    #print(i)
    if(df[,gene.col][i] == prev.id)
      df[,sym.col][i] = prev.name
    else
    {
      gene.sym = map.genes.ids(df[,gene.col][i], map.df)
      if(length(gene.sym) == 0)
        next
      df[,sym.col][i] = gene.sym
    }
    prev.id = df[,gene.col][i]
    prev.name = df[,sym.col][i]
  }
  return(df[,sym.col])
}

find.sub.int <- function(df, gene.col, int.genes)
{
  df[,gene.col] = as.character(df[,gene.col])
  req.indexes = c()
  for(i in seq(length(df[,gene.col])))
  {
    if(df[,gene.col][i] %in% int.genes)
      req.indexes = c(req.indexes, i)
  }
  print(req.indexes)
  return(df[req.indexes,])
}