get.imp.genes <- function(method, assay, number)
{
  if(method == 1) ###most expression
    return(order(rowMeans(assay), decreasing = T)[1:number])
  if(method == 2)
  {
    RowVar <- rowSums((assay - rowMeans(assay))^2)/(dim(assay)[2] - 1)
    print(RowVar)
    return(order(RowVar, decreasing = T)[1:number])
  }
}