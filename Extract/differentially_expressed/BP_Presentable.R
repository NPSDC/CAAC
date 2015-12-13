convert.bpterm <- function(bp.term)
{
  l = length(bp.term)
  bp = data.frame(t(data.frame(strsplit(bp.term, '\t'))))
  colnames(bp) <- c("Category","Term", "Count", "Percentage", "Pvalue", "Genes", "List.Total","Pop.Hits", 
                    "Pop.Total", "Fold.Enrichment", "Bonferroni", "Benjamini", "FDR")
  rownames(bp) <- seq_len(l)
  return (bp)
  }

write.bp.csv <- function(bp.csv)
{
  temp = read.delim(bp.csv, sep = '')
  te <- sapply(temp, as.character)
  bp <- convert.bpterm(te)
  write.csv(bp, bp.csv)
}

write.bp.list.csv <- function(list.csv)
{
  null <- sapply(list.csv, write.bp.csv)
}
  
csv.files = list.files(pattern = '.csv')
write.bp.list.csv(csv.files)
