########DEC 15##############################
###Contains the  data of all 6 patients along with all genes
hepatocytes = read.csv('hepatocytes.csv')
hcc.data.complete = read.csv('~/Dropbox/honours/Extract/hcc_data/msb145122-sup-0008-Dataset5/Dataset 5.csv', stringsAsFactors = F)
indexes.match.hep =  match(hcc.data.complete$GENES, hepatocytes$Gene) #Contains the match of indexes of hcc in hepatocytes
hcc.data.complete = transform(hcc.data.complete, HPA.50 = as.character(hepatocytes$Canc.Level.50[indexes.match.hep]))

levels(hcc.data.complete$HPA.50) = c(levels(hcc.data.complete$HPA.50), 'None')
hcc.data.complete$HPA.50[hcc.data.complete$HPA.50 == 'Not detected'] = 'None'
match.among = list() #contains match among hpa data and patient data
match.among = sapply(hcc.data.complete$GENES, function(x)
  {})
for(i in seq(length(hcc.data.complete$GENES)))
{
  vals = as.character(hcc.data.complete[i,])
  ans = c()
  row.length = length(hcc.data.complete)
  if(is.na(vals[row.length]))
    ans = c(NA,NA,NA,NA,NA,NA)
  else
  {  
    for(j in seq(2,row.length - 1))
    {
      if(as.character(vals[j]) == as.character(vals[row.length]))
               ans = c(ans, colnames(hcc.data.complete)[j])
    }
  }
  match.among[[i]] = ans
}
names(match.among) = hcc.data.complete$GENES

source('extract.R')
  