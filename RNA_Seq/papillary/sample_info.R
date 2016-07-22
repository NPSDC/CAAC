##Using feature_extract.R
##We use exp_prof from Feature_Extract.R or load it from our cancer directory
##Note that we had 32 normal samples
##But for 1 normal sample we dont have matched tumor its id is 182
labels = sapply(sample.ids,function(x) ##Normal or Tumor (01 or 11) 
{
  unlist(strsplit(x,split = '-', fixed = T))[4]
})
sum(is.na(labels)) == 0
type = sapply(labels, function(x)
{
  if(x == '11A' | x == '11B')
    x = 'N'
  else
    x = "T"
})
sum(is.na(type)) == 0
pats = sapply(sample.ids,function(x) ##Patients
{
  m = unlist(strsplit(x,split = '-', fixed = T))
  paste(m[1],m[2],m[3], sep = '-')
})

length(unique(pats))
normal.indexes = which(type == 'N')
tumor.indexes = which(type == 'T')
diff.ids = c() ##ids of samples to be taken for finding degs

for(i in normal.indexes)
{
  for(j in tumor.indexes)
  {
    if(pats[i] == pats[j])
      diff.ids = c(diff.ids,i,j)
  }
}
diff.ids = unique(diff.ids)
length(intersect(diff.ids, tumor.indexes))
length(diff.ids)


sample.ids = colnames(exp_prof)
sum(sample.ids == sort(sample.ids)) == length(sample.ids)
sample.info = data.frame(sample.names = sample.ids[diff.ids], type = type[diff.ids],
                           pat.ids = pats[diff.ids], pat.nums = rep(1, length(diff.ids)), stringsAsFactors = F)
sample.info$pat.comb = mapply(function(x,y){
  paste(as.character(x),y,sep = '')},
  sample.info$pat.nums, sample.info$type)


exp_prof_diff = exp_prof[,diff.ids] ##Only for tumor-match control

j = 1
pat = sample.info$pat.ids[1]
for(i in seq_along(diff.ids))
{
  if(sample.info$pat.ids[i] != pat)
  {
    j = j + 1
    pat = sample.info$pat.ids[i]
  }
  sample.info$pat.nums[i] = j
}

sample.info$type = as.factor(sample.info$type)  ##Needed for design matrix
sample.info$pat.nums = as.factor(sample.info$pat.nums)
