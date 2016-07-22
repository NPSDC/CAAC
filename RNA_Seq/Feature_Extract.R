library(TCGAbiolinks)
setwd('~/Downloads/data/mRNA_counts/')
mRNA_files = list.files()[!grepl('txt',list.files())] ##all files except manifest
library(rjson)
json_data <- fromJSON(file = '../meta_mRNA_counts.json') 
pat.ids <- sapply(seq(length(mRNA_files)), function(x)
{
  unlist(json_data[[x]][[14]])[3] ##14 column corrseponds to associated entities which contain patient id
})
length(unique(pat.ids))
file.ids <- sapply(seq(length(mRNA_files)), function(x)
{
  substr(json_data[[x]][4],1, nchar(json_data[[2]][4])-3)})
length(unique(file.ids))

sample.file.map = data.frame(pat.ids = pat.ids, file.ids = file.ids)

files.to.be.copied = c()
for(i in mRNA_files)
{
  R.utils::gunzip(paste(i,list.files(i)[1],sep = '/'))
  files.to.be.copied = c(files.to.be.copied, paste(i,list.files(i)[1], sep='/'))
}
file.copy(files.to.be.copied,'.')
unlink(mRNA_files, recursive = T) ##Deleting the directories
##So far unzipped the expression values and copied them to previous folder and deleted the original directories

##Reading expression value of each file for each gene checking whehther the genes match
##Then store them in a proper frame exp_prof
mRNA_files = list.files(pattern = 'counts')
file1 = read.delim(mRNA_files[1], header = F)
ens.ids.all = as.character(file1$V1) ##all ensembl ids including transcript
g = sapply(ens.ids.all, function(x) 
{
  unlist(strsplit(x, split = '.', fixed = T))[1]
}) ##removing the symbols after .

##For removing ids that are greater that are not entrez ids
library(biomaRt)
listMarts() 
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)
listFilters(ensembl) 
listAttributes(ensembl)
genes.entrez = getBM(attributes = c('ensembl_gene_id', 'entrezgene'), filters = 'ensembl_gene_id', values = g, mart = ensembl)
View(genes.entrez)
ens.ids = genes.entrez$ensembl_gene_id[(!is.na(genes.entrez$entrezgene))] #ens ids that are not na
##We get all the ensembl ids that have entrez genes
length(unique(ens.ids)) != length(ens.ids)
##
ens.ids = unique(ens.ids)
##Testing
##Making sure that for ens ids we have entrez ids
for(i in ens.ids)
{
  ens.ind = match(i,genes.entrez$ensembl_gene_id)
  if(is.na(genes.entrez$entrezgene[ens.ind]))
  {
    print('lut gaye')
    break
  }
}
##Making sure for all genes except ens ids there is not entrez ids
for(i in setdiff(genes.entrez$ensembl_gene_id, ens.ids))
{
  ens.ind = match(i,genes.entrez$ensembl_gene_id)
  if(!is.na(genes.entrez$entrezgene[ens.ind]))
  {
    print('lut gaye')
    break
  }
}
sum(is.na(ens.ids)) == 0

indexes.ens = match(ens.ids, g) ##Mapping ens ids in our original gene list
##Test
j = 1
for(i in indexes.ens)
{
  if(unlist(strsplit(ens.ids.all[i],'.', fixed = T))[1] != ens.ids[j])
  {
    print('lut gaye')
    break
  } 
  j = j +1
}

##Building a count matrix
ens.ids = ens.ids.all[indexes.ens] ##All ens ids which have transcripts
exp_prof = file1$V2[indexes.ens]
for(i in c(2:length(mRNA_files)))
{
  file = read.delim(mRNA_files[i], header = F)
  ids = as.character(file$V1)
  if(sum(is.na(match(ids, ens.ids.all))) != 0 | length(ids)!=length(ens.ids.all))
  {
    print('bhag')
    break
  }
  if(length(intersect(ids, ens.ids)) != length(ens.ids) | sum(is.na(match(intersect(ids, ens.ids), ens.ids))) != 0)
  {
    print('bhag')
    break
  }
  #print(indexes)
  indexes = match(ens.ids, ids)
  exp_prof = cbind(exp_prof, file$V2[indexes.ens])
  print(i)
}
exp_prof = read.csv('~/Dropbox/honours/RNA_Seq/counts.csv', header = T)
remove(file)
exp_prof = data.frame(exp_prof)
file.names.match = match(mRNA_files, sample.file.map$file.ids)
if(sum(is.na(file.names.match)) != 0)
  print('need to talk')
colnames(exp_prof) = sample.file.map$pat.ids[file.names.match]
rownames(exp_prof) = ens.ids
b = read.delim(mRNA_files[2], header = F)

ens.ids = sort(ens.ids)
ind.sort = match(ens.ids, rownames(exp_prof)) ##sorted rownames
for(i in seq_along(colnames(exp_prof)))
  exp_prof[,i] = exp_prof[,i][ind.sort]
rownames(exp_prof) = ens.ids

##Sorting samples
sample.ids = sort(colnames(exp_prof))
ind.sam.sort = match(sample.ids, colnames(exp_prof))
exp_prof = exp_prof[ind.sam.sort]
exp_prof = as.data.frame(exp_prof)
labels = sapply(sample.ids,function(x) ##Normal or Tumor 
{
  unlist(strsplit(x,split = '.', fixed = T))[4]
})
group.samples = sapply(labels, function(x)
  {
  if(x == '11A' | x == '11B')
    x = 'N'
  else
    x = "T"
})
pat.ids = sapply(sample.ids,function(x)
{
  m = unlist(strsplit(x,split = '.', fixed = T))
  paste(m[1],m[2],m[3], sep = '-')
})
normal.indexes = which(group.samples == 'N')
tumor.indexes = which(group.samples == 'T')
diff.ids = c() ##ids of samples to be taken for finding degs

for(i in normal.indexes)
{
  for(j in tumor.indexes)
  {
    if(pat.ids[i] == pat.ids[j])
      diff.ids = c(diff.ids,i,j)
  }
}
diff.ids = unique(diff.ids)
setwd('~/Dropbox/honours/RNA_Seq')
write.csv(exp_prof, 'counts.csv')
