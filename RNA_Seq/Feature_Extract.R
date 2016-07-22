##Parameters to be given wd, mRNA_file_path, total case files, file_name, proj_dir, cases

wd = '~/Downloads/data/papillary'
setwd(wd)
total_case_files = 321
cases = 289
mRNA_files_path = paste(wd,'mrna',sep = '/')
file_name = 'metadata.cart.mRna.json'
proj_dir = '~/Dropbox/honours/RNA_Seq/papillary/'

mRNA_files = list.files(mRNA_files_path)[!grepl('txt',list.files(mRNA_files_path))] ##all files except manifest
length(mRNA_files) == total_case_files
library(rjson)

json_data <- fromJSON(file = file_name) 
pat.ids <- sapply(seq(length(json_data)), function(x)
{
  unlist(json_data[[x]][[14]])[3] ##14 column corrseponds to associated entities which contain patient id
})
length(unique(pat.ids))
file.ids <- sapply(seq(length(json_data)), function(x)
{
  substr(json_data[[x]][4],1, nchar(json_data[[2]][4])-3)})
length(unique(file.ids)) == total_case_files

sample.file.map = data.frame(pat.ids = pat.ids, file.ids = file.ids)
sample.file.map = sample.file.map[order(sample.file.map$pat.ids),]
View(sample.file.map)
files.to.be.copied = c()
setwd(mRNA_files_path)
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
length(mRNA_files) == total_case_files
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
genes.entrez = genes.entrez[order(genes.entrez$ensembl_gene_id),]

ens.ids = genes.entrez$ensembl_gene_id[(!is.na(genes.entrez$entrezgene))] #ens ids that are not na

##We get all the ensembl ids that have entrez genes
length(unique(ens.ids)) != length(ens.ids)
##
ens.ids = unique(ens.ids)
length(ens.ids)
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
ens.ids = ens.ids.all[indexes.ens] ##All ens ids which have gene transcripts
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
remove(file)


exp_prof = data.frame(exp_prof)
file.names.match = match(mRNA_files, sample.file.map$file.ids)
if(sum(is.na(file.names.match)) != 0)
  print('need to talk')
colnames(exp_prof) = sample.file.map$pat.ids[file.names.match]
rownames(exp_prof) = ens.ids

sum(ens.ids == sort(ens.ids)) == length(ens.ids)
ind.sort = match(ens.ids, rownames(exp_prof)) ##sorted rownames
for(i in seq_along(colnames(exp_prof)))
  exp_prof[,i] = exp_prof[,i][ind.sort]
rownames(exp_prof) = ens.ids

##Sorting samples
sample.ids = sort(colnames(exp_prof))

ind.sam.sort = match(sample.ids, colnames(exp_prof))
exp_prof = exp_prof[,ind.sam.sort]
typeof(exp_prof)
write.csv(exp_prof, paste(proj_dir,'count.csv', sep = ''))

