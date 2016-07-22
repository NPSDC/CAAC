##Parameters to be given wd, mRNA_file_path, total case files, file_name, proj_dir, cases

library(biomaRt)
library(rjson)

get.mRNA_files = function(mRNA_files_path, total_case_files)
{
  mRNA_files = list.files(mRNA_files_path)[!grepl('txt',list.files(mRNA_files_path))] ##all files except manifest
  if(length(mRNA_files) == total_case_files)
    return(mRNA_files)
  else
    print('sucker')
}

get.json.mapping <- function(path,file_name)
{
  ##returns the file_id corresponding to patient id
  file_name = paste(path,file_name,sep='/')
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
  return(sample.file.map)
}

unzip.files <- function(mRNA_files, mRNA_files_path, proj_dir, total_case_files)
{
  ##Unzips the files from their original folders and copies them to previous directory and deleting the original directory
  files.to.be.copied = c()
  setwd(mRNA_files_path)
  for(i in mRNA_files)
  {
    R.utils::gunzip(paste(i,list.files(i)[1],sep = '/'))
    files.to.be.copied = c(files.to.be.copied, paste(i,list.files(i)[1], sep='/'))
  }
  file.copy(files.to.be.copied,'.')
  unlink(mRNA_files, recursive = T) ##Deleting the directories
  mRNA_files = list.files()[!grepl('MANIFEST.txt',list.files())]
  setwd(proj_dir)
  if(length(mRNA_files) == total_case_files)
    return(mRNA_files)
  else
  {
    print('error in zipping')
    return(-1)
  }
}
##So far unzipped the expression values and copied them to previous folder and deleted the original directories
##Reading expression value of each file for each gene checking whehther the genes match
##Then store them in a proper frame exp_prof

get.genes.files <- function(mRNA_files, mRNA_files_path, proj_dir)
{
  setwd(mRNA_files_path)
  file1 = read.delim(mRNA_files[1], header = F)
  ens.ids.all = sort(as.character(file1$V1)) ##all ensembl ids including transcript
  setwd(proj_dir)
  return(ens.ids.all)
}

remove.dots <- function(ens.ids.all)
{
  g = sapply(ens.ids.all, function(x) 
  {
    unlist(strsplit(x, split = '.', fixed = T))[1]
  }) ##removing the symbols after .
}

##For removing ids that are greater that are not entrez ids

get.entrez <- function(ens.ids.all)
{
  ensembl=useMart("ensembl")
  #listDatasets(ensembl)
  ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)
  #listFilters(ensembl) 
  #listAttributes(ensembl)
  genes.entrez = getBM(attributes = c('ensembl_gene_id', 'entrezgene'), filters = 'ensembl_gene_id', values = g, mart = ensembl)
  genes.entrez = genes.entrez[order(genes.entrez$ensembl_gene_id),]
  return(genes.entrez)
}

get.ens.ids.with.entrez <- function(genes.entrez)
{
  ens.ids = genes.entrez$ensembl_gene_id[(!is.na(genes.entrez$entrezgene))] #ens ids that are not na
  length(unique(ens.ids)) != length(ens.ids)
  ens.ids = unique(ens.ids)
  length(ens.ids)
  return(ens.ids)
}
##We get all the ensembl ids that have entrez genes

##Testing
test <- function(ens.ids, genes.entrez)
##Making sure that for ens ids we have entrez ids
{
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
}
get.indexes <- function(ens.ids, g)
{
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
  return(indexes.ens)
}
##Building a count matrix
get.count.matrix <- function(ens.ids.all, indexes.ens, mRNA_files, sample.file.map, mRNA_files_path,proj_dir)
{
  setwd(mRNA_files_path)
  file1 = read.delim(mRNA_files[1], header = F)
  file1 = file1[order(file1$V1),]
  ens.ids = ens.ids.all[indexes.ens] ##All ens ids which have gene transcripts
  exp_prof = file1$V2[indexes.ens]
  for(i in c(2:length(mRNA_files)))
  {
  file = read.delim(mRNA_files[i], header = F)
  file = file[order(file$V1),]
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
  exp_prof = cbind(exp_prof, file$V2[indexes])
  print(i)
  }
  #print(exp_prof[1:15,])
  remove(file)
  
  file.names.match = match(mRNA_files, sample.file.map$file.ids)
  if(sum(is.na(file.names.match)) != 0)
  {
    print('need to talk')
    return(-1)
  }
  #print(exp_prof[1:10,])
  colnames(exp_prof) = sample.file.map$pat.ids[file.names.match]
  rownames(exp_prof) = ens.ids
  exp_prof = data.frame(exp_prof)
  sum(ens.ids == sort(ens.ids)) == length(ens.ids)
  ind.sort = match(ens.ids, rownames(exp_prof)) ##sorted rownames
  for(i in seq_along(colnames(exp_prof)))
    exp_prof[,i] = exp_prof[,i][ind.sort]
  rownames(exp_prof) = ens.ids
  
  ##Sorting samples
  sample.ids = sort(colnames(exp_prof))
  
  ind.sam.sort = match(sample.ids, colnames(exp_prof))
  exp_prof = exp_prof[,ind.sam.sort]
  #print(exp_prof[1:10,])
  setwd(proj_dir)
  return(exp_prof)
}
do.everything <- function(wd, mRNA_file_path, total_case_files, file_name, proj_dir)
{
  sample.file.map = get.json.mapping(wd, file_name)
  files = get.mRNA_files(mRNA_files_path, total_case_files)
  files = unzip.files(files, mRNA_files_path, proj_dir, total_case_files)
  ens.ids.all = get.genes.files(tfiles, mRNA_files_path, proj_dir)
  g = remove.dots(t.ens.ids.all)
  genes.entrez <- get.entrez(ens.ids.all)
  ens.ids <- get.ens.ids.with.entrez(genes.entrez)
  test(ens.ids, genes.entrez)
  indexes.ens <- get.indexes(ens.ids, g)
  exp_fpqm <- get.count.matrix(ens.ids.all, indexes.ens, files, sample.file.map , mRNA_files_path, proj_dir)
  
}
#typeof(exp_prof)
#write.csv(exp_prof, paste(proj_dir,'count_mrna.csv', sep = ''))

