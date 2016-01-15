######DEC 14####################
types.cancer = levels(canc$Tumor)
unique.cancer.genes = unique(canc$Gene)
genes.all = c()  #Genes measured for all cancers
for(i in unique(canc$Gene))
{
  if(length(which(i == canc$Gene)) == 80)
  {
    print(i)
    genes.all = c(genes.all, i)
  }
}
genes.not.all.indexes = which(is.na(match(unique.cancer.genes, genes.all)))#Indexes of absent genes in unique cancer genes
genes.not.all = unique.cancer.genes[genes.not.all.indexes] #All genes not measured for all cancers

cancer.not.present = sapply(genes.not.all, function(x)
{
    levels(canc$Tumor)[which(is.na(match(levels(canc$Tumor), canc$Tumor[which(canc$Gene == x)])))]
}) #Associates cancer names with each gene which have not been measured for all 20 cancers
#genes.not.all = setdiff(unique.cancer.genes, genes.all)  #Alternative approach if only gene names required
cancer.all = data.frame(Genes = rep(as.character(unique.cancer.genes), each = 40), 
                        Tumor = rep(levels(canc$Tumor), each = 2, length(unique.cancer.genes)),
                        Level = rep(c('Present', 'Not detected'), 20*length(unique.cancer.genes)),
                        count = rep(c(0),40*length(unique.cancer.genes)), 
                        total = rep(c(0),40*length(unique.cancer.genes)),
                        actual.level.50 = rep(c('Present', 'Not detected'), 20*length(unique.cancer.genes))
                        )

j = 1
for(i in unique.cancer.genes)
{
   for(k in levels(canc$Tumor))
   {
     index.cancer = which(canc$Gene == i & canc$Tumor == k)
     
     if(length(index.cancer) != 0)
     {
       cancer.all$count[j] = sum(canc$Count.patients[index.cancer][1:3])
       cancer.all$count[j + 1] = canc$Count.patients[index.cancer][4]
       #a = calculate.level(0.5, cancer.all$count[j], cancer.all$count[j+1], 
      #                         c("Present", "Not detected"))
       #print(a)
       cancer.all$actual.level.50[j] = calculate.level(0.5, cancer.all$count[j], cancer.all$count[j+1], 
                                                       c("Present", "Not detected"))
       cancer.all$actual.level.50[j+1] = cancer.all$actual.level.50[j]
       cancer.all$total[j] = canc$Total.patients[index.cancer][1]
       cancer.all$total[j+1] = canc$Total.patients[index.cancer][1]
       if(cancer.all$total[j] != (cancer.all$count[j] + cancer.all$count[j+1]))
         print("TOTAL MISMATCH")
     }
     else
     {
       cancer.all$count[j] = NA
       cancer.all$count[j + 1] = NA
       cancer.all$actual.level.50[j] = NA
       cancer.all$actual.level.50[j+1] = NA
       cancer.all$total[j] = NA
       cancer.all$total[j+1] = NA
     }
     j = j + 2
  }
}

#Contains the final result gene wise with respect to each cancer and its gene
cancer.all.gene.wise = data.frame(Gene = unique.cancer.genes, 
                                  Breast.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  Carcinoid = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  Cervical.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  Colorectal.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  Endometrial.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  Giloma.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  HeadAndNeck.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  Liver.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  Lung.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  Lymphoma.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  Melanoma.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  Ovarian.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  Pancreatic.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  Prostate.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  Renal.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  Skin.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  Stomach.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  Testis.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  Thyroid.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  Urotheal.Cancer = as.character(rep(c('Present'), length(unique.cancer.genes))),
                                  stringsAsFactors = F)

for(i in seq(length(unique.cancer.genes)))
{
  gene = unique.cancer.genes[i]
  vals = c(as.character(unique.cancer.genes[i]))
  indexes = which(cancer.all$Genes == gene)
  if(length(indexes) != 40)
  {
    print("Length Not 40")
    print(indexes)
  }
  for(k in seq(1,length(indexes),2))
  {
    vals = c(vals, as.character(cancer.all$actual.level.50[indexes[k]]))
  }
  cancer.all.gene.wise[i,] = vals
}


###Starting with the mca analysis
data.cancer.gene = t(cancer.all.gene.wise)
data.cancer.gene = data.cancer.gene[-1, ]
data.cancer.gene = data.frame(data.cancer.gene)
colnames(data.cancer.gene) = cancer.all.gene.wise$Gene
not.na.cols = setdiff(cancer.all.gene.wise$Gene, genes.not.all)
data.cancer.gene = data.cancer.gene[, not.na.cols]
data.cancer.gene = data.frame(data.cancer.gene, stringsAsFactors = F)
View(data.cancer.gene)
library(FactoMineR)
library(factoextra)

res.mca.cancer <- MCA(data.cancer.gene)
fviz_screeplot(res.mca.cancer) #1st 2 dim above 1/19
plot(res.mca.cancer)
fviz_mca_biplot(res.mca.cancer)
fviz_mca_ind(res.mca.cancer)

metabolic.genes = read.delim('mca_all_cancer/metabolic_genes.txt', header = F)
metabolic.genes = as.character(metabolic_genes$V1[1:3750])
#Contains genes that are present in my data cancer
metabolic.genes.not.present = metabolic.genes[which(is.na(match(metabolic.genes, colnames(data.cancer.gene))))]
metabolic.genes.present = setdiff(metabolic.genes, metabolic.genes.not.present)
data.metabolic.genes = data.cancer.gene[, match(metabolic.genes.present, colnames(data.cancer.gene))]
View(data.metabolic.genes)
res.mca.cancer.metabolic = MCA(data.metabolic.genes)
fviz_screeplot(res.mca.cancer.metabolic)
fviz_mca_biplot(res.mca.cancer.metabolic, select.var = list(name = metabolic.genes.present.graph, cos2 = 0.6))

metabolic.genes.present.graph = sapply(metabolic.genes.present, function(x)
  {
  paste(x, '_Present', sep = '')
})
metabolic.genes.present.graph.abs = sapply(metabolic.genes.present, function(x)
{
    paste(x, '_Not detected', sep = '')
})
metabolic.genes.present.graph = c(metabolic.genes.present.graph, metabolic.genes.present.graph.abs)
remove(metabolic.genes.present.graph.abs)
fviz_mca_biplot(res.mca.cancer, label = 'var',  select.var = list(name = c('ENSG00000000419_Present', 'ENSG00000000938_Not detected')),col.var = 'red')
fviz_mca_var(res.mca.cancer, label = 'none', select.var = list(name = metabolic.genes.present.graph, cos2 = 0.0), axes = 1:2 )
fviz_mca_biplot(res.mca.cancer, map = 'rowprincipal', label = 'none', select.var = list(name = metabolic.genes.present.graph, cos2 = 0.6), axes = 1:2) + 
theme_minimal()  

fviz_mca_ind(res.mca.cancer.metabolic)
#Learning MCA (Refer important links)
#library(FactoMineR)
#library(factoextra)

########fviz_screeplot(res.mca) finds suitable dimensions
# Our data contains 13 rows and 4 columns.
#If the data were random, the expected value of the eigenvalue for each
#axis would be 1/(nrow(housetasks)-1) = 1/12 = 8.33% in terms of rows.
#Likewise, the average axis should account for 1/(ncol(housetasks)-1) = 1/3 = 33.33% 
#in terms of the 4 columns

######plot(res.mca) 
#######fviz_mca_biplot(res.mca) better representation of the above 
#######fviz_mca_biplot(res.mca) + theme_minimal() even better version
#Note in the above you can always arguments such as axes = c(1,2), col.row = "blue", col.col = "red"
#Above is a symmetric plot in which distance b/w row/row or column/column can be interpreted but not
#b/w column and row 
 
####fviz_mca_var(res.mca) plots only variables(columns) 
# library("corrplot")
######corrplot(var$contrib, is.corr = FALSE) plots the contributions of variables along different dimensions
######fviz_contrib(res.mca, choice = "var", axes = 1) contribution of variables along dimension 1 
#fviz_contrib(res.mca, choice = "var", axes = 1:2) total contribution
#fviz_contrib(res.mca, choice = "var", axes = 1, top = 10) top 10 variables

#####fviz_mca_var(res.mca, col.var="contrib")+
#####scale_color_gradient2(low="white", mid="blue", 
#                      high="red", midpoint=2)+theme_minimal() 
# plots the variables but also shows which variable has larger contribution
# 
# 
# 
# Individuals
#fviz_contrib(res.mca, choice = "ind", axes = 1:2, top = 20) contribution individual wise
#fviz_mca_ind(res.mca, col.ind="contrib")+
#scale_color_gradient2(low="white", mid="blue", 
#                      high="red", midpoint=0.85)+theme_minimal()
#fviz_mca_biplot(res.mca, map ="rowprincipal", arrow = c(TRUE, TRUE))


###Hierarchical Clustering
res.cancer.hcpc <- HCPC(res.mca.cancer)
res.cancer.metabolic.hcpc <- HCPC(res.mca.cancer.metabolic)
#res.cancer.hcpc$data.clust its last column contains the cluster to which that row belongs along with the original info
#res.cancer.hcpc$desc.var
#res.cancer.hcpc$desc.ind
#res.cancer.hcpc$call
#res.cancer.hcpc$desc.axes
plot(res.cancer.hcpc, axes = c(1,2), choice = 'tree', title = 'Clustering')
plot(res.cancer.metabolic.hcpc,  choice = 'tree', title = 'Clustering')
###MCA patient wise
hepatocytes = read.csv('hepatocytes.csv')
hcc.data.complete = read.csv('~/Dropbox/honours/Extract/hcc_data/msb145122-sup-0008-Dataset5/Dataset 5.csv', stringsAsFactors = F)
indexes.match.hep =  match(hcc.data.complete$GENES, hepatocytes$Gene) #Contains the match of indexes of hcc in hepatocytes
hcc.data.complete = transform(hcc.data.complete, HPA.50 = as.character(hepatocytes$Canc.Level.50[indexes.match.hep]))

levels(hcc.data.complete$HPA.50) = c(levels(hcc.data.complete$HPA.50), 'None')
hcc.data.complete$HPA.50[hcc.data.complete$HPA.50 == 'Not detected'] = 'None'

data.hcc.complete = t(hcc.data.complete)
data.hcc.complete = data.frame(data.hcc.complete)
View(data.hcc.complete)
data.hcc.complete = data.hcc.complete[-1,]

change.vals.rows <- function(data, orig.vals, rep.vals)
  ###This function changes the vals row wise and replaces them with rep.vals
  ###data - data.frame on which to be acted
  ###orig.vals - types of vals in rows,a vector
  ###rep.vals - vals orig.vals to be replaced with, a vector(Note orig.vals and rep.vals should
  ###have same length)
  ###eg if in data (orig vals) 'Present' and 'None' is present and in rep.vals 'N' and 'P' replace
  ###them orig vals with rep.vals
{
  no.cols = length(colnames(data))
  for(i in seq(no.cols))
  {
    levels(data[,i]) = c(levels(data[,i]), rep.vals)
    for(j in seq(length(orig.vals)))
        data[, i][which(data[, i] == orig.vals[j])] = rep.vals[j]
  }
  return(data)
}
data.hcc.complete = change.vals.rows(data.hcc.complete, c('Present', 'None'), c('P', 'N'))
res.mca.hcc = MCA(data.hcc.complete)
fviz_screeplot(res.mca.hcc)
fviz_mca_biplot(res.mca.hcc, select.var = list(cos2 = 20), geom = c('arrow', 'text'))
res.hcpc.hcc = HCPC(res.mca.hcc)
plot(res.hcpc.hcc, choice = 'tree', title = 'Clustering')

####Genes for which all cancers have same values
pos.vals = c('Present', 'Not detected')
genes.same = c()
for(i in seq(length(genes.all)))
{
  if(sum(data.cancer.gene[,i] == pos.vals[1]) == 20  || sum(data.cancer.gene[,i] == pos.vals[2]) == 20)
    genes.same = c(genes.same, genes.all[i])
}
indexes.match.same = match(genes.same, genes.all) #indexes of genes for which all cancers have
#same concentration in colnames(data.cancer.gene)
res.mca.cancer.same = MCA(data.cancer.gene[,indexes.match.same])



###CELL CYCLE
get.entrez <- function(filename)
{
  biological.process <- read.delim(filename)
  biological.process.genes.list <- as.character(unique(biological.process[,4]))#column 4 is entrez ids
  biological.process.genes.entrez <- unlist(sapply(biological.process.genes.list, function(x)
  {
    strsplit(x, split = ';')
  }))
  new_list = strsplit(filename, split = '/')
  new_list[[1]][length(new_list[[1]])] = 'entrez.genes.txt'
  write(biological.process.genes.entrez, paste(new_list[[1]], collapse = '/'))
}

get.ensembl <- function(filename)
{
  cell.cycle.ensembl <- read.csv(filename)
  cell.cycle.genes.ensembl <- unique(cell.cycle.ensembl[,1]) #column 1 is ensembl ids in data frame
  return(cell.cycle.genes.ensembl)  
}
get.entrez('mca_all_cancer/cell_cycle/uniprot-cell+cycle.xls')
cell.cycle.genes.ensembl <- get.ensembl('mca_all_cancer/cell_cycle/ensembl.csv')

###Generalised BP
biological.processes <- c('cell_cycle', 'apoptosis')
lapply(biological.processes, function(x)
  {
  file = 'mca_all_cancer/cell/uniprot.xls'
  new_path = strsplit(file, split = '/')
  new_path[[1]][2] = x
  get.entrez(paste(new_path[[1]], collapse = '/'))
})

bp.genes.ensembl <- lapply(biological.processes, function(x)
  {
  file = 'mca_all_cancer/cell/ensembl.csv'
  new_path = strsplit(file, split = '/')
  new_path[[1]][2] = x
  get.ensembl(paste(new_path[[1]], collapse = '/'))
})
names(bp.genes.ensembl) <- biological.processes

which(is.na(match(cell.cycle.genes.ensembl, cancer.all.gene.wise$Gene)))
cell.cycle.genes.absent <- cell.cycle.genes.ensembl[which(is.na(match(cell.cycle.genes.ensembl, cancer.all.gene.wise$Gene)))]


cell.cycle.genes.mca <- setdiff(cell.cycle.genes.ensembl, cell.cycle.genes.absent)
data.cancer.cell.cycle <- data.cancer.gene[, match(cell.cycle.genes.mca, colnames(data.cancer.gene))]
res.mca.cancer.cell.cycle <- MCA(data.cancer.cell.cycle)
fviz_mca_ind(res.mca.cancer.cell.cycle)
