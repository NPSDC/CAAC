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
data.cancer.gene = data.frame(data.cancer.gene)
View(data.cancer.gene)
library(FactoMineR)
mca.data.cancer <- MCA(data.cancer.gene, ind.sup = 1)

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

