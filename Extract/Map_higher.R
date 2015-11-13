  hcc_data = read.csv('C:/Users/NoorPratap/Dropbox/honours/Extract/hcc_data/msb145122-sup-0006-Dataset3/Dataset 3.csv')
lev = levels(hcc_data$Patient.2177)
genes = as.character(hcc_data$GENES[which( hcc_data$Patient.2177 == lev[3] | hcc_data$Patient.2177 == lev[1] | 
                                      hcc_data$Patient.2177 == lev[4])]) #contains genes present(high,med,low) in patient 1
write(genes, 'hcc_genes.txt') #contains all present genes in patient 1 in high conc
hepatocytes_all = normal.liver[which(normal.liver$Cell.type == "hepatocytes"),] # Will combine high,med, low
hepatocytes_all$Level = set.level(hepatocytes_all$Level, conc.hep, 2)
genes.normal.hep.present = as.character(hepatocytes_all$Gene[which(hepatocytes_all$Level == 'Present')])

write(genes.normal.hep.present, 'hep_normal_genes.txt')
mod.canc.liver.2 = data.frame(Gene = rep(unique(cancer.liver$Gene),2),  #Contains the modified levels of liver cancer
                              Level = factor(rep(c("Present", "Not detected"), l)), #and their levels w.r.t proportions
                              count = rep(0L, 2*l), total = rep(0L, 2*l))

i = 0
for (x in unique.cancer.genes) #Sets the modified cancer liver with thresholds and counts
{
  genes.indexes = cancer.liver$Gene == x
  mod.canc.liver.2$Gene[c(2*i +1 , 2*i + 2)] <- x
  mod.canc.liver.2$Level[c(2*i +1 , 2*i + 2)] <- c("Present", "Not detected")
  mod.canc.liver.2$count[c(2*i +1 , 2*i + 2)] <- c(cancer.liver$Count.patients[which(genes.indexes & high.indexes)]
                                                 + cancer.liver$Count.patients[which(genes.indexes & med.indexes)]
                                                + cancer.liver$Count.patients[which(genes.indexes & low.indexes)], 
                                                   cancer.liver$Count.patients[which(genes.indexes & not.det.indexes)] )
  mod.canc.liver.2$total[c(2*i +1 , 2*i + 2)] <- mod.canc.liver.2$count[2*i + 1] + mod.canc.liver.2$count[2*i + 2]
  levels = sapply(c(0.5, 0.75, 0.9, 1), function(x)
  {
    calculate.level(x, mod.canc.liver.2$count[c(2*i + 1)], mod.canc.liver.2$count[c(2*i + 2)], c('Present', 'Not detected'))
  })
  mod.canc.liver.2$actual.level.50[c(2*i +1 , 2*i + 2)] <- levels[1]
  mod.canc.liver.2$actual.level.75[c(2*i +1 , 2*i + 2)] <- levels[2]
  mod.canc.liver.2$actual.level.90[c(2*i +1 , 2*i + 2)] <- levels[3]
  mod.canc.liver.2$actual.level.100[c(2*i +1 , 2*i + 2)] <- levels[4]
  i = i + 1
}
m = mod.canc.liver.2
hepatocytes_all$Canc.Level.50 = m$actual.level.50[indexes.hep]
hepatocytes_all$Canc.Level.75 = m$actual.level.75[indexes.hep]
hepatocytes_all$Canc.Level.90 = m$actual.level.90[indexes.hep]
hepatocytes_all$Canc.Level.100 = m$actual.level.100[indexes.hep]

####Contains presence or absence of cancer w.r.t each gene in different levels of cancer
list.cancers = list(hepatocytes_all$Canc.Level.50, hepatocytes_all$Canc.Level.75, hepatocytes_all$Canc.Level.90, hepatocytes_all$Canc.Level.100)

genes.cancer = sapply(   #contains all genes associated with different levels of cancer (Present as well as not detected but not NA)
        list.cancers, function(x){
        #print(x)
       as.character(hepatocytes_all$Gene[is.na(x) == F])
        #print(a)
      }
)
names(genes.cancer) <- c("50", "75", "90", "100")

genes.cancer.present <- sapply(list.cancers, function(x) #genes present in all different levels
  {
   as.character(hepatocytes_all$Gene[is.na(x) == F & x == 'Present' ])
})
names(genes.cancer.present) <- c("50", "75", "90", "100")

genes.cancer.absent <- sapply(list.cancers, function(x)
  {
    hepatocytes_all$Gene[is.na(x) == F && x == 'Not detected']
})
names(genes.cancer.absent) <- c("50", "75", "90", "100")


find.gene.mismatch <- function(gene.list.1, gene.list.2)
{
  return(gene.list.1[!is.na(match(gene.list.1, gene.list.2))])
}

diff.expressed = get.mismatch(hepatocytes_all$Level, hepatocytes_all$Gene, list.cancers, length(hepatocytes_all$Gene))
diff.expressed.present = mapply(function(x,y)
      {
        find.gene.mismatch(x,y)
      },genes.cancer.present, diff.expressed)

names(diff.expressed) <- c("50", "75", "90", "100")
names(diff.expressed.present) <- c("50", "75", "90", "100")

write(diff.expressed$`50`, 'diff.expressed.50.txt')
write(diff.expressed$`75`, 'diff.expressed.75.txt')
write(diff.expressed$`90`, 'diff.expressed.90.txt')
write(diff.expressed$`100`, 'diff.expressed.100.txt')

write(diff.expressed.present$`50`, 'diff.expressed.present.50.txt')
write(diff.expressed.present$`75`, 'diff.expressed.present.75.txt')
write(diff.expressed.present$`90`, 'diff.expressed.present.90.txt')
write(diff.expressed.present$`100`, 'diff.expressed.present.100.txt')


####Comparing cancer patients
patient.list = list(hcc_data$Patient.2177, hcc_data$Patient.2280, hcc_data$Patient.2556, hcc_data$Patient.2766,
                    hcc_data$Patient.3196, hcc_data$Patient.3477)
patient.ids = colnames(hcc_data)[2:length(colnames(hcc_data))] 
patient.similarites = list()


find.patient.similarity <- function()
{
  ans = sapply(patient.list, function(x)
    {
    sapply(patient.list, function(y)
      {
      as.character(hcc_data$GENES[x == y])
    })
    
  })
  return(ans)
}
ans = find.patient.similarity()

featurise <- function(Patient.levels, vals)
{
  Vectorised.levels = rep(-2, length(Patient.levels))
  Vectorised.levels[Patient.levels == 'Weak'] = vals[2]
  Vectorised.levels[Patient.levels == 'Moderate'] = vals[3]
  Vectorised.levels[Patient.levels == 'Strong'] = vals[4]
  Vectorised.levels[Patient.levels == 'Negative'] = vals[1]
  return(Vectorised.levels)
}
get.featurised.levels <- function(patient.list, values)
{
  Levels = c()
  Levels <- sapply(patient.list, function(x)
    {
      featurise(x, values)
  })
  names(Levels) <- patient.ids
  return(Levels)
}
set.featurised.level <- function(main.data, Levels)
{
  Featurised.vec = hcc_data
  Featurised.vec$Patient.2177 = Levels[,1]
  Featurised.vec$Patient.2280 = Levels[,2]
  Featurised.vec$Patient.2556 = Levels[,3]
  Featurised.vec$Patient.2766 = Levels[,4]
  Featurised.vec$Patient.3196 = Levels[,5]
  Featurised.vec$Patient.3477 = Levels[,6]
  return(Featurised.vec)
}
get.kmeans <- function(Levels)
{
  scaled.data = scale(t(Levels), scale=F)
  no.of.clusters = c(2,3,4,5)
  fits.kmeans = sapply(no.of.clusters, function(x)
  {
    kmeans(scaled.data, x)
  }
  )
  return(fits.kmeans)
}
get.hcl <- function(Levels)
{
  scaled.data = scale(t(Levels), scale=F) 
  distance <- dist(scaled.data)
  hc <- hclust(distance)
  plot(hc)
}
Levels.1 = get.featurised.levels(patient.list, c(-1, 0, 0.5, 1))
Featurised.vec.1 = set.featurised.level(hcc_data, Levels.1)

Levels.2 <- get.featurised.levels(patient.list, c(-1, 1, 1, 1))
Featurised.vec.2 = set.featurised.level(hcc_data, Levels.2)

fit.kmeans.1 = get.kmeans(Levels.1)
fit.kmeans.2 = get.kmeans(Levels.2)

get.hcl(Levels.1)
get.hcl(Levels.2)

hcc.data.complete = read.csv('C:/Users/NoorPratap/Dropbox/honours/Extract/hcc_data/msb145122-sup-0008-Dataset5/Dataset 5.csv')
Levels.3 = get.featurised.levels(patient.list, c(-1, 1, 1, 1))
Featurised.vec.3 = set.featurised.level(hcc_data, Levels.2)
fit.kmeans.3 = get.kmeans(Levels.3)
get.hcl(Levels.3)
