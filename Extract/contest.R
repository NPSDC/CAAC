indexes.red.log <- sapply(seq(16613), function(x)
  {
  sum(cancer.all.gene.wise[x,] == 'Present') == 20 | sum(cancer.all.gene.wise[x,] == 'Not detected') == 20
})
genes.red <- cancer.all.gene.wise$Gene[indexes.same]
indexes.red <- match(genes.red, cancer.all.gene.wise$Gene)
test.data = data.frame(Gene = genes.red, 
                                  Breast.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  Carcinoid = as.character(rep(c('Present'), length(genes.red))),
                                  Cervical.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  Colorectal.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  Endometrial.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  Giloma.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  HeadAndNeck.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  Liver.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  Lung.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  Lymphoma.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  Melanoma.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  Ovarian.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  Pancreatic.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  Prostate.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  Renal.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  Skin.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  Stomach.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  Testis.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  Thyroid.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  Urotheal.Cancer = as.character(rep(c('Present'), length(genes.red))),
                                  stringsAsFactors = F)
le = length(test.data$Gene)
for(i in c(2:21))
{
  test.data[sample(1:le, sample(500:3000, 1)), i ] = 'Not detected'
}
lengths.cancers = list()
len <- function(data.set)
{
  ans2 <- lapply(data.set, length)
  names(ans2) <- names(data.set)
  return(ans2)
}
for(i in seq(20))
{
  lengths.cancers[[i]] = len(ans1[[i]])
}
rands = list()
for(i in seq(10000))
{
  rands[[i]] = sample(1:19, 5)
}
id = c(1:20)
test.max <- lapply(id, function(x)
  {
  find.sim(test1.data, 21, x, setdiff(id,x))  
})
names(test.max) = colnames(test1.data[1:20])
lengths.cancers = lapply(test.max, len)
names(lengths.cancers) = colnames(test1.data[1:20])
maximum <- lapply(lengths.cancers, function(x)
  {
  max(unlist(x))  
})
names(maximum) = colnames(test1.data[1:20])

names.max <- mapply(function(x,y)
  {
  names(which(x ==y))
}, lengths.cancers, maximum)
names(names.max) = colnames(test1.data[1:20])
