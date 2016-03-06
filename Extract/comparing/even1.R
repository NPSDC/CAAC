simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}
add.vals <- function(data,gene.index, tissue.index, count.index, total.index) #adds the 2 levels for each gene and each cancer
{
  unique.cancer.genes = unique(data[,gene.index])
  canc.add = data.frame(Gene = unique.cancer.genes, 
                                    Breast.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    Carcinoid = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    Cervical.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    Colorectal.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    Endometrial.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    Giloma.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    HeadAndNeck.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    Liver.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    Lung.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    Lymphoma.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    Melanoma.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    Ovarian.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    Pancreatic.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    Prostate.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    Renal.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    Skin.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    Stomach.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    Testis.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    Thyroid.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    Urotheal.Cancer = as.character(rep(c('NA'), length(unique.cancer.genes))),
                                    stringsAsFactors = F)
  gene.old = data[,gene.index][1]
  for(i in seq(1,length(data[,gene.index])/8/13/13,4))
  {
    gene.new = data[,gene.index][i]
    j = 1
    if(gene.new != gene.old)
    {
      gene.old = gene.new
      j = j + 1
    }
    tissue = as.character(data[i,tissue.index])
    tissue = simpleCap(tissue)
    tissue = unlist(strsplit(tissue, ' '))
    print(length(tissue))
    if(length(tissue == 2))
      tissue = paste(tissue[1], tissue[2], sep = '.')
    tissue.col = match(tissue, colnames(canc.add))
    print(tissue)
    print(tissue.col)
    if(data[i+3, count.index] != data[i+3,total.index])
      canc.add[j,tissue.col] = 'Present'
    else
      canc.add[j,tissue.col] = 'Not detected'
  }
}
canc.present.1 = add.vals(canc, 1, 2, 4, 5)
