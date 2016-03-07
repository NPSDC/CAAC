source('preprocess.R')

tissue.indexes = list(c(4), c(8), c(9,31), c(11,12), c(7), c(25,32), c(21), c(5,22), c(23), c(26),
                      c(27), c(30), c(19), c(42,43), c(44), c(45), c(47))
tissue.names = list(c('breast'), c('cervix, uterine'), c('colon','rectum'), 
                    c('endometrium 1', 'endometrium 2'), c('cerbral cortex'), 
                    c('oral mucosa', 'salivary gland'), c('liver'), c('bronchus', 'lung'), 
                    c('lymph node'), c('ovary'), c('pancreas'), c('prostate'), c('kidney'), 
                    c('stomach 1', 'stomach 2'), c('testis'), c('thyroid gland'), c('urinary bladder'))
tissue.names = unique(normal$Tissue)
colorectal.tissue = preprocess(tissue.indexes[[3]], normal, 2, 4, tissue.names)

tumor.indexes = c(c(1),c(3:10),c(12:15),c(17:20))
tumor.names = as.character(unique(canc$Tumor))
tissues.names = unlist(sapply(tumor.names[tumor.indexes], alter.name))

tissues.normal = lapply(tissue.indexes, function(x)
  {
    preprocess(x, normal, 2, 4, tissue.names)  
})
names(tissues.normal) = tissues.names

tissues.canc = lapply(tumor.indexes, function(x)
  {
  reduce.levels(x, c('Present', 'Not detected'), list('Present', 'Not detected'), tumor.names)
})
names(tissues.canc) = tissues.names

diff.expressed.all.canc = mapply(function(x,y)
  {
    find.diff.expressed(x,y,c(4,5),c(1,1))
}, tissues.normal, tissues.canc)
names(diff.expressed.all.canc) = tissues.names

diff.expressed.up.all.canc = mapply( function(x,y)
  {
  find.up.regulated(x,5,1,y)
},diff.expressed.all.canc,tissues.canc)
names(diff.expressed.up.all.canc) = tissues.names

diff.expressed.down.all.canc = find.down.regulated(diff.expressed.all.canc, diff.expressed.up.all.canc, tissues.names)
