library(rjson)
library(dplyr)
pap.cli.json <- fromJSON(file = '~/Downloads/data/papillary/clinical.cart.2016-07-21T12-04-02.150390.json')

pat_stages = lapply(seq_along(pap.cli.json), function(x)
{
  stage = pap.cli.json[[x]][[1]][[1]][[6]]
})
names(pat_stages) = sapply(seq_along(pap.cli.json), function(x)
{
  pat.id = pap.cli.json[[x]][[1]][[1]][[5]]
  pat.id = unlist(strsplit(pat.id, split = '_', fixed = T))[1]
  })

library('TCGAbiolinks')
query = TCGAquery_subtype('kirp')
pap.sub = filter(query) %>% select(patient, tumor_type.KIRP.path.)
length(intersect(pap.sub$patient, names(pat_stages)))

pap.comb = data.frame(name = names(pat_stages), stage = unlist(pat_stages), 
                      type = rep(levels(pap.sub$tumor_type.KIRP.path.)[3], length(names(pat_stages))),
                      stringsAsFactors = F)
ind.comb = match(intersect(pap.sub$patient, names(pat_stages)), pap.comb$name)
pap.comb$type[ind.comb] = pap.sub$tumor_type.KIRP.path.[match(intersect(pap.sub$patient, names(pat_stages)), pap.sub$patient)]
pap.comb = pap.comb[order(pap.comb[,1]), ]
pap.comb$type[pap.comb$type == 3] = 'Unclassified Papillary RCC'

stages = unique(pap.comb$stage)
stages = stages[1:4]
type = unique(pap.comb$type)
type = type[2:3]
