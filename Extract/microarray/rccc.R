################################################################
##Not using their script exactly since we they don't normalize the data
##For annotating purpose we use AnnotationDbi since the expression list we get does not contain
##the gene names
library("affy")
library('simpleaffy')
library(gcrma)
setwd('~/Dropbox/honours/Extract/microarray/renal/data/A')
files <- list.files(pattern = '.CEL')
MyData <- ReadAffy(filenames = files)

cell.rma <- rma(MyData)

library(RColorBrewer)
library(affyPLM)
cols <- brewer.pal(8, "Set1")
boxplot(MyData, col=cols)
boxplot(eset, col=cols)
hist(eset, col=cols)
celfiles.qc <- fitPLM(MyData)
RLE(celfiles.qc, main="RLE")
distance <- dist(t(as.matrix(eset)),method="maximum")
clusters <- hclust(distance)
plot(clusters)

#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO(filename = "GSE781-GPL96_series_matrix.txt", GSEMatrix =TRUE)
#gset <- read.delim(file = 'GSE781-GPL96_series_matrix.txt')
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))
fvarLabels(eset) <- make.names(fvarLabels(eset))

# group names for all samples
gsms <- "01011010110101010"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
eset$description <- fl
design <- model.matrix(~ description + 0, eset)
design1 <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
colnames(design1) <- levels(fl)
fit <- lmFit(eset, design)
fit.1 <- lmFit(gset, design1)
cont.matrix <- makeContrasts(G1-G0, levels=design)
cont.matrix1 <- makeContrasts(G1-G0, levels=design1)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
fit2.1<- contrasts.fit(fit.1, cont.matrix1)
fit2.1 <- eBayes(fit2.1, 0.01 )
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=50000)
tT.1 <- topTable(fit2.1, adjust="fdr", sort.by="B", number=50000)
# load NCBI platform annotation
gpl <- annotation(eset)
platf <- getGEO(gpl, AnnotGPL=F)
ncbifd <- data.frame(attr(dataTable(platf), "table"))

# replace original platform annotation
tT <- tT[setdiff(colnames(tT), setdiff(fvarLabels(eset), "ID"))]
tT <- merge(tT, ncbifd, by="ID")
tT <- tT[order(tT$P.Value), ]  # restore correct order

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

###Plotting
gsms <- "10100101001010101"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  #set group names

# order samples by group
ex <- exprs(gset)[ , order(sml)]
ex1 <- exprs(eset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("cancer","normal")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(eset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(eset)))/2),4,2,1))
title <- paste ("GSE781", '/', annotation(eset), " selected samples", sep ='')
boxplot(ex1, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")


 ##Annotating
source("http://bioconductor.org/biocLite.R")
biocLite("hgu133a.db")
library(hgu133a.db)
library(AnnotationDbi)
tT$Probe.ID <- rownames(tT)
probe.ids <- rownames(tT)
col.db <- columns(hgu133a.db)

trans.df <- select(hgu133a.db, probe.ids, c('SYMBOL', 'ENTREZID', 'ENSEMBL', 'GENENAME'))
req.dfs <- merge(tT, tT.1[,c(1,10,11,12,14,15,16),], by.x = 'Probe.ID', by.y = 'ID')

degs.renal.int.micro <- intersect(map.genes.ids(diff.expressed.all.canc.corr$Renal.Cancer, genes.ids),
                                  unique(req.dfs$Gene.Symbol))
##Testing
sum(req.dfs$adj.P.Val[match(tT$Probe.ID, req.dfs$Probe.ID)] == tT$adj.P.Val)

##Using abs p-val < 0..05
library(dplyr)
temp <- filter(req.dfs, P.Value < 0.05, abs(logFC) > 0.6 )
temp <- filter(req.dfs, adj.P.Val < 0.05, abs(logFC) > 0.6)
setwd('~/Dropbox/honours/Extract/microarray/renal/')

duplicated.ids <- which(duplicated(temp$Gene.Symbol))

  
filtered <- filter.genes(temp,2,9) 
micr.prot.int <- intersect(filtered$Gene.Symbol, map.genes.ids(diff.expressed.all.canc.corr$Renal.Cancer, genes.ids))
mean(filtered$logFC < 0)
down.reg <- filter(filtered, logFC < 0) %>% select(Gene.Symbol) %>% unlist %>% as.character
up.reg <- filter(filtered, logFC > 0) %>% select(Gene.Symbol) %>% unlist %>% as.character

micr.prot.ints = list() ##intersection of intersected degs of microarray with protein out of those
##which are up in microarray and then down and up degs of protein and then same for down
micr.prot.ints[['up']][['up']] = intersect(up.reg, map.genes.ids(diff.expressed.up.all.canc.corr$Renal.Cancer, genes.ids))
micr.prot.ints[['up']][['down']] = intersect(up.reg, map.genes.ids(diff.expressed.down.all.canc.corr$Renal.Cancer, genes.ids))
micr.prot.ints[['down']][['up']] = intersect(down.reg, map.genes.ids(diff.expressed.up.all.canc.corr$Renal.Cancer, genes.ids))
micr.prot.ints[['down']][['down']] = intersect(down.reg, map.genes.ids(diff.expressed.down.all.canc.corr$Renal.Cancer, genes.ids))
names(micr.prot.ints) = c('up', 'down')
names(micr.prot.ints$up) = c('up', 'down')
names(micr.prot.ints$down) = c('up', 'down')

##Protein data
renal.prot.norm <- filter(normal, Tissue == 'kidney')##Normal protein data with supportive reliablity
renal.prot.canc <- filter(canc, Tumor == 'renal cancer')##Canc data with counts of renal


renal.prot.norm$Gene.Sym = rep(c('NA'), length(renal.prot.norm$Gene))
renal.prot.norm$Gene.Sym = add.gene.name(renal.prot.norm, 1, 7, genes.ids)
renal.prot.canc$Gene.Sym = rep(c('NA'), length(renal.prot.canc$Gene))
renal.prot.canc$Gene.Sym = add.gene.name(renal.prot.canc, 1, 8, genes.ids)
cancer.all.gene.wise.corr$Gene.Sym = rep(c('NA'), length(cancer.all.gene.wise.corr$Gene))
cancer.all.gene.wise.corr$Gene.Sym = add.gene.name(cancer.all.gene.wise.corr, 1, 22, genes.ids)

renal.micr.up.int = list()
renal.prot.int = find.sub.int(renal.prot.norm, 7, micr.prot.int)
renal.prot.canc.int = find.sub.int(renal.prot.canc, 8, micr.prot.int)
reliable <- filter(renal.prot.int, Reliability == 'Supportive') %>%  select(Gene.Sym) %>% unlist  %>% unique
unreliable <- filter(renal.prot.int, Reliability == 'Uncertain') %>%  select(Gene.Sym) %>% unlist  %>% unique


###Intersecting degs of microarray with overall protein and not with others
###All present and absent genes in protein data for Renal Cancer
renal.up.prot <- filter(cancer.all.gene.wise.corr, Renal.Cancer == 'Present') %>% select(Gene.Sym) %>% unlist %>% as.character
renal.down.prot <-filter(cancer.all.gene.wise.corr, Renal.Cancer == 'Not detected') %>% select(Gene.Sym) %>% unlist %>% as.character

renal.micr.up.int = list() ##Intersecting up degs of microarray with proteins which are up and down in protein data
renal.micr.down.int = list() ##Intersecting down degs of microarray with proteins which are up and down in protein data
renal.micr.up.int[['up']] = intersect(up.reg, renal.up.prot)
renal.micr.up.int[['down']] = intersect(up.reg, renal.down.prot)
renal.micr.down.int[['up']] = intersect(down.reg, renal.up.prot)
renal.micr.down.int[['down']] = intersect(down.reg, renal.down.prot)
setwd('~/Dropbox/honours/Extract/microarray/renal/protein/up/')
write(renal.micr.up.int$up, 'up.txt')
write(renal.micr.up.int$down, 'down.txt')
setwd('~/Dropbox/honours/Extract/microarray/renal/protein/down/')
write(renal.micr.down.int$up, 'up.txt')
write(renal.micr.down.int$down, 'down.txt')

###testing
#testing renal.prot.norm
a= sapply(unique(renal.prot.norm$Gene), function(x){
  length(which(x == renal.prot.norm$Gene)) == 1
})
b = which(a==T)
c = match(unique(renal.prot.norm$Gene)[b], renal.prot.norm$Gene)
