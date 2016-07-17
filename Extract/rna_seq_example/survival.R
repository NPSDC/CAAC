library(GEOquery)
getGEOSuppFiles("GSE20986")
untar("GSE20986_RAW.tar", exdir="data")
cels <- list.files("data/", pattern = "[gz]")
sapply(paste("data", cels, sep="/"), gunzip)

library(affy)

source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite('simpleaffy')
biocLite('affyPLM')

library('simpleaffy')
celfiles <- read.affy(covdesc="phenodata.txt", path="data")
celfiles.gcrma <- gcrma(celfiles)

##Quality Control checks
library(RColorBrewer)
cols <- brewer.pal(8, "Set1")
boxplot(celfiles, col=cols)

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles.gcrm
library(affyPLM)
boxplot(celfiles.gcrma, col=cols)
# the boxplots are somewhat skewed by the normalisation algorithm
# and it is often more informative to look at density plots
# Plot a density vs log intensity histogram for the unnormalised data
hist(celfiles, col=cols)
# Plot a density vs log intensity histogram for the normalised data
hist(celfiles.gcrma, col=cols)

# Perform probe-level metric calculations on the CEL files:
celfiles.qc <- fitPLM(celfiles)
# Create an image of GSM24662.CEL
image(celfiles.qc, which=1, add.legend=TRUE)
image(celfiles.qc, which=4, add.legend=TRUE)
# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero. GSM524665.CEL is an outlier
RLE(celfiles.qc, main="RLE")
NUSE(celfiles.qc, main="NUSE")

eset <- exprs(celfiles.gcrma)
distance <- dist(t(eset),method="maximum")
clusters <- hclust(distance)
plot(clusters)

##Filtering
celfiles.filtered <- nsFilter(celfiles.gcrma, require.entrez=FALSE, remove.dupEntrez=FALSE)
celfiles.filtered$filter.log
samples <- celfiles.gcrma$Target              
samples <- as.factor(samples)
design <- model.matrix(~0 + samples)
colnames(design) <- c("choroid", "huvec", "iris", "retina")
library(limma)
fit <- lmFit(exprs(celfiles.filtered$eset), design)
# set up a contrast matrix to compare tissues v cell line
contrast.matrix <- makeContrasts(huvec_choroid = huvec - choroid, huvec_retina = huvec - retina, huvec_iris <- huvec - iris, levels=design)
huvec_fits <- contrasts.fit(fit, contrast.matrix)
huvec_ebFit <- eBayes(huvec_fits)
topTable(huvec_ebFit, number=10, coef=1)
biocLite("hgu133plus2.db")
library(hgu133plus2.db)
library(annotate)
probeset.list <- topTable(huvec_ebFit, coef=1, number=10000, lfc=4)
gene.symbols <- getSYMBOL(probeset.list$ID, "hgu133plus2")


###Same thing on renal
setwd('~/Dropbox/honours/Extract/microarray/renal/data/A')
cels <- list.files(pattern = "L.gz")
sapply(cels, gunzip)
cel.files <- read.affy(covdesc = 'phenodata.txt')
cel.gcrma <- gcrma(cel.files)
cel.rma <-  rma(cel.files)
cols <- brewer.pal(8, "Set1")
boxplot(cel.files, col=cols)  
boxplot(cel.gcrma, col=cols)  
boxplot(cel.rma, col=cols)  
hist(cel.files, col=cols)
hist(cel.gcrma, col=cols)
hist(cel.rma, col=cols)
cel.qc <- fitPLM(cel.files)
RLE(cel.qc, main="RLE")

#eset <- exprs(cel.gcrma)
eset1 <- exprs(cel.rma)
distance <- dist(t(eset),method="maximum")
clusters <- hclust(distance)
plot(clusters)

cel.filtered <- nsFilter(cel.gcrma, require.entrez=FALSE, remove.dupEntrez=FALSE)
samples <- cel.rma$Target

design <- model.matrix(~0 + samples)
colnames(design) <- c('N', 'C')
library(limma)
fit_3 <- lmFit(exprs(filter$eset), design)
fit_1 <- lmFit(eset1, design)
contrast.matrix_1 <- makeContrasts(N-C, levels = design)
huvec_fits <- contrasts.fit(fit_1, contrast.matrix_1)
huvec_ebFit <- eBayes(huvec_fits)
huvec_fits_2 <- contrasts.fit(fit_3, contrast.matrix_1)
huvec_ebFit_2 <- eBayes(huvec_fits_2)

tT_1 = topTable(huvec_ebFit, number=100000,  sort.by = 'B')
tT_1$ID <- rownames(tT_1)
probe.ids <- rownames(tT_1)

tT_2 = topTable(huvec_ebFit_2, number=100000,  sort.by = 'B')
tT_2$ID <- rownames(tT_2)
probe.ids <- rownames(tT_2)
library(dplyr)
req <- filter(tT, adj.P.Val < 0.05, abs(logFC) > 0.6)

##Annotating
library('hgu133a.db')
library(hgu133plus2.db)
library(annotate)
library(AnnotationDbi)
library(dplyr)
trans.df <- select(hgu133a.db, probe.ids, c('SYMBOL',  'GENENAME'))
data("Affyhgu133aExprtab")
gene.symbols <- getSYMBOL(rownames(tT_1), "hgu133plus2.db")
tT_1 <- cbind(gene.symbols, tT_1)

#####Script
library("affy")
setwd('~/Dropbox/honours/Extract/microarray/renal/data/A')
MyData <- ReadAffy()
eset <- rma(MyData)
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

fvarLabels(eset) <- make.names(fvarLabels(eset))
gsms <- "01011010110101010"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)

eset$description <- fl
des <- model.matrix(~ description + 0, eset)
colnames(des) <- levels(fl)
fit <- lmFit(eset, des)
cont.matrix <- makeContrasts(G1-G0, levels=des)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=50000)


gset <- getGEO(filename = "GSE781-GPL96_series_matrix.txt", GSEMatrix =TRUE)
fvarLabels(gset) <- make.names(fvarLabels(gset))
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }
gset$description <- fl
design1 <- model.matrix(~ description + 0, gset)
colnames(design1) <- levels(fl)
fit.1 <- lmFit(gset, design1)
cont.matrix1 <- makeContrasts(G1-G0, levels=design1)
fit.2 <- contrasts.fit(fit.1, cont.matrix1)
fit.2 <- eBayes(fit.2, 0.01)
tT.2 <- topTable(fit.2, adjust="fdr", sort.by="B", number=50000)
