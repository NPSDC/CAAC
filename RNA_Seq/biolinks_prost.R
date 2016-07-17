library(TCGAbiolinks)
query <- TCGAquery(tumor = 'PRAD', platform = 'IlluminaHiSeq_RNASeqV2') ##prad stands for prostate cancer
prad_subtype <- TCGAquery_subtype(tumor = 'prad')
TCGAdownload(query, path = "~/Downloads/data/biolinks/", type = "rsem.genes.results")
count_data <- TCGAprepare(query, dir = '~/Downloads/data/mRNA_fpqm/gdc_download_20160707_095145' )
                                            
a = read.delim('~/Downloads/data/broad_institute/PRAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt')

library('RTCGAToolbox')
getFirehoseDatasets()
stddata = getFirehoseRunningDates()
brcaData = getFirehoseData (dataset="PRAD", runDate="20151101",forceDownload = TRUE,
                            Clinic=TRUE)
data(RTCGASample)

dir = '~/Downloads/data/mRNA_fpqm/gdc_download_20160707_095145'
files <- NULL
dirs <- gsub(".tar.gz","",basename(query$deployLocation))
for (i in seq_along(dirs)) {
  aux <- list.files(file.path(dir,dirs[i]), full.names = TRUE,
                    recursive = TRUE)
  files <- c(files, aux )
}
idx <- grep("MANIFEST|README|CHANGES|DESCRIPTION|DATA_USE",files)
if (length(idx) > 0) {
  files <- files[-idx]
}
