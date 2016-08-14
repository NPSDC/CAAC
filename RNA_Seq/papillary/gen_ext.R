install.packages('~/varSelRF_0.7-5.tar.gz', repos = NULL, type = 'source')
library(varSelRF)
library(snow)

exp_fpqm_tumor_log <- log2(exp_fpqm_tumor+1)
df.stages$sample.id == colnames(exp_fpqm_tumor)
indexes.stages.reported <- which(df.stages$stage != 'not reported')
stages.levels <- df.stages$stage[indexes.stages.reported]
stages.levels <- droplevels(stages.levels)

exp_fpqm_tumor_log_reported <- exp_fpqm_tumor_log[,indexes.stages.reported]
exp_fpqm_tumor_log_reported <- t(exp_fpqm_tumor_log_reported)
rownames(exp_fpqm_tumor_log_reported) <- df.stages$short.id[indexes.stages.reported]

exp_fpqm_tumor_reported <- exp_fpqm_tumor[,indexes.stages.reported]
exp_fpqm_tumor_reported <- t(exp_fpqm_tumor_reported)
rownames(exp_fpqm_tumor_reported) <- df.stages$short.id[indexes.stages.reported]

exp_prof_tumor_reported = exp_prof_tumor[,indexes.stages.reported]
exp_prof_tumor_reported = t(exp_prof_tumor_reported)
 
cl <- makeMPIcluster(8)
tumor.fpqm.varSelRF <- varSelRF(xdata = exp_fpqm_tumor_reported, Class = stages.levels, 
                                 keep.forest = T)
tumor.fpqm.boot <- varSelRFBoot(exp_fpqm_tumor_reported, stages.levels, bootnumber = 10, usingCluster = T,
                                TheCluster = cl, srf= tumor.fpqm.varSelRF)
tumor.fpqm.log.varSelRF <- varSelRF(xdata = exp_fpqm_tumor_log_reported, Class = stages.levels, 
                                keep.forest = T)
tumor.fpqm.log.boot <- varSelRFBoot(exp_fpqm_tumor_log_reported, stages.levels, bootnumber = 10, usingCluster = T,
                                TheCluster = cl, srf= tumor.fpqm.log.varSelRF)
tumor.exp.varSelRF <- varSelRF(xdata = exp_prof_tumor_reported, Class = stages.levels, 
                                keep.forest = T)
tumor.exp.boot <- varSelRFBoot(exp_prof_tumor_reported, stages.levels, bootnumber = 10, usingCluster = T,
                                TheCluster = cl, srf= tumor.exp.varSelRF)

stopCluster(cl)

####

plot.lists <- function(varSelRF.lists, index)
{
  library(varSelRF)
  library(ggplot2)
  len <- length(varSelRF.lists)
  #print(len)
  par(mfrow=c(1,len))
  for(i in varSelRF.lists)
    plot(i,which=index)

}

plot.lists(list(tumor.fpqm.varSelRF, tumor.fpqm.log.varSelRF, tumor.exp.varSelRF), index = 1)
plot.lists(list(tumor.fpqm.varSelRF, tumor.fpqm.log.varSelRF, tumor.exp.varSelRF), index = 2)
plot(tumor.exp.boot, subject.names = df.stages$short.id)
selProbPlot(tumor.exp.boot)
summary(tumor.exp.varSelRF)

dir.copy = paste(getwd(), 'results/tumor/varSelRF', sep = '/')


write.csv(tumor.exp.varSelRF$selec.history, paste(dir.copy, 'history_exp.csv', sep = '/'))
write.csv(tumor.fpqm.varSelRF$selec.history, paste(dir.copy, 'history_fpqm.csv', sep = '/'))
write.csv(tumor.fpqm.log.varSelRF$selec.history, paste(dir.copy, 'history_fpqm_log.csv', sep = '/'))

tumor.exp.varSelRF$selected.vars
tumor.fpqm.varSelRF$selected.vars
tumor.fpqm.log.varSelRF$selected.vars

tumor.exp.varSelRF$initialOrderedImportances[1:30]
tumor.fpqm.varSelRF$initialOrderedImportances[1:30]
tumor.fpqm.log.varSelRF$initialOrderedImportances[1:30]

##BootStrap
tumor.exp.boot$bootstrap.pred.error
tumor.fpqm.boot$bootstrap.pred.error
tumor.fpqm.log.boot$bootstrap.pred.error

write.csv(tumor.exp.boot$all.data.randomForest$confusion, paste(dir.copy, 'confusion_exp.csv', sep ='/'))
write.csv(tumor.fpqm.boot$all.data.randomForest$confusion, paste(dir.copy, 'confusion_fpqm.csv', sep ='/'))
write.csv(tumor.fpqm.log.boot$all.data.randomForest$confusion, paste(dir.copy, 'confusion_fpqm_log.csv', sep ='/'))

tumor.exp.boot$all.data.randomForest$importance
tumor.fpqm.boot$all.data.randomForest$importance
tumor.fpqm.log.boot$all.data.randomForest$importance

View(tumor.exp.boot$class.predictions)

tumor.exp.boot$number.of.vars
tumor.fpqm.boot$number.of.vars
tumor.fpqm.log.boot$number.of.vars

tumor.exp.boot$overlap
tumor.fpqm.boot$overlap
tumor.fpqm.log.boot$overlap
