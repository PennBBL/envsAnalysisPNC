library(parallel)
library(mgcv)

demo <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_demographics_go1_20161212.csv")
envs <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/environment/n9498_go1_environment_factor_scores_tymoore_20150909.csv")
tracker <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_tracker_go1_20161212.csv")
tracker <- tracker[c("bblid","scanid")]

demo <- merge(tracker, demo, by=c("bblid","scanid"))
demo <- merge(demo, envs, by=c("bblid"))

model <- "~ s(ageAtScan1, k=4) + sex + race2"
path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_antsCtVol.csv"
cores = 1

#idemo <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/idemo//n1601_idemo_jlf_intersect_roivals.csv")
#idemo <- idemo[,c(2,3,grep(pattern = "_cope1_Task", names(idemo)))]
#write.csv(idemo, "/home/agarza/envsMeduAnalysis/n1601_idemo_jlfIntersect_Task.csv", row.names=F)

FunctionAnalyze <- function(model, Path, cores) {
  
  dataSubj <- read.csv(path)
  
  demo <- merge(dataSubj, demo, by=c("bblid","scanid"))
  
  dataSubj <- dataSubj[, - which(names(dataSubj) == "bblid")]
  dataSubj <- dataSubj[, - which(names(dataSubj) == "scanid")]
  
  output <- as.data.frame(matrix(NA, nrow = dim(dataSubj)[2], ncol=7))
  names(output) <- c("names","t.medu","p.medu","pfdr.medu","t.envs","p.envs","pfdr.envs")
  
  
  ##Do MEDU
  medu.model <- paste(model, "+ medu1")
  m.medu <- mclapply(1:dim(dataSubj)[2], function(x) {as.formula(paste(paste0("dataSubj[,",x,"]"), medu.model, sep=""))})
  
  model.medu <- mclapply(m.medu, function(x) {
    foo <- summary(gam(formula = x, data=demo, method="REML"))
  }, mc.cores=cores)
  
  val.medu <- t(mcmapply(function(x) {
    x$p.table[dim(x$p.table)[1], 3:4]
  }, model.medu))
  
  output[,2:3] <- val.medu
  
  #DO ENVS
  envs.model <- paste(model, "+ envSES")
  m.envs <- mclapply(1:dim(dataSubj)[2], function(x) {as.formula(paste(paste0("dataSubj[,",x,"]"), envs.model, sep=""))})
  
  model.envs <- mclapply(m.envs, function(x) {
    foo <- summary(gam(formula = x, data=demo, method="REML"))
  }, mc.cores=cores)
  
  val.envs <- t(mcmapply(function(x) {
    x$p.table[dim(x$p.table)[1], 3:4]
  }, model.envs))
  
  output[,5:6] <- val.envs
  
  output$names <- names(dataSubj)
  output$pfdr.medu <- p.adjust(output$p.medu, "fdr")
  output$pfdr.envs <- p.adjust(output$p.envs, "fdr")
  
  return(output)
}

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionVol.csv"
outputVol <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionCt.csv"
outputCt <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAtroposIntersectionGMD.csv"
outputgmd <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/asl/n1601_jlfAntsCTIntersectionPcaslValues.csv"
outputcbf <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/dti/n1601_JHUTractFA.csv"
outputfa <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/dti/n1601_JHUTractTR.csv"
outputmd <- FunctionAnalyze(model, path, cores)

path <- "/home/agarza/envsMeduAnalysis/n1601_idemo_jlfIntersect_Task.csv"
outputidemo <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/nback//n1601_jlfIntersectNbackValues.csv"
outputnback <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/rest/n1601_jlfAntsCTIntersectionAlff.csv"
outputAlff <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/rest/n1601_jlfAntsCTIntersectionReHo.csv"
outputreho <- FunctionAnalyze(model, path, cores)

tmp <- list(outputAlff,outputcbf, outputCt,outputfa, outputgmd, outputidemo, outputmd,outputnback, outputreho, outputVol)

finalOutput <- do.call(rbind, tmp)
write.csv(finalOutput, "/home/agarza/envsMeduAnalysis/n1601_wholepnc_envs_medu.csv", row.names=F)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct//n1601_freesurferSurfaceArea.csv"
outputsurfarea <- FunctionAnalyze(model, path, cores)

