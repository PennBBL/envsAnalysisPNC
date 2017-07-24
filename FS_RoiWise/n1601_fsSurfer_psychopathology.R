library(parallel)
library(mgcv)
library(e1071)
library(ggplot2)

demo <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_demographics_go1_20161212.csv")
envs <- read.csv("/data/joy/BBL/studies/pnc/n9498_dataFreeze/environment/n9498_go1_environment_factor_scores_tymoore_20150909.csv")
tracker <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_tracker_go1_20161212.csv")
tracker <- tracker[c("bblid","scanid")]
clinical <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/clinical/n1601_goassess_itemwise_bifactor_scores_20161219.csv")

demo <- merge(tracker, demo, by=c("bblid","scanid"))
demo <- merge(demo, envs, by=c("bblid"))
demo <- merge(demo, clinical, by=c("bblid","scanid"))

healthExclude <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/health/n1601_health_20161214.csv")
T1Exclude <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_t1QaData_20170306.csv")
demo <- merge(demo, healthExclude, by=c("bblid","scanid"))
demo <- merge(demo, T1Exclude, by=c("bblid","scanid"))

cores <- 1

demo <- demo[which(demo$healthExcludev2 == 0), ]
demo <- demo[which(demo$fsFinalExclude == 0), ]


gamROI <- function(demo, Path, ids, model, cores) {
  
  dataSubj <- read.csv(path)
  
  dataSubj <- merge(demo, dataSubj, by=ids)
  
  model.formula <- mclapply((dim(demo)[2] + 1):dim(dataSubj)[2], function(x) {as.formula(paste(paste0("dataSubj[,",x,"]"), model, sep=""))}, mc.cores=cores)
  
  m <- mclapply(model.formula, function(x) {
    foo <- summary(gam(formula = x, data=demo, method="REML"))
  }, mc.cores=cores)
  
  length.names.p <- length(rownames(m[[1]]$p.table))                 
  
  output <- as.data.frame(matrix(NA, nrow = length((dim(demo)[2] + 1):dim(dataSubj)[2]), ncol= (1+3*length.names.p)))
  names(output)[1] <- "names"
  
  for (i in 1:length.names.p) {
    dep.val <- rownames(m[[1]]$p.table)[i]
    names(output)[2 + (i-1)*3 ] <- paste0("tval.",dep.val)
    names(output)[3 + (i-1)*3 ] <- paste0("pval.",dep.val)
    names(output)[4 + (i-1)*3 ] <- paste0("pvalfdr.",dep.val)
    
    val.tp <- t(mcmapply(function(x) {
      x$p.table[which(rownames(x$p.table) == dep.val), 3:4]
    }, m, mc.cores=cores))
    
    output[,(2 + (i-1)*3):(3 + (i-1)*3)] <- val.tp
    output[,(4 + (i-1)*3)] <- p.adjust(output[,(3 + (i-1)*3)], "fdr")
  }
  
  output$names <- names(dataSubj)[(dim(demo)[2] + 1):dim(dataSubj)[2]]
  p.output <- output
  
  if (is.null(m[[1]]$s.table)) {
    
    return(p.output)
    
  } else {
    
    length.names.s <- length(rownames(m[[1]]$s.table))
    output <- as.data.frame(matrix(NA, nrow = length((dim(demo)[2] + 1):dim(dataSubj)[2]), ncol= (1+2*length.names.s)))
    
    names(output)[1] <- "names"
    
    for (i in 1:length.names.s) {
      
      dep.val <- rownames(m[[1]]$s.table)[i]
      names(output)[2 + (i-1)*2 ] <- paste0("pval.",dep.val)
      names(output)[3 + (i-1)*2 ] <- paste0("pvalfdr.",dep.val)
      
      val.tp <- mcmapply(function(x) {
        x$s.table[which(rownames(x$s.table) == dep.val), 4]
      }, m, mc.cores=cores)
      
      output[,(2 + (i-1)*2)] <- val.tp
      output[,(3 + (i-1)*2)] <- p.adjust(output[,(2 + (i-1)*2)], "fdr")
      output$names <- names(dataSubj)[(dim(demo)[2] + 1):dim(dataSubj)[2]]
      
      s.output <- output
      output <- merge(p.output, s.output, by="names")
      return(output)
      
    }
    
  }
}

model <- "~ s(ageAtScan1, k=4) + sex + race2 + averageManualRating + mood_4factorv2 + psychosis_4factorv2 + externalizing_4factorv2 + phobias_4factorv2 + overall_psychopathology_4factorv2"

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_freesurferSurfaceArea_20161220.csv"
tarea <- gamROI(demo, path, ids=c("bblid","scanid"), model, cores=1)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_freesurferVol_20161220.csv"
tvol <- gamROI(demo, path, ids=c("bblid","scanid"), model, cores=1)

write.csv(tarea, "~/envsMeduAnalysis/n1601_freesurfer_psychopathology_area.csv")
write.csv(tvol, "~/envsMeduAnalysis/n1601_freesurfer_psychopathology_volume.csv")
