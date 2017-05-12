library(parallel)
library(mgcv)
library(e1071)
library(ggplot2)

demo <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_demographics_go1_20161212.csv")
envs <- read.csv("/data/joy/BBL/studies/pnc/n9498_dataFreeze/environment/n9498_go1_environment_factor_scores_tymoore_20150909.csv")
tracker <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_tracker_go1_20161212.csv")
tracker <- tracker[c("bblid","scanid")]
volume <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_antsCtVol_20161006.csv")
clinical <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/clinical/n1601_goassess_itemwise_bifactor_scores_20161219.csv")

demo <- merge(tracker, demo, by=c("bblid","scanid"))
demo <- merge(demo, envs, by=c("bblid"))
demo <- merge(demo, volume, by=c("bblid","scanid"))
demo <- merge(demo, clinical, by=c("bblid","scanid"))

healthExclude <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/health/n1601_health_20161214.csv")
T1Exclude <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_t1QaData_20170306.csv")
demo <- merge(demo, healthExclude, by=c("bblid","scanid"))
demo <- merge(demo, T1Exclude, by=c("bblid","scanid"))

cores <- 1

demo <- demo[which(demo$healthExcludev2 == 0), ]
demo <- demo[which(demo$t1Exclude == 0), ]

demo$sex <- as.factor(demo$sex)
demo$race2[which(demo$race2 == 3)] <- 2
demo$race2 <- as.factor(demo$race2)



FunctionAnalyze <- function(model, Path, cores) {
  
  dataSubj <- read.csv(path)
  
  demo.ids <- demo[c("bblid","scanid")]
  
  dataSubj <- merge(demo.ids, dataSubj, by=c("bblid","scanid"))
  
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
    x$p.table[which(rownames(x$p.table) == "medu1"), 3:4]
  }, model.medu))
  
  output[,2:3] <- val.medu
  
  #DO ENVS
  envs.model <- paste(model, "+ envSES")
  m.envs <- mclapply(1:dim(dataSubj)[2], function(x) {as.formula(paste(paste0("dataSubj[,",x,"]"), envs.model, sep=""))})
  
  model.envs <- mclapply(m.envs, function(x) {
    foo <- summary(gam(formula = x, data=demo, method="REML"))
  }, mc.cores=cores)
  
  val.envs <- t(mcmapply(function(x) {
    x$p.table[which(rownames(x$p.table) == "envSES"), 3:4]
  }, model.envs))
  
  output[,5:6] <- val.envs
  
  output$names <- names(dataSubj)
  output$pfdr.medu <- p.adjust(output$p.medu, "fdr")
  output$pfdr.envs <- p.adjust(output$p.envs, "fdr")
  
  return(output)
}

model <- "~ s(ageAtScan1, k=4) + sex + race2 + averageManualRating"

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_freesurferSurfaceArea_20161220.csv"
outputarea <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_freesurferCt_20161220.csv"
outputct <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_freesurferVol_20161220.csv"
outputVol <- FunctionAnalyze(model, path, cores)




path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/clinical//n1601_goassess_itemwise_bifactor_scores_20161219.csv"
outputbifactor <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/clinical/n1601_goassess_itemwise_corrtraits_scores_20161219.csv"
outputcorrtraits <- FunctionAnalyze(model, path, cores)

write.csv(outputarea, "~/envsMeduAnalysis/n1601_FsArea_mainmodel.output.csv")
write.csv(outputbifactor, "~/envsMeduAnalysis/n1601_bifactor_mainmodel.output.csv")
write.csv(outputcorrtraits, "~/envsMeduAnalysis/n1601_corrtraits_mainmodel.output.csv")

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/cnb/n1601_cnb_factor_scores_tymoore_20151006.csv"
outputcnb <- FunctionAnalyze(model, path, cores)
write.csv(outputcnb, "~/envsMeduAnalysis/n1601_cnb_mainmodel.output.csv")

outputcnb <- outputcnb[1:13, ]

png("~/envsMeduAnalysis/envs_cnb.png", res=400, width=12, height=9, units="in")

ggplot(data=outputcnb, aes(x=names, y=t.envs)) +
  geom_bar(colour="black", stat="identity")+
  theme_bw() + 
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("envSES effect in CNB", ) + 
  geom_hline(aes(yintercept = qnorm(1 - (.025)/dim(outputcnb)[1], mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)))

dev.off()

png("~/envsMeduAnalysis/psycho_cnb.png", res=400, width=12, height=9, units="in")


ggplot(data=outputbifactor, aes(x=names, y=t.envs)) +
  geom_bar(colour="black", stat="identity")+
  theme_bw() + 
  theme(text = element_text(size=20), 
          axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5)) +
  ggtitle("envSES effect in Psychopathology Bifactor Scores", ) +
  geom_hline(aes(yintercept = qnorm(1 - (.025)/dim(outputbifactor)[1], mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE))) +
  geom_hline(aes(yintercept =  - qnorm(1 - (.025)/dim(outputbifactor)[1], mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)))

dev.off()

png("~/envsMeduAnalysis/corrtrait_cnb.png", res=400, width=12, height=9, units="in")

ggplot(data=outputcorrtraits, aes(x=names, y=t.envs)) +
  geom_bar(colour="black", stat="identity")+
  theme_bw() + 
  theme(text = element_text(size=20), 
          axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x=element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("envSES effect in CorrTrait Scores", ) + 
  geom_hline(aes(yintercept =  - qnorm(1 - (.025)/dim(outputcorrtraits)[1], mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)))

dev.off()
