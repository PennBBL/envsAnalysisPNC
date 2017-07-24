library(parallel)
library(mgcv)
library(e1071)

demo <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_demographics_go1_20161212.csv")
envs <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/environment/n9498_go1_environment_factor_scores_tymoore_20150909.csv")
tracker <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_tracker_go1_20161212.csv")
tracker <- tracker[c("bblid","scanid")]
volume <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_antsCtVol.csv")
clinical <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/clinical/n1601_goassess_itemwise_bifactor_scores_20161219.csv")

demo <- merge(tracker, demo, by=c("bblid","scanid"))
demo <- merge(demo, envs, by=c("bblid"))
demo <- merge(demo, volume, by=c("bblid","scanid"))
demo <- merge(demo, clinical, by=c("bblid","scanid"))

healthExclude <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/health/n1601_health_20161214.csv")
T1Exclude <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_t1QaData_v2.csv")
demo <- merge(demo, healthExclude, by=c("bblid","scanid"))
demo <- merge(demo, T1Exclude, by=c("bblid","scanid"))

cores <- 1

demo <- demo[which(demo$ltnExcludev2 == 0), ]
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

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionVol.csv"
outputVol <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionCt.csv"
outputCt <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAtroposIntersectionGMD.csv"
outputgmd <- FunctionAnalyze(model, path, cores)

model <- "~ s(ageAtScan1, k=4) + sex + race2 + mprage_antsCT_vol_TBV"

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionVol.csv"
outputVol2 <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionCt.csv"
outputCt2 <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAtroposIntersectionGMD.csv"
outputgmd2 <- FunctionAnalyze(model, path, cores)

model <- "~ s(ageAtScan1, k=4) + sex + race2 + mprage_antsCT_vol_TBV + averageManualRating"

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionVol.csv"
outputVol3 <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionCt.csv"
outputCt3 <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAtroposIntersectionGMD.csv"
outputgmd3 <- FunctionAnalyze(model, path, cores)

model <- "~ s(ageAtScan1, k=4) + sex + race2 + averageManualRating + medu1 + envSES"

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionVol.csv"
outputVol4 <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionCt.csv"
outputCt4 <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAtroposIntersectionGMD.csv"
outputgmd4 <- FunctionAnalyze(model, path, cores)

model <- "~ s(ageAtScan1, k=4) + sex + race2 + mprage_antsCT_vol_TBV + averageManualRating + medu1 + envSES + overall_psychopathology_4factorv2"

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionVol.csv"
outputVol5 <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionCt.csv"
outputCt5 <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAtroposIntersectionGMD.csv"
outputgmd5<- FunctionAnalyze(model, path, cores)

model <- "~ s(ageAtScan1, k=4) + sex + race2 + averageManualRating + medu1 + envSES + overall_psychopathology_4factorv2"

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionVol.csv"
outputVol6 <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionCt.csv"
outputCt6 <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAtroposIntersectionGMD.csv"
outputgmd6<- FunctionAnalyze(model, path, cores)


tmp <- list(outputCt,outputgmd,outputVol)
finalOutput <- do.call(rbind, tmp)
write.csv(finalOutput, "/home/agarza/envsMeduAnalysis/n1601_structural_envs_medu_manualrating.csv", row.names=F)

tmp <- list(outputCt2,outputgmd2,outputVol2)
finalOutput <- do.call(rbind, tmp)
write.csv(finalOutput, "/home/agarza/envsMeduAnalysis/n1601_structural_envs_medu_TBV.csv", row.names=F)

tmp <- list(outputCt3,outputgmd3,outputVol3)
finalOutput <- do.call(rbind, tmp)
write.csv(finalOutput, "/home/agarza/envsMeduAnalysis/n1601_structural_envs_medu_TBV_manualrating.csv", row.names=F)

tmp <- list(outputCt4,outputgmd4,outputVol4)
finalOutput <- do.call(rbind, tmp)
write.csv(finalOutput, "/home/agarza/envsMeduAnalysis/n1601_structural_envs_medu_on_model_TBV_manualrating.csv", row.names=F)

tmp <- list(outputCt5,outputgmd5,outputVol5)
finalOutput <- do.call(rbind, tmp)
write.csv(finalOutput, "/home/agarza/envsMeduAnalysis/n1601_structural_envs_medu_on_model_TBV_manualrating_OverallPsychopathology.csv", row.names=F)

tmp <- list(outputCt6,outputgmd6,outputVol6)
finalOutput <- do.call(rbind, tmp)
write.csv(finalOutput, "/home/agarza/envsMeduAnalysis/n1601_structural_envs_medu_on_model_manualrating_OverallPsychopathology.csv", row.names=F)


funFitModel <- function(model, Path, cores, name.output) {
  
  dataSubj <- read.csv(path)
  
  demo.ids <- demo[c("bblid","scanid")]
  
  dataSubj <- merge(demo.ids, dataSubj, by=c("bblid","scanid"))
  
  dataSubj <- dataSubj[, - which(names(dataSubj) == "bblid")]
  dataSubj <- dataSubj[, - which(names(dataSubj) == "scanid")]
  
  output <- as.data.frame(matrix(NA, nrow = dim(dataSubj)[2], ncol=4))
  names(output) <- c("names","t.val","p.val","pfdr.val")

  #Fit Model
  
  medu.model <- paste(model, "+ medu1")
  m.medu <- mclapply(1:dim(dataSubj)[2], function(x) {as.formula(paste(paste0("dataSubj[,",x,"]"), medu.model, sep=""))})
  
  model.medu <- mclapply(m.medu, function(x) {
    foo <- summary(gam(formula = x, data=demo, method="REML"))
  }, mc.cores=cores)
  
  val.medu <- t(mcmapply(function(x) {
    x$p.table[which(rownames(x$p.table) == name.output), 3:4]
  }, model.medu))
  
  output[,2:3] <- val.medu
    
  output$names <- names(dataSubj)
  output$pfdr.val <- p.adjust(output$p.val, "fdr")
  
  names(output) <- c("names",paste0("t.val.",name.output)
                     ,paste0("p.val.",name.output)
                     ,paste0("pfdr.val.",name.output))
  
  return(output)
}



model <- "~ s(ageAtScan1, k=4) + sex + race2 + averageManualRating   + envSES*sex"
name.output <- "sex2:envSES"

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionVol.csv"
outputVol7 <- funFitModel(model, path, cores, name.output)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionCt.csv"
outputCt7 <- funFitModel(model, path, cores, name.output)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAtroposIntersectionGMD.csv"
outputgmd7<- funFitModel(model, path, cores, name.output)

model <- "~  sex + race2 + averageManualRating + envSES  + envSES*ageAtScan1"
name.output <- "envSES:ageAtScan1"

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionVol.csv"
outputVol8 <- funFitModel(model, path, cores, name.output)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionCt.csv"
outputCt8 <- funFitModel(model, path, cores, name.output)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAtroposIntersectionGMD.csv"
outputgmd8<- funFitModel(model, path, cores, name.output)

model <- "~  sex + race2 + averageManualRating + envSES  + envSES*ageAtScan1"
name.output <- "envSES"

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionVol.csv"
outputVol9 <- funFitModel(model, path, cores, name.output)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionCt.csv"
outputCt9 <- funFitModel(model, path, cores, name.output)

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAtroposIntersectionGMD.csv"
outputgmd9<- funFitModel(model, path, cores, name.output)

tmp <- list(outputCt7,outputgmd7,outputVol7)
finalOutput <- do.call(rbind, tmp)
write.csv(finalOutput, "/home/agarza/envsMeduAnalysis/n1601_structural_envs_sex_age_race_Rating_sexEnvs_OutputSexEnvs.csv", row.names=F)

tmp <- list(outputCt8,outputgmd8,outputVol8)
finalOutput <- do.call(rbind, tmp)
write.csv(finalOutput, "/home/agarza/envsMeduAnalysis/n1601_structural_envs_sex_age_race_Rating_ageEnvs_OutputageEnvs.csv", row.names=F)

tmp <- list(outputCt9,outputgmd9,outputVol9)
finalOutput <- do.call(rbind, tmp)
write.csv(finalOutput, "/home/agarza/envsMeduAnalysis/n1601_structural_envs_sex_age_race_Rating_ageEnvs_OutputEnvs.csv", row.names=F)

finalDat <- demo

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionVol.csv"
tmp <- read.csv(path)
finalDat <- merge(finalDat, tmp, by=c("bblid","scanid"))

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionCt.csv"
tmp <- read.csv(path)
finalDat <- merge(finalDat, tmp, by=c("bblid","scanid"))

path <- "/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAtroposIntersectionGMD.csv"
tmp <- read.csv(path)
finalDat <- merge(finalDat, tmp, by=c("bblid","scanid"))


#Function to run SVR

set.seed(1)

svr.10fold<-function(y,features){
  nfolds=10
  n<-length(y)
  fits<-rep(NA,n)
  rand.vec<-sample( ceiling(seq(0.0001,(nfolds-0.0001), length.out=n)) )
  for(fold in 1:nfolds){
    y.sub<-y[which(rand.vec != fold) ]
    features.sub<-features[which(rand.vec!=fold),]
    svr.fold<-svm(y=y.sub, x=features.sub, cross=1, kernel='linear')
    fits[which(rand.vec==fold)]<-predict(svr.fold,features[which(rand.vec==fold),])
    cat('Done fold', paste(fold, nfolds, sep='/'), '\n')
  }
  return(fits)
}

#Brain Age analysis
y.age <- finalDat$ageAtScan1
features <- finalDat[, grep(names(finalDat) , pattern = "mprage")]


fits.bage <- svr.10fold(y.age, features)
fits.bage <- fits.bage
#rm(fits)

cor.test(fits.bage, y.age)

finalDat$brainAge <- fits.bage
gamBrainAge <- gam(brainAge ~ s(ageAtScan1, k=4) + envSES, data=finalDat, method = "REML")
summary(gamBrainAge)

#Enviroment Variables
y.env <- finalDat$envSES
features <- finalDat[, grep(names(finalDat) , pattern = "mprage")]

fits.env <- svr.10fold(y.env, features)
cor.test(fits.env, y.env)
finalDat$fittedEnv <- fits.env
model <- lm(fittedEnv ~ envSES, data=finalDat)
summary(model)

library(ggplot2)
ggplot(finalDat, aes(x=fittedEnv, y=envSES)) +
  geom_point(shape=1) + geom_abline(slope = 1, intercept = 0) + ylim(-2.4,2.2) + xlim(-2.4,2.2) +xlab("Observed SES") + ylab("Fitted SES")


##Run CT
features <- finalDat[, grep(names(finalDat) , pattern = "jlf_ct")]
fits.env <- svr.10fold(y.env, features)
cor.test(fits.env, y.env)
finalDat$fittedEnv <- fits.env
model <- lm(fittedEnv ~ envSES, data=finalDat)
summary(model)

##Run Vol
features <- finalDat[, grep(names(finalDat) , pattern = "jlf_vol")]
fits.env <- svr.10fold(y.env, features)
cor.test(fits.env, y.env)
finalDat$fittedEnv <- fits.env
model <- lm(fittedEnv ~ envSES, data=finalDat)
summary(model)

cor <- 0
for (i in 1:100) {
  
  sample <- sample(1:1200, 1200)
  fits.env <- svr.10fold(y.env, features[sample, ])
  cor <- c(cor, cor.test(fits.env, y.env)$estimate)
  print(i)
}

mean(cor[-1])

##Run GMD
features <- finalDat[, grep(names(finalDat) , pattern = "jlf_gmd")]
fits.env <- svr.10fold(y.env, features)
cor.test(fits.env, y.env)
finalDat$fittedEnv <- fits.env
model <- lm(fittedEnv ~ envSES, data=finalDat)
summary(model)

library(glmnet)

glmmod <- glmnet(y=y.env, x=as.matrix(features), alpha=1, family="gaussian")
cv.glmmod <- cv.glmnet(y=y.env, x=as.matrix(features), alpha=1, family="gaussian")
best.lambda <- cv.glmmod$lambda.min

glmmod <- glmnet(y=y.env, x=as.matrix(features), alpha=1, family="gaussian", lambda=best.lambda)

