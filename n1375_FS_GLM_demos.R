library(parallel)
library(mgcv)
library(e1071)
library(ggplot2)

demo <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_demographics_go1_20161212.csv")
envs <- read.csv("/data/joy/BBL/studies/pnc/n9498_dataFreeze/environment/n9498_go1_environment_factor_scores_tymoore_20150909.csv")
tracker <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_tracker_go1_20161212.csv")
tracker <- tracker[c("bblid","scanid")]
clinical <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/clinical/n1601_goassess_itemwise_bifactor_scores_20161219.csv")
datexscanid <- read.csv("/data/joy/BBL/projects/envsMeduAnalysis/envsAnalysisPNC/n1375_bblid_datexscanid.csv")

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

demo$sex <- as.factor(demo$sex)
demo$race2[which(demo$race2 == 3)] <- 2
demo$race2 <- as.factor(demo$race2)

demo <- subset(demo, select=c("bblid","scanid","sex","race2","ageAtScan1","medu1","fedu1","envSES","averageManualRating"))
names(demo)[4] <- paste("race")

demo$age <- demo$ageAtScan1 / 12
demo$ageSqrd <- (demo$age - mean(demo$age))^2

demo <- demo[!is.na(demo$medu1) | !is.na(demo$fedu1),]

## subset those with missing medu1
meduMissing <- demo[is.na(demo$medu1),]
meduMissing <- subset(meduMissing, select=c("bblid","scanid","fedu1"))
meduMissing$pedu <- (meduMissing$fedu1)/2

## subset those with missing fedu1
feduMissing <- demo[is.na(demo$fedu1),]
feduMissing <- subset(feduMissing, select=c("bblid","scanid","medu1"))
feduMissing$pedu <- (feduMissing$medu1)/2

## subset those with both medu1 and fedu1
both <- demo[complete.cases(demo), ]
both <- subset(both,select=c("bblid","scanid","medu1","fedu1"))
both$pedu <- (both$medu1 + both$fedu1)/2

# recombine for final demos
avgMiss <- merge(feduMissing, meduMissing, by=c("bblid","scanid","pedu"),all=TRUE)
avg <- merge(both,avgMiss,by=c("bblid","scanid","medu1","fedu1","pedu"),all=TRUE)
fullDemos <- merge(demo,avg, by=c("bblid","scanid","medu1","fedu1"),all=TRUE,row.names=FALSE)
final <- merge(datexscanid, fullDemos,by=c("bblid"))
write.csv(final,"/data/joy/BBL/projects/envsMeduAnalysis/envsAnalysisPNC/n1375_envs_demos.csv",row.names=FALSE)
