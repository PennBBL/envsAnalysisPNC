library(parallel)
library(mgcv)
library(e1071)
library(ggplot2)

demo <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_demographics_go1_20161212.csv")
envs <- read.csv("/data/joy/BBL/studies/pnc/n9498_dataFreeze/environment/n9498_go1_environment_factor_scores_tymoore_20150909.csv")
tracker <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_tracker_go1_20161212.csv")
tracker <- tracker[c("bblid","scanid")]
clinical <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/clinical/n1601_goassess_itemwise_bifactor_scores_20161219.csv")
cnb <- read.csv("/data/joy/BBL/studies/pnc/n9498_dataFreeze//cnb//n9498_cnb_factor_scores_fr_20170202.csv")

demo <- merge(tracker, demo, by=c("bblid","scanid"))
demo <- merge(demo, envs, by=c("bblid"))
demo <- merge(demo, clinical, by=c("bblid","scanid"))
demo <- merge(demo, cnb, by=c("bblid"))

healthExclude <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/health/n1601_health_20161214.csv")
T1Exclude <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_t1QaData_20170306.csv")
demo <- merge(demo, healthExclude, by=c("bblid","scanid"))
demo <- merge(demo, T1Exclude, by=c("bblid","scanid"))

area <- read.csv("/data//joy//BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionVol_20170412.csv")
area$mean_area <- rowMeans(area[,3:131])

demo <- merge(demo,area, by=c("bblid","scanid"))

demo <- demo[which(demo$healthExcludev2 == 0), ]
demo <- demo[which(demo$t1Exclude == 0), ]

t1 <- gam(overall_psychopathology_4factorv2 ~ s(ageAtScan1, k=4) + sex + race2 + averageManualRating + envSES, data=demo, method="REML")
t2 <- gam(overall_psychopathology_4factorv2 ~ s(ageAtScan1, k=4) + sex + race2 + averageManualRating + envSES + mean_area, data=demo, method="REML")

t3 <- gam(mean_area ~ s(ageAtScan1, k=4) + sex + race2 + averageManualRating + envSES, data=demo, method="REML")


alpha <- 227.48
se.alpha <- 21.533
beta <- -6.795e-05
se.beta <- 3.509e-05
SE <- sqrt((alpha^2)*(se.beta^2) + (beta^2)*(se.alpha^2))


plist <- 1

for (i in 17:21) {
  
  t2 <- gam(demo[,i] ~ s(ageAtScan1, k=4) + sex + race2 + averageManualRating + envSES + mean_area, data=demo, method="REML")
  beta <- summary(t2)$p.table[6,1]
  se.beta <- summary(t2)$p.table[6,2]  
  SE <- sqrt((alpha^2)*(se.beta^2) + (beta^2)*(se.alpha^2))
  sobel.stat <- (alpha)*(beta) /  SE
  pval <- pnorm(abs(sobel.stat), lower.tail=F) * 2
  plist <- c(plist, pval)
  
}


plist <- plist[-1]

pvalues <- cbind(names(demo)[17:21], plist, p.adjust(plist, "bonferroni"))
