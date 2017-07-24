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

area <- read.csv("/data//joy//BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_freesurferSurfaceArea_20161220.csv")
area$mean_area <- rowMeans(area[,3:72])

demo <- merge(demo,area, by=c("bblid","scanid"))

demo <- demo[which(demo$healthExcludev2 == 0), ]
demo <- demo[which(demo$fsFinalExclude == 0), ]

t1 <- gam(NAR_Overall_Accuracy ~ s(ageAtScan1, k=4) + sex + race2 + averageManualRating + envSES, data=demo, method="REML")
t2 <- gam(NAR_Overall_Accuracy ~ s(ageAtScan1, k=4) + sex + race2 + averageManualRating + envSES +overall_psychopathology_4factorv2 + mean_area, data=demo, method="REML")

t3 <- gam(mean_area ~ s(ageAtScan1, k=4) + sex + race2 + averageManualRating +overall_psychopathology_4factorv2 + envSES, data=demo, method="REML")


alpha <- 133.56
se.alpha <- 12.90
beta <- 3.126e-04
se.beta <- 4.563e-05

plist <- 1

for (i in 22:34) {
  
  t2 <- gam(demo[,i] ~ s(ageAtScan1, k=4) + sex + race2 + averageManualRating + envSES + mean_area +overall_psychopathology_4factorv2, data=demo, method="REML")
  beta <- summary(t2)$p.table[6,1]
  se.beta <- summary(t2)$p.table[6,2]  
  SE <- sqrt((alpha^2)*(se.beta^2) + (beta^2)*(se.alpha^2))
  sobel.stat <- (alpha)*(beta) /  SE
  pval <- pnorm(sobel.stat, lower.tail=F) * 2
  plist <- c(plist, pval)                    
}


plist <- plist[-1]

pvalues <- cbind(names(demo)[22:34], plist, p.adjust(plist, "bonferroni"))


zlist <- 1

for (i in 22:34) {
  
  t2 <- gam(demo[,i] ~ s(ageAtScan1, k=4) + sex + race2 + averageManualRating + envSES + mean_area, data=demo, method="REML")
  beta <- summary(t2)$p.table[6,1]
  se.beta <- summary(t2)$p.table[6,2]  
  SE <- sqrt((alpha^2)*(se.beta^2) + (beta^2)*(se.alpha^2))
  sobel.stat <- (alpha)*(beta) /  SE
  zval <- sobel.stat
  zlist <- c(zlist, zval)                    
}


zlist <- zlist[-1]


zvalues <- as.data.frame(cbind(names(demo)[22:34], zlist))
names(zvalues) <- "Factors"

zvalues$Factors <- as.character(zvalues$Factors)


for (i in 1:13) {
  zvalues$Factors[i] <- strsplit(as.character(zvalues$Factors[i]), "NAR_")[[1]][2]
}

png("~/envsMeduAnalysis/n1601Analysis/n1601_mediation_SES_Area_Cognitive.png", res=400, width=12, height=9, units="in")
# Very basic bar graph
ggplot(data=zvalues, aes(x=Factors, y=zlist)) +
  geom_bar(colour="black", stat="identity")+
  theme_bw() + 
  theme(text = element_text(size=20), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5)) +
  ggtitle("Sobel Z-Score for effect of SES on cognition mediated by Cortical Area", ) +
  ylab("Sobel Z-score")

dev.off()
