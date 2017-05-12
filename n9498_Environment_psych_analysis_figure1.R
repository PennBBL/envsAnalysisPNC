##############################################################################
################                                               ###############
################             Env Analysis Figure 1             ###############
################           Angel Garcia de la Garza            ###############
################              angelgar@upenn.edu               ###############
################                 05/10/2016                    ###############
##############################################################################




##############################################################################
################     Load libraries and data                   ###############
##############################################################################

library(parallel)
library(mgcv)
library(e1071)
library(ggplot2)

#load Data
demo <- read.csv("/data/joy/BBL/studies/pnc/n9498_dataFreeze//demographics//n9498_demographics_go1_20161212.csv")
envs <- read.csv("/data/joy/BBL/studies/pnc/n9498_dataFreeze//environment//n9498_go1_environment_factor_scores_tymoore_20150909.csv")
tracker <- read.csv("/data/joy/BBL/studies/pnc/n9498_dataFreeze//demographics//n9498_demographics_go1_20161212.csv")
tracker <- tracker[c("bblid")]
clinical <- read.csv("/data/joy/BBL/studies/pnc/n9498_dataFreeze//clinical//n9498_goassess_itemwise_bifactor_scores_20161219.csv")
health <- read.csv("/data/joy/BBL/studies/pnc/n9498_dataFreeze/health/n9498_health_20170405.csv")
health <- health[which(health$medicalrating < 3),]

#merge Data
demo <- merge(tracker, demo, by=c("bblid"))
demo <- merge(demo, envs, by=c("bblid"))
demo <- merge(demo, clinical, by=c("bblid"))
demo <- merge(demo, health, by="bblid")

#declare number of cores (leave this to one)
cores <- 1

#Clean Data
demo$sex <- as.factor(demo$sex)
demo$race2[which(demo$race2 == 3)] <- 2
demo$race2 <- as.factor(demo$race2)



##############################################################################
################     Function to analyze data                  ###############
##############################################################################


FunctionAnalyze <- function(model, Path, cores) {
  
  dataSubj <- read.csv(path)
  
  demo.ids <- demo[c("bblid")]
  
  dataSubj <- merge(demo.ids, dataSubj, by=c("bblid"))
  
  dataSubj <- dataSubj[, - which(names(dataSubj) == "bblid")]
  
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

##############################################################################
################  Analyze Bifactor and CorrTraits for figures  ###############
##############################################################################



#Declare model
model <- "~ s(ageAtClinicalAssess1, k=4) + sex + race2"


#Analyze bifactor scores 
path <- "/data/joy/BBL/studies/pnc/n9498_dataFreeze//clinical//n9498_goassess_itemwise_bifactor_scores_20161219.csv"
outputbifactor <- FunctionAnalyze(model, path, cores)

path <- "/data/joy/BBL/studies/pnc/n9498_dataFreeze//clinical//n9498_goassess_itemwise_corrtraits_scores_20161219.csv"
outputcorrtraits <- FunctionAnalyze(model, path, cores)

#Clean Names of output statistics (parameters from model)
outputbifactor$names <- c("Anxious-Misery","Psychosis","Behavioral","Fear","Overall Psychopathology")

outputbifactor$names <- factor(outputbifactor$names, levels = c("Overall Psychopathology","Anxious-Misery","Psychosis","Behavioral","Fear"))


##############################################################################
################  		Figure left panel c  	       ###############
##############################################################################

png("~/envsMeduAnalysis/n9498Analysis/psycho_envs.png", res=400, width=6.7, height=6.7, units="in")

ggplot(data=outputbifactor, aes(x=names, y=t.envs)) +
  geom_bar(aes(fill=names), stat="identity", color="black")+
  theme_classic() + 
  theme(text = element_text(size=20), 
        axis.text.x = element_text(angle = 45, hjust = 1, size=20, vjust= 1),
        axis.text.y = element_text(size=20),
        axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5)) +
  ylab("Association of SES and Psychopathology \n (z-score)") + 
  ylim(-9,2.5) +
  annotate("text",x=1,y=-9, size=6, label="***") +
  annotate("text",x=3,y=-9, size=6, label="***") +
  annotate("text",x=5,y=-9, size=6, label="*") +
  annotate("text",x=4,y=-9, size=6, label="***") + geom_hline(yintercept=0) + 
  scale_fill_manual(name = "Diagnosis",
                    values = c("darkgreen","blue", "darkorchid4","red2","gold2"),
                    guide=FALSE) +
  scale_x_discrete(labels=c("Overall \nPsychopathology","Anxious-Misery","Psychosis","Behavioral","Fear"))

  
dev.off()

##############################################################################
################  Corrtrait figure; not used in analysis       ###############
##############################################################################


png("~/envsMeduAnalysis/n9498Analysis/corrtrait_envs.png", res=400, width=6.7, height=6.7, units="in")

ggplot(data=outputcorrtraits, aes(x=names, y=t.envs)) +
  geom_bar(colour="black", stat="identity", fill="gray75")+
  theme_classic() + 
  theme(text = element_text(size=20), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x=element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("envSES effect in CorrTrait Scores", ) + 
  geom_hline(aes(yintercept =  - qnorm(1 - (.025)/dim(outputcorrtraits)[1], mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)))

dev.off()


##############################################################################
################  Creating diagnosis accross sample  	       ###############
##############################################################################


#Create diagnosis accross sample
#This was borrowed from Toni
diag <- read.csv("/data/joy/BBL/studies/pnc/n9498_dataFreeze/clinical/n9498_goassess_psych_summary_vars_20131014.csv")
diagPS <- read.csv("/data/joy/BBL/studies/pnc/n9498_dataFreeze/clinical/n9498_diagnosis_dxpmr_20161014.csv")

diag <- merge(diag, diagPS, by="bblid")

##Make variables where 1 = diagnosis.

#ADHD
diag$Add <- NA
diag$Add[which(diag$smry_add==4)] <- 1

#Agoraphobia
diag$Agr <- NA
diag$Agr[which(diag$smry_agr==4)] <- 1

#Anorexia
diag$Ano <- NA
diag$Ano[which(diag$smry_ano==4)] <- 1

#Bulimia
diag$Bul <- NA
diag$Bul[which(diag$smry_bul==4)] <- 1

#Conduct Disorder
diag$Con <- NA
diag$Con[which(diag$smry_con==4)] <- 1

#Generalized Anxiety Disorder
diag$Gad <- NA
diag$Gad[which(diag$smry_gad==4)] <- 1

#Mania
diag$Man <- NA
diag$Man[which(diag$smry_man==4)] <- 1

#Major Depressive Disorder
diag$Mdd <- NA
diag$Mdd[which(diag$smry_dep==4)] <- 1

#OCD
diag$Ocd <- NA
diag$Ocd[which(diag$smry_ocd==4)] <- 1

#Oppositional Defiant Disorder
diag$Odd <- NA
diag$Odd[which(diag$smry_odd==4)] <- 1

#Panic Disorder
diag$Pan <- NA
diag$Pan[which(diag$smry_pan==4)] <- 1

#Psychosis
diag$Ps <- NA
diag$Ps[which(diag$goassessDxpmr6 == "PS")] <- 1

#Posttraumatic Stress Disorder
diag$Ptd <- NA
diag$Ptd[which(diag$smry_ptd==4)] <- 1

#Separation Anxiety Disorder
diag$Sep <- NA
diag$Sep[which(diag$smry_sep==4)] <- 1

#Social Anxiety Disorder
diag$Soc <- NA
diag$Soc[which(diag$smry_soc==4)] <- 1

#Specific Phobia
diag$Sph <- NA
diag$Sph[which(diag$smry_phb==4)] <- 1

#Typically Developing
dxNames <- c("bblid","Add","Agr","Ano","Bul","Con","Gad","Man","Mdd","Ocd","Odd","Pan","Ps","Ptd","Sep","Soc","Sph")
dxDf <- data.matrix(diag[,dxNames])
diag$totDx <- rowSums(dxDf[,2:17], na.rm=TRUE) #This is how many people have how many diagnoses: sum(totDx==0):414, sum(totDx==1):307, sum(totDx>=2):638
diag$Td <- NA
diag$Td[which(diag$totDx==0)] <- 1


##############################################################################
################  Nonregressed figure; not used in analysis    ###############
##############################################################################




library(psych)

demo <- merge(demo,diag, by="bblid")
describeBy(demo$envSES, demo$Td)


dxNames <- c("Agr","Con","Gad","Mdd","Ocd","Odd","Ps","Ptd","Sep","Soc","Sph","Td")

t <- describeBy(demo$envSES, demo$Add, mat = T)
t$group1 <- as.character(as.factor(t$group1))
t$group1[1] <- "Add"


for (j in dxNames) {
  i <- which(names(demo) == j)
  temp <- describeBy(demo$envSES, demo[,i], mat = T)
  temp$group1 <- as.character(as.factor(temp$group1))
  temp$group1[1] <- j
  
  t <- rbind(t, temp)  
}

t$group1 <- c("ADHD","Agoraphobia","Conduct","GAD","MDD","OCD","ODD","Psychosis","PTSD","Sep. Anxiety","Social Anxiety","Specific Phobias","TD")


png("~/envsMeduAnalysis/n9498Analysis/diagnosis_environment.png", res=400, width=13.4, height=6.7, units="in")

ggplot(t, aes(x=group1, y=mean)) + 
  geom_bar(position=position_dodge(), stat="identity",fill="orange", col="black") +
  geom_errorbar(aes(ymin=mean-2*se, ymax=mean+2*se), width=.1, position=position_dodge(.9)) + 
  ylab("Neighborhood-level SES Score") +
  xlab("Diagnosis") + 
  theme_classic() +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x=element_blank(), plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept=0)

dev.off()


##############################################################################
################  Figure A  				       ###############
##############################################################################

#Calculate regressed environment for each diagnosis group

model <- gam(envSES ~ s(ageAtClinicalAssess1, k=4) + sex + race2, data=demo)
demo$resienvSES <- model$residuals

t <- describeBy(demo$resienvSES, demo$Add, mat = T)
t$group1 <- as.character(as.factor(t$group1))
t$group1[1] <- "Add"


for (j in dxNames) {
  i <- which(names(demo) == j)
  temp <- describeBy(demo$resienvSES, demo[,i], mat = T)
  temp$group1 <- as.character(as.factor(temp$group1))
  temp$group1[1] <- j
  
  t <- rbind(t, temp)  
}

t$group1 <- c("ADHD","Agoraphobia","Conduct","GAD","MDD","OCD","ODD","Psychosis","PTSD","Sep. Anxiety","Social Anxiety","Specific Phobias","TD")

png("~/envsMeduAnalysis/n9498Analysis/diagnosis_environment_regressed.png", res=400, width=13.4, height=6.7, units="in")

ggplot(t, aes(x=group1, y=mean)) + 
  geom_bar(position=position_dodge(), stat="identity", aes(fill=mean), col="black") +
  geom_errorbar(aes(ymin=mean-2*se, ymax=mean+2*se), width=.1, position=position_dodge(.9)) + 
  ylab("Neighborhood-level SES Score") +
  xlab("Diagnosis") + 
  theme_classic() +
  theme(text = element_text(size=24), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x=element_blank(), plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept=0) +
  scale_fill_gradient(low = "skyblue", high = "skyblue4") + guides(fill=FALSE)

dev.off()


##############################################################################
################Non regressed envs vs psychopathology	       ###############
################     (IGNORE this is in the figure)  	       ###############
##############################################################################



model <- gam(overall_psychopathology_4factorv2 ~ s(ageAtClinicalAssess1, k=4) + sex + race2, data=demo)
demo$resiOverall <- model$residuals

cor.test(demo$resiOverall, demo$resienvSES)

anoText <- paste0("r(",cor.test(demo$resiOverall, demo$resienvSES)$parameter,
                  ")=", round(cor.test(demo$resiOverall, demo$resienvSES)$estimate, digits = 2),
                  ", p<0.001")

ggplot(demo, aes(x=overall_psychopathology_4factorv2, y=envSES)) +
  geom_point(shape=1, color="gray39",size=1) +    # Use hollow circles
  geom_smooth(method=lm, se=T, col="gray29") +
  ylab("Neighborhood-level SES Score") +
  xlab("Overall Psychopathology Score") + 
  annotate("text", x = 2, y = 2.2, size=6, label = anoText) + 
  theme_classic() +
  theme(text = element_text(size=24), plot.title = element_text(hjust = 0.5))

##############################################################################
################ Hex plot of overall vs. environment	       ###############
################     (IGNORE this is in the figure)  	       ###############
##############################################################################



png("~/envsMeduAnalysis/n9498Analysis/overall_environment_regressed_hex.png", res=400, width=6.7, height=6.7, units="in")


ggplot(demo, aes(x=resiOverall, y=resienvSES)) +
  stat_binhex(bins=50) +    # Use hollow circles
 
  geom_smooth(method=lm, se=T, col="gray19") +
  ylab("Neighborhood-level SES Score") +
  xlab("Overall Psychopathology Score") +  
  theme_classic() +
  theme(text = element_text(size=24), plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient(limits=c(0, 25), low="gray70", high="gray30", space="Lab") +
  guides(fill=F)

dev.off()


##############################################################################
################ Hex plot of overall vs. environment	       ###############
################     (IGNORE this is in the figure)  	       ###############
##############################################################################



png("~/envsMeduAnalysis/n9498Analysis/overall_environment_regressed.png", res=400, width=6.7, height=6.7, units="in")

ggplot(demo, aes(x=resiOverall, y=resienvSES)) +
  geom_point(shape=16, aes(color=resiOverall), size=1) +    # Use hollow circles
  geom_smooth(method=lm, se=T, fill="forestgreen", color="forestgreen") +
  ylab("Neighborhood-level SES Score") +
  xlab("Overall Psychopathology Score") + 
  theme_classic() +
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) + 
  scale_colour_gradient(low = "#99FF99", high = "#003300") + guides(color=FALSE)

dev.off()


##############################################################################
##########    Plot of effect of residualized envs and corrtraits   ###########
################     (IGNORE this is in the figure)  	       ###############
##############################################################################


path <- "/data/joy/BBL/studies/pnc/n9498_dataFreeze//clinical//n9498_goassess_itemwise_corrtraits_scores_20161219.csv"
corr <- read.csv(path)

demo <- merge(demo, corr, "bblid")

corrNames <- outputcorrtraits$names
outCorr <- outputcorrtraits

t <- 1
for (i in corrNames) {
  j <- which(names(demo) == i)
  model <- gam(demo[,j] ~ s(ageAtClinicalAssess1, k=4) + sex + race2, data=demo)
  res <- model$residuals
  outCorr$t.envs[t] <- cor.test(res,demo$resienvSES)$estimate
  t <- t + 1
}

png("~/envsMeduAnalysis/n9498Analysis/corrtrait_environment_regressed.png", res=400, width=6.7, height=6.7, units="in")

ggplot(data=outCorr, aes(x=names, y=t.envs)) +
  geom_bar(colour="black", stat="identity")+
  theme_classic() + 
  theme(text = element_text(size=24), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x=element_blank(), plot.title = element_text(hjust = 0.5)) +
  ylim(-.125, 0) + ylab("Partial Correlation") +
  annotate("text",x=3,y=-.12, size=6, label="All correlations are significant at P < 0.001")
  
dev.off()


