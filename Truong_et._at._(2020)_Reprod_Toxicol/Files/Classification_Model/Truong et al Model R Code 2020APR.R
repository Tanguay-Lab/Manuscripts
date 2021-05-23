###The multi-dimensional embryonic zebrafish platform predicts flame retardant neurotoxicity
## Lisa. Truong
### May 2020 

rm(list=ls())
output.prefix <- toupper(format(Sys.Date(),"%Y%b%d"))
library(plyr); library(reshape2); library(ggplot2);
library(xlsx);
library(dplyr); library(car)
library(caret); library(mlbench);
library(skimr)

##run for all folders in the "results" directory
file.dir <-####list the directory here
setwd(file.dir)

##read in the raw data. 
devtox <- read.xlsx("DevTox Data.xlsx", sheetIndex = 1, check.names=F)
EPR <- read.xlsx("DevTox Data.xlsx", sheetIndex = 2, check.names=F)
LPR <- read.xlsx("DevTox Data.xlsx", sheetIndex = 3, check.names=F)
physchem <- read.xlsx("DevTox Data.xlsx", sheetIndex = 4, check.names=F)
class <- read.xlsx("DevTox Data.xlsx", sheetIndex = 5, check.names=F)
devtox.bmd <- read.xlsx("DevTox Data.xlsx", sheetIndex = 6, check.names=F)
key <- read.xlsx("DevTox Data.xlsx", sheetIndex = 7, check.names=F)
lpr.auc <- read.xlsx("DevTox Data.xlsx", sheetIndex = 8, check.names=F)
lpr.mov <- read.xlsx("DevTox Data.xlsx", sheetIndex = 9, check.names=F)

##combine it all togethers
zf.all <- class %>% 
  join(EPR) %>% 
  join(LPR) %>%
  join(physchem) %>%
  join(devtox.bmd) %>%
  join(devtox) %>%
  join(lpr.auc) %>%
  join(lpr.mov) %>%
  rename(EPR.LEC = CONC.EPR,
         LPR.AUC.LEC = CONC.LPR,
         Morph.LEC = any.effect.LEL,
         Morph.BMD = any.effect.BMD)

head(zf.all)


###set up modeling functions
classification <- "Class"
phys.property <- c("LogP",
                   "BoilingPt", "Density", "EnthalpyVaporization", "FlashPt",
                   "IndexRefraction", "MolarRefractivity", "MW", "NumFreeRotatingBonds", "NumHbondacceptor",
                   "NumHbonddonors", "NumRule5Violations", "VapourPressure")
zf.lec <- c("Morph.LEC", "EPR.LEC", "LPR.AUC.LEC")
zf.bmd <- c("Morph.BMD")#,"LPR.MOV.BMD", "LPR.AUC.BMD")#


##a list of the 3 models to be built
base.model <- c(classification, phys.property)
lec.model <- c(base.model, zf.lec)
bmd.model <- c(base.model, zf.bmd)
all.model <- c(base.model, zf.lec, zf.bmd)
#########################################################
##create two sets - test and train
#########################################################
row.names(zf.all) <- zf.all[,1]

##change everything to be numeric
zf.all[,c(3:ncol(zf.all))] <- sapply(zf.all[,c(3:ncol(zf.all))], function(x) as.numeric(as.character(x)))


####list the model 
model <- "bmd_mov" ### change the prefix
model.choice <- bmd.model ##the model to run. ##see lines 56-59
zf <- zf.all %>% select(all_of(model.choice)) 
head(zf)
colnames(zf)
dim(zf)


set.seed(998)

##change those with NA to a million
inTrain <- createDataPartition(y=zf$Class, p=.75, list=F)

##create the train and test set
training <- zf[inTrain,]
training1 <- training
training1$Chemical.ID <- row.names(training1)

testing <- zf[-inTrain,]
testing1 <- testing
testing1$Chemical.ID <- row.names(testing1)


ctrl2 <- trainControl(method="repeatedcv",
                      classProbs = TRUE,
                      number = 10,
                      savePredictions = 'final',
                      repeats = 5) ##cross validation, and bootstrap


sum <- data.frame()

############################################
## predict class of FRC
################################################
##algorithm spot check
metric <- "Accuracy"
# Random Forest
fit.rf <- train(Class~., data=training, method="rf", metric=metric, trControl=ctrl2)


# Table comparison
results <- resamples(list(rf=fit.rf))

summary(results)


# boxplot comparison
bwplot(results)
# Dot-plot comparison
dotplot(results)

###############################################################################################################
###now for the best models..
###rf
set.seed(999)
fitcontrol <- trainControl(
  method = 'cv',
  number=5,
  savePredictions = 'final',
  classProbs = T,
  summaryFunction = twoClassSummary
)

modelfit.rf <- train(Class ~ .,
                   data=training,
                   method="rf",
                   tuneLength=15,
                   trControl=ctrl2,
                   metric="ROC",
                   preProc=c("center", "scale"))
summary(modelfit.rf)
plot(modelfit.rf)

mClass.rf <- predict(modelfit.rf, newdata=testing)
head(mClass.rf)

cm.rf <- confusionMatrix(data=mClass.rf, testing$Class)
cm.rf

temp <- t(cm.rf$overall)
row.names(temp) <- "rf"
sum <- rbind(sum, temp)

dat1.1 <- data.frame(t(cm.rf$byClass))
dat1.1$desc <- row.names(dat1.1)



###play around model so don't screw up the real one
modelfit.rf5 <- train(Class ~ .,
                     data=training,
                     method="rf",
                     tuneLength=15,
                     trControl=ctrl2,
                     metric="ROC",
                     importance=TRUE,
                     preProc=c("center", "scale"))

##important variables
var_rf <- varImp(modelfit.rf5, scale=F)$importance
pdf(paste0("FRC_RF_Model_ImpVar_Needle_Plots_", model, "_", output.prefix, ".pdf"), width=11, height=8) ##prints out a PDF
plot(varImp(modelfit.rf5),  main="Important Variables")
dev.off()


modelfit.rf5
plot(modelfit.rf5)

mClass.rf5 <- predict(modelfit.rf5, newdata=testing)

cm.rf5 <- confusionMatrix(data=mClass.rf5, testing$Class)
cm.rf5


###writing out all the files
wb<-createWorkbook(type="xlsx")
datasheet <- createSheet(wb, sheetName = "raw_data")
addDataFrame(zf, datasheet)

sheet <- createSheet(wb, sheetName = "Overall_summary")
addDataFrame(sum, sheet)

rfsheet <- createSheet(wb, sheetName = "rf")
addDataFrame(dat1.1, rfsheet)

trainingsheet <- createSheet(wb, sheetName = "trainingset")
addDataFrame(training1, trainingsheet)

testingsheet <- createSheet(wb, sheetName = "testset")
addDataFrame(testing1, testingsheet)

rf_var <- createSheet(wb, sheetName = "rf_varimp")
addDataFrame(var_rf, rf_var)

saveWorkbook(wb, paste0("FRC_Model_Searching_Results_", model, "_", output.prefix,".xlsx"))

save(modelfit.rf,file=paste0("FRC_Toxicity_Models_", model, "_", output.prefix,".RData"))

###############################
## data visualization
data <- cm.rf$byClass
data.BA <- as.numeric(data[,11])
data.BA <- data.frame(data.BA)
data.BA$Class <- c("APE", "BP", "CPE", "Other", "OtherBrominated", "PBDE")
data.BA$Sensitivity <- as.numeric(data[,1])
data.BA$PPV <- as.numeric(data[,3])
colnames(data.BA)[1] <- c("Balance.Accuracy")
data.BA$Class <- as.factor(data.BA$Class)
data.BA$Balance.Accuracy <- data.BA$Balance.Accuracy * 100

data.BA.melt <- melt(data.BA)
data.BA.melt$value <- data.BA.melt$value*100



pdf(paste0("RF_Model_3_Parameter_results_", model, "_", output.prefix, ".pdf"), width=11, height=8) ###PDF

print(ggplot(data.BA, aes(Class,Balance.Accuracy, fill=Class)) + 
        geom_bar(stat = "identity") +
        ylab("Balance Accuracy")+
        scale_fill_grey() + theme_classic()+
        theme(axis.title.y = element_text(size=13, face="bold"),
              axis.title.x = element_text(size=13, face="bold"),
              axis.text	= element_text(size=10))+
        theme(strip.text.y = element_text(size= 14, face="bold"))
)

dev.off()
