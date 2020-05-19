#Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

#Install Biobase
BiocManager::install("Biobase")
BiocManager::install("limma")

#Load packages
library(Biobase)
library(limma)

#change working directory to one incl. necessary files
DATA.DIR <- "C:/Users/sabde/OneDrive/Documents/UM documents/BMS Year 3/Internship + Thesis/Sabrina_BMS_Bachelor-Internship/R analysis of metabolome data/"
setwd(DATA.DIR)
getwd()

# Read the metabolome data
metabolome <- metabolome_abundance.11
rownames(metabolome) <- metabolome[,1]
metabolome <- metabolome[,-1]

#load study description
desc <- Subject.data.T2DM.iHMP.5

#check order of description and normalised data columns
colnames(metabolome) == desc$SubjectID
# OK

#groups
# Two groups: IS and IR
#################################################################################################################


## VitD
InsgroupM <- factor(desc$IR_IS_classification, levels=c("IS","IR"))
class(InsgroupM)


#Translate statistical model to R code (creation of design matrix)
designM<- model.matrix(~ InsgroupM)
head(designM,2)
colSums(designM)


#Fit the model with lmFit 
fitM<- lmFit(metabolome,designM)

#Calculate the t-statistics 
fitM<- eBayes(fitM)

#Export results
resultsM <- topTable(fitM, number=20009, coef="InsgroupMIR")
write.table(resultsM, file="c:/Users/sabde/Documents/metabolome_limma_analysis.txt", sep="\t", row.names = TRUE, col.names = NA)


##############################################################################
#### With log transformation
metabolome.log <- log10(metabolome + 1)

#check order of description and normalised data columns
colnames(metabolome.log) == desc$SubjectID
# OK

#groups
# Two groups: IS and IR
#################################################################################################################


## VitD
InsgroupM <- factor(desc$IR_IS_classification, levels=c("IS","IR"))
class(InsgroupM)


#Translate statistical model to R code (creation of design matrix)
designM<- model.matrix(~ InsgroupM)
head(designM,2)
colSums(designM)


#Fit the model with lmFit 
fitM.log<- lmFit(metabolome.log,designM)

#Calculate the t-statistics 
fitM.log<- eBayes(fitM.log)

#Export results
resultsM.log <- topTable(fitM.log, number=20009, coef="InsgroupMIR")
write.table(resultsM.log, file="c:/Users/sabde/Documents/metabolome_limma_analysis.log.txt", sep="\t", row.names = TRUE, col.names = NA)

##############################################################################
#### w/ new log transformation
install.packages("metabolomics")
library(metabolomics)
