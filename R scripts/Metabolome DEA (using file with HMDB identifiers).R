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
DATA.DIR <- "C:/Users/sabde/"
setwd(DATA.DIR)
getwd()

# Read the metabolome data
metabolome <- read.delim("metabolome_abundance 13.txt", header = TRUE, sep = "\t", row.names = NULL)
dim(metabolome)
df = data.frame(matrix(1:323))
df
rownames(df) = make.names(metabolome[,1],unique=TRUE)
#
rownames(metabolome) <- df[,1]
metabolome <- metabolome[,2:61]
#metabolome <- metabolome[,-1]
dim(metabolome)
class(metabolome)
metabolomeMat <- data.matrix(metabolome)
dim(metabolomeMat)
class(metabolomeMat)

#load study description
desc <- read.delim("Subject data T2DM iHMP 5.txt", as.is=TRUE)
dim(desc)

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
