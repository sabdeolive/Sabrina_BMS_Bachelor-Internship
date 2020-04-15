#Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

#Install Biobase
BiocManager::install("Biobase")

#Assigning metabolome data to metabolome and creating matrix
metabolome <- metabolome_abundance.11
class(metabolome)
dim(metabolome)
metabolomeMat <- data.matrix(metabolome)
dim(metabolomeMat)
class(metabolomeMat)

#Assigning subject data to subject
subject <-Subject.data.T2DM.iHMP.5
class(subject)
dim(subject)

#Load package
library(Biobase)

#Creating ExpressionSet object
eset <- ExpressionSet(assayData = metabolomeMat,
                      phenoData = AnnotatedDataFrame(subject))
dim(eset)

#Install limma
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

#Translate statistical model to R code (creation of design matrix)
design<- model.matrix(~IR_IS_classification, data = pData(eset))
head(design,2)
colSums(design)

#Load limma
library(limma)

#Fit the model with lmFit 
fit<- lmFit(eset,design)

#Calculate the t-statistics 
fit<- eBayes(fit)

#Summarize results
results<- decideTests(fit[,"IR_IS_classificationIS"])
summary(results)