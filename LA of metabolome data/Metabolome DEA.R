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

library(ggplot2)
qplot(colSums(metabolome))
##############################################################################
#### With log transformation
metabolome.log <- log10(metabolome + 1)
class(metabolome.log)
dim(metabolome.log)
metabolomeMat.log <- data.matrix(metabolome.log)
dim(metabolomeMat.log)
class(metabolomeMat.log)

#Assigning subject data to subject
subject <-Subject.data.T2DM.iHMP.5
class(subject)
dim(subject)

#Load package
library(Biobase)

#Creating ExpressionSet object
eset.log <- ExpressionSet(assayData = metabolomeMat.log,
                      phenoData = AnnotatedDataFrame(subject))
# Error in validObject(.Object) : 
#   invalid class "ExpressionSet" object: 1: sampleNames differ between assayData and phenoData
# invalid class "ExpressionSet" object: 2: sampleNames differ between phenoData and protocolData

dim(eset.log)

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
fit.log<- lmFit(eset.log,design)

#Calculate the t-statistics 
fit.log<- eBayes(fit.log)

#Summarize results
results.log<- decideTests(fit.log[,"IR_IS_classificationIS"])
summary(results.log)

