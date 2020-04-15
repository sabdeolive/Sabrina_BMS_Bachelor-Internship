#Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

#Install Biobase
BiocManager::install("Biobase")

#Assigning proteome data to proteome and creating matrix
proteome <- proteome_abundance.9
class(proteome)
dim(proteome)
proteomeMat <- data.matrix(proteome)
dim(proteomeMat)
class(proteomeMat)

#Assigning subject data to subject
subject <-Subject.data.T2DM.iHMP.5
class(subject)
dim(subject)

#Load package
library(Biobase)

#Creating ExpressionSet object
eset <- ExpressionSet(assayData = proteomeMat,
                      phenoData= AnnotatedDataFrame(subject))
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

fit

#Extracting top ranked proteins from fit
TableProt <- topTable(fit,coef=2)

#Export TableProt
write.table(TableProt,
            file = "c:/Users/sabde/Documents/TopRankedProteinsFromR.txt", 
            sep="\t")