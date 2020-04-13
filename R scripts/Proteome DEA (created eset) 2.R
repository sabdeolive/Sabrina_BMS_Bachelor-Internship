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

#
