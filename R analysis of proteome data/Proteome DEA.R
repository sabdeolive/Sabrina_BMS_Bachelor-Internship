#Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

#Install Biobase
BiocManager::install("Biobase")

#Assigning proteome data to variable p
prot <- proteome_abundance_8
dim(prot)

#Assigning subject data to variable s
subj <-Subject_data_T2DM_iHMP_4[,c("SubjectID","IR_IS_classification")]
dim(subj)

#Load package
library(Biobase)

#Creating ExpressionSet object
eset <- ExpressionSet(assayData = prot,
                      phenoData= AnnotatedDataFrame(subj))

#Install limma
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")