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

DATA.DIR <- "C:/Users/Administrator/Documents/GitHub/Sabrina_BMS_Bachelor-Internship/R analysis of proteome data/"
setwd(DATA.DIR)
getwd()

# Read the proteome data
prot <- read.delim("proteome_abundance9.txt", as.is=TRUE)
rownames(prot) <- prot[,1]
prot <- prot[,-1]

#load study description
desc <- read.delim("Subjectdata _2DM_iHMP5.txt", as.is=TRUE)

#check order of description and normalised data columns
colnames(prot) == desc$SubjectID
# OK

#groups
# Two groups: IS and IR
#################################################################################################################


## VitD
Insgroup <- factor(desc$IR_IS_classification, levels=c("IS","IR"))



#Assigning proteome data to proteome and creating matrix
#proteome <- proteome_abundance.9
#class(proteome)
#dim(proteome)
#proteomeMat <- data.matrix(proteome)
#dim(proteomeMat)
#class(proteomeMat)

#Assigning subject data to subject
#subject <-Subject.data.T2DM.iHMP.5
#class(subject)
#dim(subject)



#Creating ExpressionSet object
#eset <- ExpressionSet(assayData = proteomeMat,
#                     phenoData= AnnotatedDataFrame(subject))
#dim(eset)


#Translate statistical model to R code (creation of design matrix)
design<- model.matrix(~ Insgroup)
head(design,2)
colSums(design)



#Fit the model with lmFit 
fit<- lmFit(prot,design)

#Calculate the t-statistics 
fit<- eBayes(fit)

#Summarize results
#results<- decideTests(fit[,"IR_IS_classificationIS"])
#summary(results)

results <- topTable(fit, number=20009, coef="InsgroupIR")
write.table(results, file="proteome_limma_analysis.txt", sep="\t", row.names = TRUE, col.names = NA)


fit

#Extracting top ranked proteins from fit
TableProt <- topTable(fit,coef=2)

#Export TableProt
write.table(TableProt,
            file = "c:/Users/sabde/Documents/TopRankedProteinsFromR.txt", 
            sep="\t")