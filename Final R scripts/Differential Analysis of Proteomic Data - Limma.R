######################################################################
# This script contains the code to perform the differential analysis of the iHMP T2D proteomic data using limma
######################################################################
library(BiocManager)
library(Biobase)
library(limma)

### Load necessary data
## Proteome data
prot <- proteome_abundance9 
rownames(prot) <- prot[,1] 
prot <- prot[,-1]

## Subject description data
desc <- Subjectdata._2DM_iHMP5

### Check if the order of the subject and proteome data match
colnames(prot) == desc$SubjectID


### Create a factor of the IR and IS classification to use as the explanatory variable for the linear model
Insgroup <- factor(desc$IR_IS_classification, levels=c("IS","IR"))

### Convert this factor to a design matrix in order to use with lmfit
design <- model.matrix(~ Insgroup)
# check to see which classification is given value 1 
head(design,2)

### Fit the model with lmFit 
fit <- lmFit(prot,design)

### Calculate the t-statistics 
fit <- eBayes(fit)

### Extract a table of the results
results <- topTable(fit, number=20009, coef="InsgroupIR")

