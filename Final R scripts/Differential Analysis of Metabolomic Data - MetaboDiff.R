######################################################################
# This script contains the code to perform the differential analysis of the iHMP T2D metabolomic data using MetaboDiff
######################################################################

### Install necessary packages 
# NOTE: have to have downloaded RTools for R version 3.6.3 prior to this 

# install.packages("WGCNA")
# install.packages("devtools")
# devtools::install_github("andreasmock/MetaboDiff") # Delete ellipsis in OneDrive R folder if necessary and Click YES.

### Load necessary packages
library("devtools")
library(MetaboDiff)

### Load the data
# NOTE: make first column of all datasets the rownames
assay <- metabolome_abundance.11
rowData <- `iPOP_Metablolite_Annotation.(excl..metabolites.without.HMDB.identifier).B`
colData <- Subject.data.T2DM.iHMP.4..2
colData <- colData[,-8]

### Merge all objects into a MultiAssayExperiment object
met <- create_mae(assay, rowData, colData)

### Preprocess data 
## Annotate metabolites 
met <- get_SMPDBanno(met,
                     column_kegg_id=2,
                     column_hmdb_id=3,
                     column_chebi_id=NA)

## Define cut-off for raw metabolite measurements
met <- knn_impute(met,cutoff=0.4) # i.e. only keeping those with raw measurements in more than 60% of cases

## Check if there are outliers and exclude if necessary
outlier_heatmap(met,
                group_factor="IR_IS_classification",
                label_colors=c("darkseagreen","dodgerblue"),
                k=2)

## Normalization (variance stabilizing normalization)
met.nor <- normalize_met(met)

### Hypothesis testing w/ correction for multiple testing (applies a Student's T-Test since there are only 2 groups and p-values are corrected using Benjamini-Hochberg procedure)
# (testing of the hypothesis: do the individual metabolites show differential abundance between the 2 groups?)
met.test <- diff_test(met.nor,
                      group_factors = "IR_IS_classification")

### Extracting the results
DA.results <- metadata(met.test) 
# order of the metabolites is the same as in metabolome_abundance.11



