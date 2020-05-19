# NOTE: have to have downloaded RTools for R version 3.6.3 prior to this 

#### Installation and loading of MetaboDiff
install.packages("WGCNA")

install.packages("devtools")
library("devtools")
devtools::install_github("andreasmock/MetaboDiff") # Delete ellipsis in OneDrive R folder if necessary and Click YES.

library(MetaboDiff)

#### Load data
assay <- metabolome_abundance.11
rowData <- `iPOP_Metablolite_Annotation.(excl..metabolites.without.HMDB.identifier).B`
colData <- Subject.data.T2DM.iHMP.4..2

#### Merge all objects into MultiAssayExperiment object2
met <- create_mae(assay, rowData, colData)

remove.packages("ellipsis")