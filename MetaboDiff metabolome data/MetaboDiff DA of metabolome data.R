# NOTE: have to have downloaded RTools for R version 3.6.3 prior to this 

#### Installation and loading of MetaboDiff
install.packages("WGCNA")

install.packages("devtools")
library("devtools")
devtools::install_github("andreasmock/MetaboDiff") # Delete ellipsis in OneDrive R folder if necessary and Click YES.

library(MetaboDiff)

#### Load data (NOTE: make first column row names)
assay <- metabolome_abundance.11
rowData <- `iPOP_Metablolite_Annotation.(w.1.HMDB.per.metabolite).2B`
colData <- Subject.data.T2DM.iHMP.4..2
colData <- colData[,-8]

#### Merge all objects into MultiAssayExperiment object2
met <- create_mae(assay, rowData, colData)
# A MultiAssayExperiment object of 1 listed
# experiment with a user-defined name and respective class.
# Containing an ExperimentList class object of length 1:
#   [1] raw: SummarizedExperiment with 323 rows and 60 columns
# Features:
#   experiments() - obtain the ExperimentList instance
# colData() - the primary/phenotype DFrame
# sampleMap() - the sample availability DFrame
# `$`, `[`, `[[` - extract colData columns, subset, or experiment
# *Format() - convert into a long or wide DFrame
# assays() - convert ExperimentList to a SimpleList of matrices

############################################################################
##### Data preprocessing 

#### Defining cut-off for raw metabolite measurements
met <- knn_impute(met,cutoff=0.4) # i.e. only keeping those with raw measurements in more than 60% of cases
# A MultiAssayExperiment object of 2 listed
# experiments with user-defined names and respective classes.
# Containing an ExperimentList class object of length 2:
#   [1] raw: SummarizedExperiment with 323 rows and 60 columns
# [2] imputed: SummarizedExperiment with 323 rows and 60 columns
# Features:
#   experiments() - obtain the ExperimentList instance
# colData() - the primary/phenotype DFrame
# sampleMap() - the sample availability DFrame
# `$`, `[`, `[[` - extract colData columns, subset, or experiment
# *Format() - convert into a long or wide DFrame
# assays() - convert ExperimentList to a SimpleList of matrices

#### Metabolite annotation
met <- get_SMPDBanno(met,
                     column_kegg_id=2,
                     column_hmdb_id=3,
                     column_chebi_id=NA)

#### Imputation of missing values
na_heatmap(met,
           group_factor="IR_IS_classification",
           label_colors=c("darkseagreen","dodgerblue")) 
# Error: You should have at least two distinct break values. -> seems that this error occurs if do not have NAs in input data (https://github.com/andreasmock/MetaboDiff/issues/6)

#### Outlier heatmap and removal of outliers
outlier_heatmap(met,
                group_factor="IR_IS_classification",
                label_colors=c("darkseagreen","dodgerblue"),
                k=2)
# remove outliers (CHECK if necessary)
met.out <- remove_cluster(met,cluster=2) # 60 -> 39 columns (21 samples removed)
# A MultiAssayExperiment object of 2 listed
# experiments with user-defined names and respective classes.
# Containing an ExperimentList class object of length 2:
#   [1] raw: SummarizedExperiment with 323 rows and 39 columns
# [2] imputed: SummarizedExperiment with 323 rows and 39 columns
# Features:
#   experiments() - obtain the ExperimentList instance
# colData() - the primary/phenotype DFrame
# sampleMap() - the sample availability DFrame
# `$`, `[`, `[[` - extract colData columns, subset, or experiment
# *Format() - convert into a long or wide DFrame
# assays() - convert ExperimentList to a SimpleList of matrices

#### Normalization (variance stabilizing normalization)
met.nor <- normalize_met(met.out)
# A MultiAssayExperiment object of 4 listed
# experiments with user-defined names and respective classes.
# Containing an ExperimentList class object of length 4:
#   [1] raw: SummarizedExperiment with 323 rows and 39 columns
# [2] imputed: SummarizedExperiment with 323 rows and 39 columns
# [3] norm: SummarizedExperiment with 323 rows and 39 columns
# [4] norm_imputed: SummarizedExperiment with 323 rows and 39 columns
# Features:
#   experiments() - obtain the ExperimentList instance
# colData() - the primary/phenotype DFrame
# sampleMap() - the sample availability DFrame
# `$`, `[`, `[[` - extract colData columns, subset, or experiment
# *Format() - convert into a long or wide DFrame
# assays() - convert ExperimentList to a SimpleList of matrices

#### Quality-control normalization (show distribution of raw and normalized metabolic measurements for every sample in the study set)
quality_plot(met.nor,
             group_factor="IR_IS_classification",
             label_colors=c("darkseagreen","dodgerblue")) 

############################################################################
##### Data analysis 
#### Unsupervised analysis (exploration of metabolomic profiles between samples)
source("http://peterhaschke.com/Code/multiplot.R")
multiplot(
  pca_plot(met.nor,
           group_factor="IR_IS_classification",
           label_colors=c("darkseagreen","dodgerblue")),
  tsne_plot(met.nor,
            group_factor="IR_IS_classification",
            label_colors=c("darkseagreen","dodgerblue")),
  cols=2)
# no metabolome-wide difference between the classifications.
# the tSNE plot does not reveal a distinct difference in metabolomic profiles between the different classification samples.

#### Hypothesis testing w/ correction for multiple testing (applies a Student's T-Test since there are only 2 groups and p-values are corrected using Benjamini-Hochberg procedure)
#(do the individual metabolites show differential abundance between the 2 groups?)
met.test <- diff_test(met.nor,
                group_factors = c("IR_IS_classification"))

str(metadata(met.test), max.level = 2)
# List of 1
# $ ttest_IR_IS_classification_IS_vs_IR:'data.frame':	323 obs. of  3 variables:
#   ..$ pval       : num [1:323] 0.232 0.673 0.137 0.451 0.26 ...
# ..$ adj_pval   : num [1:323] 0.758 0.939 0.655 0.91 0.805 ...
# ..$ fold_change: num [1:323] 0.1019 -0.0338 -0.2102 0.0715 0.4798 ...

DA.results <- metadata(met.test)
View(DA.results) # which rows refer to which metabolites?

### Visualization of this comparative analysis using volcano plot
par(mfrow=c(1,2))
volcano_plot(met.test, 
             group_factor="IR_IS_classification",
             label_colors=c("darkseagreen","dodgerblue"),
             p_adjust = FALSE)
volcano_plot(met.test, 
             group_factor="IR_IS_classification",
             label_colors=c("darkseagreen","dodgerblue"),
             p_adjust = TRUE)


#### Metabolic correlation network analysis 
met.net <- met.test %>%
  diss_matrix %>%
  identify_modules(min_module_size=5) %>%
  name_modules(pathway_annotation="SUB_PATHWAY") %>%
  calculate_MS(group_factors=c("IR_IS_classification"))
# ..cutHeight not given, setting it to 0.985  ===>  99% of the (truncated) height range in dendro.
# ..done.

# no alpha given 

### Association of sample traits with correlation modules (MS plot)
# (which metabolic subpathways (i.e. modeules) differ between groups?) 
MS_plot(met.net,
        group_factor="IR_IS_classification",
        p_value_cutoff=0.05,
        p_adjust=FALSE)
# no pathways mentioned and no significant module significance  

### Exploration of individual metabolites with correlation module (MOI plot)
# (which metabolite is most closely related to the individual subpathway (i.e. modules)?)
MOI_plot(met.net,
         group_factor="IR_IS_classification",
         MOI = 2,
         label_colors=c("darkseagreen","dodgerblue"),
         p_adjust = FALSE) + xlim(c(-1,8))
# DOES NOT WORK: Error in data.frame(mets = mets, x = x, y = y, fc = fc) : arguments imply differing number of rows: 0, 16
