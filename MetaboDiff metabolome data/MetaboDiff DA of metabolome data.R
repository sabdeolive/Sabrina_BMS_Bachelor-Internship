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

met.out <- remove_cluster(met,cluster=1)
# A MultiAssayExperiment object of 2 listed
# experiments with user-defined names and respective classes.
# Containing an ExperimentList class object of length 2:
#   [1] raw: SummarizedExperiment with 323 rows and 21 columns
# [2] imputed: SummarizedExperiment with 323 rows and 21 columns
# Features:
#   experiments() - obtain the ExperimentList instance
# colData() - the primary/phenotype DFrame
# sampleMap() - the sample availability DFrame
# `$`, `[`, `[[` - extract colData columns, subset, or experiment
# *Format() - convert into a long or wide DFrame
# assays() - convert ExperimentList to a SimpleList of matrices

met.out2 <- remove_cluster(met,cluster=2)
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
met.nor <- normalize_met(met)
# A MultiAssayExperiment object of 4 listed
# experiments with user-defined names and respective classes.
# Containing an ExperimentList class object of length 4:
#   [1] raw: SummarizedExperiment with 323 rows and 60 columns
# [2] imputed: SummarizedExperiment with 323 rows and 60 columns
# [3] norm: SummarizedExperiment with 323 rows and 60 columns
# [4] norm_imputed: SummarizedExperiment with 323 rows and 60 columns
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
                group_factors = "IR_IS_classification")

str(metadata(met.test), max.level = 2)
# List of 1
# $ ttest_IR_IS_classification_IS_vs_IR:'data.frame':	323 obs. of  3 variables:
#   ..$ pval       : num [1:323] 0.0302 0.08 0.8556 0.3998 0.1194 ...
# ..$ adj_pval   : num [1:323] 0.326 0.445 0.927 0.758 0.521 ...
# ..$ fold_change: num [1:323] 0.145 -0.1154 -0.0211 0.0681 0.6101 ...

DA.results <- metadata(met.test)
View(DA.results[["ttest_IR_IS_classification_IS_vs_IR"]]) # which rows refer to which metabolites?

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

dev.off()

#### Metabolic correlation network analysis 
met.net <- met.test %>%
  diss_matrix %>%
  identify_modules(min_module_size=5) %>%
  name_modules(pathway_annotation="Pathway") %>%
  calculate_MS(group_factors="IR_IS_classification")
# ..cutHeight not given, setting it to 0.99  ===>  99% of the (truncated) height range in dendro.
# ..done.

# no alpha given 

met.net1 <- met.test %>%
  diss_matrix %>%
  identify_modules(min_module_size=5) %>%
  name_modules(pathway_annotation="Chemical.Class") %>%
  calculate_MS(group_factors="IR_IS_classification")
# ..cutHeight not given, setting it to 0.99  ===>  99% of the (truncated) height range in dendro.
# ..done.

# no alpha given 

### Association of sample traits with correlation modules (MS plot)
# (which metabolic subpathways (i.e. modules) differ between groups?) 
MS_plot(met.net,
        group_factor="IR_IS_classification",
        p_value_cutoff=0.05,
        p_adjust=FALSE)
# looks weird and no significant module significance  
dev.off()

# (which chemical classes differ between groups?) 
MS_plot(met.net1,
        group_factor="IR_IS_classification",
        p_value_cutoff=0.05,
        p_adjust=FALSE)
# no significant module significance  
dev.off()

### Exploration of individual metabolites with correlation module (MOI plot)
# (which metabolite is most closely related to the individual subpathway (i.e. modules)?)
MOI_plot(met.net,
         group_factor="IR_IS_classification",
         MOI = 2,
         label_colors=c("darkseagreen","dodgerblue"),
         p_adjust = FALSE) + xlim(c(-1,8))
# DOES NOT WORK: Error in data.frame(mets = mets, x = x, y = y, fc = fc) : arguments imply differing number of rows: 0, 16

############################################################################

# Trying to figure out if there is a sample variable that the 2 clusters refer to 
metout.BMI <- as.vector(met.out@colData@listData[["BMI"]])
metout2.BMI <- as.vector(met.out2@colData@listData[["BMI"]])
mean(metout.BMI) # 28.46095
mean(metout2.BMI, na.rm = TRUE) # 29.10211

metout.SSPG <- as.vector(met.out@colData@listData[["SSPG"]])
metout2.SSPG <- as.vector(met.out2@colData@listData[["SSPG"]])
mean(metout.SSPG) # 166.0624
mean(metout2.SSPG) # 141.0321
# seems like a possibility? since study has divided patients into IR if SSPG > 150 mg/dl and IR if SSPG < 150 mg/dl (https://www.nature.com/articles/s41586-019-1236-x)
# according to 'Subject data T2DM iHMP (testing SSPG vs. class) 4.3' datafile this SSPG criteria was used for classification, hence, there may just be some outliers?

metout.Age <- as.vector(met.out@colData@listData[["Age"]])
metout2.Age <- as.vector(met.out2@colData@listData[["Age"]])
mean(metout.Age) # 55.19048
mean(metout2.Age) # 56.55615