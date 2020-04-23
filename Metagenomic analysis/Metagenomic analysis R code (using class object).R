#Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

# Check if BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# Install HMP2Data package using BiocManager
BiocManager::install("HMP2Data")

library(BiocManager)
library(HMP2Data)
library(phyloseq)
library(biomformat)
library(ggplot2)
library(readxl)
library(plyr)
library(gridExtra)
library(dplyr)

T2D <- T2D16S()
T2D

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 12062 taxa and 2208 samples ]
#sample_data() Sample Data:       [ 2208 samples by 13 sample variables ]
#tax_table()   Taxonomy Table:    [ 12062 taxa by 7 taxonomic ranks ]

#Create a contingency table of the number of taxa in each phylum
table(tax_table(T2D)[, "Phylum"]) #4 phylum showed count of only 1

# Compute prevalence of every feature
prevalence.df = apply(X = otu_table(T2D),
                      MARGIN = ifelse(taxa_are_rows(T2D), yes = 1, no = 2),
                      FUN = function(x){sum(x > 0)})

# add taxonomy and total read counts
prevalence.df = data.frame(Prevalence = prevalence.df,
                           TotalAbundance = taxa_sums(T2D),
                           tax_table(T2D))
plyr::ddply(prevalence.df, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})


# Elusimicrobia, OP11, SBR1093, Thermotogae,TM6, WS2 are taxa that are only present in 1 or 2 samples (CHECK)
# Filter the Phylum that only appear in 1 or 2 samples with low read count
filterPhyla = c("Elusimicrobia", "OP11", "SBR1093", "Thermotogae","TM6", "WS2")
T2D.f= subset_taxa(T2D, !Phylum %in% filterPhyla) # 12054 taxa
T2D.f


### prevalence filter by classification (NOTE: The taxonomic filter is included in this)
# subset samples
#IR.subset.f <- subset_samples(T2D.f,classification =="IR") 
#IS.subset.f <- subset_samples(T2D.f,classification =="IS") 



