#Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

if (!requireNamespace("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("phyloseq")

# Check if BiocManager is installed
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# Install HMP2Data package using BiocManager
#BiocManager::install("HMP2Data")

library(BiocManager)
library(phyloseq)
library(biomformat)
library(readxl)
library(plyr)
library(ggplot2)
library(gridExtra)
library(dplyr)

#splitting data into different experimental groups
