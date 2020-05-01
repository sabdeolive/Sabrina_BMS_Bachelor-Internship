library(BiocManager)
library(HMP2Data)
library(phyloseq)
library(biomformat)
library(ggplot2)
library(readxl)
library(plyr)
library(gridExtra)
library(dplyr)

sessionInfo()

T2D <- T2D16S()
subject_info<- `Subject.data.T2DM.iHMP.(with.sample.IDs.for.subjects.in.phyloseq.sam_data)`
subject_info
T2D

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 12062 taxa and 2208 samples ]
#sample_data() Sample Data:       [ 2208 samples by 13 sample variables ]
#tax_table()   Taxonomy Table:    [ 12062 taxa by 7 taxonomic ranks ]

#### Make rownames of the subject data the same as those of the phyloseq sample data 
rownames(subject_info) <- subject_info[,2] 

#### Removing sample/subject not in the subject data from the phyloseq sample data
T2D.rm <- subset_samples(T2D, !(file_name %in% "HMP2_J34982_1_NS_T0_B0_0120_ZIWTAHN-01_ANAV8"))
T2D.rm
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 12062 taxa and 2207 samples ]
# sample_data() Sample Data:       [ 2207 samples by 13 sample variables ]
# tax_table()   Taxonomy Table:    [ 12062 taxa by 7 taxonomic ranks ]
dim(subject_info) # 2207 rows = samples
# the no. of samples in the phyloseq now matches the number of samples in the subject data YAY

#### Merge phyloseq with removed sample and subject data
subject_info.sd <- sample_data(subject_info)
subject_info.sd
class(subject_info.sd)
T2D.sd <- merge_phyloseq(T2D.rm@sam_data, subject_info.sd)
T2D.sd
View(T2D.sd)
T2D.rm@sam_data <- T2D.sd
T2D.rm
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 12062 taxa and 2207 samples ]
# sample_data() Sample Data:       [ 2207 samples by 22 sample variables ]
# tax_table()   Taxonomy Table:    [ 12062 taxa by 7 taxonomic ranks ]
View(T2D.rm)
View(T2D.rm@sam_data)

#T2D.clas <- merge_phyloseq(T2D.sd,T2D.rm@otu_table, T2D.rm@tax_table) 
#T2D.clas 
#View(T2D.clas@sam_data)#didn't work, have NA values for all IR/IS classifications. 