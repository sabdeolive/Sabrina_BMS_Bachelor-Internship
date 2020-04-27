#Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

#Install phyloseq
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

sessionInfo()

T2D <- T2D16S()
subject_info<-Subject.data.T2DM.iHMP.1
subject_info
T2D

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 12062 taxa and 2208 samples ]
#sample_data() Sample Data:       [ 2208 samples by 13 sample variables ]
#tax_table()   Taxonomy Table:    [ 12062 taxa by 7 taxonomic ranks ]

#### Changing file name variable in T2D16S_samp to subject_ID
as_subject_ID <- substr(T2D@sam_data[["file_name"]],29,35)
as_subject_ID
class(as_subject_ID)


#### Filter NA values
table(tax_table(T2D)[,"Family"],exclude = NULL) #1650 NA
T2D <- subset_taxa(T2D, !is.na(Family)) 
T2D #12062 taxa - 1650 = 10412 taxa

table(tax_table(T2D)[,"Genus"], exclude = NULL) #4303 NA
T2D <- subset_taxa(T2D, !is.na(Genus)) 
T2D #10412 taxa -> 6109 taxa


table(tax_table(T2D) [,"Species"], exclude = NULL) #4824 NA
T2D <- subset_taxa(T2D, !is.na(Species))
T2D #6109 taxa -> 1285 taxa

#### Filter phylum that do not appear in a lot of samples 
##Create a contingency table of the number of taxa in each phylum
table(tax_table(T2D)[, "Phylum"]) #3 phylum (Acidobacteria, Fusobacteria and Synergistetes) showed count of only 1

## Compute prevalence of every feature
prevalence.df = apply(X = otu_table(T2D),
                      MARGIN = ifelse(taxa_are_rows(T2D), yes = 1, no = 2),
                      FUN = function(x){sum(x > 0)})

## add taxonomy and total read counts
prevalence.df = data.frame(Prevalence = prevalence.df,
                           TotalAbundance = taxa_sums(T2D),
                           tax_table(T2D))
plyr::ddply(prevalence.df, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})


## Acidobacteria, cyanobacteria and fusobacteria are taxa that are only present in 1 or 2 samples (CHECK)
# Filter the Phylum that only appear in 1 or 2 samples with low read count (might not be necessary as use prevalence threshold later on CHECK)
filterPhyla = c("Acidobacteria", "Cyanobacteria", "Fusobacteria")
T2D.fil = subset_taxa(T2D, !Phylum %in% filterPhyla) # 1285 taxa
T2D.fil #1281 taxa


#### Subsetting subject info for IR and IS
IR<-subset.data.frame(subject_info, IR_IS_classification=="IR")
IS<-subset.data.frame(subject_info, IR_IS_classification=="IS")

IR #35 individuals
IS #31 individuals


#### Remove subjects that are not in the gut 16s and metabolome data from the IR and IS subsets

###IR
##IR subjects not in gut 16s = ZPJT5MY, ZQEFRDE
##IR subjects not in metabolome = ZSOZWGV

##Finding rows of subjects not in gut 16s
which(grepl("ZPJT5MY", IR$SubjectID)) #row 13
which(grepl("ZQEFRDE", IR$SubjectID)) #row 14

##Finding rows of subjects not in metabolome
which(grepl("ZSOZWGV", IR$SubjectID)) #row 18

#Removing row 13,14 and 18 from IR data
IR_excl <- IR[-c(13,14,18),]
View(IR_excl)
dim(IR_excl) #35 - 3 = 32 subjects
which(grepl("ZPJT5MY", IR_excl$SubjectID))
which(grepl("ZQEFRDE", IR_excl$SubjectID))
which(grepl("ZSOZWGV", IR_excl$SubjectID))

#IR_to_exclude <- row.names(IR)
#IR_all_data <- IR[-c(ZPJT5MY, ZQEFRDE, ZSOZWGV),]
#IR_all_data <- IR[-c(1,2,3),]
#class(IR)

###IS
##IS subjects not in gut 16s = ZP0DXQ0, ZRGK7U8, ZS2DMX7
##IS subjects not in metabolome = none

##Finding rows of subjects not in gut 16s
which(grepl("ZP0DXQ0", IS$SubjectID)) #row 15
which(grepl("ZRGK7U8", IS$SubjectID)) #row 20
which(grepl("ZS2DMX7", IS$SubjectID)) #row 22

#Removing row 13,14 and 18 from IS data
IS_excl <- IS[-c(15,20,22),] #31 - 3 = 28 subjects
dim(IS_excl)
View(IS_excl)
which(grepl("ZP0DXQ0", IS_excl$SubjectID))
which(grepl("ZRGK7U8", IS_excl$SubjectID))
which(grepl("ZS2DMX7", IS_excl$SubjectID))


#### Make subject IDs from each classification into a vector
class(IR_excl)#data frame
class(IS_excl)#data frame

IR_v <- as.vector(IR_excl[,1])
class(IR_v)#character
IR_v

IS_v <- as.vector(IS_excl[,1])
class(IS_v)#character
IS_v


#### Subset the phyloseq class object into classification using IR_v and IS_v
IR_ps.fil <- subset_samples(T2D.fil, as_subject_ID == IR_v)
IR_ps.fil

IS_ps.fil <- subset_samples(T2D.fil, as_subject_ID == IS_v)
IS_ps.fil


#######################################################################################
  

##add as_subject_ID to T2D16S_samp data
#first change as_subject_ID to sample data
#SD_subject_ID <- sample_data(as_subject_ID, errorIfNULL = T)


#T2D16S_samp <- merge_phyloseq(T2D@sam_data,as_subject_ID)
#T2D16S_samp
#order the subject



#match_subject_ID <- sample_data(as_subject_ID, 
#ample_data(T2D16S_samp$match_sample_ID <-sample_sums()


####Remove subjects that are not in the gut 16s and metabolome data from the IR and IS subsets


###T2D
##Remove subjects that are not classified(= unknown)


#, not in gut 16s and not in metabolome data from the T2D data.
#subjects not in gut 16s = ZP0DXQ0, ZPJT5MY, ZQEFRDE, ZRGK7U8, ZS2DMX7
#subjects not in metabolome = ZSOZWGV
#T2D_all_data
  
#########################################################################################
###Filter NA values
# table(tax_table(T2D)[,"Family"],exclude = NULL) #1650 NA
# T2D <- subset_taxa(T2D, !is.na(Family)) 
# T2D #12062 taxa - 1650 = 10412 taxa
# 
# table(tax_table(T2D)[,"Genus"], exclude = NULL) #5953 NA
# T2D <- subset_taxa(T2D, !is.na(Genus)) 
# T2D #10412 taxa -> 6109 taxa
# 
# 
# table(tax_table(T2D) [,"Species"], exclude = NULL) #4824 NA
# T2D <- subset_taxa(T2D, !is.na(Species))
# T2D #6109 taxa -> 1285 taxa
# 
# ###Create a contingency table of the number of taxa in each phylum
# table(tax_table(T2D)[, "Phylum"]) #4 phylum showed count of only 1
# 
# # Compute prevalence of every feature
# prevalence.df = apply(X = otu_table(T2D),
#                       MARGIN = ifelse(taxa_are_rows(T2D), yes = 1, no = 2),
#                       FUN = function(x){sum(x > 0)})
# 
# # add taxonomy and total read counts
# prevalence.df = data.frame(Prevalence = prevalence.df,
#                            TotalAbundance = taxa_sums(T2D),
#                            tax_table(T2D))
# plyr::ddply(prevalence.df, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# 
# 
# # Elusimicrobia, OP11, SBR1093, Thermotogae,TM6, WS2 are taxa that are only present in 1 or 2 samples (CHECK)
# # Filter the Phylum that only appear in 1 or 2 samples with low read count
# filterPhyla = c("Elusimicrobia", "OP11", "SBR1093", "Thermotogae","TM6", "WS2")
# T2D.fil = subset_taxa(T2D, !Phylum %in% filterPhyla) # 12054 taxa
# T2D.fil 

#### Prevalence filter and plotting prevalence for each classification 
### Subset taxa (done before, see 'Subset the phyloseq class object into classification using IR_v and IS_v')
names.OTU <- taxa_names(T2D.fil)

## IR
keep.IR.taxa <- names.OTU[rowSums(IR_ps.fil@otu_table)>0]
IR_ps.fil <- prune_taxa(keep.IR.taxa, IR_ps.fil) #no of taxa?

## IS
keep.IS.taxa <- names.OTU[rowSums(IS_ps.fil@otu_table)>0]
IS_ps.fil <- prune_taxa(keep.IS.taxa, IS_ps.fil) #no of taxa?

### Prevalence filter IR (filtering of taxa)
##Subset the remaining phyla 
prevalence.df.IR <- subset(prevalence.df, Phylum %in% get_taxa_unique(IR_ps.fil, "Phylum"))
## Visualize prevalence and total read count in order to determine prevalence threshold
ggplot(prevalence.df.IR, aes(TotalAbundance, Prevalence / nsamples(IR_ps.fil),color=Phylum)) + 
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) + 
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") + facet_wrap(~Phylum) + 
  theme(legend.position="none")
## Define prevalence threshold as 5% of total samples (CHECK)
prevalenceThreshold.IR <- 0.05*nsamples(IR_ps.fil) #i.e. taxa have to appear in a minimum of 0.5% of samples or they will be removed
prevalenceThreshold.IR # 1.05

### Execute this prevalence filter using prune_taxa() function (i.e. keeps rows/taxa where prevalence >= 0.05 in IR group)
keepTaxa.IR <- rownames(prevalence.df.IR)[(prevalence.df.IR$Prevalence >= prevalenceThreshold.IR)]
IR_ps.fil #415 taxa
IR_ps.fil <- prune_taxa(keepTaxa.IR,IR_ps.fil) 
IR_ps.fil #415 taxa -> 413 taxa


### Prevalence filter IS
### ggplot for prevalence IS


### prevalence filter by classification (NOTE: The taxonomic filter is included in this)
# subset samples
# IR.subset.f <- subset_samples(T2D.f,classification =="IR")
# IS.subset.f <- subset_samples(T2D.f,classification =="IS")



