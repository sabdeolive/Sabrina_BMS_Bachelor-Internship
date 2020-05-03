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
View(subject_info)
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

# does not seem like the row names of the T2D.rm sample data are the file names anymore, therefore, make them the file names
#rownames(T2D.rm@sam_data) <- T2D.rm@sam_data[,"Sample.ID.from.phyloseq"]

#### Only include subjects that are in the metabolome data, that have gut 16s data and that IR/IS classification
### Make vector of subjects to include (should be 50)
subjects.to.incl_df <- `Subject.data.T2DM.iHMP.(excl..subjects.wo.class,.metabolome,.gut.16s.and.not.in.phyloseq.to.use.as.vector)`
View(subjects.to.incl_df)
subjectIDs.to.incl_vec <- as.vector(subjects.to.incl_df[,2])
subjectIDs.to.incl_vec
class(subjectIDs.to.incl_vec)
length(subjectIDs.to.incl_vec) #1778 samples

### Remove all the subjects that are not in this vector from T2D.rm using prune_samples function (only keeps those that have been defined by e.g. vector)
T2D.rm <- prune_samples(subjectIDs.to.incl_vec, T2D.rm)
T2D.rm #no. of samples matches length of vector = correct
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 12062 taxa and 1778 samples ]
# sample_data() Sample Data:       [ 1778 samples by 22 sample variables ]
# tax_table()   Taxonomy Table:    [ 12062 taxa by 7 taxonomic ranks ]
View(T2D.rm)

#### Remove all samples from nasal cavity
T2D.rm.gut <- subset_samples(T2D.rm, sample_body_site == 'feces')
T2D.rm.gut
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 12062 taxa and 869 samples ]
# sample_data() Sample Data:       [ 869 samples by 22 sample variables ]
# tax_table()   Taxonomy Table:    [ 12062 taxa by 7 taxonomic ranks ]
View(T2D.rm.gut)
View(T2D.rm.gut@sam_data) # 50 samples = correct. 
# NOTE: Final phyloseq to work with (prior to taxa filtering) = T2D.rm.gut

write.table(T2D.rm.gut@otu_table@.Data, file="c:/Users/sabde/Documents/T2D ps w filtered subjects otu_table.txt", sep="\t", row.names = TRUE, col.names = NA)
# This file was used to compare the sample IDs to those in the metabolome data file Sample IDs in the metabolome data that did not match those in the phyloseq object were removed from the metabolome data file.

#### Removal of all samples from the T2D.rm.gut that are not shared between the T2D.rm.gut and the metabolome data
# NOTE: T2D.rm.gut has 869 samples and the metabolome data has 555 (after removing samples from the metabolome data that did not match the sample IDs in the phyloseq)
### Load metabolome data
metabolomics <- `metabolome_abundance.(excl..all.samples.not.in.phyloseq).7.2`

### Make a vector of the sample IDs in the metabolomics data file 
sampleIDs.met.vec <- as.vector(metabolomics[,1])
sampleIDs.met.vec
length(sampleIDs.met.vec) # 555 samples = correct

### Remove all the samples that are not in this vector from T2D.rm.gut using prune_samples function (only keeps those that have been defined by e.g. vector) 
T2D.rm.gut
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 12062 taxa and 869 samples ]
# sample_data() Sample Data:       [ 869 samples by 22 sample variables ]
# tax_table()   Taxonomy Table:    [ 12062 taxa by 7 taxonomic ranks ]
T2D.F.sub_sam <- prune_samples(sampleIDs.met.vec, T2D.rm.gut)
T2D.F.sub_sam
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 12062 taxa and 555 samples ]
# sample_data() Sample Data:       [ 555 samples by 22 sample variables ]
# tax_table()   Taxonomy Table:    [ 12062 taxa by 7 taxonomic ranks ]


#### Filter NA values
table(tax_table(T2D.F.sub_sam)[,"Family"],exclude = NULL) #1650 NA
T2D.fil <- subset_taxa(T2D.F.sub_sam, !is.na(Family)) 
T2D.fil #12062 taxa - 1650 = 10412 taxa

table(tax_table(T2D.F.sub_sam)[,"Genus"], exclude = NULL) #5953 NA (went up from 4303 to 5953 after merging and filtering of samples+subjects?)
T2D.fil <- subset_taxa(T2D.F.sub_sam, !is.na(Genus)) 
T2D.fil #10412 taxa -> 6109 taxa

##Filtering by species = not necessary 
# table(tax_table(T2D.F.sub_sam) [,"Species"], exclude = NULL) #4824 NA
# T2D.fil <- subset_taxa(T2D.F.sub_sam, !is.na(Species))
# T2D.fil #6109 taxa -> 1285 taxa

#### Subset the phyloseq class object into classification
IR_ps.fil <- subset_samples(T2D.fil, IR_IS_classification == "IR")
IR_ps.fil # T2D.fil = 555 samples, IR_ps.fil = 226 samples.
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 6109 taxa and 226 samples ]
# sample_data() Sample Data:       [ 226 samples by 22 sample variables ]
# tax_table()   Taxonomy Table:    [ 6109 taxa by 7 taxonomic ranks ]

IS_ps.fil <- subset_samples(T2D.fil, IR_IS_classification == "IS")
IS_ps.fil # T2D = 555 samples, IS_ps.fil = 329 samples (CHECK).
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 6109 taxa and 329 samples ]
# sample_data() Sample Data:       [ 329 samples by 22 sample variables ]
# tax_table()   Taxonomy Table:    [ 6109 taxa by 7 taxonomic ranks ]

# 226 + 329 = 555 samples (correct as T2D.fil only included classified subjects)

#################################################################################

####Create a contingency table of the number of taxa in each phylum
table(tax_table(T2D.fil)[, "Phylum"]) #1 phylum (Thermotogae) showed count of only 1

### Compute prevalence of every feature(/taxa?)
prevalence.df = apply(X = otu_table(T2D.fil),
                      MARGIN = ifelse(taxa_are_rows(T2D.fil), yes = 1, no = 2),
                      FUN = function(x){sum(x > 0)})

### add taxonomy and total read counts
prevalence.df = data.frame(Prevalence = prevalence.df,
                           TotalAbundance = taxa_sums(T2D.fil),
                           tax_table(T2D.fil))
plyr::ddply(prevalence.df, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))}) #column 1 = mean prevalence, column 2 = prevalence sum

## [Thermi],  Acidobacteria, Armatimonadetes,Chloroflexi, Cyanobacteria, Gemmatimonadetes, Nitrospirae, Planctomycetes, Tenericutes and Thermotogae are taxa that are only present in 0 or 1 sample.

#################################################################################

#### Prevalence filter and plotting prevalence for each classification
### Subset taxa
names.OTU <- taxa_names(T2D.fil)

## IR
keep.IR.taxa <- names.OTU[rowSums(IR_ps.fil@otu_table)>0] #makes a character vector of all the taxa to keep (i.e. all those present in at least 1 sample) in the phyloseq for IR
IR_ps.fil #6109 taxa 
IR_ps.fil <- prune_taxa(keep.IR.taxa, IR_ps.fil)
IR_ps.fil #2313 taxa

## IS
keep.IS.taxa <- names.OTU[rowSums(IS_ps.fil@otu_table)>0]
IS_ps.fil #2522 taxa
IS_ps.fil <- prune_taxa(keep.IS.taxa, IS_ps.fil)
IS_ps.fil #2522 taxa (CHECK!!!)

### Prevalence filter IR (filtering of taxa)
##Subset the remaining phyla 
prevalence.df.IR <- subset(prevalence.df, Phylum %in% get_taxa_unique(IR_ps.fil, "Phylum")) #get_taxa_uniqe is used to determine the different taxa present for a particular taxonomic rank in a given phyloseq-class object
## Visualize prevalence and total read count in order to determine prevalence threshold
ggplot(prevalence.df.IR, aes(TotalAbundance, Prevalence / nsamples(IR_ps.fil),color=Phylum)) + 
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) + 
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") + facet_wrap(~Phylum) + 
  theme(legend.position="none")
## Define prevalence threshold as 0.25% of total samples (CHECK: chose 0.05 as it seems there is a slight natural separation at around this prevalence in the actinobacteria and proteobacteria taxa)
prevalenceThreshold.IR <- 0.0025*nsamples(IR_ps.fil) #i.e. taxa have to appear in a minimum of 0.5% of samples or they will be removed
prevalenceThreshold.IR # 0.565 (i.e. the taxa would have to be prevalent in 0.565 samples in order to be considered)

### Execute this prevalence filter using prune_taxa() function (i.e. keeps rows/taxa where prevalence >= 0.05 in IR group)
keepTaxa.IR <- rownames(prevalence.df.IR)[(prevalence.df.IR$Prevalence >= prevalenceThreshold.IR)]
IR_ps.fil #2313 taxa
IR_ps.fil <- prune_taxa(keepTaxa.IR,IR_ps.fil) 
IR_ps.fil #2313 taxa -> 2313 taxa (CHECK!!!)


### Prevalence filter IS
##Subset the remaining phyla 
prevalence.df.IS <- subset(prevalence.df, Phylum %in% get_taxa_unique(IS_ps.fil, "Phylum")) #get_taxa_uniqe is used to determine the different taxa present for a particular taxonomic rank in a given phyloseq-class object
## Visualize prevalence and total read count in order to determine prevalence threshold
ggplot(prevalence.df.IS, aes(TotalAbundance, Prevalence / nsamples(IS_ps.fil),color=Phylum)) + 
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) + 
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") + facet_wrap(~Phylum) + 
  theme(legend.position="none")
## Define prevalence threshold as 0.25% of total samples (CHECK: chose 0.05 as it seems there is a slight natural separation at around this prevalence in the actinobacteria and proteobacteria taxa)
prevalenceThreshold.IS <- 0.0025*nsamples(IS_ps.fil) #i.e. taxa have to appear in a minimum of 0.5% of samples or they will be removed
prevalenceThreshold.IS # 0.8225 (i.e. the taxa would have to be prevalent in 0.8225 samples in order to be considered)

### Execute this prevalence filter using prune_taxa() function (i.e. keeps rows/taxa where prevalence >= 0.05 in IR group)
keepTaxa.IS <- rownames(prevalence.df.IS)[(prevalence.df.IS$Prevalence >= prevalenceThreshold.IS)]
IS_ps.fil #2522 taxa
IS_ps.fil <- prune_taxa(keepTaxa.IS,IS_ps.fil) 
IS_ps.fil #2522 taxa -> 2522 taxa (CHECK!!!)


### Filter taxa of the whole T2D.rm.gut phyloseq-class object using IR and IS prevalence filtration
keepTaxa.T2D.fil <- c(keepTaxa.IR, keepTaxa.IS) 
T2D.fil #6109 taxa 
T2D.fil <- prune_taxa(keepTaxa.T2D.fil, T2D.fil)
T2D.fil #2888 taxa

########################################################################################################

#### PCoA comparing IR and IS
pslog <- transform_sample_counts(T2D.fil, function(x) log(1 + x))
out.pcoa.log <- ordinate(pslog, method = "PCoA", distance = "bray") # Anna used "PCoA" and "bray"?
evals <- out.pcoa.log$values$Eigenvalues
plot_ordination(pslog, out.pcoa.log, type = "samples", color = "IR_IS_classification") +
  labs(col = "classification") +
  coord_fixed(sqrt(evals[2] / evals[1])) # what does type = "samples" do? Anna added it to the general workflow.

# Including race as another variable
# plot_ordination(pslog, out.pcoa.log, color = "IR_IS_classification", shape = "Race") +
# labs(col = "classification") +
#   coord_fixed(sqrt(evals[2] / evals[1]))

###########################################################################################################################

#### Alpha diversity
BiocManager::install("microbiome")
library(microbiome)

alpha.div <- head(alpha(T2D.fil, index = "all"))
View(alpha.div)

###########################################################################################################################

#### Venn diagrams (DOESN'T WORK)
install.packages("VennDiagram")
library(VennDiagram)

names.OTU <- taxa_names(T2D.fil)
names.OTU.IR <- taxa_names(IR_ps.fil)
names.OTU.IS <- taxa_names(IS_ps.fil)
library("grid")
grid.newpage()
venn.diagram(names.OTU.IR,names.OTU.IS, "IR", "IS", colors= c("#e87396","#2a96a0","#9bc6ff", euler=T)) # to attempt to scale the relative volumes of the sets

###########################################################################################################################



###########################################################################################################################

#### Multitable analysis 
### Quick check
dim(metabolomics)
# 555 330
T2D.fil
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2888 taxa and 555 samples ]
# sample_data() Sample Data:       [ 555 samples by 22 sample variables ]
# tax_table()   Taxonomy Table:    [ 2888 taxa by 7 taxonomic ranks ]

# same number of samples

### Removing metabolites that are zero across many samples
keep_ix <- rowSums(metabolomics == 0) <=3
metabolomics.fil <- metabolomics[keep_ix]

### Removing microbes that are zero across many samples

