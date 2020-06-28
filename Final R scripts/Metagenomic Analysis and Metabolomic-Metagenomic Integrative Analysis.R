######################################################################
# This script contains the code to perform the analysis of the metagenomic data and the integrative analysis of the metabolomic and metagenomic data
######################################################################


############################## Metagenomic analysis ##############################

library(BiocManager)
library(HMP2Data)
library(phyloseq)
library(biomformat)
library(ggplot2)
library(readxl)
library(plyr)
library(gridExtra)
library(dplyr)

###### Load data
T2D <- T2D16S()
subject_info<- `Subject.data.T2DM.iHMP.(with.sample.IDs.for.subjects.in.phyloseq.sam_data)`

##### Phyloseq preprocessing
#### Make rownames of the subject data the same as those of the phyloseq sample data 
rownames(subject_info) <- subject_info[,2] 

#### Remove samples/subjects not in the subject data from the phyloseq sample data
T2D.rm <- subset_samples(T2D, !(file_name %in% "HMP2_J34982_1_NS_T0_B0_0120_ZIWTAHN-01_ANAV8"))
# NOTE: the number of samples in the phyloseq should now match the number of samples in the subject data

#### Merge phyloseq with removed sample and subject data
subject_info.sd <- sample_data(subject_info)
T2D.sd <- merge_phyloseq(T2D.rm@sam_data, subject_info.sd)
T2D.rm@sam_data <- T2D.sd

#### Only include classified subjects that are in the metabolome data
### Make vector of subjects to include (should be 50)
subjects.to.incl_df <- `Subject.data.T2DM.iHMP.(excl..subjects.wo.class,.metabolome,.gut.16s.and.not.in.phyloseq.to.use.as.vector)`
subjectIDs.to.incl_vec <- as.vector(subjects.to.incl_df[,2])

### Remove all the subjects that are not in this vector from T2D.rm using prune_samples function (only keeps those that have been defined by e.g. vector)
T2D.rm <- prune_samples(subjectIDs.to.incl_vec, T2D.rm)
# number of samples should match the length of the previous vector

#### Remove all samples from the nasal cavity
T2D.rm.gut <- subset_samples(T2D.rm, sample_body_site == 'feces')

#### Removal of all samples from the T2D.rm.gut that are not shared between the T2D.rm.gut and the metabolome data
# NOTE: T2D.rm.gut has 869 samples and the metabolome data has 555 (after removing samples from the metabolome data that did not match the sample IDs in the phyloseq)
### Load metabolome data
metabolomics <- `metabolome_abundance.(excl..all.sample.IDs.that.have.duplicates).7.4`
### Make a vector of the sample IDs in the metabolomics data file 
sampleIDs.met.vec <- as.vector(metabolomics[,1])
### Remove all the samples that are not in this vector from T2D.rm.gut using prune_samples function (only keeps those that have been defined by e.g. vector) 
T2D.F.sub_sam <- prune_samples(sampleIDs.met.vec, T2D.rm.gut)
# NOTE: number of samples should match number of rows in the metabolomics data

#### Remove samples that have an Axis1 value of less than -2.8 (i.e. outliers from PCA)
### Load file with samples to include
excl.out <- `List.of.samples.to.include.from.PCA.(to.use.as.vector)`
### Make a vector of the sample IDs in the excl.out data file 
excl.out.vec <- as.vector(excl.out[,1])
### Remove all the samples that are not in this vector from T2D.F.sub_sam using prune_samples function (only keeps those that have been defined by e.g. vector) 
T2D.excl.out <- prune_samples(excl.out.vec, T2D.F.sub_sam)

#### Filter any taxa with NA values for family and/or genus
### Family
table(tax_table(T2D.excl.out)[,"Family"],exclude = NULL)
T2D.fil <- subset_taxa(T2D.excl.out, !is.na(Family)) 
### Genus
table(tax_table(T2D.excl.out)[,"Genus"], exclude = NULL)
T2D.fil <- subset_taxa(T2D.excl.out, !is.na(Genus)) 

#### Subset the phyloseq class object into classification
### IR
IR_ps.fil <- subset_samples(T2D.fil, IR_IS_classification == "IR")
### IS
IS_ps.fil <- subset_samples(T2D.fil, IR_IS_classification == "IS")

#### Apply a prevalence filter
### Compute prevalence of every feature(/taxa)
prevalence.df = apply(X = otu_table(T2D.fil),
                      MARGIN = ifelse(taxa_are_rows(T2D.fil), yes = 1, no = 2),
                      FUN = function(x){sum(x > 0)})
### Add taxonomy and total read counts to this prevalence data frame
prevalence.df = data.frame(Prevalence = prevalence.df,
                           TotalAbundance = taxa_sums(T2D.fil),
                           tax_table(T2D.fil))
### Subset taxa
names.OTU <- taxa_names(T2D.fil)
## IR
keep.IR.taxa <- names.OTU[rowSums(IR_ps.fil@otu_table)>0] # makes a character vector of all the taxa to keep (i.e. all those present in at least 1 sample) in the phyloseq for IR
IR_ps.fil <- prune_taxa(keep.IR.taxa, IR_ps.fil)
## IS
keep.IS.taxa <- names.OTU[rowSums(IS_ps.fil@otu_table)>0]
IS_ps.fil <- prune_taxa(keep.IS.taxa, IS_ps.fil)
### Creating a prevalence filter for IR (to filter taxa)
## Subset the remaining phyla 
prevalence.df.IR <- subset(prevalence.df, Phylum %in% get_taxa_unique(IR_ps.fil, "Phylum"))
## Visualize prevalence and total read count in order to determine prevalence threshold
ggplot(prevalence.df.IR, aes(TotalAbundance, Prevalence / nsamples(IR_ps.fil),color=Phylum)) + 
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) + 
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") + facet_wrap(~Phylum) + 
  theme(legend.position="none")
## Define prevalence threshold as 10% of total samples 
prevalenceThreshold.IR <- 0.10*nsamples(IR_ps.fil) #i.e. taxa have to appear in a minimum of 10% of samples or they will be removed
### Apply this prevalence filter using prune_taxa() 
keepTaxa.IR <- rownames(prevalence.df.IR)[(prevalence.df.IR$Prevalence >= prevalenceThreshold.IR)]
IR_ps.fil <- prune_taxa(keepTaxa.IR,IR_ps.fil) 
### Creating a prevalence filter for IS (to filter taxa)
## Subset the remaining phyla 
prevalence.df.IS <- subset(prevalence.df, Phylum %in% get_taxa_unique(IS_ps.fil, "Phylum")) #get_taxa_uniqe is used to determine the different taxa present for a particular taxonomic rank in a given phyloseq-class object
## Visualize prevalence and total read count in order to determine prevalence threshold
ggplot(prevalence.df.IS, aes(TotalAbundance, Prevalence / nsamples(IS_ps.fil),color=Phylum)) + 
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) + 
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") + facet_wrap(~Phylum) + 
  theme(legend.position="none")
## Define prevalence threshold as 10% of total samples (CHECK: chose 0.05 as it seems there is a slight natural separation at around this prevalence in the actinobacteria and proteobacteria taxa)
prevalenceThreshold.IS <- 0.10*nsamples(IS_ps.fil) #i.e. taxa have to appear in a minimum of 10% of samples or they will be removed
### Apply this prevalence filter using prune_taxa()
keepTaxa.IS <- rownames(prevalence.df.IS)[(prevalence.df.IS$Prevalence >= prevalenceThreshold.IS)]
IS_ps.fil <- prune_taxa(keepTaxa.IS,IS_ps.fil) 
### Filter taxa of the whole T2D phyloseq object using filtered IR and IS subsets
keepTaxa.T2D.fil <- c(keepTaxa.IR, keepTaxa.IS) 
T2D.fil <- prune_taxa(keepTaxa.T2D.fil, T2D.fil)
### Further filtering to allow for sparse CCA (i.e. the integrative analysis)
library("genefilter")
T2D.fil <- prune_taxa(taxa_sums(T2D.fil) > 4, T2D.fil)
T2D.fil <- filter_taxa(T2D.fil, filterfun(kOverA(40,2)),TRUE)
Y <- otu_table(T2D.fil)
Y[Y>50] <- 50 
T2D.fil@otu_table@.Data <- Y
T2D.fil
# Final phyloseq composition:
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 362 taxa and 402 samples ]
# sample_data() Sample Data:       [ 402 samples by 22 sample variables ]
# tax_table()   Taxonomy Table:    [ 362 taxa by 7 taxonomic ranks ]

##### Data analysis
#### PCoA comparing IR and IS (investigation of metagenomic separation)
pslog <- transform_sample_counts(T2D.fil, function(x) log(1 + x))
out.pcoa.log <- ordinate(pslog, method = "PCoA", distance = "bray")
evals <- out.pcoa.log$values$Eigenvalues
plot_ordination(pslog, out.pcoa.log, type = "samples", color = "IR_IS_classification") +
  labs(col = "classification") +
  coord_fixed(sqrt(evals[2] / evals[1])) 

#### PERMANOVA significance test for PCoA (tests whether the observed separation is significant)
library(microbiome)
library(ggplot2)
library(dplyr)
### Isolate relative abundances and sample metadata
T2D.fil.rel <- microbiome::transform(T2D.fil, "compositional")
T2D.fil.otu <- abundances(T2D.fil.rel)
T2D.fil.meta <- meta(T2D.fil.rel)
### Perform PERMANOVA
library(vegan)
permanova.pcoa <- adonis(t(T2D.fil.otu)~IR_IS_classification, 
                         data = T2D.fil.meta, permutations=99, method = "bray")
print(permanova.pcoa)
## Get p-value
print(as.data.frame(permanova.pcoa$aov.tab)["IR_IS_classification", "Pr(>F)"])
#### Check homogeneity of group dispersions to see if differences in group variance could be an explanation for the metagenomic separation
dist <- vegdist(t(T2D.fil.otu))
anova(betadisper(dist, T2D.fil.meta$IR_IS_classification))
#### Ivestigation of top taxa separating the groups
coef.pcoa <- coefficients(permanova.pcoa)["IR_IS_classification1",]
top.coef.pcoa <- coef.pcoa[rev(order(abs(coef.pcoa)))[1:20]]
par(mar=c(3, 14, 2, 1)) # changing margins: c(bottom, left, top, right) -> ives the number of lines of margin to be specified on the four sides of the plot.
barplot(sort(top.coef.pcoa), horiz = T, las = 1, main = "Top taxa", xlab = "PERMANOVA coefficient")
#### Bar charts for phylum abundance visualization 
### Redo IR and IS subsets so that they are consistent with the further filtering performed
IR_ps.fil <- subset_samples(T2D.fil, IR_IS_classification == "IR")
IS_ps.fil <- subset_samples(T2D.fil, IR_IS_classification == "IS")
### Bar chart of phylum in IR samples
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pIR <- plot_bar(IR_ps.fil, x = "file_name", fill = "Phylum")
pdIR <- pIR$data
IR.bp <- ggplot(pdIR, aes(x = file_name, y = Abundance, fill = Phylum)) + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  scale_fill_manual(values=cbPalette) + xlab("IR samples") + ylab("Abundance count") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
### Bar chart of phylum in IS samples
pIS <- plot_bar(IS_ps.fil, x = "file_name", fill = "Phylum")
pdIS <- pIS$data
IS.bp <- ggplot(pdIS, aes(x = file_name, y = Abundance, fill = Phylum)) + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  scale_fill_manual(values=cbPalette) + xlab("IS samples") + ylab("Abundance count") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
### Combine bar chart into one figure
library("ggpubr")
ggarrange(IR.bp, IS.bp, labels = c("A", "B"), ncol=1, nrow =2)
#### Box plots of mean abundance for each phylum
library("microbiome")
### Box plot for Actinobacteria
Actinobacteria <- subset_taxa(T2D.fil, Phylum == "Actinobacteria") # worked
Actino.mean <- colMeans(Actinobacteria@otu_table@.Data) #could log transform these means: log2(Actino.mean)
Actino.mean <- as.data.frame(Actino.mean)
Actino.mean <- t(Actino.mean)
rownames(Actino.mean) <- "Actinobacteria abundance count"
Actinobacteria@otu_table@.Data <- Actino.mean # worked
Actino.plot <- boxplot_abundance(Actinobacteria, x = "IR_IS_classification", y = "Actinobacteria abundance count",  violin = FALSE, na.rm = FALSE, show.points = FALSE) + xlab("Classification")
print(Actino.plot)
## Calculating p-values for Actinobacteria 
Actinobacteria <- subset_taxa(T2D.fil, Phylum == "Actinobacteria")
Act.IR <- subset_samples(Actinobacteria, IR_IS_classification == "IR")
Act.IS <- subset_samples(Actinobacteria, IR_IS_classification == "IS")
Wilcox.Act <- wilcox.test(Act.IR@otu_table@.Data, Act.IS@otu_table@.Data, paired = FALSE)
Wilcox.Act
### Box plot for Bacteroidetes
Bacteroidetes <- subset_taxa(T2D.fil, Phylum == "Bacteroidetes") # worked
Bacter.mean <- colMeans(Bacteroidetes@otu_table@.Data)
Bacter.mean <- as.data.frame(Bacter.mean)
Bacter.mean <- t(Bacter.mean)
rownames(Bacter.mean) <- "Bacteroidetes abundance count"
Bacteroidetes@otu_table@.Data <- Bacter.mean # worked
Bacter.plot <- boxplot_abundance(Bacteroidetes, x = "IR_IS_classification", y = "Bacteroidetes abundance count",  violin = FALSE, na.rm = FALSE, show.points = FALSE) + xlab("Classification") + ylim(c(0,20))
print(Bacter.plot)
## Calculating p-values for Bacteroidetes
Bacteroidetes <- subset_taxa(T2D.fil, Phylum == "Bacteroidetes")
Bac.IR <- subset_samples(Bacteroidetes, IR_IS_classification == "IR")
Bac.IS <- subset_samples(Bacteroidetes, IR_IS_classification == "IS")
Wilcox.Bac <- wilcox.test(Bac.IR@otu_table@.Data, Bac.IS@otu_table@.Data, paired = FALSE)
Wilcox.Bac
### Box plot for Firmicutes
Firmicutes <- subset_taxa(T2D.fil, Phylum == "Firmicutes") # worked
Firmi.mean <- colMeans(Firmicutes@otu_table@.Data)
Firmi.mean <- as.data.frame(Firmi.mean)
Firmi.mean <- t(Firmi.mean)
rownames(Firmi.mean) <- "Firmicutes abundance count"
Firmicutes@otu_table@.Data <- Firmi.mean # worked
Firmi.plot <- boxplot_abundance(Firmicutes, x = "IR_IS_classification", y = "Firmicutes abundance count",  violin = FALSE, na.rm = FALSE, show.points = FALSE) + xlab("Classification") + ylim(c(0,20))
print(Firmi.plot)
## Calculating p-values for Firmicutes 
Firmicutes <- subset_taxa(T2D.fil, Phylum == "Firmicutes")
Firmi.IR <- subset_samples(Firmicutes, IR_IS_classification == "IR")
Firmi.IS <- subset_samples(Firmicutes, IR_IS_classification == "IS")
Wilcox.Firmi <- wilcox.test(Firmi.IR@otu_table@.Data, Firmi.IS@otu_table@.Data, paired = FALSE)
Wilcox.Firmi
### Box plot for Proteobacteria
Proteobacteria <- subset_taxa(T2D.fil, Phylum == "Proteobacteria") # worked
Proteo.mean <- colMeans(Proteobacteria@otu_table@.Data)
Proteo.mean <- as.data.frame(Proteo.mean)
Proteo.mean <- t(Proteo.mean)
rownames(Proteo.mean) <- "Proteobacteria abundance count"
Proteobacteria@otu_table@.Data <- Proteo.mean # worked
Proteo.plot <- boxplot_abundance(Proteobacteria, x = "IR_IS_classification", y = "Proteobacteria abundance count",  violin = FALSE, na.rm = FALSE, show.points = FALSE) + xlab("Classification") + ylim(c(0,20))
print(Proteo.plot)
## Calculating p-values for Proteobacteria 
Proteobacteria <- subset_taxa(T2D.fil, Phylum == "Proteobacteria")
Proteo.IR <- subset_samples(Proteobacteria, IR_IS_classification == "IR")
Proteo.IS <- subset_samples(Proteobacteria, IR_IS_classification == "IS")
Wilcox.Proteo <- wilcox.test(Proteo.IR@otu_table@.Data, Proteo.IS@otu_table@.Data, paired = FALSE)
Wilcox.Proteo
### Box plot for Verrucomicrobia
Verrucomicrobia <- subset_taxa(T2D.fil, Phylum == "Verrucomicrobia") # worked
Verruco.mean <- colMeans(Verrucomicrobia@otu_table@.Data)
Verruco.mean <- as.data.frame(Verruco.mean)
Verruco.mean <- t(Verruco.mean)
rownames(Verruco.mean) <- "Verrucomicrobia abundance count"
Verrucomicrobia@otu_table@.Data <- Verruco.mean # worked
Verruco.plot <- boxplot_abundance(Verrucomicrobia, x = "IR_IS_classification", y = "Verrucomicrobia abundance count",  violin = FALSE, na.rm = FALSE, show.points = FALSE) + xlab("Classification") 
print(Verruco.plot)
## Calculating p-values for Verrucomicrobia
Verrucomicrobia <- subset_taxa(T2D.fil, Phylum == "Verrucomicrobia")
Verruco.IR <- subset_samples(Verrucomicrobia, IR_IS_classification == "IR")
Verruco.IS <- subset_samples(Verrucomicrobia, IR_IS_classification == "IS")
Wilcox.Verruco <- wilcox.test(Verruco.IR@otu_table@.Data, Verruco.IS@otu_table@.Data, paired = FALSE)
Wilcox.Verruco
#### Create one figure with all box plots for mean phylum abundance
library(ggpubr)
ggarrange(Actino.plot, Bacter.plot, Firmi.plot, Proteo.plot, Verruco.plot, labels = c("A", "B", "C", "D", "E", "F"), ncol=3, nrow =2)

############################## Integrative analysis ##############################

##### Preprocessing of the metabolomic data
#### Exclude outlier samples (determined from PCA) from metabolomic data
metabolomics.EO <- metabolomics[metabolomics$SampleID %in% excl.out.vec,]
# Now T2D.fil and the metabolomics data have the same number of samples 
#### Make the sample IDs the rownames, remove the sample ID column and transpose
rownames(metabolomics.EO) <- metabolomics.EO[,1] 
metabolomics.EO <- metabolomics.EO[,-1]
metabolomics.EO <- t(metabolomics.EO)
#### Removing samples with low abundance across many samples and transforming
keep_ix <- rowSums(metabolomics.EO == 0) <=3
metabolomics.EO.fil <- metabolomics.EO[keep_ix,]
metabolomics.EO.fil.log <- log(1 + metabolomics.EO.fil,base = 10)
#### Isolating the OTU table from T2D.fil
X <- otu_table(T2D.fil)
##### Applying the spares CCA to determine which metabolites and microbes can best explain the covariation between the datasets
library(PMA)
cca_res <- CCA(t(X), t(metabolomics.EO.fil.log), penaltyx = .15, penaltyz = .15)
cca_res
#### Obtain a data frame of these results
library(ade4)
library(factoextra)
library(magrittr)
library(genefilter)
library(ggrepel)
feature_info <- cbind(t(X[cca_res$u != 0, ]), t(metabolomics.EO.fil.log[cca_res$v != 0, ]))
