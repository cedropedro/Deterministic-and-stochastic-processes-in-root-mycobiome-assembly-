# Load the phyloseq package
library(phyloseq)
library(reshape2)
library(reshape)
library(dplyr)
library(lme4)
library(lmerTest)
library(multcompView)
library(emmeans)
library(multcomp)



#year Comparison
otu.file<- read.csv("Year_OTU_TableNoSingle.csv")
tax.file<- read.csv("Year_Tax_NoSingle.csv")
pilot.mapping<- read.csv("Mapping_file_YearComparison.csv")
RegionMAP<- read.csv("RegionMAP_Transect_toCombine.csv") #Data with Region MAP
pilot.mapping<-left_join(pilot.mapping,RegionMAP, by= "Sample" ) #Joining data

#Calculation start here
tax.file$X<-factor(tax.file$X)
otu.file$OTU_ID<-factor(otu.file$OTU_ID)

#phyloseq only recognize matrices, this command transform data from data.frame to matrix
tax.file1<-as.matrix(tax.file, rownames.force = NA)
otu.file1<-data.matrix(otu.file, rownames.force = NA)

#phyloseq needs the same row name for tax and OTU table, make sure they are the same 
rownames(tax.file1) <- rownames(otu.file1)

#this command help you define what is a tax and OTU and then combine them into a phyloseq object 
OTU = otu_table(otu.file1, taxa_are_rows = TRUE)
TAX = tax_table(tax.file1)
pilot.phyloseq = phyloseq(OTU, TAX)
rownames(pilot.mapping) <- sample_names(pilot.phyloseq)
sampledata<-sample_data(pilot.mapping)
pilot.phyloseq= merge_phyloseq(pilot.phyloseq, sampledata)

# remove the OTU_ID row 
phylopilot = subset_samples(pilot.phyloseq, sample_names(pilot.phyloseq) != 'OTU_ID')

#Phyloseq boject
phylopilotA = subset_samples(phylopilot, Host== 'SCSC')
phylopilotYear = subset_samples(phylopilotA, X2Year== 'Yes' )
phylopilotYearmodif = subset_samples(phylopilotYear, In2YearAnalysis!= 'No' )
phylopilotYearmodif = subset_samples(phylopilotYear, Sample!= 'HPC1035' ) # remove this as it has very low number of reads

# Calculation of alpha diversity Alpha Diversity - based on phyloseq
sample.richness<-estimate_richness(phylopilotYearmodif,split=T, c("Observed", "Chao1", "Shannon"))
H <- sample.richness$Shannon
S1 <- sample.richness$Observed
S <- log(S1)
evenness <- H/S
evenness
sample.richness$Evenness = evenness
sample.richness
sample.richness$Evenness[sample.richness$Evenness == "NaN"] <- 1

# Joining the metadata to the alpha diversity calculated
sample.richness$Sample <- rownames(sample.richness)
sample.richness_2 <- left_join(sample.richness,pilot.mapping )

# LMER - Shannon
MDc1 = lmer(Shannon ~ MAP2.clim*Code*Year + (1|Transect), data = sample.richness_2)
summary(MDc1)
anova(MDc1)







