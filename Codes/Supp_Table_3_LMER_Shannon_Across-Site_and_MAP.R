# DESCRIPTION: LMER OF ALPHA DIVERSITY ACROSS A GRADIENT IN MAP

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


# Spatial data 
otu.file<- read.csv("OTU_Table_NoSingleandDoubletons.csv")
tax.file<- read.csv("Tax_table_NoSingleandDouble.csv")
pilot.mapping<- read.csv("Mapping_file_2019_2020.csv")
RegionMAP<- read.csv("RegionMAP_Transect_toCombine.csv") #Data with Region MAP
pilot.mapping<-left_join(pilot.mapping,RegionMAP, by= "Sample" ) #Joining data

# Calculation start here
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

#For the graph of 2018, I did not include the phylopilot while drawing this
phylopilotA = subset_samples(phylopilot, Host== 'SCSC')
phylopilotAmodif = subset_samples(phylopilotA, Year!= 'Y_2018') #Remove 2018 samples
phylopilotAmodif = subset_samples(phylopilotAmodif, RemoveExtra != 'No') # Remove duplicate samples used for sequencing run
phylopilotAmodif = subset_samples(phylopilotAmodif, RemoveExtra_KNZ != 'No') # Remove extra samples from KNZ


# Calculation of alpha diversity Alpha Diversity - based on phyloseq
sample.richness<-estimate_richness(phylopilotAmodif,split=T, c("Observed", "Chao1", "Shannon"))
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


#Stat
MDc1b = lmer(Shannon ~ MAP2.clim + Code +  (1|Trans), data = sample.richness_2)
summary(MDc1b)
anova(MDc1b)

library(stargazer)
stargazer(MDc1b, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

