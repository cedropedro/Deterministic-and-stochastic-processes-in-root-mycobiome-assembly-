# Load the phyloseq package
library(phyloseq)
library(reshape2)
library(reshape)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggpubr)
library(multcompView)
library(emmeans)
library(multcomp)
library(geosphere)


#Spatial
otu.file<- read.csv("OTU_Table_NoSingleandDoubletons.csv")
tax.file<- read.csv("Tax_table_NoSingleandDouble.csv")
pilot.mapping<- read.csv("Mapping_file_2019_2020.csv")
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

#For the graph of 2018, I did not include the phylopilot while drawing this
phylopilotA = subset_samples(phylopilot, Host== 'SCSC')
phylopilotAmodif = subset_samples(phylopilotA, Year!= 'Y_2018') #Remove 2018 samples
phylopilotAmodif = subset_samples(phylopilotAmodif, RemoveExtra != 'No') # Remove duplicate samples used for sequencing run
phylopilotAmodif = subset_samples(phylopilotAmodif, RemoveExtra_KNZ != 'No') # Remove extra samples from KNZ

#SPATIAL
plot_NMDS <- ordinate(
  physeq = phylopilotAmodif, 
  method = "NMDS", 
  distance = "bray",
  trymax = 1000,
  maxit = 500)
ordplot.type8<-plot_ordination(
  physeq = phylopilotAmodif,
  ordination = plot_NMDS,
  shape = "Code", color = "RegionMAP")
ordplot.type8New<- ordplot.type8 +
  theme_bw()+
  scale_shape_manual(values = 2:16)+
  labs(shape="Sites", color="MAP (mm)")+
  scale_colour_gradientn(colours= c("#D22B2B","orange","blue","blue","dark blue"),
                         values = NULL,
                         space = "Lab",
                         na.value = "grey50",
                         guide = "colourbar",
                         aesthetics = "colour")+ guides(shape = guide_legend(order = 1))
ordplot.type8New

# For wide pdf ncol should be = 3 in the  facet_wrap
ggsave(file= "Figure_2.pdf", plot=ordplot.type8New, width = 7, height = 5)
