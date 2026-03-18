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

setwd("")


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

#Year
plot_NMDS <- ordinate(
  physeq = phylopilotYearmodif, 
  method = "NMDS", 
  distance = "bray")

ordplot.type<-plot_ordination(
  physeq = phylopilotYearmodif,
  ordination = plot_NMDS,
  color = "RegionMAP", shape = "Code_Year")


ordplot.type$data$Code <- factor(ordplot.type$data$Code, 
                                 levels = c("SP", "PC", "CWR", "HSC", "KWD"))
ordplot.type$data$Code_Year <- factor(ordplot.type$data$Code_Year, 
                                      levels = c("SP_2018", "SP_2019","PC_2018","PC_2019", "CWR_2018", 
                                                 "CWR_2019","HSC_2018","HSC_2019", "KWD_2018","KWD_2019"))
ordplot.type$data$Year<- factor(ordplot.type$data$Year)
levels(ordplot.type$data$Year)<- c("2018", "2019")
ordplot.type$data$Ellipse<-ordplot.type$data$Year


NMDS_Year_BySite<-
  ordplot.type +
  scale_shape_manual(values = c(1, 19, 0, 15, 2, 17, 5, 18, 7, 3))+
  labs(shape="Sites/Year", color="MAP(mm)")+
  scale_colour_gradientn(colours= c("#D22B2B","orange","blue","dark blue"),
                         values = NULL,
                         space = "Lab",
                         na.value = "grey50",
                         guide = "colourbar",
                         aesthetics = "colour")+guides(shape = guide_legend(order = 1))+
  geom_point(size = 1.5, stroke = 0.7)+
  stat_ellipse(aes(linetype = Ellipse), type = "t", level = 0.95, size = 0.8) +  
  scale_linetype_manual(values = c("2018" = "dashed", "2019" = "solid")) + 
  guides(fill=guide_legend(title="New Legend Title"))+
  theme_bw()+
  facet_wrap(~ Code,scales = "free", ncol=3)+
  guides(
    color = guide_colorbar(order = 1),   # MAP legend first
    shape = guide_legend(order = 2),     # Sites/Year legend second
    linetype = guide_legend(order = 2)   # Year legend last
  )
NMDS_Year_BySite

# For wide pdf ncol should be = 3 in the  facet_wrap
ggsave(file= "Figure_4.pdf", plot=NMDS_Year_BySite, width = 8, height = 5.5)

