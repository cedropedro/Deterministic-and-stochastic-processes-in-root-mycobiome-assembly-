# Load the phyloseq package
library(reshape2)
library(reshape)
library(dplyr)
library(lme4)
library(lmerTest)
library(multcompView)
library(multcomp)
library(vegan)

# Between year Comparison

# OTU table
OTU.table <- read.csv("Year_OTU_TableNoSingle.csv")

# Mapping file
mapping <- read.csv("Mapping_file_YearComparison.csv")
RegionMAP<- read.csv("RegionMAP_Transect_toCombine.csv") #Data with Region MAP
mapping<-left_join(mapping,RegionMAP, by= "Sample" ) #Joining data

# Remove OTU_id column as you don't need
OTU.table<-OTU.table[,2:ncol(OTU.table)]
mapping<-mapping[2:nrow(mapping),]
head(OTU.table)

# Data transformation for Bray_curtis calculation
#hellinger-transform data converts to abundance and square root transforms it reduces impact of abundant outliers
#first step converts to relative abundance
OTU.table<- t(t(OTU.table)/colSums(OTU.table))
head(OTU.table)
#this step square root transforms
OTU.table<-sqrt(OTU.table)
head(OTU.table)

# distance matrix
## distance metric=bray
data.dist<-vegdist(t(OTU.table),method = "bray", binary = FALSE)
data.dist
data.dist2<-as.matrix(data.dist)

# Adonis function
adonis_Year <- adonis2(data.dist ~ mapping$MAP2.clim*mapping$Location*mapping$Year+mapping$Transect:mapping$Location, permutations=9999)
adonis_Year



