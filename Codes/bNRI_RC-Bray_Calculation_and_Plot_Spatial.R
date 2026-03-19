#DESCRIPTION: "bNRI and RC-Bray Calculation By hand

# Load packages
library(ape)
library(picante)
library(iCAMP)
library(tidyverse)
library(dplyr)
library(reshape2)
library(phylocomr)
library(phylobase)
library(phylotools)
library(phytools)


# Calling the tree
Mytree<- ape::read.tree(file= 'HP_Tree_ForbNRIRC.txt')

# Site - level analyses
## Calculating bNRI and RC-Bray


# Define function to ensure only unique pairwise comparisons are kept
filter_unique_comparisons <- function(df) {
  # Remove self-comparisons
  df <- df[df$Name != df$variable, ]
  
  # Remove duplicated pairwise comparisons ("a, b" vs. "b, a")
  df <- df[!duplicated(t(apply(df[, c('Name', 'variable')], 1, sort))), ]
  
  return(df)
}

# SITE LEVEL ANALYSIS ;  Comaprison Among site

## calling the data and transforming it to be used in teh analysis 

community <- read.csv("OTU_Table_Site_icamp_NoSingletons_2019.csv")
Community<- community[,-1]
ComB <- t(Community)
colnames(ComB) = community[,1]

### For transformed calculation of MPD
Communityb<- t(t(Community)/colSums(Community)) # this for the relative abundance
Communityb<-sqrt(Communityb)# this step square root transforms
ComBtransformed <- t(Communityb)
colnames(ComBtransformed) = community[,1]

# For presence/absence Analysis - Unweighted data
ComC<- as.data.frame(Community) %>% mutate_if(is.numeric, ~1 * (. > 0))
ComD <- t(ComC)
colnames(ComD) = community[,1]


# Generate tree and phylogenetic distance matrix
Tip <- colnames(ComD)
Tree <- keep.tip(Mytree, Tip) # ro only keeps the tips corresponding to community matrix
Phylodist2 <- cophenetic(Tree)

# Run bNRI analysis
bNRI <- bNRIn.p(comm = ComB, dis = Phylodist2, nworker = 4, memo.size.GB = 50,
                weighted = TRUE, rand = 1000, output.bMPD = TRUE, 
                sig.index = 'bNRI', unit.sum = NULL, 
                correct.special = TRUE, detail.null = FALSE, special.method = 'MPD',
                dirichlet = TRUE)

# Format and filter bNRI output
NRIBeta <- bNRI$index
BetaMP <- bNRI$betaMPD.obs
NRIBeta_df <- as.data.frame(NRIBeta)
NRIBeta_df$Name <- rownames(NRIBeta_df)
NRIBeta_melt <- melt(NRIBeta_df, id = 'Name')
NRIBeta_melt <- filter_unique_comparisons(NRIBeta_melt)

# Run RCPC analysis
RCMSB2 <- RC.pc(comm = ComB, rand = 1000, nworker = 4, weighted = TRUE, sig.index = 'RC' )
RCBC_df <- as.data.frame(RCMSB2$index)
RCBC_df$Name <- rownames(RCBC_df)
RCBC_melt <- melt(RCBC_df, id = 'Name')
RCBC_melt <- filter_unique_comparisons(RCBC_melt)

# Combine the filtered data frames
Combined_df <- full_join(NRIBeta_melt, RCBC_melt, by = c('Name', 'variable'))
colnames(Combined_df)<- c("sample1","sample2","bNRI","RC")


# Add 'Process' column based on conditions
Combined_df <- Combined_df %>%
  mutate(Process = case_when(
    bNRI > 2 ~ "Competition",
    bNRI < -2 ~ "Environmental filtering",
    bNRI >= -2 & bNRI <= 2 & RC > 0.95 ~ "Dispersal limitation",
    bNRI >= -2 & bNRI <= 2 & RC < -0.95 ~ "Homogenizing dispersal",
    bNRI >= -2 & bNRI <= 2 & RC >= -0.95 & RC <= 0.95 ~ "Ecological Drift"))

# Add contribution column to make plotting easier
Combined_df <- Combined_df %>%
  mutate(Process_contribution = 1)

# Export final combined result (optional)
#write.csv(Combined_df, file = 'SITE_Unweighted_Combined_bNRI_RCBC_Output.csv')

cat('Combined file successfully saved as Combined_bNRI_RCBC_Output.csv\n')



## Plotting

ComUnweighSite<- read.csv("SITE_Unweighted_Combined_bNRI_RCBC_Output.csv")

Sample1<- read.csv("Sites_Mapping_Sample1.csv")
Sample2<- read.csv("Sites_Mapping_Sample2.csv")

##Joining dataset
Site1 <-left_join(ComUnweighSite,Sample1, by="sample1")
Site2 <-left_join(Site1,Sample2, by="sample2")
Site2<- subset(Site2, Hostsamp1 == "SCSC")
Site2<- subset(Site2, Hostsamp2 == "SCSC")
any(is.na(Site2$Pooling)) # This should be false so that there are no row missing data

# Order the processes
Site2$Process <- factor(Site2$Process, 
                        levels = c("Competition", "Environmental filtering", 
                                   "Dispersal limitation", "Homogenizing dispersal", 
                                   "Ecological Drift"))

# Plot
Sites1<- ggplot(Site2, aes(fill=Process, y=Process_contribution, x=Pooling)) + 
  geom_bar(position="fill", stat="identity")+
  theme(axis.title.x= element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))
Sites1Unweight <- Sites1 + 
  scale_fill_manual(values = c("#231F20","dark cyan","#800080","#B5651D")) +
  scale_y_continuous(labels = scales::percent)+
  ylab("Percent pairwise comparison")+ theme (strip.text = element_text(size = 11))+
  theme(legend.title = element_blank())
Sites1Unweight

#ggsave(filename="Plots_bNRIByHand/Site_Region_asPool/bNRI_RC-Sites_2019asPool.pdf", plot=Sites1Unweight,height= 4, width= 4, unit= "in")


# BLOCK LEVEL ANALYSIS- AMONG BLOCKS WITHIN SITES

## Site as pool
### Calculating bNRI and RC_Bray

library(dplyr)
library(reshape2)

# Define filter function to handle a vs b = b vs a and a vs a cases
filter_unique_comparisons <- function(df) {
  df <- df[df$Name != df$variable, ]  # Remove self-comparisons (a vs a)
  df <- df[!duplicated(t(apply(df[, c('Name', 'variable')], 1, sort))), ]
  return(df)
}

# Read plant mapping file
#Block_mapping <- read.csv('Block_Mapping_file_Rmv.csv')

# Define sites
sites <- c("CWR", "DMR", "HSC", "HTM", "KNZ", "KWD", "MP", "PC", "RR", "SP")


# Initialize list to store results
AmongBlock_list <- list()

# Load community OTU table once
#community <- read.csv("OTU_Table_Transect_New_NoSingletons_2019.csv")


# Transpose the file for merging with plant metadata
community_t <- as.data.frame(t(community[, -1]))
colnames(community_t) <- community[, 1]
community_t$Transectsample1 <- rownames(community_t)

# Join with plant mapping data
community_t_joined <- left_join(community_t, Block_mapping[, c(3,11)], by = "Transectsample1")

community_t_joined$Transectsample1 # this should go to 87 for analysis that consider 2018 and 2019
# Ensure row names are preserved
rownames(community_t_joined) <- community_t_joined$Transectsample1

# Loop through each site
for (site in sites) {
  Site_com <- subset(community_t_joined, community_t_joined$Codesample1 == site)
  
  if (nrow(Site_com) < 2) {
    cat("Skipping site:", site, "due to insufficient samples.\n")
    next
  }
  
  Site_com <- as.data.frame(t(Site_com[, -c((ncol(Site_com) - 1), ncol(Site_com))]))
  Site_com_sub <- Site_com[rowSums(Site_com) > 0, ]
  ComD <- t(Site_com_sub)
  
  SiteTip <- rownames(Site_com_sub)
  SiteTree <- keep.tip(Mytree, SiteTip)
  Phylodist2 <- cophenetic(SiteTree)
  
  # Run bNRI analysis
  bNRI <- bNRIn.p(comm = ComD, dis = Phylodist2, nworker = 4, memo.size.GB = 50,
                  weighted = FALSE, rand = 1000, output.bMPD = TRUE, 
                  sig.index = 'bNRI', unit.sum = NULL, 
                  correct.special = TRUE, detail.null = FALSE, special.method = 'MPD',
                  dirichlet = TRUE)
  
  # Process bNRI output
  NRIBeta_df <- as.data.frame(bNRI$index)
  NRIBeta_df$Name <- rownames(NRIBeta_df)
  NRIBeta_melt <- melt(NRIBeta_df, id = 'Name')
  NRIBeta_melt <- filter_unique_comparisons(NRIBeta_melt)
  
  # Run RCPC analysis
  RCMSB2 <- RC.pc(comm = ComD, rand = 1000, nworker = 4, weighted = FALSE)
  RCBC_df <- as.data.frame(RCMSB2$index)
  RCBC_df$Name <- rownames(RCBC_df)
  RCBC_melt <- melt(RCBC_df, id = 'Name')
  RCBC_melt <- filter_unique_comparisons(RCBC_melt)
  
  # Combine the filtered data frames
  Combined_df <- full_join(NRIBeta_melt, RCBC_melt, by = c('Name', 'variable'))
  colnames(Combined_df) <- c("sample1", "sample2", "bNRI", "RC")
  
  # Add 'Process' column
  Combined_df <- Combined_df %>%
    mutate(Process = case_when(
      bNRI > 2 ~ "Competition",
      bNRI < -2 ~ "Environmental filtering",
      bNRI >= -2 & bNRI <= 2 & RC > 0.95 ~ "Dispersal limitation",
      bNRI >= -2 & bNRI <= 2 & RC < -0.95 ~ "Homogenizing dispersal",
      TRUE ~ "Ecological Drift"
    )) %>%
    mutate(Process_contribution = 1)  # Add contribution column
  
  # Store results in the list
  AmongBlock_list[[site]] <- Combined_df
  
  cat("Processed and saved:", site, "\n")
}

# Combine all data frames into one
All_AmongBlocks <- do.call(rbind, AmongBlock_list)

# Save final results to a single file
#write.csv(All_AmongBlocks, file = 'Block_EachSite_Individually_Unweighted_Combined_bNRI_RCBC_Output.csv')

cat("All Sites combined into a single dataframe: All_AmongBlocks\n")



### Plotting

#### Within sites - between block

##Unweigthed
ComUnweightBlock<- read.csv("Block_EachSite_Individually_Unweighted_Combined_bNRI_RCBC_Output.csv")

#Data
Sample1<- read.csv("Block_Mapping_Sample1_Mar2024.csv")
Sample2<- read.csv("Block_Mapping_Sample2_Mar2024.csv")

##Joining dataset
Block1 <-left_join(ComUnweightBlock,Sample1, by="sample1")
Block2 <-left_join(Block1,Sample2, by="sample2")
Block2<- subset(Block2, Hostsamp1 == "SCSC")
Block2<- subset(Block2, Hostsamp2 == "SCSC")
Block2<- subset(Block2, Yearsample1 != "Y_2018")
Block2<- subset(Block2, Yearsample2 != "Y_2018")
output.df <- Block2
output.df['com'] = rep(NA, nrow(Block2)) #Make combination line
output.df$com[which(as.character(output.df$Codesamp1) == as.character(output.df$Codesamp2))] = 'Within'
output.df$com[which(as.character(output.df$Codesamp1) != as.character(output.df$Codesamp2))] = 'Between'
output.df['comTrans'] = rep(NA, nrow(Block2)) #Make combination line
output.df$comTrans[which(as.character(output.df$Transectsamp1) == as.character(output.df$Transectsamp2))] = 'WithinTrans'
output.df$comTrans[which(as.character(output.df$Transectsamp1) != as.character(output.df$Transectsamp2))] = 'BetweenTrans'

output.df1<- subset(output.df, com == "Within")

# Order the processes
output.df1$Process <- factor(output.df1$Process, 
                             levels = c("Competition", "Environmental filtering", 
                                        "Dispersal limitation", "Homogenizing dispersal", 
                                        "Ecological Drift"))


#Plot
Blocks1b<- ggplot(output.df1, aes(fill=Process, y=Process_contribution, x=Codesamp1)) + 
  geom_bar(position="fill", stat="identity")+
  theme(axis.title.x= element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))
Block2Unweightb <- Blocks1b + 
  scale_fill_manual(values = c("#231F20","dark cyan","#800080","#B5651D")) +
  scale_y_continuous(labels = scales::percent)+
  ylab("Relative Importance")+ theme (strip.text = element_text(size = 11))+
  theme(legend.title = element_blank())+
  scale_x_discrete(limits=c("SP","HTM", "PC","CWR","DMR","RR","KNZ","HSC","KWD","MP"))
Block2Unweightb


#ggsave(filename="Plots_bNRIByHand/Blocks/SiteAsPool/BlockWithinSitesbetweenBlocks_Unweighted.pdf", plot=Block2Unweightb,height= 3.5, width= 6.5, unit= "in")


# PLANT LEVEL ANALYSIS: Among Plant within sites

## Site as pool
### Calculating bNRI-RCBray


filter_unique_comparisons <- function(df) {
  # Ensure columns exist before proceeding
  if (!all(c('Name', 'variable') %in% colnames(df))) {
    stop('Dataframe must contain columns: Name and variable')
  }
  
  # Remove self-comparisons (a vs a)
  df <- df[df$Name != df$variable, ]
  
  # Remove duplicated comparisons (a vs b = b vs a)
  df <- df[!duplicated(t(apply(df[, c('Name', 'variable')], 1, sort))), ]
  
  return(df)
}


# Read plant mapping file
Plant_mapping <- read.csv('Plant_Mapping_file_Rmv.csv')

# Define sites
sites <- c("CWR", "DMR", "HSC", "HTM", "KNZ", "KWD", "MP", "PC", "RR", "SP")

# Initialize lists to store data
Transect_list <- list()

# Loop through each site
for (site in sites) {
  
  # Load the site OTU table
  site_path <- paste0('OTU-Table-with-Actual-Tip-Label-', site, 'Nosingleton_2019.csv')
  community <- read.csv(site_path)
  
  # Transpose and prepare for merging with plant metadata
  community_t <- as.data.frame(t(community[, -1]))
  colnames(community_t) <- community[, 1]
  community_t$Sample <- rownames(community_t)
  
  # Join with plant mapping data
  community_t_joined <- dplyr::left_join(community_t, Plant_mapping[, 1:2], by = "Sample")
  rownames(community_t_joined) <- community_t_joined$Sample
  
  if (ncol(community_t_joined) <= 2) {
    cat('Skipping site', site, 'due to missing or insufficient data after join.\n')
    next
  }
  
  # Remove metadata columns and transpose
  Site_com_sub <- as.data.frame(t(community_t_joined[, -c((ncol(community_t_joined) - 1), ncol(community_t_joined))]))
  Site_com_sub <- Site_com_sub[rowSums(Site_com_sub) > 0, ]  
  ComD <- t(Site_com_sub)
  
  if (nrow(ComD) < 2) {
    cat('Skipping site', site, 'due to insufficient samples.\n')
    next
  }
  
  SiteTip <- rownames(Site_com_sub)
  SiteTree <- keep.tip(Mytree, SiteTip)
  Phylodist2 <- cophenetic(SiteTree)
  
  # Calculate bNRI
  bNRI <- bNRIn.p(comm = ComD, dis = Phylodist2, nworker = 4, memo.size.GB = 50,
                  weighted = TRUE, rand = 1000, output.bMPD = TRUE, 
                  sig.index = 'bNRI', unit.sum = NULL, 
                  correct.special = TRUE, detail.null = FALSE, special.method = 'MPD',
                  dirichlet = TRUE)
  
  NRIBeta_df <- as.data.frame(bNRI$index)
  NRIBeta_df$Name <- rownames(NRIBeta_df)
  NRIBeta_melt <- reshape2::melt(NRIBeta_df, id = 'Name')
  NRIBeta_melt <- filter_unique_comparisons(NRIBeta_melt)
  
  # Calculate RCBC
  RCMSB2 <- RC.pc(comm=ComD, rand=1000, nworker=4, weighted=FALSE)
  RCBC_df <- as.data.frame(RCMSB2$index)
  RCBC_df$Name <- rownames(RCBC_df)
  RCBC_melt <- reshape2::melt(RCBC_df, id = 'Name')
  RCBC_melt <- filter_unique_comparisons(RCBC_melt)
  
  # Combine dataframes
  Combined_df <- dplyr::full_join(NRIBeta_melt, RCBC_melt, by = c('Name', 'variable'))
  colnames(Combined_df) <- c('sample1', 'sample2', 'bNRI', 'RC')
  
  # Define ecological processes
  Combined_df <- Combined_df %>%
    dplyr::mutate(Process = dplyr::case_when(
      bNRI > 2 ~ "Competition",
      bNRI < -2 ~ "Environmental filtering",
      bNRI >= -2 & bNRI <= 2 & RC > 0.95 ~ "Dispersal limitation",
      bNRI >= -2 & bNRI <= 2 & RC < -0.95 ~ "Homogenizing dispersal",
      TRUE ~ "Ecological Drift"
    )) %>%
    dplyr::mutate(Process_contribution = 1)
  
  Transect_list[[site]] <- Combined_df
  
  cat('Processed site:', site, '\n')
}

# Combine all sites into a single dataframe
All_Transect <- do.call(rbind, Transect_list)

cat('All sites combined into a single dataframe\n')
#write.csv(All_Transect, file='All-Plant-WithinSites-Level-Analysis_SiteAsPool.csv')


### Plotting
#### Data

##Unweigthed
ComUnweightPlant<- read.csv("All-Plant-WithinSites-Level-Analysis_SiteAsPool_Unweighted.csv")

Sample1<- read.csv("Plant_Mapping_Sample1.csv")
Sample2<- read.csv("Plant_Mapping_Sample2.csv")

##Joining dataset
Plant1 <-left_join(ComUnweightPlant,Sample1, by="sample1")
Plant2 <-left_join(Plant1,Sample2, by="sample2")
Plant2<- subset(Plant2, Hostsample1 == "SCSC")
Plant2<- subset(Plant2, Hostsample2 == "SCSC")
output.df <- Plant2
output.df['com'] = rep(NA, nrow(Plant2)) #Make combination line
output.df$com[which(as.character(output.df$Codesample1) == as.character(output.df$Codesample2))] = 'Within'
output.df$com[which(as.character(output.df$Codesample1) != as.character(output.df$Codesample2))] = 'Site'
output.df['comTrans'] = rep(NA, nrow(Plant2)) #Make combination line
output.df$comTrans[which(as.character(output.df$Transectsample1) == as.character(output.df$Transectsample2))] = 'WithinTrans'
output.df$comTrans[which(as.character(output.df$Transectsample1) != as.character(output.df$Transectsample2))] = 'BetweenTrans'




#### Within sites - between block
output.df1<- subset(output.df, com == "Within")
output.df2b<- subset(output.df1, comTrans == "BetweenTrans")

# Order the processes
output.df2b$Process <- factor(output.df2b$Process, 
                              levels = c("Competition", "Environmental filtering", 
                                         "Dispersal limitation", "Homogenizing dispersal", 
                                         "Ecological Drift"))


#Plot
Plants1b<- ggplot(output.df2b, aes(fill=Process, y=Process_contribution, x=Codesample1)) + 
  geom_bar(position="fill", stat="identity")+
  theme(axis.title.x= element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))
Plant2Unweightb <- Plants1b + 
  scale_fill_manual(values = c("light blue","#231F20","dark cyan","#800080","#B5651D")) +
  scale_y_continuous(labels = scales::percent)+
  ylab("Relative Importance")+ theme (strip.text = element_text(size = 11))+
  theme(legend.title = element_blank())+
  scale_x_discrete(limits=c("SP","HTM", "PC","CWR","DMR","RR","KNZ","HSC","KWD","MP"))
Plant2Unweightb

#ggsave(filename="Plots_bNRIByHand/Plants/SiteAsPool/PlantWithinSitesbetweenBlocks_Unweighted.pdf", plot=Plant2Unweightb,height= 3.5, width= 6.5, unit= "in")



#### Within sites - within block

output.df1<- subset(output.df, com == "Within")
output.df2a<- subset(output.df1, comTrans == "WithinTrans")

# Order the processes
output.df2a$Process <- factor(output.df2a$Process, 
                              levels = c("Competition", "Environmental filtering", 
                                         "Dispersal limitation", "Homogenizing dispersal", 
                                         "Ecological Drift"))


unique(output.df2a$Process)
#Plot
Plants1a<- ggplot(output.df2a, aes(fill=Process, y=Process_contribution, x=Codesample1)) + 
  geom_bar(position="fill", stat="identity")+
  theme(axis.title.x= element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))
Plant2Unweighta <- Plants1a + 
  scale_fill_manual(values = c("#231F20","dark cyan","#800080","#B5651D")) +
  scale_y_continuous(labels = scales::percent)+
  ylab("Relative Importance")+ theme (strip.text = element_text(size = 11))+
  theme(legend.title = element_blank())+
  scale_x_discrete(limits=c("SP","HTM", "PC","CWR","DMR","RR","KNZ","HSC","KWD","MP"))
Plant2Unweighta

#ggsave(filename="Plots_bNRIByHand/Plants/SiteAsPool/PlantWithinSitesWithinBlocks_Unweighted.pdf", plot=Plant2Unweighta,height= 3.5, width= 6.5, unit= "in")
