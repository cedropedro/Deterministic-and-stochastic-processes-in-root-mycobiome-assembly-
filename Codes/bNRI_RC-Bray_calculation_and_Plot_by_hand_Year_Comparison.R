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
Mytree<- ape::read.tree(file= 'HP_Tree_ForbNTIRC.txt')

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


# Block-level analyses

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

Block_mapping <- read.csv('Block_Mapping_file_Rmv_2018_2019.csv') # This for using sites of both year as a species pool 
                                                                  # assuming the species are often available and the environment of mediate plant recruitment

# Define sites
sites <- c("CWR", "DMR", "HSC", "HTM", "KNZ", "KWD", "MP", "PC", "RR", "SP")


# Initialize list to store results
AmongBlock_list <- list()

# Load community OTU table once
# community <- read.csv("OTU_Table_Transect_New_NoSingletons_2018_2019.csv")


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

# Save final results to a single file - Optional

#write.csv(All_AmongBlocks, file = 'Block_EachSite_Individually_Unweighted_Combined_bNRI_RCBC_Output_2018_2019.csv')

cat("All Sites combined into a single dataframe: All_AmongBlocks\n")



### Plotting

#### Within sites - between block

##Unweigthed
ComUnweightBlock<- read.csv("Block_EachSite_Individually_Unweighted_Combined_bNRI_RCBC_Output_2018_2019.csv")

#Data
Sample1<- read.csv("Block_Mapping_Sample1_Mar2024.csv")
Sample2<- read.csv("Block_Mapping_Sample2_Mar2024.csv")

##Joining dataset
Block1 <-left_join(ComUnweightBlock,Sample1, by="sample1")
Block2 <-left_join(Block1,Sample2, by="sample2")
Block2<- subset(Block2, Hostsamp1 == "SCSC")
Block2<- subset(Block2, Hostsamp2 == "SCSC")
Block2<- subset(Block2, X2Yearsample1 == "Yes") # only when using both years data
Block2<- subset(Block2, X2Yearsample1 == "Yes") # only when using both years data
output.df <- Block2
output.df['com'] = rep(NA, nrow(Block2)) #Make combination line
output.df$com[which(as.character(output.df$Codesamp1) == as.character(output.df$Codesamp2))] = 'Within'
output.df$com[which(as.character(output.df$Codesamp1) != as.character(output.df$Codesamp2))] = 'Between'
output.df['Year'] = rep(NA, nrow(Block2)) #Make combination line
output.df$Year[which(as.character(output.df$Yearsample1) == as.character(output.df$Yearsample2))] = 'WithinYear'
output.df$Year[which(as.character(output.df$Yearsample1) != as.character(output.df$Yearsample2))] = 'BetweenYear'

# Only comapring within the same year
output.df<- subset(output.df, Year == "WithinYear")

# Only comapring within the same site
output.df1<- subset(output.df, com == "Within")

# Order the processes
output.df1$Process <- factor(output.df1$Process, 
                             levels = c("Competition", "Environmental filtering", 
                                        "Dispersal limitation", "Homogenizing dispersal", 
                                        "Ecological Drift"))


# Plot
# 2018
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
  scale_x_discrete(limits=c("SP", "PC","CWR","HSC","KWD"))+
  facet_wrap(vars(Yearsample1), ncol = 1)
Block2Unweightb

#ggsave(filename="Plots_bNRIByHand/Blocks/SiteAsPool/BlockWithinSitesbetweenBlocks_Unweighted_2018.pdf", plot=Block2Unweightb,height= 3.5, width= 6.5, unit= "in")

