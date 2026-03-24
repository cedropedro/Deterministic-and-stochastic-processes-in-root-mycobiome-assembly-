# iCAMP Approach: Year Comparison

## uploading tree

library(ape)
library(picante)
library(iCAMP)
library(tidyverse)
library(dplyr)
library(reshape2)
library(dplyr)

# Calling the tree
Mytree<- ape::read.tree(file= 'HP_Tree_ForbNRIRC.txt')

## BLOCK LEVEL ANALYSIS- AMONG BLOCKS WITHIN SITES
# We are using Site as a species pool

# Read plant mapping file
#Plant_mapping <- read.csv('Block_Mapping_file_Rmv_2018_2019.csv') # This for using sites of both year as a species pool assuming the species are often available and the environment of mediate plant recruitment and to compare among years as well


# Define sites
sites <- c("CWR", "DMR", "HSC", "HTM", "KNZ", "KWD", "MP", "PC", "RR", "SP")


# Initialize lists to store data
AmongBlock_Bin18_list <- list()

# Load community OTU table once (instead of reloading it for each site)
community <- read.csv("OTU_Table_Transect_New_NoSingletons_2018_2019.csv")

# Transpose the file for merging with plant metadata
community_t <- as.data.frame(t(community[, -1]))  # Transpose
colnames(community_t) <- community[, 1]  # Set column names
community_t$Transectsample1 <- rownames(community_t)  # Add Sample column

# Join with plant mapping data
community_t_joined <- left_join(community_t, Plant_mapping[, c(3,11)], by = "Transectsample1")

community_t_joined$Transectsample1

# Ensure row names are preserved
rownames(community_t_joined) <- community_t_joined$Transectsample1
unique(community_t_joined$Codesample1)

# Loop through each site
for (site in sites) {
  
  # Subset by site
  Site_com <- subset(community_t_joined, community_t_joined$Codesample1 == site) 
  
  # Check if there are at least 2 samples
  if (nrow(Site_com) < 2) {
    cat("Skipping site:", site, "due to insufficient samples.\n")
    next
  }
  
  # Remove last two columns and transpose
  Site_com <- as.data.frame(t(Site_com[, -c((ncol(Site_com) - 1), ncol(Site_com))]))  
  
  # Remove OTUs that are absent in all samples
  Site_com_sub <- Site_com[rowSums(Site_com) > 0, ]  
  
  # Convert to presence/absence data
  ComC <- Site_com_sub %>% mutate_if(is.numeric, ~1 * (. > 0))
  ComD <- t(ComC)
  
  # Define the phylogenetic tree for the block
  SiteTip <- rownames(Site_com_sub)
  SiteTree <- keep.tip(Mytree, SiteTip)
  
  # Run iCAMP with max bin 18 - site as pool
  icamp.outUnwSite <- icamp.big(comm = ComD, tree = SiteTree, pd.wd = getwd(),
                                rand = 1000, nworker = 4, bin.size.limit = 18, 
                                phylo.metric = c("bNRI", "bNTI"), sig.index = 'SES.RC')
  # Remove unwanted files
  files_to_remove <- c("iCAMP.process.bNRIiRCa.csv", "iCAMP.iCAMP.SES.RC.detail.rda",
                       "path.rda", "pd.bin", "pd.desc", "pd.taxon.name.csv")
  
  # Construct full paths and remove only existing files
  file_paths <- file.path(getwd(), files_to_remove)
  existing_files <- file_paths[file.exists(file_paths)]
  if (length(existing_files) > 0) file.remove(existing_files)
  
  
  # Store results
  AmongBlock_Bin18_list[[paste0("icampBin18_", site)]] <- icamp.outUnwSite$bNRIiRCa
  
  # Export results to the global environment **after** the loop
  list2env(AmongBlock_Bin18_list, envir = .GlobalEnv)
  
  # Remove unwanted files
  files_to_remove <- c("iCAMP.process.bNRIiRCa.csv", "iCAMP.iCAMP.SES.RC.detail.rda",
                       "path.rda", "pd.bin", "pd.desc", "pd.taxon.name.csv")
  
  # Construct full paths and remove only existing files
  file_paths <- file.path(getwd(), files_to_remove)
  existing_files <- file_paths[file.exists(file_paths)]
  if (length(existing_files) > 0) file.remove(existing_files)
  
  # Print progress
  cat("Processed and saved:", site, "\n")
}


# Combine all elements of Transect_Bin18_list into a single data frame
All_AmongBlocks_Bin18 <- do.call(rbind, AmongBlock_Bin18_list)

# Print completion message
cat("All Sites combined into a single dataframe: All_Transect_Bin18\n")

# Save results
#write.csv(All_AmongBlocks_Bin18, file='Block_EachSite_Individually_Unweighted-icamp_output-bNRI-Bin18-NoSingletons_BothYear.csv')


# PLOTS
### Block Year comparison

## Unweigthed
ComUnweightBlock<- read.csv("Block_EachSite_Individually_Unweighted-icamp_output-bNRI-Bin18-NoSingletons_BothYear.csv")

## Mapping file
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
output.df<- subset(output.df, Year == "WithinYear")


#Within
output.df2<- subset(output.df, com == "Within")

Blockforrun <- output.df2[,c(2:8,11,15)]
any(is.na(Blockforrun$Codesamp1)) # This should be false so that there are no row missing data

BlockFinal = reshape2::melt(Blockforrun, id = c("sample1", "sample2", "Codesamp1", "Yearsample1"))
levels(BlockFinal$variable)<- c("Competition","Environmental filtering", "Dispersal limitation","Homogenizing dispersal", "Ecological Drift")
BlockFinal$Yearsample1 <- factor(BlockFinal$Yearsample1)
levels(BlockFinal$Yearsample1)<- c("2018", "2019")

Block1<- ggplot(BlockFinal, aes(fill=variable, y=value, x=Codesamp1)) + 
  geom_bar(position="fill", stat="identity")+
  theme(axis.title.x= element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))
Block1Unweight <- Block1 + 
  scale_fill_manual(values = c("light blue","#231F20","dark cyan","#800080","#B5651D")) +
  scale_y_continuous(labels = scales::percent)+
  ylab("Relative Importance")+ theme (strip.text = element_text(size = 11))+
  theme(legend.title = element_blank())+
  scale_x_discrete(limits=c("SP", "PC","CWR","HSC","KWD"))+
  facet_wrap(vars(Yearsample1), ncol = 1)
Block1Unweight

#ggsave(filename="Plots/Block/iCAMP-Block_Bin18_IndivSitesasPool_Year_Comparison_BothYearPool.pdf", plot=Block1Unweight,height= 7.5, width= 5.5, unit= "in")



