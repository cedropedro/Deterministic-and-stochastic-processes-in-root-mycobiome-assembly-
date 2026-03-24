# iCAMP Approach: Spatial Analysis

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

## SITE LEVEL ANALYSIS: AMONG SITES ## 

# Remove unwanted files - always start with this if you have previously run iCAMP in your working directory
files_to_remove <- c("iCAMP.process.bNRIiRCa.csv", "iCAMP.iCAMP.SES.RC.detail.rda",
                     "path.rda", "pd.bin", "pd.desc", "pd.taxon.name.csv")

# Construct full paths and remove only existing files
file_paths <- file.path(getwd(), files_to_remove)
existing_files <- file_paths[file.exists(file_paths)]
if (length(existing_files) > 0) file.remove(existing_files)

#Site level
community <- read.csv("OTU_Table_Site_icamp_NoSingletons_2019.csv")
Community<- community[,-1]

# Weighted data
ComB <- t(Community)
colnames(ComB) = community[,1]

# Unweighted data
ComC<- Community %>% mutate_if(is.numeric, ~1 * (. > 0))
ComD <- t(ComC)
colnames(ComD) = community[,1]

# changing some stuff
rand.time=1000
nworker=4

# Define the phylogenetic tree for the block
Tip <- colnames(ComD)
Tree <- keep.tip(Mytree, Tip)


#Unweigthed
icamp.outUnw=icamp.big(comm=ComD,tree=Tree,pd.wd=getwd(),
                       rand=1000, nworker=4,bin.size.limit=18, 
                       phylo.metric= c("bNRI", "bNTI"),sig.index= 'SES.RC')

#write.csv(icamp.outUnw$bNRIiRCa, file=' Site_RegionAsPool_Unweighted-combined-icamp_output-bNRI-Bin18-NoSingletons.csv')


## BLOCK LEVEL ANALYSIS- AMONG BLOCKS WITHIN SITES
# We are using Site as a species pool

# Read plant mapping file
Plant_mapping <- read.csv('Block_Mapping_file_Rmv_2018_2019.csv') # This for using sites of both year as a species pool assuming the species are often available and the environment of mediate plant recruitment and to compare among years as well

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


# PLANT LEVEL ANALYSIS: Among Plant within sites

## Making plant dataset for each site individually using the singleton file
Plant_mapping <- read.csv('Plant_Mapping_file.csv') # Plant mapping file
Otu_table <- read.csv('OTU-Table-with-Actual-Tip-Label_NoSingletons.csv') # Site OTU table

# Transposing the file so it could be joind to the plant metadata and Site subset
Otu_table_t<-as.data.frame(t(Otu_table[,-1])) # transposing the file
colnames(Otu_table_t) <- Otu_table[,1] 
Otu_table_t$Sample<- rownames(Otu_table_t)

# Joining
Otu_table_t_joined<- left_join(Otu_table_t,Plant_mapping[,c(1,7,10:11)], by= "Sample" )
rownames(Otu_table_t_joined) <- rownames(Otu_table_t)
Otu_table_t_joined<-subset(Otu_table_t_joined, Otu_table_t_joined$Hostsample1 != "Corn") # To remove corn sample
Otu_table_t_joined<-subset(Otu_table_t_joined, Otu_table_t_joined$Yearsample1 != "Y_2018") # To remove sample from

unique(Otu_table_t_joined$Codesample1)


#Now making subset
Sites <- c("CWR", "DMR", "HSC", "HTM", "KNZ", "KWD","MP", "PC", "RR", "SP")
Site_tables <- list()  # Initialize list to store results


# Sub-setting
for (Site in Sites) {
  
  # Get the site OTU dynamically
  Site_com <- subset(Otu_table_t_joined, Otu_table_t_joined$Codesample1 == Site)
  
  Site_com <- as.data.frame(t(Site_com[, -c((ncol(Site_com) - 3),(ncol(Site_com) - 2),(ncol(Site_com) - 1), ncol(Site_com))]))  # Transpose and remove last 4 columns: Sample, Site, Host, Year
  Site_com_sub <- Site_com[rowSums(Site_com) > 0, ] # Remove OTUs absent in all samples
  
  #Store results
  Site_tables[[paste0(Site)]] <- Site_com_sub
  
  #Save each file separately
  write.csv(Site_com_sub, file = paste0('OTU-Table-with-Actual-Tip-Label-', Site, 'Nosingleton.csv') )
  
  # Print progress
  cat("Data stored for site:", "Otu_table",Site, "\n")
}

#Convert list to object in the environment
list2env(Site_tables, envir = .GlobalEnv)



## running the loop for plant level within site - Because the loop does not work, each time, I have change what sites was. 
 ### Clear code could be  be made. However in the main directory we have to delete:
"iCAMP.process.bNRIiRCa.csv", "iCAMP.iCAMP.SES.RC.detail.rda", "path.rda", "pd.bin","pd.desc", "pd.taxon.name.csv" 

## Define each site
sites <- c("CWR", "DMR", "HSC", "HTM", "KNZ", "KWD","MP", "PC", "RR", "SP")

# Initialize lists to store data for each site
Site_Bin18_list <- list()

# Loop through each site
for (site in sites) {
  
  # Get the site OTU dynamically
  site_path <- paste0('OTU-Table-with-Actual-Tip-Label-', site, 'Nosingleton_2019.csv')
  community <- read.csv(site_path)  # Remove the quotes around site_path
  Community <- community[, -1]  # Remove the first column
  Community <- Community[rowSums(Community) > 0, ]  # Remove OTUs absent in all samples
  ComB <- t(Community)  # Transpose the matrix
  colnames(ComB) <- community[rowSums(community[, -1]) > 0, 1]  # Set the first column as column names
  ComC<- Community %>% mutate_if(is.numeric, ~1 * (. > 0))
  ComD <- t(ComC)
  colnames(ComD) <- community[rowSums(community[, -1]) > 0, 1]  # Set the first column as column names
  
  # Define the tree
  SiteTip <- community[rowSums(community[, -1]) > 0, 1]
  SiteTree <- keep.tip(Mytree, SiteTip)
  
  # Set parameters
  rand.time = 1000
  nworker = 4
  
  # Run iCAMP with max bin 18 
  icamp.outUnw = icamp.big(comm=ComD, tree=Mytree, pd.wd=getwd(),
                           rand=1000, nworker=4, bin.size.limit=18, 
                           phylo.metric= c("bNRI", "bNTI"), sig.index= 'SES.RC')
  
  # Store results
  Site_Bin18_list[[paste0("icampBin18", site)]] <- icamp.outUnw$bNRIiRCa
  
  # Convert each element in Bin18 and Bin100 list to individual objects in the global environment
  list2env(Site_Bin18_list, envir = .GlobalEnv)
  
  # Define a list of files to remove - this 
  files_to_remove <- c("iCAMP.process.bNRIiRCa.csv", "iCAMP.iCAMP.SES.RC.detail.rda", "path.rda", 
                       "pd.bin", "pd.desc", "pd.taxon.name.csv")
  
  # Construct full paths
  file_paths <- file.path(getwd(), files_to_remove)
  
  # Remove files if they exist
  file.remove(file_paths[file.exists(file_paths)])
  
  
  # Print progress
  cat("Data stored for site:", site, "\n")
}


Combining the data

# Combining the data
Bin18_Together <- rbind(icampBin18CWR,icampBin18DMR, icampBin18SP, icampBin18PC, icampBin18HSC, icampBin18KNZ, icampBin18KWD, icampBin18HTM, icampBin18MP, icampBin18RR)

#write.csv(Bin18_Together, file='All- EachSiteIndividually-Plants-Unweighted-icamp_output-bNRI-Bin18-NoSingletons_2019_BothYears.csv')



# PLOTS

## Sites

#Site
##Unweigthed
ComUnweighSite<- read.csv(" Site_RegionAsPool_Unweighted-combined-icamp_output-bNRI-Bin18-NoSingletons.csv")

# Mapping file
Sample1<- read.csv("Sites_Mapping_Sample1_copy.csv")
Sample2<- read.csv("Sites_Mapping_Sample2_copy.csv")

##Joining dataset
Site1 <-left_join(ComUnweighSite,Sample1, by="sample1")
Site2 <-left_join(Site1,Sample2, by="sample2")
Site2<- subset(Site2, Hostsamp1 == "SCSC")
Site2<- subset(Site2, Hostsamp2 == "SCSC")
Siteforrun <- Site2[,c(2:8,22)]
#Siteforrun <- Site2[,c(1:7,21)]
any(is.na(Siteforrun$Pooling)) # This should be false so that there are no row missing data

Sites = melt(Siteforrun, id = c("sample1", "sample2", "Pooling"))
levels(Sites$variable)<- c("Competition","Environmental filtering", "Dispersal limitation","Homogenizing dispersal", "Ecological Drift")

Sites1<- ggplot(Sites, aes(fill=variable, y=value, x=Pooling)) + 
  geom_bar(position="fill", stat="identity")+
  theme(axis.title.x= element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))
Sites1Unweight <- Sites1 + 
  scale_fill_manual(values = c("light blue","#231F20","dark cyan","#800080","#B5651D")) +
  scale_y_continuous(labels = scales::percent)+
  ylab("Percent pairwise comparison")+ theme (strip.text = element_text(size = 11))+
  theme(legend.title = element_blank())
Sites1Unweight


#ggsave(filename="iCAMP-Sites_Bin18_2019_IndivSite_asPool.pdf", plot=Sites1Unweight,height= 4, width= 4, unit= "in")


## Block 

##Unweigthed

ComUnweightBlock<- read.csv("Block_EachSite_Individually_Unweighted-icamp_output-bNRI-Bin18-NoSingletons_BothYear.csv")

# Mapping file
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
output.df['Year'] = rep(NA, nrow(Block2)) #Make combination line
output.df$Year[which(as.character(output.df$Yearsample1) == as.character(output.df$Yearsample2))] = 'WithinYear'
output.df$Year[which(as.character(output.df$Yearsample1) != as.character(output.df$Yearsample2))] = 'BetweenYear'
output.df<- subset(output.df, Year == "WithinYear")


#Within
output.df2<- subset(output.df, com == "Within")

Blockforrun <- output.df2[,c(2:8,15)]
any(is.na(Blockforrun$Codesamp1)) # This should be false so that there are no row missing data

BlockFinal = reshape2::melt(Blockforrun, id = c("sample1", "sample2", "Codesamp1"))
levels(BlockFinal$variable)<- c("Competition","Environmental filtering", "Dispersal limitation","Homogenizing dispersal", "Ecological Drift")

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
  scale_x_discrete(limits=c("SP","HTM", "PC","CWR","DMR","RR","KNZ","HSC","KWD","MP"))
Block1Unweight
#ggsave(filename="Plots/Block/iCAMP-Block_Bin18_2019IndivSitesasPool_Both_Year.pdf", plot=Block1Unweight,height= 3.5, width= 6.5, unit= "in")


## Individual plants

##Unweigthed
# Site as species pool
ComUnweightPlant<- read.csv("All- EachSiteIndividually-Plants-Unweighted-icamp_output-bNRI-Bin18-NoSingletons_2019_BothYears.csv")

# Mapping files
Sample1<- read.csv("Plant_Mapping_Sample1_copy.csv")
Sample2<- read.csv("Plant_Mapping_Sample2_copy.csv")

##Joining dataset
Plant1 <-left_join(ComUnweightPlant,Sample1, by="sample1")
Plant2 <-left_join(Plant1,Sample2, by="sample2")
Plant2<- subset(Plant2, Hostsample1 == "SCSC")
Plant2<- subset(Plant2, Hostsample2 == "SCSC")
output.df <- Plant2
output.df['com'] = rep(NA, nrow(Plant2)) #Make combination line
output.df$com[which(as.character(output.df$Codesample1) == as.character(output.df$Codesample2))] = 'Within'
output.df$com[which(as.character(output.df$Codesample1) != as.character(output.df$Codesample2))] = 'Between'
output.df['comTrans'] = rep(NA, nrow(Plant2)) #Make combination line
output.df$comTrans[which(as.character(output.df$Transectsample1) == as.character(output.df$Transectsample2))] = 'WithinTrans'
output.df$comTrans[which(as.character(output.df$Transectsample1) != as.character(output.df$Transectsample2))] = 'BetweenTrans'

## For the Within_Block data jump straight to Within sites - within block
## For the All- EachSiteIndividually data jump straight to Within sites - between block
unique(output.df$com)
unique(output.df$comTrans)


## Within sites - within block
output.df1<- subset(output.df, com == "Within")
output.df2a<- subset(output.df1, comTrans == "WithinTrans")
Plantforrun <- output.df2a[,c(2:8,17)]
any(is.na(Plantforrun$Codesample1)) # This should be FALSE so that there are no row missing data
PlantFinal2 = reshape2::melt(Plantforrun, id = c("sample1", "sample2", "Codesample1"))
levels(PlantFinal2$variable)<- c("Competition","Environmental filtering", "Dispersal limitation","Homogenizing dispersal", "Ecological Drift")


# Plot
Plants1<- ggplot(PlantFinal2, aes(fill=variable, y=value, x=Codesample1)) + 
  geom_bar(position="fill", stat="identity")+
  theme(axis.title.x= element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))
Plant2Unweightw <- Plants1 + 
  scale_fill_manual(values = c("light blue","#231F20","dark cyan","#800080","#B5651D")) +
  scale_y_continuous(labels = scales::percent)+
  ylab("Relative Importance")+ theme (strip.text = element_text(size = 11))+
  theme(legend.title = element_blank())+
  scale_x_discrete(limits=c("SP","HTM", "PC","CWR","DMR","RR","KNZ","HSC","KWD","MP"))
Plant2Unweightw


#ggsave(filename="Plots/Plants/iCAMP-Plant_Bin100_2019IndivSitesasPool_Both_Year.pdf", plot=Plant2Unweightw,height= 3.5, width= 6.5, unit= "in")
