
# Description" Comparison of Bray-Curtis similarity between and among blocks

# Load the required package
library(reshape2)
library(reshape)
library(dplyr)
library(ggpubr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(multcompView)
library(emmeans)


# Drawing

#Calling the data set - These are bray-curtis distances obtain from vegan : see PERMANOVA codes 
dist.abund3<- read.csv("Analysis_without_SingletonsAndDoubletons/Plant_Bray_curtis_Matrix_NoSingle_Double.csv")
dist.long= melt(dist.abund3, id= "Name")
dist.long= subset(dist.long, value != "0") # removing the comparisons of the same sample (diagonal) 
dist.long$Bray = 1-dist.long$value # To get similarity


# Now bring the mapping file
Sample1<- read.csv("Plant_Mapping_Sample1_BC_Comp_B.csv")
Sample2<- read.csv("Plant_Mapping_Sample2_BC_Comp.csv")

## Joining dataset
Plant1 <-left_join(dist.long,Sample1, by="Name")
Plant2 <-left_join(Plant1,Sample2, by="variable")
Plant2<- subset(Plant2, Hostsample1 == "SCSC")
Plant2<- subset(Plant2, Hostsample2 == "SCSC")
Plant2b<- subset(Plant2, RemoveExtraSample1 == "Yes")
Plant2b<- subset(Plant2b, RemoveExtraSample2 == "Yes")

# Adding new variables based conditions of comparisons
output.df <- Plant2b

## Make combination line Within vs Between sites comparisonoutput.df['com'] = rep(NA, nrow(Plant2b)) 
output.df$com[which(as.character(output.df$Codesample1) == as.character(output.df$Codesample2))] = 'Within Sites'
output.df$com[which(as.character(output.df$Codesample1) != as.character(output.df$Codesample2))] = 'Between Sites'

## Make combination line Within vs Between blocks comparisons
output.df['comTrans'] = rep(NA, nrow(Plant2b)) 
output.df$comTrans[which(as.character(output.df$Transectsample1) == as.character(output.df$Transectsample2))] = 'Within Block'
output.df$comTrans[which(as.character(output.df$Transectsample1) != as.character(output.df$Transectsample2))] = 'Between Block'

##  Make combination line Within vs Between year comparisons
output.df['comYear'] = rep(NA, nrow(Plant2b)) 
output.df$comYear[which(as.character(output.df$Yearsample1) == as.character(output.df$Yearsample2))] = 'Within Year'
output.df$comYear[which(as.character(output.df$Yearsample1) != as.character(output.df$Yearsample2))] = 'Between Year'
output.df = transform(output.df, comTrans=factor(comTrans,levels = c("Within Block", "Between Block")))
output.df = transform(output.df, com=factor(com,levels = c("Within Sites", "Between Sites")))
output.df = transform(output.df, comYear=factor(comYear,levels = c("Within Year", "Between Year")))


## Within vs Between Block of the same site
output.df2<- subset(output.df, com == "Within Sites") # Considering within sites pairwise comparisons only
output.df3<-subset(output.df2, comYear == "Within Year") # Considering within year only
output.df3<-subset(output.df2, comYear == "Within Year") # Considering within year only
output.df3 <- subset(output.df3,Yearsample1 != 'Y_2018' ) # removing 2018 dataset  as it was only used compaisons between year
output.df3 <- subset(output.df3,Yearsample2 != 'Y_2018' )


# Plot 
## Value are synonym to bray curtis distances
Plants4 <- ggboxplot(output.df3, "Codesample1", "value", ylab="Community dissimilarity",fill = "comTrans")
RCBC1a <- Plants4 + theme(axis.title.x = element_blank())+ 
  scale_fill_manual(values = c("Dark green","#d8cbc4"))+
  scale_x_discrete(limits=c("SP","HTM", "PC","CWR","DMR","RR","KNZ","HSC","KWD","MP"))+
  theme(legend.title=element_blank())
RCBC1a 
ggsave(filename="suppl_Figure2.pdf", plot=RCBC1a,height= 4, width= 6, unit= "in")

# Statistical Analysis

MDcBY = lmer(value ~ comTrans * Categorical_MAPsample1 + (1|Transectsample1), data = output.df3)
summary(MDcBY)
anova(MDcBY) # perform ANOVA


