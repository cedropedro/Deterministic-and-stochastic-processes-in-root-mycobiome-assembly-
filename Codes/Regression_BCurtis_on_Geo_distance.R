#DESCRIPTION: Regression of Bray-Curtis distances on geographic distances in km

library(ggplot2)
library(geosphere)

# Load your data with columns as OTUs/variables, and rows as samples

## Insert the data
df <- read.csv("Mantel_Spatial_2018_2019_NoSingletons.csv")
dF <- df[,41:ncol(df)]

# abundance data frame
abund = df[,41:ncol(df)]

# Subset your samples into smaller dataframe
## environmental vector
Rainfall = df$Rainfall

# longitude and latitude 
geo = data.frame(df$Longitude, df$Latitude)

# Transform them into a distance matrix
## abundance data frame - bray curtis dissimilarity
dist.abund = vegdist(abund, method = "bray", na.rm = TRUE)

## geographic data frame - haversine distance 
d.geo = distm(geo, fun = distHaversine)
dist.geo = as.dist(d.geo)

# Now we can run the mantel command
## abundance vs geographic 
abund_geo  = mantel(dist.abund, dist.geo, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_geo

lm= lm(dist.abund~dist.geo)
summary(lm)

# Now we can do the plots
## First convert the distance matrix (multiple vector) into one vector, and combine matrix into one dataframe mat

aa = as.vector(dist.abund)
gg = as.vector(dist.geo)
c2=as.data.frame(gg)
#write.csv(c2, file= "distance_data.csv")

#new data frame with vectorized distance matrices
mat = data.frame(aa,gg)

#Community composition vs geographic distance

dd = ggplot(mat, aes(y = aa, x = gg/1000)) + 
  geom_point(size = 4, alpha = 0.75, colour = "black",shape = 21) + 
  geom_smooth(method = "lm", colour = "Red", alpha = 0.2) + 
  labs(x = "Geographic distance (km)", y = "Bray-Curtis Dissimilarity", fill = "Difference in rainfall") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"),
         legend.position = "top",
         legend.text = element_text(size = 10, face = "bold"),
         legend.title = element_text(size = 11, face = "bold")) +
  scale_fill_continuous(high = "navy", low = "skyblue")
dd
ggsave(filename="BCvsGeodist-CH2_Dissimilarity_NoSingletons.jpeg", plot=dd,height= 4, width= 6, unit= "in")

