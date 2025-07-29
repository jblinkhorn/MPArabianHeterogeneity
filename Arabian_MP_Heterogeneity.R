library(psych)
library(ggplot2)
library(ggpubr)
library(GGally)
library(ROCR)

####load the data and prep####

levallois_data = read.csv("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/Levallois_flake_data_2025_02_12.csv") #jb
levallois_data$SpCode = 0 
levallois_data[which(levallois_data$Site == "Tor Faraj"),"SpCode"] <- 1 #Neanderthal
levallois_data[which(levallois_data$Site == "Kebara X"),"SpCode"] <- 1 #Neanderthal
levallois_data[which(levallois_data$Site == "Ksar Akil 28A"),"SpCode"] <- 1 #Neanderthal
levallois_data[which(levallois_data$Site == "JKF-1"),"SpCode"] <- 2 #unknown
levallois_data[which(levallois_data$Site == "ALM-3"),"SpCode"] <- 2 #unknown

levallois_data[which(levallois_data$SpCode == 0),"Sp"] <- "sapiens"
levallois_data[which(levallois_data$SpCode == 1),"Sp"] <- "neanderthal"
levallois_data[which(levallois_data$SpCode == 2),"Sp"] <- "unknown"

levallois_nona = levallois_data[-which(apply(levallois_data, 1, function(x)any(is.na(x)))),]
var_cols = 3:10 #set var_cols for PCA as raw measurements
colnames(levallois_nona)[which(names(levallois_nona) == "DorsalSA_TwoTrapezoids")] <- "FlakeArea" #

# Reorder 'Site' based on 'Region' for plotting
levallois_nona$Site <- factor(levallois_nona$Site, levels = unique(levallois_nona$Site[order(levallois_nona$Region)]))

####Figure 3 ####
# Updated ggplot code
Elong <- ggplot(data = levallois_nona, mapping = aes(Site, Elongation, group = Site)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 0.5) +
  geom_boxplot(aes(colour = Region, fill = Region), alpha = 0.8, outlier.shape = NA) +
  theme_minimal() +
  theme(text = element_text(family = "sans", size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5),
        axis.title.x = element_blank())

Conv <- ggplot(data = levallois_nona, mapping = aes(Site, Convergence, group = Site)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 0.5) +
  geom_boxplot(aes(colour = Region, fill = Region), alpha = 0.8, outlier.shape = NA) +
  theme_minimal() +
  theme(text = element_text(family = "sans", size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5),
        axis.title.x = element_blank())

Platform <- ggplot(data = levallois_nona, mapping = aes(Site, PlatformArea, group = Site)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 0.5) +
  geom_boxplot(aes(colour = Region, fill = Region), alpha = 0.8, outlier.shape = NA) +
  theme_minimal() +
  theme(text = element_text(family = "sans", size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5),
        axis.title.x = element_blank())

Flake <- ggplot(data = levallois_nona, mapping = aes(Site, FlakeArea, group = Site)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 0.5) +
  geom_boxplot(aes(colour = Region, fill = Region), alpha = 0.8, outlier.shape = NA) +
  theme_minimal() +
  theme(text = element_text(family = "sans", size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5),
        axis.title.x = element_blank())

# Load necessary packages
library(ggplot2)
library(dplyr)

# Calculate percentages within each Site
dorsal_summary <- levallois_nona %>%
  group_by(Region, Site, Dorsal.scar.pattern) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = count / sum(count))  # Normalize within Site

# Combine Region and Site into a single ordered factor
dorsal_summary$Region_Site <- factor(paste(dorsal_summary$Region, dorsal_summary$Site, sep = " - "),
                                     levels = unique(paste(dorsal_summary$Region, dorsal_summary$Site, sep = " - ")))
library(viridis)
# Create bar chart with properly spaced and ordered sites
Scar <- ggplot(dorsal_summary, aes(x = Region_Site, y = percentage, fill = Dorsal.scar.pattern)) +
  geom_bar(stat = "identity", position = "fill") +  # Ensures each bar totals 100%
  scale_fill_viridis(discrete = TRUE) +
  scale_y_continuous(labels = scales::percent_format()) +  # Format as percentages
  scale_x_discrete(labels = unique(dorsal_summary$Site)) +  # Display only 'Site' names
  theme_minimal() +
  theme(text = element_text(family = "sans", size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank()) +
  labs(y = "Dorsal Scar %", fill = "Dorsal Scar Pattern")

rawmat_summary <- levallois_nona %>%
  group_by(Region, Site, raw_mat_simplified) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = count / sum(count))  # Normalize within Site

# Combine Region and Site into a single ordered factor
rawmat_summary$Region_Site <- factor(paste(rawmat_summary$Region, rawmat_summary$Site, sep = " - "),
                                     levels = unique(paste(rawmat_summary$Region, rawmat_summary$Site, sep = " - ")))
library(viridis)
# Create bar chart with properly spaced and ordered sites
RawMat <- ggplot(rawmat_summary, aes(x = Region_Site, y = percentage, fill = raw_mat_simplified)) +
  geom_bar(stat = "identity", position = "fill") +  # Ensures each bar totals 100%
  scale_fill_viridis(discrete = TRUE, option="turbo") +
  scale_x_discrete(labels = unique(rawmat_summary$Site)) + 
  scale_y_continuous(labels = scales::percent_format()) +  # Format as percentages
  theme_minimal() +
  theme(text = element_text(family = "sans", size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank()) +
  labs(y = "Raw Material %", fill = "Raw Material")


library(ggplot2)
library(cowplot)
library(gridExtra)

# Extract legends for Flake, RawMat, and Scar
gtable_Flake <- ggplotGrob(Flake + theme(legend.position = "right"))
legend_Flake <- gtable_Flake$grobs[[which(sapply(gtable_Flake$grobs, function(x) x$name) == "guide-box")]]

gtable_RawMat <- ggplotGrob(RawMat + theme(legend.position = "bottom", legend.title = element_blank())+
                              guides(fill = guide_legend(ncol = 2)))
legend_RawMat <- gtable_RawMat$grobs[[which(sapply(gtable_RawMat$grobs, function(x) x$name) == "guide-box")]]

gtable_Scar <- ggplotGrob(Scar + theme(legend.position = "bottom", legend.title = element_blank())+
                            guides(fill = guide_legend(ncol = 2)))
legend_Scar <- gtable_Scar$grobs[[which(sapply(gtable_Scar$grobs, function(x) x$name) == "guide-box")]]

# Remove legends from individual plots
Flake <- Flake + theme(legend.position = "none", axis.text.x = element_blank())
Elong <- Elong + theme(legend.position = "none", axis.text.x = element_blank())
Platform <- Platform + theme(legend.position = "none", axis.text.x = element_blank())
Conv <- Conv + theme(legend.position = "none", axis.text.x = element_blank())
RawMat <- RawMat + theme(legend.position = "none")
Scar <- Scar + theme(legend.position = "none")

# Arrange main plots in a grid
main_plots <- plot_grid(Flake, Platform, Elong, Conv, RawMat, Scar, ncol = 2, align="v")

legend_right <- plot_grid(NULL, legend_Flake, ncol = 2, rel_widths = c(0, 0.3))

legend_bottom <- plot_grid(legend_RawMat, legend_Scar, NULL, ncol = 3, rel_widths = c(0.4, 0.4, 0.1))  # Centers the Flake legend


# First, arrange the main plots with the right-hand legends
main_right <- plot_grid(
  main_plots, legend_right, 
  ncol = 2, 
  rel_widths = c(6, 1) # Ensures right-side legends don’t encroach on bottom space
)

# Now, add the bottom legend separately beneath everything
final_plot <- plot_grid(
  main_right,
  legend_bottom,
  ncol = 1,
  rel_heights = c(1, 0.1)  # Prevents right-hand legends from overlapping below
)

# Display final combined plot
final_plot

ggsave("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/final_plot.pdf", plot = final_plot, width = 12, height = 12, units = "in")

ggsave("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/final_plot.jpg", plot = final_plot, width = 12, height = 12, units = "in")


####PCA analysis####

## pre pca
KMO(levallois_nona[,var_cols])

cortest.bartlett(cor(levallois_nona[,var_cols]), n=nrow(levallois_nona[,var_cols]))

## separate out training data
set.seed(909)
scaled_vars = apply(levallois_nona[, var_cols], 2, scale)
levallois_nona[, var_cols] <- scaled_vars
unknown_sp <- which(levallois_nona$SpCode==2)

## run pca
pca_levallois_train = prcomp(levallois_nona[-unknown_sp,var_cols], retx = T, scale = F, center = F)

##rotate JKF-1 and ALM-1 or orientation of comparative dataset
unknown_proj = as.matrix(levallois_nona[unknown_sp,var_cols]) %*% pca_levallois_train$rotation
unknown_proj = cbind(levallois_nona[unknown_sp,], unknown_proj)
levallois_proj_train = cbind(levallois_nona[-unknown_sp,], pca_levallois_train$x)
levallois_proj = rbind(levallois_proj_train, unknown_proj)

pca_levallois_variance = summary(pca_levallois_train)
pc_var_data = t(pca_levallois_variance$importance[c(2,3),])
pc_var_data = cbind(1:nrow(pc_var_data), round(pc_var_data, 2) * 100)
pc_var_data <- as.data.frame(pc_var_data)
names(pc_var_data) <- c("PC", "Variance", "Cumulative")

ggplot(data = pc_var_data) +
  geom_path(mapping = aes(x = PC, y = Variance), colour="darkgrey") +
  geom_point(mapping = aes(x = PC, y = Variance)) +
  coord_cartesian(ylim=c(0,60),
                  xlim=c(1,8),
                  clip="off") +
  scale_x_continuous(breaks=c(1:8)) +
  labs(title="Principle Components", y = "% Variance") +
  annotate("text",
           x = pc_var_data$PC,
           y = -5,
           label = pc_var_data$Cumulative,
           size = 3) +
  annotate("text",
           x = 0.25,
           y = -5,
           label = "% Variance",
           size = 3) +
  theme_minimal() +
  theme(text = element_text(family="Times", size=12),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 10, unit = "mm"),
        plot.title = element_text(face="bold",hjust=0.5,size=15),
        axis.title.x.bottom = element_text(margin=margin(t=5,r=1,b=1,l=1, unit ="mm")),
        axis.text.x = element_text(margin=margin(t=8,r=1,b=1,l=1, unit ="mm")))

nsites = length(unique(levallois_proj$Site))
unknown_sp <- which(levallois_proj$SpCode==2)

levallois_proj$Site <- factor(levallois_proj$Site, 
                              levels = unique(levallois_proj$Site[order(levallois_proj$Region, levallois_proj$Site)]))

####Figure 4####

p1 <- ggplot(data = levallois_proj, 
             mapping = aes(Site, PC1, group = Site)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 0.5) +
  geom_boxplot(aes(colour = Region, fill = Region), alpha = 0.8, outlier.shape = NA) +
  theme_minimal() +
  theme(text = element_text(family="Times", size=12),
        plot.title = element_text(face="bold", hjust=0.5, size=15),
        axis.text.x = element_text(size=10, angle=45, hjust=1),  # Optional: Rotate labels for readability
        axis.title.x = element_blank())


p2 <- ggplot(data = levallois_proj, 
             mapping = aes(Site, PC2, group = Site)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 0.5) +
  geom_boxplot(aes(colour = Region, fill = Region), alpha = 0.8, outlier.shape = NA) +
  theme_minimal() +
  theme(text = element_text(family="Times", size=12),
        plot.title = element_text(face="bold", hjust=0.5, size=15),
        axis.text.x = element_text(size=10, angle=45, hjust=1),  # Optional: Rotate labels for readability
        axis.title.x = element_blank())

p3 <- ggplot(data = levallois_proj, 
             mapping = aes(Site, PC3, group = Site)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 0.5) +
  geom_boxplot(aes(colour = Region, fill = Region), alpha = 0.8, outlier.shape = NA) +
  theme_minimal() +
  theme(text = element_text(family="Times", size=12),
        plot.title = element_text(face="bold", hjust=0.5, size=15),
        axis.text.x = element_text(size=10, angle=45, hjust=1),  # Optional: Rotate labels for readability
        axis.title.x = element_blank())

p4 <- ggplot(data = levallois_proj, 
             mapping = aes(Site, PC4, group = Site)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 0.5) +
  geom_boxplot(aes(colour = Region, fill = Region), alpha = 0.8, outlier.shape = NA) +
  theme_minimal() +
  theme(text = element_text(family="sans", size=12),
        plot.title = element_text(face="bold", hjust=0.5, size=15),
        axis.text.x = element_text(size=10, angle=45, hjust=1),  # Optional: Rotate labels for readability
        axis.title.x = element_blank())

library(ggplot2)
library(cowplot)

# Ensure no legends and no x-axis labels (except p4)
p1 <- p1 + theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank())
p2 <- p2 + theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank())
p3 <- p3 + theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank())
p4 <- p4 + theme(legend.position = "none")  # Keep x-axis labels for p4

# Convert p4 to a gtable object and extract the legend
p4_gtable <- ggplotGrob(p4 + theme(legend.position = "bottom"))
legend <- p4_gtable$grobs[[which(sapply(p4_gtable$grobs, function(x) x$name) == "guide-box")]]

# Stack plots in a single column
plots_grid <- plot_grid(p1, p2, p3, p4, ncol = 1, align = "v")

# Add legend below the stacked plots
final_plot <- plot_grid(plots_grid, legend, ncol = 1, rel_heights = c(1, 0.1))

final_plot

ggsave("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/final_plot2.pdf", plot = final_plot, width = 8, height = 12, units = "in")
ggsave("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/final_plot2.jpg", plot = final_plot, width = 8, height = 12, units = "in")

####Figure 8 ####

library(psych)
library(ggplot2)
library(ggpubr)
library(GGally)
library(ROCR)

p5 <- ggplot(data = levallois_proj, 
             mapping = aes(Site, PC2, group = Site)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 0.5) +
  geom_boxplot(aes(colour = as.factor(Type), fill = as.factor(Type)), alpha = 0.5, outlier.shape = NA) +
  scale_fill_manual(values = c("purple", "green"), 
                    name = "Site Type", 
                    labels = c("Open Air", "Cave or Rockshelter")) +  # Fill colors for boxplots
  scale_color_manual(values = c("black", "black"), 
                     name = "Site Type", 
                     labels = c("Open Air", "Cave or Rockshelter")) +  # Outline colors for boxplots
  theme_minimal() +
  theme(text = element_text(family="Times", size=12),
        plot.title = element_text(face="bold", hjust=0.5, size=15),
        axis.text.x = element_text(size=10, angle=45, hjust=1),  # Optional: Rotate labels for readability
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

p5a <- ggplot() +
  # Original density plot for SpCode != 2
  geom_density(data = subset(levallois_proj, SpCode != 2),  
               aes(x = PC2, fill = as.factor(Type)), alpha = 0.6) +
  
  # Overlay density plot for SpCode == 2
  geom_density(data = subset(levallois_proj, SpCode == 2),  
               aes(x = PC2, color = Site, linetype = Site), alpha = 0.1, size = 1) +
  
  # Set manual colors & rename legend labels
  scale_fill_manual(values = c("purple", "green"), 
                    name = "Site Type", 
                    labels = c("Open Air", "Cave or Rockshelter")) +
  
  scale_color_manual(values = c("Black", "black")) +
  
  scale_linetype_manual(values = c("dotted", "longdash")) +
  
  theme_minimal() +
  theme(legend.position = "bottom")

p6 <- ggplot(data = levallois_proj, 
             mapping = aes(Site, PC2, group = Site)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 0.5) +
  geom_boxplot(aes(colour = Sp, fill = Sp), alpha = 0.5, outlier.shape = NA) +
  scale_fill_manual(values = c("red", "blue", "grey50"), 
                    name = "Species", 
                    labels = c("Neanderthal", "Homo sapiens", "Unknown")) +  # Fill colors for boxplots
  scale_color_manual(values = c("black", "black", "black"), 
                     name = "Species", 
                     labels = c("Neanderthal", "Homo sapiens", "Unknown")) +  # Outline colors for boxplots
  theme_minimal() +
  theme(text = element_text(family="Times", size=12),
        plot.title = element_text(face="bold", hjust=0.5, size=15),
        axis.text.x = element_text(size=10, angle=45, hjust=1),  # Optional: Rotate labels for readability
        axis.title.x = element_blank())+
  theme(legend.position = "none")

p6a <- ggplot() +
  # Original density plot for SpCode != 2
  geom_density(data = subset(levallois_proj, SpCode != 2),  
               aes(x = PC2, fill = Sp), alpha = 0.6) +
  
  # Overlay density plot for SpCode == 2
  geom_density(data = subset(levallois_proj, SpCode == 2),  
               aes(x = PC2, color = Site, linetype = Site), alpha = 0.1, size = 1) +
  
  # Set manual colors & rename legend labels
  scale_fill_manual(values = c("red", "blue"), 
                    name = "Species", 
                    labels = c("Neanderthal", "Homo sapiens")) +
  
  scale_color_manual(values = c("Black", "black")) +
  
  scale_linetype_manual(values = c("dotted", "longdash")) +
  
  theme_minimal() +
  theme(legend.position = "bottom")

library(patchwork)

(p5 + p6) / (p5a + p6a)

#need to add a line to save this

#####UPDATED OLD FIGURE of JKF1 and ALM3 ON PC2 and PC£3####


#JKF1_hull_index <- chull(subset(levallois_proj,Site=="JKF-1")[,c("PC3","PC2")])
#JKF1_hull_coords <- subset(levallois_proj,Site=="JKF-1")[JKF1_hull_index,c("PC3","PC2")]
#JKF1_hull_centroid <- colMeans(JKF1_hull_coords)
#JKF1_hull_centroid <- data.frame(PC3=JKF1_hull_centroid[1],PC2=JKF1_hull_centroid[2])

#ALM3_hull_index <- chull(subset(levallois_proj,Site=="ALM-3")[,c("PC3","PC2")])
#ALM3_hull_coords <- subset(levallois_proj,Site=="ALM-3")[ALM3_hull_index,c("PC3","PC2")]
#ALM3_hull_centroid <- colMeans(ALM3_hull_coords)
#ALM3_hull_centroid <- data.frame(PC3=ALM3_hull_centroid[1],PC2=ALM3_hull_centroid[2])

#ggplot(data = subset(levallois_proj, SpCode != 2), aes(x = PC3, y = PC2)) +
#  stat_density_2d(aes(fill=Sp, alpha=stat(level)),
#                  geom="polygon") +
#  geom_polygon(data=subset(levallois_proj,Site=="JKF-1")[JKF1_hull_index,c("PC3","PC2")],
#               aes(x = PC3, y = PC2),
#               fill=NA,
#               colour="black") +
#  geom_point(data=JKF1_hull_centroid,
#             aes(x = PC3, y = PC2),
#             shape=21,
#             colour="black") +
#  geom_polygon(data=subset(levallois_proj,Site=="ALM-3")[ALM3_hull_index,c("PC2","PC3")],
#               aes(x = PC3, y = PC2),
#               fill=NA,
#               colour="blue") +
#  geom_point(data=ALM3_hull_centroid,
#             aes(x = PC3, y = PC2),
#             shape=21,
#             colour="blue") +
#  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
#  scale_color_manual(values=c("#F8766D", "#00BA38")) +
#  scale_x_continuous(breaks=c(-5:10)) +
#  scale_y_continuous(breaks=c(-4:3)) +
#  labs(title="PCA Scores") +
#  theme_minimal() +
#  theme(text = element_text(family="Times", size=12),
#        plot.title = element_text(face="bold",hjust=0.5,size=15),
#        legend.position = "top")
#
#ggsave(filename="./Images/pca_23_JKF1_ALM3.pdf",
#       device = "pdf")

#####Generating variables for analysis####
library(pastclim)

##extracting pastclim data##

train <- levallois_proj_train
full <- levallois_proj
train_sites <- unique(full[, c("Site", "N", "E", "Min", "Mid", "Max")])
train_sites$lat <- train_sites$N
train_sites$lon <- train_sites$E
train_sites$reps <- table(as.character(full$Site))

geog_list <- list()
time_series_list <- list()
for (i in 1:nrow(train_sites)){#extract time series climatic data for BIO01, BIO04, BIO12, BIO15 for each site location spanning timeframe of occupation
  
  time_series_list[[i]] <- pastclim::location_series(x = data.frame(train_sites[i,]), 
                                                     time_bp = seq(train_sites$Min[[i]]*-1000, train_sites$Max[[i]]*-1000, by=-1000),
                                                     dataset = "Krapp2021", 
                                                     bio_variables = c("bio01", "bio04", "bio12", "bio15"))
}

for (i in 1:nrow(train_sites)){#extract altitude and rugosity data for each site
  geog_list[[i]] <- pastclim::location_series(x = data.frame(train_sites[i,]), 
                                              time_bp = 0,
                                              dataset = "Krapp2021", 
                                              bio_variables = c("altitude", "rugosity"))}

geog <- dplyr::bind_rows(geog_list)
geog$name <- as.character(train_sites$Site)
geog2 <- geog[order(geog$name),]

####Lithics####

PC1_dist <- dist(train$PC1, method = "euclidean")
PC2_dist <- dist(train$PC2, method = "euclidean")
PC3_dist <- dist(train$PC3, method = "euclidean")
PC4_dist <- dist(train$PC4, method = "euclidean")

#####Cost-Paths#####

library(magrittr)
library(raster)
library(dendextend)
library(NbClust)
library(sf)

msalsa <- SpatialPointsDataFrame(train_sites[,c("E", "N")], train_sites, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs  ") )
msalsa <- spTransform(msalsa, CRSobj = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

caf <- data.frame(msalsa@coords)
names(caf) <- c("E", "N")
points_proj <- extent(round(min(caf$E)-200000), round(max(caf$E)+200000), round(min(caf$N)-200000), round(max(caf$N)+200000))
points <- extent(round(min(train_sites$E)-2), round(max(train_sites$E)+2), round(min(train_sites$N)-2), round(max(train_sites$N)+2))

srtm <- terra::rast("C:/Users/jblink/OneDrive - The University of Liverpool/CB MSA_LSA/SRTM_1km.tif")
srtm_2 <- crop(srtm, points)
this_crs <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "
srtm_3 <- srtm_2 %>% terra::project(this_crs)
srtm_3 <- raster(srtm_3)

library(gdistance)
altDiff <- function(x){x[2] - x[1]}
hd <- transition(srtm_3, altDiff, 8, symm=FALSE) #nb based on queen movement; nb gives negative values
slope <- geoCorrection(hd)
adj <- adjacent(srtm_3, cells=1:ncell(srtm_3), pairs=TRUE, directions=8) #nb based on queen movement
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05)) # based on speed of movement from toblers function
Conductance <- geoCorrection(speed)
xc <- costDistance(Conductance, msalsa) # values are seconds (?) using toblers function for m/s
costpath_dist <- as.dist(xc)

####permutations of age-dependent variables####

climate_combos <- vector("list", 1000)
for (i in 1:1000) {random_rows <- lapply(time_series_list, function(df) {df[sample(nrow(df), 1), ]  # Select a random row from each dataframe
  })
  # Combine the selected random rows into a single dataframe
  climate_combos[[i]] <- do.call(rbind, random_rows)
}

####Mantel Tests and Multiple Matrix Regressions####

##Our analysis was run using a distributed computing system; the code here does run, but will take a long time####
#the distributed approach split this into 1000 discrete indexed analyses; as each analysis is on a slightly different dataset, there is no interdependence#

scale_fun <- function(x){(x-min(x))/(max(x)-min(x))} #scales to 0:1
PC_list <- list(PC1_dist, PC2_dist, PC3_dist, PC4_dist)

start_date <- Sys.Date()

single_k_list <- list()
multiple_k_list <- list()
for(k in c(1:1000)){
  site_level <- data.frame(matrix(ncol = 0, nrow = nrow(geog)))
  site_level$altitude <- geog$altitude
  site_level$rugosity <- geog$rugosity
  site_level$age <- climate_combos[[k]]$time_bp
  site_level$bio01 <- climate_combos[[k]]$bio01
  site_level$bio04 <- climate_combos[[k]]$bio04
  site_level$bio12 <- climate_combos[[k]]$bio12
  site_level$bio15 <- climate_combos[[k]]$bio15
  site_level <- site_level[rep(seq_len(nrow(site_level)), times = train_sites$reps), ]
  
  repeated_rows <- xc[rep(seq_len(nrow(xc)), times = train_sites$reps), ]
  repeated_matrix <- repeated_rows[, rep(seq_len(ncol(xc)), times = train_sites$reps)]
  
  dist_list <- list(costpath = dist(repeated_matrix),
                    altitude = dist(site_level$altitude),
                    rugosity = dist(site_level$rugosity), 
                    age = dist(site_level$age),
                    bio01 = dist(site_level$bio01),
                    bio04 = dist(site_level$bio04),
                    bio12 = dist(site_level$bio12),
                    bio15 = dist(site_level$bio15),
                    species = dist(train$SpCode, method = "euclidean") ,
                    rawmat = dist(train[,13:16], method="binary"),
                    type = dist(train$Type, method="binary"),
                    prop = dist(train$Prop_RT_Core),
                    count = dist(train$Count_RT_Core))
  
  dist_list2 <- lapply(dist_list, function(x) scale_fun(x))
  PC_list2 <- lapply(PC_list, function(x) scale_fun(x))
  
  simple_mantel_list <- list()
  multi_mantel_list <- list()
  for(i in 1:4){
    mantel_out <- lapply(dist_list2, function(x) vegan::mantel(PC_list2[[i]], x, method="pearson", permutations = 100))
    mantelresults2 <- as.data.frame(lapply(mantel_out, function(x) c(statistic = x$statistic, pvalue = x$signif)))
    mantelresults2[3,] <- p.adjust(t(mantelresults2[2,]), method = "BH")                     
    simple_mantel_list[[i]] <- t(mantelresults2)
    
    xs <- mantelresults2[2,]<0.05 #up to here
    #all_dists_sig <- dist_list2
    all_dists_sig <- dist_list2[xs==T]
    multi1 <- phytools::multi.mantel(PC_list2[[i]], all_dists_sig, nperm=100)
    names(multi1$coefficients) <- c("intercept", names(all_dists_sig))
    multi_mantel_list[[i]] <- capture.output(multi1)
    
  }
  single_k_list[[k]] <- simple_mantel_list
  multiple_k_list[[k]] <- multi_mantel_list
  
  sheet_names <- c("PC1", "PC2", "PC3", "PC4")
  
  for(i in 1:length(simple_mantel_list)){xlsx::write.xlsx(simple_mantel_list[[i]], paste0("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/Simple Mantel Tests_", start_date, "_", k, ".xlsx"), sheetName = sheet_names[i], append = T)}
  for(i in 1:length(multi_mantel_list)){xlsx::write.xlsx(multi_mantel_list[[i]], paste0("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/Multiple Mantel Tests_", start_date, "_",k, ".xlsx"), sheetName = sheet_names[i], append = T)}
  
  print(k)
}

#####Sample Code for distributed analysis####

scale_fun <- function(x){(x-min(x))/(max(x)-min(x))} #scales to 0:1
PC_list <- list(PC1_dist, PC2_dist, PC3_dist, PC4_dist)

site_level_list <- list()

for(k in c(1:1000)){
site_level <- data.frame(matrix(ncol = 0, nrow = nrow(geog)))
site_level$altitude <- geog$altitude
site_level$rugosity <- geog$rugosity
site_level$age <- climate_combos[[k]]$time_bp
site_level$bio01 <- climate_combos[[k]]$bio01
site_level$bio04 <- climate_combos[[k]]$bio04
site_level$bio12 <- climate_combos[[k]]$bio12
site_level$bio15 <- climate_combos[[k]]$bio15
site_level <- site_level[rep(seq_len(nrow(site_level)), times = train_sites$reps), ]
repeated_rows <- xc[rep(seq_len(nrow(xc)), times = train_sites$reps), ]
repeated_matrix <- repeated_rows[, rep(seq_len(ncol(xc)), times = train_sites$reps)]
site_level$costpath <- repeated_matrix
site_level <- lapply(site_level, function(x) scale_fun(x))
site_level_list[[k]] <- site_level
}

PC_list2 <- lapply(PC_list, function(x) scale_fun(x))

setwd("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/distance matrices/indexed")

#export each stable distance matrix

as.matrix(dist(site_level$altitude)) %>% write.csv("altitude.csv")
as.matrix(site_level$costpath) %>% write.csv("costpath.csv")
as.matrix(dist(site_level$rugosity)) %>% write.csv("rugosity.csv")
as.matrix(dist(train$SpCode, method="binary")) %>% write.csv("species.csv")
as.matrix(dist(train[,13:16], method="binary")) %>% write.csv("rawmat.csv")
as.matrix(dist(train$Type, method="binary")) %>% write.csv("type.csv")
as.matrix(dist(scale_fun(train$Prop_RT_Core))) %>% write.csv("prop.csv")
as.matrix(dist(scale_fun(train$Count_RT_Core))) %>% write.csv("count.csv")

as.matrix(dist(PC_list2[[1]])) %>% write.csv("PC1.csv")
as.matrix(dist(PC_list2[[2]])) %>% write.csv("PC2.csv")
as.matrix(dist(PC_list2[[3]])) %>% write.csv("PC3.csv")
as.matrix(dist(PC_list2[[4]])) %>% write.csv("PC4.csv")

#export each indexed distance matrix

for(k in 1:1000){
  as.matrix(dist(site_level_list[[k]]$age)) %>% write.csv(paste0("age", k, ".csv", sep=""))
  as.matrix(dist(site_level_list[[k]]$bio01)) %>% write.csv(paste0("MAT", k, ".csv", sep=""))
  as.matrix(dist(site_level_list[[k]]$bio04)) %>% write.csv(paste0("MATV", k, ".csv", sep=""))
  as.matrix(dist(site_level_list[[k]]$bio12)) %>% write.csv(paste0("MAP", k, ".csv", sep=""))
  as.matrix(dist(site_level_list[[k]]$bio15)) %>% write.csv(paste0("MAPV", k, ".csv", sep=""))
  
}

#create output files#

for (i in 0:999) {
  file_name <- paste0("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/distance matrices/indexed/PC1_Simple", i, ".csv")
  file.create(file_name)}

for (i in 0:999) {
  file_name <- paste0("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/distance matrices/indexed/PC1_Multiple", i, ".csv")
  file.create(file_name)}

for (i in 0:999) {
  file_name <- paste0("C:/Users/jblink/OneDrive - The University of Liverpool/NITD/distance matrices/indexed/PC2_Simple", i, ".csv")
  file.create(file_name)}

for (i in 0:999) {
  file_name <- paste0("C:/Users/jblink/OneDrive - The University of Liverpool/NITD/distance matrices/indexed/PC2_Multiple", i, ".csv")
  file.create(file_name)}

for (i in 0:999) {
  file_name <- paste0("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/distance matrices/indexed/PC3_Simple", i, ".csv")
  file.create(file_name)}

for (i in 0:999) {
  file_name <- paste0("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/distance matrices/indexed/PC3_Multiple", i, ".csv")
  file.create(file_name)}

for (i in 0:999) {
  file_name <- paste0("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/distance matrices/indexed/PC4_Simple", i, ".csv")
  file.create(file_name)}

for (i in 0:999) {
  file_name <- paste0("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/distance matrices/indexed/PC4_Multiple", i, ".csv")
  file.create(file_name)}

#UPLOAD DATA TO DISTRIBUTED COMPUTING SYSTEM#

#Analysis set up for distributed computing#

dist_list <- list(costpath = read.csv("costpath.csv", row.names=1),
                  altitude = read.csv("altitude.csv", row.names=1),
                  rugosity = read.csv("rugosity.csv", row.names=1), 
                  age = read.csv("age.csv", row.names=1),
                  bio01 = read.csv("MAT.csv", row.names=1),
                  bio04 = read.csv("MATV.csv", row.names=1),
                  bio12 = read.csv("MAP.csv", row.names=1),
                  bio15 = read.csv("MAPV.csv", row.names=1),
                  species = read.csv("species.csv", row.names=1),
                  rawmat = read.csv("rawmat.csv", row.names=1),
                  type = read.csv("type.csv", row.names=1),
                  prop = read.csv("prop.csv", row.names=1),
                  count = read.csv("count.csv", row.names=1))

PC1 <- read.csv("PC1.csv", row.names=1) #edit for each PC

mantel_out <- lapply(dist_list, function(x) vegan::mantel(PC1, x, method="pearson", permutations = 1000))
mantelresults2 <- as.data.frame(lapply(mantel_out, function(x) c(statistic = x$statistic, pvalue = x$signif)))

xs <- mantelresults2[2,]<0.05
all_dists_sig <- dist_list[xs==T]
multi1 <- phytools::multi.mantel(PC1, all_dists_sig, nperm=1000)
names(multi1$coefficients) <- c("intercept", names(all_dists_sig))
write.csv(t(mantelresults2), "Simple_PC1.csv")#edit for each PC
write.csv(capture.output(multi1), "Multi_PC1.csv")#edit for each PC

#####Figure if run on for loop with XLSX outputs; need to edit locations etc####

library(readxl)
#d <- read_excel("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/MM/Multiple Mantel Tests_2025-02-12_1.xlsx")

library(ggplot2)
library(dplyr)
library(cowplot)

# Sample data preparation
patterns <- c("age", "altitude", "bio01", "bio04", "bio12", "bio15", "costpath", "count", "prop", "rawmat", "rugosity", "species", "type")
all_ids <- c(patterns, "R-Squared")
plots <- list()
plot_titles <- c("Box Plot of Model Estimates for each variable and R-Squared values for PC1",
                 "Box Plot of Model Estimates for each variable and R-Squared values for PC2",
                 "Box Plot of Model Estimates for each variable and R-Squared values for PC3",
                 "Box Plot of Model Estimates for each variable and R-Squared values for PC4")


plot_titles <- c("PC1", "PC2", "PC3", "PC4")

total_runs <- 500


data_combo_list <- list()
model_combo_list <- list()

for(j in 1:4) {
#  j <- 2
  data_combo <- data.frame()
  model_combo <- data.frame()
# j <- 2
  for(k in 1:100) {
    #k <- 1
    #a <- as.data.frame(multiple_k_list[[k]][[j]]) #from completed analysis
    d <- read_excel(paste0("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/MM2/Multiple Mantel Tests_2025-03-12_", k, ".xlsx"), sheet=paste0("PC", j)) #from saved excel
    a <- data.frame(d[[2]]) #from saved excel
    colnames(a) <- "data"
    
    filtered_data <- data.frame()
    for (pattern in patterns) {
      matches <- a[grepl(pattern, a$data, ignore.case = TRUE), ]
      filtered_data <- rbind(filtered_data, matches)
    }
    filtered_data$data <- as.character(filtered_data[[1]])
    
    data_string <- paste(filtered_data$data, collapse = "\n")
    split_data <- read.table(text = data_string, fill = TRUE, header = FALSE, stringsAsFactors = FALSE)
    split_data_df <- data.frame(split_data, stringsAsFactors = FALSE)
    split_data_df <- split_data_df[,1:4]
    split_data_df$Model <- k
    colnames(split_data_df) <- c("ID", "Value1", "Value2", "Significance", "Model")
    
    data_combo <- rbind(data_combo, split_data_df)
    data_combo <- na.omit(data_combo)
    
    
    stats <- c("R-squared", "F-statistic")
    filtered_data2 <- data.frame()
    for (pattern in stats) {
      matches <- a[grepl(pattern, a$data, ignore.case = TRUE), ]
      filtered_data2 <- rbind(filtered_data2, matches)
    }
    filtered_data2$data <- as.character(filtered_data2[[1]])
    
    # Corrected regular expression extraction
    model_stats <- c(k, as.numeric(unlist(regmatches(filtered_data2[[1]], gregexpr("[[:digit:]]+\\.?[[:digit:]]*", filtered_data2[[1]])))))
    model_combo <- rbind(model_combo, model_stats)
  }
  
  for(k in 101:300) {
    #k <- 1
    #a <- as.data.frame(multiple_k_list[[k]][[j]]) #from completed analysis
    d <- read_excel(paste0("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/MM2/Multiple Mantel Tests_2025-05-06_", k, ".xlsx"), sheet=paste0("PC", j)) #from saved excel
    a <- data.frame(d[[2]]) #from saved excel
    colnames(a) <- "data"
    
    filtered_data <- data.frame()
    for (pattern in patterns) {
      matches <- a[grepl(pattern, a$data, ignore.case = TRUE), ]
      filtered_data <- rbind(filtered_data, matches)
    }
    filtered_data$data <- as.character(filtered_data[[1]])
    
    data_string <- paste(filtered_data$data, collapse = "\n")
    split_data <- read.table(text = data_string, fill = TRUE, header = FALSE, stringsAsFactors = FALSE)
    split_data_df <- data.frame(split_data, stringsAsFactors = FALSE)
    split_data_df <- split_data_df[,1:4]
    split_data_df$Model <- k
    colnames(split_data_df) <- c("ID", "Value1", "Value2", "Significance", "Model")
    
    data_combo <- rbind(data_combo, split_data_df)
    data_combo <- na.omit(data_combo)
    
    
    stats <- c("R-squared", "F-statistic")
    filtered_data2 <- data.frame()
    for (pattern in stats) {
      matches <- a[grepl(pattern, a$data, ignore.case = TRUE), ]
      filtered_data2 <- rbind(filtered_data2, matches)
    }
    filtered_data2$data <- as.character(filtered_data2[[1]])
    
    # Corrected regular expression extraction
    model_stats <- c(k, as.numeric(unlist(regmatches(filtered_data2[[1]], gregexpr("[[:digit:]]+\\.?[[:digit:]]*", filtered_data2[[1]])))))
    model_combo <- rbind(model_combo, model_stats)
  }
  
  for(k in 301:500) {
    #k <- 1
    #a <- as.data.frame(multiple_k_list[[k]][[j]]) #from completed analysis
    d <- read_excel(paste0("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/MM2/Multiple Mantel Tests_2025-05-08_", k, ".xlsx"), sheet=paste0("PC", j)) #from saved excel
    a <- data.frame(d[[2]]) #from saved excel
    colnames(a) <- "data"
    
    filtered_data <- data.frame()
    for (pattern in patterns) {
      matches <- a[grepl(pattern, a$data, ignore.case = TRUE), ]
      filtered_data <- rbind(filtered_data, matches)
    }
    filtered_data$data <- as.character(filtered_data[[1]])
    
    data_string <- paste(filtered_data$data, collapse = "\n")
    split_data <- read.table(text = data_string, fill = TRUE, header = FALSE, stringsAsFactors = FALSE)
    split_data_df <- data.frame(split_data, stringsAsFactors = FALSE)
    split_data_df <- split_data_df[,1:4]
    split_data_df$Model <- k
    colnames(split_data_df) <- c("ID", "Value1", "Value2", "Significance", "Model")
    
    data_combo <- rbind(data_combo, split_data_df)
    data_combo <- na.omit(data_combo)
    
    
    stats <- c("R-squared", "F-statistic")
    filtered_data2 <- data.frame()
    for (pattern in stats) {
      matches <- a[grepl(pattern, a$data, ignore.case = TRUE), ]
      filtered_data2 <- rbind(filtered_data2, matches)
    }
    filtered_data2$data <- as.character(filtered_data2[[1]])
    
    # Corrected regular expression extraction
    model_stats <- c(k, as.numeric(unlist(regmatches(filtered_data2[[1]], gregexpr("[[:digit:]]+\\.?[[:digit:]]*", filtered_data2[[1]])))))
    model_combo <- rbind(model_combo, model_stats)
  }
  
  names(model_combo) <- c("Model", "R_Squared", "F_Statistic", "Permutations", "Significance")
  
  data_combo_list[[j]] <- data_combo
  model_combo_list[[j]] <- model_combo
}  
  
for(j in 1:4) {  
  
  #j <- 1
  
  data_combo <- data_combo_list[[j]]
  model_combo <- model_combo_list[[j]]
  
  df <- data_combo
  df$Value1 <- as.numeric(df$Value1)
  
  significant_models <- model_combo$Model[model_combo$Significance < 0.05]
  model_combo_subset <- model_combo[model_combo$Model %in% significant_models, ]
  df_subset <- df[df$Model %in% significant_models, ]
  
  excluded <- vector()
  # Ensure all variables are present
  for (id in all_ids) {
    if (!id %in% df_subset$ID) {
      excluded <- c(excluded,id)
      df_subset <- rbind(df_subset, data.frame(ID = id, Value1 = NA, Value2 = NA, Significance = NA, Model = NA))
    }
  }
  
  # Calculate the total counts and significant counts
  total_counts <- table(df_subset$ID)
  total_counts2 <- as.vector(total_counts)
  names(total_counts2) <- names(total_counts)
  
  total_counts2 <- total_counts2[names(total_counts2)!="R-Squared"]
  total_counts2
  
  significant_counts <- table(df_subset$ID[df_subset$Significance < 0.05])
  significant_counts
  
  # Create a vector of zero counts for IDs with no significant results
  significant_counts_full <- rep(0, length(all_ids))
  names(significant_counts_full) <- all_ids
  significant_counts_full[names(significant_counts)] <- significant_counts
  
  significant_counts_full <-  significant_counts_full[names( significant_counts_full)!="R-Squared"]
  significant_counts_full
  
  
  # Calculate the percentages
  percentages <- significant_counts_full / total_runs * 100
  percentages
  # Create labels with total counts and percentages on separate lines
  #labels <- c(paste(patterns, "\n(n=", total_counts[patterns], " / ", round(percentages[patterns], 2), "%)", sep = ""), "R-Squared")
labels <- vector()
  
  # Loop through the patterns
  for (item in patterns) {
    if (item %in% excluded) {
      labels <- c(labels, item)
    } else {
      labels <- c(labels, paste(item, "\n(n=", total_counts2[item], " / ", round(percentages[item], 2), "%)", sep = ""))
    }
  }
  
  # Add "R-Squared" to labels
  labels <- c(labels, "R-Squared")
  
  labels
  
  # Store the generated plots
  plots[[j]] <- ggplot() +
    # First boxplot from df_subset
    geom_boxplot(data = df_subset, aes(x = factor(ID, levels = all_ids), y = Value1, group = ID), fill = "grey", alpha = 0.5) +
    geom_jitter(data = df_subset, aes(x = factor(ID, levels = all_ids), y = Value1, color = Significance < 0.05), width = 0.2) +
    # Second boxplot from model_combo_subset
    geom_boxplot(data = model_combo_subset, aes(x = factor("R-Squared", levels = all_ids), y = as.numeric(R_Squared)), fill = "lightblue", alpha = 0.5) +
    labs(title = plot_titles[j],
         #x = "Variable (count of models in which the variable is present / proportion of models where variable is present that it is significant)",
         x = "",
         y = "R-squared / Model Estimate") +
    scale_x_discrete(labels = c(labels, "R-Squared")) +
    scale_color_manual(name = "Significance < 0.05", values = c("TRUE" = "red", "FALSE" = "blue")) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10, angle = 0, hjust = 1),
      axis.title.y = element_text(size = 10, face = "italic", hjust = 0.5)
    )
}

# Combine all plots into a single figure
combined_plot <- plot_grid(plotlist = plots, ncol = 1, align = "v")

# Print the combined plot
print(combined_plot)






####Simple Mantel: Figure 6 from individual csvs####

library(readxl)
library(ggplot2)
library(dplyr)
library(cowplot)

# Sample data preparation


#this is only necessary if you're running an unequal number of analyses - e.g. hold problems taking ages on HTCondor
file_numbers <- list()
files <- list.files("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/PC1", 
                    pattern = "^PC1_Simple[0-9]+\\.csv$", full.names = TRUE)
file_numbers[[1]] <- as.numeric(gsub("^(PC1_Simple)(\\d+)\\.csv$", "\\2", basename(files)))

files <- list.files("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/PC2", 
                    pattern = "^PC2_Simple[0-9]+\\.csv$", full.names = TRUE)
file_numbers[[2]] <- as.numeric(gsub("^(PC2_Simple)(\\d+)\\.csv$", "\\2", basename(files)))

files <- list.files("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/PC3", 
                    pattern = "^PC3_Simple[0-9]+\\.csv$", full.names = TRUE)
file_numbers[[3]] <- as.numeric(gsub("^(PC3_Simple)(\\d+)\\.csv$", "\\2", basename(files)))

files <- list.files("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/PC4", 
                    pattern = "^PC4_Simple[0-9]+\\.csv$", full.names = TRUE)
file_numbers[[4]] <- as.numeric(gsub("^(PC4_Simple)(\\d+)\\.csv$", "\\2", basename(files)))


big_df_list <- list()
for(j in 1:4) {
  df_list <- list()
  i1 <- file_numbers[[j]]
  i2 <- order(i1)
  for(k in file_numbers[[j]]) {
    df_list[[i2[[which(i1==k)]]]] <- read.csv(paste0("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/PC",j,"/PC",j,"_Simple",k,".csv"))

    }
    big_df_list[[j]] <- df_list}

mantel_means <- list()
r2_means <- list()
for(j in 1:4){
  # Assuming your list of data frames is called 'df_list'
  combined_df <- bind_rows(big_df_list[[j]], .id = "source")  # Adds an ID column for traceability
  
  # Box plot of p-value grouped by variable name
  ggplot(combined_df, aes(x = X, y = pvalue)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x = "Variable Name", y = "P-value", title = "Box Plot of P-value by Variable Name") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
  
  
  mean_pvalue_by_var <- combined_df %>%
    group_by(`X`) %>%
    summarize(mean_pvalue = mean(pvalue, na.rm = TRUE))
  
  mean_r2_by_var <- combined_df %>%
    group_by(`X`) %>%
    summarize(mean_r2 = mean(statistic, na.rm = TRUE))
  
  print(mean_pvalue_by_var)
  
  mantel_means[[j]] <- mean_pvalue_by_var
  r2_means[[j]] <- mean_r2_by_var
}

combined <- cbind(mantel_means[[1]][1], do.call(cbind, lapply(mantel_means, function(df) df[-1])))
colnames(combined) <- c("Variable", "PC1", "PC2", "PC3", "PC4")
write.csv(combined, "C:/Users/jblin/OneDrive - The University of Liverpool/NITD/Table6.csv")

combined2 <- cbind(r2_means[[1]][1], do.call(cbind, lapply(r2_means, function(df) df[-1])))
colnames(combined2) <- c("Variable", "PC1", "PC2", "PC3", "PC4")
write.csv(combined2, "C:/Users/jblin/OneDrive - The University of Liverpool/NITD/Table6b.csv")

library(ggplot2)
library(dplyr)

list_of_lists <- big_df_list

# Assuming your list of lists is named 'list_of_lists'
combined_df <- bind_rows(
  lapply(seq_along(list_of_lists), function(i) {
    bind_rows(list_of_lists[[i]]) %>%
      mutate(PC = paste0("PC", i))  # Assigning group labels
  })
)

library(dplyr)

# Identify significance and check if ALL p-values within each PC-Variable group are ≤ 0.05
group_flags <- combined_df %>%
  group_by(PC, X) %>%
  summarize(has_sig_pvalue = any(pvalue <= 0.05),  # At least one significant p-value
            all_significant = all(pvalue <= 0.05),  # All p-values are significant
            .groups = "drop")

# Merge flags back into the main dataframe
combined_df <- combined_df %>%
  left_join(group_flags, by = c("PC", "X"))

# Assign colors based on significance conditions
combined_df <- combined_df %>%
  mutate(color_category = case_when(
    all_significant  ~ "dark red",   # ALL values ≤ 0.05
    has_sig_pvalue  ~ "red",         # Some values ≤ 0.05, but not all
    !all_significant ~ "dark blue",  # No values ≤ 0.05, but keeping distinction
    TRUE            ~ "blue"         # Default for varied non-significant cases
  ))

combined_df$X <- gsub("prop", "P.I.", combined_df$X)

# Create the box plot with updated color logic
simple_mantel <- ggplot(combined_df, aes(x = X, y = pvalue, color = color_category)) +
  geom_boxplot() +
  scale_color_manual(values = c("dark red" = "#8B0000",
                                "red" = "red",
                                "dark blue" = "#00008B",
                                "blue" = "blue"),
                     labels = c("dark red" = "All Significant",
                                "red" = "Some Significant",
                                "dark blue" = "None Significant",
                                "blue" = "None Significant (Varied)")) +
  theme_minimal() +
  labs(x = "Variable Name", y = "P-value", title = "P-values by Variable Across PCs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +  # Move legend to the bottom
  guides(color = guide_legend(title = NULL)) +  # Remove legend title
  facet_wrap(~ PC, ncol = 1)

ggsave("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/simple_mantel.pdf", plot = simple_mantel, width = 8.27, height = 11.69, units = "in")

####Multiple Matrix Regression: Figure 7 from individual csvs####

library(readxl)
library(ggplot2)
library(dplyr)
library(cowplot)

# Sample data preparation
patterns <- c("age", "altitude", "bio01", "bio04", "bio12", "bio15", "costpath", "count", "prop", "rawmat", "rugosity", "species", "type")
all_ids <- c(patterns, "R-Squared")
plots <- list()
plot_titles <- c("Box Plot of Model Estimates for each variable and R-Squared values for PC1",
                 "Box Plot of Model Estimates for each variable and R-Squared values for PC2",
                 "Box Plot of Model Estimates for each variable and R-Squared values for PC3",
                 "Box Plot of Model Estimates for each variable and R-Squared values for PC4")


plot_titles <- c("PC1", "PC2", "PC3", "PC4")

data_combo_list <- list()
model_combo_list <- list()

#this is only necessary if you're running an unequal number of analyses - e.g. hold problems taking ages on HTCondor
file_numbers <- list()
files <- list.files("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/PC1", 
                    pattern = "^PC1_Multiple[0-9]+\\.csv$", full.names = TRUE)
file_numbers[[1]] <- as.numeric(gsub("^(PC1_Multiple)(\\d+)\\.csv$", "\\2", basename(files)))

files <- list.files("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/PC2", 
                    pattern = "^PC2_Multiple[0-9]+\\.csv$", full.names = TRUE)
file_numbers[[2]] <- as.numeric(gsub("^(PC2_Multiple)(\\d+)\\.csv$", "\\2", basename(files)))

files <- list.files("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/PC3", 
                    pattern = "^PC3_Multiple[0-9]+\\.csv$", full.names = TRUE)
file_numbers[[3]] <- as.numeric(gsub("^(PC3_Multiple)(\\d+)\\.csv$", "\\2", basename(files)))

files <- list.files("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/PC4", 
                    pattern = "^PC4_Multiple[0-9]+\\.csv$", full.names = TRUE)
file_numbers[[4]] <- as.numeric(gsub("^(PC4_Multiple)(\\d+)\\.csv$", "\\2", basename(files)))

for(j in 1:4) {
  data_combo <- data.frame()
  model_combo <- data.frame()

  for(k in file_numbers[[j]]) {
    d <- read.csv(paste0("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/PC",j,"/PC",j,"_Multiple",k,".csv"))
    a <- data.frame(d[[2]]) #from saved excel
    colnames(a) <- "data"
    
    filtered_data <- data.frame()
    for (pattern in patterns) {
      matches <- a[grepl(pattern, a$data, ignore.case = TRUE), ]
      filtered_data <- rbind(filtered_data, matches)
    }
    filtered_data$data <- as.character(filtered_data[[1]])
    
    data_string <- paste(filtered_data$data, collapse = "\n")
    split_data <- read.table(text = data_string, fill = TRUE, header = FALSE, stringsAsFactors = FALSE)
    split_data_df <- data.frame(split_data, stringsAsFactors = FALSE)
    split_data_df <- split_data_df[,1:4]
    split_data_df$Model <- k
    colnames(split_data_df) <- c("ID", "Value1", "Value2", "Significance", "Model")
    
    data_combo <- rbind(data_combo, split_data_df)
    data_combo <- na.omit(data_combo)
    
    stats <- c("R-squared", "F-statistic")
    filtered_data2 <- data.frame()
    for (pattern in stats) {
      matches <- a[grepl(pattern, a$data, ignore.case = TRUE), ]
      filtered_data2 <- rbind(filtered_data2, matches)
    }
    filtered_data2$data <- as.character(filtered_data2[[1]])
    
    model_stats <- c(k, as.numeric(unlist(regmatches(filtered_data2[[1]], gregexpr("[[:digit:]]+\\.?[[:digit:]]*", filtered_data2[[1]])))))
    model_combo <- rbind(model_combo, model_stats)
  }

  names(model_combo) <- c("Model", "R_Squared", "F_Statistic", "Permutations", "Significance")
  
  data_combo_list[[j]] <- data_combo
  model_combo_list[[j]] <- model_combo
}  

# Compute means for each dataframe
summary_table <- do.call(rbind, lapply(seq_along(model_combo_list), function(i) {
  df <- model_combo_list[[i]]
  data.frame(
    Model = paste0("PC", i),
    R_Squared = mean(df$R_Squared, na.rm = TRUE),
    F_Statistic = mean(df$F_Statistic, na.rm = TRUE),
    Significance = mean(df$Significance, na.rm = TRUE)
  )
}))

write.csv(summary_table, "C:/Users/jblin/OneDrive - The University of Liverpool/NITD/MMR_Summary_table.csv")

table8 <- list()
table9 <- list()

# Load dplyr for clean grouping and summarizing
library(dplyr)

# Group by ID and compute means for each dataframe
grouped_summary <- lapply(seq_along(data_combo_list), function(i) {
  data_combo_list[[i]] %>%
    group_by(ID) %>%
    summarise(
      Value1 = mean(Value1, na.rm = TRUE),
      Value2 = mean(Value2, na.rm = TRUE),
      Significance = mean(Significance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Model = paste0("PC", i)) %>%
    dplyr::select(Model, everything())
})

# Combine all summaries into one neat table
final_table <- bind_rows(grouped_summary)

library(tidyr)

wide_table <- final_table %>%
  pivot_wider(
    names_from = Model,
    values_from = c(Value1, Value2, Significance),
    names_sep = "_"
  )

write.csv(wide_table, "C:/Users/jblin/OneDrive - The University of Liverpool/NITD/MMR_Summary_data_table.csv")

for(j in 1:4) {
  data_combo <- data_combo_list[[j]]
  model_combo <- model_combo_list[[j]]
  
  df <- data_combo
  df$Value1 <- as.numeric(df$Value1)
  
  significant_models <- model_combo$Model[model_combo$Significance < 0.05]
  model_combo_subset <- model_combo[model_combo$Model %in% significant_models, ]
  df_subset <- df[df$Model %in% significant_models, ]
  
  excluded <- vector()
  for (id in all_ids) {
    if (!id %in% df_subset$ID) {
      excluded <- c(excluded,id)
      df_subset <- rbind(df_subset, data.frame(ID = id, Value1 = NA, Value2 = NA, Significance = NA, Model = NA))
    }
  }
  
  # Calculate the total counts and significant counts
  total_counts <- table(df_subset$ID)
  total_counts2 <- as.vector(total_counts)
  names(total_counts2) <- names(total_counts)
  total_counts2 <- total_counts2[names(total_counts2)!="R-Squared"]
  significant_counts <- table(df_subset$ID[df_subset$Significance < 0.05])
  
  # Create a vector of zero counts for IDs with no significant results
  significant_counts_full <- rep(0, length(all_ids))
  names(significant_counts_full) <- all_ids
  significant_counts_full[names(significant_counts)] <- significant_counts
  
  significant_counts_full <-  significant_counts_full[names( significant_counts_full)!="R-Squared"]
  
  total_runs <- length(file_numbers[[j]])
  # Calculate the percentages
  percentages <- significant_counts_full / total_runs * 100
  percentages
  labels <- vector()
  
  table8[[j]] <- significant_counts_full
  table9[[j]] <- total_counts2
  
  # Loop through the patterns
  #for (item in patterns) {
  #  if (item %in% excluded) {
  #    labels <- c(labels, item)
  #  } else {
  #    labels <- c(labels, paste(item, "\n(n=", total_counts2[item], " / ", round(percentages[item], 2), "%)", sep = ""))
  #  }
  #}
  
  #labels <- c(labels, "R-Squared")
  labels <- c("Age","Altitude", "Bio01","Bio04","Bio12","Bio15","Costpath","Count", "P.I.","Rawmat","Rugosity","Species","Type", "R-squared")
  # Store the generated plots
  plots[[j]] <- ggplot() +
    # First boxplot from df_subset
    geom_boxplot(data = df_subset, aes(x = factor(ID, levels = all_ids), y = Value1, group = ID), fill = "grey", alpha = 0.5) +
    geom_jitter(data = df_subset, aes(x = factor(ID, levels = all_ids), y = Value1, color = Significance < 0.05), width = 0.2) +
    # Second boxplot from model_combo_subset
    geom_boxplot(data = model_combo_subset, aes(x = factor("R-Squared", levels = all_ids), y = as.numeric(R_Squared)), fill = "lightblue", alpha = 0.5) +
    labs(title = plot_titles[j],
         #x = "Variable (count of models in which the variable is present / proportion of models where variable is present that it is significant)",
         x = "",
         y = "R-squared / Model Estimate") +
    scale_x_discrete(labels = c(labels)) +#, "R-Squared")) +
    scale_color_manual(name = "Significance < 0.05", values = c("TRUE" = "red", "FALSE" = "blue")) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10, angle = 0, hjust = 1),
      axis.title.y = element_text(size = 10, face = "italic", hjust = 0.5)
    )
}

# Combine all plots into a single figure
combined_plot <- plot_grid(plotlist = plots, ncol = 1, align = "v")

# Print the combined plot
print(combined_plot)

ggsave("C:/Users/jblin/OneDrive - The University of Liverpool/NITD/multi_mantel.pdf", plot = combined_plot, width = 11.69, height = 8.79, units = "in")

dfx <- data.frame(table8)
names(dfx) <- plot_titles

dfy <- data.frame(table9)
names(dfy) <- plot_titles

write.csv(cbind(dfx, dfy), "C:/Users/jblin/OneDrive - The University of Liverpool/NITD/proptab.csv")

####PC1-4 SpCode Logistic regressions####

library(psych)
library(ggplot2)
library(ggpubr)
library(GGally)
library(ROCR)

#use AIC to select variables and PC to examine

AICs <- c()
xs <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8")
for(j in 1:8){
  f <- as.formula(paste("SpCode~",paste(xs[1:j],collapse="+")))
  logit_levallois <- glm(f,
                         data=levallois_proj[-unknown_sp,],
                         family=binomial(link="logit"))
  AICs <- c(AICs,logit_levallois$aic)}

LogitAICs <- data.frame(PC=1:8,AIC=round(AICs,0))

#quick cross validation for the logit model to create a meaningful ROC curve
d <- levallois_proj[-unknown_sp,]
n <- nrow(d)
##repeat many times to get confidence envelope for the ROC curve
##the following code is inefficient, but the rocr ouput needs to be captured
##for use with ggplot instead of rocr's default plots

#PC1-4

rocr_cv_perf <- data.frame(r=numeric(),tpr=numeric(),fpr=numeric())
rocr_cv_perf_auc <- data.frame(r=numeric(),tpr=numeric(),fpr=numeric())
for( j in 1:1000){
  random_sample_index <- base::sample(1:n,n/2)
  d_train <- d[random_sample_index,]
  d_cv <- d[-random_sample_index,]
  logit_levallois_train <- glm(SpCode~PC1+PC2+PC3+PC4,
                               data=d_train,
                               family=binomial(link="logit")
  )
  logit_levallois_cv_pred <- predict(logit_levallois_train,type="response",newdata=d_cv)
  rocr_pred <- prediction(logit_levallois_cv_pred,d_cv$SpCode)
  rocr_perf <- performance(rocr_pred,measure="tpr",x.measure="fpr")
  nn <- length(rocr_perf@x.values[[1]])
  rocr_cv_perf <- rbind(rocr_cv_perf,data.frame(r=rep(j,nn),tpr=rocr_perf@y.values[[1]],fpr=rocr_perf@x.values[[1]]))
  rocr_perf <- performance(rocr_pred,measure="auc")
  nn <- length(rocr_perf@y.values[[1]])
  rocr_cv_perf_auc <- rbind(rocr_cv_perf_auc,data.frame(r=rep(j,nn),auc=rocr_perf@y.values[[1]]))
}

#plotting multiple ROC curves from CV analysis
PC1_4_curve <- ggplot(data = rocr_cv_perf,
       aes(y=tpr,x=fpr,group=r)) +
  geom_line(alpha=0.01) +
  geom_line(data=data.frame(x=seq(0,1,0.001),y=seq(0,1,0.001)),
            aes(x=x,y=y),
            linetype=2,
            colour="grey",
            inherit.aes=F) +
  labs(title="ROC Curves from Logit Cross Validation (PC1-4)",
       y="TPR",
       x="FPR") +
  theme_minimal() +
  theme(text = element_text(family="Times", size=12),
        plot.title = element_text(face="bold",hjust=0.5,size=15))

PC1_4_curve

ggsave(filename="C:/Users/jblin/OneDrive - The University of Liverpool/NITD/PC1-4logit_cv_roc.pdf",
       device = "pdf")

##
####PC1-4 SPCode AUC histogram####
n <- length(rocr_cv_perf_auc$auc)
nbins <- 1 + (3.222 * log(n))

nbins <-23
PC1_4_hist <- ggplot(data = data.frame(AUC=rocr_cv_perf_auc$auc),
       mapping = aes(x = AUC, stat(density))) +
  geom_histogram(bins=nbins,
                 color="white") +
  labs(x="AUC",
       y="Density",
       title="AUC from Logit Cross Validation (PC1-4)") +
  theme_minimal() +
  theme(#text = element_text(family="Times", size=12),
    plot.title = element_text(face="bold",hjust=0.5,size=15))


PC1_4_hist
ggsave(filename="C:/Users/jblin/OneDrive - The University of Liverpool/NITD/PC1-4auc_distribution.pdf",
       device = "pdf")


logit_levallois <- glm(SpCode~PC1+PC2+PC3+PC4,
                       data=levallois_proj[-unknown_sp,],
                       family=binomial(link="logit"))

logit_predict_JKF1 <- predict(logit_levallois,
                              newdata=levallois_proj[which(levallois_proj$Site=="JKF-1"),],
                              type="response",
                              se.fit=T)

summary(logit_predict_JKF1$fit)

logit_predict_ALM3 <- predict(logit_levallois,
                              newdata=levallois_proj[which(levallois_proj$Site=="ALM-3"),],
                              type="response",
                              se.fit=T)

summary(logit_predict_ALM3$fit)



####PC2-3 SpCode Logistic regressions####

#USE MMR OUTPUTS TO SELECT PC" & PC3 FOR ANALYSIS

rocr_cv_perf <- data.frame(r=numeric(),tpr=numeric(),fpr=numeric())
rocr_cv_perf_auc <- data.frame(r=numeric(),tpr=numeric(),fpr=numeric())
for( j in 1:1000){
  random_sample_index <- base::sample(1:n,n/2)
  d_train <- d[random_sample_index,]
  d_cv <- d[-random_sample_index,]
  logit_levallois_train <- glm(SpCode~PC2+PC3,
                               data=d_train,
                               family=binomial(link="logit")
  )
  logit_levallois_cv_pred <- predict(logit_levallois_train,type="response",newdata=d_cv)
  rocr_pred <- prediction(logit_levallois_cv_pred,d_cv$SpCode)
  rocr_perf <- performance(rocr_pred,measure="tpr",x.measure="fpr")
  nn <- length(rocr_perf@x.values[[1]])
  rocr_cv_perf <- rbind(rocr_cv_perf,data.frame(r=rep(j,nn),tpr=rocr_perf@y.values[[1]],fpr=rocr_perf@x.values[[1]]))
  rocr_perf <- performance(rocr_pred,measure="auc")
  nn <- length(rocr_perf@y.values[[1]])
  rocr_cv_perf_auc <- rbind(rocr_cv_perf_auc,data.frame(r=rep(j,nn),auc=rocr_perf@y.values[[1]]))
}

#plotting multiple ROC curves from CV analysis
PC2_3_curve <- ggplot(data = rocr_cv_perf,
                      aes(y=tpr,x=fpr,group=r)) +
  geom_line(alpha=0.01) +
  geom_line(data=data.frame(x=seq(0,1,0.001),y=seq(0,1,0.001)),
            aes(x=x,y=y),
            linetype=2,
            colour="grey",
            inherit.aes=F) +
  labs(title="ROC Curves from Logit Cross Validation (PC2-3)",
       y="TPR",
       x="FPR") +
  theme_minimal() +
  theme(text = element_text(family="Times", size=12),
        plot.title = element_text(face="bold",hjust=0.5,size=15))

PC2_3_curve
ggsave(filename="C:/Users/jblin/OneDrive - The University of Liverpool/NITD/PC2-3logit_cv_roc.pdf",
       device = "pdf")


####PC2-3 SPCode AUC histogram####
n <- length(rocr_cv_perf_auc$auc)
nbins <- 1 + (3.222 * log(n))

nbins <-23
PC2_3_hist <- ggplot(data = data.frame(AUC=rocr_cv_perf_auc$auc),
                     mapping = aes(x = AUC, stat(density))) +
  geom_histogram(bins=nbins,
                 color="white") +
  labs(x="AUC",
       y="Density",
       title="AUC from Logit Cross Validation (PC2-3)") +
  theme_minimal() +
  theme(#text = element_text(family="Times", size=12),
    plot.title = element_text(face="bold",hjust=0.5,size=15))

PC2_3_hist

ggsave(filename="C:/Users/jblin/OneDrive - The University of Liverpool/NITD/PC2-3auc_distribution.pdf",
       device = "pdf")


logit_levallois_PC2_3 <- glm(SpCode~PC2+PC3,
                       data=levallois_proj[-unknown_sp,],
                       family=binomial(link="logit"))

logit_predict_JKF1_PC2_3 <- predict(logit_levallois_PC2_3,
                              newdata=levallois_proj[which(levallois_proj$Site=="JKF-1"),],
                              type="response",
                              se.fit=T)

summary(logit_predict_JKF1$fit)
summary(logit_predict_JKF1_PC2_3$fit)

logit_predict_ALM3_PC2_3 <- predict(logit_levallois_PC2_3,
                              newdata=levallois_proj[which(levallois_proj$Site=="ALM-3"),],
                              type="response",
                              se.fit=T)

summary(logit_predict_ALM3$fit)
summary(logit_predict_ALM3_PC2_3$fit)


logit_df <- as.data.frame(rbind(summary(logit_predict_JKF1$fit), 
                       summary(logit_predict_JKF1_PC2_3$fit),
                        summary(logit_predict_ALM3$fit), 
                        summary(logit_predict_ALM3_PC2_3$fit)))

row.names(logit_df) <- c("JKF-1 (PC1-4)",
                         "JKF-1 (PC2-3)",
                         "ALM-3 (PC1-4)",
                         "ALM-3 (PC2-3)")

write.csv(logit_df, "C:/Users/jblin/OneDrive - The University of Liverpool/NITD/logit_df.csv")


plot(logit_predict_JKF1_PC2_3$fit)

sort(logit_predict_JKF1_PC2_3$fit)
sort(logit_predict_ALM3_PC2_3$fit)

library(patchwork)

(PC1_4_curve + PC2_3_curve) / (PC1_4_hist + PC2_3_hist)

####End####