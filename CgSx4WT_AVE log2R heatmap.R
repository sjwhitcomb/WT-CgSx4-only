# line AVE log2R heat map
# R = plant/AVE in WT at same timepoint 
# log2R
# mean(log2R) for each CgSx4 x timepoint

# log2R data of SOLUBLE metabolites-----------------------------------
# N = 4
# AVE log2R of CgSx4 x timepoint
# rows = CgSx4 x timepoint
# columns = metabolites

df123s.sol <- read.delim("CgSx4WT_AVE log2R_soluble met x Tgline x timepoint.txt", header = TRUE) # leaf (1,2,3 month harvests) and mature seed harvest
df123.sol <- subset(df123s.sol, df123s.sol$tissue == "leaf") # leaf only df
dfs.sol <- subset(df123s.sol, df123s.sol$tissue == "seed") # seed only df

##  matricies of metabolite data for clustering and heatmaps

# leaf and seed matrix
m123s.sol <- as.matrix(df123s.sol[,c(5:dim(df123s.sol)[2])]) # only metabolite columns
rownames(m123s.sol) <- df123s.sol$group
t.m123s.sol <- t(m123s.sol)
rownames(t.m123s.sol) # metabolites
colnames(t.m123s.sol) # CgSx4 x timepoint groups

# leaf only martrix 
m123.sol <- as.matrix(df123.sol[,c(5:dim(df123.sol)[2])]) # only metabolite columns
rownames(m123.sol) <- df123.sol$group
t.m123.sol <- t(m123.sol)
rownames(t.m123.sol) # metabolites
colnames(t.m123.sol) # CgSx4 x timepoint groups

# seed only martrix
ms.sol <- as.matrix(dfs.sol[,c(7:13, 15:28)]) # only metabolite columns
# only NAs for Hse in seeds, so remove from matrix!
rownames(ms.sol) <- dfs.sol$group
t.ms.sol <- t(ms.sol)
rownames(t.ms.sol) # metabolites
colnames(t.ms.sol) # CgSx4 x timepoint groups

# log2R data of PROTEIN-INCORPORATED amino acids-----------------------------------
# N varies!
# N = 6 for WT @ 1, 2, 3 months
# N = 3 CgSx4 in leaves
# EXCEPT: 2month CgSx4 = 2
# N = 3 for WT, CgSx4 in seeds
# AVE log2R of CgSx4 x timepoint
# rows = CgSx4 x timepoint
# columns = metabolites


df123s.prot <- read.delim("CgSx4WT_AVE log2R_protein-incorporated aa x Tgline x timepoint.txt", header = TRUE) # leaf (1,2,3 month harvests) and mature seed harvest
df123.prot <- subset(df123s.prot, df123s.prot$tissue == "leaf") # leaf only df
dfs.prot <- subset(df123s.prot, df123s.prot$tissue == "seed") # seed only df

##  matricies of metabolite data for clustering and heatmaps

# leaf and seed matrix
m123s.prot <- as.matrix(df123s.prot[,c(5:dim(df123s.prot)[2])]) # only metabolite columns
rownames(m123s.prot) <- df123s.prot$group
t.m123s.prot <- t(m123s.prot)
rownames(t.m123s.prot) # metabolites
colnames(t.m123s.prot) # Tg x timepoint groups

# leaf only martrix 
m123.prot <- as.matrix(df123.prot[,c(5:dim(df123.prot)[2])]) # only metabolite columns
rownames(m123.prot) <- df123.prot$group
t.m123.prot <- t(m123.prot)
rownames(t.m123.prot) # metabolites
colnames(t.m123.prot) # Tg x timepoint groups

# seed only martrix
ms.prot <- as.matrix(dfs.prot[,c(5:dim(dfs.prot)[2])])
rownames(ms.prot) <- dfs.prot$group
t.ms.prot <- t(ms.prot)
rownames(t.ms.prot) # metabolites
colnames(t.ms.prot) # Tg x timepoint groups


# log2R data of SOLUBLE metabolites and PROTEIN-INCORPORATED amino acids------------
# N varies!
# AVE log2R of CgSx4 x timepoint
# rows = CgSx4 x timepoint
# columns = metabolites

df123s.sol.prot <- read.delim("CgSx4WT_AVE log2R_soluble met and protein-incorporated aa x Tgline x timepoint.txt", header = TRUE) # leaf (1,2,3 month harvests) and mature seed harvest

df123.sol.prot <- subset(df123s.sol.prot, df123s.sol.prot$tissue == "leaf") # leaf only df
dfs.sol.prot <- subset(df123s.sol.prot, df123s.sol.prot$tissue == "seed") # seed only df

# seed only martrix 
ms.sol.prot <- as.matrix(dfs.sol.prot[,c(8:16, 18:27, 29, 30)]) # only metabolite columns
rownames(ms.sol.prot) <- dfs.sol.prot$group
t.ms.sol.prot <- t(ms.sol.prot)
rownames(t.ms.sol.prot) # metabolites
colnames(t.ms.sol.prot) # Tg x timepoint groups


# packages for heirarchical clustering and heatmaps --------------------------------

library(ComplexHeatmap)
# Distance calculation methods
##"euclidean" = default
## "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
# Clustering methods
## "complete" = default?
## "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
library(dendextend)
library(circlize)
## colorRamp2(c(-3, 0, 3), c("green", "white", "red") ## values between -3 and 3 are linearly interpolated to obtain corresponding colors, values larger than 3 are all mapped to red and values less than -3 are all mapped to green (so the color mapping function demonstrated here is robust to outliers)
library(viridis)
## creates a vector of n equally spaced colors, designed in such a way that it will analytically be perfectly perceptually-uniform, both in regular form and also when converted to black-and-white. It is also designed to be perceived by readers with the most common form of color blindness

# make heatmaps -----------------

####### 123 leaves only

# SOLUBLE metablites only
# vertical matrix (metabolite rows, Tg_timepoint columns)
leaves.sol <- Heatmap(t.m123.sol, colorRamp2(c(-3, 0, 3), c("purple", "white", "red")),
                      name = "Log2R",
                      cluster_columns = TRUE,
                      column_dend_reorder = FALSE,
                      column_names_side = "top",
                      column_dend_height = unit(1.5, "cm"),
                      column_names_gp=gpar(fontsize = 9),
                      cluster_rows = TRUE,
                      clustering_distance_rows = "euclidean",
                      clustering_method_rows = "complete",
                      row_names_side = "left",
                      row_dend_width = unit(2, "cm"),
                      row_names_gp=gpar(fontsize = 9),
                      heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  color_bar = "continuous", 
                                                  legend_direction = "horizontal",
                                                  legend_width = unit(30, "mm")))


pdf("CgSx4WT_heatmap_AVElog2R_leaves x soluble metabolites.pdf",
     width= 6.30/2,
     height = 6.30*1.1)
draw(leaves.sol, show_heatmap_legend = TRUE, heatmap_legend_side = "bottom")
dev.off()

####### seed only

# SOLUBLE metabolites and PROTEIN-incorporated amino acids
# vertical matrix (metabolite rows, Tg_timepoint columns)
seeds <- Heatmap(t.ms.sol.prot, colorRamp2(c(-1.5, 0, 1.5), c("darkgreen", "white", "orange")),
                 name = "Log2R",
                 cluster_columns = FALSE,
                 column_names_side = "top",
                 column_names_gp=gpar(fontsize = 9.4),
                 cluster_rows = TRUE,
                 clustering_distance_rows = "euclidean",
                 clustering_method_rows = "complete",
                 row_names_side = "left",
                 row_dend_width = unit(2, "cm"),
                 row_names_gp=gpar(fontsize = 9.4),
                 heatmap_legend_param = list(title_gp = gpar(fontsize = 8.3),
                                             labels_gp = gpar(fontsize = 7.3),
                                             color_bar = "continuous", 
                                             legend_direction = "horizontal", 
                                             legend_width = unit(30, "mm")))

pdf("CgSx4WT_heatmap_AVElog2R_seeds x soluble and prot.pdf",
    width= 6.30/3,
    height = 6.30)
draw(seeds, show_heatmap_legend = TRUE, heatmap_legend_side = "bottom")
dev.off()
