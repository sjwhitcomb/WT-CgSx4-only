# log2(AVE 1monthWT normalized concentration of soluble metabolite) for CgSx4 x leaf timepoint

# AVE log2_WT1m normalizeddata of SOLUBLE metabolites-----------------------------------
# rows = CgSx4 x timepoint
# columns = metabolites

df123.sol <- read.delim("CgSx4WT_log2 AVE Norm WT1m_leaf_soluble metabolites.txt", header = TRUE)

# matrix (WT + CgSx4) x 3 leaf timepoints = 6
m123.sol6 <- as.matrix(df123.sol[,c(6:dim(df123.sol)[2])]) # only metabolite columns
rownames(m123.sol6) <- df123.sol$group
t.m123.sol6 <- t(m123.sol6)
rownames(t.m123.sol6) # metabolites
colnames(t.m123.sol6) # Tg x timepoint groups

# matrix (CgSx4) x 1m + (WT + CgSx4) x 2m and 3m = 5
m123.sol5 <- m123.sol6[-1,]
t.m123.sol5 <- t(m123.sol5)

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

# make heatmaps -------------------------------------------
# cluster rows only matrix.5
# viridis color scale        
Heatmap(t.m123.sol5, col = colorRamp2(c(-3, 0, 3), viridis(3, opt = "D")),
        name = "Log2 \n(WT 1m\nnorm.)",
        cluster_rows = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "complete",
        row_names_side = "left",
        row_names_gp=gpar(cex=1),
        row_dend_side = "left",
        row_dend_width = unit(3, "cm"),
        cluster_columns = FALSE,
        column_order = c("CgSx4_leaf.1m","WT_leaf.2m","CgSx4_leaf.2m", "WT_leaf.3m","CgSx4_leaf.3m"),
        column_names_side = "top",
        column_names_gp=gpar(cex=1))


# cluster rows only matrix.12
# viridis color scale        
hm6 <- Heatmap(t.m123.sol6, col = colorRamp2(c(-3, 0, 3), viridis(3, opt = "D")),
                name = "Log2 \n(WT 1m\nnorm.)",
                cluster_rows = TRUE,
                clustering_distance_rows = "euclidean",
                clustering_method_rows = "complete",
                row_names_side = "left",
                row_names_gp=gpar(fontsize = 9),
                row_dend_side = "left",
                row_dend_width = unit(2, "cm"),
                cluster_columns = FALSE,
                column_order = c("WT_leaf.1m", "CgSx4_leaf.1m","WT_leaf.2m","CgSx4_leaf.2m","WT_leaf.3m","CgSx4_leaf.3m"),
                column_names_side = "top",
                column_names_gp=gpar(fontsize = 9),
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                            labels_gp = gpar(fontsize = 7),
                                            color_bar = "continuous", 
                                            legend_direction = "vertical", 
                                            legend_height = unit(40, "mm")))
hm6

pdf("CgSx4WT SupFig2.pdf", width = 4, height= 6.5)
draw(hm6,heatmap_legend_side = "left")
dev.off()


