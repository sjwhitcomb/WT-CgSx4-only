## NOTE!! data from CgS vigor experiment (not Cuong's experiment); SPRING-240214-9

# DATA ------------------------

#### seed yield 
# FW in g
# N = 6
# all plants! including those that had black fungus on their panicles
# collected seeds only from "mature" panicles, then dried and cleaned but not dehulled
df.yield<-read.delim("CgSx4WT_seed yield_all plants_line.txt", header = TRUE)
df.yield$line<-factor(df.yield$line, levels = c("WT", "CgSx4"))

#### 1000 seed DW
# DW in g
# N = 4 TECHNICAL replicates! (pooled seeds from multiple plants of same genotype)
# collected seeds only from "mature" panicles, then dried and cleaned but not dehulled

df.seedweight<-read.delim("CgSx4WT_1000seedDW_SPRING-240214-9_technical replicates.txt", header = TRUE)
df.seedweight$line<-factor(df.seedweight$line, levels = c("WT", "CgSx4"))

#### tillers
# number of tillers 40 days after transplantation
# WT N = 5
# CgSx4 N = 6
df.tillers<-read.delim("CgSx4WT_tillers d40_line.txt", header = TRUE)
df.tillers$line<-factor(df.tillers$line, levels = c("WT", "CgSx4"))

#### final vetetative DW
# DW in g
# WT N = 5
# CgSx4 N = 6
df.veg<-read.delim("CgSx4WT_final vegetative DW_line.txt", header = TRUE)
df.veg$line<-factor(df.veg$line, levels = c("WT", "CgSx4"))


# create ggplot objects -------------------
library(ggplot2)
library(grid)

# function agronomic.univscatter() to create univariate scatter plot OBJECT
# for AGRONOMIC data
# no legend

# df = df
# var = df$var
# title = ""
# units = "", for y axis title
# jw = jitter width
# jh = jitter height
# ylim = number or NA
# return is ggplot OBJECT


agronomic.univscatter <- function (df, var, title, units, ylim, jw, jh){
  pj <- position_jitter(width = jw, 
                        height = jh)
  yLabel <- units
  
  usp <- ggplot(df, aes(x = line, y = var)) +
    geom_point(aes(fill = line), shape = 25, size = 1.75, stroke = 0.25, alpha = 0.8, position = pj) +
    scale_fill_manual(values=c("black", "grey")) +
    xlab("") +
    ylab(yLabel) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 8)) +
    theme(axis.title.x=element_text(size=7)) +
    theme(axis.title.y=element_text(size=7)) +
    theme(axis.text.x=element_text(color="black", size=6, angle = 90, hjust = 0.8, vjust = 0.5)) +
    theme(axis.text.y=element_text(color="black", size=6)) +
    theme(panel.background = element_rect(fill = "white", color = "black")) +
    theme(panel.border = element_rect(colour="black", fill=NA)) +
    scale_y_continuous(limits = c(0, ylim)) +
    theme(legend.position="none")
  return(usp)
}


## FOR FIGURE 2
y <- agronomic.univscatter(df.yield, df.yield$seedyield, "seed yield", "g DW", NA, 0.1, 0.05)
sw <- agronomic.univscatter(df.seedweight, df.seedweight$seed.1000, "1000 seed weight", "g DW", 40, 0.1, 0.05)
t <- agronomic.univscatter(df.tillers, df.tillers$tillers, "tillers", "number at 40 days", 15, 0.1, 0.05)
vw <- agronomic.univscatter(df.veg, df.veg$finalVegDW, "final vegetative DW", "g DW", 75, 0.1, 0.05)


# equalize widths and arrange graphs in grid -------------------
library(grid)    # for pmax
library(gridExtra) # to arrange the plots
library(ggplot2)   # to construct the plots
library(gtable)   # to add columns to gtables of plots without legends


equalize <- function(plots) {
  g <- lapply(plots, ggplotGrob) # Convert to gtables
  g.widths = lapply(g, function(x) grid:::unit.list(x$widths)) # Apply the un-exported unit.list() from grid to each plot
  g3.widths <- lapply(g.widths, function(x) x[1:3]) # Get first three widths from each plot
  g3max.widths <- do.call(unit.pmax, g3.widths) # Get maximum widths for first three widths across the plots
  for(i in 1:length(plots)) g[[i]]$widths[1:3] = g3max.widths # Apply the maximum widths to each plot
  return(g)
}

## 4 graphs for simplified Fig.2
plots <- list(y, sw, t, vw) # in order that want that arranged (filled by row)
g <- equalize(plots)
pdf("CgSx4 Fig2.pdf", width = 6.3, height= 1.5*6.30/4)
do.call(grid.arrange, c(g, ncol = 4)) # Draw it
dev.off()



dpi <- 300
tiff("F2_1r x 4c_width160mm.tif", width = 6.30*dpi, height= 1.5*6.30/4*dpi, res = dpi)
do.call(grid.arrange, c(g, ncol = 4)) # Draw it
dev.off()


