# univariate scatter plot w/jitter (each plant), raw concentration/activity data

#### soluble metabolite data --------------------------------------------------------
# N = 4 for all soluble metabolites x 2 rice lines (WT + CgSx4) x 4 harvests (3 leaf + 1 seed)
# ug/mg FW ## ions units 
# nmol/mg FW ## thiols units
# nmol/g FW ## amino acids units

df123s.sol <- read.delim("CgSx4WT_raw_soluble metabolites.txt", header = TRUE)
df123s.sol$line<-factor(df123s.sol$line, levels = c("WT", "CgSx4"))

df123.sol <- subset(df123s.sol, df123s.sol$tissue == "leaf")
dfs.sol <- subset(df123s.sol, df123s.sol$tissue == "seed")

#### protein-incorporated amino acid data ------------------------------------------
# N varies!
# N = 6 for WT @ 1, 2, 3 months
# N = 3 for CgSx4 in leaves EXCEPT @ 2months N = 2
# N = 3 in seeds
#umol/g DW!! ## protein-incorporated amino acid units

df123s.prot <- read.delim("CgSx4WT_raw_protein incorporated aa.txt", header = TRUE)
df123s.prot$line<-factor(df123s.prot$line, levels = c("WT", "CgSx4"))

df123.prot <- subset(df123s.prot, df123s.prot$tissue == "leaf")
dfs.prot <- subset(df123s.prot, df123s.prot$tissue == "seed") 


# functions to create ggplot objects -------------------
library(ggplot2)
library(grid)


# function leaf.univscatter() to create univariate scatter plot OBJECT
# for LEAF soluble and protein-incorporated DATA
# no legend

# df = df
# var = df$var
# metN = "", for graph title
# type = "", NULL, for graph title
# units = "", for y axis title
# jw = jitter width
# jh = jitter height
# dw = dodge width
# return is ggplot OBJECT

leaf.univscatter <- function (df, var, metN, type, units, jw, jh, dw){
  pjd <- position_jitterdodge(jitter.width = jw, 
                              jitter.height = jh, 
                              dodge.width = dw)
  title <- paste(type, metN, "in leaves")
  yLabel <- units
  
  usp <- ggplot(df, aes(x = harvest, y = var)) +
    geom_point(aes(fill = line), shape = 21, size = 1.75, stroke = 0.25, alpha = 0.9, position = pjd) +
    scale_fill_manual(values=c("black", "grey")) +
    xlab("plant age") +
    ylab(yLabel) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 8)) +
    theme(axis.title.x=element_text(size=7)) +
    theme(axis.title.y=element_text(size=7)) +
    theme(axis.text.x=element_text(color="black", size=6)) +
    theme(axis.text.y=element_text(color="black", size=6)) +
    theme(panel.background = element_rect(fill = "white", color = "black")) +
    theme(panel.border = element_rect(colour="black", fill=NA)) +
    scale_y_continuous(limits = c(0, NA)) +
    theme(legend.position="none")
  return(usp)
}



# function seed.univscatter() to create univariate scatter plot OBJECT
# for SEED soluble and protein-incorporated DATA
# no legend

# df = df
# var = df$var
# metN = "", for graph title
# type = "", NULL, for graph title
# units = "", for y axis title
# jw = jitter width
# jh = jitter height
# dw = dodge width
# return is ggplot OBJECT

seed.univscatter <- function (df, var, metN, type, units, jw, jh){
  pj <- position_jitter(width = jw, 
                        height = jh)
  title <- paste(type, metN, "in seeds")
  yLabel <- units
  
  usp <- ggplot(df, aes(x = line, y = var)) +
    geom_point(aes(fill = line), shape = 23, size = 1.75, stroke = 0.25, alpha = 0.9, position = pj) +
    scale_fill_manual(values=c("black", "grey")) +
    xlab("") +
    ylab(yLabel) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 8)) +
    theme(axis.title.x=element_text(size=7)) +
    theme(axis.title.y=element_text(size=7)) +
    theme(axis.text.x=element_text(color="black", size=6,angle = 90, hjust = 0.8, vjust = 0.5)) +
    theme(axis.text.y=element_text(color="black", size=6)) +
    theme(panel.background = element_rect(fill = "white", color = "black")) +
    theme(panel.border = element_rect(colour="black", fill=NA)) +
    scale_y_continuous(limits = c(0, NA)) +
    theme(legend.position="none")
  return(usp)
}



## ggplot objects ------

## for simplified Fig. 1
Met.leaves.sol <- leaf.univscatter(df123.sol, df123.sol$Met, "Met", "soluble", "nmol g-1 FW", 0.2, 0.05, 0.6)
Met.leaves.prot <- leaf.univscatter(df123.prot, df123.prot$Met, "Met", "protein", "µmol g-1 DW", 0.2, 0.05, 0.6)
Met.seeds.sol <- seed.univscatter(dfs.sol, dfs.sol$Met, "Met", "soluble", "nmol g-1 FW", 0.05, 0.05)
Met.seeds.prot <- seed.univscatter(dfs.prot, dfs.prot$Met, "Met", "protein", "µmol g-1 DW", 0.05, 0.05)


## for simplified Fig. 4
Arg.leaves.sol <- leaf.univscatter(df123.sol, df123.sol$Arg, "Arg", "soluble", "nmol g-1 FW", 0.2, 0.05, 0.6)
Gln.leaves.sol <- leaf.univscatter(df123.sol, df123.sol$Gln, "Gln", "soluble", "nmol g-1 FW", 0.2, 0.05, 0.6)
His.leaves.sol <- leaf.univscatter(df123.sol, df123.sol$His, "His", "soluble", "nmol g-1 FW", 0.2, 0.05, 0.6)
Hser.leaves.sol <- leaf.univscatter(df123.sol, df123.sol$Hser, "Hser", "soluble", "nmol g-1 FW", 0.2, 0.05, 0.6)
Leu.leaves.sol <- leaf.univscatter(df123.sol, df123.sol$Leu, "Leu", "soluble", "nmol g-1 FW", 0.2, 0.05, 0.6)
Lys.leaves.sol <- leaf.univscatter(df123.sol, df123.sol$Lys, "Lys", "soluble", "nmol g-1 FW", 0.2, 0.05, 0.6)
SO4.leaves <- leaf.univscatter(df123.sol, df123.sol$SO4, "SO4", NULL, "µg mg-1 FW", 0.2, 0.05, 0.6)
Cys.leaves.sol <- leaf.univscatter(df123.sol, df123.sol$Cys, "Cys", "soluble", "nmol mg-1 FW", 0.2, 0.05, 0.6)
GSH.leaves.sol <- leaf.univscatter(df123.sol, df123.sol$GSH, "GSH", "soluble", "nmol mg-1 FW", 0.2, 0.05, 0.6)
Hcys.leaves.sol <- leaf.univscatter(df123.sol, df123.sol$Hcys, "Hcys", "soluble", "nmol mg-1 FW", 0.2, 0.05, 0.6)



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


# assemble panels save PDF ------------

## 2 leaf graphs + 2 seed graphs (simplified CGS Fig.1)
plots <- list(Met.leaves.sol, Met.leaves.prot, Met.seeds.sol, Met.seeds.prot) # in order that want that arranged (filled by row)
g <- equalize(plots)
pdf("CgSx4 Fig1.pdf", width = 3.54, height= 1.5*2*3.54/2)
do.call(grid.arrange, c(g, ncol = 2)) # Draw it
dev.off()

## 10 leaf graphs (simplified CGS Fig.4)
plots <- list(Arg.leaves.sol, Cys.leaves.sol, Gln.leaves.sol, GSH.leaves.sol, Hcys.leaves.sol, His.leaves.sol, Hser.leaves.sol, Leu.leaves.sol, Lys.leaves.sol, SO4.leaves) # in order that want that arranged (filled by row)
g <- equalize(plots)
pdf("CgSx4 Fig4.pdf", width = 7.4, height= 3*1.25*7.4/4)
do.call(grid.arrange, c(g, ncol = 4)) # Draw it
dev.off()