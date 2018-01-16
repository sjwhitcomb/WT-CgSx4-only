# univariate scatter plot w/jitter (each plant), raw concentration/activity data

#### soluble metabolite data -----------------
# N = 4 for all soluble metabolites x 2 rice lines (WT + CgSx4) x 4 harvests (3 leaf + 1 seed)
# ug/mg FW ## ions units 
# nmol/mg FW ## thiols units
# nmol/g FW ## amino acids units

df123s.sol <- read.delim("CgSx4WT_raw_soluble metabolites.txt", header = TRUE)
df123s.sol$line<-factor(df123s.sol$line, levels = c("WT", "CgSx4"))
df123.sol <- subset(df123s.sol, df123s.sol$tissue == "leaf")
dfs.sol <- subset(df123s.sol, df123s.sol$tissue == "seed") # NOTE: no seed data for CGS, SMM, Hse, NO3, Hcys

#### protein-incorporated amino acid data -----------------
# N varies!
# N = 6 for WT @ 1, 2, 3 months
# N = 3 for CgSx4 in leaves EXCEPT @ 2month CgSx4 = 2
# N = 3 for WT and CgSx4
# umol/g DW!! ## protein-incorporated amino acid units

df123s.prot <- read.delim("CgSx4WT_raw_protein incorporated aa.txt", header = TRUE)
df123s.prot$line<-factor(df123s.prot$line, levels = c("WT", "CgSx4"))
df123.prot <- subset(df123s.prot, df123s.prot$tissue == "leaf")
dfs.prot <- subset(df123s.prot, df123s.prot$tissue == "seed") 

# packages ------
library(ggplot2)
library(grid)


#### iterate over all soluble LEAF metabolites for supplemental figure ----------
# iterate -> named! list of ggplot objects AND tiff files

# ions and thiols and amino acids
output <- vector(mode = "list", length = length(df123.sol))
names(output) <- names(df123.sol)

for(i in c(9:32)) {
  metabolite <- names(df123.sol[i])
  saveName <-paste("univariate scatter", metabolite, "soluble", "leaves" ,"noLegend.tiff", sep = "_")
  
  pjd <- position_jitterdodge(jitter.width = 0.2, 
                              jitter.height = 0.05, 
                              dodge.width = 0.6)
  
  if(i >= 28 & i <= 29) {
    yLabel <- "µg mg-1 FW" 
    title <- paste(metabolite, "in leaves")
    
    # ions, no legend
    sp <- ggplot(df123.sol, aes_string(x = "harvest", y = names(df123.sol)[i])) +
      geom_point(aes(fill = line), shape = 21, size = 1.75, stroke = 0.25, alpha = 0.9, position = pjd) +
      scale_fill_manual(values=c("black", "grey")) +
      xlab("plant age") +
      ylab(yLabel) +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size = 8)) +
      theme(axis.title = element_text(size = 7)) +
      theme(axis.text.x=element_text(color="black", size=6)) +
      theme(axis.text.y=element_text(color="black", size=6)) +
      theme(panel.background = element_rect(fill = "white", color = "black")) +
      theme(panel.border = element_rect(colour="black", fill=NA)) +
      scale_y_continuous(limits = c(0, NA)) +
      theme(legend.position="none")
    
    output[[i]] <- sp
    ggsave(saveName, width=1.6, height=2, unit="in", dpi=300)
    
  } else if (i >= 30 & i <= 32){
    yLabel <- "µg mg-1 FW" 
    title <- paste("soluble", metabolite, "in leaves")
    
    # thiols, no legend
    sp <- ggplot(df123.sol, aes_string(x = "harvest", y = names(df123.sol)[i])) +
      geom_point(aes(fill = line), shape = 21, size = 1.75, stroke = 0.25, alpha = 0.9, position = pjd) +
      scale_fill_manual(values=c("black", "grey")) +
      xlab("plant age") +
      ylab(yLabel) +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size = 8)) +
      theme(axis.title = element_text(size = 7)) +
      theme(axis.text.x=element_text(color="black", size=6)) +
      theme(axis.text.y=element_text(color="black", size=6)) +
      theme(panel.background = element_rect(fill = "white", color = "black")) +
      theme(panel.border = element_rect(colour="black", fill=NA)) +
      scale_y_continuous(limits = c(0, NA)) +
      theme(legend.position="none")
    
    output[[i]] <- sp
    ggsave(saveName, width=1.6, height=2, unit="in", dpi=300)
    
  } else if (i >= 9 & i <= 27){
    yLabel <- "nmol g-1 FW"
    title <- paste("soluble", metabolite, "in leaves")
    
    # amino acids, no legend
    sp <- ggplot(df123.sol, aes_string(x = "harvest", y = names(df123.sol)[i])) +
      geom_point(aes(fill = line), shape = 21, size = 1.75, stroke = 0.25, alpha = 0.9, position = pjd) +
      scale_fill_manual(values=c("black", "grey")) +
      xlab("plant age") +
      ylab(yLabel) +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size = 8)) +
      theme(axis.title = element_text(size = 7)) +
      theme(axis.text.x=element_text(color="black", size=6)) +
      theme(axis.text.y=element_text(color="black", size=6)) +
      theme(panel.background = element_rect(fill = "white", color = "black")) +
      theme(panel.border = element_rect(colour="black", fill=NA)) +
      scale_y_continuous(limits = c(0, NA)) +
      theme(legend.position="none")
    
    output[[i]] <- sp
    ggsave(saveName, width= 1.6, height=2, unit="in", dpi=300)
  } else
    print("indexing not working")
}

objects <- output[c(9:32)] ## remove empty elements of the output list
alph <- objects[order(names(objects))] ## alphabetize the elements

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

## for simplified Supplemental fig 1
plots <- alph # in alphabetical order, note, grid will be filled by row
g <- equalize(plots)
pdf("CgSx4 SupFig1.pdf", width = 1.6*6, height= 2*4)
do.call(grid.arrange, c(g, ncol = 6)) # Draw it
dev.off()


#### iterate over all soluble SEED metabolites for supplemental figure ----------
# iterate > named! list of ggplot objects AND tiff files

# ions and thiols and amino acids
output <- vector(mode = "list", length = length(dfs.sol))
names(output) <- names(dfs.sol)

for(i in c(9:32)) {
  metabolite <- names(dfs.sol[i])
  saveName <-paste("univariate scatter", metabolite, "soluble", "seeds" ,"noLegend.tiff", sep = "_")
  
  pj <- position_jitter(width = 0.05, 
                        height = 0.05) 
  
  if(i >= 28 & i <= 29) {
    yLabel <- "µg mg-1 FW" 
    title <- paste(metabolite, "in seeds")
    
    # ions, no legend
    sp <- ggplot(dfs.sol, aes_string(x = "line", y = names(dfs.sol)[i])) +
      geom_point(aes(fill = line), shape = 23, size = 1.75, stroke = 0.25, alpha = 0.9, position = pj) +
      scale_fill_manual(values=c("black", "grey")) +
      xlab("") +
      ylab(yLabel) +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size = 8)) +
      theme(axis.title = element_text(size = 6)) +
      theme(axis.text.x=element_text(color="black", size=6,angle = 90, hjust = 0.8, vjust = 0.5)) +
      theme(axis.text.y=element_text(color="black", size=6)) +
      theme(panel.background = element_rect(fill = "white", color = "black")) +
      theme(panel.border = element_rect(colour="black", fill=NA)) +
      scale_y_continuous(limits = c(0, NA)) +
      theme(legend.position="none")
    
    output[[i]] <- sp
    ggsave(saveName, width=1.6, height=2, unit="in", dpi=300)
    
  } else if (i >= 30 & i <= 32){
    yLabel <- "µg mg-1 FW" 
    title <- paste("soluble", metabolite, "in seeds")
    
    # thiols, no legend
    sp <- ggplot(dfs.sol, aes_string(x = "line", y = names(dfs.sol)[i])) +
      geom_point(aes(fill = line), shape = 23, size = 1.75, stroke = 0.25, alpha = 0.9, position = pj) +
      scale_fill_manual(values=c("black", "grey")) +
      xlab("") +
      ylab(yLabel) +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size = 8)) +
      theme(axis.title = element_text(size = 6)) +
      theme(axis.text.x=element_text(color="black", size=6,angle = 90, hjust = 0.8, vjust = 0.5)) +
      theme(axis.text.y=element_text(color="black", size=6)) +
      theme(panel.background = element_rect(fill = "white", color = "black")) +
      theme(panel.border = element_rect(colour="black", fill=NA)) +
      scale_y_continuous(limits = c(0, NA)) +
      theme(legend.position="none")
    
    output[[i]] <- sp
    ggsave(saveName, width=1.6, height=2, unit="in", dpi=300)
    
  } else if (i >= 9 & i <= 27){
    yLabel <- "nmol g-1 FW"
    title <- paste("soluble", metabolite, "in seeds")
    
    # amino acids, no legend
    sp <- ggplot(dfs.sol, aes_string(x = "line", y = names(dfs.sol)[i])) +
      geom_point(aes(fill = line), shape = 23, size = 1.75, stroke = 0.25, alpha = 0.9, position = pj) +
      scale_fill_manual(values=c("black", "grey")) +
      xlab("") +
      ylab(yLabel) +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size = 8)) +
      theme(axis.title = element_text(size = 6)) +
      theme(axis.text.x=element_text(color="black", size=6,angle = 90, hjust = 0.8, vjust = 0.5)) +
      theme(axis.text.y=element_text(color="black", size=6)) +
      theme(panel.background = element_rect(fill = "white", color = "black")) +
      theme(panel.border = element_rect(colour="black", fill=NA)) +
      scale_y_continuous(limits = c(0, NA)) +
      theme(legend.position="none")
    
    output[[i]] <- sp
    ggsave(saveName, width=1.6, height=2, unit="in", dpi=300)
  } else
    print("indexing not working")
}

objects <- output[c(9:17, 19:28, 30:31)] ## remove empty elements of the output list ## there are also empty graphs in seeds!
alph <- objects[order(names(objects))] ## alphabetize the elements

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

## for simplified Supplemental figure 4 - 21 graphs
plots <- alph # in alphabetical order, note, grid will be filled by row
g <- equalize(plots)
pdf("CgSx4 SupFig4.pdf", width = 1.6*6, height= 2*4)
do.call(grid.arrange, c(g, ncol = 6)) # Draw it
dev.off()


#### iterate over all protein-incorporated amino acids in SEED for supplemental figure ----------
# iterate > named! list of ggplot objects AND tiff files

# amino acids
output <- vector(mode = "list", length = length(dfs.prot))
names(output) <- names(dfs.prot)

for(i in c(6:19)) {
  metabolite <- names(dfs.prot[i])
  saveName <-paste("univariate scatter", metabolite, "prot", "seeds" ,"noLegend.tiff", sep = "_")
  yLabel <- "µmol g-1 DW"
  title <- paste("protein", metabolite, "in seeds")
  pj <- position_jitter(width = 0.05, 
                        height = 0.05)
  
  # no legend
  sp <- ggplot(dfs.prot, aes_string(x = "line", y = names(dfs.prot)[i])) +
    geom_point(aes(fill = line), shape = 23, size = 1.75, stroke = 0.25, alpha = 0.9, position = pj) +
    scale_fill_manual(values=c("black", "grey")) +
    xlab("") +
    ylab(yLabel) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 8)) +
    theme(axis.title = element_text(size = 6)) +
    theme(axis.text.x=element_text(color="black", size=6,angle = 90, hjust = 0.8, vjust = 0.5)) +
    theme(axis.text.y=element_text(color="black", size=6)) +
    theme(panel.background = element_rect(fill = "white", color = "black")) +
    theme(panel.border = element_rect(colour="black", fill=NA)) +
    scale_y_continuous(limits = c(0, NA)) +
    theme(legend.position="none")
  
  output[[i]] <- sp
  ggsave(saveName, width=1.6, height=2, unit="in", dpi=300)
}

objects <- output[c(6:19)] ## remove empty elements of the output list
alph <- objects[order(names(objects))] ## alphabetize the elements

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

## for simplified supplemental 5
plots <- alph # in alphabetical order, note, grid will be filled by row
g <- equalize(plots)
pdf("CgSx4WT SupFig5.pdf", width = 1.6*5, height= 2*3)
do.call(grid.arrange, c(g, ncol = 5)) # Draw it
dev.off()
