#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# This code calculates BC 
#--------------------------------------------------------------------------------------

#This is the final step.

#1. DownloadPbdb.R
#2. Dataprep.R
#3. Hex_sim.R
#4  Survival_PT.R
#5. BC.R

library(plyr)
library(dplyr)
library(deeptime)
library(ggplot2)

# Create function to calculate Biogeographic Connectedness.
#--------------------------------------------------------------------------------------
BC <- function(genera, localities) {
  O <- length(genera)
  N <- length(unique(genera))
  L <- length(unique(localities))
  return((O-N)/((L*N)-N))
}

table(Sim.c$x.cell)
table(Sim.c$y.cell)#find cells used so we can subset to those cells
table(Sim.i$y.cell)
table(Sim.i$x.cell)

cells20 <- as.data.frame(c(125, 135, 144, 246, 247, 251, 252, 253, 265, 381, 594, 601, 
             602, 614, 652, 653, 661, 662, 663, 667, 672, 673, 693))  
names(cells20) <- "cell"
                           
data.bc <- merge(data.j, cells20, by = "cell", all =FALSE) #dataset for cells with 20 or more occs
data.bc <- unique(subset(data.bc, select = -c(X, lat, lng))) #get unique genus occurrences per cell, remove unnecessary columns

#By survival status:
bc.all <-ddply(data.bc, c("interval.ma", "extinct"), function(df)c(BC(df$accepted_name, df$cell),
                                                                   length(df$accepted_name),
                                                                   length(unique(df$accepted_name))))
names(bc.all) <- c("interval.ma", "extinct", "BC", "occurrences", "genera")
bc.all$labels <- c("Survivors - After", "Survivors - Before", "Victims")

#General P-T:
bc.pt <-ddply(data.bc, c("interval.ma"), function(df)c(BC(df$accepted_name, df$cell),
                                                                   length(df$accepted_name),
                                                                   length(unique(df$accepted_name))))
names(bc.pt) <- c("interval.ma", "BC", "occurrences", "genera")
bc.pt$labels <-c("Changhsingian", "Induan")


#using all possible cells, and not just ones with 20+ occ's:
bc.allcells<-ddply(data.j, c("interval.ma", "extinct"), function(df)c(BC(df$accepted_name, df$cell),
                                                                   length(df$accepted_name),
                                                                   length(unique(df$accepted_name))))
names(bc.allcells) <- c("interval.ma", "extinct", "BC", "occurrences", "genera")
bc.allcells$labels <- c("Survivors - After", "Survivors - Before", "Victims")

#all cells + general P-T:
bc.ptall <-ddply(data.j, c("interval.ma"), function(df)c(BC(df$accepted_name, df$cell),
                                                       length(df$accepted_name),
                                                       length(unique(df$accepted_name))))
names(bc.ptall) <- c("interval.ma", "BC", "occurrences", "genera")
bc.pt$labels <-c("Changhsingian", "Induan")


table(data.j$interval.ma)
table(data.bc$interval.ma)
