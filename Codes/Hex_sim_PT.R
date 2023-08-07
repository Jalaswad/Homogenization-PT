#--------------------------------------------------------------------------------------
# This code analyzes similarity of taxa
#--------------------------------------------------------------------------------------

#This is step 3.

#1. DownloadPbdb.R
#2. Dataprep.R
#3. Hex_sim.R
#4  Survival_PT.R
#5. BC.R

#setup
#--------------------------------------------------------------------------------------
setwd("")

library(plyr)
library(dplyr)
library(tidyverse)
library(fossil)
library(gridExtra)

# For each cell, get a list of genera and numbers of occurrences
#--------------------------------------------------------------------------------------
Genus_occs_pt <-ddply(data.j, c("early_interval", "cell", "accepted_name"),
                       function(df)c(length(df$accepted_name)))
names(Genus_occs_pt) <- c("early_interval", "cell", "accepted_name", "occs")

# Get coordinates of cell centers to compare distances 
#---------------------------------------------------------------------------------
cells.pt <-(sort(unique(Genus_occs_pt$cell))) #all unique cells in the data frame
cell.pt.matrix <- data.frame(merge(cells.pt, cells.pt)) #create two columns with all pairs of values
cell.pt.matrix <-subset(cell.pt.matrix, x <y) #prevent duplicates and self-comparisons

split.genusoccs   <- split(Genus_occs_pt, Genus_occs_pt$early_interval) #split for each interval
list2env(split.genusoccs,envir=.GlobalEnv) #save as separate dataframes

#Calculate Great Circle Distance:
#---------------------------------------------------------------------------------
gcd.slc <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  long1.r <- long1*pi/180 #convert from degrees to radians for lats and longs
  long2.r <- long2*pi/180
  lat1.r <- lat1*pi/180
  lat2.r <- lat2*pi/180
  d <- acos(sin(lat1.r)*sin(lat2.r) + cos(lat1.r)*cos(lat2.r) * cos(long2.r-long1.r)) * R
  return(d) }# Distance in km

#---------------------------------------------------------------------------------
#Create dataframe with analysis for cells, coordinates, coef's, distances, # genera, etc'
#---------------------------------------------------------------------------------

#Changhsingian:
#---------------------------------------------------------------------------------

Sim.c             <- data.frame(matrix(ncol=13, nrow=length(cell.pt.matrix$x)))

names(Sim.c)      <- c("x.cell", "y.cell","x.lon", "x.lat", "y.lon", 
                       "y.lat", "Czek", "Jacc", "Dist", "Occs.x","Occs.y", "Gen.x", "Gen.y")

for (i in 1:length(cell.pt.matrix$x)) {
  
  Sim.c[i,1]     <- cell.pt.matrix$x[i] #record x cell number
  Sim.c[i,2]     <- cell.pt.matrix$y[i] #record y cell number
  
  j <- cell.pt.matrix$x[i]
  k <- cell.pt.matrix$y[i]
  
  Sim.c[i,3]     <- cells1$long[cells1$cells==j] 
  Sim.c[i,4]     <- cells1$lat[cells1$cells==j]
  Sim.c[i,5]     <- cells1$long[cells1$cells==k]
  Sim.c[i,6]     <- cells1$lat[cells1$cells==k] 
  
  data.1        <- subset(Changhsingian, cell==cell.pt.matrix$x[i]) 
  data.2        <- subset(Changhsingian, cell==cell.pt.matrix$y[i])
  data.m        <- merge(data.1, data.2, by=c("accepted_name"), all=TRUE) 
  data.m[is.na(data.m)] = 0 
  data.m$min    <- pmin(data.m$occs.x, data.m$occs.y)
  
  
  Sim.c[i,7]      <-2*abs(sum(data.m$min))/((sum(data.m$occs.x)+sum(data.m$occs.y))) #Czek
  Sim.c[i,8]     <- (length(data.m$accepted_name[data.m$occs.x>0 & data.m$occs.y>0])/ 
                       (length(data.m$accepted_name))) #calculate Jaccard coefficient
  Sim.c[i,9]      <- (gcd.slc(Sim.c[i,3], Sim.c[i,4], Sim.c[i,5], Sim.c[i,6]))
  
  Sim.c[i,10]    <- sum(data.m$occs.x)#sum of all occs in that cell for x
  Sim.c[i,11]    <- sum(data.m$occs.y) #sum for all in cell y
  Sim.c[i,12]    <- length(data.1$accepted_name) # num of genera in x
  Sim.c[i,13]    <- length(data.2$accepted_name) # num of genera in y 
  
  print(i) }

Sim.c$CutDist <- cut(Sim.c$Dist, breaks = c(0, 2000, 4000, 6000, 8000, 10000, 12000, 
                                            14000, 16000, 18000, 20000), 
                     labels = c("0", "2000", "4000", "6000",
                                "8000", "10000", "12000","14000", 
                                "16000", "18000"), include.lowest = TRUE) 

#Induan:
#---------------------------------------------------------------------------------

Sim.i             <- data.frame(matrix(ncol=13, nrow=length(cell.pt.matrix$x)))

names(Sim.i)      <- c("x.cell", "y.cell","x.lon", "x.lat", "y.lon", 
                       "y.lat", "Czek", "Jacc", "Dist",
                       "Occs.x", "Occs.y", "Gen.x", "Gen.y")

for (i in 1:length(cell.pt.matrix$x)) {
  #for (i in 1:100) {
  Sim.i[i,1]     <- cell.pt.matrix$x[i] #record x cell number
  Sim.i[i,2]     <- cell.pt.matrix$y[i] #record y cell number
  
  j <- cell.pt.matrix$x[i]
  k <- cell.pt.matrix$y[i]
  
  Sim.i[i,3]     <- cells1$long[cells1$cell==j]
  Sim.i[i,4]     <- cells1$lat[cells1$cell==j]
  Sim.i[i,5]     <- cells1$long[cells1$cell==k]
  Sim.i[i,6]     <- cells1$lat[cells1$cell==k]
  
  data.3        <- subset(Induan, cell==cell.pt.matrix$x[i]) 
  data.4        <- subset(Induan, cell==cell.pt.matrix$y[i])
  data.mi        <- merge(data.3, data.4, by=c("accepted_name"), all=TRUE) 
  data.mi[is.na(data.mi)] = 0 
  data.mi$mini    <- pmin(data.mi$occs.x, data.mi$occs.y)
  
  Sim.i[i,7]      <-2*abs(sum(data.mi$mini))/ ((sum(data.mi$occs.x)+sum(data.mi$occs.y))) 
  Sim.i[i,8]     <- (length(data.mi$accepted_name[data.mi$occs.x>0 & data.mi$occs.y>0])/ 
                       (length(data.mi$accepted_name))) #calculate Jaccard coefficient
  Sim.i[i,9]      <- (gcd.slc(Sim.i[i,3], Sim.i[i,4], Sim.i[i,5], Sim.i[i,6]))
  
  Sim.i[i,10]    <- sum(data.mi$occs.x)#sum of all occs in that cell for x
  Sim.i[i,11]    <- sum(data.mi$occs.y) #sum for all in cell y
  Sim.i[i,12]    <- length(data.3$accepted_name) # num of genera in x
  # length(data.m$genus[data.m$occs.x >0]) 
  Sim.i[i,13]    <- length(data.4$accepted_name) # num of genera in y 
  
  print(i) 
}


Sim.i$CutDist <- cut(Sim.i$Dist, breaks = c(0, 2000, 4000, 6000, 8000, 10000, 12000, 
                                            14000, 16000, 18000, 20000), 
                     labels = c("0", "2000", "4000", "6000",
                                "8000", "10000", "12000","14000", 
                                "16000", "18000"), include.lowest = TRUE) 

#Subset for minimium no. of genera to compare:
#---------------------------------------------------------------------------------
Sim.c <- subset(Sim.c, Occs.x >20 & Occs.y >20)
Sim.i<- subset(Sim.i, Occs.x >20 & Occs.y >120)

Sim.cr <- subset(Sim.c, Occs.x >20 & Occs.y >20)
Sim.ir<- subset(Sim.i, Occs.x >20 & Occs.y >20)

#Find the average and sd
#--------------------------------------------------------------------------------
#Changhsingian:
Sim.c.na <- Sim.c %>% drop_na()
Sim.i.na <- Sim.i %>% drop_na()

#detach(package:plyr) #sometimes required to get sd to work

Sum.Cz<- Sim.c.na %>% 
  group_by(CutDist) %>% 
  summarize(avg =mean(Czek), sdev = sd(Czek))
Sum.Cz$label <- 'Changhsingian'
Sum.Cz <-as.data.frame(Sum.Cz)

Sum.Ja <- Sim.c.na%>% 
  group_by(CutDist) %>% 
  summarize(avg =mean(Jacc),sdev = sd(Jacc))
Sum.Ja$label <- 'Changhsingian'
Sum.Ja <-as.data.frame(Sum.Ja)

table(Sim.c$CutDist)
table(Sim.i$CutDist)
Sum.Cz$bubbles <- c(7,4,8,7,9,4,1,3,2) 
Sum.Ja$bubbles <- c(7,4,8,7,9,4,1,3,2) 
#bubbles = number of cells compared. found in Sim dataframes using table() above

#Induan:
Sum.i.Cz <- Sim.i.na %>% 
  group_by(CutDist) %>% 
  summarize(avg =mean(Czek),sdev = sd(Czek))
Sum.i.Cz$label <- 'Induan'
Sum.i.Cz <-as.data.frame(Sum.i.Cz)

Sum.i.Ja <-Sim.i.na %>%
  group_by(CutDist) %>%
  summarize(avg =mean(Jacc),sdev = sd(Jacc))
Sum.i.Ja$label <- 'Induan'
Sum.i.Ja <-as.data.frame(Sum.i.Ja)

Sum.i.Cz$bubbles <- c(10,2,14,5,8,17,9,1) 
Sum.i.Ja$bubbles <- c(10,2,14,5,8,17,9,1)

#Concatenate Cha and Ind:
Sum.Ja.2 <- rbind(Sum.Ja, Sum.i.Ja)
Sum.Cz.2 <-rbind(Sum.Cz, Sum.i.Cz)

#plots:
#---------------------------------------------------------------------------------
library(ggplot2)

yb =c("deepskyblue","dodgerblue4" )

Sum.Cz.2$CutDist <- as.numeric(Sum.Cz.2$CutDist)
plot1.g <- Sum.Cz.2 %>%
  ggplot( aes(x=CutDist, y=avg, group=as.factor(label),
              colour=as.factor(label),size = bubbles^2))+
  scale_colour_manual(name="Age", values=yb)+
  theme_classic()+
  geom_point(alpha=0.8)+
  scale_size_area()+
  theme(text = element_text(size=25),
        axis.text.x = element_text(angle=0, hjust=0.5, size=25),
        axis.text.y = element_text(size=25)) +
  ylim(0,0.65)+
  scale_x_continuous(limits=c(1,9), breaks=c(1, 2, 3, 4, 5, 6,
                                                  7, 8, 9),
                     labels= as.character(c("2000", "4000", "6000",
                                            "8000", "10000", "12000",
                                            "14000", "16000", "18000")))+
  geom_line(size=2)+
  geom_errorbar(aes(ymin=pmax(avg-sdev,0), ymax=avg+sdev), width=0.05,size=1, alpha=0.7) + 
  theme(plot.margin = unit(c( t= 0.25, r=0.5, b=0.75,l=.25), "cm"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        plot.title = element_text(size=20), legend.position= c(.8, .9))+ theme(legend.title = element_blank())+
  guides(size = "none")+ #remove legend for bubbles
  labs (title = "Czekanowski")
plot1.g

Sum.Ja.2$CutDist <- as.numeric(Sum.Ja.2$CutDist)
plot2.g <- Sum.Ja.2 %>%
  ggplot( aes(x=CutDist, y=avg, group=as.factor(label),
              colour=as.factor(label),size = bubbles^2))+
  scale_colour_manual(name="Age", values=yb)+
  theme_classic()+
  geom_point(alpha=0.8)+
  scale_size_area()+
  theme(text = element_text(size=25),
        axis.text.x = element_text(angle=0, hjust=0.5, size=25),
        axis.text.y = element_text(size=25)) +
  ylim(0,0.65)+
  scale_x_continuous(limits=c(1,9), breaks=c(1, 2, 3, 4, 5, 6,
                                             7, 8, 9),
                     labels= as.character(c("2000", "4000", "6000",
                                            "8000", "10000", "12000",
                                            "14000", "16000", "18000")))+
  geom_line(size=2)+
  geom_errorbar(aes(ymin=pmax(avg-sdev,0), ymax=avg+sdev), width=0.05,size=1, alpha=0.7) + 
  theme(plot.margin = unit(c( t= 0.25, r=0.5, b=0.75,l=.25), "cm"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        plot.title = element_text(size=20), legend.position= "none")+
  guides(size = "none")+ #remove legend for bubbles
  labs (title = "Jaccard")
plot2.g

#Below, we will arrange both plots 1 and 2 side-by-side in one image so we can see
#how the Czekanowski and Jaccard methods differ:
#_______________________________________________________________________________
library(gridExtra) #allows you to arrange plots side by side in a grid
library(grid)

bg1<-grid.arrange(plot1.g,plot2.g, nrow=2, 
                    bottom = textGrob("", vjust = 0,  
                                      gp = gpar(cex = 0.95)),
                    left = textGrob("",  
                                    rot = 0, vjust = 0, 
                                    gp = gpar( cex = 0.95)),
                    top = textGrob("", vjust =0, 
                                   gp = gpar(fontface = "bold", fontsize =16, cex = 0.95))) 
bg1


#Below, we create the grid and save it directly to our working directory:
#___________________________________________________________________________________________________________
ggsave("bg.pt.occs20.pdf", plot = grid.arrange(plot1.g,plot2.g, nrow=2, 
                                           bottom = textGrob("Great Circle Distance", vjust = -1.5, 
                                                             gp = gpar(cex = 0.95)),
                                           left = textGrob("", 
                                                           rot = 90, vjust = 4,
                                                           gp = gpar( cex = 0.95)),
                                           top = textGrob("", vjust =0.5,
                                                          gp = gpar(fontface = "bold", fontsize =40, cex = 0.95))),
       height=11, width=18, device="pdf", dpi=500) 
