#--------------------------------------------------------------------------------------
# This code will separate genera by survival status: Survived, Exctinct, Originated and
# then run separate similarity analyses 
#--------------------------------------------------------------------------------------

#This is step 4.

#1. DownloadPbdb.R
#2. Dataprep.R
#3. Hex_sim
#4. Survival_PT.R 
#5. BC.R

library(plyr)
library(dplyr)
library(tidyverse)
library(fossil)
library(ggplot2)
library(grid)
library(gridExtra)

setwd("")

#Subset by survival status:
#--------------------------------------------------------------------------------------
#Survivors:
survivors <- subset(data.j, extinct ==0)
survivors$label <- 'S'
#Extinct:
extinct    <- subset(data.j, extinct=1)  
extinct$label <- 'E'


#Calculate Great Circle Distance:
#-------------------------------------------------------------------------------
gcd.slc <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  long1.r <- long1*pi/180 #convert from degrees to radians for lats and longs
  long2.r <- long2*pi/180
  lat1.r <- lat1*pi/180
  lat2.r <- lat2*pi/180
  d <- acos(sin(lat1.r)*sin(lat2.r) + cos(lat1.r)*cos(lat2.r) 
            * cos(long2.r-long1.r)) * R
  return(d) }# Distance in km

#Sim for each survival type:
#-------------------------------------------------------------------------------
 
#Survived
#------------------------------------------------------------------------------
Genus_occs_pt.s <-plyr::ddply(survivors, c("interval.ma", "cell", "accepted_name"),
                           function(df)c(length(df$accepted_name)))
names(Genus_occs_pt.s) <- c("interval.ma", "cell", "accepted_name", "occs")
Genus_occs_pt.s$stage <- NA
Genus_occs_pt.s$stage <- ifelse(Genus_occs_pt.s$interval.ma == 252.17, "Changhsingian", "Induan")

cells.pt.s <-(sort(unique(Genus_occs_pt.s$cell))) #all unique cells in the data frame
cell.pt.matrix.s <- data.frame(merge(cells.pt.s, cells.pt.s)) #create two columns with all pairs of values
cell.pt.matrix.s <-subset(cell.pt.matrix.s, x <y) #prevent duplicates and self-comparisons

split.genusoccs   <- split(Genus_occs_pt.s, Genus_occs_pt.s$stage)  #split for each interval
list2env(split.genusoccs,envir=.GlobalEnv) #save as separate dataframes

#Dataframe for similarity of Survivors:
Sim.s           <- data.frame(matrix(ncol=13, nrow=length(cell.pt.matrix.s$x)))

names(Sim.s)    <- c("x.cell", "y.cell","x.lon", "x.lat", "y.lon", 
                     "y.lat", "Czek", "Jacc", "Dist", "Occs.x",
                     "Occs.y", "Gen.x", "Gen.y")

for (i in 1:length(cell.pt.matrix.s$x)) {
  # for  (i in 1:100) {
  
  Sim.s[i,1]     <- cell.pt.matrix.s$x[i] #record x cell number
  Sim.s[i,2]     <- cell.pt.matrix.s$y[i] #record y cell number
  
  j <- cell.pt.matrix.s$x[i]
  k <- cell.pt.matrix.s$y[i]
  
  Sim.s[i,3]     <- cells1$long[cells1$cells==j] #record x cell lat
  Sim.s[i,4]     <- cells1$lat[cells1$cells==j] #record x cell lng
  Sim.s[i,5]     <- cells1$long[cells1$cells==k] #record y cell lat
  Sim.s[i,6]     <- cells1$lat[cells1$cells==k] #record y cell lng
  
  data.1        <- subset(Changhsingian, cell==cell.pt.matrix.s$x[i]) 
  data.2        <- subset(Changhsingian, cell==cell.pt.matrix.s$y[i])
  data.m        <- merge(data.1, data.2, by=c("accepted_name"), all=TRUE) 
  data.m[is.na(data.m)] = 0 
  data.m$min    <- pmin(data.m$occs.x, data.m$occs.y)
  
  Sim.s[i,7]      <-2*abs(sum(data.m$min))/((sum(data.m$occs.x)+sum(data.m$occs.y))) 
  Sim.s[i,8]     <- (length(data.m$accepted_name[data.m$occs.x>0 & data.m$occs.y>0])/ 
                       (length(data.m$accepted_name))) #calculate Jaccard coefficient
  Sim.s[i,9]      <- (gcd.slc(Sim.s[i,3], Sim.s[i,4], Sim.s[i,5], Sim.s[i,6]))
  
  
  Sim.s[i,10]    <- sum(data.m$occs.x)#sum of all occs in that cell for x
  Sim.s[i,11]    <- sum(data.m$occs.y) #sum for all in cell y
  Sim.s[i,12]    <- length(data.1$accepted_name) # num of genera in x
  Sim.s[i,13]    <- length(data.2$accepted_name) # num of genera in y 
  
  print(i) }

Sim.s$CutDist <- cut(Sim.s$Dist, breaks = c(0, 2000, 4000, 6000, 8000, 10000, 12000, 
                                            14000, 16000, 18000, 20000),
                     labels = c("0", "2000", "4000", "6000",
                                "8000", "10000", "12000","14000", 
                                "16000", "18000"), include.lowest = TRUE) 
Sim.s <- drop_na(Sim.s)

#Induan:
Sim.is             <- data.frame(matrix(ncol=13, nrow=length(cell.pt.matrix.s$x)))

names(Sim.is)      <- c("x.cell", "y.cell","x.lon", "x.lat", "y.lon", 
                        "y.lat", "Czek", "Jacc", "Dist", "Occs.x",
                        "Occs.y", "Gen.x", "Gen.y")

for (i in 1:length(cell.pt.matrix.s$x)) {
  # for  (i in 1:100) {
  
  Sim.is[i,1]     <- cell.pt.matrix.s$x[i] #record x cell number
  Sim.is[i,2]     <- cell.pt.matrix.s$y[i] #record y cell number
  
  j <- cell.pt.matrix.s$x[i]
  k <- cell.pt.matrix.s$y[i]
  
  Sim.is[i,3]     <- cells1$long[cells1$cells==j] #record x cell lat
  Sim.is[i,4]     <- cells1$lat[cells1$cells==j] #record x cell lng
  Sim.is[i,5]     <- cells1$long[cells1$cells==k] #record y cell lat
  Sim.is[i,6]     <- cells1$lat[cells1$cells==k] #record y cell lng
  
  data.11        <- subset(Induan, cell==cell.pt.matrix.s$x[i]) 
  data.22        <- subset(Induan, cell==cell.pt.matrix.s$y[i])
  data.mm        <- merge(data.11, data.22, by=c("accepted_name"), all=TRUE) 
  data.mm[is.na(data.mm)] = 0 
  data.mm$min    <- pmin(data.mm$occs.x, data.mm$occs.y)
  
  Sim.is[i,7]      <-2*abs(sum(data.mm$min))/((sum(data.mm$occs.x)+sum(data.mm$occs.y))) 
  Sim.is[i,8]     <- (length(data.mm$accepted_name[data.mm$occs.x>0 & data.mm$occs.y>0])/ 
                        (length(data.mm$accepted_name))) #calculate Jaccard coefficient
  Sim.is[i,9]      <- (gcd.slc(Sim.is[i,3], Sim.is[i,4], Sim.is[i,5], Sim.is[i,6]))
  
  
  Sim.is[i,10]    <- sum(data.mm$occs.x)#sum of all occs in that cell for x
  Sim.is[i,11]    <- sum(data.mm$occs.y) #sum for all in cell y
  Sim.is[i,12]    <- length(data.11$accepted_name) # num of genera in x
  Sim.is[i,13]    <- length(data.22$accepted_name) # num of genera in y 
  
  print(i) }
Sim.is$CutDist <- cut(Sim.is$Dist, breaks = c(0, 2000, 4000, 6000, 8000, 10000, 12000, 
                                              14000, 16000, 18000, 20000),
                      labels = c("0", "2000", "4000", "6000",
                                 "8000", "10000", "12000","14000", 
                                 "16000", "18000"), include.lowest = TRUE) 
Sim.is <- drop_na(Sim.is)


#Victims:
#------------------------------------------------------------------------------

Genus_occs_pt.v <-plyr::ddply(extinct, c("interval.ma", "cell", "accepted_name"),
                         function(df)c(length(df$accepted_name)))
names(Genus_occs_pt.v) <- c("interval.ma", "cell", "accepted_name", "occs")
Genus_occs_pt.v$stage <- NA
Genus_occs_pt.v$stage <- ifelse(Genus_occs_pt.v$interval.ma == 252.17, "Changhsingian", "Induan")

cells.pt.v <-(sort(unique(Genus_occs_pt.v$cell)))
cell.pt.matrix.v <- data.frame(merge(cells.pt.v, cells.pt.v))
cell.pt.matrix.v <-subset(cell.pt.matrix.v, x <y)


split.genusoccs   <- split(Genus_occs_pt.v, Genus_occs_pt.v$stage) 
list2env(split.genusoccs,envir=.GlobalEnv) 


Sim.v           <- data.frame(matrix(ncol=13, nrow=length(cell.pt.matrix.v$x)))

names(Sim.v)    <- c("x.cell", "y.cell","x.lon", "x.lat", "y.lon", 
                     "y.lat", "Czek", "Jacc", "Dist", "Occs.x",
                     "Occs.y", "Gen.x", "Gen.y")

for (i in 1:length(cell.pt.matrix.v$x)) {
  # for  (i in 1:100) {
  
  Sim.v[i,1]     <- cell.pt.matrix.v$x[i] 
  Sim.v[i,2]     <- cell.pt.matrix.v$y[i] 
  
  j <- cell.pt.matrix.v$x[i]
  k <- cell.pt.matrix.v$y[i]
  
  Sim.v[i,3]     <- cells1$long[cells1$cells==j] 
  Sim.v[i,4]     <- cells1$lat[cells1$cells==j] 
  Sim.v[i,5]     <- cells1$long[cells1$cells==k]
  Sim.v[i,6]     <- cells1$lat[cells1$cells==k] 
  
  data.1        <- subset(Changhsingian, cell==cell.pt.matrix.v$x[i]) 
  data.2        <- subset(Changhsingian, cell==cell.pt.matrix.v$y[i])
  data.m        <- merge(data.1, data.2, by=c("accepted_name"), all=TRUE) 
  data.m[is.na(data.m)] = 0 
  data.m$min    <- pmin(data.m$occs.x, data.m$occs.y)
  
  Sim.v[i,7]      <-2*abs(sum(data.m$min))/((sum(data.m$occs.x)+sum(data.m$occs.y))) 
  Sim.v[i,8]     <- (length(data.m$accepted_name[data.m$occs.x>0 & data.m$occs.y>0])/ 
                       (length(data.m$accepted_name))) 
  Sim.v[i,9]      <- (gcd.slc(Sim.v[i,3], Sim.v[i,4], Sim.v[i,5], Sim.v[i,6]))
  
  
  Sim.v[i,10]    <- sum(data.m$occs.x)
  Sim.v[i,11]    <- sum(data.m$occs.y)
  Sim.v[i,12]    <- length(data.1$accepted_name) 
  Sim.v[i,13]    <- length(data.2$accepted_name) 
  
  print(i) }

Sim.v$CutDist <- cut(Sim.v$Dist, breaks = c(0, 2000, 4000, 6000, 8000, 10000, 12000, 
                                            14000, 16000, 18000, 20000),
                     labels = c("0", "2000", "4000", "6000",
                                "8000", "10000", "12000","14000", 
                                "16000", "18000"), include.lowest = TRUE) 
Sim.v <- drop_na(Sim.v)


#Subset for minimium no. of occs' to compare:
#---------------------------------------------------------------------------------
Sim.s<- subset(Sim.s, Occs.x >20 & Occs.y> 20) #survivors in Changhsingian
Sim.is<- subset(Sim.is, Occs.x >20 & Occs.y >20)#survivors in Induan
Sim.v<- subset(Sim.v, Occs.x >20 & Occs.y>20) #victims

#Find the average and sd
#---------------------------------------------------------------------------------
detach(package:plyr)
#Survivors - Changhsingian
Sum.s.Cz<- Sim.s %>% 
  group_by(CutDist) %>% 
  summarize(avg =mean(Czek), sdev = sd(Czek))
Sum.s.Cz$label <- 'Cha'
Sum.s.Cz <-as.data.frame(Sum.s.Cz)

Sum.s.Ja <- Sim.s%>% 
  group_by(CutDist) %>% 
  summarize(avg =mean(Jacc),sdev = sd(Jacc))
Sum.s.Ja$label <- 'Cha'
Sum.s.Ja <-as.data.frame(Sum.s.Ja)

#S - Induan:
Sum.is.Cz <- Sim.is %>% 
  group_by(CutDist) %>% 
  summarize(avg =mean(Czek),sdev = sd(Czek))
Sum.is.Cz$label <- 'Ind'
Sum.is.Cz <-as.data.frame(Sum.is.Cz)

Sum.is.Ja <-Sim.is%>%
  group_by(CutDist) %>%
  summarize(avg =mean(Jacc),sdev = sd(Jacc))
Sum.is.Ja$label <- 'Ind'
Sum.is.Ja <-as.data.frame(Sum.is.Ja)

#Victims:
Sum.v.Cz <- Sim.v %>% 
  group_by(CutDist) %>% 
  summarize(avg =mean(Czek),sdev = sd(Czek))
Sum.v.Cz$label <- 'Cha'
Sum.v.Cz <-as.data.frame(Sum.v.Cz)

Sum.v.Ja <-Sim.v%>%
  group_by(CutDist) %>%
  summarize(avg =mean(Jacc),sdev = sd(Jacc))
Sum.v.Ja$label <- 'Cha'
Sum.v.Ja <-as.data.frame(Sum.v.Ja)

#Create plots to visualize data
#---------------------------------------------------------------------------------
library(ggplot2)

#Czekanowski versus Jaccard:

#Survivors:
Czek.s <- rbind(Sum.s.Cz, Sum.is.Cz)
Jacc.s<- rbind(Sum.s.Ja, Sum.is.Ja)

#Victims:
Czek.v <- Sum.v.Cz
Jacc.v <- Sum.v.Ja

#Create plot for comparing V vs. S. before extinction:
Czek.v$label2 <- "Victims"
Czek.s$label2 <-"Survivors"

Czek.scha <-subset(Czek.s, label == "Cha")
Czek.scha$label2 <- "Survivors - Before"
Czek.sind <-subset(Czek.s, label =="Ind")
Czek.sind$label2 <- "Survivors - After"

#add bubbles weighted by number of cells compared to one another:
table(Sim.s$CutDist)
table(Sim.is$CutDist)
table(Sim.v$CutDist)

Czek.scha$bubble <- c(5,3,7,8,7,2,2,2)
Czek.sind$bubble <-c(18,11,25,23,21,24,2,10,2)
Czek.v$bubble <-c(12,5,12,12,13,5,1,4,2)

Czek.sv <-rbind(Czek.scha, Czek.sind, Czek.v)

Czek.sv$CutDist=as.numeric(levels(Czek.sv$CutDist))[Czek.sv$CutDist] #make this numeric for x-axis scale
obf <- c("slateblue3","lightseagreen", "violetred")
plot10 <- Czek.sv %>%
  ggplot( aes(x=CutDist, y=avg, group=as.factor(label2),
              colour=as.factor(label2),size = bubble))+
  scale_colour_manual(name="Age", values=obf)+
  theme_classic()+
  geom_point (alpha=0.8)+
  scale_size_area()+
  geom_line(size=2)+
  theme(text = element_text(size=25),
        axis.text.x = element_text(angle=0, hjust=0.5, size=25),
        axis.text.y = element_text(size=25)) +
  ylim(0,0.65)+
  scale_x_continuous(limits=c(0, 16000), breaks=c(0, 2000, 4000, 6000, 8000, 10000,
                                                  12000, 14000, 16000),
                     labels= as.character(c("2000", "4000", "6000",
                                            "8000", "10000", "12000",
                                            "14000", "16000", "18000")))+ 
  geom_errorbar(aes(ymin=pmax(avg-sdev,0), ymax=avg+sdev), width=1,size=1.5, alpha=0.7) + 
  theme(plot.margin = unit(c( t= 0.25, r=0.5, b=0.75,l=.25), "cm"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        plot.title = element_text(size=20), legend.position= c(.7, .9))+ theme(legend.title = element_blank())+
  guides(size = "none")+ #remove legend for bubbles
  labs (title = "By Survival Status", fontface = "bold")
plot10

bg1<-grid.arrange(plot1.g,plot2.g,plot10, ncol=3, 
                  bottom = textGrob("Great Circle Distance", vjust = -1, 
                                    gp = gpar(cex = 0.95, fontface = "bold", fontsize=16)),
                  left = textGrob("",  
                                  rot = 0, vjust = 0, 
                                  gp = gpar( cex = 0.95)),
                  top = textGrob("", vjust =0, 
                                 gp = gpar(fontface = "bold", fontsize =16, cex = 0.95))) 
bg1

ggsave("sim.pt.jul2023.pdf", plot = grid.arrange(plot1.g,plot2.g,plot10, ncol=3, 
                                               bottom = textGrob("Great Circle Distance (km)", vjust = 0.6, 
                                                                 gp = gpar(cex = 2, fontface="bold", fontsize=15)),
                                               left = textGrob("Similarity", 
                                                               rot = 90, vjust = 0.6,
                                                               gp = gpar( cex = 2, fontface ="bold", fontsize=15))),
                                             #  top = textGrob("Similarity Measures", vjust =0.6,
                                               #               gp = gpar(fontface = "bold", cex = 3))),
                                                height=10, width=35, device="pdf", dpi=500) 



