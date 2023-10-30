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
victims.cha <- subset(data.j, data.j$ex==1 & data.j$interval.ma==252.17) #victims
survivors.cha <- subset(data.j, data.j$ex==0 &data.j$interval.ma==252.17) #survivors in cha

# we do not want to look at originators so we must make sure the fad is not during induan
survivors.ind <- subset(data.j,interval.ma==251.2 & fad != 251.2) #survivors in ind 

survivors <- subset(data.j,interval.ma==251.2) #survivors in ind 


originators <- subset(data.j,interval.ma==251.2 & fad == 251.2)

#Sim for each survival category:
#-------------------------------------------------------------------------------
Genus_occs_pt.s <-plyr::ddply(survivors.cha, c("interval.ma", "cell", "accepted_name"),
                           function(df)c(length(df$accepted_name)))
names(Genus_occs_pt.s) <- c("interval.ma", "cell", "accepted_name", "occs")

cells.pt.s <-(sort(unique(Genus_occs_pt.s$cell))) #all unique cells in the data frame
cell.pt.matrix.s <- data.frame(merge(cells.pt.s, cells.pt.s)) #create two columns with all pairs of values
cell.pt.matrix.s <-subset(cell.pt.matrix.s, x <y) #prevent duplicates and self-comparisons

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
  
  data.1        <- subset(Genus_occs_pt.s, cell==cell.pt.matrix.s$x[i]) 
  data.2        <- subset(Genus_occs_pt.s, cell==cell.pt.matrix.s$y[i])
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
#_________________________________________________________________________
Genus_occs_pt.is <-plyr::ddply(survivors.ind, c("interval.ma", "cell", "accepted_name"),
                               function(df)c(length(df$accepted_name)))
names(Genus_occs_pt.is) <- c("interval.ma", "cell", "accepted_name", "occs")

cells.pt.is <-(sort(unique(Genus_occs_pt.is$cell))) #all unique cells in the data frame
cell.pt.matrix.is <- data.frame(merge(cells.pt.is, cells.pt.is)) #create two columns with all pairs of values
cell.pt.matrix.is <-subset(cell.pt.matrix.is, x <y) #prevent duplicates and self-comparisons

Sim.is             <- data.frame(matrix(ncol=13, nrow=length(cell.pt.matrix.is$x)))

names(Sim.is)      <- c("x.cell", "y.cell","x.lon", "x.lat", "y.lon", 
                        "y.lat", "Czek", "Jacc", "Dist", "Occs.x",
                        "Occs.y", "Gen.x", "Gen.y")

for (i in 1:length(cell.pt.matrix.is$x)) {
  
  Sim.is[i,1]     <- cell.pt.matrix.is$x[i] #record x cell number
  Sim.is[i,2]     <- cell.pt.matrix.is$y[i] #record y cell number
  
  j <- cell.pt.matrix.is$x[i]
  k <- cell.pt.matrix.is$y[i]
  
  Sim.is[i,3]     <- cells1$long[cells1$cells==j] #record x cell lat
  Sim.is[i,4]     <- cells1$lat[cells1$cells==j] #record x cell lng
  Sim.is[i,5]     <- cells1$long[cells1$cells==k] #record y cell lat
  Sim.is[i,6]     <- cells1$lat[cells1$cells==k] #record y cell lng
  
  data.11        <- subset(Genus_occs_pt.is, cell==cell.pt.matrix.is$x[i]) 
  data.22       <- subset(Genus_occs_pt.is, cell==cell.pt.matrix.is$y[i])
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
Genus_occs_pt.v <-plyr::ddply(victims.cha, c("interval.ma", "cell", "accepted_name"),
                         function(df)c(length(df$accepted_name)))
names(Genus_occs_pt.v) <- c("interval.ma", "cell", "accepted_name", "occs")

cells.pt.v <-(sort(unique(Genus_occs_pt.v$cell)))
cell.pt.matrix.v <- data.frame(merge(cells.pt.v, cells.pt.v))
cell.pt.matrix.v <-subset(cell.pt.matrix.v, x <y)

Changhsingian <- Genus_occs_pt.v

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


#Create plot for comparing V vs. S. before victims.chaion:
Czek.v$label2 <- "Victims"
Czek.s$label2 <-"Survivors"

Czek.scha <-subset(Czek.s, label == "Cha")
Czek.scha$label2 <- "Survivors - Cha"
Czek.sind <-subset(Czek.s, label =="Ind")
Czek.sind$label2 <- "Survivors - Ind"

#add bubbles weighted by number of cells compared to one another:
table(Sim.s$CutDist)
table(Sim.is$CutDist)
table(Sim.v$CutDist)

Czek.scha$bubble <-c(3,3,3,3,4,1,2,2)
Czek.sind$bubble <-c(10,2,10,5,8,12,7,1)
Czek.v$bubble <-c(4,5,5,2,2,2,1)

Czek.sv <-rbind(Czek.scha, Czek.sind, Czek.v)

Jacc.scha <-subset(Jacc.s, label == "Cha")
Jacc.scha$label2 <- "Survivors - Cha"
Jacc.sind <-subset(Jacc.s, label =="Ind")
Jacc.sind$label2 <- "Survivors - Ind"
Jacc.v$label2 <- "Victims"

Jacc.scha$bubble <- c(3,3,3,3,4,1,2,2)
Jacc.sind$bubble <-c(10,2,10,5,8,12,7,1)
Jacc.v$bubble <-c(4,5,5,2,2,2,1)

Jacc.sv <-rbind(Jacc.scha, Jacc.sind, Jacc.v)

Jacc.sv$CutDist=as.numeric(levels(Jacc.sv$CutDist))[Jacc.sv$CutDist] #make this numeric for x-axis scale
obf <- c("lightseagreen","slateblue3", "violetred")
plot10j <- Jacc.sv %>%
  ggplot( aes(x=CutDist, y=avg, group=as.factor(label2),
              colour=as.factor(label2),size = bubble))+
  scale_colour_manual(name="Age", values=obf)+
  theme_classic()+
  geom_point (alpha=0.8)+
  scale_size_area()+
  geom_line(size=2)+
  theme(text = element_text(size=45),
        axis.text.x = element_text(angle=0, hjust=0.5, size=30),
        axis.text.y = element_text(size=30)) +
  ylim(0,0.6)+
  scale_x_continuous(limits=c(0, 16000), breaks=c(0, 2000, 4000, 6000, 8000, 10000,
                                                  12000, 14000, 16000),
                     labels= as.character(c("2000", "4000", "6000",
                                            "8000", "10000", "12000",
                                            "14000", "16000", "18000")))+ 
  geom_errorbar(aes(ymin=pmax(avg-sdev,0), ymax=avg+sdev), width=1,size=1.5, alpha=0.7) + 
  theme(plot.margin = unit(c( t= 1, r=1, b=2,l=2), "cm"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        plot.title = element_text(size=35), legend.position= "none")+
  theme(legend.title = element_blank())+
  guides(size = "none")+ #remove legend for bubbles
  labs(title = bquote(~bold('C')~    ' Jaccard (By Survival Status)'))
plot10j

Czek.sv$CutDist=as.numeric(levels(Czek.sv$CutDist))[Czek.sv$CutDist] #make this numeric for x-axis scale
obf <- c("lightseagreen","slateblue3", "violetred")
plot10 <- Czek.sv %>%
  ggplot( aes(x=CutDist, y=avg, group=as.factor(label2),
              colour=as.factor(label2),size =bubble))+
  scale_colour_manual(name="Age", values=obf)+
  theme_classic()+
  geom_point (alpha=0.8)+
  scale_size_area()+
  geom_line(size=2)+
  theme(text = element_text(size=45),
        axis.text.x = element_text(angle=0, hjust=0.5, size=30),
        axis.text.y = element_text(size=30)) +
  ylim(0,0.6)+
  scale_x_continuous(limits=c(0, 16000), breaks=c(0, 2000, 4000, 6000, 8000, 10000,
                                                  12000, 14000, 16000),
                     labels= as.character(c("2000", "4000", "6000",
                                            "8000", "10000", "12000",
                                            "14000", "16000", "18000")))+ 
  geom_errorbar(aes(ymin=pmax(avg-sdev,0), ymax=avg+sdev), width=1,size=1.5, alpha=0.7) + 
  theme(plot.margin = unit(c( t= 1, r=2, b=2,l=1), "cm"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        plot.title = element_text(size=35), legend.position= c(0.8,0.9))+
  theme(legend.title = element_blank())+
  guides(size = "none")+ #remove legend for bubbles
  labs(title = bquote(~bold('D')~    ' Czekanowski (By Survival Status)'))
plot10

bg1<-grid.arrange(plot2.g, plot1.g, plot10j,plot10, ncol=2, 
                  bottom = textGrob("Great Circle Distance", vjust = -1, 
                                    gp = gpar(cex = 0.95, fontface = "bold", fontsize=16)),
                  left = textGrob("",  
                                  rot = 0, vjust = 0, 
                                  gp = gpar( cex = 0.95)),
                  top = textGrob("", vjust =0, 
                                 gp = gpar(fontsize =25, cex = 0.95)))
bg1

ggsave("sim.pt.oct2023.pdf", plot = grid.arrange(plot1.g, plot2.g, plot10j,plot10, ncol=2, 
                                               bottom = textGrob("Great Circle Distance (km)", vjust =0, 
                                                                 gp = gpar(cex = 2,  fontsize=20)),
                                               left = textGrob("Similarity Measures", vjust =0.6,
                                                               rot = 90, gp = gpar( cex =2, fontsize=20))),
                                                height=20, width=27, device="pdf", dpi=500)

write.csv(data.j, file ="~/Desktop/PT/Codes/data.j.oct2023.csv")

