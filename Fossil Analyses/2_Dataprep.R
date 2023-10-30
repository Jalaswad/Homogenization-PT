#--------------------------------------------------------------------------------------
# This code creates shapefiles in R for use in paleogeographic mapping on R
# It will also partition pbdb data for specified time intervals
#--------------------------------------------------------------------------------------

#This is step 2.

#1. DownloadPbdb.R
#2. Dataprep.R
#3. Hex_sim.R
#4  Survival_PT.R
#5. BC.R

#setup
#--------------------------------------------------------------------------------------
#install.packages("sf", "tidyverse","broom","rgdal", "rgeos")
#install.packages("dggridR", "maps", "geosphere", "reshape2", "fossil")
#install.packages("ms_simplify", "raster")
library(sf)
library(tidyverse)
library(ggplot2)
library(rgdal)
library(rgeos)
library(raster)
library(plyr)
library(dplyr)
library(dggridR)
library(fossil)
library(viridis)
library(png)


setwd("")
pbdb <- read.csv(file='~/Codes/pbdb.aug2023.csv')


#adjust radiometric ages
interval.ma    <- ddply(pbdb,c("early_interval"), function(df)c(min(df$min_ma)))
names(interval.ma) <-c("early_interval", "interval.ma")
pbdb       <- merge(pbdb, interval.ma, by=c("early_interval"))


#find first and last occurrences and merge back into data frame, using min_ma column
fadlad <- ddply(pbdb, c("accepted_name"), function(df)c(max(df$interval.ma), #fad
                                                            min(df$interval.ma))) #lad
names(fadlad) <- c("accepted_name", "fad", "lad")

#merge fad and lad information into data frame
pbdb <- merge(pbdb, fadlad, by=c("accepted_name"))

#add extinction/survivor binary variable
pbdb$ex <- 0
pbdb$ex[pbdb$interval.ma==pbdb$lad] <- 1


# Read in shapefiles using sf package, then transform to use with tidyverse
#--------------------------------------------------------------------------------------
ma.252 <- readOGR(dsn=".", layer="252.scotese_polygon")
ma.252 <- gSimplify(ma.252, #gsimplify is a function from rdgal that simplifies shapefiles for R
                    tol = 0.05, # tol(erance): the bigger the value, the coarser the geometry will be
                    topologyPreserve = TRUE) 
ma.252 <- fortify(ma.252)
save(ma.252, file="/Users/Neville/Desktop/PT/Codes/ma.252.Rda")
paleogeo     <- ma.252


#--------------------------------------------------------------------------------------
pbdb       <- unique(subset(pbdb, select = c(reference_no, early_interval, interval.ma,
                                             accepted_name, genus, class,
                                             phylum, paleolat, paleolng, ex, fad)))

pbdb<- subset(pbdb,  (class == "Bivalvia" | class == "Gastropoda") & 
                (interval.ma==252.17 | interval.ma==251.2))
pbdb <- drop_na(pbdb, paleolat, paleolng)

data.cha <- subset(pbdb, interval.ma==252.17)
data.ind <- subset(pbdb, interval.ma == 251.2)
data.cha <-drop_na(data.cha, paleolat,paleolng)  
data.ind <- drop_na(data.ind, paleolat,paleolng)

#Filter to remove wastebin taxa
#These taxa are determined to be 'wastebin genera' in the  2006 Plotnick and Wagner Paper
#--------------------------------------------------------------------------------------
data.cha= filter(data.cha, !(accepted_name %in% c("Inoceramus", "Ostrea","Chlamys","Nuculana",
                                                  "Corbula", "Nucula","Modiolus", "Plagiostoma",
                                                  "Anomia", "Pholadomya", "Dacrydium", "Loxo", "Unio", "Perna (Perna)",
                                                  "Turritella", "Natica","Platyceras","Polinices",
                                                  "Euspira","Gyrodes","Bellerophon","Hormotoma","Conus",
                                                  "Composita", "Atrypa", "Leptaena","Lingula","Chonetes",
                                                  "Derbyia","Crurithyris","Spirifer", "Schizophoria","Camarotoechia")))

data.ind= filter(data.ind, !(accepted_name %in% c("Inoceramus", "Ostrea","Chlamys","Nuculana",
                                                  "Corbula", "Nucula","Modiolus", "Plagiostoma",
                                                  "Anomia", "Pholadomya", "Dacrydium", "Loxo", "Unio", "Perna (Perna)",
                                                  "Turritella", "Natica","Platyceras","Polinices",
                                                  "Euspira","Gyrodes","Bellerophon","Hormotoma","Conus",
                                                  "Composita", "Atrypa", "Leptaena","Lingula","Chonetes",
                                                  "Derbyia","Crurithyris","Spirifer", "Schizophoria","Camarotoechia")))

#Hexgrid processing 
#--------------------------------------------------------------------------------------
data <- subset(data.cha, select = c(interval.ma, early_interval, accepted_name,paleolat, paleolng, ex, fad))
names(data) <- c("interval.ma","early_interval", "accepted_name","lat", "lng", "ex", "fad") 
data.2 <- subset(data.ind, select = c(interval.ma, early_interval, accepted_name, paleolat, paleolng, ex, fad))
names(data.2) <- c("interval.ma","early_interval", "accepted_name","lat","lng", "ex", "fad")
data <- rbind(data, data.2) 

dggs <- dgconstruct(area = 500000) #construct grid of cells with specified area
dgmaxcell(dggs) #returns max number of cells

data.j <- data.frame(matrix(ncol=8, nrow=length(data$accepted_name))) 
names(data.j)      <- c("interval.ma"," early_interval", "accepted_name","lat", "lng", "cell", "ex", "fad")

for (j in 1:length(data.j$interval.ma)){ 
  data.j[j,1]       <- data$interval.ma[j]
  data.j[j,2]       <- data$early_interval[j] 
  data.j[j,3]       <- data$accepted_name[j]  
  data.j[j,4]       <- data$lat[j]
  data.j[j,5]       <- data$lng[j]
  data.j[j,6]       <- dgGEO_to_SEQNUM(dggs, data$lng[j], data$lat[j])$seqnum 
  # assigns coordinates to cells
  data.j[j,7]       <-data$ex[j]
  data.j[j,8]       <-data$fad[j]
  }

cells         <- sort(unique(data.j$cell)) 
cells1 <- as.data.frame(cells) #keeps a log of cells and their locations
cellcenter     <- dgSEQNUM_to_GEO(dggs,cells1$cells) #create a dataframe for the output
#^ will not run properly if you have a newer computer/ intel processor or any M chips.
cells1$lat    <- cellcenter$lat_deg #add latitude for center of each cell
cells1$long   <- cellcenter$lon_deg #add longitude for center of each cell

occs.count<-  data.j %>% group_by(accepted_name, cell) %>% dplyr::summarise(count=n())
#^ calculates how many occurrences are in each cell

grid           <- dgcellstogrid(dggs, occs.count$cell,frame=TRUE,wrapcells=TRUE)
#now we create the grid using the cells we created and the number of genera per cell
grid          <- merge(grid, occs.count,by.x="cell",by.y="cell")
cutoff        <- quantile(grid$count,0.95)
grid          <- grid %>% mutate(count=ifelse(count>cutoff,cutoff,count))

grid20        <- subset(grid, count>20) #subset to 20 or more occ's

#Plots:
#--------------------------------------------------------------------------------------
grid.stage   <- split(grid20, grid20$early_interval) #split for each interval
list2env(grid.stage,envir=.GlobalEnv) #save as separate dataframes

#Chang:
p  <- ggplot() + 
  geom_polygon(data=paleogeo, aes(x=long, y=lat, group=group), fill="gray81", color ="gray58")   +
  geom_polygon(data=Changhsingian,      aes(x=long, y=lat, group=group, fill=count), alpha=0.9)    +
  geom_path   (data=Changhsingian,      aes(x=long, y=lat, group=group),  color="black",alpha=0.8) +
  theme(panel.background = element_rect()) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        plot.title = element_text(size=35), legend.position= "none")+
  scale_fill_viridis(name = "Number of\nOccurrences",
                     limits = c(20, 150), direction=-1, oob = scales::squish)+
  labs(title = "Changhsingian")+
  theme_light()
p 

ggsave(filename="cells_chang20.pdf", plot=q2,
       height=6, width=10, units="in", device="pdf", dpi=500)
q  <- ggplot() + 
  geom_polygon(data=paleogeo, aes(x=long, y=lat, group=group), fill="gray81", color ="gray58")   +
  geom_polygon(data=Induan,      aes(x=long, y=lat, group=group, fill=count), alpha=0.9)    +
  geom_path   (data=Induan,      aes(x=long, y=lat, group=group),  color="black",alpha=0.8) +
  theme(panel.background = element_rect()) +
  xlim(-300,300)+
  scale_fill_viridis(name = "Number of\nOccurrences",
                     limits = c(20, 150), direction=-1, oob = scales::squish)+
  labs(title = "Induan")+
  theme_light()
q 
ggsave(filename="cells_ind20.pdf", plot=q,
       height=6, width=10, units="in", device="pdf", dpi=500)


#Replot on a spherical projection: (This will show that the cells are equal area and not warped)
#--------------------------------------------------------------------------------------
q1 <- q+coord_map("ortho", orientation = c(10, 40, 0))+
  # qq <- q+coord_map("ortho", orientation = c(-40, 80, 0))+ #for second set of photos
  xlab('')+ylab('')+
  theme(axis.ticks.x=element_blank())+
  theme(axis.ticks.y=element_blank())+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_blank())+
  theme_light()+
  theme(plot.title = element_text(size = 20, face="bold"),   plot.subtitle = element_text(size = 18))+
  labs(title = "                          Orientation 1", subtitle=  "Induan")
q1
ggsave(filename="sphere-induan20-Tethys.pdf", plot=q1, bg='transparent',
       height=6, width=10, units="in", device="pdf", dpi=500)

q2 <- q+coord_map("ortho", orientation = c(-40, 80, 0))+ #for second set of photos
  xlab('')+ylab('')+
  theme(axis.ticks.x=element_blank())+
  theme(axis.ticks.y=element_blank())+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_blank())+
  theme_light()+
  theme(plot.title = element_text(size = 20,face="bold"),   plot.subtitle = element_text(size = 18))+
  labs(title ="                          Orientation 2",  subtitle= "Induan")
q2 
ggsave(filename="sphere-induan20-Tethys.pdf", plot=q2, bg='transparent',
       height=6, width=10, units="in", device="pdf", dpi=500)

p1 <- p+coord_map("ortho", orientation = c(10, 40, 0))+
  xlab('')+ylab('')+
  theme(axis.ticks.x=element_blank())+
  theme(axis.ticks.y=element_blank())+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_blank())+
  theme_light()+
  theme(plot.title = element_text(size = 20),   plot.subtitle = element_text(size = 18))+
  labs(title="", subtitle= "Changhsingian")
p1
ggsave(filename="sphere-chang20-Tethys.pdf", plot=p1, bg='transparent',
       height=6, width=10, units="in", device="pdf", dpi=500)

p2 <- p+coord_map("ortho", orientation = c(-40, 80, 0))+ #for second set of photos
  xlab('')+ylab('')+
  theme(axis.ticks.x=element_blank())+
  theme(axis.ticks.y=element_blank())+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_blank())+
  theme_light()+
  theme(plot.title = element_text(size = 20),   plot.subtitle = element_text(size = 18))+
  labs(title="",subtitle= "Changhsingian")
p2
ggsave(filename="sphere-chang20-Tethys.pdf", plot=p2, bg='transparent',
       height=6, width=10, units="in", device="pdf", dpi=500)

cells20fig<- grid.arrange(q1,q2,p1,p2, nrow=2)

ggsave("orientations.pdf", plot = grid.arrange(q1,q2,p1,p2, nrow=2, 
                                                 top= textGrob("Spherical Projection of Cell Locations", vjust = 0.6, 
                                                                   gp = gpar(cex = 2,  fontsize=12))), 
                                               height=15, width=15, device="pdf", dpi=500) 

