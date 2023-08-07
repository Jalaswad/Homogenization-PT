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
#install.packages("dggridR", "maps", "geosphere", "reshape2", "fossil", "grid","gridExtra")
#install.packages("ms_simplify", "raster")
#remotes::install_github("SebKrantz/dggridR") # ^this version of dggridR is compatible with macbooks with M1/M2 chip
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
library(grid)
library(gridExtra)

setwd("")
pbdb <- read.csv(file='~/Desktop/PT/Codes/pbdb.jul2022.csv')


# Read in shapefiles using sf package, then transform to use with tidyverse
#--------------------------------------------------------------------------------------
ma.252 <- readOGR(dsn=".", layer="252.scotese_polygon")
ma.252 <- gSimplify(ma.252, #gsimplify is a function from rdgal that simplifies shapefiles for R
                    tol = 0.05, # tol(erance): the bigger the value, the coarser the geometry will be
                    topologyPreserve = TRUE) 
ma.252 <- fortify(ma.252)
save(ma.252, file="/Users/Neville/Desktop/PT/Codes/ma.252.Rda")
paleogeo     <- ma.252

#Rework data for correct time allocations to ages, + radiometric ages
#--------------------------------------------------------------------------------------
pbdb       <- unique(subset(pbdb, select = c(reference_no, early_interval,
                                             accepted_name, genus, class,
                                             phylum, paleolat, paleolng)))

pbdb<- subset(pbdb,  (class == "Bivalvia" | class == "Gastropoda") & 
                (early_interval== "Changhsingian" | early_interval== "Induan"))
pbdb <- drop_na(pbdb, paleolat, paleolng)

data.cha <- subset(pbdb, early_interval=="Changhsingian")
data.ind <- subset(pbdb, early_interval == "Induan")
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
data <- subset(data.cha, select = c(early_interval, accepted_name,paleolat, paleolng))
names(data) <- c("early_interval", "accepted_name","lat", "lng") 
data.2 <- subset(data.ind, select = c(early_interval, accepted_name, paleolat, paleolng))
names(data.2) <- c("early_interval", "accepted_name","lat","lng")
data <- rbind(data, data.2) 

dggs <- dgconstruct(area = 500000) #construct grid of cells with specified area
dgmaxcell(dggs) #returns max number of cells

data.j <- data.frame(matrix(ncol=5, nrow=length(data$accepted_name))) 

names(data.j)      <- c("early_interval", "accepted_name","lat", "lng", "cell")
for (j in 1:length(data.j$early_interval)){ 
  data.j[j,1]       <- data$early_interval[j] 
  data.j[j,2]       <- data$accepted_name[j]  
  data.j[j,3]       <- data$lat[j]
  data.j[j,4]       <- data$lng[j]
  data.j[j,5]       <- dgGEO_to_SEQNUM(dggs, data$lng[j], data$lat[j])$seqnum 
  # ^ this is a geographic code which is used to convert numbers to geographic coordinates in the dggridR package.
}

cells         <- sort(unique(data.j$cell)) 
cells1 <- as.data.frame(cells) #keeps a log of cells and their locations
cellcenter     <- dgSEQNUM_to_GEO(dggs,cells1$cells) #create a dataframe for the output
cells1$lat    <- cellcenter$lat_deg #add latitude for center of each cell
cells1$long   <- cellcenter$lon_deg #add longitude for center of each cell

occs.count<-  data.j %>% group_by(early_interval, cell) %>% dplyr::summarise(count=n())
#^ calculates how many occurrences are in each cell

grid           <- dgcellstogrid(dggs, occs.count$cell,frame=TRUE,wrapcells=TRUE)
#now we create the grid using the cells we created and the number of genera per cell
grid          <- merge(grid, gen.count,by.x="cell",by.y="cell")
cutoff        <- quantile(grid$count,0.95)
grid          <- grid %>% mutate(count=ifelse(count>cutoff,cutoff,count))

grid20        <- subset(grid, count>20) #subset to 20 or more occ's

#Plots:
#--------------------------------------------------------------------------------------
grid.stage   <- split(grid20, grid20$early_interval) #split for each interval
list2env(grid.stage,envir=.GlobalEnv) #save as separate dataframes

q  <- ggplot() + 
  geom_polygon(data=paleogeo, aes(x=long, y=lat, group=group), fill="gray81", color ="gray58")   +
  geom_polygon(data=Induan,      aes(x=long, y=lat, group=group, fill=count), alpha=0.9)    +
  geom_path   (data=Induan,      aes(x=long, y=lat, group=group),  color="black",alpha=0.8) +
  theme(panel.background = element_rect()) +
  scale_fill_viridis(name = "Number of\nOccurrences",
                     limits = c(20, 150), direction=-1, oob = scales::squish)+
  labs(title = "Induan")+
  theme_light()
q 
ggsave(filename="cells_ind20.pdf", plot=q,
       height=6, width=10, units="in", device="pdf", dpi=500)

#Chang:
q2  <- ggplot() + 
  geom_polygon(data=paleogeo, aes(x=long, y=lat, group=group), fill="gray81", color ="gray58")   +
  geom_polygon(data=Changhsingian,      aes(x=long, y=lat, group=group, fill=count), alpha=0.9, wrapcells=TRUE)    +
  geom_path   (data=Changhsingian,      aes(x=long, y=lat, group=group),  color="black",alpha=0.8, wrapcells=TRUE) +
  theme(panel.background = element_rect()) +
  scale_fill_viridis(name = "Number of\nOccurrences",
                     limits = c(20, 150), direction=-1, oob = scales::squish)+
  labs(title = "Changhsingian")+
  theme_light()
q2 
ggsave(filename="cells_chang20.pdf", plot=q2,
       height=6, width=10, units="in", device="pdf", dpi=500)

cells20fig<- grid.arrange(q2,q,nrow=2)

#Replot on a spherical projection: (This will show that the cells are equal area and not warped)
qq <- q2+coord_map("ortho", orientation = c(-50, 50, 50))+
  xlab('')+ylab('')+
  theme(axis.ticks.x=element_blank())+
  theme(axis.ticks.y=element_blank())+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_blank())+
  theme_light()+
  labs(title= "600k km cells of genera count",
       subtitle =  "Changhsingian")
qq
ggsave(filename="sphere-changhsingian20-Tethys.pdf", plot=qq, bg='transparent',
       height=6, width=10, units="in", device="pdf", dpi=500)
