
#Author: Alice
#Date: 12/5/22

#library(maptools)
library(ggplot2)
library(sf)
library(ggrepel) #for labeling points
library(ggspatial) #for adding scale bars


#Rnatural earth is a good package for map data

#import map data from shape file
california <- st_read("./ca-state-boundary/CA_State_TIGER2016.shp") #st_read is from sf

#can use st_transform to change CRS 
#CRS is the projection 
#4326 to get to lat-long 
california <- st_transform(california, 4326)

#import coordinates of sample sites
sites_remote <- read.csv("Helianthus map data - annuus remote.csv")
#contains coordinates for only Pumpkin, Manteca, and Mountain

#make a rectangle to add
rectangle <- data.frame(
  x = c(-119.15, -119.8, -119.8, -119.15), #same long as zoomed map
  y = c(36.5, 36.5, 36.9, 36.9) #same lat as zoomed map
)

california_map <- ggplot(data=california)+
  geom_sf()+
  geom_point(data=sites_remote, aes(x=long, y=lat, color=sp), pch=19, cex=2, show.legend = FALSE, alpha=0.6)+
  theme_void(base_family = "Courier", base_size = 8)+
  scale_color_manual(values=c("darkgreen", "blue"))+
     geom_polygon(data = rectangle, aes(x, y, group = 1), fill="lightgray", color="black")+ #adds rectangle around dense sites
  geom_text_repel(data=sites_remote, aes(x=long, y=lat, label=name, color=sp))+
  theme(legend.position = "none")+
  annotation_scale()+
  geom_point(data=dense_sites, aes(x=long, y=lat, color = sp), pch=19, size= 0.5, cex=2, show.legend = FALSE, alpha=0.6)
california_map 


#making a zoomed in-map####

#we need to look at just the dense sites, otherwise we get a map with the other sites plotted outside the rectange
#could have also made a subset of the original table in dplyr
#dense sites are Academy, Boyd, Canal, DelRey, Valley, and Yokuts
dense_sites <- read.csv("Helianthus map data - all dense sites.csv")

#make a subset of the data that's just the area around the sample sites
california_cropped <- st_crop(california, xmin=-119.8, xmax=-119.15,
                              ymin=36.5, ymax=36.9)
#min & max lat of relevant site
#36.55 to 36.75
#min and max long of relevant site
#-119.31 to -119.59

california_zoom <- ggplot(data=california_cropped)+
  geom_sf()+
  geom_point(data=dense_sites, aes(x=long, y=lat,color=sp), pch=16, cex=2, show.legend = FALSE, alpha=0.6)+
  theme_void(base_family = "Courier", base_size = 8)+
  scale_color_manual(values=c("darkgreen", "blue"))+
  geom_text_repel(data=dense_sites, aes(x=long, y=lat, label=name, color=sp), box.padding=0.3)+ #, nudge_x=0.03, nudge_y=0.03)
  theme(legend.position = "none")+
  annotation_scale()
california_zoom

#combine the two maps ####

library(cowplot) #for making inset maps

california_combined = ggdraw() + 
  draw_plot(california_map)+
  draw_plot(california_zoom, x = 0.51, y = 0.50, width = 0.4, height = 0.4)+
  theme_light()

california_combined

#write out file
#for Mol Ecol paper
ggsave("california_map_paper_new.pdf", width=80, units="mm")
#for poster
ggsave("california_map_void_newer.png", width=4.5, height=7, units="in")

