# Map plot
# 1 Dec 2022

library(maps)
library(mapdata)
library(maptools) #for shapefiles
library(scales) #for transparency
library(rgdal)
library(maps)
library(mapdata)
library(maptools) #for shapefiles
library(scales) #for transparency
library(ggmap)
library(ggsn)
library(dplyr)
library(ggmap)
library(cowplot)
library(sp)
library(rgeos)
library(raster)
library(maptools)
library(rgdal)
library(wesanderson)
library(dplyr)
library(png)
library(patchwork)

source('load_format_data.R')
colnames(data)
#Map
my.colors<-c('#f1bb7b', '#d67237', '#fd6467', '#5b1a18')
mycols3<-c('#7394d4', '#d8a499','#c6cdf7', '#e6a0c4', '#d67237','#f1bb7b','#fd6467', '#5b1a18')
my.colors2<-c('#d8a499', '#d67237', '#fd6467','#5b1a18' )

a<-ggplot(data, aes(newregion))+
  geom_bar(aes(fill=salmon.species), position="dodge")+
  xlab("Region")+ylab("Count")+
  theme_bw()+scale_fill_manual(values=my.colors2)+
  theme(legend.position = "none", axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text=element_text(size=11),axis.title=element_text(size=14))

bounds<-c(left=-170.35, bottom=53, right=-129.5, top=64.0)
map<-get_stamenmap(bounds, zoom=5, maptype = "terrain-background") %>% ggmap()+
  geom_point(data = data, aes(x=jitter(long,factor=500),y=jitter(lat, factor=500),fill=salmon.species),size=2.5,shape=21,
             alpha=0.8)+
  scale_fill_manual(values = my.colors2)+
  xlab("Longitude (°W)")+
  ylab("Latitude (°N)")+
  guides(fill=guide_legend("Salmon Species"))+
  theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text=element_text(size=11),axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  annotate(geom = "rect", ymax = 61.9, ymin = 59.3, xmax = -153, xmin = -144, colour = "black", fill = NA)+
  annotate(geom = "rect", ymax = 59.3, ymin = 54.3, xmax = -137.2, xmin = -130, colour = "black", fill = NA)+
  annotate(geom = "rect", ymax = 59.6, ymin = 57.5, xmax = -163, xmin = -155.5, colour = "black", fill = NA)+
  annotate(geom="text", x=-145.7, y=59.7, label="Central",color="black", size=4.5)+
  annotate(geom="text", x=-135, y=54.8, label="Southeast",color="black", size=4.5)+
  annotate(geom="text", x=-160.5, y=57.9, label="Bristol Bay",color="black", size=4.5)+
  annotate(geom="text", x=-157.5, y=55.5, label="Western",color="black", size=4.5)+
  annotate(geom="text", x=-145, y=55, label="Gulf of Alaska",color="black", size=6,fontface = 'italic')
map
ggsave("mapfig.png", last_plot())  

ggdraw(map)+ draw_plot(
  {   a
  }, x=0.645, y =0.62, 
  width=0.35, height=0.25
)
ggsave("mapfigdeluxe.png", last_plot())
