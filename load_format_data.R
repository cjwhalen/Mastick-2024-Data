#Load format data

library(tidyverse)
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(lme4)
library(lmtest)
library(nlme)
library(MASS)
library(geoR)
library(spacetime)
library(gstat)
library(sf)
library(sp)
library(patchwork)
library(gdm)
library(vegan)
library(wesanderson)
library(mapview)
#library(faux)
library(emmeans)
library(multcomp)
library(MuMIn)

cans<-read.csv("Can Data_Jan23.csv")
colnames(cans) #weird extra columns, get rid of those
cans<-cans[,1:19]

#this will be used to check for spatial autocorrelation
cans$latjitt<-jitter(cans$lat, factor=0.1, amount=0)
cans$longjitt<--jitter(cans$long, factor=0.1, amount=0)
#cans$latlongjitt<-cbind(cans$longjitt, cans$latjitt)

#make a dataframe that excludes the practice cans
allsalmon<-cans%>%
  filter(box.no. != "practice")%>%
  filter(salmon.species != "")
head(allsalmon)

str(cans)
#Tell R to consider can size a number and not a factor
allsalmon$can.size<-as.numeric(allsalmon$can.size)
allsalmon$anis.no<-as.numeric(allsalmon$anis.no)
allsalmon$cansgramall<-(allsalmon$can.size*28.3495)
allsalmon$year<-as.numeric(allsalmon$year)

### Adding a "region" column ####
#We need to change lat/long into something categorical for the models
unique(allsalmon$lat)
unique(allsalmon$long)
allsalmon$long<- -allsalmon$long
mapview(allsalmon, xcol = "longjitt", ycol = "latjitt", zcol="anis.no", crs = 4269, grid = FALSE)
#aleutians, bristol bay, southeast, kodiak island, anchorage -- check sample sizes of these for each species, see if you need to go to lower region 
allsalmon$region<-NA
allsalmon$region[ c(which(allsalmon$cannery.location=="Dillingham"|
                            allsalmon$cannery.location=="Nushagak"|
                            allsalmon$cannery.location=="Naknek"|
                            allsalmon$cannery.location=="South Naknek"|
                            allsalmon$cannery.location=="Egegik"))]<-"Bristol Bay"
allsalmon$region[ c(which(allsalmon$cannery.location=="King Cove"|
                            allsalmon$cannery.location=="Chignik"))]<-"Aleutians"
allsalmon$region[ c(which(allsalmon$cannery.location=="Alitak"|
                            allsalmon$cannery.location=="Larsen Bay"|
                            allsalmon$cannery.location=="Kodiak"))]<-"Western minus Aleutians"
allsalmon$region[ c(which(allsalmon$cannery.location=="Kenai"|
                            allsalmon$cannery.location=="Seward"|
                            allsalmon$cannery.location=="Valdez"|
                            allsalmon$cannery.location=="Cordova"))]<-"Central"
allsalmon$region[ c(which(allsalmon$cannery.location=="Excursion Inlet"|
                            allsalmon$cannery.location=="Sitka"|
                            allsalmon$cannery.location=="Petersberg"|
                            allsalmon$cannery.location=="Ketchikan"|
                            allsalmon$cannery.location=="Metlakatla"))]<-"Southeast"

#getting rid of any rows that fall outside of this region (Bellingham)
allsalmon<-allsalmon[-which(is.na(allsalmon$region)==T),]
allsalmon$newregion<-allsalmon$region
allsalmon$newregion[ c(which( allsalmon$region=="Western minus Aleutians"|
                                allsalmon$region=="Aleutians"))]<-"Western"
#adding in chilling practices to the model
allsalmon$chill<-NA
allsalmon$chill[c(which(allsalmon$year<1995))]<-"No chill"
allsalmon$chill[c(which(allsalmon$year>=1995))]<-"Chill"

#renaming red to sockeye
allsalmon$salmon.species[ c(which(allsalmon$salmon.species=="red"))]<-"sockeye"
data<-allsalmon[c(2,5,9, 13:14, 16:25)]
