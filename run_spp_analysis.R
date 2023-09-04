#MS Analysis - by species

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
library(emmeans)
library(multcomp)
library(cowplot)
library(MuMIn)

### Inputting data ####
source('load_format_data.R')

#check to see if each company is represented in each region once, some appear more than once 

data<-data %>% 
  filter(!is.na(year))
data %>% 
  count(newregion, cannery.location, company)
data %>% 
  count(salmon.species, anis.no)
data %>% 
  count(salmon.species)
data.2.use<-data[c(3, 5, 10, 11)]

data.trunc<-data %>% 
  filter(salmon.species!= "sockeye") 
length(data.trunc$can.id) #126 cans

data %>% 
  count(newregion, salmon.species)

data$company<-as.character(data$company)

### Checking for spatial and temporal autocorrelation ####
#check autocorrelation for full dataset
#are data spatially autocorrelated?
spatial<-glm(anis.no~year, data=data.trunc)
simspatial<-simulateResiduals(fittedModel = spatial)
length(simspatial)
spatialtest<-testSpatialAutocorrelation(simulationOutput = simspatial,  x = data.trunc$longjitt, y = data.trunc$latjitt)
spatialtest # NOT spatially autocorrelated

#now, check if data are temporally autocorrelated?
temporal<-glm(anis.no~year, data= data.trunc)
temporaltest<-dwtest(temporal, order.by = NULL, alternative = 'two.sided', exact = FALSE, tol = 1e-10)
temporaltest #allsalmon is not temporally autocorrelated
dwtest(temporal)

#red
red<-allsalmon%>%
  filter(salmon.species == "sockeye")
cansgramred<-(red$can.size*28.3495)

spatial.red<-glm(anis.no~year, data=red)
simspatial.red<-simulateResiduals(fittedModel = spatial.red)
spatialtest.red<-testSpatialAutocorrelation(simulationOutput = simspatial.red,  x = red$longjitt, y = red$latjitt)
spatialtest.red  # red IS spatially autocorrelated

red$timegroup<-red$year
temporal.red<-glm(anis.no~year, data= red)
temporaltest.red<-dwtest(temporal.red, order.by = NULL, alternative = 'two.sided', exact = FALSE, tol = 1e-10)
temporaltest.red #red is not temporally autocorrelated --> use glmmPQL for temporal autocorrelation
dwtest(temporal.red) #red is NOT temporally autocorrelated


### Normality checks ####
#Before, we run any analyses, we have to check the distribution of our data....
#lets look at the shape of our data, ideally we would have a bell shape distribution
 
#normal checks, All species together
hist(data.trunc$anis.no) #let's use negative binomial or poission distribution
#This very not normal! But, for good measure, lets look at the residuals
normalityanis<-glm(anis.no~year, data = data.trunc)
resanis<-simulateResiduals(fittedModel = normalityanis, n = 250)
plot(resanis) #residual vs predicted quantile deviations detected
resanis$scaledResiduals
testUniformity(resanis) #deviation significant

#normal checks, red
hist(red$anis.no) #nb or poisson
normalityred<-glm(anis.no~year, data = red)
resred<-simulateResiduals(fittedModel = normalityred, n = 250)
plot(resred) #deviations
resred$scaledResiduals
testUniformity(resred) #deviation significant

# Clearly not normal, as expected.
# For distributions like parasite count data, that are heavily skewed towards 0, 
#and reach 0 very quickly, we usually use a negative binomial distribution or a 
#poisson distribution.


### Models ####

#How does anisakid abundance change over time with species combined?
#Our models likely fit a poisson or negative binomial distribution. 
#we'll 'compete' the models, and see which has the lower AIC (noted in the summary). 
#The lower AIC wins. AIC#s within 2 units of each each 'tie'.

#Aside from spatial and temporal autocorrelation, we should also consider "industrial" autocorrelation - does canning company impact results?
data.trunc$syear<-scale(data.trunc$year)
data.trunc$scansgramall<-scale(data.trunc$cansgramall)

#### Big model selection ####
#First model (2/10/2023) with coho, chum, and pink 
#all caught nearby their cannery
#region should be accurate

NBMod2<-glmmTMB(anis.no~scansgramall +syear*salmon.species + chill+ (1|newregion/factory.code) + (1|newregion/company), data=data.trunc, family="nbinom2")
summary(NBMod2) #within 2 AIC units, keep this as the model to use
#test without outlier
data.trunc2<-data.trunc[-which(data.trunc$anis.no==115),]
NBMod2b<-glmmTMB(anis.no~scansgramall +syear*salmon.species + chill+ (1|newregion/factory.code) + (1|newregion/company), data=data.trunc2, family="nbinom2")

#change reference position
data.trunc$salmon.species[which(data.trunc$salmon.species=="chum")]<-"zchum"
NBMod3<-glmmTMB(anis.no~scansgramall +syear*salmon.species + chill+ (1|newregion/factory.code) + (1|newregion/company), data=data.trunc, family="nbinom2")
summary(NBMod3)
#change ref again
data.trunc$salmon.species[which(data.trunc$salmon.species=="coho")]<-"zcoho"
NBMod4<-glmmTMB(anis.no~scansgramall +syear*salmon.species + chill+ (1|newregion/factory.code) + (1|newregion/company), data=data.trunc, family="nbinom2")
summary(NBMod4)

#changing it back for plotting
data.trunc$salmon.species[which(data.trunc$salmon.species=="zchum")]<-"chum"
data.trunc$salmon.species[which(data.trunc$salmon.species=="zcoho")]<-"coho"

#### Permutation Test for Chum outlier ####
chum<-data.trunc%>%
  filter(salmon.species == "chum")
year<-chum$year
year
anis.no<-chum$anis.no
nsims <- 5000
n <- length(chum$year)
sample.dist <- vector(length=nsims)
test.statistic <- summary(lm(anis.no~year))$coefficients["year","Estimate"]


for (i in 1:nsims) {
  #permute the years
  perm.years <- sample(year, size=n, replace=FALSE)
  #calculate and store the trend estimate
  sample.dist[i] <- summary(lm(anis.no~ perm.years, family="nbinom2"))$coefficients["perm.years","Estimate"]
}
#this is the mean of a bunch of T and F, internally TRUE = 1, FALSE = 0
#so this calculates the proportion greater than the actual test statistic
pvalue <- mean(sample.dist >= test.statistic)


### Plotting ####

my.colors<-c( '#d67237', '#fd6467','#f1bb7b', '#5b1a18')
my.colors2<-c('#d8a499', '#d67237','#991932','#fd6467','#5b1a18' )
my.colors3<-c('#b39fcc', '#991932', '#9b5580', '#641529')

##### Chum, Coho, Pink Combo Plot ####
# need to make predictions based on all three of the fish species
data.trunc %>% 
  count(newregion)
###### Chum plot ####
data.trunc$cansgramall<-as.character(round(data.trunc$cansgramall,0))
data.trunc %>% 
  count(newregion)

#subset data by salmon species for plots
chum<-data.trunc%>%
  filter(salmon.species == "chum")
coho<-data.trunc%>%
  filter(salmon.species=="coho")
pink<-data.trunc%>%
  filter(salmon.species=="pink")

### MAKE glmmTMB model for chum, use this to derive the curve
rescope9 = simulateResiduals(NBMod2)
plot(rescope9)
testDispersion(rescope9)
unique(data.trunc$syear)
unique(data.trunc$scansgramall)
unique(data.trunc$company)
year_pred_df <- data.frame(year = 1979:2020,
                           scansgramall=-1.4795267,
                           newregion="Southeast",
                           salmon.species = "chum",
                           chill="Chill",
                           factory.code="C1",
                           company="2"
)
year_pred_df$syear<-scale(year_pred_df$year)
year_pred <- predict(NBMod2, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(NBMod2, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1979:2020
era<-1979:2020
chum.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() + theme(axis.line = element_line(colour = "black"),
                          axis.text = element_text(color = "black", size=14),
                          axis.title = element_text(size=14),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank()) +
  ylab("Anisakids per can") +  
  xlab(" ")+
  geom_line(aes(y = fit, x=era),
            data = year_pred_df) +
  geom_ribbon(data = year_quantiles, aes(x=era, ymin = lower, ymax = upper), fill='#d8a499', alpha=0.2 ) + 
  theme(text=element_text(size=12,color="black")) +
  geom_point(data = chum,
              aes(y = anis.no, x= year, shape=as.factor(cansgramall)), color='#d8a499', size=2.5, position = position_jitter(width = 0.2, height=0.2))+
  scale_color_manual(values = my.colors2)+
  scale_shape_manual(values=c(17, 15))+
  theme(legend.position = "none")+
  ggtitle("Chum")
chum.plot
minichum<-ggplot(mapping = aes(x = year)) +
  theme_classic() + theme(axis.line = element_line(colour = "black"),
                          axis.text = element_text(color = "black", size =14),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank()) +
  ylab(" ") +  
  xlab(" ")+
  geom_line(aes(y = fit, x=era),
            data = year_pred_df) +
  geom_ribbon(data = year_quantiles, aes(x=era, ymin = lower, ymax = upper), fill='#d8a499', alpha=0.2 ) + 
  theme(text=element_text(size=12,color="black")) +
  geom_point(data = chum,
             aes(y = anis.no, x= year, shape=as.factor(cansgramall)), color='#d8a499', size=2.5, position = position_jitter(width = 0.2, height=0.2))+
  scale_color_manual(values = my.colors2)+
  scale_shape_manual(values=c(17, 15))+
   ylim(0,11)+
  xlim(1980, 2013)+
  theme(legend.position = "none")
minichum

chum<-ggdraw(chum.plot)+ draw_plot(
  {   minichum
  }, x=0.15, y =0.4, 
  width=0.65, height=0.45
)
chum

###### Coho plot ####

coho.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() + theme(axis.line = element_line(colour = "black"),
                          axis.text = element_text(color = "black", size=14),
                          axis.title = element_text(size=14),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank()) +
  ylab("Anisakids per can") + xlab("Year")+ 
  theme(text=element_text(size=12,color="black")) +
  geom_point(data = coho,
              aes(y = anis.no, shape=as.factor(cansgramall)), size=2.5,color='#d67237', position = position_jitter(width = 0.1, height=0.2))+
  scale_color_manual(values = my.colors2[2])+
  scale_shape_manual(values=c(17, 15))+
  theme(legend.position = "none")+
  ggtitle("Coho")
coho.plot

###### Pink Plot #####
unique(data.trunc$syear)
unique(data.trunc$scansgramall)
unique(data.trunc$company)
year_pred_df <- data.frame(year = 1979:2020,
                           scansgramall=-1.4795267,
                           newregion="Southeast",
                           salmon.species = "pink",
                           chill="Chill",
                           factory.code="C1",
                           company="2"
)
year_pred_df$syear<-scale(year_pred_df$year)
year_pred <- predict(NBMod2, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(NBMod2, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1979:2020
era<-1979:2020
pink.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() + theme(axis.line = element_line(colour = "black"),
                          axis.text = element_text(color = "black", size=14),
                          axis.title = element_text(size=14),
                          legend.text = element_text(size=14),
                          legend.title = element_text(size=14),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank()) +
  ylab(" ") +  
  geom_line(aes(y = fit, x=era),
            data = year_pred_df) +
  geom_ribbon(data = year_quantiles, aes(x=era, ymin = lower, ymax = upper), fill='#fd6467', alpha=0.2) + 
  theme(text=element_text(size=12,color="black")) +
  geom_point(data = pink,
              aes(y = anis.no, shape=as.factor(cansgramall)), color = '#fd6467', size=2.5, position = position_jitter(width = 0.1, height=0.2))+
  scale_color_manual(values = my.colors2[4])+
  labs(shape="Can size (g)")+
  xlab(" ")+
  ylim(0,11)+
  ggtitle("Pink")
pink.plot


#### Red model selection ####
#editing code so that we have one model for all the "normal" salmon, and one model to account for sockeye, since they're collection was not as regular
#need to account for changes in time periods

red$syear<-scale(red$year)
red$scansgramall<-scale(red$cansgramall)
red$timeperiod<-NA
red$timeperiod[ c(which(red$year<2000))]<-"period1"
red$timeperiod[ c(which(red$year<2018 & red$year>1999))]<-"period2"
red$timeperiod[ c(which(red$year>=2018))]<-"period3"

#new model, to account for differences in fishing location over time
RedNBMod<-glmmTMB(anis.no~scansgramall +syear + chill+ (1|newregion*timeperiod)+ (1|factory.code) + (1|company), data=red, family="nbinom2")
summary(RedNBMod)
options(na.action="na.fail")
out.put<-MuMIn::dredge(RedNBMod)
#best model:
RedNBMod2<-glmmTMB(anis.no~scansgramall + (1|newregion*timeperiod)+ (1|factory.code) + (1|company), data=red, family="nbinom2")
summary(RedNBMod2)

##### Red plot ####
red$cansgramall2<-as.character(round(red$cansgramall,0))
redcols<-c('#5b1a18', '#AC4743')
#not a significant relationship with time, so we don't need to do the CI or predicted fit, just the raw data
red.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() + theme(axis.line = element_line(colour = "black"),
                          axis.text = element_text(color = "black", size=14),
                          axis.title = element_text(size=14),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank()) +
  ylab(" ") + xlab("Year")+ 
  geom_point(data = red,
              aes(y = anis.no, shape=cansgramall2), color = '#5b1a18', size=2.5, position=position_jitter(width = 0.2, height =0.2))+
  ggtitle("Sockeye")+
  scale_shape_manual(values = c(17, 15))+
  theme(legend.position = "none")



(chum +
    pink.plot)/(coho.plot+red.plot)

ggsave("Allplot.png", height=4.5, width=7)

rescope9 = simulateResiduals(RedNBMod2)
plot(rescope9)
testDispersion(rescope9)
unique(red$year)
unique(red$scansgramall)
unique(red$company)
red %>%
  count(company, factory.code)
year_pred_df <- data.frame(year = 1979:2019,
                           scansgramall=1.0697,
                           newregion="Bristol Bay",
                           chill="Chill",
                           factory.code="C13",
                           company="6",
                           timeperiod ="period2"
)
year_pred_df$syear<-scale(year_pred_df$year)
year_pred <- predict(RedNBMod2, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(RedNBMod2, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1979:2019
era<-1979:2019
red.plot2<-ggplot(mapping = aes(x = year)) +
  theme_classic() + theme(axis.line = element_line(colour = "black"),
                          axis.text = element_text(color = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank()) +
  ylab("Anisakids per can") +
  geom_line(aes(y = fit, x=era),
            data = year_pred_df) +
  geom_ribbon(data = year_quantiles, aes(x=era, ymin = lower, ymax = upper), fill="gray", alpha=0.4 ) +
  theme(text=element_text(size=12,color="black")) +
  geom_jitter(data = red,
              aes(y = anis.no, color=as.factor(cansgramall2)), size=2.5, position = position_jitter(width = 0.1, height=0))+
  labs(color="Can size (g)")+
  ggtitle("Sockeye")
red.plot2



