###This r-script record the statistical analysis for wood decomposition project in Mengsong
###Looking at fungal succession on tropical woody debris published in SBB in 2021
###Natural field component
###

#source("http://bioconductor.org/biocLite.R")
#biocLite("phyloseq")


###Install bioconductor from BiocManager
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("phyloseq")

#

##Clean R'sbrain
###Ordered factor understanding
##https://stats.stackexchange.com/questions/339382/results-of-lm-function-with-a-dependent-ordered-categorical-variable

rm(list=ls())
options(digits=6)



##Load needed libraries
library(phyloseq)
library(ggplot2)
library(cowplot)
library(emmeans)
library(MuMIn)
library(car)
library(multcomp)
library(multcompView)
library(car)
library(arm)
library(nlme)
library(rms)
library(AICcmodavg)
library(LMERConvenienceFunctions)
library(plyr)
library(dplyr)
library(tidyr)


####Forest characteristics

forest_data<-read.csv(file="Forests characteristics.csv", header=T, sep=",")
head(forest_data)
str(forest_data)
forest_data$Forest_type<-as.factor(forest_data$Forest_type)
forest_characteristics<-forest_data%>%
  group_by(Forest_type)%>%
  summarize(mean_BA=mean(Basal_area_m2_per_ha),
            sd_BA=sd(Basal_area_m2_per_ha),
            se_BA=sd(Basal_area_m2_per_ha)/sqrt(n()),
            mean_can_openess=mean(Canopy_openness_percent),
            sd_can_openess=sd(Canopy_openness_percent),
            se_can_openess=sd(Canopy_openness_percent)/sqrt(n()),
            mean_forest_age=mean(Forest_age_years, na.rm=T),
            sd_forest_age=sd(Forest_age_years, na.rm=T),
            se_forest_age=sd(Forest_age_years, na.rm=T)/sqrt(n())##n() counts rows in each group
  )
forest_characteristics

###anova on characteristics
###Basal area
characteristic_mod<-aov(data=forest_data,Basal_area_m2_per_ha~Forest_type)
summary(characteristic_mod)
TukeyHSD(characteristic_mod)

##Relevel Forest type
levels(forest_data$Forest_type)
forest_data$Forest_type<-relevel(forest_data$Forest_type, ref="Open_land")
levels(forest_data$Forest_type)

characteristic_mod<-aov(data=forest_data,Basal_area_m2_per_ha~Forest_type)
summary(characteristic_mod)
TukeyHSD(characteristic_mod)

##back to original ref  level
forest_data$Forest_type<-relevel(forest_data$Forest_type, ref="Mature_forest")
levels(forest_data$Forest_type)

###anova on characteristics
####Forest canopy openess
characteristic_mod3<-aov(data=forest_data,Canopy_openness_percent~Forest_type)
summary(characteristic_mod3)
TukeyHSD(characteristic_mod3)
emmeans(characteristic_mod3, specs="Forest_type")

cld(emmeans(characteristic_mod3, specs="Forest_type"))

##Relevel Forest type
levels(forest_data$Forest_type)
forest_data$Forest_type<-relevel(forest_data$Forest_type, ref="Open_land")
levels(forest_data$Forest_type)

characteristic_mod4<-aov(data=forest_data,Basal_area_m2_per_ha~Forest_type)
summary(characteristic_mod4)
TukeyHSD(characteristic_mod4)
cld(emmeans(characteristic_mod4, specs="Forest_type"))



##back to original ref  level
forest_data$Forest_type<-relevel(forest_data$Forest_type, ref="Mature_forest")
levels(forest_data$Forest_type)


###anova on characteristics
####Forest age
characteristic_mod5<-aov(data=forest_data,Forest_age_years~Forest_type)
summary(characteristic_mod5)
TukeyHSD(characteristic_mod5)
emmeans(characteristic_mod5, specs="Forest_type")
cld(emmeans(characteristic_mod5, specs="Forest_type"))

##Relevel Forest type
levels(forest_data$Forest_type)
forest_data$Forest_type<-relevel(forest_data$Forest_type, ref="Open_land")
levels(forest_data$Forest_type)

characteristic_mod6<-aov(data=forest_data,Forest_age_years~Forest_type)
summary(characteristic_mod6)
TukeyHSD(characteristic_mod6)
cld(emmeans(characteristic_mod4, specs="Forest_type"))

### The above analayses results are intergrated with microclimatic information in Dossa et al 2020
### and both integrated in Table 1 of the current MS Dossa et al 2021 SBB

###Import OTUs data
mengsong_97pick <- read.table("mengsong_97closed_guilds_r.csv", sep=",", row.names=1,header=T, check.names=F,blank.lines.skip = FALSE)   
dim(mengsong_97pick)

####Import map data
map <- read.table("mengsong_map_r_chem_SBB.csv",sep=",", header=T, check.names=F,blank.lines.skip = FALSE)
rownames(map) <- map[,1]
dim(map)



#####create subset vectors
subset_sample <- rownames(map)[which(map$control=="No")]
length(subset_sample)
sample_map <- map[subset_sample,]
sample_map <- sample_map[order(rownames(sample_map)),]
sample_map$Forest <- factor(sample_map$Forest)
levels(sample_map$Forest) ###check the forest level
###check the data to make sure the otu table has the correct number of samples
mengsong_97pick[1:2,834:835]

#####create map for separating mengsong samples
Down_map <- sample_map[which(sample_map$SampleType=="Down" |sample_map$SampleType=="Initial"),]
Down_map[] <- lapply(Down_map, function(x) if(is.factor(x)) factor(x) else x)##make sure the factor levels are consistent with the subset dataframe
levels(Down_map$Forest)
dim(Down_map)
summary(Down_map)

## This summary shows 485 samples (initial: 34, down:384, up:67) 
##So actually at time 36 mo because because of fragmentation of down part we only can collect up

##Pilot samples (which are the samples from 12 mo, 36mo)
Pilot_map <- sample_map[sample_map$SampleType=="PilotSample",]

####exclude pilot map
mengsong_map <- sample_map[!sample_map$SampleType %in% "PilotSample",]
UpDown_all_time_map <- mengsong_map[mengsong_map$Updown_pooling_18mo_36mo=="Yes",]
UpDown_all_time_map[] <- lapply(UpDown_all_time_map, function(x) if(is.factor(x)) factor(x) else x)
OpenMature_map <- mengsong_map[mengsong_map$Openmature_pooling=="Yes",]
OpenMature_map[] <- lapply(OpenMature_map, function(x) if(is.factor(x)) factor(x) else x)
dim(OpenMature_map)
summary(OpenMature_map)

###########filter out non fungi otus
subset_fungi <- rownames(mengsong_97pick)[which(mengsong_97pick$Kingdom=="Fungi")]
fungi_97pick <- mengsong_97pick[subset_fungi,]
dim(fungi_97pick)

fungi_97pick_sample <- fungi_97pick[,subset_sample]####only get sample otus
dim(fungi_97pick_sample)
#####exclude less than 10 otus 
fungi_97pick_sample_11 <- subset(fungi_97pick_sample,rowSums(fungi_97pick_sample[,]) >10)
dim(fungi_97pick_sample_11)
fungi_97pick_sample_11 <- fungi_97pick_sample_11[,order(colnames(fungi_97pick_sample_11))]###arrange the colnames in order

###import into phyloseq
otu_phylo <- otu_table(fungi_97pick_sample_11,taxa_are_rows = T)
otu_phylo[1:3,1:2]
sample_design <- sample_data(sample_map)

############put together
otu_together <- phyloseq(otu_phylo,sample_design)
otu_together

####separating the data and calculate alpha diversity index separately
############################################################################Down samples
otu_Down <- subset_samples(otu_together,SampleType=="Down" | SampleType=="Initial")
otu_Down
otu_together

otu_Down_filt <- prune_taxa(taxa_sums(otu_Down)>0,otu_Down)
otu_Down_filt
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 9643 taxa and 485 samples ]
#sample_data() Sample Data:       [ 485 samples by 48 sample variables ]

otu_Down_richness <- estimate_richness(otu_Down_filt,measures=c("Chao1","Shannon"))
qqnorm(otu_Down_richness[,1]) ###check the normalization of data, only transform chao1 
qqnorm(otu_Down_richness[,3])

Chao1_sqrt <- sqrt(otu_Down_richness[,1])
Down_richness_final <- cbind(otu_Down_richness,Chao1_sqrt)
dim(Down_richness_final)
qqnorm(Down_richness_final[,4])

head(Down_richness_final)
head(Down_map)

otu_Down_richness_map <- cbind(Down_richness_final,Down_map)##combine the map file for stat analysis
levels(otu_Down_richness_map$Forest)###check the factor level

###calculate the index and could import into excel for making the table
dim(otu_Down_richness_map)
str(otu_Down_richness_map)
otu_Down_mean_chao <- with(otu_Down_richness_map,aggregate(Chao1,by=list(Forest=Forest,Time=Time,Species=Species,Position=Position),function(x) mean(x,na.rm=T)))
otu_Down_sd_chao <- with(otu_Down_richness_map,aggregate(Chao1,by=list(Forest=Forest,Time=Time,Species=Species,Position=Position),function(x) sd(x,na.rm=T)))

##Calculate SE of initial samples
### Sample size for Litsea (initial samples)= 22
### Sample size for Castanopisis (initial samples)=23 
###Compute SE
mean_otus<-otu_Down_mean_chao[otu_Down_mean_chao$Time=="0mo",]
sd_otus<-otu_Down_sd_chao[otu_Down_sd_chao$Time=="0mo",]
se_otus<-sd_otus$x/(c(sqrt(23),sqrt(22)))
# Castanopsis 69.17099/sqrt(23) # 14.42315 # Therefore mean(SE): 178.6013(14.42315)
# Litsea 154.68577/sqrt(22) #  32.97912 # Therefore mean(SE): 275.3904(32.97912)

otu_Down_mean_shannon <- with(otu_Down_richness_map,aggregate(Shannon,by=list(Forest=Forest,Time=Time,Species=Species,Position=Position),function(x) mean(x,na.rm=T)))

otu_Down_mean_chao2 <- with(otu_Down_richness_map,aggregate(Chao1,by=list(Forest=Forest,Time=Time,Species=Species),function(x) mean(x,na.rm=T)))



###############fit chao1 of Combined Down samples (including the 67 up samples) into model

##Data visualization
####
ggplot(data=otu_Down_richness_map, aes(x=Time, y=Chao1_sqrt,group=Species, colour=Forest)) + geom_point() +
  geom_smooth(method='lm')+
  facet_wrap(~Forest+Species)

ggplot(data=otu_Down_richness_map, aes(x=Time, y=Chao1_sqrt,colour=Species)) + geom_point(position = position_dodge(width = 0.6)) +
  geom_smooth(method='lm')+
  facet_wrap(~Forest)

##########fit chao1 of Down samples (strict down not including the 67 up samples) into model
###Sort the strict down sample
otu_Down <- subset_samples(otu_together,SampleType=="Down" | SampleType=="Initial")
otu_Down

otu_Down_strict <- subset_samples(otu_Down,Position=="down")

otu_Down_strict
otu_Down_strict_filt <- prune_taxa(taxa_sums(otu_Down_strict)>0,otu_Down_strict)

otu_Down_strict_richness <- estimate_richness(otu_Down_strict_filt,measures=c("Chao1","Shannon"))
qqnorm(otu_Down_strict_richness[,1]) ###check the normalization of data,only transform chao1 
qqnorm(otu_Down_strict_richness[,3])

Chao1_strict_sqrt <- sqrt(otu_Down_strict_richness[,1])
Down_richness_strict_final <- cbind(otu_Down_strict_richness,Chao1_strict_sqrt)
dim(Down_richness_strict_final)
qqnorm(Down_richness_strict_final[,4])
str(Down_richness_strict_final)

##Down strict map file
Down_map
levels(Down_map$Position)
head(Down_map)
dim(Down_map)
colnames(Down_map)

Down_Init_map <- Down_map[!Down_map$Position %in% "up",]## Only removes up
dim(Down_Init_map)
Down_strict_map <- Down_Init_map[!Down_Init_map$Position %in% "Initial",]## this removes initial
dim(Down_strict_map)

otu_Down_richness_strict_map <- cbind(Down_richness_strict_final,Down_strict_map)##combine the map file for stat analysis
levels(otu_Down_richness_strict_map$Forest)###check the factor level
str(otu_Down_richness_strict_map)

dim(otu_Down_richness_strict_map)


###calculate the index and could import into excel for making the table
otu_Down_strict_mean_chao <- with(otu_Down_richness_strict_map,aggregate(Chao1,by=list(Forest=Forest,Time=Time,Species=Species,Position=Position,Subplot=Subplot),function(x) mean(x,na.rm=T)))
otu_Down_strict_mean_shannon <- with(otu_Down_richness_strict_map,aggregate(Shannon,by=list(Forest=Forest,Time=Time,Species=Species,Position=Position,Subplot=Subplot),function(x) mean(x,na.rm=T)))

##Plot raw data
ggplot(data=otu_Down_richness_strict_map, aes(x=Time, y=Chao1_strict_sqrt,group=Species, colour=Forest)) + geom_point() +
  geom_smooth(method='lm')+
  facet_wrap(~Forest+Species)

ggplot(data=otu_Down_richness_strict_map, aes(x=Time, y=Chao1_strict_sqrt,colour=Species)) + geom_point(position = position_dodge(width = 0.6)) +
  geom_smooth(method='lm')+
  facet_wrap(~Forest)

###Transform termite to factor
otu_Down_richness_strict_map$Termi.assum
otu_Down_richness_strict_map$Termi.assum<-as.factor(otu_Down_richness_strict_map$Termi.assum)
otu_Down_richness_strict_map$Termi.assum


### Analaysis of wood specific gravity (WSG) loss 
dim(otu_Down_richness_strict_map)
summary(otu_Down_richness_strict_map)

### to only retain the three levels Mature, Regenerating and Open land 
otu_Down_richness_strict_map$Forest
otu_Down_richness_strict_map[] <- lapply(otu_Down_richness_strict_map, function(x) if(is.factor(x)) factor(x) else x)
otu_Down_richness_strict_map$Forest
otu_Down_richness_strict_map$Plot<-as.factor(otu_Down_richness_strict_map$Plot)

otu_Down_richness_strict_map$Subplot<-as.factor(otu_Down_richness_strict_map$Subplot)



########Modeling with three level interaction included
###Import otu table tha contains otus, other information including rot and trophic modes information
otu_Down_richness_strict_map<-read.csv(file="otu_Down_richness_strict_map.csv", header=T, sep=",",stringsAsFactors=T)
##stringsAsFactors=T this helps to put any variable with string as as factor while reading the file
dim(otu_Down_richness_strict_map)
str(otu_Down_richness_strict_map)




###First need to order Time as ordered categorical factor
otu_Down_richness_strict_map$Time<-factor(otu_Down_richness_strict_map$Time, levels=c("18mo","36mo"),ordered=TRUE)
str(otu_Down_richness_strict_map)

###Check homogeneity of variance
bartlett.test(Per_WSG_loss~Forest,data=otu_Down_richness_strict_map)

####Interpretation
#### We found homogeneity of variance among the forest types (Bartlett's K-squared= 0.6646, df=2, P-value= 0.7173).

#### This also confirms that we do not need to include Forest as random factor
####Modeling 
WSG_lme_mod<-lme(data= otu_Down_richness_strict_map, fixed=Per_WSG_loss~Time*Chao1_strict_sqrt*Species*Forest*Termi.assum,random=~1|Plot/Subplot/Wood_log_ID, na.action=na.omit)
summary(WSG_lme_mod)
anova(WSG_lme_mod)


## there is no five way interaction Time chao, species, forest, termites then we update interaction 
###So -Time：Chao1_strict_sqrt:Species:Forest:Termi.assum 

WSG_lme_mod2<-update(WSG_lme_mod, .~.-Time:Chao1_strict_sqrt:Species:Forest:Termi.assum)
summary(WSG_lme_mod2)
anova(WSG_lme_mod2) 
anova(WSG_lme_mod,WSG_lme_mod2) 

##Can't remove five way but we are limiting ourselves to 3 three anyway 
## No Chao1_strict_sqrt:Species:Forest:Termi.assum interaction
##so remove -Chao1_strict_sqrt:Species:Forest:Termi.assum
WSG_lme_mod3<-update(WSG_lme_mod2, .~.-Chao1_strict_sqrt:Species:Forest:Termi.assum)
summary(WSG_lme_mod3)
anova(WSG_lme_mod3) 
anova(WSG_lme_mod2,WSG_lme_mod3) 

##Can't remove that term but we are limiting ourselves to 3 three anyway 
## No Time:Species:Forest:Termi.assum interaction
##so remove -Time:Species:Forest:Termi.assum
WSG_lme_mod4<-update(WSG_lme_mod3, .~.-Time:Species:Forest:Termi.assum)
summary(WSG_lme_mod4)
anova(WSG_lme_mod4) 
anova(WSG_lme_mod3,WSG_lme_mod4) 

##Can't remove that term but we are limiting ourselves to 3 three anyway 

## No Time:Chao1_strict_sqrt:Forest:Termi.assum interaction
##so remove -Time:Chao1_strict_sqrt:Forest:Termi.assum 
WSG_lme_mod5<-update(WSG_lme_mod4, .~.-Time:Chao1_strict_sqrt:Forest:Termi.assum )
summary(WSG_lme_mod5)
anova(WSG_lme_mod5) 
anova(WSG_lme_mod4,WSG_lme_mod5) 

##Can't remove that term but we are limiting ourselves to 3 three anyway 

## No Time:Chao1_strict_sqrt:Species:Termi.assum interaction
##so remove -Time:Chao1_strict_sqrt:Species:Termi.assum 
WSG_lme_mod6<-update(WSG_lme_mod5, .~.-Time:Chao1_strict_sqrt:Species:Termi.assum )
summary(WSG_lme_mod6)
anova(WSG_lme_mod6) 
anova(WSG_lme_mod5,WSG_lme_mod6) 

##Can't remove that term but we are limiting ourselves to 3 three anyway 
## No Time:Chao1_strict_sqrt:Species:Forest interaction 
##so remove -Time:Chao1_strict_sqrt:Species:Forest 
WSG_lme_mod7<-update(WSG_lme_mod6, .~.-Time:Chao1_strict_sqrt:Species:Forest )
summary(WSG_lme_mod7)
anova(WSG_lme_mod7) 
anova(WSG_lme_mod6,WSG_lme_mod7) 


##Can't remove that term but we are limiting ourselves to 3 three anyway 
## No Species:Forest:Termi.assum interaction
##so remove -Species:Forest:Termi.assum 
WSG_lme_mod8<-update(WSG_lme_mod7, .~.-Species:Forest:Termi.assum )
summary(WSG_lme_mod8)
anova(WSG_lme_mod8) 
anova(WSG_lme_mod7,WSG_lme_mod8) 

###Go back to WSG_lme_mod7

##so remove -Chao1_strict_sqrt:Forest:Termi.assum 
WSG_lme_mod9<-update(WSG_lme_mod7, .~.-Chao1_strict_sqrt:Forest:Termi.assum )
summary(WSG_lme_mod9)
anova(WSG_lme_mod9) 
anova(WSG_lme_mod8,WSG_lme_mod9)

## No Chao1_strict_sqrt:Species:Termi.assum interaction
##so remove -Chao1_strict_sqrt:Species:Termi.assum 
WSG_lme_mod10<-update(WSG_lme_mod9, .~.-Chao1_strict_sqrt:Species:Termi.assum )
summary(WSG_lme_mod10)
anova(WSG_lme_mod10) 
anova(WSG_lme_mod9,WSG_lme_mod10)

## No Time:Species:Termi.assum interaction
##so remove -Time:Species:Termi.assum 
WSG_lme_mod11<-update(WSG_lme_mod10, .~.-Time:Species:Termi.assum )
summary(WSG_lme_mod11)
anova(WSG_lme_mod11) 
anova(WSG_lme_mod10,WSG_lme_mod11)

###Go back to WSG_lme_mod10

## No Time:Chao1_strict_sqrt:Termi.assum interaction
##so remove -Time:Chao1_strict_sqrt:Termi.assum 
WSG_lme_mod12<-update(WSG_lme_mod10, .~.-Time:Chao1_strict_sqrt:Termi.assum )
summary(WSG_lme_mod12)
anova(WSG_lme_mod12) 
anova(WSG_lme_mod11,WSG_lme_mod12)

## No Chao1_strict_sqrt:Species:Forest interaction
##so remove -Chao1_strict_sqrt:Species:Forest 
WSG_lme_mod13<-update(WSG_lme_mod12, .~.-Chao1_strict_sqrt:Species:Forest )
summary(WSG_lme_mod13)
anova(WSG_lme_mod13) 
anova(WSG_lme_mod12,WSG_lme_mod13)

## No Time:Species:Forest  interaction
##so remove -Time:Species:Forest 
WSG_lme_mod14<-update(WSG_lme_mod13, .~.-Time:Species:Forest )
summary(WSG_lme_mod14)
anova(WSG_lme_mod14) 
anova(WSG_lme_mod13,WSG_lme_mod14)

###Go back to WSG_lme_mod13
## No Time:Chao1_strict_sqrt:Forest interaction
##so remove -Time:Chao1_strict_sqrt:Forest 
WSG_lme_mod15<-update(WSG_lme_mod13, .~.-Time:Chao1_strict_sqrt:Forest )
summary(WSG_lme_mod15)
anova(WSG_lme_mod15) 
anova(WSG_lme_mod14,WSG_lme_mod15)
 
## No Time:Chao1_strict_sqrt:Species interaction
##so remove -Time:Chao1_strict_sqrt:Species 
WSG_lme_mod16<-update(WSG_lme_mod15, .~.-Time:Chao1_strict_sqrt:Species )
summary(WSG_lme_mod16)
anova(WSG_lme_mod16) 
anova(WSG_lme_mod15,WSG_lme_mod16)

## No Species:Forest:Termi.assum interaction
##so remove -Species:Forest:Termi.assum 
WSG_lme_mod17<-update(WSG_lme_mod16, .~.-Species:Forest:Termi.assum )
summary(WSG_lme_mod17)
anova(WSG_lme_mod17) 
anova(WSG_lme_mod16,WSG_lme_mod17)

###Go back to WSG_lme_mod16
## No Chao1_strict_sqrt:Termi.assum interaction
##so remove -Chao1_strict_sqrt:Termi.assum 
WSG_lme_mod18<-update(WSG_lme_mod16, .~.-Chao1_strict_sqrt:Termi.assum )
summary(WSG_lme_mod18)
anova(WSG_lme_mod18) 
anova(WSG_lme_mod17,WSG_lme_mod18)

## No Chao1_strict_sqrt:Forest interaction
##so remove -Chao1_strict_sqrt:Forest 
WSG_lme_mod19<-update(WSG_lme_mod18, .~.-Chao1_strict_sqrt:Forest )
summary(WSG_lme_mod19)
anova(WSG_lme_mod19) 
anova(WSG_lme_mod18,WSG_lme_mod19)

## No Chao1_strict_sqrt:Species interaction
##so remove -Chao1_strict_sqrt:Species 
WSG_lme_mod20<-update(WSG_lme_mod19, .~.-Chao1_strict_sqrt:Species )
summary(WSG_lme_mod20)
anova(WSG_lme_mod20) 

anova(WSG_lme_mod19,WSG_lme_mod20)


###Diagnostics
plot(WSG_lme_mod20, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Increasing tend
plot(WSG_lme_mod20, abs(resid(., type='pearson'))~fitted(.)|Species, type=c('p', 'r'), 
     abline=0)##

####Perhaps use error correlation weight to correct this

###Allow differnt variace per species
WSG_lme_mod20a<- update(WSG_lme_mod20, weights=varIdent(form=~1|Species))
summary(WSG_lme_mod20a)

anova(WSG_lme_mod20a)

plot(WSG_lme_mod20a, abs(resid(., type='pearson'))~fitted(.)|Species, type=c('p', 'r'), 
     abline=0)##increasing trend

plot(WSG_lme_mod20a, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Still increasing trend


###Use VarExp
WSG_lme_mod20b<-update(WSG_lme_mod20, weights=varExp(form=~fitted(.)))

summary(WSG_lme_mod20b)
anova(WSG_lme_mod20b)

plot(WSG_lme_mod20b, abs(resid(., type='pearson'))~fitted(.)|Species, type=c('p', 'r'), 
     abline=0)##
plot(WSG_lme_mod20b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Increasing trend


###Use VarPower
WSG_lme_mod20c<-update(WSG_lme_mod20, weights=varPower(form=~fitted(.)))

summary(WSG_lme_mod20c)
anova(WSG_lme_mod20c)

plot(WSG_lme_mod20c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(WSG_lme_mod20c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 

AIC(WSG_lme_mod20,WSG_lme_mod20a,WSG_lme_mod20b,WSG_lme_mod20c)

##Model WSG_lme_mod20b or WSG_lme_mod20c best based on AIC， diagnostics and compare with WSG_lme_all

###Best model 
## We chose WSG_lme_mod20b

#RsquareAdj()
r.squaredGLMM(WSG_lme_mod20b)###
##   R2m      R2c
###  0.5327679 0.7803119
summary(WSG_lme_mod20b)

anova(WSG_lme_mod20b)


####Change level to Open land
####otu_Down_richness_strict_map
otu_Down_richness_strict_map$Forest
otu_Down_richness_strict_map$Forest<-relevel(otu_Down_richness_strict_map$Forest,ref="Open_land")
levels(otu_Down_richness_strict_map$Forest)


####Run model again
WSG_lme_mod20b_open<-update(WSG_lme_mod20, weights=varExp(form=~fitted(.)))

summary(WSG_lme_mod20b_open)
anova(WSG_lme_mod20b_open)
r.squaredGLMM(WSG_lme_mod20b_open)

plot(WSG_lme_mod20b_open, abs(resid(., type='pearson'))~fitted(.)|Species, type=c('p', 'r'), 
     abline=0)##
plot(WSG_lme_mod20b_open, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Increasing trend


####Change back to Mature forest as baseline level
####otu_Down_richness_strict_map
otu_Down_richness_strict_map$Forest
otu_Down_richness_strict_map$Forest<-relevel(otu_Down_richness_strict_map$Forest,ref="Mature_forest")
levels(otu_Down_richness_strict_map$Forest)


###
###Remove all terms involving Chao1
WSG_lme_mod21<-update(WSG_lme_mod20,.~.-Chao1_strict_sqrt-Time:Chao1_strict_sqrt)

summary(WSG_lme_mod21)
anova(WSG_lme_mod21)
r.squaredGLMM(WSG_lme_mod21)

###Diagnostics
plot(WSG_lme_mod21, abs(resid(., type='pearson'))~fitted(.)|Species, type=c('p', 'r'), 
     abline=0)##
plot(WSG_lme_mod21, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Increasing trend


####Use weights
###VarIdent
WSG_lme_mod21a<-update(WSG_lme_mod21,weights=varIdent(form=~1|Species))

summary(WSG_lme_mod21a)
anova(WSG_lme_mod21a)

plot(WSG_lme_mod21a, abs(resid(., type='pearson'))~fitted(.)|Species, type=c('p', 'r'), 
     abline=0)##
plot(WSG_lme_mod21a, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Increasing trend


###VarExp
WSG_lme_mod21b<-update(WSG_lme_mod21,weights=varExp(form=~fitted(.)))

summary(WSG_lme_mod21b)
anova(WSG_lme_mod21b)

plot(WSG_lme_mod21b, abs(resid(., type='pearson'))~fitted(.)|Species, type=c('p', 'r'), 
     abline=0)##
plot(WSG_lme_mod21b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Better



###VarPower
WSG_lme_mod21c<-update(WSG_lme_mod21,weights=varPower(form=~fitted(.)))

summary(WSG_lme_mod21c)
anova(WSG_lme_mod21c)

plot(WSG_lme_mod21c, abs(resid(., type='pearson'))~fitted(.)|Species, type=c('p', 'r'), 
     abline=0)##
plot(WSG_lme_mod21c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Better

AIC(WSG_lme_mod21,WSG_lme_mod21a,WSG_lme_mod21b,WSG_lme_mod21c)

##Based on AIC we select WSG_lme_mod21b or WSG_lme_mod21c are equally the best

#### Best model 
### We choose WSG_lme_mod21b
r.squaredGLMM(WSG_lme_mod21b)

summary(WSG_lme_mod21b)
#> r.squaredGLMM(WSG_lme_mod21b)
#       R2m = fixed R2c=fixed+random
#[1,] 0.547206 0.7557998

summary(WSG_lme_mod21b)
anova(WSG_lme_mod21b)

###We can conclude that Chao1 explains (R2c of model with Chao- R2c model without Chao1) 0.0246= 0.7803-0.7557 
r.squaredGLMM(WSG_lme_mod21b)## Without Chao1
r.squaredGLMM(WSG_lme_mod20b)##With Chao1

###Relevels forest ref="Open_land"
otu_Down_richness_strict_map$Forest
levels(otu_Down_richness_strict_map$Forest)
otu_Down_richness_strict_map$Forest<-relevel(otu_Down_richness_strict_map$Forest, ref="Open_land")
levels(otu_Down_richness_strict_map$Forest)

WSG_lme_mod21d<-update(WSG_lme_mod21,weights=varExp(form=~fitted(.)))

summary(WSG_lme_mod21d)

anova(WSG_lme_mod21d)


#> anova(WSG_lme_mod21d)
#numDF denDF   F-value p-value
#(Intercept)                    1   147 1036.2505  <.0001
#Time                           1   147  164.6399  <.0001
#Species                        1   102    2.1104  0.1494
#Forest                         2    24    0.2263  0.7992
#Termi.assum                    1   147    6.6633  0.0108
#Time:Species                   1   147    5.5725  0.0196
#Time:Forest                    2   147    1.1805  0.3100
#Species:Forest                 2   102    0.8213  0.4428
#Time:Termi.assum               1   147    1.4664  0.2279
#Species:Termi.assum            1   147    4.1721  0.0429
#Forest:Termi.assum             2   147    0.0706  0.9319
#Time:Species:Forest            2   147    1.4699  0.2333
#Time:Species:Termi.assum       1   147    3.1091  0.0799
#Time:Forest:Termi.assum        2   147    4.1211  0.0181
#Species:Forest:Termi.assum     2   147    1.3890  0.2526

#anova(WSG_lme_mod22d,type="marginal")

plot(WSG_lme_mod21d, abs(resid(., type='pearson'))~fitted(.)|Species, type=c('p', 'r'), 
     abline=0)##
plot(WSG_lme_mod21d, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Stlight increasing trend


r.squaredGLMM(WSG_lme_mod21d)

###Relevel back to forest ref="Mature forest"

otu_Down_richness_strict_map$Forest
levels(otu_Down_richness_strict_map$Forest)
otu_Down_richness_strict_map$Forest<-relevel(otu_Down_richness_strict_map$Forest, ref="Mature_forest")
levels(otu_Down_richness_strict_map$Forest)

####Create a new data for prediction 


###Best model WSG_lme_mod


##prediction
str(otu_Down_richness_strict_map)
otu_Down_richness_strict_map$Plot<-as.factor(otu_Down_richness_strict_map$Plot)
otu_Down_richness_strict_map$Subplot<-as.factor(otu_Down_richness_strict_map$Subplot)

preddat_18_36mo <- expand.grid(Species = c("Castanopsis_mekongensis","Litsea_cubeba"),
                       Time=c("18mo","36mo"),
                       Plot=unique(otu_Down_richness_strict_map$Plot),
                       Subplot=unique(otu_Down_richness_strict_map$Subplot),
                       Termi.assum = c("0", "1"))

dim(preddat_18_36mo)
head(preddat_18_36mo)
summary(preddat_18_36mo)

str(otu_Down_richness_strict_map)
str(preddat_18_36mo)
preddat_18_36mo$Plot<-as.factor(preddat_18_36mo$Plot)
preddat_18_36mo$Subplot<-as.factor(preddat_18_36mo$Subplot)

#write.csv(preddat_18_36mo,file="preddat_18_36mo.csv")
###We added a new column about Forest type accordingly 

preddat <- read.csv("preddat_18_36mo.csv")
head(preddat)
dim(preddat)
str(preddat)

###transform plot, subplot, forest, and termi.assum in factors

preddat$Plot<-as.factor(preddat$Plot)
preddat$Subplot<-as.factor(preddat$Subplot)
preddat$Termi.assum<-as.factor(preddat$Termi.assum)
str(preddat)


#Best model for WSG without chao1

WSG_lme_mod21b
summary(WSG_lme_mod21b)
anova(WSG_lme_mod21b)
r.squaredGLMM(WSG_lme_mod21b)


hist(otu_Down_richness_strict_map$Per_WSG_loss)



###
preddat$pred <- predict(WSG_lme_mod21b, newdata=preddat, level=0)

summary(preddat$pred)
###Prediction interval 
## [-2] drops response from formula

Designmat <- model.matrix(formula(WSG_lme_mod21b)[-2], preddat)## [-2] drops response from formula
predvar <- diag(Designmat %*% vcov(WSG_lme_mod21b) %*% t(Designmat)) 
preddat$SE <- sqrt(predvar) 
preddat$SE2 <- sqrt(predvar+WSG_lme_mod21b$sigma^2) ## sigma= residual 

preddat$SE_uc<-preddat$pred+1.96*preddat$SE
preddat$SE_lc<-preddat$pred-1.96*preddat$SE
head(preddat)

###Convert variables to categorical
preddat$Species
preddat$Species<-as.factor(preddat$Species)
levels(preddat$Species)
preddat$Time
preddat$Time<-as.factor(preddat$Time)
levels(preddat$Time)
preddat$Termi.assum
preddat$Termi.assum<-as.factor(preddat$Termi.assum)
levels(preddat$Termi.assum)
preddat$Forest
preddat$Forest<-as.factor(preddat$Forest)
levels(preddat$Forest)
preddat$Forest<-factor(preddat$Forest, levels = c("Mature_forest","Regenerating_forest","Open_land"))
levels(preddat$Forest)
head(preddat)

library(ggplot2)
library(gridExtra)
library(reshape)
library(grid)

str(preddat)


summary(preddat)



##graphing
##Relevel Species as Litsea cubeba and Castanopsis mekongensis
levels(preddat$Species)
levels(preddat$Species) <- c("Castanopsis mekongensis","Litsea cubeba")
###Make sure Litsea is taken as reference level
preddat$Species<-relevel(preddat$Species, ref="Litsea cubeba")
levels(preddat$Species)

levels(preddat$Forest)
levels(preddat$Forest) <- c("Mature forest", "Regenerating forest","Open land")

levels(preddat$Termi.assum)
levels(preddat$Termi.assum)<-c("Absence", "Presence")
levels(preddat$Termi.assum)

summary(preddat)
WSG_loss_18_36mo_graph<-ggplot(data=preddat, aes(y=pred, x=Time, colour= Forest,shape=Species, fill=Termi.assum))+
  geom_rect(data=NULL,aes(xmin=1.5,xmax=2.7,ymin=-Inf,ymax=Inf),
            fill="lightgray", colour = NA, alpha = 0.1) +
  scale_shape_manual(values=c(24,21)) +
  scale_fill_manual(values=c(NA, "black"),guide=guide_legend(override.aes=list(shape=21)))+
  geom_point(position=position_dodge(1),stat="identity", size=6)+ 
  geom_errorbar(aes(ymax = preddat$SE_uc, ymin = preddat$SE_lc), width=0, size = 0.8, position=position_dodge(1))+
  scale_y_continuous(limits=c(5,72))+
  scale_x_discrete(limits=levels(preddat$Time))+
  xlab(" Exposure duration") +
  ylab("Percentage wood specific gravity loss") + 
  labs(colour="Habitat type")+
  labs(shape="Woody species")+
  labs(fill="Termites status")+
  annotate("text", x=2.3, y=72, label= "Time p-value<0.001          ", size=5)+
  annotate("text", x=2.3, y=70, label= "Termites p-value=0.017",size=5) +
  annotate("text", x=2.3, y=68, label= "Species:Time p-value=0.020",size=5)+
  annotate("text", x=2.3, y=66, label= "      Species:Termites p-value=0.043      ",size=5)+
  annotate("text", x=2.3, y=64, label= "      Time:Habitat:Termites p-value=0.018      ",size=5)+
  #annotate("text", x=2.3, y=62, label= "      Time:Species:Termites p-value=0.070      ",size=5)+
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        strip.text = element_text(colour = "black", size = 15),
        #legend.position = c("bottom"),
        legend.position=c(0,1), 
        legend.justification=c(0,1),
        #legend.position = c(0.1,0.7),
        legend.text = element_text(colour="black",size=15,face="italic"),
        legend.title=element_text(size=15),
        #legend.text=element_text(size=12,face="italic"),
        #legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm")) 

WSG_loss_18_36mo_graph

  
  
###Figure 2 WSG dynamics 
tiff(filename="Figure 2 WSG 18 36mo with termites status.tiff", res = 400, width=5000, height=3500, compression = "lzw")
WSG_loss_18_36mo_graph
dev.off()

####WSG modelling considering rot type and trophic type
str(otu_Down_richness_strict_map)
WSG_lme_all<-lme(data= otu_Down_richness_strict_map, fixed=Per_WSG_loss~Time*Chao1_strict_sqrt*Species*Forest*Termi.assum+Time*White_Rot+Chao1_strict_sqrt*White_Rot+Species*White_Rot+Forest*White_Rot+Termi.assum*White_Rot+
                   Time*Brown_Rot+Chao1_strict_sqrt*Brown_Rot+Species*Brown_Rot+Forest*Brown_Rot+Termi.assum*Brown_Rot+Time*Soft_Rot+Chao1_strict_sqrt*Soft_Rot+Species*Soft_Rot+Forest*Soft_Rot+Termi.assum*Soft_Rot+
                   Time*Saprotroph+Chao1_strict_sqrt*Saprotroph+Species*Saprotroph+Forest*Saprotroph+Termi.assum*Saprotroph, random=~1|Plot/Subplot/Wood_log_ID, na.action=na.omit)
summary(WSG_lme_all)
anova(WSG_lme_all)


#### Simplification
WSG_lme_all

## there is no five levels chao, species, forest, termites then we update interaction 
###So -Time：Chao1_strict_sqrt:Species:Forest:Termi.assum 

WSG_lme_all2<-update(WSG_lme_all, .~.-Time:Chao1_strict_sqrt:Species:Forest:Termi.assum)
summary(WSG_lme_all2)
anova(WSG_lme_all2) 
anova(WSG_lme_all,WSG_lme_all2) 

##Can't remove but we are not considering five ways interaction so we drop that term

## remove four ways interaction
##so remove -Chao1_strict_sqrt:Species:Forest:Termi.assum
WSG_lme_all3<-update(WSG_lme_all2, .~.-Chao1_strict_sqrt:Species:Forest:Termi.assum)
summary(WSG_lme_all3)
anova(WSG_lme_all2,WSG_lme_all3) 

##Can't remove but we are not considering five ways interaction so we drop that term

## remove four ways interaction
##so remove -Time:Species:Forest:Termi.assum
WSG_lme_all4<-update(WSG_lme_all3, .~.-Time:Species:Forest:Termi.assum)
summary(WSG_lme_all4)
anova(WSG_lme_all3, WSG_lme_all4) 

##Can't remove but we are not considering five ways interaction so we drop that term

## remove four ways interaction
##so remove -Time:Chao1_strict_sqrt:Forest:Termi.assum 
WSG_lme_all5<-update(WSG_lme_all4, .~.-Time:Chao1_strict_sqrt:Forest:Termi.assum )
summary(WSG_lme_all5)
anova(WSG_lme_all5) 
anova(WSG_lme_all4, WSG_lme_all5)

##Can't remove but we are not considering five ways interaction so we drop that term

## Remove four ways interaction
##so remove -Time:Chao1_strict_sqrt:Species:Termi.assum 
WSG_lme_all6<-update(WSG_lme_all5, .~.-Time:Chao1_strict_sqrt:Species:Termi.assum )
summary(WSG_lme_all6)
anova(WSG_lme_all6) 

##Can't remove but we are not considering five ways interaction so we drop that term

## Remove four ways interaction
##so remove -Time:Chao1_strict_sqrt:Species:Forest 
WSG_lme_all7<-update(WSG_lme_all6, .~.-Time:Chao1_strict_sqrt:Species:Forest )
summary(WSG_lme_all7)
anova(WSG_lme_all7) 
anova(WSG_lme_all6, WSG_lme_all7)

##Can't remove but we are not considering five ways interaction so we drop that term

## Remove three ways interaction
##so remove -Species:Forest:Termi.assum 
WSG_lme_all8<-update(WSG_lme_all7, .~.-Species:Forest:Termi.assum )
summary(WSG_lme_all8)
anova(WSG_lme_all8) 
anova(WSG_lme_all7, WSG_lme_all8)

###Go back to WSG_lme_all7

## Remove three ways interaction
##so remove -Chao1_strict_sqrt:Forest:Termi.assum 
WSG_lme_all9<-update(WSG_lme_all7, .~.-Chao1_strict_sqrt:Forest:Termi.assum )
summary(WSG_lme_all9)
anova(WSG_lme_all9) 
anova(WSG_lme_all7,WSG_lme_all9)

## Remove three ways interaction
##so remove -Chao1_strict_sqrt:Species:Termi.assum 
WSG_lme_all10<-update(WSG_lme_all9, .~.-Chao1_strict_sqrt:Species:Termi.assum )
summary(WSG_lme_all10)
anova(WSG_lme_all10) 
anova(WSG_lme_all9,WSG_lme_all10)

## Remove three ways interaction
##so remove -Time:Species:Termi.assum 
WSG_lme_all11<-update(WSG_lme_all10, .~.-Time:Species:Termi.assum )
summary(WSG_lme_all11)
anova(WSG_lme_all11) 
anova(WSG_lme_all10,WSG_lme_all11)

##Go back to WSG_lme_all10

## Remove three ways interaction
##so remove -Time:Chao1_strict_sqrt:Termi.assum 
WSG_lme_all12<-update(WSG_lme_all10, .~.-Time:Chao1_strict_sqrt:Termi.assum )
summary(WSG_lme_all12)
anova(WSG_lme_all12) 
anova(WSG_lme_all11,WSG_lme_all12)

## Remove three ways interaction
##so remove -Chao1_strict_sqrt:Species:Forest 
WSG_lme_all13<-update(WSG_lme_all12, .~.-Chao1_strict_sqrt:Species:Forest )
summary(WSG_lme_all13)
anova(WSG_lme_all13) 
anova(WSG_lme_all12,WSG_lme_all13)

## Remove three ways interaction
##so remove -Time:Species:Forest 
WSG_lme_all14<-update(WSG_lme_all13, .~.-Time:Species:Forest )
summary(WSG_lme_all14)
anova(WSG_lme_all14) 
anova(WSG_lme_all13,WSG_lme_all14)

##Go back to WSG_lme_13

## Remove three ways interaction
##so remove -Time:Chao1_strict_sqrt:Forest 
WSG_lme_all15<-update(WSG_lme_all13, .~.-Time:Chao1_strict_sqrt:Forest )
summary(WSG_lme_all15)
anova(WSG_lme_all15) 
anova(WSG_lme_all13,WSG_lme_all15)


## Remove three ways interaction
##so remove -Time:Chao1_strict_sqrt:Species 
WSG_lme_all16<-update(WSG_lme_all15, .~.-Time:Chao1_strict_sqrt:Species )
summary(WSG_lme_all16)
anova(WSG_lme_all16) 
anova(WSG_lme_all15,WSG_lme_all16)


## Remove two ways interaction
##so remove -Chao1_strict_sqrt:Termi.assum 
WSG_lme_all17<-update(WSG_lme_all16, .~.-Chao1_strict_sqrt:Termi.assum )
summary(WSG_lme_all17)
anova(WSG_lme_all17) 
anova(WSG_lme_all16,WSG_lme_all17)

## Remove two ways interaction
##so remove -Chao1_strict_sqrt:Forest 
WSG_lme_all18<-update(WSG_lme_all17, .~.-Chao1_strict_sqrt:Forest )
summary(WSG_lme_all18)
anova(WSG_lme_all18) 
anova(WSG_lme_all17,WSG_lme_all18)


## Remove three ways interaction
##so remove -Chao1_strict_sqrt:Species 
WSG_lme_all19<-update(WSG_lme_all18, .~.-Chao1_strict_sqrt:Species )
summary(WSG_lme_all19)
anova(WSG_lme_all19) 
anova(WSG_lme_all18,WSG_lme_all19)


####So far we can't remove Time:Forest,  Time:Termites, Forest:Termites
#### Time:Species, Species:Forest because when we drop Time:Species:Forest because AIC increases

## No two ways interaction
##so remove -Species:Forest 
#WSG_lme_all20<-update(WSG_lme_all19, .~.-Species:Forest )
#summary(WSG_lme_all20)
#anova(WSG_lme_all20) 
#anova(WSG_lme_all19,WSG_lme_all20)

###Remove two ways interaction
##remove Termi.assum:Saprotroph 

WSG_lme_all21<-update(WSG_lme_all19, .~.-Termi.assum:Saprotroph )
summary(WSG_lme_all21)
anova(WSG_lme_all21) 
anova(WSG_lme_all19,WSG_lme_all21)

###Remove two ways interaction
##remove Chao1_strict_sqrt:Saprotroph 

WSG_lme_all22<-update(WSG_lme_all21, .~.-Chao1_strict_sqrt:Saprotroph )
summary(WSG_lme_all22)
anova(WSG_lme_all22) 
anova(WSG_lme_all21,WSG_lme_all22)


###Remove two ways interaction
##remove Time:Saprotroph 

WSG_lme_all23<-update(WSG_lme_all22, .~.-Time:Saprotroph )
summary(WSG_lme_all23)
anova(WSG_lme_all23) 
anova(WSG_lme_all22,WSG_lme_all23)

###Remove two ways interaction
##remove Termi.assum:Soft_Rot 

WSG_lme_all24<-update(WSG_lme_all23, .~.-Termi.assum:Soft_Rot )
summary(WSG_lme_all24)
anova(WSG_lme_all24) 
anova(WSG_lme_all23,WSG_lme_all24)

###Remove two ways interaction
##remove Forest:Soft_Rot 

WSG_lme_all25<-update(WSG_lme_all24, .~.-Forest:Soft_Rot )
summary(WSG_lme_all25)
anova(WSG_lme_all25) 
anova(WSG_lme_all24,WSG_lme_all25)

###Remove two ways interaction
##remove Species:Soft_Rot 

WSG_lme_all26<-update(WSG_lme_all25, .~.-Species:Soft_Rot )
summary(WSG_lme_all26)
anova(WSG_lme_all26) 
anova(WSG_lme_all25,WSG_lme_all26)

###Remove two ways interaction
##remove Time:Soft_Rot 

WSG_lme_all27<-update(WSG_lme_all26, .~.-Time:Soft_Rot )
summary(WSG_lme_all27)
anova(WSG_lme_all27) 
anova(WSG_lme_all26,WSG_lme_all27)

###Remove two ways interaction
##remove Termi.assum:Brown_Rot 

WSG_lme_all28<-update(WSG_lme_all27, .~.-Termi.assum:Brown_Rot )
summary(WSG_lme_all28)
anova(WSG_lme_all28) 
anova(WSG_lme_all27,WSG_lme_all28)

###Remove two ways interaction
##remove Species:Brown_Rot 

WSG_lme_all29<-update(WSG_lme_all28, .~.-Species:Brown_Rot )
summary(WSG_lme_all29)
anova(WSG_lme_all29) 
anova(WSG_lme_all28,WSG_lme_all29)

###Remove two ways interaction
##remove Chao1_strict_sqrt:Brown_Rot 

WSG_lme_all30<-update(WSG_lme_all29, .~.-Chao1_strict_sqrt:Brown_Rot )
summary(WSG_lme_all30)
anova(WSG_lme_all30) 
anova(WSG_lme_all29,WSG_lme_all30)

###Remove two ways interaction
##remove Time:Brown_Rot 

WSG_lme_all31<-update(WSG_lme_all30, .~.-Time:Brown_Rot )
summary(WSG_lme_all31)
anova(WSG_lme_all31) 
anova(WSG_lme_all30,WSG_lme_all31)

###Remove two ways interaction
##remove Time:Brown_Rot 

WSG_lme_all32<-update(WSG_lme_all31, .~.-Termi.assum:White_Rot )
summary(WSG_lme_all32)
anova(WSG_lme_all32) 
anova(WSG_lme_all31,WSG_lme_all32)

###Remove two ways interaction
##remove Forest:White_Rot 

WSG_lme_all33<-update(WSG_lme_all32, .~.-Forest:White_Rot)
summary(WSG_lme_all33)
anova(WSG_lme_all33) 
anova(WSG_lme_all32,WSG_lme_all33)

###Remove two ways interaction
##remove Chao1_strict_sqrt:White_Rot 

WSG_lme_all34<-update(WSG_lme_all33, .~.-Chao1_strict_sqrt:White_Rot )
summary(WSG_lme_all34)
anova(WSG_lme_all34) 
anova(WSG_lme_all33,WSG_lme_all34)

###Remove two ways interaction
##remove Species:Termi.assum 

#WSG_lme_all35<-update(WSG_lme_all34, .~.-Species:Termi.assum )
#summary(WSG_lme_all35)
#anova(WSG_lme_all35) 
#anova(WSG_lme_all34,WSG_lme_all35)


###No two ways interaction
##remove Time:Species 

#WSG_lme_all36<-update(WSG_lme_all35, .~.-Time:Species )
#summary(WSG_lme_all36)
#anova(WSG_lme_all36) 
#anova(WSG_lme_all35,WSG_lme_all36)

#AIC(WSG_lme_all35,WSG_lme_all36)
###go back to WSG_lme_all35

###No two ways interaction
##remove Time:White_Rot 

WSG_lme_all37<-update(WSG_lme_all34, .~.-Time:White_Rot )
summary(WSG_lme_all37)
anova(WSG_lme_all37) 
anova(WSG_lme_all34,WSG_lme_all37) 

###Remove two ways interaction
##remove Species:Saprotroph 

WSG_lme_all38<-update(WSG_lme_all37, .~.-Species:Saprotroph )
summary(WSG_lme_all38)
anova(WSG_lme_all38) 
anova(WSG_lme_all37,WSG_lme_all38) 


###Diagnostics

###Diagnostics
plot(WSG_lme_all38, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Increasing tend
plot(WSG_lme_all38, abs(resid(., type='pearson'))~fitted(.)|Species, type=c('p', 'r'), 
     abline=0)##

####Perhaps use error correlation weight to correct

###Allow differnt variace per species
WSG_lme_all38a<- update(WSG_lme_all38, weights=varIdent(form=~1|Species))
summary(WSG_lme_all38a)

anova(WSG_lme_all38a)

plot(WSG_lme_all38a, abs(resid(., type='pearson'))~fitted(.)|Species, type=c('p', 'r'), 
     abline=0)##increasing trend

plot(WSG_lme_all38a, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Still increasing trend


###Use VarExp
WSG_lme_all38b<-update(WSG_lme_all38, weights=varExp(form=~fitted(.)))


summary(WSG_lme_all38b)
anova(WSG_lme_all38b)

plot(WSG_lme_all38b, abs(resid(., type='pearson'))~fitted(.)|Species, type=c('p', 'r'), 
     abline=0)##
plot(WSG_lme_all38b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Slight increasing trend


###Use VarPower
WSG_lme_all38c<-update(WSG_lme_all38, weights=varPower(form=~fitted(.)))

summary(WSG_lme_all38c)
anova(WSG_lme_all38c)

plot(WSG_lme_all38c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##
plot(WSG_lme_all38c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Slight increasing trend

AIC(WSG_lme_all38,WSG_lme_all38a,WSG_lme_all38b,WSG_lme_all38c)

###

##Model WSG_lme_all38b or WSG_lme_all38c best based on AIC and diagnostics
## We chose WSG_lme_all38c

#RsquareAdj()
r.squaredGLMM(WSG_lme_all38c)###
##   R2m      R2c
### 0.678 0.968

###Best model WSG_lme_all38c
anova(WSG_lme_all38c)###Chao1 explains from anova table 
###Table S5
summary(WSG_lme_all38c)

r.squaredGLMM(WSG_lme_all38c)

# r.squaredGLMM(WSG_lme_all38c)
#        R2m       R2c
# [1,]  0.678 0.968

##(2.5953 +15.3281+4.3268)/(1072.8846+212.5221+ 2.5953+1.8666+0.1499+6.2778+2.0479+3.6559+0.0269+3.1406+15.3281+ 1.8989+1.1043+0.9144+2.0130+0.0815+8.5233+3.0764+4.3268+3.7519+0.7709+3.7602+6.0110+1.9980)=0.01638

###For WSG best model is WSG_lme_all38c

#####Change level to Open
###Relevel  to forest ref="Open_land"

otu_Down_richness_strict_map$Forest
levels(otu_Down_richness_strict_map$Forest)
otu_Down_richness_strict_map$Forest<-relevel(otu_Down_richness_strict_map$Forest, ref="Open_land")
levels(otu_Down_richness_strict_map$Forest)


####
WSG_lme_all38d<-update(WSG_lme_all38, weights=varPower(form=~fitted(.)))
###Table S6
summary(WSG_lme_all38d)
anova(WSG_lme_all38d)

plot(WSG_lme_all38d, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##slight increasing trend
plot(WSG_lme_all38d, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## slight increasing trend

r.squaredGLMM(WSG_lme_all38d)

# r.squaredGLMM(WSG_lme_all38d)
#        R2m       R2c
# [1,] 0.678 0.968

#####Change level to Mature
###Relevel  to forest ref="Mature_forest"

otu_Down_richness_strict_map$Forest
levels(otu_Down_richness_strict_map$Forest)
otu_Down_richness_strict_map$Forest<-relevel(otu_Down_richness_strict_map$Forest, ref="Mature_forest")
levels(otu_Down_richness_strict_map$Forest)


####WSG modelling
WSG_lme_all38c


##############Alpha diversity (Chao1) modeling
###The four interaction Time:Forest:Species:Termi.assum

alpha_lme1 <- lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum+Time:Forest+Time:Species+Forest:Species+Time:Termi.assum+Forest:Termi.assum+Species:Termi.assum+Time:Forest:Species+
                       Time:Forest:Termi.assum+Time:Species:Termi.assum+Forest:Species:Termi.assum+Time:Forest:Species:Termi.assum, random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map)


summary(alpha_lme1)
anova(alpha_lme1)


###Remove four way interaction
###Remove Time:Forest:Species:Termi.assum

alpha_lme2<-update(alpha_lme1, .~.-Time:Forest:Species:Termi.assum )
summary(alpha_lme2)
anova(alpha_lme2) 
AIC(alpha_lme2)

###No three way interaction
###Remove Time:Forest:Termi.assum

alpha_lme3<-update(alpha_lme2, .~.-Time:Forest:Termi.assum )
summary(alpha_lme3)
anova(alpha_lme3) 
AIC(alpha_lme3)

###No three way interaction
###Remove Time:Species:Termi.assum

alpha_lme4<-update(alpha_lme3, .~.-Time:Species:Termi.assum )
summary(alpha_lme4)
anova(alpha_lme4) 
AIC(alpha_lme4)

##Go back to model 3
summary(alpha_lme3)
AIC(alpha_lme3)

###No three way interaction
###Remove Time:Forest:Species

alpha_lme5<-update(alpha_lme3, .~.-Time:Forest:Species )
summary(alpha_lme5)
anova(alpha_lme5) 
AIC(alpha_lme5)

###Go back to alpha_lme3 

###Diagnostics
plot(alpha_lme3, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(alpha_lme3, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 


###Use VarExp
alpha_lme3b<-update(alpha_lme3, weights=varExp(form=~fitted(.)))

summary(alpha_lme3b)
anova(alpha_lme3b)

plot(alpha_lme3b, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Slight increasing trend
plot(alpha_lme3b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 



###Use VarPower
alpha_lme3c<-update(alpha_lme3, weights=varPower(form=~fitted(.)))

summary(alpha_lme3c)
anova(alpha_lme3c)

plot(alpha_lme3c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Slight increasing trend
plot(alpha_lme3c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 

AIC(alpha_lme3b,alpha_lme3c)

###Best alpha_lme3b
r.squaredGLMM(alpha_lme3b)

##> r.squaredGLMM(alpha_lme3b)
######R2m       R2c
##[1,] 0.3927463 0.6724915

#### Putting the rot types and Saprotrophs in the model as the dSep suggested


##############Chao1 modeling
###The four interaction Time:Forest:Species:Termi.assum

alpha_rot_lme1 <- lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum+Time:Forest+Time:Species+Forest:Species+Time:Termi.assum+Forest:Termi.assum+Species:Termi.assum+Time:Forest:Species+
                    Time:Forest:Termi.assum+Time:Species:Termi.assum+Forest:Species:Termi.assum+Time:Forest:Species:Termi.assum+White_Rot+Brown_Rot+Soft_Rot+Saprotroph, random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map)


summary(alpha_rot_lme1)
anova(alpha_rot_lme1)


###Can't remove this term but as we limit on hree ways interacion we remove
###Remove Time:Forest:Species:Termi.assum

alpha_rot_lme2<-update(alpha_rot_lme1, .~.-Time:Forest:Species:Termi.assum )
summary(alpha_rot_lme2)
anova(alpha_rot_lme2)
anova(alpha_rot_lme1, alpha_lme2)
AIC(alpha_rot_lme2)

##Can't remove this term but as we limit on hree ways interacion we remove
###No three way interaction
###Remove Time:Forest:Termi.assum

alpha_rot_lme3<-update(alpha_rot_lme2, .~.-Time:Forest:Termi.assum )
summary(alpha_rot_lme3)
anova(alpha_rot_lme3) 
anova(alpha_rot_lme2,alpha_rot_lme3)
AIC(alpha_lme3)

###No three way interaction
###Remove Time:Species:Termi.assum

alpha_rot_lme4<-update(alpha_rot_lme3, .~.-Time:Species:Termi.assum )
summary(alpha_rot_lme4)
anova(alpha_rot_lme4) 
anova(alpha_rot_lme3,alpha_rot_lme4)
AIC(alpha_rot_lme4)

##Go back to model 3
summary(alpha_rot_lme3)
AIC(alpha_rot_lme3)

###No three way interaction
###Remove Time:Forest:Species

alpha_rot_lme5<-update(alpha_rot_lme3, .~.-Time:Forest:Species )
summary(alpha_rot_lme5)
anova(alpha_rot_lme5) 
anova(alpha_rot_lme3,alpha_rot_lme5)
AIC(alpha_rot_lme5)

###Go back to alpha_lme3 

###
###No main effect of Soft Rot then remove
###Remove Soft_Rot

alpha_rot_lme6<-update(alpha_rot_lme3, .~.-Soft_Rot )
summary(alpha_rot_lme6)
anova(alpha_rot_lme6) 
anova(alpha_rot_lme3,alpha_rot_lme6)
AIC(alpha_rot_lme6)

###No main effect of Brown Rot then remove
###Remove Brown_Rot

alpha_rot_lme7<-update(alpha_rot_lme6, .~.-Brown_Rot )
summary(alpha_rot_lme7)
anova(alpha_rot_lme7) 
anova(alpha_rot_lme6,alpha_rot_lme7)
AIC(alpha_rot_lme7)


###Diagnostics
plot(alpha_rot_lme7, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(alpha_rot_lme7, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 


###Use VarExp
alpha_rot_lme7b<-update(alpha_rot_lme7, weights=varExp(form=~fitted(.)))

summary(alpha_rot_lme7b)
anova(alpha_rot_lme7b)

plot(alpha_rot_lme7b, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Slight increasing trend
plot(alpha_rot_lme7b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Not bad



###Use VarPower
alpha_rot_lme7c<-update(alpha_rot_lme7, weights=varPower(form=~fitted(.)))

summary(alpha_rot_lme7c)
anova(alpha_rot_lme7c)

plot(alpha_rot_lme7c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Slight increasing trend
plot(alpha_rot_lme7c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 

AIC(alpha_rot_lme7b,alpha_rot_lme7c)

###Best alpha_lme7b
summary(alpha_rot_lme7b)
r.squaredGLMM(alpha_rot_lme7b)

##> r.squaredGLMM(alpha_lme7b)
######R2m       R2c
##[1,] 0.4384972 0.6975767


####Change forest level

#####Change level to Open
###Relevel  to forest ref="Open_land"

otu_Down_richness_strict_map$Forest
levels(otu_Down_richness_strict_map$Forest)
otu_Down_richness_strict_map$Forest<-relevel(otu_Down_richness_strict_map$Forest, ref="Open_land")
levels(otu_Down_richness_strict_map$Forest)

###Use VarExp
alpha_rot_lme7c<-update(alpha_rot_lme7, weights=varExp(form=~fitted(.)))

##Table S17
summary(alpha_rot_lme7c)
anova(alpha_rot_lme7c)

plot(alpha_rot_lme7c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Slight increasing trend
plot(alpha_rot_lme7c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Not bad

r.squaredGLMM(alpha_rot_lme7c)


####
#####Change level to Mature
###Relevel  to forest ref="Mature_forest"

otu_Down_richness_strict_map$Forest
levels(otu_Down_richness_strict_map$Forest)
otu_Down_richness_strict_map$Forest<-relevel(otu_Down_richness_strict_map$Forest, ref="Mature_forest")
levels(otu_Down_richness_strict_map$Forest)



#####Reduce the model to its simple form including time, forest species and termites

alpha_rot_lme7b


####Remove White Rot
alpha_rot_lme8<-update(alpha_rot_lme7, .~.-White_Rot )
summary(alpha_rot_lme8)
anova(alpha_rot_lme8) 

AIC(alpha_rot_lme8)

####Remove Saprotroph
alpha_rot_lme9<-update(alpha_rot_lme8, .~.-Saprotroph )
summary(alpha_rot_lme9)
anova(alpha_rot_lme9) 
anova(alpha_rot_lme8,alpha_rot_lme9) 
AIC(alpha_rot_lme9)

####Remove Time:Forest:Species
alpha_rot_lme10<-update(alpha_rot_lme9, .~.-Time:Forest:Species )
summary(alpha_rot_lme10)
anova(alpha_rot_lme10) 
anova(alpha_rot_lme9,alpha_rot_lme10) 


###
alpha_rot_lme11<-update(alpha_rot_lme10, .~.-Time:Species:Termi.assum)
summary(alpha_rot_lme11)
anova(alpha_rot_lme11) 
anova(alpha_rot_lme9,alpha_rot_lme10)
##We can't remove Time:Forest:Species
###Go back to alpha_rot_lme9


####
alpha_rot_lme12<-update(alpha_rot_lme11, .~.-Time:Termi.assum)
summary(alpha_rot_lme12)
anova(alpha_rot_lme12) 
anova(alpha_rot_lme11,alpha_rot_lme12)

####
alpha_rot_lme13<-update(alpha_rot_lme12, .~.-Time:Species)
summary(alpha_rot_lme13)
anova(alpha_rot_lme13) 
anova(alpha_rot_lme12,alpha_rot_lme13)

####
alpha_rot_lme14<-update(alpha_rot_lme13, .~.-Time:Forest)
summary(alpha_rot_lme14)
anova(alpha_rot_lme14) 
anova(alpha_rot_lme13,alpha_rot_lme14)

#alpha_rot_lme9<-alpha_rot_lme14
###Diagnotics

plot(alpha_rot_lme9, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(alpha_rot_lme9, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 


###Use VarExp
alpha_rot_lme9b<-update(alpha_rot_lme9, weights=varExp(form=~fitted(.)))

summary(alpha_rot_lme9b)
anova(alpha_rot_lme9b)

plot(alpha_rot_lme9b, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Slight increasing trend
plot(alpha_rot_lme9b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Not bad



###Use VarPower
alpha_rot_lme9c<-update(alpha_rot_lme9, weights=varPower(form=~fitted(.)))

summary(alpha_rot_lme9c)
anova(alpha_rot_lme9c)

plot(alpha_rot_lme9c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Slight increasing trend
plot(alpha_rot_lme9c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 

AIC(alpha_rot_lme9b,alpha_rot_lme9c)

###Best alpha_lme9b
r.squaredGLMM(alpha_rot_lme9b)

summary(alpha_rot_lme9b)

##> r.squaredGLMM(alpha_lme9b)
######R2m       R2c
##[1,] 0.3927463 0.6724915

###Change level to Open land
###Relevel for Forest to ref="Open_land"

otu_Down_richness_strict_map$Forest
levels(otu_Down_richness_strict_map$Forest)
otu_Down_richness_strict_map$Forest<-relevel(otu_Down_richness_strict_map$Forest, ref="Open_land")
levels(otu_Down_richness_strict_map$Forest)
##run again the model

alpha_rot_lme9e<-update(alpha_rot_lme9, weights=varExp(form=~fitted(.)))

summary(alpha_rot_lme9e)
anova(alpha_rot_lme9e)

plot(alpha_rot_lme9e, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Slight increasing trend
plot(alpha_rot_lme9e, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Not bad

r.squaredGLMM(alpha_rot_lme9e)




############
######Rot types modeling

#####Brown rot
hist(otu_Down_richness_strict_map$Brown_Rot)
hist(log(otu_Down_richness_strict_map$Brown_Rot+1))
otu_Down_richness_strict_map$brown_rot_ln<-log(otu_Down_richness_strict_map$Brown_Rot+1)

brown_ln_mod1b<-lme(brown_rot_ln~Time*Species*Forest*Termi.assum*Species, random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map)
summary(brown_ln_mod1b)



###Simplification of most interaction except Time:Forest and Time:Species,Time:Termi.assum, Time:Species, Species:Forest, Time:Species:Forest
brown_ln_mod2b<-lme(brown_rot_ln~Time*Species+Time*Forest+Time*Species*Forest+Forest*Termi.assum,random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map)
summary(brown_ln_mod2b)
anova(brown_ln_mod2b)
### Simplifications
###Diagnostics
plot(brown_ln_mod2b, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(brown_ln_mod2b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 


###Use VarPower
brown_ln_mod2c<-update(brown_ln_mod2b, weights=varPower(form=~fitted(.)))

summary(brown_ln_mod2c)
anova(brown_ln_mod2c)

plot(brown_ln_mod2c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(brown_ln_mod2c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 

AIC(brown_ln_mod2b,brown_ln_mod2c)

###Use VarExp
brown_ln_mod2d<-update(brown_ln_mod2b, weights=varExp(form=~fitted(.)))

summary(brown_ln_mod2d)
anova(brown_ln_mod2d)

plot(brown_ln_mod2d, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend

plot(brown_ln_mod2d, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 

AIC(brown_ln_mod2b,brown_ln_mod2c,brown_ln_mod2d)

###Best model 

brown_ln_mod2d

summary(brown_ln_mod2d)
r.squaredGLMM(brown_ln_mod2d)

#> r.squaredGLMM(brown_ln_mod2d)
#  R2m       R2c
# [1,]  0.4201719 0.5168239

###Use emmeans for pairewise comparison
brown_ln_mod2d

#####Soft rot
#####soft rot
hist(otu_Down_richness_strict_map$Soft_Rot)
hist(log(otu_Down_richness_strict_map$Soft_Rot+1))
otu_Down_richness_strict_map$soft_rot_ln<-log(otu_Down_richness_strict_map$Soft_Rot+1)

soft_ln_mod1b<-lme(soft_rot_ln~Time*Species*Forest*Termi.assum*Species, random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map)
summary(soft_ln_mod1b)



###Simplification of most interaction except Time:Forest and Time:Species
soft_ln_mod2b<-lme(soft_rot_ln~Time+Species+Forest+Termi.assum,random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map)
summary(soft_ln_mod2b)
anova(soft_ln_mod2b)

###Diagnostics
plot(soft_ln_mod2b, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##decreasing trend
plot(soft_ln_mod2b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 


###Use VarPower
soft_ln_mod2c<-update(soft_ln_mod2b, weights=varPower(form=~fitted(.)))

summary(soft_ln_mod2c)
anova(soft_ln_mod2c)

plot(soft_ln_mod2c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)#
plot(soft_ln_mod2c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Perfect trend 

AIC(soft_ln_mod2b,soft_ln_mod2c)

###Use VarExp
soft_ln_mod2d<-update(soft_ln_mod2b, weights=varExp(form=~fitted(.)))

summary(soft_ln_mod2d)
anova(soft_ln_mod2d)

plot(soft_ln_mod2d, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(soft_ln_mod2d, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 

AIC(soft_ln_mod2b,soft_ln_mod2c, soft_ln_mod2d)
###Best soft_ln_mod2d VarExp

summary(soft_ln_mod2d)
r.squaredGLMM(soft_ln_mod2d)

#r.squaredGLMM(soft_ln_mod2d)
######R2m       R2c
#[1,] 0.07281921 0.1353411




###Use linear mixed efffect modeling
####White rot 

hist(otu_Down_richness_strict_map$White_Rot)
hist(log(otu_Down_richness_strict_map$White_Rot+1))
otu_Down_richness_strict_map$white_rot_ln<-log(otu_Down_richness_strict_map$White_Rot+1)

white_ln_mod1b<-lme(white_rot_ln~Time*Species*Forest*Termi.assum*Species, random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map)
summary(white_ln_mod1b)



###Simplification of most interaction except Time:Forest and Time:Species
white_ln_mod2b<-lme(white_rot_ln~Time*Species+Time*Forest+Termi.assum+Species,random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map)
summary(white_ln_mod2b)
anova(white_ln_mod2b)
anova(white_ln_mod1b,white_ln_mod2b)

### Can't remove but since we are limiting to three way interaction

###remove termites
#-Termi.assum
white_ln_mod3b<-update(white_ln_mod2b,.~.-Termi.assum)
summary(white_ln_mod3b)
anova(white_ln_mod3b)
anova(white_ln_mod2b, white_ln_mod3b)


###Diagnostics
plot(white_ln_mod3b, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(white_ln_mod3b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 


###Use VarPower
white_ln_mod3c<-update(white_ln_mod3b, weights=varPower(form=~fitted(.)))

summary(white_ln_mod3c)
anova(white_ln_mod3c)

plot(white_ln_mod3c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(white_ln_mod3c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 

AIC(white_ln_mod3b,white_ln_mod3c)

###Use VarExp
white_ln_mod3d<-update(white_ln_mod3b, weights=varExp(form=~fitted(.)))

summary(white_ln_mod3d)
anova(white_ln_mod3d)

plot(white_ln_mod3d, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(white_ln_mod3d, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 

AIC(white_ln_mod3b,white_ln_mod3c,white_ln_mod3d)

r.squaredGLMM(white_ln_mod3c)

#> r.squaredGLMM(white_ln_mod3c)
#R2m       R2c
#[1,] 0.1003066 0.3678254

summary(white_ln_mod3c)



#########



######Saprotroph abundance
#####saprotroph rot
hist(otu_Down_richness_strict_map$Saprotroph)
hist(log(otu_Down_richness_strict_map$Saprotroph+1))
otu_Down_richness_strict_map$saprotroph_ln<-log(otu_Down_richness_strict_map$Saprotroph+1)

saprotroph_ln_mod1<-lmer(saprotroph_ln~Time*Species*Forest*Termi.assum*Species+(1|Plot/Subplot/Wood_log_ID),data=otu_Down_richness_strict_map)
summary(saprotroph_ln_mod1)



###Simplification of most interaction except Time:Forest and Time:Species, Time:Termites,Time:Species:Termites

saprotroph_ln_mod2<-lmer(saprotroph_ln~Time*Species+Time*Forest+Time*Species*Termi.assum+Time*Termi.assum+(1|Plot/Subplot/Wood_log_ID),data=otu_Down_richness_strict_map)
summary(saprotroph_ln_mod2)
anova(saprotroph_ln_mod1,saprotroph_ln_mod2)


###Simplification Time:Species:Termites

saprotroph_ln_mod3<-update(saprotroph_ln_mod2, .~.-Time:Species:Termi.assum)
summary(saprotroph_ln_mod3)
anova(saprotroph_ln_mod2,saprotroph_ln_mod3 )
AIC(saprotroph_ln_mod2,saprotroph_ln_mod3)

### Remove Species:Termi.assum
saprotroph_ln_mod4<-update(saprotroph_ln_mod3, .~.-Species:Termi.assum)
summary(saprotroph_ln_mod4)
anova(saprotroph_ln_mod3,saprotroph_ln_mod4 )
AIC(saprotroph_ln_mod2,saprotroph_ln_mod3,saprotroph_ln_mod4)

### Remove Time:Termi.assum
saprotroph_ln_mod5<-update(saprotroph_ln_mod4, .~.-Time:Termi.assum)
summary(saprotroph_ln_mod5)
anova(saprotroph_ln_mod3,saprotroph_ln_mod4,saprotroph_ln_mod5 )
AIC(saprotroph_ln_mod2,saprotroph_ln_mod3,saprotroph_ln_mod4,saprotroph_ln_mod5 )

### Remove Termi.assum
saprotroph_ln_mod6<-update(saprotroph_ln_mod5, .~.-Termi.assum)
summary(saprotroph_ln_mod6)
anova(saprotroph_ln_mod3,saprotroph_ln_mod4,saprotroph_ln_mod5,saprotroph_ln_mod6 )
AIC(saprotroph_ln_mod2,saprotroph_ln_mod3,saprotroph_ln_mod4,saprotroph_ln_mod5,saprotroph_ln_mod6 )

####Diagnostic

## Do they show any trend (evidence of nonlinearity)?
## plot residuals against fitted values
strict_resids <- resid(saprotroph_ln_mod6, type='pearson')
plot(strict_resids~fitted(saprotroph_ln_mod6))
lines(lowess(strict_resids~fitted(saprotroph_ln_mod6)), col='red')##A bit trendy (increasing trend)

## are they comparable across the range of predictors?
#boxplot(resids~$Time)
length(strict_resids)
dim(otu_Down_richness_strict_map)

## Are the residuals homoscedastic?
## plot the sqrt of the absolute residuals against fitted values
plot(sqrt(abs(strict_resids))~ fitted(saprotroph_ln_mod6))
lines(lowess(sqrt(abs(strict_resids))~
               fitted(saprotroph_ln_mod6)), col='red')###

## and normally distributed?
qqPlot(strict_resids) ##Need to worry

## But what about the random effects - are they normally distributed?
qqPlot(ranef(saprotroph_ln_mod6)$Plot$'(Intercept)')## This looks fine

## to get an R^2 type statistic (not equivalent to OLS R-square)
library(arm)
library(MuMIn)

r.squaredGLMM(saprotroph_ln_mod3)

#> r.squaredGLMM(saprotroph_ln_mod3)
#R2m       R2c
#[1,] 0.09643187 0.1833403
Anova(saprotroph_ln_mod3)

####Because of heteroscedasticity we used lme
###What if i use lme

saprotroph_ln_mod1b<-lme(saprotroph_ln~Time*Species*Forest*Termi.assum*Species, random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map)
summary(saprotroph_ln_mod1b)



###Simplification of most of the interactions except Time:Forest and Time:Species, Time:Termites,Time:Species:Termites

saprotroph_ln_mod2b<-lme(saprotroph_ln~Time*Species+Time*Forest+Time*Forest*Termi.assum+Time*Termi.assum, random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map)
summary(saprotroph_ln_mod2b)
anova(saprotroph_ln_mod2b)
AIC(saprotroph_ln_mod1b,saprotroph_ln_mod2b)
anova(saprotroph_ln_mod1b,saprotroph_ln_mod2b)

###Can't simplify but we limit to three way interaction
###Simplification Time:Forest:Termi.assum 

saprotroph_ln_mod3b<-update(saprotroph_ln_mod2b, .~.-Time:Forest:Termi.assum )
summary(saprotroph_ln_mod3b)
anova(saprotroph_ln_mod3b )
AIC(saprotroph_ln_mod3b,saprotroph_ln_mod1b,saprotroph_ln_mod2b)
anova(saprotroph_ln_mod2b,saprotroph_ln_mod3b)

### Remove Forest:Termi.assum
saprotroph_ln_mod4b<-update(saprotroph_ln_mod3b, .~.-Forest:Termi.assum)
summary(saprotroph_ln_mod4b)
anova(saprotroph_ln_mod4b )
AIC(saprotroph_ln_mod3b,saprotroph_ln_mod1b,saprotroph_ln_mod2b,saprotroph_ln_mod4b)
anova(saprotroph_ln_mod3b,saprotroph_ln_mod4b)

### Remove Time:Termi.assum
saprotroph_ln_mod5b<-update(saprotroph_ln_mod4b, .~.-Time:Termi.assum)
summary(saprotroph_ln_mod5b)
anova(saprotroph_ln_mod5b)
AIC(saprotroph_ln_mod3b,saprotroph_ln_mod1b,saprotroph_ln_mod2b,saprotroph_ln_mod4b,saprotroph_ln_mod5b)
anova(saprotroph_ln_mod4b,saprotroph_ln_mod5b)

### Remove Time:Forest
saprotroph_ln_mod6b<-update(saprotroph_ln_mod5b, .~.-Time:Forest)
summary(saprotroph_ln_mod6b)
anova(saprotroph_ln_mod6b)
AIC(saprotroph_ln_mod3b,saprotroph_ln_mod1b,saprotroph_ln_mod2b,saprotroph_ln_mod4b,saprotroph_ln_mod5b,saprotroph_ln_mod6b)
anova(saprotroph_ln_mod5b,saprotroph_ln_mod6b)

###Remove Termi.assum

saprotroph_ln_mod7b<-update(saprotroph_ln_mod6b, .~.-Termi.assum)
summary(saprotroph_ln_mod7b)
anova(saprotroph_ln_mod7b)
AIC(saprotroph_ln_mod3b,saprotroph_ln_mod1b,saprotroph_ln_mod2b,saprotroph_ln_mod4b,saprotroph_ln_mod5b,saprotroph_ln_mod6b,saprotroph_ln_mod7b)
anova(saprotroph_ln_mod6b,saprotroph_ln_mod7b)

###Remove Forest

saprotroph_ln_mod8b<-update(saprotroph_ln_mod7b, .~.-Forest)
summary(saprotroph_ln_mod8b)
anova(saprotroph_ln_mod8b)
AIC(saprotroph_ln_mod3b,saprotroph_ln_mod1b,saprotroph_ln_mod2b,saprotroph_ln_mod4b,saprotroph_ln_mod5b,saprotroph_ln_mod6b,saprotroph_ln_mod7b,saprotroph_ln_mod8b)
anova(saprotroph_ln_mod7b,saprotroph_ln_mod8b)

###Diagnostics
plot(saprotroph_ln_mod8b, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Decreasing trend
plot(saprotroph_ln_mod8b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 


###Use VarPower
saprotroph_ln_mod8c<-update(saprotroph_ln_mod8b, weights=varPower(form=~fitted(.)))

summary(saprotroph_ln_mod8c)
anova(saprotroph_ln_mod8c)

plot(saprotroph_ln_mod8c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(saprotroph_ln_mod8c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Very nice one

AIC(saprotroph_ln_mod3b,saprotroph_ln_mod1b,saprotroph_ln_mod2b,saprotroph_ln_mod4b,saprotroph_ln_mod5b,saprotroph_ln_mod6b,saprotroph_ln_mod7b,saprotroph_ln_mod8b,saprotroph_ln_mod8c)
anova(saprotroph_ln_mod8b,saprotroph_ln_mod8c)

###Use VarExp
saprotroph_ln_mod8d<-update(saprotroph_ln_mod8b, weights=varExp(form=~fitted(.)))

summary(saprotroph_ln_mod8d)
anova(saprotroph_ln_mod8d)

plot(saprotroph_ln_mod8d, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(saprotroph_ln_mod8d, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## nice one

AIC(saprotroph_ln_mod8b,saprotroph_ln_mod8c,saprotroph_ln_mod8d)

###Best model saprotroph_ln_mod8d VarExp

saprotroph_ln_mod8d
summary(saprotroph_ln_mod8d)
r.squaredGLMM(saprotroph_ln_mod8d)

#> r.squaredGLMM(saprotroph_ln_mod8d)
#R2m          R2c
#[1,] 0.0001469831 0.0002698264

summary(saprotroph_ln_mod8d)



#####Saprotroph with rot types
####Simple addition of white, brown, soft,rot
saprotroph_ln_mod9b<- lme(saprotroph_ln~Time*Species+Time*Forest+Time*Forest*Termi.assum+Time*Termi.assum+Time*White_Rot+Time*Brown_Rot+Time*Soft_Rot  +Species*White_Rot+Species*Brown_Rot+Species*Soft_Rot
                            +Forest*White_Rot+Forest*Brown_Rot+Forest*Soft_Rot+Termi.assum*White_Rot+Termi.assum*Brown_Rot+Termi.assum*Soft_Rot,random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map)
summary(saprotroph_ln_mod9b)
anova(saprotroph_ln_mod9b)




##Simplification
####Three way ineraction Time:Forest:Termi.assum

saprotroph_ln_mod10b<- update(saprotroph_ln_mod9b,.~.-Time:Forest:Termi.assum)

summary(saprotroph_ln_mod10b)
anova(saprotroph_ln_mod10b)
anova(saprotroph_ln_mod9b,saprotroph_ln_mod10b)

#####Termi.assum:White_Rot
saprotroph_ln_mod11b<- update(saprotroph_ln_mod9b,.~.-Termi.assum:White_Rot)

summary(saprotroph_ln_mod11b)
anova(saprotroph_ln_mod11b)
anova(saprotroph_ln_mod10b,saprotroph_ln_mod11b)

#####Termi.assum:Brown_Rot
saprotroph_ln_mod12b<- update(saprotroph_ln_mod11b,.~.-Termi.assum:Brown_Rot)

summary(saprotroph_ln_mod12b)
anova(saprotroph_ln_mod12b)
anova(saprotroph_ln_mod11b,saprotroph_ln_mod12b)

#####Termi.assum1:Soft_Rot
saprotroph_ln_mod13b<- update(saprotroph_ln_mod12b,.~.-Termi.assum:Soft_Rot)

summary(saprotroph_ln_mod13b)
anova(saprotroph_ln_mod13b)
anova(saprotroph_ln_mod12b,saprotroph_ln_mod13b)


##Forest:Brown_Rot 
saprotroph_ln_mod14b<- update(saprotroph_ln_mod13b,.~.-Forest:Brown_Rot)

summary(saprotroph_ln_mod14b)
anova(saprotroph_ln_mod14b)
anova(saprotroph_ln_mod13b,saprotroph_ln_mod14b)

###Forest:Soft_Rot
saprotroph_ln_mod15b<- update(saprotroph_ln_mod14b,.~.-Forest:Soft_Rot)

summary(saprotroph_ln_mod15b)
anova(saprotroph_ln_mod15b)
anova(saprotroph_ln_mod14b,saprotroph_ln_mod15b)


###Forest:White_Rot
saprotroph_ln_mod16b<- update(saprotroph_ln_mod15b,.~.-Forest:White_Rot)

summary(saprotroph_ln_mod16b)
anova(saprotroph_ln_mod16b)
anova(saprotroph_ln_mod15b,saprotroph_ln_mod16b)


###Species:Brown_Rot 
saprotroph_ln_mod17b<- update(saprotroph_ln_mod16b,.~.-Species:Brown_Rot)

summary(saprotroph_ln_mod17b)
anova(saprotroph_ln_mod17b)
anova(saprotroph_ln_mod16b,saprotroph_ln_mod17b)

###Species:Soft_Rot 
saprotroph_ln_mod18b<- update(saprotroph_ln_mod17b,.~.-Species:Soft_Rot)

summary(saprotroph_ln_mod18b)
anova(saprotroph_ln_mod18b)
anova(saprotroph_ln_mod17b,saprotroph_ln_mod18b)



###Time:Brown_Rot  
saprotroph_ln_mod19b<- update(saprotroph_ln_mod18b,.~.-Time:Brown_Rot)

summary(saprotroph_ln_mod19b)
anova(saprotroph_ln_mod19b)
anova(saprotroph_ln_mod18b,saprotroph_ln_mod19b)

###Time:Soft_Rot  
saprotroph_ln_mod20b<- update(saprotroph_ln_mod19b,.~.-Time:Soft_Rot)

summary(saprotroph_ln_mod20b)
anova(saprotroph_ln_mod20b)
anova(saprotroph_ln_mod19b,saprotroph_ln_mod20b)


###Brown_Rot
saprotroph_ln_mod21b<- update(saprotroph_ln_mod20b,.~.-Brown_Rot)

summary(saprotroph_ln_mod21b)
anova(saprotroph_ln_mod21b)
anova(saprotroph_ln_mod20b,saprotroph_ln_mod21b)

###Diagnostics

plot(saprotroph_ln_mod21b, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Decreasing trend
plot(saprotroph_ln_mod21b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 


###Use VarPower
saprotroph_ln_mod21c<-update(saprotroph_ln_mod21b, weights=varPower(form=~fitted(.)))

summary(saprotroph_ln_mod21c)
anova(saprotroph_ln_mod21c)

plot(saprotroph_ln_mod21c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Slight decreasing trend
plot(saprotroph_ln_mod21c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Very nice one

AIC(saprotroph_ln_mod21b,saprotroph_ln_mod21c)
anova(saprotroph_ln_mod21b,saprotroph_ln_mod21c)

###Use VarExp
saprotroph_ln_mod21d<-update(saprotroph_ln_mod21b, weights=varExp(form=~fitted(.)))

summary(saprotroph_ln_mod21d)
anova(saprotroph_ln_mod21d)

plot(saprotroph_ln_mod21d, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Slight decreasing trend
plot(saprotroph_ln_mod21d, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## nice one

AIC(saprotroph_ln_mod21b,saprotroph_ln_mod21c, saprotroph_ln_mod21d)

###Best model saprotroph_ln_mod7d VarExp

saprotroph_ln_mod21d
###Table S9
summary(saprotroph_ln_mod21d)
r.squaredGLMM(saprotroph_ln_mod21d)

#> r.squaredGLMM(saprotroph_ln_mod21d)
#R2m          R2c
#[1,] 0.0001072938 0.0001239554


###Relevel for Forest to ref="Open_land"

otu_Down_richness_strict_map$Forest
levels(otu_Down_richness_strict_map$Forest)
otu_Down_richness_strict_map$Forest<-relevel(otu_Down_richness_strict_map$Forest, ref="Open_land")
levels(otu_Down_richness_strict_map$Forest)
##run again the model
saprotroph_ln_mod22d<-update(saprotroph_ln_mod21b, weights=varExp(form=~fitted(.)))

summary(saprotroph_ln_mod22d)
anova(saprotroph_ln_mod22d)

plot(saprotroph_ln_mod22d, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Slight decreasing trend
plot(saprotroph_ln_mod22d, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## nice one

###Best model saprotroph_ln_mod22d VarExp

saprotroph_ln_mod22d
###Table S10
summary(saprotroph_ln_mod22d)
r.squaredGLMM(saprotroph_ln_mod22d)


###Relevel back for Forest to ref="Mature_forest"

otu_Down_richness_strict_map$Forest
levels(otu_Down_richness_strict_map$Forest)
otu_Down_richness_strict_map$Forest<-relevel(otu_Down_richness_strict_map$Forest, ref="Mature_forest")
levels(otu_Down_richness_strict_map$Forest)



####Simplify by removing rot type 
summary(saprotroph_ln_mod21d)

saprotroph_ln_mod23b<- update(saprotroph_ln_mod21b,.~.-Soft_Rot-White_Rot-Time:White_Rot-Species:White_Rot)

summary(saprotroph_ln_mod23b)
anova(saprotroph_ln_mod23b)
anova(saprotroph_ln_mod21b,saprotroph_ln_mod23b)


###Remove three ways interaction Time:Forest:Termi.assum
saprotroph_ln_mod24b<- update(saprotroph_ln_mod23b,.~.-Time:Forest:Termi.assum)

summary(saprotroph_ln_mod24b)
anova(saprotroph_ln_mod24b)
anova(saprotroph_ln_mod23b,saprotroph_ln_mod24b)

## Remove Two ways interaction Time:Termi.assum

saprotroph_ln_mod25b<- update(saprotroph_ln_mod24b,.~.-Time:Termi.assum)

summary(saprotroph_ln_mod25b)
anova(saprotroph_ln_mod25b)
anova(saprotroph_ln_mod24b,saprotroph_ln_mod25b)

###Diagnostics

plot(saprotroph_ln_mod25b, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Decreasing trend
plot(saprotroph_ln_mod25b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 


###Use VarPower
saprotroph_ln_mod25c<-update(saprotroph_ln_mod25b, weights=varPower(form=~fitted(.)))

summary(saprotroph_ln_mod25c)
anova(saprotroph_ln_mod25c)

plot(saprotroph_ln_mod25c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Slight decreasing trend
plot(saprotroph_ln_mod25c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Very nice one

AIC(saprotroph_ln_mod25b,saprotroph_ln_mod25c)
anova(saprotroph_ln_mod25b,saprotroph_ln_mod25c)

###Use VarExp
saprotroph_ln_mod25d<-update(saprotroph_ln_mod25b, weights=varExp(form=~fitted(.)))

summary(saprotroph_ln_mod25d)
anova(saprotroph_ln_mod25d)

plot(saprotroph_ln_mod25d, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##
plot(saprotroph_ln_mod25d, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## nice one

AIC(saprotroph_ln_mod25b,saprotroph_ln_mod25c, saprotroph_ln_mod25d)

###Best model saprotroph_ln_mod25d VarExp


####Graphing
preddat_saprotroph <- read.csv("preddat_18_36mo.csv")
head(preddat_saprotroph)
dim(preddat_saprotroph)
str(preddat_saprotroph)
summary(preddat_saprotroph)

summary(saprotroph_ln_mod25d)

anova(saprotroph_ln_mod25d)
preddat_saprotroph
dim(preddat_saprotroph)
###Subset so the should not be termite effect here
###
summary(saprotroph_ln_mod25d)
anova(saprotroph_ln_mod25d)
preddat_saprotroph$pred <- predict(saprotroph_ln_mod25d, newdata=preddat_saprotroph, level=0)

summary(preddat_saprotroph$pred)
###Prediction interval 
## [-2] drops response from formula

Designmat_saprotroph <- model.matrix(formula(saprotroph_ln_mod25d)[-2], preddat_saprotroph)
predvar_saprotroph <- diag(Designmat_saprotroph %*% vcov(saprotroph_ln_mod25d) %*% t(Designmat_saprotroph)) 
preddat_saprotroph$SE <- sqrt(predvar_saprotroph) 
preddat_saprotroph$SE2 <- sqrt(predvar_saprotroph+saprotroph_ln_mod25d$sigma^2) ## sigma= residual 

preddat_saprotroph$SE_uc<-preddat_saprotroph$pred+1.96*preddat_saprotroph$SE
preddat_saprotroph$SE_lc<-preddat_saprotroph$pred-1.96*preddat_saprotroph$SE
head(preddat_saprotroph)


preddat_saprotroph$Species
preddat_saprotroph$Species<-as.factor(preddat_saprotroph$Species)
levels(preddat_saprotroph$Species)
preddat_saprotroph$Time
preddat_saprotroph$Time<-as.factor(preddat_saprotroph$Time)
levels(preddat_saprotroph$Time)
preddat_saprotroph$Time<-factor(preddat_saprotroph$Time, ordered=TRUE)

preddat_saprotroph$Termi.assum
preddat_saprotroph$Termi.assum<-as.factor(preddat_saprotroph$Termi.assum)
levels(preddat_saprotroph$Termi.assum)

levels(preddat_saprotroph$Termi.assum)<- c("Absence","Presence")
levels(preddat_saprotroph$Termi.assum)


preddat_saprotroph$Forest
preddat_saprotroph$Forest<-as.factor(preddat_saprotroph$Forest)
levels(preddat_saprotroph$Forest)
preddat_saprotroph$Forest<-factor(preddat_saprotroph$Forest, levels = c("Mature_forest","Regenerating_forest","Open_land"))
levels(preddat_saprotroph$Forest)
levels(preddat_saprotroph$Forest)<-c("Mature forest","Regenerating forest","Open land")
levels(preddat_saprotroph$Forest)
levels(preddat_saprotroph$Species)
preddat_saprotroph$Species<-factor(preddat_saprotroph$Species,levels=c("Litsea_cubeba","Castanopsis_mekongensis"))
levels(preddat_saprotroph$Species)
levels(preddat_saprotroph$Species)<-c("Litsea cubeba","Castanopsis mekongensis")
levels(preddat_saprotroph$Species)
head(preddat_saprotroph)
summary(preddat_saprotroph)

###Graphing

saprotroph_graph<-ggplot(data=preddat_saprotroph, aes(y=pred, x=Time, colour= Forest,shape=Species, fill=Termi.assum))+
  geom_rect(data=NULL,aes(xmin=1.5,xmax=2.7,ymin=-Inf,ymax=Inf),
            fill="lightgray", colour = NA, alpha = 0.1) +
  scale_shape_manual(values=c(24,21)) +
  scale_fill_manual(values=c(NA, "black"),guide=guide_legend(override.aes=list(shape=21)))+
  geom_point(position=position_dodge(1),stat="identity", size=6)+ 
  geom_errorbar(aes(ymax = SE_uc, ymin = SE_lc), width=0, size = 0.8, position=position_dodge(1))+
  scale_y_continuous(limits=c(7,11))+
  scale_x_discrete(limits=levels(preddat$Time))+
  xlab(" Exposure duration") +
  ylab("Saprotroph fungi abundance (log scale)") + 
  labs(colour="Habitat type")+
  labs(shape="Woody species")+
  labs(fill="Termites status")+
  annotate("text", x=2.3, y=11, label= "Time p-value<0.003          ", size=5)+
  annotate("text", x=2.3, y=10.8, label= "Species p-value=0.026",size=5) +
  annotate("text", x=2.3, y=10.6, label= "Time:Habitat p-value=0.047",size=5)+
  annotate("text", x=2.3, y=10.4, label= " Habitat:Termites p-value=0.036      ",size=5)+
  #annotate("text", x=2.3, y=64, label= "      Time:Habitat:Termites p-value=0.018      ",size=5)+
  #annotate("text", x=2.3, y=62, label= "      Time:Species:Termites p-value=0.070      ",size=5)+
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 20),
        axis.title = element_text(colour = "black", size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        strip.text = element_text(colour = "black", size = 15),
        #legend.position = c("bottom"),
        legend.position=c(0,1), 
        legend.justification=c(0,1),
        #legend.position = c(0.1,0.7),
        legend.text = element_text(colour="black",size=15,face="italic"),
        legend.title=element_text(size=15),
        #legend.text=element_text(size=12,face="italic"),
        #legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm")) 

saprotroph_graph




###
###Model Rot types with other rot type as well

####Brown rot

brown_ln_mod5b<-lme(brown_rot_ln~Time*Species+Time*Forest+Time*Species*Forest+Time*Termi.assum+
                      Time*Soft_Rot+Time*White_Rot +Species*Soft_Rot+Species*White_Rot+
                      Forest*Soft_Rot+Forest*White_Rot+Termi.assum*Soft_Rot+Termi.assum*White_Rot,
                      random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map)
summary(brown_ln_mod5b)
anova(brown_ln_mod5b)

###Remove Forest:White_Rot
brown_ln_mod6b<-update(brown_ln_mod5b,.~.-Forest:White_Rot)
summary(brown_ln_mod6b)
anova(brown_ln_mod6b)
anova(brown_ln_mod5b,brown_ln_mod6b)

###Remove Forest:Soft_Rot
brown_ln_mod7b<-update(brown_ln_mod6b,.~.-Forest:Soft_Rot)
summary(brown_ln_mod7b)
anova(brown_ln_mod7b)
anova(brown_ln_mod6b,brown_ln_mod7b)

###Remove Species:Soft_Rot
brown_ln_mod8b<-update(brown_ln_mod7b,.~.-Species:Soft_Rot)
summary(brown_ln_mod8b)
anova(brown_ln_mod8b)
anova(brown_ln_mod7b,brown_ln_mod8b)

###Remove Species:White_Rot
brown_ln_mod9b<-update(brown_ln_mod8b,.~.-Species:White_Rot)
summary(brown_ln_mod9b)
anova(brown_ln_mod9b)
anova(brown_ln_mod8b,brown_ln_mod9b)


###Remove Time:White_Rot
brown_ln_mod10b<-update(brown_ln_mod9b,.~.-Time:White_Rot)
summary(brown_ln_mod10b)
anova(brown_ln_mod10b)
anova(brown_ln_mod9b,brown_ln_mod10b)

###Remove Time:Soft_Rot
brown_ln_mod11b<-update(brown_ln_mod10b,.~.-Time:Soft_Rot)
summary(brown_ln_mod11b)
anova(brown_ln_mod11b)
anova(brown_ln_mod10b,brown_ln_mod11b)

###Diagnostics
plot(brown_ln_mod11b, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Decreasing trend
plot(brown_ln_mod11b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 


###Use VarPower
#brown_ln_mod11c<-update(brown_ln_mod11b, weights=varPower(form=~fitted(.)))
###Convergene issue
#summary(brown_ln_mod11c)
#anova(brown_ln_mod11c)

#plot(brown_ln_mod11c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
#     abline=0)##
#plot(brown_ln_mod11c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
#     abline=0)##

#AIC(brown_ln_mod11b,brown_ln_mod11c)

###Use VarExp
#brown_ln_mod11d<-update(brown_ln_mod11b, weights=varExp(form=~fitted(.)))

#summary(brown_ln_mod11d)
#anova(brown_ln_mod11d)

#plot(brown_ln_mod11d, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
#     abline=0)#
#plot(brown_ln_mod11d, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
#     abline=0)## #Not bad

#AIC(brown_ln_mod11b,brown_ln_mod11c,brown_ln_mod11d)


###USe VarIdent
###Use VarExp
brown_ln_mod11e<-update(brown_ln_mod11b, weights=varIdent(form=~fitted(.)|Species))

###Table S15
summary(brown_ln_mod11e)
anova(brown_ln_mod11e)

plot(brown_ln_mod11e, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)#
plot(brown_ln_mod11e, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 


###Best model 
r.squaredGLMM(brown_ln_mod11e)

#> r.squaredGLMM(brown_ln_mod11e)
#R2m   R2c
#[1,] 0.158 0.264
summary(brown_ln_mod11e)

###Diagnostics are not that good
###Change level to Open land

otu_Down_richness_strict_map$Forest
levels(otu_Down_richness_strict_map$Forest)
otu_Down_richness_strict_map$Forest<-relevel(otu_Down_richness_strict_map$Forest, ref="Open_land")
levels(otu_Down_richness_strict_map$Forest)

##Run again the model
brown_ln_mod11f<-update(brown_ln_mod11b, weights=varIdent(form=~fitted(.)|Species))

###Table S16
summary(brown_ln_mod11f)
anova(brown_ln_mod11f)

plot(brown_ln_mod11f, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)#
plot(brown_ln_mod11f, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 


###Best model 
r.squaredGLMM(brown_ln_mod11f)

#> r.squaredGLMM(brown_ln_mod11f)
#R2m   R2c
#[1,] 0.158 0.264

###Relevel back for Forest to ref="Mature_forest"

otu_Down_richness_strict_map$Forest
levels(otu_Down_richness_strict_map$Forest)
otu_Down_richness_strict_map$Forest<-relevel(otu_Down_richness_strict_map$Forest, ref="Mature_forest")
levels(otu_Down_richness_strict_map$Forest)




###Simplify to the simle model without other rot type
summary(brown_ln_mod11b)
###Remove Time:White_Rot
brown_ln_mod12b<-update(brown_ln_mod11b,.~.-Soft_Rot-White_Rot-Termi.assum:Soft_Rot-Termi.assum:White_Rot)
summary(brown_ln_mod12b)
anova(brown_ln_mod12b)
anova(brown_ln_mod11b,brown_ln_mod12b)

###Diagnostics
plot(brown_ln_mod12b, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Decreasing trend
plot(brown_ln_mod12b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 


###Use VarPower
brown_ln_mod12c<-update(brown_ln_mod12b, weights=varPower(form=~fitted(.)))
summary(brown_ln_mod12c)
anova(brown_ln_mod12c)

plot(brown_ln_mod12c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##
plot(brown_ln_mod12c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##

AIC(brown_ln_mod12b,brown_ln_mod12c)

###Use VarExp
brown_ln_mod12d<-update(brown_ln_mod12b, weights=varExp(form=~fitted(.)))

summary(brown_ln_mod12d)
anova(brown_ln_mod12d)

plot(brown_ln_mod12d, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)#
plot(brown_ln_mod12d, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## #Not bad

AIC(brown_ln_mod12b,brown_ln_mod12c,brown_ln_mod12d)


#####Graph

#####Prediction

preddat_brown <- read.csv("preddat_18_36mo.csv")
head(preddat_brown)
dim(preddat_brown)
str(preddat_brown)
summary(preddat_brown)

str(preddat_brown)
preddat_brown

###
summary(brown_ln_mod12c)
anova(brown_ln_mod12c)

preddat_brown$pred <- predict(brown_ln_mod12c, newdata=preddat_brown, level=0)

summary(preddat_brown$pred)
###Prediction interval 
## [-2] drops response from formula

Designmat_brown <- model.matrix(formula(brown_ln_mod12c)[-2], preddat_brown)
predvar_brown <- diag(Designmat_brown %*% vcov(brown_ln_mod12c) %*% t(Designmat_brown)) 
preddat_brown$SE <- sqrt(predvar_brown) 
preddat_brown$SE2 <- sqrt(predvar_brown+brown_ln_mod12c$sigma^2) ## sigma= residual 

preddat_brown$SE_uc<-preddat_brown$pred+1.96*preddat_brown$SE
preddat_brown$SE_lc<-preddat_brown$pred-1.96*preddat_brown$SE
head(preddat_brown)


preddat_brown$Species
preddat_brown$Species<-as.factor(preddat_brown$Species)
levels(preddat_brown$Species)
levels(preddat_brown$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(preddat_brown$Species)
preddat_brown$Species<-relevel(preddat_brown$Species, ref="Litsea cubeba")
levels(preddat_brown$Species)
preddat_brown$Time
preddat_brown$Time<-as.factor(preddat_brown$Time)
levels(preddat_brown$Time)
preddat_brown$Time<-factor(preddat_brown$Time, ordered=TRUE)

preddat_brown$Termi.assum
preddat_brown$Termi.assum<-as.factor(preddat_brown$Termi.assum)
levels(preddat_brown$Termi.assum)
levels(preddat_brown$Termi.assum)<-c("Absence", "Presence")
levels(preddat_brown$Termi.assum)
preddat_brown$Forest
preddat_brown$Forest<-as.factor(preddat_brown$Forest)
levels(preddat_brown$Forest)
preddat_brown$Forest<-factor(preddat_brown$Forest, levels = c("Mature_forest","Regenerating_forest","Open_land"))
levels(preddat_brown$Forest)
levels(preddat_brown$Forest)<-c("Mature forest","Regenerating forest","Open land")
levels(preddat_brown$Forest)
head(preddat_brown)
summary(preddat_brown)
###Graphing

brown_graph<-ggplot(data=preddat_brown, aes(y=pred, x=Time, colour= Forest,shape=Species, fill=Termi.assum))+
  geom_rect(data=NULL,aes(xmin=1.5,xmax=2.7,ymin=-Inf,ymax=Inf),
            fill="lightgray", colour = NA, alpha = 0.1) +
  scale_shape_manual(values=c(24,21)) +
  scale_fill_manual(values=c(NA, "black"),guide=guide_legend(override.aes=list(shape=21)))+
  geom_point(position=position_dodge(1),stat="identity", size=6)+ 
  geom_errorbar(aes(ymax = SE_uc, ymin = SE_lc), width=0, size = 0.8, position=position_dodge(1))+
  scale_y_continuous(limits=c(-0.5,7))+
  scale_x_discrete(limits=levels(preddat$Time))+
  xlab(" Exposure duration") +
  ylab("Brown rot fungi abundance (log scale)") + 
  labs(colour="Habitat type")+
  labs(shape="Woody species")+
  labs(fill="Termites status")+
  annotate("text", x=2.3, y=7, label= "Time p-value<0.001          ", size=5)+
  annotate("text", x=2.3, y=6.6, label= "Habitat p-value=0.012",size=5) +
  annotate("text", x=2.3, y=6.2, label= "Termites p-value=0.020",size=5)+
  annotate("text", x=2.3, y=5.8, label= "   Time:Termites p-value=0.043      ",size=5)+
  #annotate("text", x=2.3, y=64, label= "      Time:Habitat:Termites p-value=0.018      ",size=5)+
  #annotate("text", x=2.3, y=62, label= "      Time:Species:Termites p-value=0.070      ",size=5)+
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(colour = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black", size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        strip.text = element_text(colour = "black", size = 15),
        #legend.position = c("bottom"),
        legend.position=c(0,1), 
        legend.justification=c(0,1),
        #legend.position = c(0.1,0.7),
        legend.text = element_text(colour="black",size=15,face="italic"),
        legend.title=element_text(size=15),
        #legend.text=element_text(size=12,face="italic"),
        #legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm")) 

brown_graph




###Soft rot 
soft_ln_mod4b<-lme(soft_rot_ln~Time+Species+Forest+Termi.assum+Time*Brown_Rot+Time*White_Rot +Species*Brown_Rot+Species*White_Rot+
                     Forest*Brown_Rot+Forest*White_Rot+Termi.assum*Brown_Rot+Termi.assum*White_Rot,
                   random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map)
summary(soft_ln_mod4b)
anova(soft_ln_mod4b)
###Remove Forest:White_rot
soft_ln_mod5b<-update(soft_ln_mod4b,.~.-Forest:White_Rot)
summary(soft_ln_mod5b)
anova(soft_ln_mod5b)
anova(soft_ln_mod4b,soft_ln_mod5b)

###Remove Forest:Brown_rot
soft_ln_mod6b<-update(soft_ln_mod5b,.~.-Forest:Brown_Rot)
summary(soft_ln_mod6b)
anova(soft_ln_mod6b)
anova(soft_ln_mod5b,soft_ln_mod6b)


###Remove Species:White_Rot
soft_ln_mod7b<-update(soft_ln_mod6b,.~.-Species:White_Rot)
summary(soft_ln_mod7b)
anova(soft_ln_mod7b)
anova(soft_ln_mod6b,soft_ln_mod7b)

###Remove Species:Brown_Rot
soft_ln_mod8b<-update(soft_ln_mod7b,.~.-Species:Brown_Rot)
summary(soft_ln_mod8b)
anova(soft_ln_mod8b)
anova(soft_ln_mod7b,soft_ln_mod8b)

###Remove Time:White_Rot
soft_ln_mod9b<-update(soft_ln_mod8b,.~.-Time:White_Rot)
summary(soft_ln_mod9b)
anova(soft_ln_mod9b)
anova(soft_ln_mod8b,soft_ln_mod9b)

###Remove Time:Brown_Rot
soft_ln_mod10b<-update(soft_ln_mod9b,.~.-Time:Brown_Rot)
summary(soft_ln_mod10b)
anova(soft_ln_mod10b)
anova(soft_ln_mod9b,soft_ln_mod10b)

###
###Remove Termi.assum:Brown_Rot
soft_ln_mod11b<-update(soft_ln_mod10b,.~.-Termi.assum:Brown_Rot)
summary(soft_ln_mod11b)
anova(soft_ln_mod11b)
anova(soft_ln_mod10b,soft_ln_mod11b)

###Remove Termi.assum:White_Rot
#soft_ln_mod12b<-update(soft_ln_mod11b,.~.-Termi.assum:White_Rot)
#summary(soft_ln_mod12b)
#anova(soft_ln_mod12b)
#anova(soft_ln_mod11b,soft_ln_mod12b)


#####Diagnostics
###Diagnostics
plot(soft_ln_mod11b, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Decreasing trend
plot(soft_ln_mod11b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 


###Use VarPower
soft_ln_mod11c<-update(soft_ln_mod11b, weights=varPower(form=~fitted(.)))

summary(soft_ln_mod11c)
anova(soft_ln_mod11c)

plot(soft_ln_mod11c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##
plot(soft_ln_mod11c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Perfect

AIC(soft_ln_mod11b,soft_ln_mod11c)

###Use VarExp
soft_ln_mod11d<-update(soft_ln_mod11b, weights=varExp(form=~fitted(.)))

summary(soft_ln_mod11d)
anova(soft_ln_mod11d)

plot(soft_ln_mod11d, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)#
plot(soft_ln_mod11d, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## #Not bad

AIC(soft_ln_mod11b,soft_ln_mod11c,soft_ln_mod11d)

###Best model VarExp soft_ln_mod11c
r.squaredGLMM(soft_ln_mod11c)

##Table S13
summary(soft_ln_mod11c)
#> r.squaredGLMM(soft_ln_mod11c)
####R2m       R2c
####[1,] 0.09140351 0.1541454
 

####Change level to open_land
otu_Down_richness_strict_map$Forest
levels(otu_Down_richness_strict_map$Forest)
otu_Down_richness_strict_map$Forest<-relevel(otu_Down_richness_strict_map$Forest, ref="Open_land")
levels(otu_Down_richness_strict_map$Forest)

##Run again the model

soft_ln_mod11e<-update(soft_ln_mod11b, weights=varPower(form=~fitted(.)))
###Table S14
summary(soft_ln_mod11e)
anova(soft_ln_mod11e)

plot(soft_ln_mod11e, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##
plot(soft_ln_mod11e, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Perfect
r.squaredGLMM(soft_ln_mod11e)

####Change back level to Mature forest
otu_Down_richness_strict_map$Forest
levels(otu_Down_richness_strict_map$Forest)
otu_Down_richness_strict_map$Forest<-relevel(otu_Down_richness_strict_map$Forest, ref="Mature_forest")
levels(otu_Down_richness_strict_map$Forest)

####Simplify the other rot type out from the model

###Remove Brown_Rot-White_Rot-Termi.assum:White_Rot
soft_ln_mod12b<-update(soft_ln_mod11b,.~.-Brown_Rot-White_Rot-Termi.assum:White_Rot)
summary(soft_ln_mod12b)
anova(soft_ln_mod12b)
anova(soft_ln_mod11b,soft_ln_mod12b)

#####Diagnostics
###Diagnostics
plot(soft_ln_mod12b, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Decreasing trend
plot(soft_ln_mod12b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 

###Use VarPower
soft_ln_mod12c<-update(soft_ln_mod12b, weights=varPower(form=~fitted(.)))

summary(soft_ln_mod12c)
anova(soft_ln_mod12c)

plot(soft_ln_mod12c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##
plot(soft_ln_mod12c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Almost perfect

AIC(soft_ln_mod12b,soft_ln_mod12c)

###Use VarExp
soft_ln_mod12d<-update(soft_ln_mod12b, weights=varExp(form=~fitted(.)))

summary(soft_ln_mod12d)
anova(soft_ln_mod12d)

plot(soft_ln_mod12d, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)#
plot(soft_ln_mod12d, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## #Not bad

AIC(soft_ln_mod12b,soft_ln_mod12c,soft_ln_mod12d)

#Best model soft_ln_mod12d


#####Graph

#####Prediction

preddat_soft <- read.csv("preddat_18_36mo.csv")
head(preddat_soft)
dim(preddat_soft)
str(preddat_soft)
summary(preddat_soft)

str(preddat_soft)
preddat_soft

###
summary(soft_ln_mod12d)
anova(soft_ln_mod12d)
preddat_soft$pred <- predict(soft_ln_mod12d, newdata=preddat_soft, level=0)

summary(preddat_soft$pred)
###Prediction interval 
## [-2] drops response from formula

Designmat_soft <- model.matrix(formula(soft_ln_mod12d)[-2], preddat_soft)
predvar_soft <- diag(Designmat_soft %*% vcov(soft_ln_mod12d) %*% t(Designmat_soft)) 
preddat_soft$SE <- sqrt(predvar_soft) 
preddat_soft$SE2 <- sqrt(predvar_soft+soft_ln_mod12d$sigma^2) ## sigma= residual 

preddat_soft$SE_uc<-preddat_soft$pred+1.96*preddat_soft$SE
preddat_soft$SE_lc<-preddat_soft$pred-1.96*preddat_soft$SE
head(preddat_soft)


preddat_soft$Species
preddat_soft$Species<-as.factor(preddat_soft$Species)
levels(preddat_soft$Species)
levels(preddat_soft$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
preddat_soft$Species<-relevel(preddat_soft$Species, ref="Litsea cubeba")
levels(preddat_soft$Species)
preddat_soft$Time
preddat_soft$Time<-as.factor(preddat_soft$Time)
levels(preddat_soft$Time)
preddat_soft$Time<-factor(preddat_soft$Time, ordered=TRUE)

preddat_soft$Termi.assum
preddat_soft$Termi.assum<-as.factor(preddat_soft$Termi.assum)
levels(preddat_soft$Termi.assum)
levels(preddat_soft$Termi.assum)<-c("Absence","Presence")
levels(preddat_soft$Termi.assum)
preddat_soft$Forest
preddat_soft$Forest<-as.factor(preddat_soft$Forest)
levels(preddat_soft$Forest)
preddat_soft$Forest<-factor(preddat_soft$Forest, levels = c("Mature_forest","Regenerating_forest","Open_land"))
levels(preddat_soft$Forest)
levels(preddat_soft$Forest)<-c("Mature forest","Regenerating forest","Open land")
levels(preddat_soft$Forest)
head(preddat_soft)
summary(preddat_soft)

###Graphing
summary(preddat_soft)
soft_graph<-ggplot(data=preddat_soft, aes(y=pred, x=Time, colour= Forest,shape=Species, fill=Termi.assum))+
  geom_rect(data=NULL,aes(xmin=1.5,xmax=2.7,ymin=-Inf,ymax=Inf),
            fill="lightgray", colour = NA, alpha = 0.1) +
  scale_shape_manual(values=c(24,21)) +
  scale_fill_manual(values=c(NA, "black"),guide=guide_legend(override.aes=list(shape=21)))+
  geom_point(position=position_dodge(1),stat="identity", size=6)+ 
  geom_errorbar(aes(ymax = SE_uc, ymin = SE_lc), width=0, size = 0.8, position=position_dodge(1))+
  scale_y_continuous(limits=c(2,9))+
  scale_x_discrete(limits=levels(preddat$Time))+
  xlab(" Exposure duration") +
  ylab("Soft rot fungi abundance (log scale)") + 
  labs(colour="Habitat type")+
  labs(shape="Woody species")+
  labs(fill="Termites status")+
  annotate("text", x=2.3, y=9, label= "Time p-value<0.0001          ", size=5)+
  annotate("text", x=2.3, y=8.6, label= "Species p-value<0.0001",size=5) +
  annotate("text", x=2.3, y=8.2, label= "Habitat type p-value=0.020",size=5)+
  annotate("text", x=2.3, y=7.8, label= " Termites p-value<0.0001      ",size=5)+
  #annotate("text", x=2.3, y=64, label= "      Time:Habitat:Termites p-value=0.018      ",size=5)+
  #annotate("text", x=2.3, y=62, label= "      Time:Species:Termites p-value=0.070      ",size=5)+
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(colour = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black", size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        strip.text = element_text(colour = "black", size = 15),
        #legend.position = c("bottom"),
        legend.position=c(0,1), 
        legend.justification=c(0,1),
        #legend.position = c(0.1,0.7),
        legend.text = element_text(colour="black",size=15,face="italic"),
        legend.title=element_text(size=15),
        #legend.text=element_text(size=12,face="italic"),
        #legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm")) 

soft_graph


####White rot
###Simplification of most of interactions except Time:Forest and Time:Species

white_ln_mod5b<-lme(white_rot_ln~Time*Species+Time*Forest+Termi.assum+Species+Time*Brown_Rot+Time*Soft_Rot +Species*Brown_Rot+Species*Soft_Rot+
                      Forest*Brown_Rot+Forest*Soft_Rot+Termi.assum*Brown_Rot+Termi.assum*Soft_Rot,random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map)
summary(white_ln_mod5b)
anova(white_ln_mod5b)

###Remove Forest:Brown_Rot
white_ln_mod6b<-update(white_ln_mod5b,.~.-Forest:Brown_Rot)
summary(white_ln_mod6b)
anova(white_ln_mod6b)

anova(white_ln_mod5b,white_ln_mod6b)

###Remove Forest:Soft_Rot
white_ln_mod7b<-update(white_ln_mod6b,.~.-Forest:Soft_Rot)
summary(white_ln_mod7b)
anova(white_ln_mod7b)

anova(white_ln_mod6b,white_ln_mod7b)

###Remove Species:Brown_Rot
white_ln_mod8b<-update(white_ln_mod7b,.~.-Species:Brown_Rot)
summary(white_ln_mod8b)
anova(white_ln_mod8b)

anova(white_ln_mod7b,white_ln_mod8b)



###Remove Species:Soft_Rot
white_ln_mod9b<-update(white_ln_mod8b,.~.-Species:Soft_Rot)
summary(white_ln_mod9b)
anova(white_ln_mod9b)

anova(white_ln_mod8b,white_ln_mod9b)



###Remove Time:Soft_Rot
white_ln_mod10b<-update(white_ln_mod9b,.~.-Time:Soft_Rot)
summary(white_ln_mod10b)
anova(white_ln_mod10b)

anova(white_ln_mod9b,white_ln_mod10b)

###Diagnostics


###Diagnostics
plot(white_ln_mod10b, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(white_ln_mod10b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 


###Use VarPower
white_ln_mod10c<-update(white_ln_mod10b, weights=varPower(form=~fitted(.)))

summary(white_ln_mod10c)
anova(white_ln_mod10c)

plot(white_ln_mod10c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(white_ln_mod10c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 

AIC(white_ln_mod10b,white_ln_mod10c)

###Use VarExp
white_ln_mod10d<-update(white_ln_mod10b, weights=varExp(form=~fitted(.)))

summary(white_ln_mod10d)
anova(white_ln_mod10d)

plot(white_ln_mod10d, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(white_ln_mod10d, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 

AIC(white_ln_mod10b,white_ln_mod10c,white_ln_mod10d)
###Best model VarExpo
##Table S11
summary(white_ln_mod10d)
r.squaredGLMM(white_ln_mod10d)

#> r.squaredGLMM(white_ln_mod10d)
#R2m       R2c
#[1,] 0.1096027 0.3038954

summary(white_ln_mod10d)

###Relevel for Forest to ref="Open_land"

otu_Down_richness_strict_map$Forest
levels(otu_Down_richness_strict_map$Forest)
otu_Down_richness_strict_map$Forest<-relevel(otu_Down_richness_strict_map$Forest, ref="Open_land")
levels(otu_Down_richness_strict_map$Forest)

###Run the best model
###Use VarPower
white_ln_mod11c<-update(white_ln_mod10b, weights=varPower(form=~fitted(.)))

###Table S12
summary(white_ln_mod11c)
anova(white_ln_mod11c)

plot(white_ln_mod11c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(white_ln_mod11c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 

r.squaredGLMM(white_ln_mod11c)

#> r.squaredGLMM(white_ln_mod11c)
#   R2m       R2c
###[1,] 0.1082098 0.3004967


### Relevel back for FOrest to ref= "Mature_forest"
####Change level to Mature
###Relevel  to forest ref="Mature_forest"

otu_Down_richness_strict_map$Forest
levels(otu_Down_richness_strict_map$Forest)
otu_Down_richness_strict_map$Forest<-relevel(otu_Down_richness_strict_map$Forest, ref="Mature_forest")
levels(otu_Down_richness_strict_map$Forest)


#########Simplify by removing the other rot type

###Remove -Brown_Rot-White_Rot-Time:Brown_Rot
white_ln_mod11b<-update(white_ln_mod10b,.~.-Brown_Rot-Soft_Rot-Termi.assum:Brown_Rot-Termi.assum:Soft_Rot-Time:Brown_Rot)
summary(white_ln_mod11b)
anova(white_ln_mod11b)

anova(white_ln_mod10b,white_ln_mod11b)

###Remove Termi.assum
white_ln_mod12b<-update(white_ln_mod11b,.~.-Termi.assum)
summary(white_ln_mod12b)
anova(white_ln_mod12b)

anova(white_ln_mod11b,white_ln_mod12b)


###Diagnostics


###Diagnostics
plot(white_ln_mod12b, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(white_ln_mod12b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 


###Use VarPower
white_ln_mod12c<-update(white_ln_mod12b, weights=varPower(form=~fitted(.)))

summary(white_ln_mod12c)
anova(white_ln_mod12c)

plot(white_ln_mod12c, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(white_ln_mod12c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 

AIC(white_ln_mod12b,white_ln_mod12c)

###Use VarExp
white_ln_mod12d<-update(white_ln_mod12b, weights=varExp(form=~fitted(.)))

summary(white_ln_mod12d)
anova(white_ln_mod12d)

plot(white_ln_mod12d, abs(resid(., type='pearson'))~fitted(.)|Species+Termi.assum, type=c('p', 'r'), 
     abline=0)##Increasing trend
plot(white_ln_mod12d, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## 

AIC(white_ln_mod12b,white_ln_mod12c,white_ln_mod12d)
#Best based on diagnostics 
white_ln_mod12c

####Graphing
preddat_white <- read.csv("preddat_18_36mo.csv")
head(preddat_white)
dim(preddat_white)
str(preddat_white)
summary(preddat_white)

dim(preddat_white)

summary(white_ln_mod12c)
anova(white_ln_mod12c)
preddat_white$pred <- predict(white_ln_mod12c, newdata=preddat_white, level=0)

summary(preddat_white$pred)
###Prediction interval 
## [-2] drops response from formula

Designmat_white <- model.matrix(formula(white_ln_mod12c)[-2], preddat_white)
predvar_white <- diag(Designmat_white %*% vcov(white_ln_mod12c) %*% t(Designmat_white)) 
preddat_white$SE <- sqrt(predvar_white) 
preddat_white$SE2 <- sqrt(predvar_white+white_ln_mod12c$sigma^2) ## sigma= residual 

preddat_white$SE_uc<-preddat_white$pred+1.96*preddat_white$SE
preddat_white$SE_lc<-preddat_white$pred-1.96*preddat_white$SE
head(preddat_white)


preddat_white$Species
preddat_white$Species<-as.factor(preddat_white$Species)
levels(preddat_white$Species)
preddat_white$Time
preddat_white$Time<-as.factor(preddat_white$Time)
levels(preddat_white$Time)
preddat_white$Time<-factor(preddat_white$Time, ordered=TRUE)

preddat_white$Termi.assum
preddat_white$Termi.assum<-as.factor(preddat_white$Termi.assum)
levels(preddat_white$Termi.assum)
levels(preddat_white$Termi.assum)<-c("Absence", "Presence")
levels(preddat_white$Termi.assum)

preddat_white$Forest
preddat_white$Forest<-as.factor(preddat_white$Forest)
levels(preddat_white$Forest)
preddat_white$Forest<-factor(preddat_white$Forest, levels = c("Mature_forest","Regenerating_forest","Open_land"))
levels(preddat_white$Forest)
levels(preddat_white$Forest)<-c("Mature forest","Regenerating forest","Open land")
levels(preddat_white$Forest)
levels(preddat_white$Species)
preddat_white$Species<-factor(preddat_white$Species,levels=c("Litsea_cubeba","Castanopsis_mekongensis"))
levels(preddat_white$Species)
levels(preddat_white$Species)<-c("Litsea cubeba","Castanopsis mekongensis")
levels(preddat_white$Species)
head(preddat_white)
summary(preddat_white)

###Graphing

white_graph<-ggplot(data=preddat_white, aes(y=pred, x=Time, colour= Forest,shape=Species, fill=Termi.assum))+
  geom_rect(data=NULL,aes(xmin=1.5,xmax=2.7,ymin=-Inf,ymax=Inf),
            fill="lightgray", colour = NA, alpha = 0.1) +
  scale_shape_manual(values=c(24,21)) +
  scale_fill_manual(values=c(NA, "black"),guide=guide_legend(override.aes=list(shape=21)))+
  geom_point(position=position_dodge(1),stat="identity", size=6)+ 
  geom_errorbar(aes(ymax = SE_uc, ymin = SE_lc), width=0, size = 0.8, position=position_dodge(1))+
  scale_y_continuous(limits=c(3,9))+
  scale_x_discrete(limits=levels(preddat$Time))+
  xlab(" Exposure duration") +
  ylab("White rot fungi abundance (log scale)") + 
  labs(colour="Habitat type")+
  labs(shape="Woody species")+
  labs(fill="Termites status")+
  annotate("text", x=2.3, y=9, label= "Time p-value<0.0001          ", size=5)+
  annotate("text", x=2.3, y=8.6, label= "Time:Species p-value=0.002",size=5) +
  annotate("text", x=2.3, y=8.2, label= "Time:Habitat p-value=0.01",size=5)+
  #annotate("text", x=2.3, y=66, label= "      Species:Termites p-value=0.043      ",size=5)+
  #annotate("text", x=2.3, y=64, label= "      Time:Habitat:Termites p-value=0.018      ",size=5)+
  #annotate("text", x=2.3, y=62, label= "      Time:Species:Termites p-value=0.070      ",size=5)+
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 20),
        axis.title = element_text(colour = "black", size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        strip.text = element_text(colour = "black", size = 15),
        #legend.position = c("bottom"),
        legend.position=c(0,1), 
        legend.justification=c(0,1),
        #legend.position = c(0.1,0.7),
        legend.text = element_text(colour="black",size=15,face="italic"),
        legend.title=element_text(size=15),
        #legend.text=element_text(size=12,face="italic"),
        #legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm")) 

white_graph


####Put all together Brown, Soft, White and Saprotroph graphs togeter 4 in 1
tiff(filename="Figure 4 rot saprotroph ABCD.tiff", res = 400, width=8000, height=6500, compression = "lzw")
grid.arrange(arrangeGrob(brown_graph + theme(plot.margin = unit( c(1,0.5,0,0) , units = "lines")),
                         soft_graph + theme(legend.position="none", plot.margin = unit( c(1,0.5,0,0) , units = "lines")),
                         white_graph + theme(legend.position="none", plot.margin = unit( c(1,0.5,0,0) , units = "lines")),
                         saprotroph_graph + theme(legend.position="none", plot.margin = unit( c(1,0.5,0,0) , units = "lines")),
                         nrow=2, ncol=2))

dev.off()



######Graphing for alpha diversity
###Best model is alpha_rot_lme

alpha_rot_lme9b
summary(alpha_rot_lme9b)
anova(alpha_rot_lme9b)


##### Plot with predicted values and confidence interval 
##results
set.seed(2018)
alph_Comb_Down_resultsc <- list()
alph_Comb_Down_resultsc$confint <- intervals(alpha_rot_lme9b,method = 'boot',nsim=999, parallel='snow',ncpus=4,which = "fixed")
alph_Comb_Down_resultsc$r.squared <- r.squaredGLMM(alpha_rot_lme9b)

preddatc <- expand.grid(Species = c("Castanopsis_mekongensis","Litsea_cubeba"),
                        Time = factor(c("18mo","36mo"), ordered=TRUE),
                        Termi.assum=c("0","1"),
                        Forest = c("Mature_forest","Open_land","Regenerating_forest")
                        )

class(preddatc)
dim(preddatc)
str(preddatc)

preddatc<-read.csv(file="preddatc.csv", header=TRUE, sep=",")
dim(preddatc)
str(preddatc)

#####Graph

preddat_alpha_div<-preddatc
str(preddat_alpha_div)

preddat_alpha_div



###
summary(alpha_rot_lme9b)
anova(alpha_rot_lme9b)


preddat_alpha_div$pred <- predict(alpha_rot_lme9b, newdata=preddat_alpha_div, level=0)

summary(preddat_alpha_div$pred)
###Prediction interval 
## [-2] drops response from formula

Designmat_alpha_div <- model.matrix(formula(alpha_rot_lme9b)[-2], preddat_alpha_div)
predvar_alpha_div <- diag(Designmat_alpha_div %*% vcov(alpha_rot_lme9b) %*% t(Designmat_alpha_div)) 
preddat_alpha_div$SE <- sqrt(predvar_alpha_div) 
preddat_alpha_div$SE2 <- sqrt(predvar_alpha_div+alpha_rot_lme9b$sigma^2) ## sigma= residual 

preddat_alpha_div$SE_uc<-preddat_alpha_div$pred+1.96*preddat_alpha_div$SE
preddat_alpha_div$SE_lc<-preddat_alpha_div$pred-1.96*preddat_alpha_div$SE
head(preddat_alpha_div)
tail(preddat_alpha_div)

#convert to original scale from sqrt
preddat_alpha_div$pred_orig <- (preddat_alpha_div$pred)^2
preddat_alpha_div$SE_uc_all_orig <- (preddat_alpha_div$SE_uc)^2
preddat_alpha_div$SE_lc_all_orig <- (preddat_alpha_div$SE_lc)^2
head(preddat_alpha_div)

preddat_alpha_div$Species
preddat_alpha_div$Species<-as.factor(preddat_alpha_div$Species)
levels(preddat_alpha_div$Species)
preddat_alpha_div$Time
preddat_alpha_div$Time<-as.factor(preddat_alpha_div$Time)
levels(preddat_alpha_div$Time)
preddat_alpha_div$Time<-factor(preddat_alpha_div$Time, ordered=TRUE)

preddat_alpha_div$Termi.assum

preddat_alpha_div$Forest
preddat_alpha_div$Forest<-as.factor(preddat_alpha_div$Forest)
levels(preddat_alpha_div$Forest)
preddat_alpha_div$Forest<-factor(preddat_alpha_div$Forest, levels = c("Mature_forest","Regenerating_forest","Open_land"))
levels(preddat_alpha_div$Forest)
levels(preddat_alpha_div$Forest)<- c("Mature_forest","Regenerating_forest","Open_land") 
levels(preddat_alpha_div$Forest)
head(preddat_alpha_div)
summary(preddat_alpha_div)
dim(preddat_alpha_div)

alph_Comb_Down_resultsc$model.predictions <- preddat_alpha_div


###Write a csv file of the model predictions and add the time 0mo information
###mean(SE)
# Castanopsis 69.17099/sqrt(23) # 14.42315 # Therefore mean(SE): 178.6013(14.42315)
# Litsea 154.68577/sqrt(22) #  32.97912 # Therefore mean(SE): 275.3904(32.97912)
#time0mo<-data.frame(Time=rep("0mo",6),Forest= rep(c("Mature_forest","Regenerating_forest","Open_land"),2),Species=c(rep("Castanopsis_mekongensis",3),rep("Litsea_cubeba",3)),
#                    Termi.assum=c(rep("0",3),rep("1",3)),Down_strict_alpha_model8b.pred=c(rep(178.6013,3),rep(275.3904,3)),Down_strict_alpha_model8b.uc.all=c(rep(178.6013+14.42315,3),rep(275.3904+32.97912,3)), 
#                    Down_strict_alpha_model8b.lc.all=c(rep(178.6013-14.42315,3),rep(275.3904-32.97912,3)))
#dim(time0mo)
#time0mo
#dim(alph_Comb_Down_resultsb$model.predictions)

time0mob<-data.frame(Time=rep("0mo",12),Forest= rep(rep(c("Mature_forest","Regenerating_forest","Open_land"),2),2),Species=rep(c(rep("Castanopsis_mekongensis",3),rep("Litsea_cubeba",3)),2),
                     Termi.assum=c(rep("0",6),rep("1",6)),pred_orig=rep(c(rep(178.6013,3),rep(275.3904,3)),2),SE_uc_all_orig=rep(c(rep(178.6013+14.42315,3),rep(275.3904+32.97912,3)),2), 
                     SE_lc_all_orig=rep(c(rep(178.6013-14.42315,3),rep(275.3904-32.97912,3)),2))
dim(time0mob)
time0mob
head(time0mob)


####Because we don't have cases where termites were present at time 0 then we can't input this information
time0mob<-time0mob[time0mob$Termi.assum=="0",]
dim(alph_Comb_Down_resultsc$model.predictions)
head(alph_Comb_Down_resultsc$model.predictions)
head(time0mob)

head(alph_Comb_Down_resultsc$model.predictions)
summary(alph_Comb_Down_resultsc$model.predictions)
head(time0mob)
alph_Comb_Down_resultsc$model.predictions[,c(2,4,1,3,10:12)]

alph_Comb_Down_resultsc$model.predictions<-rbind(time0mob,alph_Comb_Down_resultsc$model.predictions[,c(2,4,1,3,10:12)])


alph_Comb_Down_resultsc
alph_Comb_Down_resultsc$confint
alph_Comb_Down_resultsc$r.squared
alph_Comb_Down_resultsc$model.predictions
summary(alph_Comb_Down_resultsc$model.predictions)

#### Graphing 
levels(alph_Comb_Down_resultsc$model.predictions$Time)
alph_Comb_Down_resultsc$model.predictions$Time<-as.factor(alph_Comb_Down_resultsc$model.predictions$Time)
levels(alph_Comb_Down_resultsc$model.predictions$Time)
levels(alph_Comb_Down_resultsc$model.predictions$Time)<-c("0 mo", "18 mo", "36 mo")
levels(alph_Comb_Down_resultsc$model.predictions$Forest)
alph_Comb_Down_resultsc$model.predictions$Forest<-as.factor(alph_Comb_Down_resultsc$model.predictions$Forest)
levels(alph_Comb_Down_resultsc$model.predictions$Forest)
levels(alph_Comb_Down_resultsc$model.predictions$Forest)<-c("Mature forest", "Open land", "Regenerating forest")
alph_Comb_Down_resultsc$model.predictions$Forest<-factor(alph_Comb_Down_resultsc$model.predictions$Forest, levels=c("Mature forest","Regenerating forest","Open land"))
levels(alph_Comb_Down_resultsc$model.predictions$Forest)
levels(alph_Comb_Down_resultsc$model.predictions$Termi.assum)
alph_Comb_Down_resultsc$model.predictions$Termi.assum<-as.factor(alph_Comb_Down_resultsc$model.predictions$Termi.assum)
levels(alph_Comb_Down_resultsc$model.predictions$Termi.assum)
levels(alph_Comb_Down_resultsc$model.predictions$Termi.assum)<-c("Absence", "Presence")
levels(alph_Comb_Down_resultsc$model.predictions$Termi.assum)

levels(alph_Comb_Down_resultsc$model.predictions$Species)
alph_Comb_Down_resultsc$model.predictions$Species<-as.factor(alph_Comb_Down_resultsc$model.predictions$Species)
levels(alph_Comb_Down_resultsc$model.predictions$Species)
alph_Comb_Down_resultsc$model.predictions$Species<-factor(alph_Comb_Down_resultsc$model.predictions$Species, levels=c("Litsea_cubeba","Castanopsis_mekongensis"))
levels(alph_Comb_Down_resultsc$model.predictions$Species)
levels(alph_Comb_Down_resultsc$model.predictions$Species)<-c("Litsea cubeba","Castanopsis mekongensis")
levels(alph_Comb_Down_resultsc$model.predictions$Species)



alpha_diversity_graph_December_2020<-ggplot(data=alph_Comb_Down_resultsc$model.predictions, aes(x=Time, y=pred_orig, colour= Forest,shape=Species, fill=Termi.assum))+
  geom_point(position=position_dodge(1),stat="identity", size=6)+ 
  geom_rect(data=NULL,aes(xmin=1.5,xmax=2.5,ymin=-Inf,ymax=Inf),
            fill="lightgray", colour = NA, alpha = 0.025) +
  scale_shape_manual(values=c(24,21)) +
  scale_fill_manual(values=c(NA, "black"),guide=guide_legend(override.aes=list(shape=21)))+
  geom_errorbar(aes(ymax = SE_uc_all_orig, ymin = SE_lc_all_orig), width=0, size = 0.8, position=position_dodge(1))+
  scale_y_continuous(limits=c(20,520))+
  labs(colour="Habitat type")+
  labs(shape="Woody species")+
  labs(fill="Termites status")+
  annotate("text", x=3.1, y=500, label= "Time p-value<0.001", size=5)+
  annotate("text", x=3.1, y=480, label= "Species p-value<0.001",size=5)+
  annotate("text", x=3, y=460, label= "  Habitat:Species:Termites p-value=0.016      ",size=5)+
  annotate("text", x=3, y=440, label= "Time:Species:Termites p-value=0.032      ",size=5)+
  ylab(" Alpha diversity (fungal OTUs, Chao 1)  ")+
  xlab("Exposure time (months)")+
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        strip.text = element_text(colour = "black", size = 15),
        #legend.position = c("bottom"),
        legend.position=c(0,1), 
        legend.justification=c(0,1),
        legend.text = element_text(colour="black",size=15,face="italic"),
        legend.title=element_text(size=15),
        #legend.text=element_text(size=12,face="italic"),
        #legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm"))

alpha_diversity_graph_December_2020



tiff(filename="Figure 5 Alpha diversity.tiff", res = 800, width=8750, height=7600, compression = "lzw")
alpha_diversity_graph_December_2020
dev.off()


#####Put all lme in the psem
library(piecewiseSEM)
###

#WSG_lme_all38b

#alpha_lme3b

#brown_ln_mod2c

#soft_ln_mod2d

#white_ln_mod3c

#saprotroph_ln_mod21d



####SEM


otu_Down_richness_strict_map2<-mutate(otu_Down_richness_strict_map,brown_rot_ln=log(Brown_Rot+1),soft_rot_ln=log(Soft_Rot+1),white_rot_ln=log(White_Rot+1),
                                      saprotroph_ln=log(Saprotroph+1))

head(otu_Down_richness_strict_map2)
str(otu_Down_richness_strict_map2)
summary(otu_Down_richness_strict_map2)

otu_Down_richness_strict_map2$Forest
otu_Down_richness_strict_map2$Time
otu_Down_richness_strict_map2$Time<-factor(otu_Down_richness_strict_map2$Time, ordered=FALSE)
otu_Down_richness_strict_map2$Time
otu_Down_richness_strict_map2$Termi.assum
otu_Down_richness_strict_map2$Species
otu_Down_richness_strict_map2$Species<-as.factor(otu_Down_richness_strict_map2$Species)
otu_Down_richness_strict_map2$Species
#<-as.factor(otu_Down_richness_strict_map2$Time)

#write.csv(otu_Down_richness_strict_map2,file='psem.csv')

#otu_Down_richness_strict_map2<-read.csv(file="psem.csv", header=T, sep=",")
##WSG #WSG_lme_all38b
library(piecewiseSEM)

dim(otu_Down_richness_strict_map2)

Dossa_psem2<- psem(lme(Per_WSG_loss ~ Time + Chao1_strict_sqrt + Species + Forest +Termi.assum + white_rot_ln + brown_rot_ln + soft_rot_ln + saprotroph_ln, 
                      random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                  ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                  lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                  ##Brown rot brown_ln_mod2c
                  lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                  ##Softrot soft_ln_mod2d
                  lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                  ##White rot white_ln_mod3c
                  lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                  ###Saprotrophs saprotroph_ln_mod7d
                  lme(saprotroph_ln ~ Time + Species + Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                  data=otu_Down_richness_strict_map2
)


dSep(Dossa_psem2)

###Suggesting to include some missing links
###Base on biological meaningfulness we included
### rot types on Chao1
### saprotroph on Chao1
### rot types on saprotroph
### rot type on other rot type

Dossa_psem3<- psem(lme(Per_WSG_loss ~ Time + Chao1_strict_sqrt + Species + Forest +Termi.assum + white_rot_ln + brown_rot_ln + soft_rot_ln + saprotroph_ln, 
                       random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                   lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum+white_rot_ln + brown_rot_ln + soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Brown rot brown_ln_mod2c
                   lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Softrot soft_ln_mod2d
                   lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ##White rot white_ln_mod3c
                   lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ###Saprotrophs saprotroph_ln_mod7d
                   lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                   data=otu_Down_richness_strict_map2
)


dSep(Dossa_psem3)
### There are other suggestion but none of them is showing a significant P value so it is fine
#Dossa_psem3
summary(Dossa_psem3)

anova(Dossa_psem2,Dossa_psem3)

###Model simplification on WSG loss
####Remove Chao1_strict_sqrt on WSG


Dossa_psem4<- psem(lme(Per_WSG_loss ~ Time + Species + Forest +Termi.assum + white_rot_ln + brown_rot_ln + soft_rot_ln + saprotroph_ln, 
                       random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                   lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum+white_rot_ln + brown_rot_ln + soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Brown rot brown_ln_mod2c
                   lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Softrot soft_ln_mod2d
                   lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ##White rot white_ln_mod3c
                   lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ###Saprotrophs saprotroph_ln_mod7d
                   lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                   data=otu_Down_richness_strict_map2
)



summary(Dossa_psem4)

###Remove white_rot_ln on WSG loss
Dossa_psem5<- psem(lme(Per_WSG_loss ~ Time + Species + Forest +Termi.assum + brown_rot_ln + soft_rot_ln + saprotroph_ln, 
                       random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                   lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum+white_rot_ln + brown_rot_ln + soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Brown rot brown_ln_mod2c
                   lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Softrot soft_ln_mod2d
                   lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ##White rot white_ln_mod3c
                   lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ###Saprotrophs saprotroph_ln_mod7d
                   lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                   data=otu_Down_richness_strict_map2
)



summary(Dossa_psem5)




###Remove brown_rot_ln on WSG loss

Dossa_psem6<- psem(lme(Per_WSG_loss ~ Time + Species + Forest +Termi.assum + soft_rot_ln + saprotroph_ln, 
                       random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                   lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum+white_rot_ln + brown_rot_ln + soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Brown rot brown_ln_mod2c
                   lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Softrot soft_ln_mod2d
                   lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ##White rot white_ln_mod3c
                   lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ###Saprotrophs saprotroph_ln_mod7d
                   lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                   data=otu_Down_richness_strict_map2
)



summary(Dossa_psem6)



####Remove Forest on WSG loss


Dossa_psem7<- psem(lme(Per_WSG_loss ~ Time + Species +Termi.assum + soft_rot_ln + saprotroph_ln, 
                       random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                   lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum+white_rot_ln + brown_rot_ln + soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Brown rot brown_ln_mod2c
                   lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Softrot soft_ln_mod2d
                   lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ##White rot white_ln_mod3c
                   lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ###Saprotrophs saprotroph_ln_mod7d
                   lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                   data=otu_Down_richness_strict_map2
)



summary(Dossa_psem7)

####Remove Termi.assum

Dossa_psem8<- psem(lme(Per_WSG_loss ~ Time + Species + soft_rot_ln + saprotroph_ln, 
                       random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                   lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum+white_rot_ln + brown_rot_ln + soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Brown rot brown_ln_mod2c
                   lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Softrot soft_ln_mod2d
                   lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ##White rot white_ln_mod3c
                   lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ###Saprotrophs saprotroph_ln_mod7d
                   lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                   data=otu_Down_richness_strict_map2
)



summary(Dossa_psem8)
anova(Dossa_psem7,Dossa_psem8)

###Remove Forest on Chao 
Dossa_psem9<- psem(lme(Per_WSG_loss ~ Time + Species + soft_rot_ln + saprotroph_ln, 
                       random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                   lme(Chao1_strict_sqrt ~ Time+Species+Termi.assum+white_rot_ln + brown_rot_ln + soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Brown rot brown_ln_mod2c
                   lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Softrot soft_ln_mod2d
                   lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ##White rot white_ln_mod3c
                   lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ###Saprotrophs saprotroph_ln_mod7d
                   lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                   data=otu_Down_richness_strict_map2
)

summary(Dossa_psem9)

###Remove white_rot_ln on Chao 
Dossa_psem10<- psem(lme(Per_WSG_loss ~ Time + Species + soft_rot_ln + saprotroph_ln, 
                       random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                   lme(Chao1_strict_sqrt ~ Time+Species+Termi.assum + brown_rot_ln + soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Brown rot brown_ln_mod2c
                   lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Softrot soft_ln_mod2d
                   lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ##White rot white_ln_mod3c
                   lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ###Saprotrophs saprotroph_ln_mod7d
                   lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                   data=otu_Down_richness_strict_map2
)

summary(Dossa_psem10)
anova(Dossa_psem9,Dossa_psem10)



###Remove white_rot_ln on Chao 
Dossa_psem11<- psem(lme(Per_WSG_loss ~ Time + Species + soft_rot_ln + saprotroph_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+ brown_rot_ln + soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem11)
anova(Dossa_psem10,Dossa_psem11)


###Remove brown_rot_ln on Chao 
Dossa_psem12<- psem(lme(Per_WSG_loss ~ Time + Species + soft_rot_ln + saprotroph_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+  soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem12)
anova(Dossa_psem11,Dossa_psem12)


###Remove Termi.assum on Chao 
Dossa_psem13<- psem(lme(Per_WSG_loss ~ Time + Species + soft_rot_ln + saprotroph_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+  soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem13)
anova(Dossa_psem12,Dossa_psem13)


###Remove Termi.assum on brown_rot_ln 
Dossa_psem14<- psem(lme(Per_WSG_loss ~ Time + Species + soft_rot_ln + saprotroph_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+  soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem14)
anova(Dossa_psem13,Dossa_psem14)




###Remove Species on white_rot_ln
Dossa_psem15<- psem(lme(Per_WSG_loss ~ Time + Species + soft_rot_ln + saprotroph_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+  soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem15)
anova(Dossa_psem14,Dossa_psem15)

###Remove Forest on white_rot_ln
Dossa_psem16<- psem(lme(Per_WSG_loss ~ Time + Species + soft_rot_ln + saprotroph_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+  soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem16)
anova(Dossa_psem15,Dossa_psem16)

#Remove Forest on saprotroph_ln
Dossa_psem17<- psem(lme(Per_WSG_loss ~ Time + Species + soft_rot_ln + saprotroph_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+  soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time + Species +white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem17)
anova(Dossa_psem16,Dossa_psem17)

#Remove Species on saprotroph_ln
Dossa_psem18<- psem(lme(Per_WSG_loss ~ Time + Species + soft_rot_ln + saprotroph_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+  soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time +white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem18)
anova(Dossa_psem17,Dossa_psem18)


#Remove Brown rot on saprotroph_ln
Dossa_psem19<- psem(lme(Per_WSG_loss ~ Time + Species + soft_rot_ln + saprotroph_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+  soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time +white_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem19)
anova(Dossa_psem18,Dossa_psem19)


#Include the two missing links fro dSep


Dossa_psem20<- psem(lme(Per_WSG_loss ~ Time + Species + soft_rot_ln + saprotroph_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+  soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time +white_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem20)
SEM<-summary(Dossa_psem20)
str(summary(Dossa_psem20))

###Table S18
write.csv(SEM$coefficients, file="SEM piecewise SEM model.csv")
anova(Dossa_psem19,Dossa_psem20)
###

#####Do marginal means
library(emmeans)

library(multcomp)
### WSG loss
emmeans(Dossa_psem20[[1]],specs=pairwise~"Time")
#> emmeans(Dossa_psem20[[1]],specs=pairwise~"Time")
#$emmeans
#Time emmean   SE df lower.CL upper.CL
#18mo   25.6 1.15 26     23.2     27.9
#36mo   42.3 1.26 26     39.7     44.9

#Results are averaged over the levels of: Species 
#Degrees-of-freedom method: containment 
#Confidence level used: 0.95 

#$contrasts
#contrast    estimate   SE  df t.ratio p.value
#18mo - 36mo    -16.7 1.39 160 -12.043 <.0001 

#Results are averaged over the levels of: Species 
#Degrees-of-freedom method: containment 



emmeans(Dossa_psem20[[1]],specs=pairwise~"Species")

#> emmeans(Dossa_psem20[[1]],specs=pairwise~"Species")
#$emmeans
#Species                 emmean   SE df lower.CL upper.CL
#Castanopsis_mekongensis   34.9 1.24 26     32.4     37.5
#Litsea_cubeba             33.0 1.29 26     30.3     35.6

#Results are averaged over the levels of: Time 
#Degrees-of-freedom method: containment 
#Confidence level used: 0.95 

#$contrasts
#contrast                                estimate   SE  df t.ratio p.value
#Castanopsis_mekongensis - Litsea_cubeba     1.95 1.56 104 1.250   0.2141 

#Results are averaged over the levels of: Time 
#Degrees-of-freedom method: containment 



###Chao squareroot
emmeans(Dossa_psem20[[2]],specs=pairwise~"Time")

#> emmeans(Dossa_psem20[[2]],specs=pairwise~"Time")
#$emmeans
#Time emmean    SE df lower.CL upper.CL
#18mo   16.2 0.248 26     15.6     16.7
#36mo   14.4 0.279 26     13.8     14.9

#Results are averaged over the levels of: Species 
#Degrees-of-freedom method: containment 
#Confidence level used: 0.95 

#$contrasts
#contrast    estimate    SE  df t.ratio p.value
#18mo - 36mo     1.78 0.355 160 5.013   <.0001 

#Results are averaged over the levels of: Species 
#Degrees-of-freedom method: containment 

#emmeans(Dossa_psem20[[2]],specs=pairwise~"Species")

#> emmeans(Dossa_psem20[[2]],specs=pairwise~"Species")
#$emmeans
#Species                 emmean    SE df lower.CL upper.CL
#Castanopsis_mekongensis   14.6 0.254 26     14.1     15.1
#Litsea_cubeba             15.9 0.267 26     15.4     16.5

#Results are averaged over the levels of: Time 
#Degrees-of-freedom method: containment 
#Confidence level used: 0.95 

#$contrasts
#contrast                                estimate    SE  df t.ratio p.value
#Castanopsis_mekongensis - Litsea_cubeba    -1.35 0.347 104 -3.893  0.0002 

#Results are averaged over the levels of: Time 
#Degrees-of-freedom method: containment



###Brown rot
emmeans(Dossa_psem20[[3]],specs=pairwise~"Time")
#> emmeans(Dossa_psem20[[3]],specs=pairwise~"Time")
#$emmeans
#Time emmean    SE df lower.CL upper.CL
#18mo   1.14 0.161 24    0.803     1.47
#36mo   1.78 0.175 24    1.417     2.14

#Results are averaged over the levels of: Species, Forest 
#Degrees-of-freedom method: containment 
#Confidence level used: 0.95 

#$contrasts
#contrast    estimate    SE  df t.ratio p.value
#18mo - 36mo   -0.642 0.192 160 -3.350  0.0010 

#Results are averaged over the levels of: Species, Forest 
#Degrees-of-freedom method: containment 


emmeans(Dossa_psem20[[3]],specs=pairwise~"Forest")
#> emmeans(Dossa_psem20[[3]],specs=pairwise~"Forest")
#$emmeans
#Forest              emmean    SE df lower.CL upper.CL
#Mature_forest         1.34 0.212 26    0.908     1.78
#Open_land             1.98 0.303 24    1.359     2.61
#Regenerating_forest   1.04 0.190 24    0.651     1.44

#Results are averaged over the levels of: Time, Species 
#Degrees-of-freedom method: containment 
#Confidence level used: 0.95 

#$contrasts
#contrast                            estimate    SE df t.ratio p.value
#Mature_forest - Open_land             -0.641 0.372 24 -1.724  0.2170 
#Mature_forest - Regenerating_forest    0.299 0.284 24  1.051  0.5530 
#Open_land - Regenerating_forest        0.940 0.357 24  2.633  0.0374 

#Results are averaged over the levels of: Time, Species 
#Degrees-of-freedom method: containment 
#P value adjustment: tukey method for comparing a family of 3 estimates 


emmeans(Dossa_psem20[[3]],specs=pairwise~"Species")

#> emmeans(Dossa_psem20[[3]],specs=pairwise~"Species")
#$emmeans
#Species                 emmean    SE df lower.CL upper.CL
#Castanopsis_mekongensis   1.68 0.163 24    1.341     2.01
#Litsea_cubeba             1.24 0.171 24    0.885     1.59

#Results are averaged over the levels of: Time, Forest 
#Degrees-of-freedom method: containment 
#Confidence level used: 0.95 

#$contrasts
#contrast                                estimate    SE  df t.ratio p.value
#Castanopsis_mekongensis - Litsea_cubeba     0.44 0.187 104 2.357   0.0203 

#Results are averaged over the levels of: Time, Forest 
#Degrees-of-freedom method: containment 

###Soft rot
emmeans(Dossa_psem20[[4]],specs=pairwise~"Time")

#> emmeans(Dossa_psem20[[4]],specs=pairwise~"Time")
#$emmeans
#Time emmean    SE df lower.CL upper.CL
#18mo   5.54 0.169 24     5.19     5.89
#36mo   4.76 0.178 24     4.40     5.13

#Results are averaged over the levels of: Species, Forest, Termi.assum 
#Degrees-of-freedom method: containment 
#Confidence level used: 0.95 

#$contrasts
#contrast    estimate    SE  df t.ratio p.value
#18mo - 36mo    0.776 0.177 161 4.390   <.0001 

#Results are averaged over the levels of: Species, Forest, Termi.assum 
#Degrees-of-freedom method: containment 


emmeans(Dossa_psem20[[4]],specs=pairwise~"Species")
#> emmeans(Dossa_psem20[[4]],specs=pairwise~"Species")
#$emmeans
#Species                 emmean    SE df lower.CL upper.CL
#Castanopsis_mekongensis   4.79 0.173 24     4.44     5.15
#Litsea_cubeba             5.51 0.179 24     5.14     5.88

#Results are averaged over the levels of: Time, Forest, Termi.assum 
#Degrees-of-freedom method: containment 
#Confidence level used: 0.95 

#$contrasts
#contrast                                estimate    SE  df t.ratio p.value
#Castanopsis_mekongensis - Litsea_cubeba   -0.718 0.185 104 -3.878  0.0002 

#Results are averaged over the levels of: Time, Forest, Termi.assum 
#Degrees-of-freedom method: containment 



emmeans(Dossa_psem20[[4]],specs=pairwise~"Forest")
#> emmeans(Dossa_psem20[[4]],specs=pairwise~"Forest")
#$emmeans
#Forest              emmean    SE df lower.CL upper.CL
#Mature_forest         4.65 0.228 26     4.18     5.12
#Open_land             5.77 0.321 24     5.10     6.43
#Regenerating_forest   5.03 0.204 24     4.61     5.45

#Results are averaged over the levels of: Time, Species, Termi.assum 
#Degrees-of-freedom method: containment 
#Confidence level used: 0.95 

#$contrasts
#contrast                            estimate    SE df t.ratio p.value
#Mature_forest - Open_land             -1.113 0.391 24 -2.846  0.0234 
#Mature_forest - Regenerating_forest   -0.380 0.304 24 -1.252  0.4355 
#Open_land - Regenerating_forest        0.733 0.379 24  1.936  0.1505 

#Results are averaged over the levels of: Time, Species, Termi.assum 
#Degrees-of-freedom method: containment 
#P value adjustment: tukey method for comparing a family of 3 estimates 

emmeans(Dossa_psem20[[4]],specs=pairwise~"Termi.assum")


#> emmeans(Dossa_psem20[[4]],specs=pairwise~"Termi.assum")
#$emmeans
#Termi.assum emmean    SE df lower.CL upper.CL
#0             4.73 0.170 24     4.38     5.09
#1             5.57 0.198 24     5.16     5.98

#Results are averaged over the levels of: Time, Species, Forest 
#Degrees-of-freedom method: containment 
#Confidence level used: 0.95 

#$contrasts
#contrast estimate    SE  df t.ratio p.value
#0 - 1      -0.833 0.218 161 -3.825  0.0002 

#Results are averaged over the levels of: Time, Species, Forest 
#Degrees-of-freedom method: containment


###White rot
emmeans(Dossa_psem20[[5]],specs=pairwise~"Time")
#> emmeans(Dossa_psem20[[5]],specs=pairwise~"Time")
#$emmeans
#Time emmean    SE df lower.CL upper.CL
#18mo   4.84 0.179 26     4.47     5.20
#36mo   5.80 0.197 26     5.40     6.21

#Degrees-of-freedom method: containment 
#Confidence level used: 0.95 

#$contrasts
#contrast    estimate    SE  df t.ratio p.value
#18mo - 36mo   -0.969 0.222 162 -4.356  <.0001 

#Degrees-of-freedom method: containment



###Saprotroph
emmeans(Dossa_psem20[[6]],specs=pairwise~"Time")

#> emmeans(Dossa_psem20[[6]],specs=pairwise~"Time")
#$emmeans
#Time emmean     SE df lower.CL upper.CL
#18mo   8.42 0.0857 26     8.24     8.59
#36mo   8.79 0.0934 26     8.60     8.98

#Degrees-of-freedom method: containment 
#Confidence level used: 0.95 

#$contrasts
#contrast    estimate    SE  df t.ratio p.value
#18mo - 36mo   -0.376 0.107 160 -3.531  0.0005 

#Degrees-of-freedom method: containment

###Check how much of variation soft and saprotroph explain in WSG loss

###Remove Time on WSG loss
Dossa_psem_soft_sapro<- psem(lme(Per_WSG_loss ~ soft_rot_ln + saprotroph_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+  soft_rot_ln + saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time +white_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem_soft_sapro)

#Individual R-squared:
#Response          method Marginal Conditional
#Per_WSG_loss      none     0.11        0.13
#Chao1_strict_sqrt none     0.28        0.35
#brown_rot_ln      none     0.11        0.22
#soft_rot_ln       none     0.15        0.28
#white_rot_ln      none     0.04        0.24
#saprotroph_ln     none     0.24        0.29

###Soft and saprotroph explain 11% Marginal and 13 % conditional.

####September

Dossa_psem22<- psem(lme(Per_WSG_loss ~ Time + Chao1_strict_sqrt + Species + Forest +Termi.assum + white_rot_ln + brown_rot_ln + soft_rot_ln, 
                       random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                   lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Brown rot brown_ln_mod2c
                   lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Softrot soft_ln_mod2d
                   lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ##White rot white_ln_mod3c
                   lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ###Saprotrophs saprotroph_ln_mod7d
                   #lme(saprotroph_ln ~ Time + Species + Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                   data=otu_Down_richness_strict_map2
)


dSep(Dossa_psem22)

###Suggesting to include some missing links
###Base on biological meaningfulness we included
### rot types on Chao1
### rot types on saprotroph
### rot type on other rot type

Dossa_psem23<- psem(lme(Per_WSG_loss ~ Time + Chao1_strict_sqrt + Species + Forest +Termi.assum + white_rot_ln + brown_rot_ln + soft_rot_ln, 
                       random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                   lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Brown rot brown_ln_mod2c
                   lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                   ##Softrot soft_ln_mod2d
                   lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum+white_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ##White rot white_ln_mod3c
                   lme(white_rot_ln ~ Time+Species+Forest+ Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                   ###Saprotrophs saprotroph_ln_mod7d
                   #lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                   data=otu_Down_richness_strict_map2
)


#+white_rot_ln + brown_rot_ln
#+ brown_rot_ln + soft_rot_ln

dSep(Dossa_psem23)
summary(Dossa_psem23)

anova(Dossa_psem22,Dossa_psem23)



###Model simplification on WSG loss
####Remove Chao1_strict_sqrt on WSG

Dossa_psem24<- psem(lme(Per_WSG_loss ~ Time + Species + Forest +Termi.assum + white_rot_ln + brown_rot_ln + soft_rot_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum+white_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time+Species+Forest+ Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    #lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)
summary(Dossa_psem24)

anova(Dossa_psem23,Dossa_psem24)



###Remove white_rot_ln on WSG loss
Dossa_psem25<- psem(lme(Per_WSG_loss ~ Time + Species + Forest +Termi.assum  + brown_rot_ln + soft_rot_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum+white_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time+Species+Forest+ Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    #lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem25)




###Remove brown_rot_ln on WSG loss

Dossa_psem26<- psem(lme(Per_WSG_loss ~ Time + Species + Forest +Termi.assum  + soft_rot_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum+white_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time+Species+Forest+ Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    #lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)
summary(Dossa_psem26)



####Remove Forest on WSG loss


Dossa_psem27<- psem(lme(Per_WSG_loss ~ Time + Species + Termi.assum  + soft_rot_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum+white_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time+Species+Forest+ Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    #lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)


summary(Dossa_psem27)

###Remove Forest on Chao

Dossa_psem28<- psem(lme(Per_WSG_loss ~ Time + Species + Termi.assum  + soft_rot_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+Termi.assum+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum+white_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time+Species+Forest+ Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    #lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)


summary(Dossa_psem28)
anova(Dossa_psem27,Dossa_psem28)

 
###Remove white_rot_ln on Chao 
Dossa_psem29<- psem(lme(Per_WSG_loss ~ Time + Species + Termi.assum  + soft_rot_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+Termi.assum+ brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum+white_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time+Species+Forest+ Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    #lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)


summary(Dossa_psem29)
anova(Dossa_psem28,Dossa_psem29)


###Remove Termi.assumm on Chao 
Dossa_psem30<- psem(lme(Per_WSG_loss ~ Time + Species + Termi.assum  + soft_rot_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+ brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum+white_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time+Species+Forest+ Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    #lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem30)
anova(Dossa_psem29,Dossa_psem30)



###Remove Termi.assum on brown_rot_ln 
Dossa_psem31<- psem(lme(Per_WSG_loss ~ Time + Species + Termi.assum  + soft_rot_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+ brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum+white_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time+Species+Forest+ Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    #lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem31)
anova(Dossa_psem30,Dossa_psem31)

###Remove white_rot_ln on soft_rot_ln


Dossa_psem32<- psem(lme(Per_WSG_loss ~ Time + Species + Termi.assum  + soft_rot_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+ brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time+Species+Forest+ Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    #lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem32)
anova(Dossa_psem31,Dossa_psem32)

###Remove Forest on white_rot_ln
Dossa_psem33<- psem(lme(Per_WSG_loss ~ Time + Species + Termi.assum  + soft_rot_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+ brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time+Species+ Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    #lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem33)
anova(Dossa_psem32,Dossa_psem33)

###Remove Species on white_rot_ln
Dossa_psem34<- psem(lme(Per_WSG_loss ~ Time + Species + Termi.assum  + soft_rot_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+ brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time+ Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    #lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem34)
anova(Dossa_psem33,Dossa_psem34)

###Remove Termi.assum on white_rot_ln
Dossa_psem35<- psem(lme(Per_WSG_loss ~ Time + Species + Termi.assum  + soft_rot_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+ brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    lme(brown_rot_ln ~ Time+Species+Forest+white_rot_ln + soft_rot_ln,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    lme(white_rot_ln ~ Time, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    #lme(saprotroph_ln ~ Time + Species + Forest+white_rot_ln + brown_rot_ln + soft_rot_ln, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem35)


###### Saprotroph  no rot type
Dossa_psem40<- psem(lme(Per_WSG_loss ~ Time + Chao1_strict_sqrt + Species + Forest +Termi.assum+saprotroph_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum+saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    #lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    #lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    #lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time + Species + Forest+Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)


dSep(Dossa_psem40)

summary(Dossa_psem40)
###Suggesting to include some missing links
####No missing links

###Simplification
###Remove Chao1_strict_sqrt on WSG

Dossa_psem41<- psem(lme(Per_WSG_loss ~ Time + Species + Forest +Termi.assum+saprotroph_ln, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum+saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    #lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    #lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    #lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time + Species + Forest+Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem41)


#Remove saprotroph_ln on WSG loss
Dossa_psem42<- psem(lme(Per_WSG_loss ~ Time + Species + Forest +Termi.assum, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum+saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    #lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    #lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    #lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time + Species + Forest+Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem42)


###Remove  Forest on WSG loss
Dossa_psem43<- psem(lme(Per_WSG_loss ~ Time + Species +Termi.assum, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Forest+Species+Termi.assum+saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    #lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    #lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    #lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time + Species + Forest+Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem43)


###Remove  Forest on Chao
Dossa_psem44<- psem(lme(Per_WSG_loss ~ Time + Species +Termi.assum, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+Termi.assum+saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    #lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    #lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    #lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time + Species + Forest+Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem44)

###Remove Termi.assum on chao
Dossa_psem45<- psem(lme(Per_WSG_loss ~ Time + Species +Termi.assum, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    #lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    #lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    #lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time + Species + Forest+Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem45)

####Remove Forest on Saprotroph
Dossa_psem46<- psem(lme(Per_WSG_loss ~ Time + Species +Termi.assum, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    #lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    #lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    #lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time + Species +Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem46)

###Remove Termi.assum on saprotroph

Dossa_psem47<- psem(lme(Per_WSG_loss ~ Time + Species +Termi.assum, 
                        random=~1|Plot/Subplot/Wood_log_ID,data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Alpha alpha_lme3b ###/Wood_log_ID because of convergence issue
                    lme(Chao1_strict_sqrt ~ Time+Species+saprotroph_ln, random=~1|Plot/Subplot/Wood_log_ID,  data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Brown rot brown_ln_mod2c
                    #lme(brown_rot_ln ~ Time+Species+Forest+Termi.assum,random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),#weights=VarExp(form=~fitted(.)),
                    ##Softrot soft_ln_mod2d
                    #lme(soft_rot_ln ~ Time + Species + Forest + Termi.assum, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ##White rot white_ln_mod3c
                    #lme(white_rot_ln ~ Time+Species+Forest, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2), #weights=VarExp(form=~fitted(.)),
                    ###Saprotrophs saprotroph_ln_mod7d
                    lme(saprotroph_ln ~ Time + Species, random=~1|Plot/Subplot/Wood_log_ID, data=otu_Down_richness_strict_map2),# weights=VarExp(form=~fitted(.)),
                    data=otu_Down_richness_strict_map2
)

summary(Dossa_psem47)


#####Beta diversity

##This part of the script records the analysis of the fungi DNA sequencing on a project in Mengsong
##that looks at the fungal diversity turn over (beta diversity) along a disturbance gradient.
##Sampling were done every 6 month interval for 36 mo but only initial, 18 and 36 mo samples were 
##analyzed.

library(vegan)
library(vegan3d)
library(ggplot2)

##Clean R's brain

rm(list = ls())

#ADONIS check pairwise for comparing categorical variable at each level

mengsong_97pick <- read.table("mengsong_97closed_guilds_r.csv", sep=",", row.names=1,header=T, check.names=F,blank.lines.skip = FALSE)   
dim(mengsong_97pick)

###import map file
map <- read.table("mengsong_map_r_chem_SBB.csv",sep=",", header=T, check.names=F,blank.lines.skip = FALSE)
rownames(map) <- map[,1]
dim(map)
#####create subset vectors
subset_sample <- rownames(map)[which(map$control=="No")]
subset_control <- rownames(map)[which(map$control=="Yes")]

sample_map <- map[subset_sample,]
control_map <- map[subset_control,]
control_map[] <- lapply(control_map, function(x) if(is.factor(x)) factor(x) else x)

sample_map <- sample_map[order(rownames(sample_map)),]
levels(sample_map$Forest) ###check the forest level
###check the data to make sure the otu table has the correct number of samples
mengsong_97pick[1:2,834:835]
#####create map for separating mengsong samples
Down_map <- sample_map[which(sample_map$SampleType=="Down" |sample_map$SampleType=="Initial"),]
Down_map[] <- lapply(Down_map, function(x) if(is.factor(x)) factor(x) else x)##make sure the factor levels are cosistent with the subset dataframe
levels(Down_map$Forest)
levels(Down_map$Position)
head(Down_map)
dim(Down_map)

dim(sample_map)
Down_strict_map<-Down_map[which(Down_map$Position=="down"),]
dim(Down_strict_map)

Down_strict_map[] <- lapply(Down_strict_map, function(x) if(is.factor(x)) factor(x) else x)##make sure the factor levels are cosistent with the subset dataframe
dim(Down_strict_map)
levels(Down_strict_map$Forest)
Down_strict_map$Forest<-as.factor(Down_strict_map$Forest)
levels(Down_strict_map$Forest)
Down_strict_map$Forest<-as.factor(Down_strict_map$Forest)
levels(Down_strict_map$Forest)
levels(Down_strict_map$Species)
Down_strict_map$Species<-as.factor(Down_strict_map$Species)
levels(Down_strict_map$Species)
levels(Down_strict_map$Position)
Down_strict_map$Position<-as.factor(Down_strict_map$Position)
levels(Down_strict_map$Position)
levels(Down_strict_map$Time)
Down_strict_map$Time<-as.factor(Down_strict_map$Time)
levels(Down_strict_map$Time)
Down_strict_map$Time<-factor(Down_strict_map$Time, ordered=TRUE)
levels(Down_strict_map$Time)
Down_strict_map$Time

Pilot_map <- sample_map[sample_map$SampleType=="PilotSample",]
Pilot_map[] <- lapply(Pilot_map, function(x) if(is.factor(x)) factor(x) else x)
####exclude pilot map
mengsong_map <- sample_map[!sample_map$SampleType %in% "PilotSample",]

UpDown_all_time_map <- mengsong_map[mengsong_map$Updown_pooling_18mo_36mo=="Yes",]
UpDown_all_time_map[] <- lapply(UpDown_all_time_map, function(x) if(is.factor(x)) factor(x) else x)
levels(UpDown_all_time_map$Position)

OpenMature_map <- mengsong_map[mengsong_map$Openmature_pooling=="Yes",]
OpenMature_map[] <- lapply(OpenMature_map, function(x) if(is.factor(x)) factor(x) else x)

####Remove initial time zero from UpDown_all-time_map, OpenMature_map
no_init_UpDown_all_time_map<-mengsong_map[mengsong_map$Updown_pooling_18mo_36mo=="Yes" & !mengsong_map$SampleType=="Initial" ,]
dim(no_init_UpDown_all_time_map)  
no_init_UpDown_all_time_map[] <- lapply(no_init_UpDown_all_time_map, function(x) if(is.factor(x)) factor(x) else x)

no_init_OpenMature_map<-mengsong_map[mengsong_map$Openmature_pooling=="Yes" & !mengsong_map$SampleType=="Initial" ,]
dim(no_init_OpenMature_map)  
no_init_OpenMature_map[] <- lapply(no_init_OpenMature_map, function(x) if(is.factor(x)) factor(x) else x)

###########filter out non fungi otus
subset_fungi <- rownames(mengsong_97pick)[which(mengsong_97pick$Kingdom=="Fungi")]
fungi_97pick <- mengsong_97pick[subset_fungi,]
dim(fungi_97pick)

fungi_control <- fungi_97pick[,subset_control] ###do not filter out 10 because samples have very few otus
fungi_97pick_sample <- fungi_97pick[,subset_sample]####only get sample otus

#####exclude less than 10 otus 
fungi_97pick_sample_11 <- subset(fungi_97pick_sample,rowSums(fungi_97pick_sample[,]) >10)
fungi_97pick_sample_11 <- fungi_97pick_sample_11[,order(colnames(fungi_97pick_sample_11))]###arrange the colnames in order


##########transform data to relative abundance
sort(colSums(fungi_97pick_sample_11), dec=T)
fungi_rare_norm <- t(t(fungi_97pick_sample_11)/colSums(fungi_97pick_sample_11))*100
fungi_dat_log <- log2(fungi_rare_norm + 1)


############################################calculate distance and do adonis on Down samples

#####Down strict adonis analysis and modeling

###Down strict data subsetting
Down_strict_log <- fungi_dat_log[,Down_strict_map$SampleID]
dim(Down_strict_log)
Down_strict_log[1:3,1:3]
Down_strict_log <- Down_strict_log[rowSums(Down_strict_log[,]) >0,]###filter out 0 otus because of subseting
dim(Down_strict_log)
Down_strict_map$Forest
Down_strict_map$Forest<-as.factor(Down_strict_map$Forest)
levels(Down_strict_map$Forest)

##Compute dissimilarity distance based on Bray Curtis method
library(vegan)
Down_strict_dis_bray <- vegdist(t(Down_strict_log),method="bray")
length(Down_strict_dis_bray)

############  Check whether through the repeated measurements time has any any effect on diversity

###Repeated measurement and random factors in adonis
#http://thebiobucket.blogspot.com/2011/04/repeat-measure-adonis-lately-i-had-to.html#more 
## Load packages
require(vegan)
### number of perms
B <- 999

####Trying to have the repeated measurment with random factor in Adonis
#### We have contacted Jari (the author of vegan). He said at the moment adonis..
#### .. does not deal with random factor. However, for testing time effect in 
####.. repeated measurments the above blog could help. Below, we adjust our data to
####.. above blog script.

##Beta diversity for Mengsong down 
##Data 
Down_strict_log
class(Down_strict_log)
dossa<-(t(Down_strict_log))
dossa[1:6,1:6]
dim(dossa)
head(Down_strict_map)
dim(Down_strict_map)

### Data:
sp <- dossa

### add time effect
### this will effect will be tested by adonis():
class(Down_strict_map)
str(Down_strict_map)
Down_strict_map[,c(3,4,5,6,7,8,9)]
a<-data.frame(Down_strict_map[,c(3,4,5,6,7,8,9)])
str(a)
levels(a$Forest)
levels(a$Forest)<-c(1,2,3)
levels(a$Time)
levels(a$Time)<-c(1,2)
a$Time<-as.numeric(a$Time)

levels(a$Species)
levels(a$Species)<-c(1,2)
a$Species<-as.numeric(a$Species)
levels(a$Position)
levels(a$Position)<-1

levels(a$Plot)
levels(a$Subplot)
sp_1 <- data.frame(Time=a$Time,Plot=a$Plot,Subplot=a$Subplot,Forest=a$Forest,sp, row.names = NULL)

### choose which species set to test:
test_sp <- sp_1

dim(test_sp)
str(test_sp)
test_sp[1:6,1:6]
###Hypothesis to be tested
#H0: Species composition is the same across time. so no time effect
#H1: Species composition differs between time points


### computing the true R2-value

### (btw, using dist() defaults to euclidean distance):
print(fit <- adonis(test_sp[,-c(4:length(ncol(test_sp)))] ~ Time, permutations=1, data=test_sp)) # "NO MORE dist() is needed"
### number of permutations
B <- 999

### setting up frame which will be populated by
### random r2 values:
pop <- rep(NA, B + 1)

### the first entry will be the true r2 value computed above:
pop[1] <- fit$aov.tab[1, 5]
pop[1]
## [1] 0.0274402

### set up a "permControl" object so it reflects our study design: Plots are like blocks and Subplots are like units within blocks
### we turn off mirroring as time should only flow in one direction
ctrl <- how(blocks = test_sp$Plot, plots= Plots(test_sp$Subplot), within = Within(type = "series", mirror = FALSE)) #"permControl" is not used at all

### Number of observations:
nobs <- nrow(test_sp)

### check permutation (...rows represent the sample id):
### ..they are ok!
### within in each repeated sample (= sites) time points are shuffled,
### with keeping the sequence intact (e.g., for site 1: 1,2,3 - 2,3,1 - 3,2,1)
shuffle(nobs, control = ctrl)

### loop:
### in adonis(...) you need to put permutations = 1, otherwise 
### adonis will not run
set.seed(2018)

###The following loop may take hours to run depending on your computer properties
###Line ## to ## We put the # sign to enable running it.
#for(i in 2:(B+1)){
#  idx <- shuffle(nobs, control = ctrl)
#  fit.rand <- adonis(test_sp[,-c(4:length(ncol(test_sp)))] ~ Time[idx],permutations = 1, data=test_sp) # "NO MORE dist() is needed"
#  pop[i] <- fit.rand$aov.tab[1, 5]
#}

### get the p-value of populations that have random r2 more than the true R2, 
### if pval less that 0.05 then there time effect:
#print(pval <- sum(pop >= pop[1]) / (B + 1))
### [1] 0.001

### the sign. p-value supports the H1 (->there is a time effect).
### ..and the fact that samples are not iid is allowed by
### the customized perms - so this p-value is trustworthy as opposed
### to tests not acknowledging dependency of data points..

## make a histogram to see random R2-values and the true one:
#hist(pop, xlab = "Population R2")
#abline(v = pop[1], col = 2, lty = 3)
#text(0.007, 60, paste("true R2,\np = ", pval, sep = ""))

##End of not run as it consumes times


###Getting random factor included in adonis
length(Down_strict_map$Plot)
Down_strict_dis_bray
length(Down_strict_dis_bray)
dim(Down_strict_map)
length(Down_strict_map$Wood_log_ID)


## Let's put ctrl to control for the way the permutations will  be done.
## We use strata to account for random factor
## strata="Down_strict_map$Wood_log_ID"

ctrl <- how(blocks = Down_strict_map$Plot, plots= Plots(Down_strict_map$Subplot), within = Within(type = "series", mirror = FALSE), nperm=999) #"permControl" is not used at all

head(Down_strict_map)
dim(Down_strict_map)

### Beta disperse http://thebiobucket.blogspot.com/2011/04/assumptions-for-permanova-with-adonis.html

dim(Down_strict_map)
length(Down_strict_dis_bray)

head(Down_strict_dis_bray)
str(Down_strict_dis_bray)

####Test for dispersion
##Beta disperse for time
betadisper(Down_strict_dis_bray, Down_strict_map$Time)
anova(betadisper(Down_strict_dis_bray, Down_strict_map$Time))

##### try out bias correction; compare with preedent model
########This bias matters most when comparing diversity among treatments with small, unequal numbers of samples. 
###########Setting bias.adjust=TRUE when using betadisper imposes a sqrt(n/(n-1)) correction (Stier et al. 2013).
betadisper(Down_strict_dis_bray, Down_strict_map$Time,bias.adjust=TRUE)
anova(betadisper(Down_strict_dis_bray, Down_strict_map$Time,bias.adjust=TRUE))

plot(betadisper(Down_strict_dis_bray, Down_strict_map$Time,bias.adjust=TRUE), conf=0.95)


##Beta disperse for Species
betadisper(Down_strict_dis_bray, Down_strict_map$Species)
anova(betadisper(Down_strict_dis_bray, Down_strict_map$Species))

##### try out bias correction; compare with precedent model
########This bias matters most when comparing diversity among treatments with small, unequal numbers of samples. 
###########Setting bias.adjust=TRUE when using betadisper imposes a sqrt(n/(n-1)) correction (Stier et al. 2013).
betadisper(Down_strict_dis_bray, Down_strict_map$Species,bias.adjust=TRUE)
anova(betadisper(Down_strict_dis_bray, Down_strict_map$Species,bias.adjust=TRUE))
plot(betadisper(Down_strict_dis_bray, Down_strict_map$Species,bias.adjust=TRUE), conf=0.95)

##Beta disperse for Forest
betadisper(Down_strict_dis_bray, Down_strict_map$Forest)
betadisp_Forest<-betadisper(Down_strict_dis_bray, Down_strict_map$Forest)
plot(betadisp_Forest)
betadisp_Forest$eig
length(betadisp_Forest$eig)
betadisp_Forest$distances
betadisp_Forest$call

str(betadisper(Down_strict_dis_bray, Down_strict_map$Forest))
anova(betadisper(Down_strict_dis_bray, Down_strict_map$Forest))

##### try out bias correction; compare with preedent model
########This bias matters most when comparing diversity among treatments with small, unequal numbers of samples. 
###########Setting bias.adjust=TRUE when using betadisper imposes a sqrt(n/(n-1)) correction (Stier et al. 2013).
betadisper(Down_strict_dis_bray, Down_strict_map$Forest,bias.adjust=TRUE)
anova(betadisper(Down_strict_dis_bray, Down_strict_map$Forest,bias.adjust=TRUE))
plot(betadisper(Down_strict_dis_bray, Down_strict_map$Forest,bias.adjust=TRUE), conf=0.95)

##Post hoc Tukey HSD
TukeyHSD(betadisper(Down_strict_dis_bray, Down_strict_map$Forest))
(mod.Forest.HSD <- TukeyHSD(betadisper(Down_strict_dis_bray, Down_strict_map$Forest,bias.adjust=TRUE)))
plot(mod.Forest.HSD)


###Re  draw the plot
## Results
#> (mod.Forest.HSD <- TukeyHSD(betadisper(Down_strict_dis_bray, Down_strict_map$Forest,bias.adjust=TRUE)))
#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = distances ~ group, data = df)

#$group
#                                       diff         lwr         upr    p adj
#Open_land-Mature_forest           -0.01132795 -0.02010711 -0.00254880 0.007214
#Regenerating_forest-Mature_forest  0.00177016 -0.00495896  0.00849929 0.809831
#Regenerating_forest-Open_land      0.01309811  0.00460096  0.02159526 0.000949


str(mod.Forest.HSD)
betadis_tuke<-as.data.frame(mod.Forest.HSD$group)
colnames(betadis_tuke)
betadis_tuke$Pairs
betadis_tuke$Pairs<-c("Open_Mature", "Renegerating_Mature", "Regenerating_Open")

betadis_tuke$Pairs1<-c("Open vs. Mature forest", "Renegerating vs. Mature forest", "Regenerating vs. Open forest")
str(betadis_tuke)

library(ggplot2)

betadisper_tukey<-ggplot(data=betadis_tuke, aes(x=diff, y=Pairs))+
  geom_point()+  
  geom_errorbarh(aes(xmax= betadis_tuke$upr, xmin=betadis_tuke$lwr, height =0))+
  geom_vline(xintercept = 0, lty=2)+
  xlab("Differences in centroids means") +
  ylab("Pairs comparison")+ 
  theme(plot.title=element_text(size=rel(1.5),colour="black"))+
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text.y=element_text(colour = "black", size = 8),
        axis.ticks.y=element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        strip.text = element_text(colour = "black", size = 12),
        legend.position = "none")
betadisper_tukey

betadisper_tukey1<-ggplot(data=betadis_tuke, aes(x=diff, y=Pairs1))+
  geom_point()+  
  geom_errorbarh(aes(xmax= betadis_tuke$upr, xmin=betadis_tuke$lwr, height =0))+
  geom_vline(xintercept = 0, lty=2)+
  xlab("Differences in centroids means") +
  ylab("Pairs comparison")+ 
  theme(plot.title=element_text(size=rel(1.5),colour="black"))+
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text=element_text(colour = "black", size = 10),
        axis.ticks=element_line(colour = "black"),
        axis.title = element_text(colour = "black", size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        strip.text = element_text(colour = "black", size = 12),
        legend.position = "none")
betadisper_tukey1

tiff(filename="Figure S6 Beta diversity beta disperser forest.tiff", res = 600, width=4000, height=2000, compression = "lzw")
betadisper_tukey1
dev.off()


###PerMANOVA/adonis

## Termites assum
###Full model

summary(Down_strict_map)

Down_strict_map$Termi.assum
Down_strict_map$Termi.assum<-as.factor(Down_strict_map$Termi.assum)
Down_strict_map$Termi.assum

Down_strict_dis_bray_ado1b <- adonis(Down_strict_dis_bray ~ Time+Forest+Species+Termi.assum+Time:Forest+Time:Species+Forest:Species+Time:Termi.assum+Forest:Termi.assum+Species:Termi.assum+Time:Forest:Species+
                                   Time:Forest:Termi.assum+Time:Species:Termi.assum+Forest:Species:Termi.assum+Time:Forest:Species:Termi.assum,permutations=ctrl,strata= "Down_strict_map$Wood_log_ID", data=Down_strict_map)
Down_strict_dis_bray_ado1b

####Model simplification
##Remove non-significant variables
###The four interaction Time:Forest:Species:Termites

Down_strict_dis_bray_ado2b <- adonis(Down_strict_dis_bray ~ Time+Forest+Species+Termi.assum+Time:Forest+Time:Species+Forest:Species+Time:Termi.assum+Forest:Termi.assum+Species:Termi.assum+Time:Forest:Species+
                                       Time:Forest:Termi.assum+Time:Species:Termi.assum+Forest:Species:Termi.assum,permutations=ctrl,strata= "Down_strict_map$Wood_log_ID", data=Down_strict_map)
Down_strict_dis_bray_ado2b


####Remove -Forest:Species:Termites
Down_strict_dis_bray_ado3b <- adonis(Down_strict_dis_bray ~ Time+Forest+Species+Termi.assum+Time:Forest+Time:Species+Forest:Species+Time:Termi.assum+Forest:Termi.assum+Species:Termi.assum+Time:Forest:Species+
                                       Time:Forest:Termi.assum+Time:Species:Termi.assum,permutations=ctrl,strata= "Down_strict_map$Wood_log_ID", data=Down_strict_map)
Down_strict_dis_bray_ado3b

####Remove -Time:Species:Termites
Down_strict_dis_bray_ado4b <- adonis(Down_strict_dis_bray ~ Time+Forest+Species+Termi.assum+Time:Forest+Time:Species+Forest:Species+Time:Termi.assum+Forest:Termi.assum+Species:Termi.assum+Time:Forest:Species+
                                       Time:Forest:Termi.assum,permutations=ctrl,strata= "Down_strict_map$Wood_log_ID", data=Down_strict_map)
Down_strict_dis_bray_ado4b

####Remove -Time:Forest:Termites
Down_strict_dis_bray_ado5b <- adonis(Down_strict_dis_bray ~ Time+Forest+Species+Termi.assum+Time:Forest+Time:Species+Forest:Species+Time:Termi.assum+Forest:Termi.assum+Species:Termi.assum+Time:Forest:Species
                                     ,permutations=ctrl,strata= "Down_strict_map$Wood_log_ID", data=Down_strict_map)
Down_strict_dis_bray_ado5b

####Remove -Time:Forest:Species
Down_strict_dis_bray_ado6b <- adonis(Down_strict_dis_bray ~ Time+Forest+Species+Termi.assum+Time:Forest+Time:Species+Forest:Species+Time:Termi.assum+Forest:Termi.assum+Species:Termi.assum
                                     ,permutations=ctrl,strata= "Down_strict_map$Wood_log_ID", data=Down_strict_map)
Down_strict_dis_bray_ado6b



#####Stop simplification here
###Best model Down_strict_dis_bray_ado6b
Down_strict_dis_bray_ado6b

##### Modeling of beta diversity with adonis analysis 

##Table 3 
Down_strict_dis_bray_ado6b$aov.tab


#Blocks:  Down_strict_map$Plot 
#Plots: Down_strict_map$Subplot, plot permutation: none
#Permutation: series
#Number of permutations: 999

#Terms added sequentially (first to last)

#Terms added sequentially (first to last)

#                     Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
#Time                  1      4.69   4.695  11.342 0.0274  0.001 ***
#Forest                2      2.21   1.105   2.669 0.0129  0.001 ***
#Species               1      4.11   4.113   9.938 0.0240  0.001 ***
#Termi.assum           1      1.03   1.033   2.496 0.0060  0.422    
#Time:Forest           2      1.11   0.555   1.341 0.0065  0.001 ***
#Time:Species          1      1.74   1.739   4.202 0.0102  0.001 ***
#Forest:Species        2      1.13   0.565   1.366 0.0066  0.009 ** 
#Time:Termi.assum      1      0.52   0.524   1.266 0.0031  0.017 *  
#Forest:Termi.assum    2      1.10   0.550   1.328 0.0064  0.011 *  
#Species:Termi.assum   1      0.69   0.693   1.674 0.0040  0.001 ***
#Residuals           369    152.74   0.414         0.8928           
#Total               383    171.08                 1.0000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



######Interpretation
######Fungal species community turn over is affected by time, species, forest type and interaction between Time and Termites, Time and forest,Times:Species
#########################################################################Forest: Termites, SPecies:Termites.
######Together this model factor explains (1-0.8928)*100=10.72% 
###### The beta disperse which is also known as beta diversity significance suggesting that the adonis result could be due to either



###Conduct 

##db-RDA Distance based RDA

### All factors
capscale(Down_strict_dis_bray~Time+Forest+Species+Termi.assum+Time:Forest+Time:Species+Forest:Species+Time:Termi.assum+Forest:Termi.assum+Species:Termi.assum, data=Down_strict_map)
Dossa_db_rda<-dbrda(Down_strict_dis_bray~Time+Forest+Termi.assum+Time:Termi.assum, data=Down_strict_map)
Dossa_db_rda$CA
Dossa_db_rda$Ybar
Dossa_db_rda$Ybar[1:6,1:6]
Dossa_db_rda$CCA

dim(Dossa_db_rda$CCA$wa)

dim(Dossa_db_rda$CCA$u)

head(Dossa_db_rda$CCA$u)
a<-summary(dbrda(Down_strict_dis_bray~Time+Forest+Termi.assum+Time:Termi.assum, data=Down_strict_map))
dim(a$constraints)
dim(a$sites)

#write.csv(a$constraints, file="db_rda_constraints.csv")
Down_strict_map[1:6,1:6]
Down_strict_map[380:384,1:6]


anova(dbrda(Down_strict_dis_bray~Time+Forest+Species+Termi.assum+Time:Forest+Time:Species+Forest:Species+Time:Termi.assum+Forest:Termi.assum+Species:Termi.assum, data=Down_strict_map))

###Make sure understanding how to run for eac parameter
### Only time
capscale(Down_strict_dis_bray~Time, data=Down_strict_map)
dbrda(Down_strict_dis_bray~Time, data=Down_strict_map)
summary(dbrda(Down_strict_dis_bray~Time, data=Down_strict_map))
anova(dbrda(Down_strict_dis_bray~Time, data=Down_strict_map))

#> anova(dbrda(Down_strict_dis_bray~Time, data=Down_strict_map))
#Permutation test for dbrda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: dbrda(formula = Down_strict_dis_bray ~ Time, data = Down_strict_map)
#           Df SumOfSqs     F Pr(>F)    
#Model      1     4.69 10.78  0.001 ***
#Residual 382   166.39                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

### Only Forest
capscale(Down_strict_dis_bray~Forest, data=Down_strict_map)
dbrda(Down_strict_dis_bray~Forest, data=Down_strict_map)
summary(dbrda(Down_strict_dis_bray~Forest, data=Down_strict_map))
anova(dbrda(Down_strict_dis_bray~Forest, data=Down_strict_map))

###Only time and termites

capscale(Down_strict_dis_bray~Time*Termi.assum, data=Down_strict_map)
dbrda(Down_strict_dis_bray~Time*Termi.assum, data=Down_strict_map)
summary(dbrda(Down_strict_dis_bray~Time*Termi.assum, data=Down_strict_map))
anova(dbrda(Down_strict_dis_bray~Time*Termi.assum, data=Down_strict_map))



#######################do procrustes on the down samples
####Procrustes test the null hypothesis is that there is no relationship between... 
####...the position of the points in the first ordination and the second.
####create matrix by time
####Do it for only strict!!!!!!
###MAP of DOWN STRICT 

##Clean R's brain
rm(list = ls())

mengsong_97pick <- read.table("mengsong_97closed_guilds_r.csv", sep=",", row.names=1,header=T, check.names=F,blank.lines.skip = FALSE)   
dim(mengsong_97pick)

###import map file
map <- read.table("mengsong_map_r_chem_SBB.csv",sep=",", header=T, check.names=F,blank.lines.skip = FALSE)
rownames(map) <- map[,1]
dim(map)
#####create subset vectors
subset_sample <- rownames(map)[which(map$control=="No")]
subset_control <- rownames(map)[which(map$control=="Yes")]

sample_map <- map[subset_sample,]
control_map <- map[subset_control,]
control_map[] <- lapply(control_map, function(x) if(is.factor(x)) factor(x) else x)

sample_map <- sample_map[order(rownames(sample_map)),]
levels(sample_map$Forest) ###check the forest level

###check the data to make sure the otu table has the correct number of samples
mengsong_97pick[1:2,834:835]
#####create map for separating mengsong samples
Down_map <- sample_map[which(sample_map$SampleType=="Down" |sample_map$SampleType=="Initial"),]
#######make sure the factor levels are cosistent with the subset dataframe
Down_map[] <- lapply(Down_map, function(x) if(is.factor(x)) factor(x) else x)
levels(Down_map$Forest)
levels(Down_map$Position)
head(Down_map)
dim(Down_map)

dim(sample_map)
Down_strict_map<-Down_map[which(Down_map$Position=="down"),]
dim(Down_strict_map)
str(Down_strict_map)
Down_strict_map[] <- lapply(Down_strict_map, function(x) if(is.factor(x)) factor(x) else x)##make sure the factor levels are cosistent with the subset dataframe
str(Down_strict_map)
dim(Down_strict_map)
dim(Down_strict_map)
levels(Down_strict_map$Forest)
Down_strict_map$Forest<-as.factor(Down_strict_map$Forest)
levels(Down_strict_map$Forest)
Down_strict_map$Forest<-as.factor(Down_strict_map$Forest)
levels(Down_strict_map$Forest)
levels(Down_strict_map$Species)
Down_strict_map$Species<-as.factor(Down_strict_map$Species)
levels(Down_strict_map$Species)
levels(Down_strict_map$Position)
Down_strict_map$Position<-as.factor(Down_strict_map$Position)
levels(Down_strict_map$Position)
levels(Down_strict_map$Time)
Down_strict_map$Time<-as.factor(Down_strict_map$Time)
levels(Down_strict_map$Time)
Down_strict_map$Time<-factor(Down_strict_map$Time, ordered=TRUE)
levels(Down_strict_map$Time)
Down_strict_map$Time

######Creating a map for initial samples
Initial_samp<-sample_map[which(sample_map$SampleType=="Initial"),]
Initial_samp[] <- lapply(Initial_samp, function(x) if(is.factor(x)) factor(x) else x)
levels(Initial_samp$Forest)
levels(Initial_samp$Position)
head(Initial_samp)
dim(Initial_samp)


Pilot_map <- sample_map[sample_map$SampleType=="PilotSample",]
Pilot_map[] <- lapply(Pilot_map, function(x) if(is.factor(x)) factor(x) else x)

####exclude pilot map
mengsong_map <- sample_map[!sample_map$SampleType %in% "PilotSample",]

UpDown_all_time_map <- mengsong_map[mengsong_map$Updown_pooling_18mo_36mo=="Yes",]
UpDown_all_time_map[] <- lapply(UpDown_all_time_map, function(x) if(is.factor(x)) factor(x) else x)
levels(UpDown_all_time_map$Position)

OpenMature_map <- mengsong_map[mengsong_map$Openmature_pooling=="Yes",]
OpenMature_map[] <- lapply(OpenMature_map, function(x) if(is.factor(x)) factor(x) else x)

####Remove initial time zero from UpDown_all-time_map, OpenMature_map
no_init_UpDown_all_time_map<-mengsong_map[mengsong_map$Updown_pooling_18mo_36mo=="Yes" & !mengsong_map$SampleType=="Initial" ,]
dim(no_init_UpDown_all_time_map)  
no_init_UpDown_all_time_map[] <- lapply(no_init_UpDown_all_time_map, function(x) if(is.factor(x)) factor(x) else x)

no_init_OpenMature_map<-mengsong_map[mengsong_map$Openmature_pooling=="Yes" & !mengsong_map$SampleType=="Initial" ,]
dim(no_init_OpenMature_map)  
no_init_OpenMature_map[] <- lapply(no_init_OpenMature_map, function(x) if(is.factor(x)) factor(x) else x)

####



###########filter out non fungi otus
subset_fungi <- rownames(mengsong_97pick)[which(mengsong_97pick$Kingdom=="Fungi")]
fungi_97pick <- mengsong_97pick[subset_fungi,]
dim(fungi_97pick)

fungi_control <- fungi_97pick[,subset_control] ###do not filter out 10 because samples have very few otus
fungi_97pick_sample <- fungi_97pick[,subset_sample]####only get sample otus

#####exclude less than 10 otus 
fungi_97pick_sample_11 <- subset(fungi_97pick_sample,rowSums(fungi_97pick_sample[,]) >10)
fungi_97pick_sample_11 <- fungi_97pick_sample_11[,order(colnames(fungi_97pick_sample_11))]###arrange the colnames in order


##########transform data to relative abundance
sort(colSums(fungi_97pick_sample_11), dec=T)
fungi_rare_norm <- t(t(fungi_97pick_sample_11)/colSums(fungi_97pick_sample_11))*100
fungi_dat_log <- log2(fungi_rare_norm + 1)

str(fungi_dat_log)
dim(fungi_dat_log)
head(fungi_dat_log,2)

###Down strict data subsetting
fungi_dat_log<-as.data.frame(fungi_dat_log)

str(fungi_dat_log)
Down_strict_log <- fungi_dat_log[,as.character(Down_strict_map$SampleID)]
dim(Down_strict_log)
str(fungi_dat_log)
str(Down_strict_log)

Down_strict_log <- Down_strict_log[rowSums(Down_strict_log[,]) >0,]###filter out 0 otus because of subseting
dim(Down_strict_log)
levels(Down_strict_map$Forest)

require(vegan)

Down_strict_dis_bray <- vegdist(t(Down_strict_log),method="bray")

###########the number of samples are not the same, need to subset and then do on 18mo and 36mo
###Need to create a 18 mo and 36 mo having same number ofsamples and identical ID of samples
###Since 36 mo has less sampled log than 18 mo, we will match the log ID of 36 mo to the same of
###18 mo and get the same number of logs. We matched them and found common 163 samples

##Do transpose to have samples x otus matrix

dossa<-(t(Down_strict_log))
dim(dossa)
colnames(dossa)[1:6]
rownames(dossa)[1:6]
a<-data.frame(Down_strict_map[,c(3,4,5,6,7,8,9)])
str(a)
dim(a)
head(a,6,6)
head(dossa,6,6)
str(dossa)
sp_1 <- cbind(a,dossa)
str(sp_1)
sp_1<-as.data.frame(sp_1)
dim(sp_1)
head(sp_1)
str(sp_1)
###

Down_strict_log_18mo <- sp_1[sp_1$Time=="18mo",]
Down_strict_log_18mo$Time
str(Down_strict_log_18mo)
Down_strict_log_18mo$Wood_log_ID
Down_strict_log_36mo <- sp_1[sp_1$Time=="36mo",]
Down_strict_log_36mo$Time
Down_strict_log_36mo$Wood_log_ID

str(Down_strict_log_36mo)
dim(Down_strict_log_18mo)
dim(Down_strict_log_36mo)

c<-merge(Down_strict_log_36mo,Down_strict_log_18mo, by.x="Wood_log_ID", by.y = "Wood_log_ID")
dim(c)
str(c)
names(c)
d<-as.vector(c$Wood_log_ID)
length(d)
rum<-colnames(c)
rum[9449:9455]
length(colnames(c))

##Now retrieve same log tag  information for 18 mo and 36 mo
down_common_36mo<-c[,1:9080]
str(down_common_36mo)

down_common_18mo<-c[,c(1:7,9081:length(colnames(c)))]
str(down_common_18mo)
dim(down_common_18mo)
aDown_strict_log_18mo <- c[c$Wood_log_ID==d,]
aDown_strict_log_18mo <- c[,c$Wood_log_ID==d]
str(aDown_strict_log_18mo)
dim(aDown_strict_log_18mo)

bDown_strict_log_18mo <- aDown_strict_log_18mo[aDown_strict_log_18mo$Wood_log_ID==c$Wood_log_ID,]
dim(bDown_strict_log_18mo)
cDown_strict_log_36mo <- c[,c$Wood_log_ID==d]
dDown_strict_log_36mo<-cDown_strict_log_36mo[cDown_strict_log_36mo$Wood_log_ID==c$Wood_log_ID,]


dim(aDown_strict_log_18mo)
dim(bDown_strict_log_18mo)
dim(cDown_strict_log_36mo)
dim(dDown_strict_log_36mo)
str(aDown_strict_log_18mo)

str(bDown_strict_log_18mo)
str(dDown_strict_log_36mo)

##Get this dataframes back to their original length with original columns names
bDown_strict_log_18mo<-bDown_strict_log_18mo[,c(1,9081:dim(bDown_strict_log_18mo)[2])]
dim(bDown_strict_log_18mo)
str(bDown_strict_log_18mo)
colnames(bDown_strict_log_18mo)<-c("Wood_log_ID", "Forest", "Plot", "Subplot", "Time", "Species", "Position", colnames(dossa))
str(bDown_strict_log_18mo)

dDown_strict_log_36mo<-dDown_strict_log_36mo[,1:9080]
dim(dDown_strict_log_36mo)
str(dDown_strict_log_36mo)
colnames(dDown_strict_log_36mo)<-c("Wood_log_ID", "Forest", "Plot", "Subplot", "Time", "Species", "Position", colnames(dossa))
str(dDown_strict_log_36mo)

bDown_strict_log_18mo$Wood_log_ID
dDown_strict_log_36mo$Wood_log_ID

head(bDown_strict_log_18mo)
head(dDown_strict_log_36mo)

###Sort Wood_log_ID so we know for sure procrustes is working normally
sort_log_18mo<-bDown_strict_log_18mo[order(bDown_strict_log_18mo$Wood_log_ID),]
sort_log_36mo<-dDown_strict_log_36mo[order(dDown_strict_log_36mo$Wood_log_ID),]

is.matrix(sort_log_18mo)
is.matrix(sort_log_36mo)

##Create procrustes matrices
sort_log_18mo<-sort_log_18mo[,c(1,8:length(colnames(sort_log_18mo)))]
sort_log_36mo<-sort_log_36mo[,c(1,8:length(colnames(sort_log_36mo)))]
dim(sort_log_18mo)
dim(sort_log_36mo)
sort_log_18mo[1:3,1:3]
row.names(sort_log_18mo)<-sort_log_18mo[,1]
sort_log_18mo<-sort_log_18mo[,-1]
sort_log_18mo[1:3,1:3]
proc_1<-as.matrix(sort_log_18mo)

is.matrix(proc_1)
str(proc_1)

sort_log_36mo[1:3,1:3]
row.names(sort_log_36mo)<-sort_log_36mo[,1]
sort_log_36mo<-sort_log_36mo[,-1]
sort_log_36mo[1:3,1:3]
proc_2<-as.matrix(sort_log_36mo)

is.matrix(proc_2)
str(proc_2)


###Disimilarity distance bray curtis
log_18mo.mds2 = metaMDS(proc_1,distance='bray',k=3)
log_36mo.mds2 = metaMDS(proc_2,distance='bray',k=3)


# Procrustes analysis
set.seed(2018)
log.pro2 = procrustes(log_18mo.mds2$points,log_36mo.mds2$points)

protest(log_18mo.mds2$points,log_36mo.mds2$points)

#Call:
#  protest(X = log_18mo.mds2$points, Y = log_36mo.mds2$points) 

#Procrustes Sum of Squares (m12 squared):         0.8644  
#Correlation in a symmetric Procrustes rotation:  0.3683 
#Significance:  0.001 

#Permutation: free
#Number of permutations: 999


##NULL Ho: The degree of concordance between two (or more) matrices is no greater than 
###expected given random inter-matrix associations.

##So we reject Ho and acept H1  
##H1:The degree of concordance between two (or more) matrices is greater than 
###expected given random inter-matrix associations.

## Significance:  0.001 We accept the H1 hypothesis that there is a relationship among points...
## between...the position of the points in the first ordination (here 18mo) and the second (36mo).
## Priority effect

par(mar=c(5,5,2,2))
par(mfrow=c(1,1))
par(pty='s')
str(log.pro2)
plot(log.pro2,choices=c(1,2),main='',cex.axis=2,cex.lab=2,cex=2,lwd=6)
plot(log.pro2,choices=c(1,3),main='',cex.axis=2,cex.lab=2,cex=2,lwd=6)


dim(Down_strict_log)


###Procrustes between time 0 mo and 18 mo
####Trying to have a map with 0, 18 and 36 mo together
##Clean R's brain
#rm(list = ls())

mengsong_97pick <- read.table("mengsong_97closed_guilds_r.csv", sep=",", row.names=1,header=T, check.names=F,blank.lines.skip = FALSE)   
dim(mengsong_97pick)

###import map file
map <- read.table("mengsong_map_r_chem_SBB.csv",sep=",", header=T, check.names=F,blank.lines.skip = FALSE)
rownames(map) <- map[,1]
dim(map)
#####create subset vectors
subset_sample <- rownames(map)[which(map$control=="No")]
subset_control <- rownames(map)[which(map$control=="Yes")]

sample_map <- map[subset_sample,]
control_map <- map[subset_control,]
control_map[] <- lapply(control_map, function(x) if(is.factor(x)) factor(x) else x)

sample_map <- sample_map[order(rownames(sample_map)),]
levels(sample_map$Forest) ###check the forest level

###check the data to make sure the otu table has the correct number of samples
mengsong_97pick[1:2,834:835]
#####create map for separating mengsong samples
Down_map <- sample_map[which(sample_map$SampleType=="Down" |sample_map$SampleType=="Initial"),]
#######make sure the factor levels are cosistent with the subset dataframe
Down_map[] <- lapply(Down_map, function(x) if(is.factor(x)) factor(x) else x)
levels(Down_map$Forest)
Down_map$Forest<-as.factor(Down_map$Forest)
levels(Down_map$Forest)
levels(Down_map$Position)
Down_map$Position<-as.factor(Down_map$Position)
levels(Down_map$Position)
head(Down_map)
dim(Down_map)

dim(sample_map)
Down_strict_initial_map<-Down_map[which(Down_map$Position=="down"|Down_map$Position=="Initial"),]
dim(Down_strict_initial_map)
str(Down_strict_initial_map)
Down_strict_initial_map[] <- lapply(Down_strict_initial_map, function(x) if(is.factor(x)) factor(x) else x)##make sure the factor levels are cosistent with the subset dataframe
str(Down_strict_initial_map)
dim(Down_strict_initial_map)
levels(Down_strict_initial_map$Forest)
Down_strict_initial_map$Forest<-as.factor(Down_strict_initial_map$Forest)
levels(Down_strict_initial_map$Forest)
levels(Down_strict_initial_map$Species)
Down_strict_initial_map$Species<-as.factor(Down_strict_initial_map$Species)
levels(Down_strict_initial_map$Species)

levels(Down_strict_initial_map$Position)
Down_strict_initial_map$Position<-as.factor(Down_strict_initial_map$Position)
levels(Down_strict_initial_map$Position)

######Creating a map for initial samples
Initial_samp<-sample_map[which(sample_map$SampleType=="Initial"),]
Initial_samp[] <- lapply(Initial_samp, function(x) if(is.factor(x)) factor(x) else x)
levels(Initial_samp$Forest)
levels(Initial_samp$Position)
head(Initial_samp)
dim(Initial_samp)


Pilot_map <- sample_map[sample_map$SampleType=="PilotSample",]
Pilot_map[] <- lapply(Pilot_map, function(x) if(is.factor(x)) factor(x) else x)

####exclude pilot map
mengsong_map <- sample_map[!sample_map$SampleType %in% "PilotSample",]

UpDown_all_time_map <- mengsong_map[mengsong_map$Updown_pooling_18mo_36mo=="Yes",]
UpDown_all_time_map[] <- lapply(UpDown_all_time_map, function(x) if(is.factor(x)) factor(x) else x)
levels(UpDown_all_time_map$Position)

OpenMature_map <- mengsong_map[mengsong_map$Openmature_pooling=="Yes",]
OpenMature_map[] <- lapply(OpenMature_map, function(x) if(is.factor(x)) factor(x) else x)

####Remove initial time zero from UpDown_all-time_map, OpenMature_map
no_init_UpDown_all_time_map<-mengsong_map[mengsong_map$Updown_pooling_18mo_36mo=="Yes" & !mengsong_map$SampleType=="Initial" ,]
dim(no_init_UpDown_all_time_map)  
no_init_UpDown_all_time_map[] <- lapply(no_init_UpDown_all_time_map, function(x) if(is.factor(x)) factor(x) else x)

no_init_OpenMature_map<-mengsong_map[mengsong_map$Openmature_pooling=="Yes" & !mengsong_map$SampleType=="Initial" ,]
dim(no_init_OpenMature_map)  
no_init_OpenMature_map[] <- lapply(no_init_OpenMature_map, function(x) if(is.factor(x)) factor(x) else x)

###########filter out non fungi otus
subset_fungi <- rownames(mengsong_97pick)[which(mengsong_97pick$Kingdom=="Fungi")]
fungi_97pick <- mengsong_97pick[subset_fungi,]
dim(fungi_97pick)

fungi_control <- fungi_97pick[,subset_control] ###do not filter out 10 because samples have very few otus
fungi_97pick_sample <- fungi_97pick[,subset_sample]####only get sample otus

#####exclude less than 10 otus 
fungi_97pick_sample_11 <- subset(fungi_97pick_sample,rowSums(fungi_97pick_sample[,]) >10)
fungi_97pick_sample_11 <- fungi_97pick_sample_11[,order(colnames(fungi_97pick_sample_11))]###arrange the colnames in order


##########transform data to relative abundance
sort(colSums(fungi_97pick_sample_11), dec=T)
fungi_rare_norm <- t(t(fungi_97pick_sample_11)/colSums(fungi_97pick_sample_11))*100
###Log 2 transformation
fungi_dat_log <- log2(fungi_rare_norm + 1)

str(fungi_dat_log)
dim(fungi_dat_log)
head(fungi_dat_log,2)

###Down strict and initial data subsetting
fungi_dat_log<-as.data.frame(fungi_dat_log)

str(fungi_dat_log)
Down_strict_initial_log <- fungi_dat_log[,as.character(Down_strict_initial_map$SampleID)]
dim(Down_strict_initial_log)
str(Down_strict_initial_log)
str(Down_strict_initial_log)

###filter out 0 otus because of subseting
Down_strict_initial_log <- Down_strict_initial_log[rowSums(Down_strict_initial_log[,]) >0,]
dim(Down_strict_initial_log)
levels(Down_strict_initial_map$Forest)

require(vegan)

Down_strict_initial_dis_bray <- vegdist(t(Down_strict_initial_log),method="bray")

###########the number of samples are not the same, we need to subset 0 mo and 18 mo
###Need to create a 0 mo, 18 mo having same number of samples and identical ID of samples

dossa2<-(t(Down_strict_initial_log))
dim(dossa2)
colnames(dossa2)[1:6]
rownames(dossa2)[1:6]
head(Down_strict_initial_map)
tail(Down_strict_initial_map)

a<-data.frame(Down_strict_initial_map[,c(2,3,4,5,6,7,8,9)])
str(a)
dim(a)
head(a,6,6)
head(dossa2,6,6)
str(dossa2)

###Combine a and dossa 2
sp_2 <- cbind(a,dossa2)
str(sp_2)
sp_2<-as.data.frame(sp_2)
dim(sp_2)
head(sp_2)
str(sp_2)

###
Down_strict_initial_log_18mo <- sp_2[sp_2$Time=="18mo",]
Down_strict_initial_log_18mo$Time
Down_strict_initial_log_18mo[1:6,1:9]
str(Down_strict_initial_log_18mo)
Down_strict_initial_log_18mo$Wood_log_ID
dim(Down_strict_initial_log_18mo)


#### Need to create time 0mo 
Down_strict_initial_log_0mo <- sp_2[sp_2$Time=="0mo",]
Down_strict_initial_log_0mo$Time
str(Down_strict_initial_log_0mo)
Down_strict_initial_log_0mo$Wood_log_ID


####This is based on tree id, we need to repeat time 0 community information of all logs from same tree
####We need to create a dataframe of log_ID and tree of 18 mo length, then merge this with the data on time 0mo based on tree

Tree_log_id<-Down_strict_initial_log_18mo[,c(1,5)]

dim(Tree_log_id)

#colnames(Tree_log_id)<-c("Tree","Log_ID")
head((Tree_log_id))

###Create a merged data on time 0mo that contains log id and repeat info ofcommunity
merge_0mo<-merge(Tree_log_id,Down_strict_initial_log_0mo,  by.x="Tree", by.y="Tree", all.x=T)

dim(Down_strict_initial_log_0mo)
dim(merge_0mo)
str(merge_0mo)
head(merge_0mo)
merge_0mo[1:6,1:6]
#Removed the second identical Wood_log_ID.y 
merge_0mo<-merge_0mo[,-6]

head((Tree_log_id))
#rownames(order(merge_0mo$Wood_log_ID.x))

#merge_0mo$Wood_log_ID<-merge_0mo$Wood_log_ID.x
Down_strict_initial_log_0mo<-merge_0mo
dim(Down_strict_initial_log_0mo)
head(Down_strict_initial_log_0mo)
Down_strict_initial_log_0mo[1:6,1:9]



dim(Down_strict_initial_log_18mo)

dim(Down_strict_initial_log_0mo)

###Sort Wood_log_ID so we know for sure procrustes is working normally
sort_Down_strict_initial_log_0mo<-Down_strict_initial_log_0mo[order(Down_strict_initial_log_0mo$Wood_log_ID.x),]
sort_Down_strict_initial_log_18mo<-Down_strict_initial_log_18mo[order(Down_strict_initial_log_18mo$Wood_log_ID),]

is.matrix(sort_Down_strict_initial_log_0mo)
is.matrix(sort_Down_strict_initial_log_18mo)

Log_not_needed<-sort_Down_strict_initial_log_0mo[which(is.na(sort_Down_strict_initial_log_0mo$"13377")),]
dim(Log_not_needed)                          

##Remove those missing log from 18mo
sort_Down_strict_initial_log_18mo<-sort_Down_strict_initial_log_18mo[-which(is.na(sort_Down_strict_initial_log_0mo$"13377")),]
dim(sort_Down_strict_initial_log_18mo)

###Now remove from 0mo
sort_Down_strict_initial_log_0mo<-sort_Down_strict_initial_log_0mo[-which(is.na(sort_Down_strict_initial_log_0mo$"13377")),]
dim(sort_Down_strict_initial_log_0mo)

##Create procrustes matrices
sort_Down_strict_initial_log_0mo[1:6,1:10]
sort_Down_strict_initial_log_0mo<-sort_Down_strict_initial_log_0mo[,c(2,9:length(colnames(sort_Down_strict_initial_log_0mo)))]

sort_Down_strict_initial_log_18mo[1:6,1:10]
sort_Down_strict_initial_log_18mo<-sort_Down_strict_initial_log_18mo[,c(5,9:length(sort_Down_strict_initial_log_18mo))]

dim(sort_Down_strict_initial_log_0mo)
dim(sort_Down_strict_initial_log_18mo)
sort_Down_strict_initial_log_0mo[1:3,1:3]
sort_Down_strict_initial_log_18mo[1:3,1:3]

row.names(sort_Down_strict_initial_log_0mo)<-sort_Down_strict_initial_log_0mo[,1]
sort_Down_strict_initial_log_0mo[1:3,1:3]
sort_Down_strict_initial_log_0mo<-sort_Down_strict_initial_log_0mo[,-1]
sort_Down_strict_initial_log_0mo[1:3,1:3]
proc_time0mo<-as.matrix(sort_Down_strict_initial_log_0mo)

is.matrix(proc_time0mo)
str(proc_time0mo)

sort_Down_strict_initial_log_18mo[1:3,1:3]
row.names(sort_Down_strict_initial_log_18mo)<-sort_Down_strict_initial_log_18mo[,1]
sort_Down_strict_initial_log_18mo[1:3,1:3]
sort_Down_strict_initial_log_18mo<-sort_Down_strict_initial_log_18mo[,-1]
sort_Down_strict_initial_log_18mo[1:3,1:3]
proc_time18mo<-as.matrix(sort_Down_strict_initial_log_18mo)

is.matrix(proc_time18mo)
str(proc_time18mo)


###Disimilarity distance bray curtis
set.seed(2018)
log_0mo.mds2 = metaMDS(proc_time0mo,distance='bray',k=3)
log_18mo.mds2 = metaMDS(proc_time18mo,distance='bray',k=3)


# Procrustes analysis
initial_log.pro2 = procrustes(log_0mo.mds2$points,log_18mo.mds2$points)
protest(log_0mo.mds2$points,log_18mo.mds2$points)

#Call:
#  protest(X = log_0mo.mds2$points, Y = log_18mo.mds2$points) 

#Procrustes Sum of Squares (m12 squared):        0.8818
#Correlation in a symmetric Procrustes rotation: 0.3438 
#Significance:  0.001 

#Permutation: free
#Number of permutations: 999

##NULL Ho: The degree of concordance between two (or more) matrices is no greater than 
###expected given random inter-matrix associations.

##So we reject Ho and accept H1  
##H1:The degree of concordance between two (or more) matrices is greater than 
###expected given random inter-matrix associations.

## Significance:  0.001 We accept the H1 hypothesis that there is a relationship among points...
## between...the position of the points in the first ordination (here 0mo) and the second (18mo).
## Priority effect

par(mar=c(5,5,2,2))
par(mfrow=c(1,1))
par(pty='s')
str(initial_log.pro2)
plot(initial_log.pro2,choices=c(1,2),main='',cex.axis=2,cex.lab=2,cex=2,lwd=6)
plot(initial_log.pro2,choices=c(1,3),main='',cex.axis=2,cex.lab=2,cex=2,lwd=6)



#par(pty='s')
#plot(initial_log.pro2, main="Procrustes test between communities at 0 mo and 18 mo")

#dev.off()



####Put the two procrustes (0 to 18 mo, 18mo to 36mo)
###Figure S5
tiff(filename="Figure S5 of beta diversity procrustes between 0 mo and 18 mo and 18 mo 36 mo.tiff", res = 600, width=7500, height=3400, compression = "lzw")
par(mfrow=c(1,2))
par(pty='s')

plot(initial_log.pro2, main="")
text(-1,1,"A")


plot(log.pro2, main="")## down.proc
text(-0.04,0.02,"B")
dev.off()



##############
####PHyla, trophic mode, rot types analyses
#####

###For only the down strict 384 and initial samples 34 so in total 418 samples
########## 418 samples (384 and initial samples 34)

rm(list=ls())

library(reshape2)
library(ggplot2)

mengsong_97pick <- read.table("mengsong_97closed_guilds_r.csv", sep=",", row.names=1,header=T, check.names=F,blank.lines.skip = FALSE)   
dim(mengsong_97pick)

###import map file
map <- read.table("mengsong_map_r_chem_SBB.csv",sep=",", header=T, check.names=F,blank.lines.skip = FALSE)
rownames(map) <- map[,1]
dim(map)
#####create subset vectors
subset_sample <- rownames(map)[which(map$control=="No")]
subset_control <- rownames(map)[which(map$control=="Yes")]

sample_map <- map[subset_sample,]
dim(sample_map)
sample_map <- sample_map[order(rownames(sample_map)),]
dim(sample_map)
sample_map$Forest <- factor(sample_map$Forest)
levels(sample_map$Forest) ###check the forest level
###check the data to make sure the otu table has the correct number of samples
mengsong_97pick[1:2,834:835]

#####create map for separating mengsong samples
Down_map <- sample_map[which(sample_map$SampleType=="Down" |sample_map$SampleType=="Initial"),]
dim(Down_map)

Down_map[] <- lapply(Down_map, function(x) if(is.factor(x)) factor(x) else x)##make sure the factor levels are consistent with the subset dataframe
levels(Down_map$Forest)

Pilot_map <- sample_map[sample_map$SampleType=="PilotSample",]

####exclude pilot map
mengsong_map <- sample_map[!sample_map$SampleType %in% "PilotSample",]
dim(mengsong_map)

mengsong_map[] <- lapply(mengsong_map, function(x) if(is.factor(x)) factor(x) else x)
mengsong_map_order <- mengsong_map[order(rownames(mengsong_map)),]
dim(mengsong_map_order)

###########filter out non fungi otus
subset_fungi <- rownames(mengsong_97pick)[which(mengsong_97pick$Kingdom=="Fungi")]
fungi_97pick <- mengsong_97pick[subset_fungi,]
dim(fungi_97pick)

fungi_97pick_mengsong <- fungi_97pick[,!colnames(fungi_97pick) %in% subset_control & !colnames(fungi_97pick) %in% Pilot_map[,1]]####only get mengsong sample otus
dim(fungi_97pick_mengsong)


fungi_97pick_mengsong[1:2,636:637]
sample_l0 <- length(mengsong_map[,1])
sample_l0

#####exclude less than 10 otus 
fungi_97pick_mengsong_11 <- subset(fungi_97pick_mengsong,rowSums(fungi_97pick_mengsong[,1:sample_l0]) >10)
fungi_97pick_mengsong_11[1:2,636:637]

dim(fungi_97pick_mengsong_11)
fungi_97pick_mengsong_11[1:2,637:652]
fungi_taxo<-fungi_97pick_mengsong_11[,637:652]

#fungi_97pick_mengsong_11 <- fungi_97pick_mengsong_11[,colnames(fungi_97pick_mengsong_11) %in% rownames(Down_map)|colnames(fungi_97pick_mengsong_11) %in% "taxonomy"]
#dim(fungi_97pick_mengsong_11)
#fungi_97pick_mengsong_11[1:2,636:486]
###################work on strict down only for publication
dim(Down_map)
Down_strict_initial_map <- Down_map[!Down_map$Position %in%"up",]
dim(Down_strict_initial_map)
remove_names <- Down_map[Down_map$Position %in%"up",]
dim(remove_names)

Down_strict_initial_map[] <- lapply(Down_strict_initial_map, function(x) if(is.factor(x)) factor(x) else x)##make sure the factor levels are cosistent with the subset dataframe
dim(Down_strict_initial_map)
str(Down_strict_initial_map)
Down_strict_initial_map_order <- Down_strict_initial_map[order(rownames(Down_strict_initial_map)),]
dim(Down_strict_initial_map_order)


####Make sure fungi_97pick_mengsong_11 is 434== 418 samples+ 16 taxonomy columns

###
dim(mengsong_map)
dim(fungi_97pick_mengsong_11)

fungi_97pick_down_strict_initial_11 <- fungi_97pick_mengsong_11[,colnames(fungi_97pick_mengsong_11) %in% rownames(Down_map)]
dim(fungi_97pick_down_strict_initial_11)

fungi_97pick_down_strict_initial_11 <- fungi_97pick_down_strict_initial_11[,!colnames(fungi_97pick_down_strict_initial_11) %in% rownames(remove_names)]
dim(fungi_97pick_down_strict_initial_11)


fungi_97pick_down_strict_initial_11<-cbind(fungi_97pick_down_strict_initial_11,fungi_taxo)
dim(fungi_97pick_down_strict_initial_11)

sample_l <- length(Down_strict_initial_map[,1])
sample_l


Phylum3_agg <- aggregate(fungi_97pick_down_strict_initial_11[,1:sample_l], by=list(fungi_97pick_down_strict_initial_11$Phylum), FUN=sum)
head(Phylum3_agg)
dim(Phylum3_agg)
names(Phylum3_agg)[1]<- paste("Phylum")
rownames(Phylum3_agg) <- Phylum3_agg[,1]
head(Phylum3_agg)
rownames(Phylum3_agg)

dim(Phylum3_agg)
head(Phylum3_agg)
Phylum3_agg[1:3,1:3]

#### NMDS  on phylum

####################################################################
####Do NMDS for 18 and 36 mo for phylum
###########################
library(vegan)
dim(Phylum3_agg)[2]
set.seed(2018)
Down_initial_phylum<-Phylum3_agg[,2:dim(Phylum3_agg)[2]]
dim(Down_initial_phylum)

##Order
Down_initial_phylum_order<-Down_initial_phylum[,order(colnames(Down_initial_phylum))]
head(Down_initial_phylum_order)

#########do adonis on strict down otu
Down_strict_dis_bray <- vegdist(t(Down_initial_phylum_order),method="bray")

#######forest and time are significant based on this model
Down_strict_NMDS<- metaMDS(t(Down_initial_phylum_order),distance="bray",k=2,trymax=10,autotransform=F,wascores=T,expand=T,trace=1,plot=F)

str(Down_strict_NMDS)
summary(Down_strict_NMDS)



Down_strict_NMDS

dim(Down_initial_phylum)

##Convergence at k=2 stress=0.1045
str(Down_strict_NMDS)
Down_strict_nmds1 <- data.frame(MDS1=Down_strict_NMDS$points[,1],MDS2=Down_strict_NMDS$points[,2])


match(rownames(Down_strict_nmds1),rownames(Down_strict_initial_map))

Down_strict_nmds_map <-cbind(Down_strict_nmds1,Down_strict_initial_map)

Down_strict_nmds_map$Position
Down_strict_nmds_map$Position <- factor(Down_strict_nmds_map$Position,levels=c("Initial","down"))

###Need to create centroids
###Subsetting
#Litsea_down_18_36mo <- subset(preddat_orig, Species == "Lit_cub"& Pos_soil=="down"&(Number_mo=="18 months"|Number_mo=="36 months"))


str(Down_strict_nmds_map)
levels(Down_strict_nmds_map$Species)
Down_strict_nmds_map$Species<-as.factor(Down_strict_nmds_map$Species)
levels(Down_strict_nmds_map$Species)
levels(Down_strict_nmds_map$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(Down_strict_nmds_map$Forest)
Down_strict_nmds_map$Forest<-as.factor(Down_strict_nmds_map$Forest)
levels(Down_strict_nmds_map$Forest)
levels(Down_strict_nmds_map$Forest)<-c("Initial","Mature forest", "Open land", "Regenerating forest")
Down_strict_nmds_map$Forest<-factor(Down_strict_nmds_map$Forest, levels=c("Initial","Mature forest", "Regenerating forest","Open land"))
levels(Down_strict_nmds_map$Forest)



###Cleaned figures
####By time

NMDS_by_time<-ggplot(data=Down_strict_nmds_map,aes(x=MDS1,y=MDS2,color=Time))+
  ylab("NMDS2")+xlab("NMDS1")+geom_point(size=2)+
  facet_wrap(~ Species,ncol=2)+theme(strip.text.x = element_text(size = 20, colour = "black"))+
  theme_bw()+theme(panel.border = element_rect(fill=NA,color="black"))+
  theme(panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 20),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', size=20),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.direction="vertical",legend.title=element_blank())+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  theme(strip.text = element_text(colour = "black", size = 20))+
  scale_color_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo Exposure    ", "18 mo Exposure   ", "36 mo Exposure    "))+
  stat_ellipse(data=Down_strict_nmds_map,aes(x=MDS1,y=MDS2,group=Time),linetype=2, size=0.8)+
  annotate("text", x=0.1, y=-1.8, label="2D-Stress = 0.10,\n F (1,382) = 11.34, P = 0.001 (18 vs. 36 mo)", size = 5, fontface="plain", colour="black")

NMDS_by_time

####By forest type


NMDS_by_forest<-ggplot(data=Down_strict_nmds_map,aes(x=MDS1,y=MDS2,shape=Forest, col=Forest))+
  ylab("NMDS2")+xlab("NMDS1")+geom_point(size=2)+
  facet_wrap(~ Species,ncol=2)+theme(strip.text.x = element_text(size = 20, colour = "black"))+
  theme_bw()+theme(panel.border = element_rect(fill=NA,color="black"))+
  theme(panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 20),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', size=20),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.direction="vertical",legend.title=element_blank())+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"), legend.position = c(0.1,0.1))+
  theme(strip.text = element_text(colour = "black", size = 20))+
  stat_ellipse(data=Down_strict_nmds_map,aes(x=MDS1,y=MDS2,group=Forest),linetype=2, size=0.8)+
  annotate("text", x=0.1, y=-1.8, label="2D-Stress = 0.10,\n F (2, 381) = 2.669, P = 0.001 (18 vs. 36 mo)", size = 5, fontface="plain", colour="black")
NMDS_by_forest

####END Phylum


#### NMDS  on Class from FUNGuild

####Based on Class

Class_agg <- aggregate(fungi_97pick_down_strict_initial_11[,1:sample_l], by=list(fungi_97pick_down_strict_initial_11$Class), FUN=sum)
head(Class_agg)
dim(Class_agg)
tail(Class_agg)
Class_agg[1:3,1:3]
str(Class_agg)
names(Class_agg)[1]<- paste("Class")
rownames(Class_agg) <- Class_agg[,1]
head(Class_agg)
rownames(Class_agg)

dim(Class_agg)

####################################################################
####Do NMDS for 18 and 36 mo for class
###########################

dim(Class_agg)[2]

Down_initial_Class<-Class_agg[,2:dim(Class_agg)[2]]
dim(Down_initial_Class)

##Order
Down_initial_Class_order<-Down_initial_Class[,order(colnames(Down_initial_Class))]
head(Down_initial_Class_order)

#########do adonis on strict down otu
Down_strict_dis_bray <- vegdist(t(Down_initial_Class_order),method="bray")

#######forest and time are significant based on this model
Down_strict_NMDS<- metaMDS(t(Down_initial_Class_order),distance="bray",k=2,trymax=10,autotransform=F,wascores=T,expand=T,trace=1,plot=F)

str(Down_strict_NMDS)
summary(Down_strict_NMDS)

## No convergence for k=2

## No convergence at k=3 as well

## Move to k=4

Down_strict_NMDS<- metaMDS(t(Down_initial_Class_order),distance="bray",k=4,trymax=200,autotransform=F,wascores=T,expand=T,trace=1,plot=F)
Down_strict_NMDS

##Convergence at k=4 stress=0.11
str(Down_strict_NMDS)

Down_strict_nmds1 <- data.frame(MDS1=Down_strict_NMDS$points[,1],MDS2=Down_strict_NMDS$points[,2])


match(rownames(Down_strict_nmds1),rownames(Down_strict_initial_map))

Down_strict_nmds_map <-cbind(Down_strict_nmds1,Down_strict_initial_map) ##Down_strict_map

Down_strict_nmds_map$Position
Down_strict_nmds_map$Position <- factor(Down_strict_nmds_map$Position,levels=c("Initial","down"))

###Need to create centroids
###Subsetting
#Litsea_down_18_36mo <- subset(preddat_orig, Class == "Lit_cub"& Pos_soil=="down"&(Number_mo=="18 months"|Number_mo=="36 months"))


str(Down_strict_nmds_map)
levels(Down_strict_nmds_map$Species)
Down_strict_nmds_map$Species<-as.factor(Down_strict_nmds_map$Species)
levels(Down_strict_nmds_map$Species)
levels(Down_strict_nmds_map$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(Down_strict_nmds_map$Species)

levels(Down_strict_nmds_map$Forest)
levels(Down_strict_nmds_map$Forest)<-c("Initial","Mature forest", "Open land", "Regenerating forest")
Down_strict_nmds_map$Forest<-factor(Down_strict_nmds_map$Forest, levels=c("Initial","Mature forest", "Regenerating forest", "Open land"))
levels(Down_strict_nmds_map$Forest)



###Cleaned figures for Figure 6
####By time

NMDS_by_time<-ggplot(data=Down_strict_nmds_map,aes(x=MDS1,y=MDS2,color=Time))+
  ylab("NMDS2")+xlab("NMDS1")+geom_point(size=2)+
  facet_wrap(~ Species,ncol=2)+theme(strip.text.x = element_text(size = 20, colour = "black"))+
  theme_bw()+theme(panel.border = element_rect(fill=NA,color="black"))+
  theme(panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 20),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.direction="vertical",legend.title=element_blank())+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"), legend.position = c(0.9,0.18))+
  theme(strip.text = element_text(colour = "black", size = 20))+
  scale_color_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo Exposure    ", "18 mo Exposure   ", "36 mo Exposure    "))+
  stat_ellipse(data=Down_strict_nmds_map,aes(x=MDS1,y=MDS2,group=Time),linetype=2, size=0.8)+
  annotate("text", x=0, y=-2.3, label="4D-Stress = 0.11, F (1,382) = 11.34, P = 0.001 (18 vs. 36mo)", size = 5, fontface="plain", colour="black")

NMDS_by_time

####By forest type


NMDS_by_forest<-ggplot(data=Down_strict_nmds_map,aes(x=MDS1,y=MDS2,shape=Forest, col=Forest))+
  ylab("NMDS2")+xlab("NMDS1")+geom_point(size=2)+
  facet_wrap(~ Species,ncol=2)+theme(strip.text.x = element_text(size = 20, colour = "black"))+
  theme_bw()+theme(panel.border = element_rect(fill=NA,color="black"))+
  theme(panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 20),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', size=20),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.direction="vertical",legend.title=element_blank())+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"), legend.position = c(0.9,0.18))+
  theme(strip.text = element_text(colour = "black", size = 20))+
  stat_ellipse(data=Down_strict_nmds_map,aes(x=MDS1,y=MDS2,group=Forest),linetype=2, size=0.8)+
  annotate("text", x=0, y=-2.3, label="4D-Stress = 0.11, F (2, 381) = 2.669, P = 0.001 (18 vs. 36mo)", size = 5, fontface="plain", colour="black")
  ##\n this is to make a new line
NMDS_by_forest


library(gridExtra)
library(reshape)
library(grid)

##Figure 6
tiff(filename="Figure 6 of Classes nmds with ellipses with k 4.tiff", res = 600, width=8200, height=8500, compression = "lzw")
grid.arrange(arrangeGrob(NMDS_by_time + theme(plot.margin = unit( c(0,1,0,0) , units = "lines")), ##top,right,down,left
                         NMDS_by_forest + theme(plot.margin = unit( c(0,1,0,0) , units = "lines")),
                         nrow=2))


dev.off()

####END Class


#### NMDS  on Order from FUNGuild

####Based on Order

Order_agg <- aggregate(fungi_97pick_down_strict_initial_11[,1:sample_l], by=list(fungi_97pick_down_strict_initial_11$Order), FUN=sum)
head(Order_agg)
dim(Order_agg)
tail(Order_agg)
Order_agg[1:3,1:3]
str(Order_agg)
names(Order_agg)[1]<- paste("Order")
rownames(Order_agg) <- Order_agg[,1]
head(Order_agg)
rownames(Order_agg)

dim(Order_agg)

####################################################################
####Do NMDS for 18 and 36 mo for order
###########################

dim(Order_agg)[2]

Down_initial_Order<-Order_agg[,2:dim(Order_agg)[2]]
dim(Down_initial_Order)

##Order
Down_initial_Order_order<-Down_initial_Order[,order(colnames(Down_initial_Order))]
head(Down_initial_Order_order)

#########do adonis on strict down otu
Down_strict_dis_bray <- vegdist(t(Down_initial_Order_order),method="bray")

#######forest and time are significant based on this model
Down_strict_NMDS<- metaMDS(t(Down_initial_Order_order),distance="bray",k=2,trymax=10,autotransform=F,wascores=T,expand=T,trace=1,plot=F)

str(Down_strict_NMDS)
summary(Down_strict_NMDS)
Down_strict_NMDS

## No convergence for k=2


## Move to k=4

Down_strict_NMDS<- metaMDS(t(Down_initial_Order_order),distance="bray",k=4,trymax=200,autotransform=F,wascores=T,expand=T,trace=1,plot=F)
Down_strict_NMDS
###No convergence after k=4
####END



#### NMDS  on Family from FUNGuild

####Based on Family

Family_agg <- aggregate(fungi_97pick_down_strict_initial_11[,1:sample_l], by=list(fungi_97pick_down_strict_initial_11$Family), FUN=sum)
head(Family_agg)
dim(Family_agg)
tail(Family_agg)
Family_agg[1:3,1:3]
str(Family_agg)
names(Family_agg)[1]<- paste("Family")
rownames(Family_agg) <- Family_agg[,1]
head(Family_agg)
rownames(Family_agg)

dim(Family_agg)

####################################################################
####Do NMDS for 18 and 36 mo for family
###########################

dim(Family_agg)[2]

Down_initial_Family<-Family_agg[,2:dim(Family_agg)[2]]
dim(Down_initial_Family)

##Order
Down_initial_Family_order<-Down_initial_Family[,order(colnames(Down_initial_Family))]
head(Down_initial_Family_order)

#########do adonis on strict down otu
Down_strict_dis_bray <- vegdist(t(Down_initial_Family_order),method="bray")

#######forest and time are significant based on this model
Down_strict_NMDS<- metaMDS(t(Down_initial_Family_order),distance="bray",k=2,trymax=10,autotransform=F,wascores=T,expand=T,trace=1,plot=F)

str(Down_strict_NMDS)
summary(Down_strict_NMDS)
Down_strict_NMDS

## No convergence for k=2

## Move to k=4

Down_strict_NMDS<- metaMDS(t(Down_initial_Family_order),distance="bray",k=4,trymax=200,autotransform=F,wascores=T,expand=T,trace=1,plot=F)
Down_strict_NMDS
##No convergence reached at k=4 so we stopped


#### NMDS  on Genus from FUNGuild

####Based on Genus

Genus_agg <- aggregate(fungi_97pick_down_strict_initial_11[,1:sample_l], by=list(fungi_97pick_down_strict_initial_11$Genus), FUN=sum)
head(Genus_agg)
dim(Genus_agg)
tail(Genus_agg)
Genus_agg[1:3,1:3]
str(Genus_agg)
names(Genus_agg)[1]<- paste("Genus")
rownames(Genus_agg) <- Genus_agg[,1]
head(Genus_agg)
rownames(Genus_agg)

dim(Genus_agg)

####################################################################
####Do NMDS for 18 and 36 mo for genus
###########################

dim(Genus_agg)[2]

Down_initial_Genus<-Genus_agg[,2:dim(Genus_agg)[2]]
dim(Down_initial_Genus)

##Order
Down_initial_Genus_order<-Down_initial_Genus[,order(colnames(Down_initial_Genus))]
head(Down_initial_Genus_order)

#########do adonis on strict down otu
Down_strict_dis_bray <- vegdist(t(Down_initial_Genus_order),method="bray")

#######forest and time are significant based on this model
Down_strict_NMDS<- metaMDS(t(Down_initial_Genus_order),distance="bray",k=2,trymax=10,autotransform=F,wascores=T,expand=T,trace=1,plot=F)

str(Down_strict_NMDS)
summary(Down_strict_NMDS)
Down_strict_NMDS
## No convergence for k=2

#Down_strict_nmds1 <- data.frame(MDS1=Down_strict_NMDS$points[,1],MDS2=Down_strict_NMDS$points[,2])
#Down_strict_nmds_map <-cbind(Down_strict_nmds1,Down_strict_map)


## Move to k=4

Down_strict_NMDS<- metaMDS(t(Down_initial_Genus_order),distance="bray",k=4,trymax=200,autotransform=F,wascores=T,expand=T,trace=1,plot=T)
Down_strict_NMDS
###No convergence t k=4 so w ended here
####END Genus



#####NMDS on 
####Based on Species

Sp_agg <- aggregate(fungi_97pick_down_strict_initial_11[,1:sample_l], by=list(fungi_97pick_down_strict_initial_11$Species), FUN=sum)
head(Sp_agg)
dim(Sp_agg)
tail(Sp_agg)
Sp_agg[1:3,1:3]
str(Sp_agg)
names(Sp_agg)[1]<- paste("species")
rownames(Sp_agg) <- Sp_agg[,1]
head(Sp_agg)
rownames(Sp_agg)

dim(Sp_agg)

####################################################################
####Do NMDS for 18 and 36 mo for species
###########################

dim(Sp_agg)[2]

Down_initial_Sp<-Sp_agg[,2:dim(Sp_agg)[2]]
dim(Down_initial_Sp)

##Order
Down_initial_Sp_order<-Down_initial_Sp[,order(colnames(Down_initial_Sp))]
head(Down_initial_Sp_order)

#########do adonis on strict down otu
Down_strict_dis_bray <- vegdist(t(Down_initial_Sp_order),method="bray")

#######forest and time are significant based on this model
Down_strict_NMDS<- metaMDS(t(Down_initial_Sp_order),distance="bray",k=2,trymax=10,autotransform=F,wascores=T,expand=T,trace=1,plot=F)

str(Down_strict_NMDS)
summary(Down_strict_NMDS)

## No convergence for k=2

#Down_strict_nmds1 <- data.frame(MDS1=Down_strict_NMDS$points[,1],MDS2=Down_strict_NMDS$points[,2])
#Down_strict_nmds_map <-cbind(Down_strict_nmds1,Down_strict_map)


## Move to k=4

Down_strict_NMDS<- metaMDS(t(Down_initial_Sp_order),distance="bray",k=4,trymax=200,autotransform=F,wascores=T,expand=T,trace=1,plot=T)
Down_strict_NMDS

####END Species

#####
####Based on OTUs

fungi_97pick_down_strict_initial_11[1:2, 1:4]
fungi_97pick_down_strict_initial_11$OTU_ID<-rownames(fungi_97pick_down_strict_initial_11)
otus_agg <- aggregate(fungi_97pick_down_strict_initial_11[,1:sample_l], by=list(fungi_97pick_down_strict_initial_11$OTU_ID), FUN=sum)
head(otus_agg)
dim(otus_agg)
names(otus_agg)[1]<- paste("otu")
rownames(otus_agg) <- otus_agg[,1]
head(otus_agg)
rownames(otus_agg)

dim(otus_agg)
head(otus_agg)
otus_agg[1:3,1:3]

#### NMDS  on otus

####################################################################
####Do NMDS for 18 and 36 mo for OTUs
###########################

dim(otus_agg)[2]

Down_initial_otus<-otus_agg[,2:dim(otus_agg)[2]]
dim(Down_initial_otus)

##Order
Down_initial_otus_order<-Down_initial_otus[,order(colnames(Down_initial_otus))]
head(Down_initial_otus_order)

#########do adonis on strict down otu
Down_strict_dis_bray <- vegdist(t(Down_initial_otus_order),method="bray")

#######forest and time are significant based on this model
Down_strict_NMDS<- metaMDS(t(Down_initial_otus_order),distance="bray",k=2,trymax=200,autotransform=F,wascores=T,expand=T,trace=1,plot=F)

str(Down_strict_NMDS)
summary(Down_strict_NMDS)

## No convergence for k=2

## Move to k=4

Down_strict_NMDS<- metaMDS(t(Down_initial_otus_order),distance="bray",k=4,trymax=300,autotransform=F,wascores=T,expand=T,trace=1,plot=T)
Down_strict_NMDS

###No convergence, I stopped.

####END OTUs


######


#####Barplots 

###Proportion based on just total OTUs list

###Proportion of unidentified based on OTUs presence
PHYLUM_otus<-fungi_97pick_down_strict_initial_11[,1]
###Make OTUs abundance more than 1 to be one just as present
PHYLUM_otus[PHYLUM_otus>0]<-1
###Make OTUs abundance 0 be one just as present in the list of otus
PHYLUM_otus[PHYLUM_otus==0]<-1
PHYLUM_otus[1:10]
####Everything should be one

PHYLUM_otus<-aggregate(PHYLUM_otus, by=list(fungi_97pick_down_strict_initial_11$Phylum), FUN=sum)
head(PHYLUM_otus)
names(PHYLUM_otus)[1]<- paste("PHYLUM")
rownames(PHYLUM_otus) <- PHYLUM_otus[,1]
head(PHYLUM_otus)
rownames(PHYLUM_otus)


PHYLUM_modify <- c("Ascomycota",  "Basidiomycota",  "Chytridiomycota",  "Glomeromycota",  "Rozellomycota",  "unidentified",  "Zygomycota")
PHYLUM_otus[,1] <- PHYLUM_modify
str(PHYLUM_otus)

###
PHYLUM_otus2 <- aggregate(PHYLUM_otus[,2:ncol(PHYLUM_otus)], by=list(PHYLUM_otus$PHYLUM), FUN=sum)
head(PHYLUM_otus2)
names(PHYLUM_otus2)[1]<- paste("PHYLUM")
rownames(PHYLUM_otus2) <- PHYLUM_otus2[,1]
str(PHYLUM_otus2)

########construct data frame for calculate, use the melt function
PHYLUM_otus_sample_only <- PHYLUM_otus2[,2:ncol(PHYLUM_otus2)]
PHYLUM_otus_sample_only_order <-PHYLUM_otus_sample_only

###Percent of unidentified based on OTUs
PHYLUM_otus_type_prop<-PHYLUM_otus_sample_only_order
str(PHYLUM_otus_type_prop)
PHYLUM_otus_percent<-PHYLUM_otus_type_prop*100/sum(PHYLUM_otus_type_prop)
PHYLUM_otus_per<-cbind(PHYLUM_otus2,PHYLUM_otus_percent)

sum(PHYLUM_otus_per$PHYLUM_otus_percent)

####


###Proportion of unidentified based on OTUs presence
PHYLUM_pres_abs<-fungi_97pick_down_strict_initial_11[,1:sample_l]
PHYLUM_pres_abs[PHYLUM_pres_abs>0]<-1

PHYLUM_pres_abs[1:3]
PHYLUM<-aggregate(PHYLUM_pres_abs, by=list(fungi_97pick_down_strict_initial_11$Phylum), FUN=sum)
head(PHYLUM)
names(PHYLUM)[1]<- paste("PHYLUM")
rownames(PHYLUM) <- PHYLUM[,1]
head(PHYLUM)
rownames(PHYLUM)


PHYLUM_modify <- c("Ascomycota",  "Basidiomycota",  "Chytridiomycota",  "Glomeromycota",  "Rozellomycota",  "unidentified",  "Zygomycota")
PHYLUM[,1] <- PHYLUM_modify
str(PHYLUM)

###
PHYLUM2 <- aggregate(PHYLUM[,2:ncol(PHYLUM)], by=list(PHYLUM$PHYLUM), FUN=sum)
head(PHYLUM2)
names(PHYLUM2)[1]<- paste("PHYLUM")
rownames(PHYLUM2) <- PHYLUM2[,1]
str(PHYLUM2)

########construct data frame for calculate, use the melt function
PHYLUM_sample_only <- PHYLUM2[,2:ncol(PHYLUM2)]
PHYLUM_sample_only_order <- PHYLUM_sample_only[,order(colnames(PHYLUM_sample_only))]

###Percent of unidentified based on OTUs
PHYLUM_type_prop<-PHYLUM_sample_only_order
PHYLUM_type_prop$total<-rowSums(PHYLUM_type_prop)
PHYLUM_type_prop[1:7, 1:3]
PHYLUM_type_prop$percent<-(PHYLUM_type_prop$total/sum(PHYLUM_type_prop$total) * 100)
PHYLUM_type_prop[1:7, 418:420]
sum(PHYLUM_type_prop$percent)


#####Phylum

Phylum3_agg2 <- aggregate(Phylum3_agg[,2:ncol(Phylum3_agg)], by=list(Phylum3_agg$Phylum), FUN=sum)
head(Phylum3_agg2)
names(Phylum3_agg2)[1]<- paste("Phylum")
rownames(Phylum3_agg2) <- Phylum3_agg2[,1]

########construct data frame for calculate, use the melt function
Phylum_sample_only <- Phylum3_agg2[,2:ncol(Phylum3_agg2)]
Phylum_sample_only_order <- Phylum_sample_only[,order(colnames(Phylum_sample_only))]

#head(Phylum_sample_only_order)
Phylum_sample_only_order[1:7,1:3]

Phylum_prop<-Phylum_sample_only_order
Phylum_prop$total<-rowSums(Phylum_prop)

Phylum_prop$percent<-(Phylum_prop$total/sum(Phylum_prop$total) * 100)
sum(Phylum_prop$percent)

str(Phylum_sample_only_order)


#########combine with map file
head(t(Phylum_sample_only_order))
match(rownames(t(Phylum_sample_only_order)),Down_strict_initial_map$SampleID)
###Matching not correct

Phylum_combine <- cbind(t(Phylum_sample_only_order),Down_strict_initial_map)
dim(Phylum_combine)
str(Phylum_combine)
#write.csv(Phylum_combine, file="Phylum combine august2020.csv")
melt_id1 <- nrow(Phylum_sample_only_order)+1
melt_id1
melt_id2 <- ncol(Phylum_combine)
melt_id2

library(reshape2)
##Use melt to actually put the trait under one column and add set of rows with the ith information of ith trait 
## Basically melt is collapsing several columns to two (one for the columns and one for their values)
Phylum_melt <- reshape2::melt(Phylum_combine,id=melt_id1:melt_id2,variable.name = "Phylum",value.name = "Count")
str(Phylum_melt)
dim(Phylum_melt)
#Phylum_melt[1:3,1:3]
Phylum_time <- with(Phylum_melt,aggregate(Count,by=list(Time=Time,Phylum=Phylum,Species=Species,Forest=Forest),function(x) mean(x,na.rm=T)))
dim(Phylum_time)
Try_Phylum_total<-Phylum_time
head(Try_Phylum_total)
levels(Try_Phylum_total$Time)
Try_Phylum_total$Time<-as.factor(Try_Phylum_total$Time)
levels(Try_Phylum_total$Time)

levels(Try_Phylum_total$Species)
Try_Phylum_total$Species<-as.factor(Try_Phylum_total$Species)
levels(Try_Phylum_total$Species)

levels(Try_Phylum_total$Forest)
Try_Phylum_total$Forest<-as.factor(Try_Phylum_total$Forest)
levels(Try_Phylum_total$Forest)

levels(Try_Phylum_total$Phylum)


sum(Try_Phylum_total$x)
dim(Try_Phylum_total)

str(Try_Phylum_total)
l_try <- length(levels(Try_Phylum_total$Phylum))*2

Initial_strict_down_try <- do.call("rbind", replicate(3, Try_Phylum_total[1:l_try,], simplify = FALSE))
dim(Initial_strict_down_try)

Initial_strict_down_try$Forest <- c(paste(replicate(l_try,"Mature_forest")),paste(replicate(l_try,"Open_land")),paste(replicate(l_try,"Regenerating_forest")))
dim(Initial_strict_down_try)

strict_down_try_plot <- rbind(Initial_strict_down_try,Try_Phylum_total[(l_try+1):nrow(Try_Phylum_total),])
dim(strict_down_try_plot)
strict_down_try_plot[] <- lapply(strict_down_try_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_try_plot$Forest <- factor(strict_down_try_plot$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))

strict_down_try_plot$Phylum    #### <- factor(strict_down_try_plot$Phylum,levels=c("Ascomycota",  "Basidiomycota",  "Chytridiomycota",  
dim(strict_down_try_plot)                                   ##                                        "Glomeromycota",  "Rozellomycota",  "unidentified",  "Zygomycota"))
head(strict_down_try_plot,42)

strict_down_try_plot$x_percent<-(strict_down_try_plot$x)*100/sum(strict_down_try_plot$x)
sum(strict_down_try_plot$x_percent)


tail(strict_down_try_plot)
levels(strict_down_try_plot$Species)
levels(strict_down_try_plot$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_try_plot$Forest)
levels(strict_down_try_plot$Forest)<-c("Open land", "Regenerating forest", "Mature forest")

ggplot(strict_down_try_plot, aes(x = Time, y = x_percent, fill = Phylum)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  ylab("Relative abundance (%)")+
  xlab("Exposure time (Months)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))



library(plyr)

Phylum_perct<- ddply(Try_Phylum_total, .(Forest,Time,Species), mutate, Phylum_pct = x / sum(x) * 100)
levels(Phylum_perct$Time)
l_initial <- length(levels(Phylum_perct$Phylum))*2

Initial_strict_down_rep <- do.call("rbind", replicate(3, Phylum_perct[1:l_initial,], simplify = FALSE))
Initial_strict_down_rep$Forest <- c(paste(replicate(l_initial,"Mature_forest")),paste(replicate(l_initial,"Open_land")),paste(replicate(l_initial,"Regenerating_forest")))

strict_down_perct_plot <- rbind(Initial_strict_down_rep,Phylum_perct[(l_initial+1):nrow(Phylum_perct),])
strict_down_perct_plot[] <- lapply(strict_down_perct_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_perct_plot$Forest <- factor(strict_down_perct_plot$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))

strict_down_perct_plot$Phylum ###  #### <- factor(strict_down_try_plot$Phylum,levels=c("Ascomycota",  "Basidiomycota",  "Chytridiomycota",  
##                                        "Glomeromycota",  "Rozellomycota",  "unidentified",  "Zygomycota"))

levels(strict_down_perct_plot$Species)
levels(strict_down_perct_plot$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_perct_plot$Forest)
levels(strict_down_perct_plot$Forest)<-c("Open land", "Regenerating forest", "Mature forest")

ggplot(strict_down_perct_plot, aes(x = Time, y = Phylum_pct, fill = Phylum)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  xlab("Exposure time (Months)")+ ylab("Abundance (%)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))

###
#################exclude unidentified  from bar graph
head(Phylum3_agg2)
rownames(Phylum3_agg2)

Phylum_agg3 <- Phylum3_agg2[-c(6),]
dim(Phylum_agg3)
########construct data frame for calculate, use the melt function
Phylum_sample_only2 <- Phylum_agg3[,2:ncol(Phylum_agg3)]
Phylum_sample_only_order2 <- Phylum_sample_only2[,order(colnames(Phylum_sample_only2))]

#########combine with map file
head(t(Phylum_sample_only_order2))
match(rownames(t(Phylum_sample_only_order2)),Down_strict_initial_map$SampleID)

Phylum_combine2 <- cbind(t(Phylum_sample_only_order2),Down_strict_initial_map)
dim(Phylum_combine2)

melt_id1 <- nrow(Phylum_sample_only_order2)+1
melt_id2 <- ncol(Phylum_combine2)
Phylum_melt2 <- reshape2::melt(Phylum_combine2,id=melt_id1:melt_id2,variable.name = "Phylum",value.name = "Count")
dim(Phylum_melt2)

Phylum_time2 <- with(Phylum_melt2,aggregate(Count,by=list(Time=Time,Phylum=Phylum,Species=Species,Forest=Forest),function(x) mean(x,na.rm=T)))
str(Phylum_time2)

Try_Phylum_total2<-Phylum_time2
head(Try_Phylum_total2)
sum(Try_Phylum_total2$x)


l_try <- length(levels(Try_Phylum_total2$Phylum))*2

Try_Phylum_total2$Phylum
Initial_strict_down_try2 <- do.call("rbind", replicate(3, Try_Phylum_total2[1:l_try,], simplify = FALSE))
Initial_strict_down_try2$Forest <- c(paste(replicate(l_try,"Mature_forest")),paste(replicate(l_try,"Open_land")),paste(replicate(l_try,"Regenerating_forest")))


strict_down_try_plot2 <- rbind(Initial_strict_down_try2,Try_Phylum_total2[(l_try+1):nrow(Try_Phylum_total2),])
strict_down_try_plot2[] <- lapply(strict_down_try_plot2, function(x) if(is.factor(x)) factor(x) else x)
strict_down_try_plot2$Forest <- factor(strict_down_try_plot2$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))

strict_down_try_plot2$Phylum ###<- factor(strict_down_try_plot2$Phylum,levels=c( "Ascomycota" , "Basidiomycota",  "Chytridiomycota",
                                ####"Glomeromycota" , "Rozellomycota" , "Zygomycota"))

strict_down_try_plot2$x_percent<-(strict_down_try_plot2$x)*100/sum(strict_down_try_plot2$x)
sum(strict_down_try_plot2$x_percent)

levels(strict_down_try_plot2$Species)
strict_down_try_plot2$Species<-as.factor(strict_down_try_plot2$Species)
levels(strict_down_try_plot2$Species)
levels(strict_down_try_plot2$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_try_plot2$Species)

levels(strict_down_try_plot2$Forest)
levels(strict_down_try_plot2$Forest)<-c("Open land", "Regenerating forest", "Mature forest")

phylum_all_no_unid<-ggplot(strict_down_try_plot2, aes(x = Time, y = x_percent, fill = Phylum)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  ylab("Relative abundance (%)")+
  xlab("Exposure time (months)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))
phylum_all_no_unid


Phylum_perct2<- ddply(Phylum_time2, .(Forest,Time,Species), mutate, Phylum_pct2 = x / sum(x) * 100)
levels(Phylum_perct2$Time)
l_initial <- length(levels(Phylum_perct2$Phylum))*2

Initial_strict_down_rep2 <- do.call("rbind", replicate(3, Phylum_perct2[1:l_initial,], simplify = FALSE))
Initial_strict_down_rep2$Forest <- c(paste(replicate(l_initial,"Mature_forest")),paste(replicate(l_initial,"Open_land")),paste(replicate(l_initial,"Regenerating_forest")))

strict_down_perct_plot2 <- rbind(Initial_strict_down_rep2,Phylum_perct2[(l_initial+1):nrow(Phylum_perct2),])
strict_down_perct_plot2[] <- lapply(strict_down_perct_plot2, function(x) if(is.factor(x)) factor(x) else x)
strict_down_perct_plot2$Forest <- factor(strict_down_perct_plot2$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))

strict_down_perct_plot2$Phylum ###<- factor(strict_down_try_plot2$Phylum,levels=c( "Ascomycota" , "Basidiomycota",  "Chytridiomycota",
####"Glomeromycota" , "Rozellomycota" , "Zygomycota"))
levels(strict_down_perct_plot2$Species)
strict_down_perct_plot2$Species<-as.factor(strict_down_perct_plot2$Species)
levels(strict_down_perct_plot2$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_perct_plot2$Species)
levels(strict_down_perct_plot2$Forest)
levels(strict_down_perct_plot2$Forest)<-c("Open land", "Regenerating forest", "Mature forest")
levels(strict_down_perct_plot2$Forest)

phylum_all_no_unid_100<-ggplot(strict_down_perct_plot2, aes(x = Time, y = Phylum_pct2, fill = Phylum)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3)+
  ylab("Abundance (%)")+
  xlab("Exposure time (months)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))
phylum_all_no_unid_100



###Consider putting Asco and Basidio and anyhing else as "others"
head(Phylum3_agg2)
rownames(Phylum3_agg2)
str(Phylum3_agg2)

####modify Phylum category to fit our three phyla groups of interest "Ascomycota" , "Basidiomycota",  and "Others" 
Phylum_modify <- c( "Ascomycota" , "Basidiomycota", "Others","Others" , "Others" , " unidentified", "Others")
Phylum3_agg2[,1] <- Phylum_modify

Phylum_agg4 <- Phylum3_agg2[-c(6),]
str(Phylum_agg4)

head(Phylum_agg4,2,3)
Phylum_agg4$Phylum<-as.factor(Phylum_agg4$Phylum)
Phylum_agg4 <- aggregate(Phylum_agg4[,2:ncol(Phylum_agg4)], by=list(Phylum_agg4$Phylum), FUN=sum)
str(Phylum_agg4)
names(Phylum_agg4)[1]<- paste("Phylum")
rownames(Phylum_agg4)
rownames(Phylum_agg4) <- c("Ascomycota", "Basidiomycota","Others")#Trait_agg2[,1]

########construct data frame for calculate, use the melt function
Phylum_sample_only3 <- Phylum_agg4[,2:ncol(Phylum_agg4)]
Phylum_sample_only_order3 <- Phylum_sample_only3[,order(colnames(Phylum_sample_only3))]

#########combine with map file
head(t(Phylum_sample_only_order3))
match(rownames(t(Phylum_sample_only_order3)),Down_strict_initial_map$SampleID)
####The matching is correct

Phylum_combine3 <- cbind(t(Phylum_sample_only_order3),Down_strict_initial_map)
dim(Phylum_combine3)

melt_id1 <- nrow(Phylum_sample_only_order3)+1
melt_id1
melt_id2 <- ncol(Phylum_combine3)
melt_id2
Phylum_melt3 <- reshape2::melt(Phylum_combine3,id=melt_id1:melt_id2,variable.name = "Phylum",value.name = "Count")

Phylum_time3 <- with(Phylum_melt3,aggregate(Count,by=list(Time=Time,Phylum=Phylum,Species=Species,Forest=Forest),function(x) mean(x,na.rm=T)))
str(Phylum_time3)
Phylum_time3$Phylum

Try_Phylum_total3<-Phylum_time3
Try_Phylum_total3$Phylum

head(Try_Phylum_total3)
sum(Try_Phylum_total3$x)
dim(Try_Phylum_total3)
Try_Phylum_total3

l_try <- length(levels(Try_Phylum_total3$Phylum))*2
Initial_strict_down_try3 <-Try_Phylum_total3
Initial_strict_down_try3 <- do.call("rbind", replicate(3, Try_Phylum_total3[1:l_try,], simplify = FALSE))
Initial_strict_down_try3$Forest <- c(paste(replicate(l_try,"Mature_forest")),paste(replicate(l_try,"Open_land")),paste(replicate(l_try,"Regenerating_forest")))


strict_down_try_plot3 <- rbind(Initial_strict_down_try3,Try_Phylum_total3[(l_try+1):nrow(Try_Phylum_total3),])
strict_down_try_plot3[] <- lapply(strict_down_try_plot3, function(x) if(is.factor(x)) factor(x) else x)
strict_down_try_plot3$Forest <- factor(strict_down_try_plot3$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))

strict_down_try_plot3$x_percent<-(strict_down_try_plot3$x)*100/sum(strict_down_try_plot3$x)
sum(strict_down_try_plot3$x_percent)
dim(strict_down_try_plot3)

strict_down_try_plot3$Phylum ###<- factor(strict_down_try_plot3$Phylum,levels=c( "Ascomycota" , "Basidiomycota",  "Chytridiomycota",
####"Others"))
levels(strict_down_try_plot3$Species)
strict_down_try_plot3$Species<-as.factor(strict_down_try_plot3$Species)
levels(strict_down_try_plot3$Species)
levels(strict_down_try_plot3$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_try_plot3$Forest)
levels(strict_down_try_plot3$Forest)<-c("Open land", "Regenerating forest", "Mature forest")
levels(strict_down_try_plot3$Forest)

phylum_Asc_Bas_Oth_not100<-ggplot(strict_down_try_plot3, aes(x = Time, y = x_percent, fill = Phylum)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  ylab("Relative abundance (%)")+
  xlab("Exposure time (months)")+
  #annotate("text", x=0.5, y=5, label= "1209",size=5)+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 20),axis.text.y=element_text(colour = 'black', size = 20))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=20),axis.text.x = element_text(angle=0,colour = 'black', size=20,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=20), legend.title=element_text(size=20),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))+
  scale_fill_brewer(palette = "Paired")### With this you can change to colors needed
  
phylum_Asc_Bas_Oth_not100

strict_down_try_plot3

Phylum_perct3<- ddply(Phylum_time3, .(Forest,Time,Species), mutate, Phylum_pct3 = x / sum(x) * 100)
levels(Phylum_perct3$Time)
l_initial <- length(levels(Phylum_perct3$Phylum))*2

Initial_strict_down_rep3 <- do.call("rbind", replicate(3, Phylum_perct3[1:l_initial,], simplify = FALSE))
Initial_strict_down_rep3$Forest <- c(paste(replicate(l_initial,"Mature_forest")),paste(replicate(l_initial,"Open_land")),paste(replicate(l_initial,"Regenerating_forest")))

strict_down_perct_plot3 <- rbind(Initial_strict_down_rep3,Phylum_perct3[(l_initial+1):nrow(Phylum_perct3),])
strict_down_perct_plot3[] <- lapply(strict_down_perct_plot3, function(x) if(is.factor(x)) factor(x) else x)
strict_down_perct_plot3$Forest <- factor(strict_down_perct_plot3$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))

strict_down_perct_plot3$Phylum ###<- factor(strict_down_try_plot3$Phylum,levels=c( "Ascomycota" , "Basidiomycota",  "Chytridiomycota",
####"Others" ))
levels(strict_down_perct_plot3$Species)
levels(strict_down_perct_plot3$Species)<-c("Castanopsis mekongensis","Litsea cubeba")

levels(strict_down_perct_plot3$Forest)<-c("Open land", "Regenerating forest", "Mature forest")
levels(strict_down_perct_plot3$Forest)


phylum_Asc_Bas_Oth_100<-ggplot(strict_down_perct_plot3, aes(x = Time, y = Phylum_pct3, fill = Phylum)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3)+
  ylab("Abundance (%)")+
  xlab("Exposure time (months)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 20),axis.text.y=element_text(colour = 'black', size = 20))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=20),axis.text.x = element_text(angle=0,colour = 'black', size=20,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=20), legend.title=element_text(size=20),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))+
  scale_fill_brewer(palette = "Paired")
phylum_Asc_Bas_Oth_100



###

### Traits on rot types


###Proportion based on just total OTUs list

###Proportion of unidentified based on OTUs presence
ROT_otus<-fungi_97pick_down_strict_initial_11[,1]
###Make OTUs abundance more than 1 to be one just as present
ROT_otus[ROT_otus>0]<-1
###Make OTUs abundance 0 be one just as present in the list of otus
ROT_otus[ROT_otus==0]<-1
ROT_otus[1:10]
####Everything should be one

ROT_otus<-aggregate(ROT_otus, by=list(fungi_97pick_down_strict_initial_11$Trait), FUN=sum)
head(ROT_otus)
names(ROT_otus)[1]<- paste("ROT")
rownames(ROT_otus) <- ROT_otus[,1]
head(ROT_otus)
rownames(ROT_otus)


ROT_modify <- c("Unassigned","Blue-Staining", "Brown Rot","Hypogeous","Unassigned","Poisonous", "Soft Rot", "Soft Rot", "White Rot")
ROT_otus[,1] <- ROT_modify
str(ROT_otus)

###
ROT_otus2 <- aggregate(ROT_otus[,2:ncol(ROT_otus)], by=list(ROT_otus$ROT), FUN=sum)
head(ROT_otus2)
names(ROT_otus2)[1]<- paste("ROT")
rownames(ROT_otus2) <- ROT_otus2[,1]
str(ROT_otus2)

########construct data frame for calculate, use the melt function
ROT_otus_sample_only <- ROT_otus2[,2:ncol(ROT_otus2)]
ROT_otus_sample_only_order <-ROT_otus_sample_only

###Percent of unidentified based on OTUs
ROT_otus_type_prop<-ROT_otus_sample_only_order
str(ROT_otus_type_prop)
ROT_otus_percent<-ROT_otus_type_prop*100/sum(ROT_otus_type_prop)
ROT_otus_per<-cbind(ROT_otus2,ROT_otus_percent)
ROT_otus_per
sum(ROT_otus_per$ROT_otus_percent)


###Proportion of unidentified based on OTUs presence absence across all samples
Rot_pres_abs<-fungi_97pick_down_strict_initial_11[,1:sample_l]
Rot_pres_abs[Rot_pres_abs>0]<-1

Rot_pres_abs[1:3]
ROT<-aggregate(Rot_pres_abs, by=list(fungi_97pick_down_strict_initial_11$Trait), FUN=sum)
head(ROT)
names(ROT)[1]<- paste("Rot")
rownames(ROT) <- ROT[,1]
head(ROT)
rownames(ROT)

ROT_modify <- c("Unassigned","Blue-Staining", "Brown Rot","Hypogeous","Unassigned","Poisonous", "Soft Rot", "Soft Rot", "White Rot" )
ROT[,1] <- ROT_modify
str(ROT)

###
ROT2 <- aggregate(ROT[,2:ncol(ROT)], by=list(ROT$Rot), FUN=sum)
head(ROT2)
names(ROT2)[1]<- paste("Rot")
rownames(ROT2) <- ROT2[,1]
str(ROT2)




########construct data frame for calculate, use the melt function
ROT_sample_only <- ROT2[,2:ncol(ROT2)]
ROT_sample_only_order <- ROT_sample_only[,order(colnames(ROT_sample_only))]

###Percent of unidentified based on OTUs
ROT_type_prop<-ROT_sample_only_order
ROT_type_prop$total<-rowSums(ROT_type_prop)
ROT_type_prop[1:7, 1:3]
ROT_type_prop$percent<-(ROT_type_prop$total/sum(ROT_type_prop$total) * 100)
ROT_type_prop[1:7, 418:420]
sum(ROT_type_prop$percent)




##ROT TYPES
Rot_agg <- aggregate(fungi_97pick_down_strict_initial_11[,1:sample_l], by=list(fungi_97pick_down_strict_initial_11$Trait), FUN=sum)
head(Rot_agg)
names(Rot_agg)[1]<- paste("Rot")
rownames(Rot_agg) <- Rot_agg[,1]
head(Rot_agg)
rownames(Rot_agg)

####modify trophic 
#Trophic_modify <- c("Unassigned","Blue-Staining", "Brown Rot","Hypogeous","NULL","Poisonous", "Soft Rot", "Soft Rot fungus (Seehann et al. 1975)", "White Rot" )

##
Rot_modify <- c("Unassigned","Blue-Staining", "Brown Rot","Hypogeous","Unassigned","Poisonous", "Soft Rot", "Soft Rot", "White Rot" )
Rot_agg[,1] <- Rot_modify
str(Rot_agg)

###
Rot_agg2 <- aggregate(Rot_agg[,2:ncol(Rot_agg)], by=list(Rot_agg$Rot), FUN=sum)
head(Rot_agg2)
names(Rot_agg2)[1]<- paste("Rot")
rownames(Rot_agg2) <- Rot_agg2[,1]
str(Rot_agg2)




########construct data frame for calculate, use the melt function
Rot_sample_only <- Rot_agg2[,2:ncol(Rot_agg2)]
Rot_sample_only_order <- Rot_sample_only[,order(colnames(Rot_sample_only))]

###Percent of unidentified
Rot_type_prop<-Rot_sample_only_order
Rot_type_prop$total<-rowSums(Rot_type_prop)
Rot_type_prop[1:7, 1:3]
Rot_type_prop$percent<-(Rot_type_prop$total/sum(Rot_type_prop$total) * 100)
Rot_type_prop[1:7, 418:420]
sum(Rot_type_prop$percent)


#########combine with map file
head(t(Rot_sample_only_order))
match(rownames(t(Rot_sample_only_order)),Down_strict_initial_map$SampleID)
##Matching is corect

Rot_combine <- cbind(t(Rot_sample_only_order),Down_strict_initial_map)
dim(Rot_combine)
melt_id1 <- nrow(Rot_sample_only_order)+1
melt_id1
melt_id2 <- ncol(Rot_combine)
melt_id2

library(reshape2)
##Use melt to actually put the trait under one column and add set of rows with the ith information of ith trait 
## Basically melt is collapsing several columns to two (one for the columns and one for their values)
Rot_melt <- reshape2::melt(Rot_combine,id=melt_id1:melt_id2,variable.name = "Rot",value.name = "Count")
str(Rot_melt)
summary(Rot_melt)

Rot_time <- with(Rot_melt,aggregate(Count,by=list(Time=Time,Rot=Rot,Species=Species,Forest=Forest),function(x) mean(x,na.rm=T)))
Try_rot_total<-Rot_time

head(Try_rot_total)
sum(Try_rot_total$x)


l_try <- length(levels(Try_rot_total$Rot))*2

Initial_strict_down_try <- do.call("rbind", replicate(3, Try_rot_total[1:l_try,], simplify = FALSE))
Initial_strict_down_try$Forest <- c(paste(replicate(l_try,"Mature_forest")),paste(replicate(l_try,"Open_land")),paste(replicate(l_try,"Regenerating_forest")))


strict_down_try_plot <- rbind(Initial_strict_down_try,Try_rot_total[(l_try+1):nrow(Try_rot_total),])
strict_down_try_plot[] <- lapply(strict_down_try_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_try_plot$Forest <- factor(strict_down_try_plot$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))

strict_down_try_plot$Rot
strict_down_try_plot$Rot <- factor(strict_down_try_plot$Rot,levels=c("Unassigned","Blue-Staining", "Hypogeous", "Poisonous","Brown Rot", 
                                                                     "Soft Rot",  "White Rot"))
strict_down_try_plot$x_percent<-(strict_down_try_plot$x)*100/sum(strict_down_try_plot$x)
sum(strict_down_try_plot$x_percent)

levels(strict_down_try_plot$Species)
levels(strict_down_try_plot$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_try_plot$Forest)
levels(strict_down_try_plot$Forest)<-c("Open land", "Regenerating forest", "Mature forest")

ggplot(strict_down_try_plot, aes(x = Time, y = x_percent, fill = Rot)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  ylab("Relative abundance (%)")+
  xlab("Exposure time (Months)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))



#library(plyr)

Rot_perct<- ddply(Rot_time, .(Forest,Time,Species), mutate, Rot_pct = x / sum(x) * 100)
levels(Rot_perct$Time)
l_initial <- length(levels(Rot_perct$Rot))*2

Initial_strict_down_rep <- do.call("rbind", replicate(3, Rot_perct[1:l_initial,], simplify = FALSE))
Initial_strict_down_rep$Forest <- c(paste(replicate(l_initial,"Mature_forest")),paste(replicate(l_initial,"Open_land")),paste(replicate(l_initial,"Regenerating_forest")))

strict_down_perct_plot <- rbind(Initial_strict_down_rep,Rot_perct[(l_initial+1):nrow(Rot_perct),])
strict_down_perct_plot[] <- lapply(strict_down_perct_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_perct_plot$Forest <- factor(strict_down_perct_plot$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))
levels(strict_down_perct_plot$Forest) <- c("Open land","Regenerating forest","Mature forest")
levels(strict_down_perct_plot$Forest)

levels(strict_down_perct_plot$Species)
strict_down_perct_plot$Species<-as.factor(strict_down_perct_plot$Species)
levels(strict_down_perct_plot$Species)
levels(strict_down_perct_plot$Species) <- c("Castanopsis mekongensis","Litsea_cubeba")
levels(strict_down_perct_plot$Species)
strict_down_perct_plot$Rot <- factor(strict_down_perct_plot$Rot,levels=c("Unassigned","Blue-Staining", "Hypogeous", "Poisonous","Brown Rot", 
                                                                         "Soft Rot",  "White Rot"))

rot_all_una_not100<-ggplot(strict_down_perct_plot, aes(x = Time, y = Rot_pct, fill = Rot)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  xlab("Exposure time (Months)")+ ylab("Abundance (%)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 20),axis.text.y=element_text(colour = 'black', size = 20))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=20),axis.text.x = element_text(angle=0,colour = 'black', size=20,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=20), legend.title=element_text(size=20),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))+
  scale_fill_brewer(palette = "Paired")
rot_all_una_not100



###
#################exclude unassigned and other (not rot types) from bar graph
head(Rot_agg2)
rownames(Rot_agg2)
Rot_agg3 <- Rot_agg2[-c(1,3,4,6),]
########construct data frame for calculate, use the melt function
Rot_sample_only2 <- Rot_agg3[,2:ncol(Rot_agg3)]
Rot_sample_only_order2 <- Rot_sample_only2[,order(colnames(Rot_sample_only2))]

#########combine with map file
head(t(Rot_sample_only_order2))
match(rownames(t(Rot_sample_only_order2)),Down_strict_initial_map$SampleID)
###Matching is correct

Rot_combine2 <- cbind(t(Rot_sample_only_order2),Down_strict_initial_map)
dim(Rot_combine2)

melt_id1 <- nrow(Rot_sample_only_order2)+1
melt_id2 <- ncol(Rot_combine2)
Rot_melt2 <- reshape2::melt(Rot_combine2,id=melt_id1:melt_id2,variable.name = "Rot",value.name = "Count")
head(Rot_melt2)
Rot_time2 <- with(Rot_melt2,aggregate(Count,by=list(Time=Time,Rot=Rot,Species=Species,Forest=Forest),function(x) mean(x,na.rm=T)))

Try_rot_total2<-Rot_time2

head(Try_rot_total2)
sum(Try_rot_total2$x)

l_try <- length(levels(Try_rot_total2$Rot))*2

Try_rot_total2$Rot
Initial_strict_down_try2 <- do.call("rbind", replicate(3, Try_rot_total2[1:l_try,], simplify = FALSE))
Initial_strict_down_try2$Forest <- c(paste(replicate(l_try,"Mature_forest")),paste(replicate(l_try,"Open_land")),paste(replicate(l_try,"Regenerating_forest")))


strict_down_try_plot2 <- rbind(Initial_strict_down_try2,Try_rot_total2[(l_try+1):nrow(Try_rot_total2),])
strict_down_try_plot2[] <- lapply(strict_down_try_plot2, function(x) if(is.factor(x)) factor(x) else x)
strict_down_try_plot2$Forest <- factor(strict_down_try_plot2$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))

strict_down_try_plot2$Rot <- factor(strict_down_try_plot2$Rot,levels=c("Brown Rot", "Soft Rot", "White Rot"))
strict_down_try_plot2$Rot <- factor(strict_down_try_plot2$Rot,levels=c("White Rot", "Soft Rot", "Brown Rot"))
strict_down_try_plot2$x_percent<-(strict_down_try_plot2$x)*100/sum(strict_down_try_plot2$x)


levels(strict_down_try_plot2$Species)
strict_down_try_plot2$Species<-as.factor(strict_down_try_plot2$Species)
levels(strict_down_try_plot2$Species)
levels(strict_down_try_plot2$Species)<-c("Castanopsis mekongensis","Litsea cubeba")


levels(strict_down_try_plot2$Forest)
levels(strict_down_try_plot2$Forest)<-c("Open land", "Regenerating forest", "Mature forest")
levels(strict_down_try_plot2$Forest)

str(strict_down_try_plot2)
rot_bro_sof_whi_not100<-ggplot(strict_down_try_plot2, aes(x = Time, y = x_percent, fill = Rot)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  ylab("Relative abundance (%)")+
  xlab("Exposure time (months)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 20),axis.text.y=element_text(colour = 'black', size = 20))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=20),axis.text.x = element_text(angle=0,colour = 'black', size=20,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=20), legend.title=element_text(size=20),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))+
  scale_fill_brewer(palette = "Paired")
rot_bro_sof_whi_not100

dim(strict_down_try_plot2)



Rot_perct2<- ddply(Rot_time2, .(Forest,Time,Species), mutate, Rot_pct2 = x / sum(x) * 100)
levels(Rot_perct2$Time)
l_initial <- length(levels(Rot_perct2$Rot))*2

Initial_strict_down_rep2 <- do.call("rbind", replicate(3, Rot_perct2[1:l_initial,], simplify = FALSE))
Initial_strict_down_rep2$Forest <- c(paste(replicate(l_initial,"Mature_forest")),paste(replicate(l_initial,"Open_land")),paste(replicate(l_initial,"Regenerating_forest")))

strict_down_perct_plot2 <- rbind(Initial_strict_down_rep2,Rot_perct2[(l_initial+1):nrow(Rot_perct2),])
strict_down_perct_plot2[] <- lapply(strict_down_perct_plot2, function(x) if(is.factor(x)) factor(x) else x)
strict_down_perct_plot2$Forest <- factor(strict_down_perct_plot2$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))

strict_down_perct_plot2$Rot <- factor(strict_down_perct_plot2$Rot,levels=c("Brown Rot", "Soft Rot", "White Rot"))
strict_down_perct_plot2$Rot <- factor(strict_down_perct_plot2$Rot,levels=c("White Rot", "Soft Rot", "Brown Rot"))

levels(strict_down_perct_plot2$Species)
strict_down_perct_plot2$Species<-as.factor((strict_down_perct_plot2$Species))
levels(strict_down_perct_plot2$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_perct_plot2$Species)
levels(strict_down_perct_plot2$Forest)
levels(strict_down_perct_plot2$Forest)<-c("Open land", "Regenerating forest", "Mature forest")
levels(strict_down_perct_plot2$Forest)


rot_bro_sof_whi_100<-ggplot(strict_down_perct_plot2, aes(x = Time, y = Rot_pct2, fill = Rot)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3)+
  ylab("Abundance (%)")+
  xlab("Exposure time (months)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 20),axis.text.y=element_text(colour = 'black', size = 20))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=20),axis.text.x = element_text(angle=0,colour = 'black', size=20,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=20), legend.title=element_text(size=20),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))+
  scale_fill_brewer(palette = "Paired")
rot_bro_sof_whi_100


####

####Trophic traits

###Proportion based on just total OTUs list

###Proportion of unidentified based on OTUs presence
TROPHIC_otus<-fungi_97pick_down_strict_initial_11[,1]
###Make OTUs abundance more than 1 to be one just as present
TROPHIC_otus[TROPHIC_otus>0]<-1
###Make OTUs abundance 0 be one just as present in the list of otus
TROPHIC_otus[TROPHIC_otus==0]<-1
TROPHIC_otus[1:10]
####Everything should be one

TROPHIC_otus<-aggregate(TROPHIC_otus, by=list(fungi_97pick_down_strict_initial_11$Trophic_Mode), FUN=sum)
head(TROPHIC_otus)
names(TROPHIC_otus)[1]<- paste("TROPHIC")
rownames(TROPHIC_otus) <- TROPHIC_otus[,1]
head(TROPHIC_otus)
rownames(TROPHIC_otus)


TROPHIC_modify <- c("Unassigned","Pathogen-Saprotroph-Symbiotroph", "Pathotroph","Pathotroph-Saprotroph","Pathotroph-Saprotroph-Symbiotroph", 
                    "Pathotroph-Symbiotroph", "Saprotroph","Saprotroph-Pathotroph-Symbiotroph", "Saprotroph-Symbiotroph","Symbiotroph")
TROPHIC_otus[,1] <- TROPHIC_modify
str(TROPHIC_otus)

###
TROPHIC_otus2 <- aggregate(TROPHIC_otus[,2:ncol(TROPHIC_otus)], by=list(TROPHIC_otus$TROPHIC), FUN=sum)
head(TROPHIC_otus2)
names(TROPHIC_otus2)[1]<- paste("TROPHIC")
rownames(TROPHIC_otus2) <- TROPHIC_otus2[,1]
str(TROPHIC_otus2)

########construct data frame for calculate, use the melt function
TROPHIC_otus_sample_only <- TROPHIC_otus2[,2:ncol(TROPHIC_otus2)]
TROPHIC_otus_sample_only_order <-TROPHIC_otus_sample_only

###Percent of unidentified based on OTUs
TROPHIC_otus_type_prop<-TROPHIC_otus_sample_only_order
str(TROPHIC_otus_type_prop)
TROPHIC_otus_percent<-TROPHIC_otus_type_prop*100/sum(TROPHIC_otus_type_prop)
TROPHIC_otus_per<-cbind(TROPHIC_otus2,TROPHIC_otus_percent)
TROPHIC_otus_per
sum(TROPHIC_otus_per$TROPHIC_otus_percent)


####Proportion of unidentified based on otus presence absence across all samples
TROPHIC_pres_abs<-fungi_97pick_down_strict_initial_11[,1:sample_l]
TROPHIC_pres_abs[TROPHIC_pres_abs>0]<-1

TROPHIC_pres_abs[1:3]
TROPHIC<-aggregate(TROPHIC_pres_abs, by=list(fungi_97pick_down_strict_initial_11$Trophic_Mode), FUN=sum)
head(TROPHIC)
names(TROPHIC)[1]<- paste("TROPHIC")
rownames(TROPHIC) <- TROPHIC[,1]
head(TROPHIC)
rownames(TROPHIC)

TROPHIC_modify <- c("Unassigned","Pathogen-Saprotroph-Symbiotroph", "Pathotroph","Pathotroph-Saprotroph","Pathotroph-Saprotroph-Symbiotroph", 
                    "Pathotroph-Symbiotroph", "Saprotroph","Saprotroph-Pathotroph-Symbiotroph", "Saprotroph-Symbiotroph","Symbiotroph" )

TROPHIC[,1] <- TROPHIC_modify
str(TROPHIC)

###
TROPHIC2 <- aggregate(TROPHIC[,2:ncol(TROPHIC)], by=list(TROPHIC$TROPHIC), FUN=sum)
head(TROPHIC2)
names(TROPHIC2)[1]<- paste("TROPHIC")
rownames(TROPHIC2) <- TROPHIC2[,1]
str(TROPHIC2)

########construct data frame for calculate, use the melt function
TROPHIC_sample_only <- TROPHIC2[,2:ncol(TROPHIC2)]
TROPHIC_sample_only_order <- TROPHIC_sample_only[,order(colnames(TROPHIC_sample_only))]

###Percent of unidentified based on OTUs
TROPHIC_type_prop<-TROPHIC_sample_only_order
TROPHIC_type_prop$total<-rowSums(TROPHIC_type_prop)
TROPHIC_type_prop[1:10, 1:3]
TROPHIC_type_prop$percent<-(TROPHIC_type_prop$total/sum(TROPHIC_type_prop$total) * 100)
TROPHIC_type_prop[1:10, 418:420]
sum(TROPHIC_type_prop$percent)



#####
Trophic_agg <- aggregate(fungi_97pick_down_strict_initial_11[,1:sample_l], by=list(fungi_97pick_down_strict_initial_11$Trophic_Mode), FUN=sum)

head(Trophic_agg)
dim(Trophic_agg)
Trophic_agg[1:8,1:6]
names(Trophic_agg)[1]<- paste("Trophic")
rownames(Trophic_agg) <- Trophic_agg[,1]
head(Trophic_agg)
str(Trophic_agg)

####modify trophic 
Trophic_modify <- c("Unassigned","Pathogen-Saprotroph-Symbiotroph", "Pathotroph","Pathotroph-Saprotroph","Pathotroph-Saprotroph-Symbiotroph", 
                    "Pathotroph-Symbiotroph", "Saprotroph","Saprotroph-Pathotroph-Symbiotroph", "Saprotroph-Symbiotroph","Symbiotroph" )
Trophic_agg[,1] <- Trophic_modify
Trophic_agg[1:10,1:6]
str(Trophic_agg)
dim(Trophic_agg)
head(Trophic_agg,2)

Trophic_agg2 <- aggregate(Trophic_agg[,2:ncol(Trophic_agg)], by=list(Trophic_agg$Trophic), FUN=sum)
head(Trophic_agg2)
Trophic_agg2[1:10,1:6]
names(Trophic_agg2)[1]<- paste("Trophic")
rownames(Trophic_agg2) <- Trophic_agg2[,1]
dim(Trophic_agg2)

########construct data frame for calculate, use the melt function
Trophic_sample_only <- Trophic_agg2[,2:ncol(Trophic_agg2)]
head(Trophic_sample_only)
dim(Trophic_sample_only)

Trophic_sample_only_order <- Trophic_sample_only[,order(colnames(Trophic_sample_only))]
dim(Trophic_sample_only_order)
head(Trophic_sample_only_order)
Trophic_sample_only_order[1:6,1:6]

###Percent of unidentified based on abundance of OTUs
Trophic_sample_only_order[1:7,1:3]

Trophic_prop<-Trophic_sample_only_order
Trophic_prop$total<-rowSums(Trophic_prop)
Trophic_prop[1:10, 418:419]
Trophic_prop$percent<-(Trophic_prop$total/sum(Trophic_prop$total) * 100)
sum(Trophic_prop$percent)

Trophic_prop[1:10, 418:420]
#########combine with map file
dim(Down_strict_initial_map)
head(Down_strict_initial_map)
head(t(Trophic_sample_only_order))
match(rownames(t(Trophic_sample_only_order)),Down_strict_initial_map$SampleID)
###Matching is correct

Trophic_combine <- cbind(t(Trophic_sample_only_order),Down_strict_initial_map)
dim(Trophic_combine)

melt_id1 <- nrow(Trophic_sample_only_order)+1
melt_id1
melt_id2 <- ncol(Trophic_combine)
melt_id2

library(reshape2)

##Use melt to actually put the trait under one column and add set of rows with the ith information of ith trait 
## Basically melt is collapsing several columnS to two (one for the columns and one for their values)
#library(reshape2)
Trophic_melt <- reshape2::melt(Trophic_combine,id=melt_id1:melt_id2,variable.name="Trophic",value.name="Count")
str(Trophic_melt)

Trophic_time <- with(Trophic_melt,aggregate(Count,by=list(Time=Time,Trophic=Trophic,Species=Species,Forest=Forest),function(x) mean(x,na.rm=T)))
Try_trophic_total<-Trophic_time
head(Try_trophic_total)
sum(Try_trophic_total$x)



l_try <- length(levels(Try_trophic_total$Trophic))*2

Initial_strict_down_try <- do.call("rbind", replicate(3, Try_trophic_total[1:l_try,], simplify = FALSE))
Initial_strict_down_try$Forest <- c(paste(replicate(l_try,"Mature_forest")),paste(replicate(l_try,"Open_land")),paste(replicate(l_try,"Regenerating_forest")))


strict_down_try_plot <- rbind(Initial_strict_down_try,Try_trophic_total[(l_try+1):nrow(Try_trophic_total),])
strict_down_try_plot[] <- lapply(strict_down_try_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_try_plot$Forest
strict_down_try_plot$Forest <- factor(strict_down_try_plot$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))
strict_down_try_plot$Trophic
strict_down_try_plot$Trophic <- factor(strict_down_try_plot$Trophic,levels=c("Unassigned","Pathogen-Saprotroph-Symbiotroph", "Pathotroph","Pathotroph-Saprotroph","Pathotroph-Saprotroph-Symbiotroph", 
                                                                             "Pathotroph-Symbiotroph", "Saprotroph","Saprotroph-Pathotroph-Symbiotroph", "Saprotroph-Symbiotroph","Symbiotroph"))

strict_down_try_plot$x_percent<-(strict_down_try_plot$x)*100/sum(strict_down_try_plot$x)
sum(strict_down_try_plot$x_percent)

levels(strict_down_try_plot$Species)
strict_down_try_plot$Species<-as.factor(strict_down_try_plot$Species)
levels(strict_down_try_plot$Species)
levels(strict_down_try_plot$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_try_plot$Species)

levels(strict_down_try_plot$Forest)
levels(strict_down_try_plot$Forest)<-c("Open land", "Regenerating forest", "Mature forest")

ggplot(strict_down_try_plot, aes(x = Time, y = x_percent, fill = Trophic)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  ylab("Relative abundance (%)")+
  xlab("Exposure time (Months)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))



Trophic_perct<- ddply(Trophic_time, .(Forest,Time,Species), mutate, Trophic_pct = x / sum(x) * 100)
levels(Trophic_perct$Time)
l_initial <- length(levels(Trophic_perct$Trophic))*2

Initial_strict_down_rep <- do.call("rbind", replicate(3, Trophic_perct[1:l_initial,], simplify = FALSE))
Initial_strict_down_rep$Forest <- c(paste(replicate(l_initial,"Mature_forest")),paste(replicate(l_initial,"Open_land")),paste(replicate(l_initial,"Regenerating_forest")))

strict_down_perct_plot <- rbind(Initial_strict_down_rep,Trophic_perct[(l_initial+1):nrow(Trophic_perct),])
strict_down_perct_plot[] <- lapply(strict_down_perct_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_perct_plot$Forest <- factor(strict_down_perct_plot$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))
strict_down_perct_plot$Trophic <- factor(strict_down_perct_plot$Trophic,levels=c("Unassigned","Pathogen-Saprotroph-Symbiotroph", "Pathotroph","Pathotroph-Saprotroph","Pathotroph-Saprotroph-Symbiotroph", 
                                                                                 "Pathotroph-Symbiotroph", "Saprotroph","Saprotroph-Pathotroph-Symbiotroph", "Saprotroph-Symbiotroph","Symbiotroph"))

trophic_all_unas_100<-ggplot(strict_down_perct_plot, aes(x = Time, y = Trophic_pct, fill = Trophic)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +ylab("Percentage")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))#+

trophic_all_unas_100



###
#################exclude unassigned from bar graph
head(Trophic_agg2)


Trophic_agg3 <- Trophic_agg2[-10,]
########construct data frame for calculate, use the melt function
Trophic_sample_only2 <- Trophic_agg3[,2:ncol(Trophic_agg3)]
Trophic_sample_only_order2 <- Trophic_sample_only2[,order(colnames(Trophic_sample_only2))]


head(t(Trophic_sample_only_order2))
match(rownames(t(Trophic_sample_only_order2)),Down_strict_initial_map$SampleID)
#########combine with map file

match(rownames(t(Trophic_sample_only_order2)),Down_strict_initial_map$SampleID)

Trophic_combine2 <- cbind(t(Trophic_sample_only_order2),Down_strict_initial_map)
dim(Trophic_combine2)

melt_id1 <- nrow(Trophic_sample_only_order2)+1
melt_id2 <- ncol(Trophic_combine2)
Trophic_melt2 <- reshape2::melt(Trophic_combine2,id=melt_id1:melt_id2,variable.name = "Trophic",value.name = "Count")

Trophic_time2 <- with(Trophic_melt2,aggregate(Count,by=list(Time=Time,Trophic=Trophic,Species=Species,Forest=Forest),function(x) mean(x,na.rm=T)))
str(Trophic_time2)

Try_trophic_total2<-Trophic_time2
head(Try_trophic_total2)
sum(Try_trophic_total2$x)


l_try <- length(levels(Try_trophic_total2$Trophic))*2
Initial_strict_down_try2 <- do.call("rbind", replicate(3, Try_trophic_total2[1:l_try,], simplify = FALSE))
Initial_strict_down_try2$Forest <- c(paste(replicate(l_try,"Mature_forest")),paste(replicate(l_try,"Open_land")),paste(replicate(l_try,"Regenerating_forest")))


strict_down_try_plot2 <- rbind(Initial_strict_down_try2,Try_trophic_total2[(l_try+1):nrow(Try_trophic_total2),])
strict_down_try_plot2[] <- lapply(strict_down_try_plot2, function(x) if(is.factor(x)) factor(x) else x)
strict_down_try_plot2$Forest <- factor(strict_down_try_plot2$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))
strict_down_try_plot2$Trophic
strict_down_try_plot2$Trophic <- factor(strict_down_try_plot2$Trophic,levels=c("Pathogen-Saprotroph-Symbiotroph", "Pathotroph","Pathotroph-Saprotroph","Pathotroph-Saprotroph-Symbiotroph", 
                                                                               "Pathotroph-Symbiotroph", "Saprotroph","Saprotroph-Pathotroph-Symbiotroph", "Saprotroph-Symbiotroph","Symbiotroph"))

strict_down_try_plot2$x_percent<-(strict_down_try_plot2$x)*100/sum(strict_down_try_plot2$x)
sum(strict_down_try_plot2$x_percent)

levels(strict_down_try_plot2$Species)
strict_down_try_plot2$Species<-as.factor(strict_down_try_plot2$Species)
levels(strict_down_try_plot2$Species)
levels(strict_down_try_plot2$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_try_plot2$Species)
levels(strict_down_try_plot2$Forest)
levels(strict_down_try_plot2$Forest)<-c("Open land", "Regenerating forest", "Mature forest")


trophic_all_not100<-ggplot(strict_down_try_plot2, aes(x = Time, y = x_percent, fill = Trophic)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  ylab("Relative abundance (%)")+
  xlab("Exposure time (months)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))
trophic_all_not100


Trophic_perct2<- ddply(Trophic_time2, .(Forest,Time,Species), mutate, Trophic_pct2 = x / sum(x) * 100)
levels(Trophic_perct2$Time)
l_initial <- length(levels(Trophic_perct2$Trophic))*2

Initial_strict_down_rep2 <- do.call("rbind", replicate(3, Trophic_perct2[1:l_initial,], simplify = FALSE))
Initial_strict_down_rep2$Forest <- c(paste(replicate(l_initial,"Mature_forest")),paste(replicate(l_initial,"Open_land")),paste(replicate(l_initial,"Regenerating_forest")))

strict_down_perct_plot2 <- rbind(Initial_strict_down_rep2,Trophic_perct2[(l_initial+1):nrow(Trophic_perct2),])
strict_down_perct_plot2[] <- lapply(strict_down_perct_plot2, function(x) if(is.factor(x)) factor(x) else x)
strict_down_perct_plot2$Forest
strict_down_perct_plot2$Forest <- factor(strict_down_perct_plot2$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))
strict_down_perct_plot2$Trophic <- factor(strict_down_perct_plot2$Trophic,levels=c("Pathogen-Saprotroph-Symbiotroph", "Pathotroph","Pathotroph-Saprotroph","Pathotroph-Saprotroph-Symbiotroph", 
                                                                                   "Pathotroph-Symbiotroph", "Saprotroph","Saprotroph-Pathotroph-Symbiotroph", "Saprotroph-Symbiotroph","Symbiotroph"))
levels(strict_down_perct_plot2$Species)
strict_down_perct_plot2$Species<-as.factor(strict_down_perct_plot2$Species)
levels(strict_down_perct_plot2$Species)
levels(strict_down_perct_plot2$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_perct_plot2$Species)

levels(strict_down_perct_plot2$Forest)
levels(strict_down_perct_plot2$Forest)<-c("Open land", "Regenerating forest", "Mature forest")



trophic_all_100<-ggplot(strict_down_perct_plot2, aes(x = Time, y = Trophic_pct2, fill = Trophic)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3)+
  ylab("Abundance (%)")+
  xlab("Exposure time (months)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))
trophic_all_100


####Exclude non specific trophic such as ("Pathogen-Saprothroph-Symbiotroph", Pathotroph-Saprothroph, 
######"Pathotroph-Saprothroph-Symbiothroph", "Pathrothroph-Symbiothroph", "Saprothroph-Pathothroph-Symbiothroph","
#######Saprothroph-Symbiothroph")

#################exclude the above from the bar graph
head(Trophic_agg2)
Trophic_agg2$Trophic
dim(Trophic_agg2)

###Trophic_agg2 for 

Trophic_agg4 <- Trophic_agg2[-c(1,3,4,5,7,8),]
########construct data frame for calculate, use the melt function
Trophic_sample_only4 <- Trophic_agg4[,2:ncol(Trophic_agg4)]
Trophic_sample_only_order4 <- Trophic_sample_only4[,order(colnames(Trophic_sample_only4))]

#########combine with map file
head(t(Trophic_sample_only_order4))
match(rownames(t(Trophic_sample_only_order4)),Down_strict_initial_map$SampleID)


Trophic_combine4 <- cbind(t(Trophic_sample_only_order4),Down_strict_initial_map)
dim(Trophic_combine4)

melt_id1 <- nrow(Trophic_sample_only_order4)+1
melt_id2 <- ncol(Trophic_combine4)
Trophic_melt4 <- reshape2::melt(Trophic_combine4,id=melt_id1:melt_id2,variable.name = "Trophic",value.name = "Count")

Trophic_time4 <- with(Trophic_melt4,aggregate(Count,by=list(Time=Time,Trophic=Trophic,Species=Species,Forest=Forest),function(x) mean(x,na.rm=T)))

Try_trophic_total4<-Trophic_time4
head(Try_trophic_total4)
sum(Try_trophic_total4$x)


l_try <- length(levels(Try_trophic_total4$Trophic))*2
Initial_strict_down_try4 <- do.call("rbind", replicate(3, Try_trophic_total4[1:l_try,], simplify = FALSE))
Initial_strict_down_try4$Forest <- c(paste(replicate(l_try,"Mature_forest")),paste(replicate(l_try,"Open_land")),paste(replicate(l_try,"Regenerating_forest")))


strict_down_try_plot4 <- rbind(Initial_strict_down_try4,Try_trophic_total4[(l_try+1):nrow(Try_trophic_total4),])
strict_down_try_plot4[] <- lapply(strict_down_try_plot4, function(x) if(is.factor(x)) factor(x) else x)
strict_down_try_plot4$Forest
strict_down_try_plot4$Forest <- factor(strict_down_try_plot4$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))
strict_down_try_plot4$Trophic
strict_down_try_plot4$Trophic <- factor(strict_down_try_plot4$Trophic,levels=c("Pathotroph","Saprotroph","Symbiotroph", "Unassigned"))

strict_down_try_plot4$x_percent<-(strict_down_try_plot4$x)*100/sum(strict_down_try_plot4$x)
dim(strict_down_try_plot4)


levels(strict_down_try_plot4$Species)
strict_down_try_plot4$Species<-as.factor(strict_down_try_plot4$Species)
levels(strict_down_try_plot4$Species)
levels(strict_down_try_plot4$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_try_plot4$Species)

levels(strict_down_try_plot4$Forest)
levels(strict_down_try_plot4$Forest)<-c("Open land", "Regenerating forest", "Mature forest")
levels(strict_down_try_plot4$Forest)

trophic_Pat_Sap_Sym_not100<-ggplot(strict_down_try_plot4, aes(x = Time, y = x_percent, fill = Trophic)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  ylab("Relative abundance (%)")+
  xlab("Exposure time (months)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))
trophic_Pat_Sap_Sym_not100




Trophic_perct4<- ddply(Trophic_time4, .(Forest,Time,Species), mutate, Trophic_pct4 = x / sum(x) * 100)
levels(Trophic_perct4$Time)

l_initial <- length(levels(Trophic_perct4$Trophic))*2

Initial_strict_down_rep4 <- do.call("rbind", replicate(3, Trophic_perct4[1:l_initial,], simplify = FALSE))
Initial_strict_down_rep4$Forest <- c(paste(replicate(l_initial,"Mature_forest")),paste(replicate(l_initial,"Open_land")),paste(replicate(l_initial,"Regenerating_forest")))

strict_down_perct_plot4 <- rbind(Initial_strict_down_rep4,Trophic_perct4[(l_initial+1):nrow(Trophic_perct4),])
strict_down_perct_plot4[] <- lapply(strict_down_perct_plot4, function(x) if(is.factor(x)) factor(x) else x)
strict_down_perct_plot4$Forest <- factor(strict_down_perct_plot4$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))
strict_down_perct_plot4$Trophic <- factor(strict_down_perct_plot4$Trophic,levels=c("Pathotroph","Saprotroph","Symbiotroph","Unassigned"))

levels(strict_down_perct_plot4$Species)
strict_down_perct_plot4$Species<-as.factor(strict_down_perct_plot4$Species)
levels(strict_down_perct_plot4$Species)
levels(strict_down_perct_plot4$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_perct_plot4$Species)
levels(strict_down_perct_plot4$Forest)
levels(strict_down_perct_plot4$Forest)<-c("Open land", "Regenerating forest", "Mature forest")
levels(strict_down_perct_plot4$Forest)


trophic_Pat_Sap_Sym_una_100<-ggplot(strict_down_perct_plot4, aes(x = Time, y = Trophic_pct4, fill = Trophic)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3)+
  ylab("Abundance (%)")+
  xlab("Exposure time (months)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))

trophic_Pat_Sap_Sym_una_100


###Trophic Exclude unassigned

####Exclude non specific trophic such as ("Pathogen-Saprothroph-Symbiotroph", Pathotroph-Saprothroph, 
######"Pathotroph-Saprothroph-Symbiothroph", "Pathrothroph-Symbiothroph", "Saprothroph-Pathothroph-Symbiothroph","
#######Saprothroph-Symbiothroph"), "Unassigned"

#################exclude the above from the bar graph
head(Trophic_agg2)
Trophic_agg2$Trophic

Trophic_agg5 <- Trophic_agg2[-c(1,3,4,5,7,8,10),]
########construct data frame for calculate, use the melt function
Trophic_sample_only5 <- Trophic_agg5[,2:ncol(Trophic_agg5)]
Trophic_sample_only_order5 <- Trophic_sample_only5[,order(colnames(Trophic_sample_only5))]

#########combine with map file
head(t(Trophic_sample_only_order5))
match(rownames(t(Trophic_sample_only_order5)),Down_strict_initial_map$SampleID)
##Matching is correct

Trophic_combine5 <- cbind(t(Trophic_sample_only_order5),Down_strict_initial_map)
dim(Trophic_combine5)

melt_id1 <- nrow(Trophic_sample_only_order5)+1
melt_id2 <- ncol(Trophic_combine5)
Trophic_melt5 <- reshape2::melt(Trophic_combine5,id=melt_id1:melt_id2,variable.name = "Trophic",value.name = "Count")

Trophic_time5 <- with(Trophic_melt5,aggregate(Count,by=list(Time=Time,Trophic=Trophic,Species=Species,Forest=Forest),function(x) mean(x,na.rm=T)))

Try_trophic_total5<-Trophic_time5
head(Try_trophic_total5)
sum(Try_trophic_total5$x)



l_try <- length(levels(Try_trophic_total5$Trophic))*2
Initial_strict_down_try5 <- do.call("rbind", replicate(3, Try_trophic_total5[1:l_try,], simplify = FALSE))
Initial_strict_down_try5$Forest <- c(paste(replicate(l_try,"Mature_forest")),paste(replicate(l_try,"Open_land")),paste(replicate(l_try,"Regenerating_forest")))


strict_down_try_plot5 <- rbind(Initial_strict_down_try5,Try_trophic_total5[(l_try+1):nrow(Try_trophic_total5),])
strict_down_try_plot5[] <- lapply(strict_down_try_plot5, function(x) if(is.factor(x)) factor(x) else x)
strict_down_try_plot5$Forest
strict_down_try_plot5$Forest <- factor(strict_down_try_plot5$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))
strict_down_try_plot5$Trophic
strict_down_try_plot5$Trophic <- factor(strict_down_try_plot5$Trophic,levels=c("Pathotroph","Saprotroph","Symbiotroph"))
strict_down_try_plot5$Trophic <- factor(strict_down_try_plot5$Trophic,levels=c("Saprotroph","Symbiotroph","Pathotroph"))

strict_down_try_plot5$x_percent<-(strict_down_try_plot5$x)*100/sum(strict_down_try_plot5$x)
sum(strict_down_try_plot5$x_percent)
dim(strict_down_try_plot5)
levels(strict_down_try_plot5$Species)
strict_down_try_plot5$Species<-as.factor(strict_down_try_plot5$Species)
levels(strict_down_try_plot5$Species)
levels(strict_down_try_plot5$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_try_plot5$Forest)
strict_down_try_plot5$Forest<-as.factor(strict_down_try_plot5$Forest)
levels(strict_down_try_plot5$Forest)
levels(strict_down_try_plot5$Forest)<-c("Open land", "Regenerating forest", "Mature forest")
levels(strict_down_try_plot5$Forest)

trophic_Pat_Sap_Sym_not100<-ggplot(strict_down_try_plot5, aes(x = Time, y = x_percent, fill = Trophic)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  ylab("Relative abundance (%)")+
  xlab("Exposure time (months)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 20),axis.text.y=element_text(colour = 'black', size = 20))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=20),axis.text.x = element_text(angle=0,colour = 'black', size=20,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=20), legend.title=element_text(size=20),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))+
  scale_fill_brewer(palette = "Paired")### With this you can change to colors needed
trophic_Pat_Sap_Sym_not100

str(Trophic_time5)

Trophic_perct5<- ddply(Trophic_time5, .(Forest,Time,Species), mutate, Trophic_pct5 = x / sum(x) * 100)
levels(Trophic_perct5$Time)
l_initial <- length(levels(Trophic_perct5$Trophic))*2

Initial_strict_down_rep5 <- do.call("rbind", replicate(3, Trophic_perct5[1:l_initial,], simplify = FALSE))
Initial_strict_down_rep5$Forest <- c(paste(replicate(l_initial,"Mature_forest")),paste(replicate(l_initial,"Open_land")),paste(replicate(l_initial,"Regenerating_forest")))

strict_down_perct_plot5 <- rbind(Initial_strict_down_rep5,Trophic_perct5[(l_initial+1):nrow(Trophic_perct5),])
strict_down_perct_plot5[] <- lapply(strict_down_perct_plot5, function(x) if(is.factor(x)) factor(x) else x)
strict_down_perct_plot5$Forest <- factor(strict_down_perct_plot5$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))
strict_down_perct_plot5$Trophic <- factor(strict_down_perct_plot5$Trophic,levels=c("Pathotroph","Saprotroph","Symbiotroph"))
strict_down_perct_plot5$Trophic <- factor(strict_down_perct_plot5$Trophic,levels=c("Saprotroph","Symbiotroph","Pathotroph"))


levels(strict_down_perct_plot5$Species)
strict_down_perct_plot5$Species<-as.factor(strict_down_perct_plot5$Species)
levels(strict_down_perct_plot5$Species)
levels(strict_down_perct_plot5$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_perct_plot5$Forest)
levels(strict_down_perct_plot5$Forest)<-c("Open land", "Regenerating forest", "Mature forest")
levels(strict_down_perct_plot5$Forest)



trophic_Pat_Sap_Sym_100<-ggplot(strict_down_perct_plot5, aes(x = Time, y = Trophic_pct5, fill = Trophic)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3)+
  ylab("Abundance (%)")+
  xlab("Exposure time (months)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 20),axis.text.y=element_text(colour = 'black', size = 20))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=20),axis.text.x = element_text(angle=0,colour = 'black', size=20,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=20), legend.title=element_text(size=20),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))+
  scale_fill_brewer(palette = "Paired")### With this you can change to colors needed
trophic_Pat_Sap_Sym_100


###Computing the average of Pathotrophs, Saproptrophs  and Symbiothroph across time and
#####omputation of mean
str(strict_down_perct_plot5)

Trophic_means<-strict_down_perct_plot5%>%
  group_by(Trophic, Time)%>%
  summarise(mean_percent=mean(Trophic_pct5),
            sd_percent=sd(Trophic_pct5))
Trophic_means

####Saprotrooh peaked at 18 mo around 82.3%




# A tibble: 9 x 4
# Groups:   Trophic [3]
#Trophic     Time  mean_percent sd_percent
#<fct>       <chr>        <dbl>      <dbl>
#  1 Saprotroph  0mo          77.5       6.76 
#2 Saprotroph  18mo         82.3       2.58 
#3 Saprotroph  36mo         77.4      10.2  
#4 Symbiotroph 0mo           1.25      0.726
#5 Symbiotroph 18mo          6.12      3.44 
#6 Symbiotroph 36mo         14.2      11.0  
#7 Pathotroph  0mo          21.3       7.49 
#8 Pathotroph  18mo         11.6       3.13 
#9 Pathotroph  36mo          8.43      2.26 

###We can also go deeper by looking at specific Species at a time

Trophic_means2<-strict_down_perct_plot5%>%
  group_by(Species,Trophic, Time)%>%
  summarise(mean_percent=mean(Trophic_pct5),
            sd_percent=sd(Trophic_pct5))
Trophic_means2

####




#> Trophic_means2
# A tibble: 18 x 5
# Groups:   Species, Trophic [6]
#Species                 Trophic     Time  mean_percent sd_percent
#<fct>                   <fct>       <chr>        <dbl>      <dbl>
#1 Castanopsis mekongensis Saprotroph  0mo         83.7         0   
#2 Castanopsis mekongensis Saprotroph  18mo        83.4         2.38
#3 Castanopsis mekongensis Saprotroph  36mo        84.5         1.06
#4 Castanopsis mekongensis Symbiotroph 0mo          1.91        0   
#5 Castanopsis mekongensis Symbiotroph 18mo         5.49        2.42
#6 Castanopsis mekongensis Symbiotroph 36mo         6.34        3.50
#7 Castanopsis mekongensis Pathotroph  0mo         14.4         0   
#8 Castanopsis mekongensis Pathotroph  18mo        11.2         4.46
#9 Castanopsis mekongensis Pathotroph  36mo         9.21        3.06
#10 Litsea cubeba           Saprotroph  0mo         71.3         0   
#11 Litsea cubeba           Saprotroph  18mo        81.3         2.77
#12 Litsea cubeba           Saprotroph  36mo        70.3        10.5 
#13 Litsea cubeba           Symbiotroph 0mo          0.586       0   
#14 Litsea cubeba           Symbiotroph 18mo         6.75        4.74
#15 Litsea cubeba           Symbiotroph 36mo        22.0        10.3 
#16 Litsea cubeba           Pathotroph  0mo         28.1         0   
#17 Litsea cubeba           Pathotroph  18mo        12.0         2.03
#18 Litsea cubeba           Pathotroph  36mo         7.66        1.27


#####Figure 2
####Put al three panels together

###Put figure in panels
##load libraries
library(ggplot2)
library(gridExtra)
library(reshape)
library(grid)



tiff(filename="Figure 3 Phylum Trophic rot.tiff", res = 400, width=11000, height=7500, compression = "lzw")
grid.arrange(arrangeGrob(phylum_Asc_Bas_Oth_not100 + theme(legend.position="right",
                                                                  plot.margin = unit( c(1,0.5,0,0) , units = "lines")),
                         trophic_Pat_Sap_Sym_not100 + theme(legend.position="right",
                                                                          plot.margin = unit( c(1,0.5,0,0) , units = "lines")),
                         rot_bro_sof_whi_not100 + theme(legend.position="right",
                                                                           plot.margin = unit( c(1,0.5,0,0) , units = "lines")),
                         nrow=2, ncol=2))

dev.off()

###Note that the number otus on each bar was added in powerpoint

###Figure S3  Abundance up to 100

tiff(filename="Figure S3 Phylum Trophic rot rot with unassigned.tiff", res = 400, width=11000, height=7500, compression = "lzw")
grid.arrange(arrangeGrob(phylum_Asc_Bas_Oth_100 + theme(legend.position="right",
                                                           plot.margin = unit( c(1,0.5,0,0) , units = "lines")),
                         trophic_Pat_Sap_Sym_100 + theme(legend.position="right",
                                                            plot.margin = unit( c(1,0.5,0,0) , units = "lines")),
                         rot_bro_sof_whi_100 + theme(legend.position="right",
                                                        plot.margin = unit( c(1,0.5,0,0) , units = "lines")),
                         rot_all_una_not100 + theme(legend.position="right",
                                                        plot.margin = unit( c(1,0.5,0,0) , units = "lines")),
                         nrow=2, ncol=2))

dev.off()




### Class level

###Proportion based on just total OTUs list

###Proportion of unidentified based on OTUs presence
CLASS_otus<-fungi_97pick_down_strict_initial_11[,1]
###Make OTUs abundance more than 1 to be one just as present
CLASS_otus[CLASS_otus>0]<-1
###Make OTUs abundance 0 be one just as present in the list of otus
CLASS_otus[CLASS_otus==0]<-1
CLASS_otus[1:10]
####Everything should be one

CLASS_otus<-aggregate(CLASS_otus, by=list(fungi_97pick_down_strict_initial_11$Class), FUN=sum)
head(CLASS_otus)
names(CLASS_otus)[1]<- paste("CLASS")
rownames(CLASS_otus) <- CLASS_otus[,1]
head(CLASS_otus)
rownames(CLASS_otus)


###
CLASS_otus2 <- aggregate(CLASS_otus[,2:ncol(CLASS_otus)], by=list(CLASS_otus$CLASS), FUN=sum)
head(CLASS_otus2)
names(CLASS_otus2)[1]<- paste("CLASS")
rownames(CLASS_otus2) <- CLASS_otus2[,1]
str(CLASS_otus2)

########construct data frame for calculate, use the melt function
CLASS_otus_sample_only <- CLASS_otus2[,2:ncol(CLASS_otus2)]
CLASS_otus_sample_only_order <-CLASS_otus_sample_only

###Percent of unidentified based on OTUs
CLASS_otus_type_prop<-CLASS_otus_sample_only_order
str(CLASS_otus_type_prop)
CLASS_otus_percent<-CLASS_otus_type_prop*100/sum(CLASS_otus_type_prop)
CLASS_otus_per<-cbind(CLASS_otus2,CLASS_otus_percent)
CLASS_otus_per
sum(CLASS_otus_per$CLASS_otus_percent)


####Proportion of unidentified based on presence/absence of OTUs
CLASS_pres_abs<-fungi_97pick_down_strict_initial_11[,1:sample_l]
CLASS_pres_abs[CLASS_pres_abs>0]<-1

CLASS_pres_abs[1:3]
CLASS<-aggregate(CLASS_pres_abs, by=list(fungi_97pick_down_strict_initial_11$Class), FUN=sum)
head(CLASS)
names(CLASS)[1]<- paste("CLASS")
rownames(CLASS) <- CLASS[,1]
head(CLASS)
rownames(CLASS)


str(CLASS)

###
CLASS2 <- aggregate(CLASS[,2:ncol(CLASS)], by=list(CLASS$CLASS), FUN=sum)
head(CLASS2)
names(CLASS2)[1]<- paste("CLASS")
rownames(CLASS2) <- CLASS2[,1]
str(CLASS2)

########construct data frame for calculate, use the melt function
CLASS_sample_only <- CLASS2[,2:ncol(CLASS2)]
CLASS_sample_only_order <- CLASS_sample_only[,order(colnames(CLASS_sample_only))]

###Percent of unidentified based on OTUs
CLASS_type_prop<-CLASS_sample_only_order
CLASS_type_prop$total<-rowSums(CLASS_type_prop)
CLASS_type_prop[1:37, 1:3]
CLASS_type_prop$percent<-(CLASS_type_prop$total/sum(CLASS_type_prop$total) * 100)
CLASS_type_prop[1:37, 418:420]
sum(CLASS_type_prop$percent)


###############collape to Class level, aggregate works same as the collapse function in qiime

Class_agg <- aggregate(fungi_97pick_down_strict_initial_11[,1:sample_l], by=list(fungi_97pick_down_strict_initial_11$Class), FUN=sum)
head(Class_agg)
dim(Class_agg)
names(Class_agg)[1]<- paste("Class")
rownames(Class_agg) <- Class_agg[,1]
head(Class_agg)
rownames(Class_agg)

str(Class_agg)

Class_agg2 <- aggregate(Class_agg[,2:ncol(Class_agg)], by=list(Class_agg$Class), FUN=sum)
head(Class_agg2)
names(Class_agg2)[1]<- paste("Order")
rownames(Class_agg2) <- Class_agg2[,1]

########construct data frame for calculate, use the melt function
Class_sample_only <- Class_agg2[,2:ncol(Class_agg2)]
Class_sample_only_order <- Class_sample_only[,order(colnames(Class_sample_only))]


###Percent of unidentified based on OTUs abundance
Class_sample_only_order[1:37,1:3]

Class_prop<-Class_sample_only_order
Class_prop$total<-rowSums(Class_prop)
Class_prop[1:37, 418:419]
Class_prop$percent<-(Class_prop$total/sum(Class_prop$total) * 100)
sum(Class_prop$percent)

Class_prop[1:37, 418:420]


#########combine with map file
head(t(Class_sample_only_order))
match(rownames(t(Class_sample_only_order)),Down_strict_initial_map$SampleID)
###Matching is correct
Class_combine <- cbind(t(Class_sample_only_order),Down_strict_initial_map)
dim(Class_combine)

#write.csv(Trait_combine, file="traits combine.csv")
melt_id1 <- nrow(Class_sample_only_order)+1
melt_id1
melt_id2 <- ncol(Class_combine)
melt_id2

library(reshape2)
##Use melt to actually put the trait under one column and add set of rows with the ith information of ith trait 
## Basically melt is collapsing several columns to two (one for the columns and one for their values)
Class_melt <- reshape2::melt(Class_combine,id=melt_id1:melt_id2,variable.name = "Class",value.name = "Count")
str(Class_melt)
summary(Class_melt)

Class_time <- with(Class_melt,aggregate(Count,by=list(Time=Time,Class=Class,Species=Species,Forest=Forest),function(x) mean(x,na.rm=T)))

Try_Class_total<-Class_time
head(Try_Class_total)
sum(Try_Class_total$x)




l_try <- length(levels(Try_Class_total$Class))*2

Initial_strict_down_try <- do.call("rbind", replicate(3, Try_Class_total[1:l_try,], simplify = FALSE))
Initial_strict_down_try$Forest <- c(paste(replicate(l_try,"Mature_forest")),paste(replicate(l_try,"Open_land")),paste(replicate(l_try,"Regenerating_forest")))


strict_down_try_plot <- rbind(Initial_strict_down_try,Try_Class_total[(l_try+1):nrow(Try_Class_total),])
strict_down_try_plot[] <- lapply(strict_down_try_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_try_plot$Forest <- factor(strict_down_try_plot$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))

strict_down_try_plot$Class
strict_down_try_plot$Class 

strict_down_try_plot$x_percent<-(strict_down_try_plot$x)*100/sum(strict_down_try_plot$x)
dim(strict_down_try_plot)
levels(strict_down_try_plot$Species)
levels(strict_down_try_plot$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_try_plot$Forest)
levels(strict_down_try_plot$Forest)<-c("Open land", "Regenerating forest", "Mature forest")



Class_perct<- ddply(Class_time, .(Forest,Time,Species), mutate, Class_pct = x / sum(x) * 100)
summary(Class_perct)

head(Class_perct)
str(Class_perct)


levels(Class_perct$Time)
l_initial <- length(levels(Class_perct$Class))*2

Initial_strict_down_rep <- do.call("rbind", replicate(3, Class_perct[1:l_initial,], simplify = FALSE))
Initial_strict_down_rep$Forest <- c(paste(replicate(l_initial,"Mature_forest")),paste(replicate(l_initial,"Open_land")),paste(replicate(l_initial,"Regenerating_forest")))

strict_down_perct_plot <- rbind(Initial_strict_down_rep,Class_perct[(l_initial+1):nrow(Class_perct),])
strict_down_perct_plot[] <- lapply(strict_down_perct_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_perct_plot$Forest <- factor(strict_down_perct_plot$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))

strict_down_perct_plot$Class ###


###Computing the average of classes across time and

str(strict_down_try_plot)
###Overall
Class_means_all<-strict_down_try_plot%>%
  group_by(Class)%>%
  summarise(mean_percent=mean(x_percent),
            sd_percent=sd(x_percent))
Class_means_all
sort(Class_means_all$mean_percent, decreasing =TRUE)

strict_down_try_plot
####Per species per forest, and time
Class_means1<-strict_down_try_plot%>%
  group_by(Forest,Species ,Class,Time)%>%
  summarise(mean_percent=mean(x_percent),
            sd_percent=sd(x_percent))
Class_means1
summary(Class_means1)

Class_means1[Class_means1$mean_percent>1,]

levels(Class_means1$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(Class_means1$Forest)<-c("Open land", "Regenerating forest", "Mature forest")


Class_means1_higher1percent<-ggplot(Class_means1[Class_means1$mean_percent>1,], aes(x = Time, y = mean_percent, fill = Class)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  ylab("Relative abundance (%)")+
  xlab("Exposure time")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))

Class_means1_higher1percent


#####Computation of mean for abundance
str(strict_down_perct_plot)
#library(dplyr)

Class_means_all<-strict_down_perct_plot%>%
  group_by(Class)%>%
  summarise(mean_percent=mean(Class_pct),
            sd_percent=sd(Class_pct))
Class_means_all
sort(Class_means_all$mean_percent, decreasing =TRUE)

Class_means2<-strict_down_perct_plot%>%
  group_by(Forest,Species,Class,Time)%>%
  summarise(mean_percent=mean(Class_pct),
            sd_percent=sd(Class_pct))
Class_means2
summary(Class_means2)
Class_means2[Class_means2$mean_percent>10,]
Class_means2[Class_means2$mean_percent>5,]

Class_means2$Species
Class_means2$Species<-as.factor(Class_means2$Species)
levels(Class_means2$Species)
Class_means2$Species<-factor(Class_means2$Species,levels=c("Castanopsis_mekongensis","Litsea_cubeba"))
levels(Class_means2$Species)
levels(Class_means2$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(Class_means2$Forest)
levels(Class_means2$Forest)<-c("Open land", "Regenerating forest", "Mature forest")
levels(Class_means2$Forest)

####

Class_means2_higher_10_percent<-ggplot(data=Class_means2[Class_means2$mean_percent>10,],aes(x = Time, y = mean_percent, fill = Class))+
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  xlab("Exposure duration (Months)")+ 
  ylab("Abundance (%)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))
Class_means2_higher_10_percent


Class_means2_higher_5_percent<-ggplot(data=Class_means2[Class_means2$mean_percent>5,],aes(x = Time, y = mean_percent, fill = Class))+
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  xlab("Exposure duration (Months)")+ 
  ylab("Abundance (%)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))+
  scale_fill_brewer(palette="Paired")

Class_means2_higher_5_percent



###Trying to fill the remaining of 100 % with Other
Class_means2[Class_means2$mean_percent>5,]

Class_means2$Class_new<-ifelse(Class_means2$mean_percent<=5, " All classes <= 5%", as.character(Class_means2$Class))
Class_means2$Class_new<-as.factor(Class_means2$Class_new)
levels(Class_means2$Class_new)
Class_means2$Class_new<-factor(Class_means2$Class_new, levels=c(" All classes <= 5%"," unidentified"," Agaricomycetes"," Dacrymycetes", " Dothideomycetes"," Eurotiomycetes",
                                                                " Leotiomycetes"," Pezizomycotina_cls_Incertae_sedis",
                                                                " Saccharomycetes"," Sordariomycetes"," Tremellomycetes"))
levels(Class_means2$Class_new)
head(Class_means2_higher_5_percent)

Class_means2_higher_5_with_other<-ggplot(data=Class_means2,aes(x = Time, y = mean_percent, fill = Class_new))+
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  xlab("Exposure time (Months)")+ 
  ylab("Abundance (%)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  #scale_fill_discrete(name="Most abundant classes > 5%")+
  labs(fill="Most abundant classes > 5%")+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 20),axis.text.y=element_text(colour = 'black', size = 20))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=20),axis.text.x = element_text(angle=0,colour = 'black', size=20,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=20), legend.title=element_text(size=20),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))+
  scale_fill_brewer(palette="Paired")

Class_means2_higher_5_with_other



####PPT figure for classes
Class_means2_higher_5_with_other<-ggplot(data=Class_means2,aes(x = Time, y = mean_percent, fill = Class_new))+
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  xlab("Exposure time (Months)")+ 
  ylab("Abundance (%)")+
  theme(strip.text.x = element_text(size = 18, colour = "black"))+
  #scale_fill_discrete(name="Most abundant classes > 5%")+
  labs(fill="Most abundant classes > 5%")+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 23),axis.text.y=element_text(colour = 'black', size = 23))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=23),axis.text.x = element_text(angle=0,colour = 'black', size=23,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=20), legend.title=element_text(size=20),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))+
  scale_fill_brewer(palette="Paired")

Class_means2_higher_5_with_other

###END

###Order


###Proportion based on just total OTUs list

###Proportion of unidentified based on OTUs presence
ORDER_otus<-fungi_97pick_down_strict_initial_11[,1]
###Make OTUs abundance more than 1 to be one just as present
ORDER_otus[ORDER_otus>0]<-1
###Make OTUs abundance 0 be one just as present in the list of otus
ORDER_otus[ORDER_otus==0]<-1
ORDER_otus[1:10]
####Everything should be one

ORDER_otus<-aggregate(ORDER_otus, by=list(fungi_97pick_down_strict_initial_11$Order), FUN=sum)
head(ORDER_otus)
names(ORDER_otus)[1]<- paste("ORDER")
rownames(ORDER_otus) <- ORDER_otus[,1]
head(ORDER_otus)
rownames(ORDER_otus)


###
ORDER_otus2 <- aggregate(ORDER_otus[,2:ncol(ORDER_otus)], by=list(ORDER_otus$ORDER), FUN=sum)
head(ORDER_otus2)
names(ORDER_otus2)[1]<- paste("ORDER")
rownames(ORDER_otus2) <- ORDER_otus2[,1]
str(ORDER_otus2)

########construct data frame for calculate, use the melt function
ORDER_otus_sample_only <- ORDER_otus2[,2:ncol(ORDER_otus2)]
ORDER_otus_sample_only_order <-ORDER_otus_sample_only

###Percent of unidentified based on OTUs
ORDER_otus_type_prop<-ORDER_otus_sample_only_order
str(ORDER_otus_type_prop)
ORDER_otus_percent<-ORDER_otus_type_prop*100/sum(ORDER_otus_type_prop)
ORDER_otus_per<-cbind(ORDER_otus2,ORDER_otus_percent)
ORDER_otus_per
sum(ORDER_otus_per$ORDER_otus_percent)


####Proportion of unidentified based on presence/absence of OTUs
ORDER_pres_abs<-fungi_97pick_down_strict_initial_11[,1:sample_l]
ORDER_pres_abs[ORDER_pres_abs>0]<-1

ORDER_pres_abs[1:3]
ORDER<-aggregate(ORDER_pres_abs, by=list(fungi_97pick_down_strict_initial_11$Order), FUN=sum)
head(ORDER)
names(ORDER)[1]<- paste("ORDER")
rownames(ORDER) <- ORDER[,1]
head(ORDER)
rownames(ORDER)


str(ORDER)

###
ORDER2 <- aggregate(ORDER[,2:ncol(ORDER)], by=list(ORDER$ORDER), FUN=sum)
head(ORDER2)
names(ORDER2)[1]<- paste("ORDER")
rownames(ORDER2) <- ORDER2[,1]
str(ORDER2)

########construct data frame for calculate, use the melt function
ORDER_sample_only <- ORDER2[,2:ncol(ORDER2)]
ORDER_sample_only_order <- ORDER_sample_only[,order(colnames(ORDER_sample_only))]

###Percent of unidentified based on OTUs
ORDER_type_prop<-ORDER_sample_only_order
ORDER_type_prop$total<-rowSums(ORDER_type_prop)
ORDER_type_prop[120:129, 1:3]
ORDER_type_prop$percent<-(ORDER_type_prop$total/sum(ORDER_type_prop$total) * 100)
ORDER_type_prop[120:129, 418:420]
sum(ORDER_type_prop$percent)


###
Order_agg <- aggregate(fungi_97pick_down_strict_initial_11[,1:sample_l], by=list(fungi_97pick_down_strict_initial_11$Order), FUN=sum)
head(Order_agg)
head(Order_agg)
names(Order_agg)[1]<- paste("Order")
rownames(Order_agg) <- Order_agg[,1]
head(Order_agg)
rownames(Order_agg)

str(Order_agg)

Order_agg2 <- aggregate(Order_agg[,2:ncol(Order_agg)], by=list(Order_agg$Order), FUN=sum)
head(Order_agg2)
names(Order_agg2)[1]<- paste("Order")
rownames(Order_agg2) <- Order_agg2[,1]

########construct data frame for calculate, use the melt function
Order_sample_only <- Order_agg2[,2:ncol(Order_agg2)]
Order_sample_only_order <- Order_sample_only[,order(colnames(Order_sample_only))]

###Percent of unidentified
Order_sample_only_order[120:129,1:3]

Order_prop<-Order_sample_only_order
Order_prop$total<-rowSums(Order_prop)
Order_prop[120:129, 418:419]
Order_prop$percent<-(Order_prop$total/sum(Order_prop$total) * 100)
sum(Order_prop$percent)

Order_prop[120:129, 418:420]

#########combine with map file
head(t(Order_sample_only_order))
match(rownames(t(Order_sample_only_order)),Down_strict_initial_map$SampleID)
##Matching is correct

Order_combine <- cbind(t(Order_sample_only_order),Down_strict_initial_map)
dim(Order_combine)

#write.csv(Trait_combine, file="traits combine.csv")
melt_id1 <- nrow(Order_sample_only_order)+1
melt_id1
melt_id2 <- ncol(Order_combine)
melt_id2

library(reshape2)
##Use melt to actually put the trait under one column and add set of rows with the ith information of ith trait 
## Basically melt is collapsing several columns to two (one for the columns and one for their values)
Order_melt <- reshape2::melt(Order_combine,id=melt_id1:melt_id2,variable.name = "Order",value.name = "Count")
str(Order_melt)


Order_time <- with(Order_melt,aggregate(Count,by=list(Time=Time,Order=Order,Species=Species,Forest=Forest),function(x) mean(x,na.rm=T)))

Try_Order_total<-Order_time

head(Try_Order_total)
sum(Try_Order_total$x)


l_try <- length(levels(Try_Order_total$Order))*2

Initial_strict_down_try <- do.call("rbind", replicate(3, Try_Order_total[1:l_try,], simplify = FALSE))
Initial_strict_down_try$Forest <- c(paste(replicate(l_try,"Mature_forest")),paste(replicate(l_try,"Open_land")),paste(replicate(l_try,"Regenerating_forest")))


strict_down_try_plot <- rbind(Initial_strict_down_try,Try_Order_total[(l_try+1):nrow(Try_Order_total),])
strict_down_try_plot[] <- lapply(strict_down_try_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_try_plot$Forest
strict_down_try_plot$Forest <- factor(strict_down_try_plot$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))

strict_down_try_plot$Order 
strict_down_try_plot$x_percent<-(strict_down_try_plot$x)*100/sum(strict_down_try_plot$x)

levels(strict_down_try_plot$Species)
strict_down_try_plot$Species<-as.factor(strict_down_try_plot$Species)
levels(strict_down_try_plot$Species)
levels(strict_down_try_plot$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_try_plot$Species)

levels(strict_down_try_plot$Forest)
levels(strict_down_try_plot$Forest)<-c("Open land", "Regenerating forest", "Mature forest")



Order_perct<- ddply(Order_time, .(Forest,Time,Species), mutate, Order_pct = x / sum(x) * 100)
summary(Order_perct)

head(Order_perct)
str(Order_perct)


levels(Order_perct$Time)
l_initial <- length(levels(Order_perct$Order))*2

Initial_strict_down_rep <- do.call("rbind", replicate(3, Order_perct[1:l_initial,], simplify = FALSE))
Initial_strict_down_rep$Forest <- c(paste(replicate(l_initial,"Mature_forest")),paste(replicate(l_initial,"Open_land")),paste(replicate(l_initial,"Regenerating_forest")))

strict_down_perct_plot <- rbind(Initial_strict_down_rep,Order_perct[(l_initial+1):nrow(Order_perct),])
strict_down_perct_plot[] <- lapply(strict_down_perct_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_perct_plot$Forest <- factor(strict_down_perct_plot$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))

strict_down_perct_plot$Order ###



###Computing the average of orders across time and

str(strict_down_try_plot)

Order_means1<-strict_down_try_plot%>%
  group_by(Forest,Species ,Order,Time)%>%
  summarise(mean_percent=mean(x_percent),
            sd_percent=sd(x_percent))
Order_means1
summary(Order_means1)
Order_means1[Order_means1$mean_percent>mean(Order_means1$mean_percent),]

Order_means1[Order_means1$mean_percent>1,]


Orders_means1_higher_0.5_percent<-ggplot(Order_means1[Order_means1$mean_percent>0.5,], aes(x = Time, y = mean_percent, fill = Order)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  ylab("Relative abundance (%)")+
  xlab("Exposure time")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))#+
 # scale_fill_brewer(palette="Paired")

Orders_means1_higher_0.5_percent

#####Computation of mean for relative abundance
str(strict_down_perct_plot)
#library(dplyr)
##Overall
Order_means_all<-strict_down_perct_plot%>%
  group_by(Order)%>%
  summarise(mean_percent=mean(Order_pct),
            sd_percent=sd(Order_pct))
Order_means_all
sort(Order_means_all$mean_percent, decreasing =TRUE)


Order_means2<-strict_down_perct_plot%>%
  group_by(Forest,Species ,Order,Time)%>%
  summarise(mean_percent=mean(Order_pct),
            sd_percent=sd(Order_pct))
Order_means2
summary(Order_means2)
Order_means2[Order_means2$mean_percent>10,]
Order_means2[Order_means2$mean_percent>5,]

####

Order_means2_higher_8_percent<-ggplot(data=Order_means2[Order_means2$mean_percent>8,],aes(x = Time, y = mean_percent, fill = Order))+
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  xlab("Exposure time (Months)")+ 
  ylab("Abundance (%)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))+
  scale_fill_brewer(palette="Paired")
Order_means2_higher_8_percent



###Trying to fill the remaining of 100 % with Other
Order_means2[Order_means2$mean_percent>8,]

Order_means2$Order_new<-ifelse(Order_means2$mean_percent<=8, " All orders <= 8 %", as.character(Order_means2$Order))
levels(Order_means2$Order_new)

levels(Order_means2$Species)
Order_means2$Species<-as.factor(Order_means2$Species)
levels(Order_means2$Species)

table(is.na(Order_means2$Order_new))
Order_means2$Order_new<-as.factor(Order_means2$Order_new)
levels(Order_means2$Order_new)

####
Order_means2$Order_new<-factor(Order_means2$Order_new, levels=c(" All orders <= 8 %"," unidentified"," Agaricales"," Capnodiales"," Chaetothyriales"," Coniochaetales",
                                                              " Eurotiales"," Helotiales" ," Hypocreales"," Pleosporales"," Polyporales",
                                                              " Saccharomycetales"," Xylariales" ))
levels(Order_means2$Order_new)
####

##### higher 5 %

Order_means2_higher_5_percent<-ggplot(data=Order_means2[Order_means2$mean_percent>5,],aes(x = Time, y = mean_percent, fill = Order))+
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  xlab("Exposure time (Months)")+ 
  ylab("Abundance (%)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))
Order_means2_higher_5_percent



###Trying to fill the remaining of 100 % with Other
Order_means2[Order_means2$mean_percent>5,]

Order_means2$Order_new<-ifelse(Order_means2$mean_percent<=5, " All orders <= 5 %", as.character(Order_means2$Order))
levels(Order_means2$Order_new)
table(is.na(Order_means2$Order_new))
Order_means2$Order_new<-as.factor(Order_means2$Order_new)
levels(Order_means2$Order_new)


Order_means2$Order_new<-factor(Order_means2$Order_new,levels=c(" All orders <= 5 %"," unidentified"," Agaricales"," Botryosphaeriales"," Capnodiales", " Chaetothyriales"," Coniochaetales"," Dacrymycetales",                      
                                                             " Eurotiales"," Helotiales"," Hypocreales"," Pezizomycotina_ord_Incertae_sedis"," Pleosporales"," Polyporales"," Saccharomycetales"," Sordariales", 
                                                             " Sordariomycetidae_ord_Incertae_sedis", " Xylariales"))                          



levels(Order_means2$Order_new)

levels(Order_means2$Species)
Order_means2$Species<-as.factor(Order_means2$Species)
levels(Order_means2$Species)
#Order_means2$Species<-factor(Order_means2$Species,levels=c("Castanopsis_mekongensis","Litsea_cubeba"))
levels(Order_means2$Species)
levels(Order_means2$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(Order_means2$Forest)<-c("Open land", "Regenerating forest", "Mature forest")


Order_means2_higher_5_percent_with_other<-ggplot(data=Order_means2,aes(x = Time, y = mean_percent, fill = Order_new))+
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  xlab("Exposure time (Months)")+ ylab("Abundance (%)")+
  scale_fill_discrete(name="Most abundant orders > 5 %")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 20),axis.text.y=element_text(colour = 'black', size = 20))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=20),axis.text.x = element_text(angle=0,colour = 'black', size=20,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=20), legend.title=element_text(size=20),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))

Order_means2_higher_5_percent_with_other




####Family
###Proportion based on just total OTUs list

###Proportion of unidentified based on OTUs presence
FAMILY_otus<-fungi_97pick_down_strict_initial_11[,1]
###Make OTUs abundance more than 1 to be one just as present
FAMILY_otus[FAMILY_otus>0]<-1
###Make OTUs abundance 0 be one just as present in the list of otus
FAMILY_otus[FAMILY_otus==0]<-1
FAMILY_otus[1:10]
####Everything should be one

FAMILY_otus<-aggregate(FAMILY_otus, by=list(fungi_97pick_down_strict_initial_11$Family), FUN=sum)
head(FAMILY_otus)
names(FAMILY_otus)[1]<- paste("FAMILY")
rownames(FAMILY_otus) <- FAMILY_otus[,1]
head(FAMILY_otus)
rownames(FAMILY_otus)


###
FAMILY_otus2 <- aggregate(FAMILY_otus[,2:ncol(FAMILY_otus)], by=list(FAMILY_otus$FAMILY), FUN=sum)
head(FAMILY_otus2)
names(FAMILY_otus2)[1]<- paste("FAMILY")
rownames(FAMILY_otus2) <- FAMILY_otus2[,1]
str(FAMILY_otus2)

########construct data frame for calculate, use the melt function
FAMILY_otus_sample_only <- FAMILY_otus2[,2:ncol(FAMILY_otus2)]
FAMILY_otus_sample_only_order <-FAMILY_otus_sample_only

###Percent of unidentified based on OTUs
FAMILY_otus_type_prop<-FAMILY_otus_sample_only_order
str(FAMILY_otus_type_prop)
FAMILY_otus_percent<-FAMILY_otus_type_prop*100/sum(FAMILY_otus_type_prop)
FAMILY_otus_per<-cbind(FAMILY_otus2,FAMILY_otus_percent)
FAMILY_otus_per
sum(FAMILY_otus_per$FAMILY_otus_percent)


####Proportion of unidentified based on presence/absence of OTUs
FAMILY_pres_abs<-fungi_97pick_down_strict_initial_11[,1:sample_l]
FAMILY_pres_abs[FAMILY_pres_abs>0]<-1

FAMILY_pres_abs[1:3]
FAMILY<-aggregate(FAMILY_pres_abs, by=list(fungi_97pick_down_strict_initial_11$Family), FUN=sum)
head(FAMILY)
names(FAMILY)[1]<- paste("FAMILY")
rownames(FAMILY) <- FAMILY[,1]
head(FAMILY)
rownames(FAMILY)


str(FAMILY)

###
FAMILY2 <- aggregate(FAMILY[,2:ncol(FAMILY)], by=list(FAMILY$FAMILY), FUN=sum)
head(FAMILY2)
names(FAMILY2)[1]<- paste("FAMILY")
rownames(FAMILY2) <- FAMILY2[,1]
str(FAMILY2)

########construct data frame for calculate, use the melt function
FAMILY_sample_only <- FAMILY2[,2:ncol(FAMILY2)]
FAMILY_sample_only_order <- FAMILY_sample_only[,order(colnames(FAMILY_sample_only))]

###Percent of unidentified based on OTUs
FAMILY_type_prop<-FAMILY_sample_only_order
FAMILY_type_prop$total<-rowSums(FAMILY_type_prop)
FAMILY_type_prop[290:308, 1:3]
FAMILY_type_prop$percent<-(FAMILY_type_prop$total/sum(FAMILY_type_prop$total) * 100)
FAMILY_type_prop[290:308, 418:420]
sum(FAMILY_type_prop$percent)

###############collape to Family level, aggregate works same as the collapse function in qiime

Family_agg <- aggregate(fungi_97pick_down_strict_initial_11[,1:sample_l], by=list(fungi_97pick_down_strict_initial_11$Family), FUN=sum)
head(Family_agg)
head(Family_agg)
names(Family_agg)[1]<- paste("Family")
rownames(Family_agg) <- Family_agg[,1]
head(Family_agg)
rownames(Family_agg)

str(Family_agg)

Family_agg2 <- aggregate(Family_agg[,2:ncol(Family_agg)], by=list(Family_agg$Family), FUN=sum)
head(Family_agg2)
names(Family_agg2)[1]<- paste("Family")
rownames(Family_agg2) <- Family_agg2[,1]

########construct data frame for calculate, use the melt function
Family_sample_only <- Family_agg2[,2:ncol(Family_agg2)]
Family_sample_only_order <- Family_sample_only[,order(colnames(Family_sample_only))]

###Percent of unidentified
Family_sample_only_order[290:308,1:3]

Family_prop<-Family_sample_only_order
Family_prop$total<-rowSums(Family_prop)
Family_prop[290:308, 418:419]
Family_prop$percent<-(Family_prop$total/sum(Family_prop$total) * 100)
sum(Family_prop$percent)

Family_prop[290:308, 418:420]


#########combine with map file
head(t(Family_sample_only_order))
match(rownames(t(Family_sample_only_order)),Down_strict_initial_map$SampleID)

Family_combine <- cbind(t(Family_sample_only_order),Down_strict_initial_map)
dim(Family_combine)


melt_id1 <- nrow(Family_sample_only_order)+1
melt_id1
melt_id2 <- ncol(Family_combine)
melt_id2

library(reshape2)
##Use melt to actually put the trait under one column and add set of rows with the ith information of ith trait 
## Basically melt is collapsing several columns to two (one for the columns and one for their values)
Family_melt <- reshape2::melt(Family_combine,id=melt_id1:melt_id2,variable.name = "Family",value.name = "Count")
str(Family_melt)

Family_melt$Forest
Family_time <- with(Family_melt,aggregate(Count,by=list(Time=Time,Family=Family,Species=Species,Forest=Forest),function(x) mean(x,na.rm=T)))
Try_Family_total<-Family_time
head(Try_Family_total)
sum(Try_Family_total$x)


l_try <- length(levels(Try_Family_total$Family))*2

str(Try_Family_total)
Initial_strict_down_try <- do.call("rbind", replicate(3, Try_Family_total[1:l_try,], simplify = FALSE))
Initial_strict_down_try$Forest <- c(paste(replicate(l_try,"Mature_forest")),paste(replicate(l_try,"Open_land")),paste(replicate(l_try,"Regenerating_forest")))


strict_down_try_plot <- rbind(Initial_strict_down_try,Try_Family_total[(l_try+1):nrow(Try_Family_total),])
strict_down_try_plot[] <- lapply(strict_down_try_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_try_plot$Forest <- factor(strict_down_try_plot$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))
strict_down_try_plot[] <- lapply(strict_down_try_plot, function(x) if(is.factor(x)) factor(x) else x)

strict_down_try_plot$Family #
strict_down_try_plot$x_percent<-(strict_down_try_plot$x)*100/sum(strict_down_try_plot$x)

levels(strict_down_try_plot$Species)
levels(strict_down_try_plot$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_try_plot$Forest)
levels(strict_down_try_plot$Forest)<-c("Open land", "Regenerating forest", "Mature forest")


Family_perct<- ddply(Family_time, .(Forest,Time,Species), mutate, Family_pct = x / sum(x) * 100)
summary(Family_perct)

head(Family_perct)
str(Family_perct)


levels(Family_perct$Time)
l_initial <- length(levels(Family_perct$Family))*2

Initial_strict_down_rep <- do.call("rbind", replicate(3, Family_perct[1:l_initial,], simplify = FALSE))
Initial_strict_down_rep$Forest <- c(paste(replicate(l_initial,"Mature_forest")),paste(replicate(l_initial,"Open_land")),paste(replicate(l_initial,"Regenerating_forest")))

strict_down_perct_plot <- rbind(Initial_strict_down_rep,Family_perct[(l_initial+1):nrow(Family_perct),])
strict_down_perct_plot[] <- lapply(strict_down_perct_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_perct_plot$Forest <- factor(strict_down_perct_plot$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))
strict_down_perct_plot[] <- lapply(strict_down_perct_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_perct_plot$Family ###
###Computing the average of Families across time and

str(strict_down_try_plot)

#library(dplyr)
Family_means1<-strict_down_try_plot%>%
  group_by(Forest,Species,Family,Time)%>%
  summarise(mean_percent=mean(x_percent,na.rm=T),
            sd_percent=sd(x_percent))
Family_means1
summary(Family_means1)
Family_means1[Family_means1$mean_percent>mean(Family_means1$mean_percent),]

Family_means1[Family_means1$mean_percent>1,]
#Family_means1$Species<-as.factor(Family_means1$Species)

#Family_means1$Species<-factor(Family_means1$Species,levels=c("Castanopsis_mekongensis","Litsea_cubeba"))
levels(Family_means1$Species)
levels(Family_means1$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(Family_means1$Forest)
levels(Family_means1$Forest)<-c("Open land", "Regenerating forest", "Mature forest")



#####Computation of mean for abundance
str(strict_down_perct_plot)
#library(dplyr)
summary(strict_down_perct_plot)

Family_means2<-strict_down_perct_plot%>%
  group_by(Forest,Species,Family,Time)%>%
  summarise(mean_percent=mean(Family_pct),
            sd_percent=sd(Family_pct))
Family_means2
summary(Family_means2)
Family_means2[Family_means2$mean_percent>10,]
Family_means2[Family_means2$mean_percent>5,]
Family_means2[Family_means2$mean_percent>7,]

####

Family_means2_higher_7_percent<-ggplot(data=Family_means2[Family_means2$mean_percent>7,],aes(x = Time, y = mean_percent, fill = Family))+
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  xlab("Incubation duration (Months)")+ 
  ylab("Abundance (%)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))+
  scale_fill_brewer(palette="Paired")

Family_means2_higher_7_percent


###Trying to fill the remaining of 100 % with Other
Family_means2[Family_means2$mean_percent>7,]

Family_means2$Family_new<-ifelse(Family_means2$mean_percent<=7, " All families <= 7 %", as.character(Family_means2$Family))
Family_means2$Family_new<-as.factor(Family_means2$Family_new)
levels(Family_means2$Family_new)

Family_means2$Family_new<-factor(Family_means2$Family_new, levels=c(" All families <= 7 %"," unidentified"," Coniochaetaceae"," Dermateaceae",
                                                                    " Herpotrichiellaceae"," Hyaloscyphaceae"," Hypocreaceae"," Meruliaceae",
                                                                    " Saccharomycetales_fam_Incertae_sedis"," Sclerotiniaceae",
                                                                    " Strophariaceae"," Trichocomaceae", " Xylariaceae"))


levels(Family_means2$Family_new)

Family_means2$Species
Family_means2$Species<-as.factor(Family_means2$Species)

Family_means2$Species
Family_means2$Species<-factor(Family_means2$Species,levels=c("Castanopsis_mekongensis","Litsea_cubeba"))

levels(Family_means2$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(Family_means2$Species)
levels(Family_means2$Forest)
levels(Family_means2$Forest)<-c("Open land", "Regenerating forest", "Mature forest")


Family_means2_higher_7_percent_with_other<-ggplot(data=Family_means2,aes(x = Time, y = mean_percent, fill = Family_new))+
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  xlab("Exposure time (Months)")+ 
  ylab("Abundance (%)")+
  #scale_fill_discrete(name="Most abundant families > 7 %")+
  labs(fill="Most abundant families > 7 %")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 20),axis.text.y=element_text(colour = 'black', size = 20))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=20),axis.text.x = element_text(angle=0,colour = 'black', size=20,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=20), legend.title=element_text(size=20),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))#+

Family_means2_higher_7_percent_with_other


####Genus level
####Genus

###Proportion based on just total OTUs list

###Proportion of unidentified based on OTUs presence
GENUS_otus<-fungi_97pick_down_strict_initial_11[,1]
###Make OTUs abundance more than 1 to be one just as present
GENUS_otus[GENUS_otus>0]<-1
###Make OTUs abundance 0 be one just as present in the list of otus
GENUS_otus[GENUS_otus==0]<-1
GENUS_otus[1:10]
####Everything should be one

GENUS_otus<-aggregate(GENUS_otus, by=list(fungi_97pick_down_strict_initial_11$Genus), FUN=sum)
head(GENUS_otus)
names(GENUS_otus)[1]<- paste("GENUS")
rownames(GENUS_otus) <- GENUS_otus[,1]
head(GENUS_otus)
rownames(GENUS_otus)


###
GENUS_otus2 <- aggregate(GENUS_otus[,2:ncol(GENUS_otus)], by=list(GENUS_otus$GENUS), FUN=sum)
head(GENUS_otus2)
names(GENUS_otus2)[1]<- paste("GENUS")
rownames(GENUS_otus2) <- GENUS_otus2[,1]
str(GENUS_otus2)

########construct data frame for calculate, use the melt function
GENUS_otus_sample_only <- GENUS_otus2[,2:ncol(GENUS_otus2)]
GENUS_otus_sample_only_order <-GENUS_otus_sample_only

###Percent of unidentified based on OTUs
GENUS_otus_type_prop<-GENUS_otus_sample_only_order
str(GENUS_otus_type_prop)
GENUS_otus_percent<-GENUS_otus_type_prop*100/sum(GENUS_otus_type_prop)
GENUS_otus_per<-cbind(GENUS_otus2,GENUS_otus_percent)
GENUS_otus_per[848:855,]
sum(GENUS_otus_per$GENUS_otus_percent)


####Proportion of unidentified based on presence/absence of OTUs
GENUS_pres_abs<-fungi_97pick_down_strict_initial_11[,1:sample_l]
GENUS_pres_abs[GENUS_pres_abs>0]<-1

GENUS_pres_abs[1:3]
GENUS<-aggregate(GENUS_pres_abs, by=list(fungi_97pick_down_strict_initial_11$Genus), FUN=sum)
head(GENUS)
names(GENUS)[1]<- paste("GENUS")
rownames(GENUS) <- GENUS[,1]
head(GENUS)
rownames(GENUS)


str(GENUS)

###
GENUS2 <- aggregate(GENUS[,2:ncol(GENUS)], by=list(GENUS$GENUS), FUN=sum)
head(GENUS2)
names(GENUS2)[1]<- paste("GENUS")
rownames(GENUS2) <- GENUS2[,1]
str(GENUS2)

########construct data frame for calculate, use the melt function
GENUS_sample_only <- GENUS2[,2:ncol(GENUS2)]
GENUS_sample_only_order <- GENUS_sample_only[,order(colnames(GENUS_sample_only))]

###Percent of unidentified based on OTUs
GENUS_type_prop<-GENUS_sample_only_order
GENUS_type_prop$total<-rowSums(GENUS_type_prop)
GENUS_type_prop[848:855, 1:3]
GENUS_type_prop$percent<-(GENUS_type_prop$total/sum(GENUS_type_prop$total) * 100)
GENUS_type_prop[848:855, 418:420]
sum(GENUS_type_prop$percent)


###############collape to Genus level, aggregate works same as the collapse function in qiime

Genus_agg <- aggregate(fungi_97pick_down_strict_initial_11[,1:sample_l], by=list(fungi_97pick_down_strict_initial_11$Genus), FUN=sum)
head(Genus_agg)
head(Genus_agg)
names(Genus_agg)[1]<- paste("Genus")
rownames(Genus_agg) <- Genus_agg[,1]
head(Genus_agg)
rownames(Genus_agg)

str(Genus_agg)
dim(Genus_agg)

Genus_agg2 <- aggregate(Genus_agg[,2:ncol(Genus_agg)], by=list(Genus_agg$Genus), FUN=sum)
head(Genus_agg2)
names(Genus_agg2)[1]<- paste("Genus")
rownames(Genus_agg2) <- Genus_agg2[,1]

########construct data frame for calculate, use the melt function
Genus_sample_only <- Genus_agg2[,2:ncol(Genus_agg2)]
Genus_sample_only_order <- Genus_sample_only[,order(colnames(Genus_sample_only))]

###Percent of unidentified
Genus_sample_only_order[848:855,1:3]

Genus_prop<-Genus_sample_only_order
Genus_prop$total<-rowSums(Genus_prop)
Genus_prop[848:855, 418:419]
Genus_prop$percent<-(Genus_prop$total/sum(Genus_prop$total) * 100)
sum(Genus_prop$percent)

Genus_prop[848:855, 418:420]

#########combine with map file
head(t(Genus_sample_only_order))
match(rownames(t(Genus_sample_only_order)),Down_strict_initial_map$SampleID)


Genus_combine <- cbind(t(Genus_sample_only_order),Down_strict_initial_map)
dim(Genus_combine)


melt_id1 <- nrow(Genus_sample_only_order)+1
melt_id1
melt_id2 <- ncol(Genus_combine)
melt_id2

library(reshape2)
##Use melt to actually put the trait under one column and add set of rows with the ith information of ith trait 
## Basically melt is collapsing several columns to two (one for the columns and one for their values)
Genus_melt <- reshape2::melt(Genus_combine,id=melt_id1:melt_id2,variable.name = "Genus",value.name = "Count")
str(Genus_melt)

Genus_melt$Forest
Genus_time <- with(Genus_melt,aggregate(Count,by=list(Time=Time,Genus=Genus,Species=Species,Forest=Forest),function(x) mean(x,na.rm=T)))
Try_Genus_total<-Genus_time
head(Try_Genus_total)
sum(Try_Genus_total$x)



l_try <- length(levels(Try_Genus_total$Genus))*2

str(Try_Genus_total)
Initial_strict_down_try <- do.call("rbind", replicate(3, Try_Genus_total[1:l_try,], simplify = FALSE))
Initial_strict_down_try$Forest <- c(paste(replicate(l_try,"Mature_forest")),paste(replicate(l_try,"Open_land")),paste(replicate(l_try,"Regenerating_forest")))


strict_down_try_plot <- rbind(Initial_strict_down_try,Try_Genus_total[(l_try+1):nrow(Try_Genus_total),])
strict_down_try_plot[] <- lapply(strict_down_try_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_try_plot$Forest <- factor(strict_down_try_plot$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))
strict_down_try_plot[] <- lapply(strict_down_try_plot, function(x) if(is.factor(x)) factor(x) else x)

strict_down_try_plot$Genus 
strict_down_try_plot$x_percent<-(strict_down_try_plot$x)*100/sum(strict_down_try_plot$x)
sum(strict_down_try_plot$x_percent)

levels(strict_down_try_plot$Species)
levels(strict_down_try_plot$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_try_plot$Forest)
levels(strict_down_try_plot$Forest)<-c("Open land", "Regenerating forest", "Mature forest")




Genus_perct<- ddply(Genus_time, .(Forest,Time,Species), mutate, Genus_pct = x / sum(x) * 100)
summary(Genus_perct)

head(Genus_perct)
str(Genus_perct)


levels(Genus_perct$Time)
l_initial <- length(levels(Genus_perct$Genus))*2

Initial_strict_down_rep <- do.call("rbind", replicate(3, Genus_perct[1:l_initial,], simplify = FALSE))
Initial_strict_down_rep$Forest <- c(paste(replicate(l_initial,"Mature_forest")),paste(replicate(l_initial,"Open_land")),paste(replicate(l_initial,"Regenerating_forest")))

strict_down_perct_plot <- rbind(Initial_strict_down_rep,Genus_perct[(l_initial+1):nrow(Genus_perct),])
strict_down_perct_plot[] <- lapply(strict_down_perct_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_perct_plot$Forest <- factor(strict_down_perct_plot$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))
strict_down_perct_plot[] <- lapply(strict_down_perct_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_perct_plot$Genus ###

summary(strict_down_perct_plot)
strict_down_perct_plot$Time
###Computing the average of Genera across time and

str(strict_down_try_plot)

#library(dplyr)
Genus_means1<-strict_down_try_plot%>%
  group_by(Species,Forest,Genus,Time)%>%
  summarise(mean_percent=mean(x_percent,na.rm=T),
            sd_percent=sd(x_percent))
Genus_means1
summary(Genus_means1)
Genus_means1[Genus_means1$mean_percent>mean(Genus_means1$mean_percent),]

Genus_means1[Genus_means1$mean_percent>1,]


Genus_means1_higher0.5percent<-ggplot(Genus_means1[Genus_means1$mean_percent>0.5,], aes(x = Time, y = mean_percent, fill = Genus)) +
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  ylab("Relative abundance (%)")+
  xlab("Incubation time")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))#+
  #scale_fill_brewer(palette="Paired")

Genus_means1_higher0.5percent


#####Computation of mean for relative abundance
str(strict_down_perct_plot)
#library(dplyr)
summary(strict_down_perct_plot)


Genus_means2<-strict_down_perct_plot%>%
  group_by(Forest,Species ,Genus,Time)%>%
  summarise(mean_percent=mean(Genus_pct),
            sd_percent=sd(Genus_pct))

Genus_means2
summary(Genus_means2)
Genus_means2[Genus_means2$mean_percent>10,]
Genus_means2[Genus_means2$mean_percent>5,]




####

Genus_means2_higher_5_percent<-ggplot(data=Genus_means2[Genus_means2$mean_percent>5,],aes(x = Time, y = mean_percent, fill = Genus))+
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  xlab("Incubation duration (Months)")+ 
  ylab("Abundance (%)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))#+
#scale_fill_brewer(palette="Paired")
Genus_means2_higher_5_percent



###Trying to fill the remaining of 100 % with Other
Genus_means2[Genus_means2$mean_percent>5,]
summary(Genus_means2[Genus_means2$mean_percent>5,])
Genus_means2$Genus_new<-ifelse(Genus_means2$mean_percent<=5, " All genera <= 5 %", as.character(Genus_means2$Genus))
Genus_means2$Genus_new<-as.factor(Genus_means2$Genus_new)
levels(Genus_means2$Genus_new)
Genus_means2$Genus_new<-factor(Genus_means2$Genus_new, levels=c(" All genera <= 5 %"," unidentified"," Botrytis"," Candida"," Coniochaeta"," Exophiala",
                                                                " Hyaloscypha"," Hypholoma"," Lecythophora"," Mollisia"," Nemania"," Neofusicoccum"," Penicillium",
                                                                " Pezicula", " Phlebia"," Sugiyamaella"," Trichoderma"," Xylomelasma"))

levels(Genus_means2$Genus_new)
Genus_means2$Species<-as.factor(Genus_means2$Species)
levels(Genus_means2$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(Genus_means2$Species)

Genus_means2$Forest<-as.factor(Genus_means2$Forest)
levels(Genus_means2$Forest)
levels(Genus_means2$Forest)<-c("Open land", "Regenerating forest", "Mature forest")


Genus_means2_higher_5_with_other<-ggplot(data=Genus_means2,aes(x = Time, y = mean_percent, fill = Genus_new))+
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
  xlab("Exposure time (Months)")+ 
  ylab("Abundance (%)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  #scale_fill_discrete(name="Most abundant genera > 5 %")+
  labs(fill="Most abundant genera > 5 %")+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 20),axis.text.y=element_text(colour = 'black', size = 20))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=20),axis.text.x = element_text(angle=0,colour = 'black', size=20,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=20), legend.title=element_text(size=20),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))#+
  #scale_fill_brewer(palette="Paired")

Genus_means2_higher_5_with_other


summary(Genus_means2)

####Putting figures for Order, Class, Family and Genera together

###Put figure in panels
##load libraries
library(ggplot2)
library(gridExtra)
library(reshape)
library(grid)


####Figure S4
tiff(filename="Figure S4 abundant taxa from Order to genera.tiff", res = 400, width=11000, height=7500, compression = "lzw")
grid.arrange(arrangeGrob(Class_means2_higher_5_with_other + theme(legend.position="right",
                                                                  plot.margin = unit( c(1,0.5,0,0) , units = "lines")),
                         Order_means2_higher_5_percent_with_other + theme(legend.position="right",
                                                                         plot.margin = unit( c(1,0.5,0,0) , units = "lines")),
                         Family_means2_higher_7_percent_with_other + theme(legend.position="right",
                                                                           plot.margin = unit( c(1,0.5,0,0) , units = "lines")),
                         Genus_means2_higher_5_with_other + theme(legend.position="right",
                                                                  plot.margin = unit( c(1,0.5,0,0) , units = "lines")),
                         nrow=2, ncol=2))
dev.off()

###Note A, B, C, D were added in powerpoint

#Class_means2_higher_5_with_other
#Order_means2_higher_5_percent_with_other
#Family_means2_higher_7_percent_with_other
#Genus_means2_higher_5_with_other


###Species

####Species level
####Species
###############collape to Class level, aggregate works same as the collapse function in qiime

Species_agg <- aggregate(fungi_97pick_down_strict_initial_11[,1:sample_l], by=list(fungi_97pick_down_strict_initial_11$Species), FUN=sum)
head(Species_agg)
head(Species_agg)
names(Species_agg)[1]<- paste("Species")
rownames(Species_agg) <- Species_agg[,1]
head(Species_agg)
rownames(Species_agg)

str(Species_agg)

Species_agg2 <- aggregate(Species_agg[,2:ncol(Species_agg)], by=list(Species_agg$Species), FUN=sum)
head(Species_agg2)
names(Species_agg2)[1]<- paste("Species")
rownames(Species_agg2) <- Species_agg2[,1]

########construct data frame for calculate, use the melt function
Species_sample_only <- Species_agg2[,2:ncol(Species_agg2)]
Species_sample_only_order <- Species_sample_only[,order(colnames(Species_sample_only))]

#########combine with map file
Species_combine <- cbind(t(Species_sample_only_order),Down_strict_initial_map)
dim(Species_combine)

#write.csv(Trait_combine, file="traits combine.csv")
melt_id1 <- nrow(Species_sample_only_order)+1
melt_id1
melt_id2 <- ncol(Species_combine)
melt_id2

library(reshape2)
##Use melt to actually put the trait under one column and add set of rows with the ith information of ith trait 
## Basically melt is collapsing several columns to two (one for the columns and one for their values)
Species_melt <- reshape2::melt(Species_combine,id=melt_id1:melt_id2,variable.name = "Species_F",value.name = "Count")
str(Species_melt)

Species_melt$Forest
Species_time <- with(Species_melt,aggregate(Count,by=list(Time=Time,Species_F=Species_F,Species=Species,Forest=Forest),function(x) mean(x,na.rm=T)))

Try_Species_total<-Species_time

head(Try_Species_total)
sum(Try_Species_total$x)

l_try <- length(levels(Try_Species_total$Species_F))*2

str(Try_Species_total)
Initial_strict_down_try <- do.call("rbind", replicate(3, Try_Species_total[1:l_try,], simplify = FALSE))
Initial_strict_down_try$Forest <- c(paste(replicate(l_try,"Mature_forest")),paste(replicate(l_try,"Open_land")),paste(replicate(l_try,"Regenerating_forest")))


strict_down_try_plot <- rbind(Initial_strict_down_try,Try_Species_total[(l_try+1):nrow(Try_Species_total),])
strict_down_try_plot[] <- lapply(strict_down_try_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_try_plot$Forest <- factor(strict_down_try_plot$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))
strict_down_try_plot[] <- lapply(strict_down_try_plot, function(x) if(is.factor(x)) factor(x) else x)

strict_down_try_plot$Species_F 

strict_down_try_plot$x_percent<-(strict_down_try_plot$x)*100/sum(strict_down_try_plot$x)

levels(strict_down_try_plot$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(strict_down_try_plot$Forest)<-c("Open land", "Regenerating forest", "Mature forest")

#ggplot(strict_down_try_plot, aes(x = Time, y = x_percent, fill = Species_F)) +
#  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +
#  ylab("Abundance (%)")+
#  xlab("Incubation time")+
#  theme(strip.text.x = element_text(size = 15, colour = "black"))+
#  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
#  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
#  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
#  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
#  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))



#library(plyr)

Species_perct<- ddply(Species_time, .(Forest,Time,Species), mutate, Species_pct = x / sum(x) * 100)
summary(Species_perct)

head(Species_perct)
#tail()
str(Species_perct)


levels(Species_perct$Time)
l_initial <- length(levels(Species_perct$Species_F))*2

Initial_strict_down_rep <- do.call("rbind", replicate(3, Species_perct[1:l_initial,], simplify = FALSE))
Initial_strict_down_rep$Forest <- c(paste(replicate(l_initial,"Mature_forest")),paste(replicate(l_initial,"Open_land")),paste(replicate(l_initial,"Regenerating_forest")))

strict_down_perct_plot <- rbind(Initial_strict_down_rep,Species_perct[(l_initial+1):nrow(Species_perct),])
strict_down_perct_plot[] <- lapply(strict_down_perct_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_perct_plot$Forest <- factor(strict_down_perct_plot$Forest,levels=c("Open_land","Regenerating_forest","Mature_forest"))
strict_down_perct_plot[] <- lapply(strict_down_perct_plot, function(x) if(is.factor(x)) factor(x) else x)
strict_down_perct_plot$Species_F ###


###Computing the average of Species across time and

str(strict_down_try_plot)
levels(strict_down_try_plot$Forest)

#library(dplyr)
Species_means1<-strict_down_try_plot%>%
  group_by(Forest,Species_F,Species,Time)%>%
  summarise(mean_percent=mean(x_percent,na.rm=T),
            sd_percent=sd(x_percent))
Species_means1
summary(Species_means1)
Species_means1[Species_means1$mean_percent>mean(Species_means1$mean_percent),]

Species_means1[Species_means1$mean_percent>0.5,]

levels(Species_means1$Forest)


#####Computation of mean for relative abundance
str(strict_down_perct_plot)
#library(dplyr)
summary(strict_down_perct_plot)

Species_means2<-strict_down_perct_plot%>%
  group_by(Forest,Species ,Species_F,Time)%>%
  summarise(mean_percent=mean(Species_pct),
            sd_percent=sd(Species_pct))
Species_means2
summary(Species_means2)
Species_means2[Species_means2$mean_percent>10,]
Species_means2[Species_means2$mean_percent>5,]

####

Species_means2_higher_5_percent<-ggplot(data=Species_means2[Species_means2$mean_percent>5,],aes(x = Time, y = mean_percent, fill = Species_F))+
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +xlab("Incubation duration (Months)")+ ylab("Relative abundance (%)")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))

Species_means2_higher_5_percent

#dev.off()




###Trying to fill the remaining of 100 % with Other
Species_means2[Species_means2$mean_percent>5,]
Species_means2$Species_new<-ifelse(Species_means2$mean_percent<=5, " Others", as.character(Species_means2$Species_F))
Species_means2$Species_new<-as.factor(Species_means2$Species_new)
levels(Species_means2$Species_new)
#Species_means2$Species_new<-factor(Species_means2$Species_new, levels=c(" Others"," Agaricales_sp" ," Agaricomycetes_sp"," Ascomycota_sp"," Beauveria_kipukae",
#                                                                        " Botrytis_caroliniana"," Candida_sp"," Capnodiales_sp"," Mollisia_dextrinospora"," Phlebia_sp"," Pleosporales_sp",
#                                                                        " Psathyrellaceae_sp"," Scytalidium_album"," Sordariomycetes_sp"," Trichoderma_neokoningii"," Xylariaceae_sp" ))

Species_means2$Species<-as.factor(Species_means2$Species)
levels(Species_means2$Species)<-c("Castanopsis mekongensis","Litsea cubeba")
levels(Species_means2$Species)

Species_means2$Forest<-as.factor(Species_means2$Forest)
levels(Species_means2$Forest)
levels(Species_means2$Forest)<-c("Open land", "Regenerating forest", "Mature forest")


Species_means2_higher_5_percent<-ggplot(data=Species_means2,aes(x = Time, y = mean_percent, fill = Species_new))+
  geom_bar(stat = 'identity',  position = 'stack')+ facet_wrap(~ Species+Forest,ncol=3) +xlab("Incubation duration (Months)")+ ylab("Relative abundance (%)")+
  scale_fill_discrete(name="Species")+
  theme(strip.text.x = element_text(size = 15, colour = "black"))+
  theme(axis.title.y = element_text(vjust=0.8,hjust=0.5,colour = 'black', size = 15),axis.text.y=element_text(colour = 'black', size = 15))+
  theme(axis.title.x = element_text(colour = 'black', angle=0,size=15),axis.text.x = element_text(angle=0,colour = 'black', size=15,vjust=1,hjust=0.6))+
  theme(legend.background = element_rect(fill="white",colour=NULL, size=.5),legend.text=element_text(size=15), legend.title=element_text(size=15),legend.direction="vertical")+
  theme(legend.key=element_rect(colour="white",size=0.1),legend.key.size=unit(.5,"cm"))+
  scale_x_discrete(breaks=c("0mo","18mo","36mo"),labels=c("0 mo","18 mo","36 mo"))

Species_means2_higher_5_percent


####VENN DIAGRAM 

getwd()
####import otu table with guild and taxonomy level split into each column
mengsong_97pick <- read.table("mengsong_97closed_guilds_r.csv", sep=",", row.names=1,header=T, check.names=F,blank.lines.skip = FALSE)   
dim(mengsong_97pick)

###import map file
map <- read.table("mengsong_map_r_chem_SBB.csv",sep=",", header=T, check.names=F,blank.lines.skip = FALSE)
rownames(map) <- map[,1]
dim(map)
#####create subset vectors
subset_sample <- rownames(map)[which(map$control=="No")]
subset_control <- rownames(map)[which(map$control=="Yes")]

sample_map <- map[subset_sample,]
sample_map <- sample_map[order(rownames(sample_map)),]
sample_map$Forest <- factor(sample_map$Forest)
sample_map$Forest

sample_map$Species
sample_map$Species<-as.factor(sample_map$Species)
levels(sample_map$Species) ###check the species level

###check the data to make sure the otu table has the correct number of samples
mengsong_97pick[1:2,834:835]
#####create map for separating mengsong samples
Down_map <- sample_map[which(sample_map$SampleType=="Down" |sample_map$SampleType=="Initial"),]
Down_map[] <- lapply(Down_map, function(x) if(is.factor(x)) factor(x) else x)##make sure the factor levels are cosistent with the subset dataframe
levels(Down_map$Forest)

Pilot_map <- sample_map[sample_map$SampleType=="PilotSample",]
####exclude pilot map
mengsong_map <- sample_map[!sample_map$SampleType %in% "PilotSample",]
mengsong_map[] <- lapply(mengsong_map, function(x) if(is.factor(x)) factor(x) else x)

###########filter out non fungi otus
subset_fungi <- rownames(mengsong_97pick)[which(mengsong_97pick$Kingdom=="Fungi")]
fungi_97pick <- mengsong_97pick[subset_fungi,]
dim(fungi_97pick)

fungi_97pick_mengsong <- fungi_97pick[,!colnames(fungi_97pick) %in% subset_control & !colnames(fungi_97pick) %in% Pilot_map[,1]]####only get mengsong sample otus
fungi_97pick_mengsong[1:2,636:637]
sample_l <- length(mengsong_map[,1])
#####exclude less than 10 otus 
fungi_97pick_mengsong_11 <- subset(fungi_97pick_mengsong,rowSums(fungi_97pick_mengsong[,1:sample_l]) >10)
fungi_97pick_mengsong_11[1:2,636:637]

#########################################work on strict down
Down_strict_initial_map <- Down_map[!Down_map$Position %in%"up",]
Down_strict_initial_map[] <- lapply(Down_strict_initial_map, function(x) if(is.factor(x)) factor(x) else x)##make sure the factor levels are cosistent with the subset dataframe

fungi_97pick_down_strict_initial_11 <- fungi_97pick_mengsong_11[,rownames(Down_strict_initial_map)]
dim(fungi_97pick_down_strict_initial_11)

#####create subset sampleID for time point OTU calculation
subset_down_initial <- rownames(Down_strict_initial_map)[which(Down_strict_initial_map$Time=="0mo")]
subset_down_18 <- rownames(Down_strict_initial_map)[which(Down_strict_initial_map$Time=="18mo")]
subset_down_36 <- rownames(Down_strict_initial_map)[which(Down_strict_initial_map$Time=="36mo")]


#############subset otu table to initial time point
subset_down_initial_otu <- fungi_97pick_down_strict_initial_11[,subset_down_initial]
subset_down_initial_otu[,"Total"] <- rowSums(subset_down_initial_otu)###add the total summ as a column
subset_down_initial_otu_order <- subset_down_initial_otu[order(-subset_down_initial_otu$Total),]###sort by descending order

###exclude 0 otus
biovenn_down_initial_otu_all <- subset_down_initial_otu_order[subset_down_initial_otu_order$Total>0,]


#############subset otu table to 18mo time point
subset_down_18mo_otu <- fungi_97pick_down_strict_initial_11[,subset_down_18]
subset_down_18mo_otu[,"Total"] <- rowSums(subset_down_18mo_otu)###add the total summ as a column
subset_down_18mo_otu_order <- subset_down_18mo_otu[order(-subset_down_18mo_otu$Total),]###sort by descending order

###exclude 0 otus
biovenn_down_18mo_otu_all <- subset_down_18mo_otu_order[subset_down_18mo_otu_order$Total>0,]
dim(biovenn_down_18mo_otu_all)


#############subset otu table to 36 mo time point
subset_down_36mo_otu <- fungi_97pick_down_strict_initial_11[,subset_down_36]
subset_down_36mo_otu[,"Total"] <- rowSums(subset_down_36mo_otu)
subset_down_36mo_otu_order <- subset_down_36mo_otu[order(-subset_down_36mo_otu$Total),]###sort by descending order
#dim(biovenn_down_36mo_otu)
###exclude 0 otus
biovenn_down_36mo_otu_all <- subset_down_36mo_otu_order[subset_down_36mo_otu_order$Total>0,]
dim(biovenn_down_36mo_otu_all)



#######separate species
################litsea
Down_strict_initial_Litsea_map <-Down_strict_initial_map[Down_strict_initial_map$Species=="Litsea_cubeba",]
fungi_97pick_down_strict_initial_Litsea_11 <- fungi_97pick_mengsong_11[,rownames(Down_strict_initial_Litsea_map)]
dim(fungi_97pick_down_strict_initial_Litsea_11)

#########################################work on strick down Litsea
#####create subset sampleID for time point OTU calculation
subset_down_initial_Litsea <- rownames(Down_strict_initial_Litsea_map)[which(Down_strict_initial_Litsea_map$Time=="0mo")]
subset_down_18_LC <- rownames(Down_strict_initial_Litsea_map)[which(Down_strict_initial_Litsea_map$Time=="18mo")]
length(subset_down_18_LC)
subset_down_36_LC <- rownames(Down_strict_initial_Litsea_map)[which(Down_strict_initial_Litsea_map$Time=="36mo")]

###subset by habitat type
subset_down_18_mature_LC <- rownames(Down_strict_initial_Litsea_map)[which(Down_strict_initial_Litsea_map$Time=="18mo"&Down_strict_initial_Litsea_map$Forest=="Mature_forest")]
length(subset_down_18_mature_LC)
subset_down_18_regenerating_LC <- rownames(Down_strict_initial_Litsea_map)[which(Down_strict_initial_Litsea_map$Time=="18mo"&Down_strict_initial_Litsea_map$Forest=="Regenerating_forest")]
length(subset_down_18_regenerating_LC)
subset_down_18_open_LC <- rownames(Down_strict_initial_Litsea_map)[which(Down_strict_initial_Litsea_map$Time=="18mo"&Down_strict_initial_Litsea_map$Forest=="Open_land")]
length(subset_down_18_open_LC)


subset_down_36 <- rownames(Down_strict_initial_Litsea_map)[which(Down_strict_initial_Litsea_map$Time=="36mo")]
length(subset_down_36)

subset_down_36_mature_LC <- rownames(Down_strict_initial_Litsea_map)[which(Down_strict_initial_Litsea_map$Time=="36mo"&Down_strict_initial_Litsea_map$Forest=="Mature_forest")]
length(subset_down_36_mature_LC)
subset_down_36_regenerating_LC <- rownames(Down_strict_initial_Litsea_map)[which(Down_strict_initial_Litsea_map$Time=="36mo"&Down_strict_initial_Litsea_map$Forest=="Regenerating_forest")]
length(subset_down_36_regenerating_LC)
subset_down_36_open_LC <- rownames(Down_strict_initial_Litsea_map)[which(Down_strict_initial_Litsea_map$Time=="36mo"&Down_strict_initial_Litsea_map$Forest=="Open_land")]
length(subset_down_36_open_LC)

#############subset otu table to initial_Litsea time point
subset_down_initial_Litsea_otu <- fungi_97pick_down_strict_initial_Litsea_11[,subset_down_initial_Litsea]
subset_down_initial_Litsea_otu[,"Total"] <- rowSums(subset_down_initial_Litsea_otu)###add the total summ as a column
subset_down_initial_Litsea_otu_order <- subset_down_initial_Litsea_otu[order(-subset_down_initial_Litsea_otu$Total),]###sort by descending order
###exclude 0 otus
biovenn_down_initial_Litsea_otu_all <- subset_down_initial_Litsea_otu_order[subset_down_initial_Litsea_otu_order$Total>0,]
###export for biovenn plotting
write.csv(biovenn_down_initial_Litsea_otu_all,file="biovenn_down_initial_Litsea_otu_all.csv",row.names = T)

#############subset otu table to 18mo timepoint
subset_down_18mo_otu <- fungi_97pick_down_strict_initial_Litsea_11[,subset_down_18_LC]
subset_down_18mo_otu[,"Total"] <- rowSums(subset_down_18mo_otu)###add the total summ as a column
subset_down_18mo_otu_order <- subset_down_18mo_otu[order(-subset_down_18mo_otu$Total),]###sort by descending order

###exclude 0 otus
biovenn_down_18mo_otu_all <- subset_down_18mo_otu_order[subset_down_18mo_otu_order$Total>0,]
###export for biovenn plotting
write.csv(biovenn_down_18mo_otu_all,file="biovenn_down_18mo_Litsea_otu_all.csv",row.names = TRUE)

#########subdivide by habitat type at 18mo

###Mature forest
subset_down_18mo_otu_litsea_mature <- fungi_97pick_down_strict_initial_Litsea_11[,subset_down_18_mature_LC]
subset_down_18mo_otu_litsea_mature[,"Total"] <- rowSums(subset_down_18mo_otu_litsea_mature)###add the total summ as a column
subset_down_18mo_otu_litsea_mature_order <- subset_down_18mo_otu_litsea_mature[order(-subset_down_18mo_otu_litsea_mature$Total),]###sort by descending order
dim(subset_down_18mo_otu_litsea_mature_order)

###exclude 0 otus
litsea_down_18mo_otu_mature_all <- subset_down_18mo_otu_litsea_mature_order[subset_down_18mo_otu_litsea_mature_order$Total>0,]
dim(litsea_down_18mo_otu_mature_all)
###export for biovenn plotting
write.csv(litsea_down_18mo_otu_mature_all,file="litsea_down_18mo_otu_mature_all.csv",row.names = TRUE)

###Regenerating forest
subset_down_18mo_otu_litsea_regenerating <- fungi_97pick_down_strict_initial_Litsea_11[,subset_down_18_regenerating_LC]
subset_down_18mo_otu_litsea_regenerating[,"Total"] <- rowSums(subset_down_18mo_otu_litsea_regenerating)###add the total summ as a column
subset_down_18mo_otu_litsea_regenerating_order <- subset_down_18mo_otu_litsea_regenerating[order(-subset_down_18mo_otu_litsea_regenerating$Total),]###sort by descending order
dim(subset_down_18mo_otu_litsea_regenerating_order)


###exclude 0 otus
litsea_down_18mo_otu_regenerating_all <- subset_down_18mo_otu_litsea_regenerating_order[subset_down_18mo_otu_litsea_regenerating_order$Total>0,]
dim(litsea_down_18mo_otu_regenerating_all)
###export for biovenn plotting
write.csv(litsea_down_18mo_otu_regenerating_all,file="litsea_down_18mo_otu_regenerating_all.csv",row.names = TRUE)

###Open land

subset_down_18mo_otu_litsea_open <- fungi_97pick_down_strict_initial_Litsea_11[,subset_down_18_open_LC]
subset_down_18mo_otu_litsea_open[,"Total"] <- rowSums(subset_down_18mo_otu_litsea_open)###add the total summ as a column
subset_down_18mo_otu_litsea_open_order <- subset_down_18mo_otu_litsea_open[order(-subset_down_18mo_otu_litsea_open$Total),]###sort by descending order
dim(subset_down_18mo_otu_litsea_open_order)


###exclude 0 otus
litsea_down_18mo_otu_open_all <- subset_down_18mo_otu_litsea_open_order[subset_down_18mo_otu_litsea_open_order$Total>0,]
dim(litsea_down_18mo_otu_open_all)
###export for biovenn plotting
write.csv(litsea_down_18mo_otu_open_all,file="litsea_down_18mo_otu_open_all.csv",row.names = TRUE)




#############subset otu table to 36mo time point
subset_down_36mo_otu <- fungi_97pick_down_strict_initial_Litsea_11[,subset_down_36_LC]
subset_down_36mo_otu[,"Total"] <- rowSums(subset_down_36mo_otu)###add the total summ as a column
subset_down_36mo_otu_order <- subset_down_36mo_otu[order(-subset_down_36mo_otu$Total),]###sort by descending order
###exclude 0 otus
biovenn_down_36mo_otu_all <- subset_down_36mo_otu_order[subset_down_36mo_otu_order$Total>0,]

###export for biovenn plotting
write.csv(biovenn_down_36mo_otu_all,file="biovenn_down_36mo_Litsea_otu_all.csv",row.names = TRUE)

#########subdivide by habitat type at 36mo

###Mature forest
subset_down_36mo_otu_litsea_mature <- fungi_97pick_down_strict_initial_Litsea_11[,subset_down_36_mature_LC]
subset_down_36mo_otu_litsea_mature[,"Total"] <- rowSums(subset_down_36mo_otu_litsea_mature)###add the total summ as a column
subset_down_36mo_otu_litsea_mature_order <- subset_down_36mo_otu_litsea_mature[order(-subset_down_36mo_otu_litsea_mature$Total),]###sort by descending order
dim(subset_down_36mo_otu_litsea_mature_order)

###exclude 0 otus
litsea_down_36mo_otu_mature_all <- subset_down_36mo_otu_litsea_mature_order[subset_down_36mo_otu_litsea_mature_order$Total>0,]
dim(litsea_down_36mo_otu_mature_all)
###export for biovenn plotting
write.csv(litsea_down_36mo_otu_mature_all,file="litsea_down_36mo_otu_mature_all.csv",row.names = TRUE)

###Regenerating forest
subset_down_36mo_otu_litsea_regenerating <- fungi_97pick_down_strict_initial_Litsea_11[,subset_down_36_regenerating_LC]
subset_down_36mo_otu_litsea_regenerating[,"Total"] <- rowSums(subset_down_36mo_otu_litsea_regenerating)###add the total summ as a column
subset_down_36mo_otu_litsea_regenerating_order <- subset_down_36mo_otu_litsea_regenerating[order(-subset_down_36mo_otu_litsea_regenerating$Total),]###sort by descending order
dim(subset_down_36mo_otu_litsea_regenerating_order)


###exclude 0 otus
litsea_down_36mo_otu_regenerating_all <- subset_down_36mo_otu_litsea_regenerating_order[subset_down_36mo_otu_litsea_regenerating_order$Total>0,]
dim(litsea_down_36mo_otu_regenerating_all)
###export for biovenn plotting
write.csv(litsea_down_36mo_otu_regenerating_all,file="litsea_down_36mo_otu_regenerating_all.csv",row.names = TRUE)

###Open land

subset_down_36mo_otu_litsea_open <- fungi_97pick_down_strict_initial_Litsea_11[,subset_down_36_open_LC]
subset_down_36mo_otu_litsea_open[,"Total"] <- rowSums(subset_down_36mo_otu_litsea_open)###add the total summ as a column
subset_down_36mo_otu_litsea_open_order <- subset_down_36mo_otu_litsea_open[order(-subset_down_36mo_otu_litsea_open$Total),]###sort by descending order
dim(subset_down_36mo_otu_litsea_open_order)


###exclude 0 otus
litsea_down_36mo_otu_open_all <- subset_down_36mo_otu_litsea_open_order[subset_down_36mo_otu_litsea_open_order$Total>0,]
dim(litsea_down_36mo_otu_open_all)
###export for biovenn plotting
write.csv(litsea_down_36mo_otu_open_all,file="litsea_down_36mo_otu_open_all.csv",row.names = TRUE)


###################subset to castanopsis, prepare dataset for biovenn
Down_strict_initial_Castanopsis_map <-Down_strict_initial_map[Down_strict_initial_map$Species=="Castanopsis_mekongensis",]
fungi_97pick_down_strict_initial_Castanopsis_11 <- fungi_97pick_mengsong_11[,rownames(Down_strict_initial_Castanopsis_map)]
dim(fungi_97pick_down_strict_initial_Castanopsis_11)



#########################################work on strict down 
#####create subset sampleID for time point OTU calculation
subset_down_initial_Castanopsis <- rownames(Down_strict_initial_Castanopsis_map)[which(Down_strict_initial_Castanopsis_map$Time=="0mo")]
subset_down_18_CM <- rownames(Down_strict_initial_Castanopsis_map)[which(Down_strict_initial_Castanopsis_map$Time=="18mo")]
length(subset_down_18_CM)
subset_down_36_CM <- rownames(Down_strict_initial_Castanopsis_map)[which(Down_strict_initial_Castanopsis_map$Time=="36mo")]
length(subset_down_36_CM)

###subset by habitat type
subset_down_18_mature_CM <- rownames(Down_strict_initial_Castanopsis_map)[which(Down_strict_initial_Castanopsis_map$Time=="18mo"&Down_strict_initial_Castanopsis_map$Forest=="Mature_forest")]
length(subset_down_18_mature_CM)
subset_down_18_regenerating_CM <- rownames(Down_strict_initial_Castanopsis_map)[which(Down_strict_initial_Castanopsis_map$Time=="18mo"&Down_strict_initial_Castanopsis_map$Forest=="Regenerating_forest")]
length(subset_down_18_regenerating_CM)
subset_down_18_open_CM <- rownames(Down_strict_initial_Castanopsis_map)[which(Down_strict_initial_Castanopsis_map$Time=="18mo"&Down_strict_initial_Castanopsis_map$Forest=="Open_land")]
length(subset_down_18_open_CM)


subset_down_36_mature_CM <- rownames(Down_strict_initial_Castanopsis_map)[which(Down_strict_initial_Castanopsis_map$Time=="36mo"&Down_strict_initial_Castanopsis_map$Forest=="Mature_forest")]
length(subset_down_36_mature_CM)

subset_down_36_regenerating_CM <- rownames(Down_strict_initial_Castanopsis_map)[which(Down_strict_initial_Castanopsis_map$Time=="36mo"&Down_strict_initial_Castanopsis_map$Forest=="Regenerating_forest")]
length(subset_down_36_regenerating_CM)
subset_down_36_open_CM <- rownames(Down_strict_initial_Castanopsis_map)[which(Down_strict_initial_Castanopsis_map$Time=="36mo"&Down_strict_initial_Castanopsis_map$Forest=="Open_land")]
length(subset_down_36_open_CM)



#############subset otu table to initial_Castanopsis timepoint
subset_down_initial_Castanopsis_otu <- fungi_97pick_down_strict_initial_Castanopsis_11[,subset_down_initial_Castanopsis]
subset_down_initial_Castanopsis_otu[,"Total"] <- rowSums(subset_down_initial_Castanopsis_otu)###add the total summ as a column
subset_down_initial_Castanopsis_otu_order <- subset_down_initial_Castanopsis_otu[order(-subset_down_initial_Castanopsis_otu$Total),]###sort by descending order
###exclude 0 otus
biovenn_down_initial_Castanopsis_otu_all <- subset_down_initial_Castanopsis_otu_order[subset_down_initial_Castanopsis_otu_order$Total>0,]
###export for biovenn plotting
write.csv(biovenn_down_initial_Castanopsis_otu_all,file="biovenn_down_initial_Castanopsis_otu_all.csv",row.names = TRUE)

#############subset otu table to 18mo timepoint
subset_down_18mo_otu <- fungi_97pick_down_strict_initial_Castanopsis_11[,subset_down_18_CM]

subset_down_18mo_otu[,"Total"] <- rowSums(subset_down_18mo_otu)###add the total summ as a column
subset_down_18mo_otu_order <- subset_down_18mo_otu[order(-subset_down_18mo_otu$Total),]###sort by descending order

###exclude 0 otus
biovenn_down_18mo_otu_all <- subset_down_18mo_otu_order[subset_down_18mo_otu_order$Total>0,]
###export for biovenn plotting
write.csv(biovenn_down_18mo_otu_all,file="biovenn_down_18mo_Castanopsis_otu_all.csv",row.names = TRUE)

#########subdivide by habitat type at 18mo 

###Mature forest
subset_down_18mo_otu_Castanopsis_mature <- fungi_97pick_down_strict_initial_Castanopsis_11[,subset_down_18_mature_CM]
subset_down_18mo_otu_Castanopsis_mature[,"Total"] <- rowSums(subset_down_18mo_otu_Castanopsis_mature)###add the total summ as a column
subset_down_18mo_otu_Castanopsis_mature_order <- subset_down_18mo_otu_Castanopsis_mature[order(-subset_down_18mo_otu_Castanopsis_mature$Total),]###sort by descending order
dim(subset_down_18mo_otu_Castanopsis_mature_order)

###exclude 0 otus
Castanopsis_down_18mo_otu_mature_all <- subset_down_18mo_otu_Castanopsis_mature_order[subset_down_18mo_otu_Castanopsis_mature_order$Total>0,]
dim(Castanopsis_down_18mo_otu_mature_all)
###export for biovenn plotting
write.csv(Castanopsis_down_18mo_otu_mature_all,file="Castanopsis_down_18mo_otu_mature_all.csv",row.names = TRUE)

###Regenerating forest
subset_down_18mo_otu_Castanopsis_regenerating <- fungi_97pick_down_strict_initial_Castanopsis_11[,subset_down_18_regenerating_CM]
subset_down_18mo_otu_Castanopsis_regenerating[,"Total"] <- rowSums(subset_down_18mo_otu_Castanopsis_regenerating)###add the total summ as a column
subset_down_18mo_otu_Castanopsis_regenerating_order <- subset_down_18mo_otu_Castanopsis_regenerating[order(-subset_down_18mo_otu_Castanopsis_regenerating$Total),]###sort by descending order
dim(subset_down_18mo_otu_Castanopsis_regenerating_order)


###exclude 0 otus
Castanopsis_down_18mo_otu_regenerating_all <- subset_down_18mo_otu_Castanopsis_regenerating_order[subset_down_18mo_otu_Castanopsis_regenerating_order$Total>0,]
dim(Castanopsis_down_18mo_otu_regenerating_all)
###export for biovenn plotting
write.csv(Castanopsis_down_18mo_otu_regenerating_all,file="Castanopsis_down_18mo_otu_regenerating_all.csv",row.names = TRUE)

###Open land

subset_down_18mo_otu_Castanopsis_open <- fungi_97pick_down_strict_initial_Castanopsis_11[,subset_down_18_open_CM]
subset_down_18mo_otu_Castanopsis_open[,"Total"] <- rowSums(subset_down_18mo_otu_Castanopsis_open)###add the total summ as a column
subset_down_18mo_otu_Castanopsis_open_order <- subset_down_18mo_otu_Castanopsis_open[order(-subset_down_18mo_otu_Castanopsis_open$Total),]###sort by descending order
dim(subset_down_18mo_otu_Castanopsis_open_order)


###exclude 0 otus
Castanopsis_down_18mo_otu_open_all <- subset_down_18mo_otu_Castanopsis_open_order[subset_down_18mo_otu_Castanopsis_open_order$Total>0,]
dim(Castanopsis_down_18mo_otu_open_all)
###export for biovenn plotting
write.csv(Castanopsis_down_18mo_otu_open_all,file="Castanopsis_down_18mo_otu_open_all.csv",row.names = TRUE)





#############subset otu table to 36mo timepoint
subset_down_36mo_otu <- fungi_97pick_down_strict_initial_Castanopsis_11[,subset_down_36_CM]
subset_down_36mo_otu[,"Total"] <- rowSums(subset_down_36mo_otu)###add the total summ as a column
subset_down_36mo_otu_order <- subset_down_36mo_otu[order(-subset_down_36mo_otu$Total),]###sort by descending order
###exclude 0 otus
biovenn_down_36mo_otu_all <- subset_down_36mo_otu_order[subset_down_36mo_otu_order$Total>0,]
dim(biovenn_down_36mo_otu_all)
###export for biovenn plotting
write.csv(biovenn_down_36mo_otu_all,file="biovenn_36mo_Castanopsis_otu.csv",row.names = TRUE)



#########subdivide by habitat type at 36mo

###Mature forest
subset_down_36mo_otu_Castanopsis_mature <- fungi_97pick_down_strict_initial_Castanopsis_11[,subset_down_36_mature_CM]
subset_down_36mo_otu_Castanopsis_mature[,"Total"] <- rowSums(subset_down_36mo_otu_Castanopsis_mature)###add the total summ as a column
subset_down_36mo_otu_Castanopsis_mature_order <- subset_down_36mo_otu_Castanopsis_mature[order(-subset_down_36mo_otu_Castanopsis_mature$Total),]###sort by descending order
dim(subset_down_36mo_otu_Castanopsis_mature_order)

###exclude 0 otus
Castanopsis_down_36mo_otu_mature_all <- subset_down_36mo_otu_Castanopsis_mature_order[subset_down_36mo_otu_Castanopsis_mature_order$Total>0,]
dim(Castanopsis_down_36mo_otu_mature_all)
###export for venny plotting
write.csv(Castanopsis_down_36mo_otu_mature_all,file="Castanopsis_down_36mo_otu_mature_all.csv",row.names = TRUE)

###Regenerating forest
subset_down_36mo_otu_Castanopsis_regenerating <- fungi_97pick_down_strict_initial_Castanopsis_11[,subset_down_36_regenerating_CM]
subset_down_36mo_otu_Castanopsis_regenerating[,"Total"] <- rowSums(subset_down_36mo_otu_Castanopsis_regenerating)###add the total summ as a column
subset_down_36mo_otu_Castanopsis_regenerating_order <- subset_down_36mo_otu_Castanopsis_regenerating[order(-subset_down_36mo_otu_Castanopsis_regenerating$Total),]###sort by descending order
dim(subset_down_36mo_otu_Castanopsis_regenerating_order)


###exclude 0 otus
Castanopsis_down_36mo_otu_regenerating_all <- subset_down_36mo_otu_Castanopsis_regenerating_order[subset_down_36mo_otu_Castanopsis_regenerating_order$Total>0,]
dim(Castanopsis_down_36mo_otu_regenerating_all)
###export for venny plotting
write.csv(Castanopsis_down_36mo_otu_regenerating_all,file="Castanopsis_down_36mo_otu_regenerating_all.csv",row.names = TRUE)

###Open land

subset_down_36mo_otu_Castanopsis_open <- fungi_97pick_down_strict_initial_Castanopsis_11[,subset_down_36_open_CM]
subset_down_36mo_otu_Castanopsis_open[,"Total"] <- rowSums(subset_down_36mo_otu_Castanopsis_open)###add the total summ as a column
subset_down_36mo_otu_Castanopsis_open_order <- subset_down_36mo_otu_Castanopsis_open[order(-subset_down_36mo_otu_Castanopsis_open$Total),]###sort by descending order
dim(subset_down_36mo_otu_Castanopsis_open_order)


###exclude 0 otus
Castanopsis_down_36mo_otu_open_all <- subset_down_36mo_otu_Castanopsis_open_order[subset_down_36mo_otu_Castanopsis_open_order$Total>0,]
dim(Castanopsis_down_36mo_otu_open_all)
###export for venny plotting
write.csv(Castanopsis_down_36mo_otu_open_all,file="Castanopsis_down_36mo_otu_open_all.csv",row.names = TRUE)


##Venn diagram is possible to draw using package VennDiagram 
###However, we use the https://www.biovenn.nl/ where we input files we constructed above.
