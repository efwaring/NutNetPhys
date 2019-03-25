library(ggplot2)
library(dplyr)
library(tidyr)
library(nlme)


full_biomass<-read.csv("full-biomass-nutrients-06-December-2017.csv")
comb<- read.csv("comb-by-plot-09-April-2018.csv")
CRU <-read.csv("CRU-annual_2018-07-06.csv")

comb$proportion_parN <- as.numeric(comb$proportion_par)

full_biomass <- full_biomass_nutrients_06_December_2017

full_biomass$pct_N <- as.numeric(full_biomass$pct_N)

full_biomass$AGN <- full_biomass$mass*(full_biomass$pct_N*0.01)

sites <- unique(full_biomass$site_code)
year <- unique(full_biomass$year)

comb1 <- comb %>% filter(site_code %in% sites) 
cru1 <- CRU %>% filter(site_code %in% sites)




all <- left_join(full_biomass, comb1, match="first")
all <- left_join(all, cru1, match="first")




all <- all %>% filter(category!="LITTER") %>% filter(category!="LEGUME") %>%
  filter(category!="WOODY")

allN <- all %>% filter(grepl("N", trt))  %>% group_by(site_code, category) %>%
  summarise(AGN=mean(AGN),
            pct_N=mean(pct_N),
            proportion_par=mean(proportion_par),
            precip.mean=mean(precip.mean),
            PET.mean=mean(PET.mean))
allN$N1="N"

allnotN <- all %>% filter(!grepl("N", trt)) %>% group_by(site_code, category) %>%
  summarise(AGN=mean(AGN),
            pct_N=mean(pct_N),
            proportion_par=mean(proportion_par))
allnotN$notN="Not N"


test <- (allN$AGN-allnotN$AGN)/allN$AGN
test2 <-(allN$pct_N-allnotN$pct_N)/allN$pct_N


allN$deltaN <- test
allN$deltapctN <- test2
allN$Wetness <- allN$precip.mean/allN$PET.mean
  


ggplot(allN, aes((precip.mean/PET.mean), deltaN, color=category))+ geom_point()+
  scale_y_continuous(limits=c(-1,1)) +
  geom_smooth(method="lm", se=F)

ggplot(allN, aes(proportion_par, deltaN, color=Wetness))+
  scale_y_continuous(limits=c(-1,1)) +
  geom_point(size=3)+
  scale_color_viridis(begin=1, end=0)+
  theme_bw()

ggplot(allN, aes(proportion_par, deltapctN, color=Wetness))+
  geom_point(size=3)+
  scale_y_continuous(limits=c(-1,1)) +
  scale_color_viridis(begin=1, end=0)+
  theme_bw()



  
