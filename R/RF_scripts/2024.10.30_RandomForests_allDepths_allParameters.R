library(pacman)
p_load(tidyverse, here, reshape2, randomForest, paletteer, iml, patchwork, ggh4x, tidytext, cowplot, corrplot, Metrics)

iData <- read.csv(here('Data', 'cp_trait_soil_data.csv'))


# RFs with combined soil depths
combTargs <- iData %>%
  select(., c(species_full, ldmc, sla, height_cm, biochar, mulch, vwc, pH, total_veg_cover,
              twi, Aspect, bare, mulch, biochar, nitrifiers,
              simpson_bacteria_15, simpson_bacteria_5, simpson_fungi_5, simpson_fungi_15, EMF_0_5, EMF_5_15,  
              percent_clay_0to5, percent_clay_5to15, sand_0_5, percent_sand_5to15, 
              total_c_0_5, total_c_5_15, total_n_0_5, total_n_5_15,
              DOC_0_5, DOC_5_15, TDN_0_5, TDN_5_15, 
              ammonium_0_5, ammonium_5_15, nitrate_0_5, nitrate_5_15, potassium_0_5, potassium_5_15, phosphate_0_5, phosphate_5_15, 
              cations_0_5, cations_5_15,
              CN_0_5, DIN_0_5, DON_0_5)) #CN_5_15, 


traitsRFdata <- combTargs %>% 
  data.frame() %>%
  mutate(species_full = case_when(
    str_detect(species_full, regex('Oreo')) ~ 'O. parryi',
    str_detect(species_full, regex('Vacc')) ~ 'V. scoparium'
  )) %>%
  rename(., LDMC = ldmc,
         SLA = sla,
         Height = height_cm)

subFunc <- function(species, trait){
  species = 'V. scoparium'
  trait = 'SLA'
  
  predictors <- traitsRFdata %>% subset(., species_full == species) %>%
    select(., -c(1:4))
  
  depVar <- traitsRFdata %>%
    subset(., species_full == species) %>%
    select(., trait)
  
  cbind(depVar, predictors) %>% na.omit
}    

subFunc('V. scoparium', SLA)

# SLA
sla_Vaccinium <- randomForest(SLA~., data = subFunc('V. scoparium', 'SLA'), mtry = 4, ntree = 1000, nodesize = 50, maxnodes = 50)
#sla_Carex <- randomForest(sla~., data = subFunc('Carex sp', 'sla'), mtry = 4, ntree = 1000, nodesize = 100, maxnodes = 50) 
sla_Oreochrysum <- randomForest(SLA~., data = subFunc('O. parryi', 'SLA'), mtry = 4, ntree = 1000, nodesize = 100, maxnodes = 50) 

# LDMC
ldmc_Vaccinium <- randomForest(LDMC~., data = subFunc('V. scoparium', 'LDMC'), mtry = 4, ntree = 1000, nodesize = 50, maxnodes = 50)
#ldmc_Carex <- randomForest(ldmc~., data = subFunc('Carex sp', 'ldmc'), mtry = 4, ntree = 1000, nodesize = 100, maxnodes = 50) 
ldmc_Oreochrysum <- randomForest(LDMC~., data = subFunc('O. parryi', 'LDMC'), mtry = 4, ntree = 1000, nodesize = 100, maxnodes = 50) 

# Height
height_Vaccinium <- randomForest(Height~., data = subFunc('V. scoparium', 'Height'), mtry = 4, ntree = 1000, nodesize = 50, maxnodes = 50)
#height_Carex <- randomForest(height_cm~., data = subFunc('Carex sp', 'height_cm'), mtry = 4, ntree = 1000, nodesize = 100, maxnodes = 50) 
height_Oreochrysum <- randomForest(Height~., data = subFunc('O. parryi', 'Height'), mtry = 4, ntree = 1000, nodesize = 100, maxnodes = 50) 


# regressions to assess model performance
predData <- iData %>%
  mutate(slaVacPred = predict(sla_Vaccinium, iData),
         slaOreoPred = predict(sla_Oreochrysum, iData),
         ldmcVacPred = predict(ldmc_Vaccinium, iData),
         ldmcOreoPred = predict(ldmc_Oreochrysum, iData),
         heightVacPred = predict(height_Vaccinium, iData),
         heightOreoPred = predict(height_Oreochrysum, iData))
        
vaccSla <- lm(sla ~ slaVacPred, data = predData %>% subset(species_full == 'Vaccinium sp')) %>% summary                       
oreoSla <- lm(sla ~ slaOreoPred, data = predData %>% subset(species_full == 'Oreochrysum parryi')) %>% summary                                

vaccLdmc <- lm(sla ~ ldmcVacPred, data = predData %>% subset(species_full == 'Vaccinium sp')) %>% summary                       
oreoLdmc <- lm(sla ~ ldmcOreoPred, data = predData %>% subset(species_full == 'Oreochrysum parryi')) %>% summary                                

vaccHeight <- lm(height_cm ~ heightVacPred, data = predData %>% subset(species_full == 'Vaccinium sp')) %>% summary                       
oreoHeight <- lm(height_cm ~ heightOreoPred, data = predData %>% subset(species_full == 'Oreochrysum parryi')) %>% summary                                


vaccSub <- predData %>%
  subset(species_full == 'Vaccinium sp')

attach(vaccSub)
bias(sla, slaVacPred)
bias(ldmc, ldmcVacPred)
bias(height_cm, heightVacPred)

rmse(sla, slaVacPred)
rmse(ldmc, ldmcVacPred)
rmse(height_cm, heightVacPred)
detach(vaccSub)

oreoSub <- predData %>%
  subset(species_full == 'Oreochrysum parryi')
attach(oreoSub)
bias(sla, slaOreoPred)
bias(ldmc, ldmcOreoPred)
bias(height_cm, heightOreoPred)

rmse(sla, slaOreoPred)
rmse(ldmc, ldmcOreoPred)
rmse(height_cm, heightOreoPred)
detach(oreoSub)

viTraits <- function(model, speciesName){
  #model = sla_Oreochrysum
  #speciesName = 'Oreochrysum'
  #trait = 'SLA'
  df = model$importance %>% 
    as.data.frame %>%
    rownames_to_column
  df %>% 
    mutate(normVarImp = IncNodePurity/max(IncNodePurity)) %>%
    select(., -IncNodePurity) %>%
    setNames(., nm = c('Predictor', speciesName))
}


# All SLA data variable importances
allSLA <- left_join(viTraits(sla_Vaccinium, 'V. scoparium'), viTraits(sla_Oreochrysum, 'O. parryi')) %>%
  mutate(Trait = 'SLA')

allLDMC <- left_join(viTraits(ldmc_Vaccinium, 'V. scoparium'), viTraits(ldmc_Oreochrysum, 'O. parryi')) %>%
  mutate(Trait = 'LDMC')

allHeight <- left_join(viTraits(height_Vaccinium, 'V. scoparium'), viTraits(height_Oreochrysum, 'O. parryi')) %>%
  mutate(Trait = 'Height')

allModels <- bind_rows(allSLA, allLDMC, allHeight) %>%
  mutate(Group = case_when(
    str_detect(Predictor, regex('veg|biochar|mulch|bare')) ~ 'Cover',
    str_detect(Predictor, regex('simpson|richness|EMF|nitrifiers')) ~ 'Microbes',
    str_detect(Predictor, regex('vwc|twi|sand|clay|Aspect')) ~ 'Moisture',
    str_detect(Predictor, regex('DOC|TDN|nitrate|ammonium|CN|DIN|DON|total')) ~ 'C & N',
    str_detect(Predictor, regex('phosphate|potassium|pH|cations')) ~ 'Other Chem'))




orderVarImp <- allModels %>%
  melt(id = c('Predictor', 'Trait', 'Group')) %>%
  group_by(Group, Predictor, variable) %>% 
  summarise(meanVarImp = mean(value)) %>%
  left_join(., allModels %>% melt(id = c('Predictor', 'Trait', 'Group'))) %>%
  mutate(value = round(value, 2)) 

ggplot(orderVarImp, aes(x=Trait, y=reorder_within(x=Predictor, by=meanVarImp, within=list(variable, Group)), fill = value, label = value)) +
  scale_y_reordered() +
  theme_bw(base_size = 20) +
  labs(fill = 'Relative variable\nimportance') +
  geom_tile(color = 'white', size = 0.5) +
  paletteer::scale_fill_paletteer_c("ggthemes::Classic Orange-Blue") +
  geom_text(color = 'white',fontface='bold', size = 5) +
  facet_grid2(Group~variable, scales = 'free_y', independent = 'y') +
  force_panelsizes(rows = c(16*0.0213, 4*0.0213, 7*0.0213, 7*0.0213, 7*0.0213)) +
  theme(axis.title = element_blank(),
        #axis.text.x = element_text(angle = 35, vjust = 1, hjust=1),
        legend.position="bottom",
        legend.key.width=unit(2,"cm"))

ggsave(here('Figures', 'FigSx_RF var imp all depths.png'), height = 40, width = 35, dpi = 300, unit = 'cm')

