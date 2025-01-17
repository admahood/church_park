library(pacman)
p_load(tidyverse, here, reshape2, randomForest, paletteer, iml, patchwork, ggh4x, tidytext, cowplot, corrplot, magrittr, pdp)

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
              CN_0_5, CN_5_15, DIN_0_5, DON_0_5))

sumarSoils <- combTargs %>%
  transmute(BA.div = (simpson_bacteria_15 + simpson_bacteria_5)/2,
            Fung.div = (simpson_fungi_5 + simpson_fungi_15)/2,
            EMF = (EMF_0_5 + EMF_5_15)/2,
            Nitrifiers = nitrifiers,
            TWI = twi,
            VWC = vwc,
            Vegetation = total_veg_cover,
            Mulch = mulch,
            Biochar = biochar,
            Bare.soil = bare,
            Clay = (percent_clay_0to5 + percent_clay_5to15)/2,
            Sand = (sand_0_5 + percent_sand_5to15)/2,
            #$totalC = (total_c_0_5 + total_c_5_15)/2,
            #totalN = (total_n_0_5 + total_n_5_15)/2,
            DOC = (DOC_0_5 + DOC_5_15)/2,
            TDN = (TDN_0_5 + TDN_5_15)/2,
            NH4 = (ammonium_0_5 + ammonium_5_15)/2,
            NO3 = (nitrate_0_5 + nitrate_5_15)/2,
            #K = (potassium_0_5 + potassium_5_15)/2,
            PO4 = (phosphate_0_5 + phosphate_5_15)/2,
            pH = pH,
            Cations = (cations_0_5 + cations_5_15)/2,
            CN = (CN_0_5 + CN_5_15)/2,
            DIN = DIN_0_5,
            DON = DON_0_5)


ggplot(iData, aes(x=treatment, y=pH, fill = treatment)) + geom_boxplot()


traitsRFdata <- data.frame(
  combTargs %>% select(., c(1:4)),
  sumarSoils) %>%
  mutate(species_full = case_when(
    str_detect(species_full, regex('Oreo')) ~ 'O. parryi',
    str_detect(species_full, regex('Vacc')) ~ 'V. scoparium'
  )) %>%
  rename(., LDMC = ldmc,
         SLA = sla,
         Height = height_cm)

subFunc <- function(species, trait){
  #species = 'V. scoparium'
  #trait = 'ldmc'
  
  predictors <- traitsRFdata %>% subset(., species_full == species) %>%
    select(., -c(1:4))
  
  depVar <- traitsRFdata %>%
    subset(., species_full == species) %>%
    select(., trait)
  
  out <- cbind(depVar, predictors) %>% 
    na.omit
  
  out[-which(duplicated(out$DON)),]
  
}    



# SLA
sla_Vaccinium <- randomForest(SLA~., data = subFunc('V. scoparium', 'SLA'), mtry = 4, ntree = 1000, nodesize = 5, maxnodes = 5)
#sla_Carex <- randomForest(sla~., data = subFunc('Carex sp', 'sla'), mtry = 4, ntree = 1000, nodesize = 100, maxnodes = 50) 
sla_Oreochrysum <- randomForest(SLA~., data = subFunc('O. parryi', 'SLA'), mtry = 4, ntree = 1000, nodesize = 5, maxnodes = 5) 

# LDMC
ldmc_Vaccinium <- randomForest(LDMC~., data = subFunc('V. scoparium', 'LDMC'), mtry = 4, ntree = 1000, nodesize = 5, maxnodes = 5)
#ldmc_Carex <- randomForest(ldmc~., data = subFunc('Carex sp', 'ldmc'), mtry = 4, ntree = 1000, nodesize = 100, maxnodes = 50) 
ldmc_Oreochrysum <- randomForest(LDMC~., data = subFunc('O. parryi', 'LDMC'), mtry = 4, ntree = 1000, nodesize = 5, maxnodes = 5) 

# Height
height_Vaccinium <- randomForest(Height~., data = subFunc('V. scoparium', 'Height'), mtry = 4, ntree = 1000, nodesize = 5, maxnodes = 5)
#height_Carex <- randomForest(height_cm~., data = subFunc('Carex sp', 'height_cm'), mtry = 4, ntree = 1000, nodesize = 100, maxnodes = 50) 
height_Oreochrysum <- randomForest(Height~., data = subFunc('O. parryi', 'Height'), mtry = 4, ntree = 1000, nodesize = 5, maxnodes = 5) 




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
    str_detect(Predictor, regex('Vegetation|Biochar|Mulch|Bare')) ~ 'Cover',
    str_detect(Predictor, regex('.div|EMF|Nitrifiers')) ~ 'Microbes',
    str_detect(Predictor, regex('VWC|TWI|Sand|Clay|Aspect')) ~ 'Moisture',
    str_detect(Predictor, regex('DOC|TDN|NH4|NO3|CN|DIN|DON|PO4|pH|Cations')) ~ 'Chemistry'))
    #str_detect(Predictor, regex('PO4|pH|Cations')) ~ 'Chem'))





# stats for results
allModels %>%
  melt(id = c('Predictor', 'Trait', 'Group')) %>%
  group_by(Predictor) %>%
  summarise(varImp = mean(value)) %>%
  arrange(desc(varImp)) %>%
  data.frame



# partial dependence plots

pdpFunc <- function(species, trait, model, predictor){
  #species <- "V. scoparium"
  #trait <- 'SLA'
  #model <- sla_Vaccinium
  #predictor = 'NH4'
  data <- subFunc(species, trait)
  Pred <- Predictor$new(model, data)
  pdp <- FeatureEffect$new(Pred, feature = predictor, method = 'pdp')
  pdp$results %>%
    mutate(Species = species,
           Trait = trait,
           Predictor = predictor) %>%
    setNames(., nm = c('x', 'y', 'type', 'Species', 'Trait', 'Predictor')) %>%
    mutate(normY = y/max(y))
}
  
allPDP <- rbind(
  pdpFunc('V. scoparium', 'SLA', sla_Vaccinium, 'NH4'),
  pdpFunc('V. scoparium', 'SLA', sla_Vaccinium, 'NO3'),
  pdpFunc('V. scoparium', 'SLA', sla_Vaccinium, 'EMF'),
  pdpFunc('V. scoparium', 'SLA', sla_Vaccinium, 'Nitrifiers'),
  
  pdpFunc('V. scoparium', 'LDMC', ldmc_Vaccinium, 'NH4'),
  pdpFunc('V. scoparium', 'LDMC', ldmc_Vaccinium, 'NO3'),
  pdpFunc('V. scoparium', 'LDMC', ldmc_Vaccinium, 'EMF'),
  pdpFunc('V. scoparium', 'LDMC', ldmc_Vaccinium, 'Nitrifiers'),
  
  pdpFunc('V. scoparium', 'Height', height_Vaccinium, 'NH4'),
  pdpFunc('V. scoparium', 'Height', height_Vaccinium, 'NO3'),
  pdpFunc('V. scoparium', 'Height', height_Vaccinium, 'EMF'),
  pdpFunc('V. scoparium', 'Height', height_Vaccinium, 'Nitrifiers'),
  
  
  pdpFunc('O. parryi', 'SLA', sla_Oreochrysum, 'NH4'),
  pdpFunc('O. parryi', 'SLA', sla_Oreochrysum, 'NO3'),
  pdpFunc('O. parryi', 'SLA', sla_Oreochrysum, 'EMF'),
  pdpFunc('O. parryi', 'SLA', sla_Oreochrysum, 'Nitrifiers'),
  
  pdpFunc('O. parryi', 'LDMC', ldmc_Oreochrysum, 'NH4'),
  pdpFunc('O. parryi', 'LDMC', ldmc_Oreochrysum, 'NO3'),
  pdpFunc('O. parryi', 'LDMC', ldmc_Oreochrysum, 'EMF'),
  pdpFunc('O. parryi', 'LDMC', ldmc_Oreochrysum, 'Nitrifiers'),
  
  pdpFunc('O. parryi', 'Height', height_Oreochrysum, 'NH4'),
  pdpFunc('O. parryi', 'Height', height_Oreochrysum, 'NO3'),
  pdpFunc('O. parryi', 'Height', height_Oreochrysum, 'EMF'),
  pdpFunc('O. parryi', 'Height', height_Oreochrysum, 'Nitrifiers')
  
)

orderVarImp <- allModels %>%
  melt(id = c('Predictor', 'Trait', 'Group')) %>%
  group_by(Group, Predictor, variable) %>% 
  summarise(meanVarImp = mean(value)) %>%
  left_join(., allModels %>% melt(id = c('Predictor', 'Trait', 'Group'))) %>%
  mutate(value = round(value, 2)) 
  

rfVarImp <- ggplot(orderVarImp, aes(x=Trait, y=reorder_within(x=Predictor, by=meanVarImp, within=list(variable, Group)), fill = value, label = value)) +
  scale_y_reordered() +
  theme_bw(base_size = 20) +
  labs(fill = 'Standardized variable\nimportance') +
  geom_tile(color = 'white', size = 0.5) +
  paletteer::scale_fill_paletteer_c("ggthemes::Classic Orange-Blue") +
  geom_text(color = 'white',fontface='bold', size = 5) +
  facet_grid2(Group~variable, scales = 'free_y', independent = 'y') +
  force_panelsizes(rows = c(10*4.5454, 4*4.5454, 4*4.5454, 4*4.5454)) +
  theme(axis.title = element_blank(),
        #axis.text.x = element_text(angle = 35, vjust = 1, hjust=1),
        legend.position="bottom",
        legend.key.width=unit(2,"cm"))
  

pdp <- ggplot(allPDP, aes(x=x, y=normY, color = Species)) + theme_bw(base_size = 20) +
  geom_line(size = 2) +
  labs(y = 'Standardized partial dependence') +
  scale_color_brewer(palette = 'Set1') +
  scale_x_continuous(n.breaks = 4, breaks = waiver()) +
  facet_grid2(Trait~Predictor, scales = 'free_x') +
  theme(axis.title.x = element_blank(),
        legend.position = 'bottom')


plot_grid(rfVarImp, pdp, align = 'vh', axis = 'tb', rel_widths = c(0.55, 0.45))

ggsave(here('Figures', 'v2_Fig4_RF on traits.png'), plot = last_plot(), height = 30, width = 50, unit = 'cm')

