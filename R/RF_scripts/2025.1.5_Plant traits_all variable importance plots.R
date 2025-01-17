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





traitsRFdata <- data.frame(
  combTargs %>% select(., c(1:4)),
  sumarSoils) %>%
  mutate(species_full = case_when(
    str_detect(species_full, regex('Oreo')) ~ 'O. parryi',
    str_detect(species_full, regex('Vacc')) ~ 'V. scoparium'
  )) %>%
  rename(., LDMC = ldmc,
         SLA = sla,
         Height = height_cm) #%>%
  #mutate_at(2:26, scale)





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
  group_by(Predictor, variable) %>%
  summarise(varImp = mean(value)) %>%
  arrange(desc(varImp)) %>%
  data.frame



# partial dependence plots

#species <- "V. scoparium"
#trait <- 'SLA'
#model <- sla_Vaccinium
#predictor = 'NH4'

pdpFunc <- function(species, trait, model, predictor){
  
  data <- subFunc(species, trait)
  Pred <- Predictor$new(model, data)
  pdp <- FeatureEffect$new(Pred, feature = predictor, method = 'pdp')
  pdp$results %>%
    mutate(Species = species,
           Trait = trait,
           Predictor = predictor) %>%
    setNames(., nm = c('x', 'y', 'type', 'Species', 'Trait', 'Predictor')) %>%
    mutate(normY = 100*((y-mean(y))/mean(y)))
}

predictors <- names(traitsRFdata[5:26])

pdpFuncMultiple <- function(species, trait, model, predictors){
  
  # Initialize an empty list to store results for each predictor
  all_results <- list()
  
  # Iterate through the vector of predictors
  for (predictor in predictors) {
    
    # Call the original pdpFunc for each predictor
    pdp_result <- pdpFunc(species, trait, model, predictor)
    
    # Append the result to the list
    all_results[[predictor]] <- pdp_result
  }
  
  # Combine all results into a single data frame
  combined_results <- bind_rows(all_results)
  
  return(combined_results)
}

orderVarImp <- allModels %>%
  melt(id = c('Predictor', 'Trait', 'Group')) %>%
  group_by(Group, Predictor, variable) %>% 
  summarise(meanVarImp = mean(value)) %>%
  left_join(., allModels %>% melt(id = c('Predictor', 'Trait', 'Group'))) %>%
  mutate(value = round(value, 2)) %>%
  setNames(., nm = c('Group', 'Predictor', 'Species', 'meanVarImp', 'Trait', 'value'))

allPDP <- rbind(
pdpFuncMultiple('V. scoparium', 'SLA', sla_Vaccinium, predictors),
pdpFuncMultiple('V. scoparium', 'LDMC', ldmc_Vaccinium, predictors),
pdpFuncMultiple('V. scoparium', 'Height', height_Vaccinium, predictors),

pdpFuncMultiple('O. parryi', 'SLA', sla_Oreochrysum, predictors),
pdpFuncMultiple('O. parryi', 'LDMC', ldmc_Oreochrysum, predictors),
pdpFuncMultiple('O. parryi', 'Height', height_Oreochrysum, predictors)
)

jnData <- left_join(orderVarImp, allPDP) 

cols <- c('CN', 'DIN', 'DON', 'NH4', 'NO3', 'Bare.soil', 'BA.div', 'EMF', 'Nitrifiers', 'Clay', 'TWI')

jnData %>% filter(., Predictor %in% cols) %>%
  subset(x > -0.35) %>%
  ggplot(., aes(x=x, y=normY, color = Species, alpha = value)) + theme_bw(base_size = 25) +
    geom_hline(yintercept = 0, color = 'black') +  
    geom_line(size = 4) +
    facet_nested(Trait ~ Group + Predictor, scales = 'free', independent = 'x') +
    scale_x_continuous(n.breaks = 4) +
    scale_color_manual(values = c('#008080FF', '#CA562CFF')) +
    scale_alpha_continuous(range = c(-0.25, 1)) +
    labs(y = 'Standardized change (%)', color = 'Species', alpha = 'Standardized\nvariable importance') +
    theme(#strip.text.y.right = element_text(angle = 360),
          axis.title.x = element_blank(),
          legend.position = 'bottom',
          axis.text = element_text(size = 13))

ggsave(here('Figures', 'Fig 4_RF pdp.png'), height = 30, width = 60, dpi = 300, units = 'cm')

