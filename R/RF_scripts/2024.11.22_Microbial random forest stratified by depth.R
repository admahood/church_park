library(pacman)
p_load(tidyverse, here, reshape2, randomForest, paletteer, iml, patchwork, ggh4x, tidytext, cowplot)

iData <- read.csv(here('Data', 'cp_trait_soil_data.csv'))


# RFs with combined soil depths
uniqueData <- iData[which(duplicated(iData$simpson_fungi_5)==F),]

shallowPreds <- uniqueData %>%
  select(., c(biochar, mulch, vwc, pH, total_veg_cover,
              twi, Aspect, bare, mulch, biochar, 
              percent_clay_0to5, sand_0_5, 
              total_c_0_5, total_n_0_5, 
              DOC_0_5, TDN_0_5, 
              ammonium_0_5, nitrate_0_5, potassium_0_5, phosphate_0_5, 
              cations_0_5, cations_5_15,
              CN_0_5, DIN_0_5, DON_0_5))
shallowDepVar <- uniqueData %>%
  select(., c(nitrifiers,EMF_0_5,simpson_bacteria_5, simpson_fungi_5))


deepPreds <- uniqueData %>%
  select(., c(biochar, mulch, vwc, pH, total_veg_cover,
              twi, Aspect, bare, mulch, biochar, 
              percent_clay_5to15, percent_sand_5to15, 
              total_c_5_15, total_n_5_15,
              DOC_5_15, TDN_5_15, 
              ammonium_5_15, nitrate_5_15, potassium_5_15, phosphate_5_15, 
              cations_5_15))
deepDepVar <- uniqueData %>%
  select(., c(nitrifiers, EMF_5_15,simpson_bacteria_15, simpson_fungi_15))

rfFunc <- function(depVar, depVars, predictors){
  #depVar <- 'nitrifiers'
  #depVars <- deepDepVar
  #predictors <- deepPreds
  formula <- as.formula(paste(depVar, '~.'))
  
  rfData <- data.frame(depVars %>% select(., depVar),
                       predictors) %>%
    na.omit
  
  randomForest(formula, data = rfData, ntrees = 5000, mtry = sqrt(ncol(rfData)))
  
}


sNitrifiers <- rfFunc('nitrifiers', shallowDepVar, shallowPreds)
sEMF <- rfFunc('EMF_0_5', shallowDepVar, shallowPreds)
sFungalDiv <- rfFunc('simpson_fungi_5', shallowDepVar, shallowPreds)
sBactDiv <- rfFunc('simpson_bacteria_5', shallowDepVar, shallowPreds)

dNitrifiers <- rfFunc('nitrifiers', deepDepVar, deepPreds)
dEMF <- rfFunc('EMF_5_15', deepDepVar, deepPreds)
dFungalDiv <- rfFunc('simpson_fungi_15', deepDepVar, deepPreds)
dBactDiv <- rfFunc('simpson_bacteria_15', deepDepVar, deepPreds)


# model performance

cor(uniqueData$nitrifiers, predict(sNitrifiers, uniqueData))^2
cor(uniqueData$EMF_0_5, predict(sEMF, uniqueData))^2
cor(uniqueData$simpson_fungi_5, predict(sFungalDiv, uniqueData))^2
cor(uniqueData$simpson_bacteria_5, predict(sBactDiv, uniqueData))^2

cor(uniqueData$nitrifiers, predict(sNitrifiers, uniqueData))^2
cor(uniqueData$EMF_5_15, predict(sEMF, uniqueData))^2
cor(uniqueData$simpson_fungi_15, predict(sFungalDiv, uniqueData))^2
cor(uniqueData$simpson_bacteria_15, predict(sBactDiv, uniqueData))^2



# VARIABLE IMPORTANCE PLOTS

varImpFunc <- function(model, depVarName, depth){
    #model <- sNitrifiers
    #depVarName <- 'Nitrifiers'
    #depth <- '0-5 cm'
    model$importance %>%
      data.frame() %>%
      rownames_to_column() %>%
      subset(., rowname != 'species_full') %>%
      mutate(RelImportance = round(IncNodePurity/max(IncNodePurity), 2), 
             DepVar = depVarName,
             soilDepth = depth) %>%
      arrange(RelImportance) %>%
      tail(n=10)
}

sVarImp <- rbind(
  varImpFunc(sNitrifiers, 'Nitrifiers', '0-5 cm'),
  varImpFunc(sEMF, 'EMF', '0-5 cm'),
  varImpFunc(sFungalDiv, 'Fungal diversity', '0-5 cm'),
  varImpFunc(sBactDiv, 'Bacterial diversity', '0-5 cm')
)  

dVarImp <- rbind(
  varImpFunc(dNitrifiers, 'Nitrifiers', '5-15 cm'),
  varImpFunc(dEMF, 'EMF', '5-15 cm'),
  varImpFunc(dFungalDiv, 'Fungal diversity', '5-15 cm'),
  varImpFunc(dBactDiv, 'Bacterial diversity', '5-15 cm')
)  

# shallow plot
ggplot(rbind(sVarImp, dVarImp), aes(x=DepVar, y=reorder_within(rowname, RelImportance, list(DepVar, soilDepth)), fill = RelImportance, label = RelImportance)) +
  geom_tile() + theme_bw(base_size = 20) +
  geom_text(color = 'white',fontface='bold', size = 5) +
  paletteer::scale_fill_paletteer_c("ggthemes::Classic Orange-Blue") +
  facet_grid2(soilDepth~DepVar, scales = 'free', independent = 'all') +
  scale_y_reordered() +
  labs(fill = 'Normalized\nvariable importance') +
  theme(legend.position = 'bottom',
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.width=unit(2.5,"cm"))

ggsave(here('Figures', 'FigSx_by depth microbe RF.png'), height = 25, width = 45, dpi = 300, units = 'cm')
