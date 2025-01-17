library(pacman)
p_load(tidyverse, here, reshape2, randomForest, paletteer, iml, patchwork, ggh4x, tidytext, cowplot, pdp)

iData <- read.csv(here('Data', 'cp_trait_soil_data.csv'))


# RFs with combined soil depths
combTargs <- iData %>%
  select(., c(species_full, sla, ldmc, biochar, mulch, vwc, pH, total_veg_cover,
              twi, Aspect, bare, mulch, biochar, nitrifiers,
              EMF_0_5, EMF_5_15,   
              simpson_bacteria_15, simpson_bacteria_5, simpson_fungi_5, simpson_fungi_15,
              percent_clay_0to5, percent_clay_5to15, sand_0_5, percent_sand_5to15, 
              total_c_0_5, total_c_5_15, total_n_0_5, total_n_5_15,
              DOC_0_5, DOC_5_15, TDN_0_5, TDN_5_15, 
              ammonium_0_5, ammonium_5_15, nitrate_0_5, nitrate_5_15, potassium_0_5, potassium_5_15, phosphate_0_5, phosphate_5_15, 
              cations_0_5, cations_5_15,
              CN_0_5, CN_5_15, DIN_0_5, DON_0_5))


sumarSoils <- combTargs %>%
  transmute(BA.div = (simpson_bacteria_15 + simpson_bacteria_5)/2,
            Fungal.div = (simpson_fungi_5 + simpson_fungi_15)/2,
            EMF = (EMF_0_5 + EMF_5_15)/2,
            Nitrifiers = nitrifiers,
            TWI = twi,
            VWC = vwc,
            SLA = sla,
            LDMC = ldmc,
            Vegetation = total_veg_cover,
            Mulch = mulch,
            Biochar = biochar,
            Bare.soil = bare,
            Clay = (percent_clay_0to5 + percent_clay_5to15)/2,
            Sand = (sand_0_5 + percent_sand_5to15)/2,
            TotalC = (total_c_0_5 + total_c_5_15)/2,
            TotalN = (total_n_0_5 + total_n_5_15)/2,
            DOC = (DOC_0_5 + DOC_5_15)/2,
            TDN = (TDN_0_5 + TDN_5_15)/2,
            NH4 = (ammonium_0_5 + ammonium_5_15)/2,
            NO3 = (nitrate_0_5 + nitrate_5_15)/2,
            K = (potassium_0_5 + potassium_5_15)/2,
            PO4 = (phosphate_0_5 + phosphate_5_15)/2,
            pH = pH,
            Cations = (cations_0_5 + cations_5_15)/2,
            #CN = (CN_0_5 + CN_5_15)/2,
            #DIN = DIN_0_5,
            DON = DON_0_5)

microRFdata <- data.frame(
  combTargs %>% select(., c(1)),
  sumarSoils) %>% 
  subset(species_full == 'Carex sp')

uniqueRFdata <- microRFdata[-which(duplicated(microRFdata$BA.div) == T),]
ncol(uniqueRFdata)



# RUN MODELS --------------------------------------------------------------
emfRF <- randomForest(EMF ~., data = microRFdata %>% select(., -c(BA.div, Fungal.div, Nitrifiers)) %>% na.omit,
                      ntree = 5000, mtry = 5, nodesize = 100, maxnodes = 10)

nitriRF <- randomForest(Nitrifiers ~., data = microRFdata %>% select(., -c(BA.div, Fungal.div, EMF)) %>% na.omit,
                        ntree = 5000, mtry = 5, nodesize = 100, maxnodes = 10)

fungDivRF <- randomForest(Fungal.div ~., data = microRFdata %>% select(., -c(BA.div, EMF, Nitrifiers)) %>% na.omit,
                        ntree = 5000, mtry = 5, nodesize = 100, maxnodes = 10)

bactDivRF <- randomForest(BA.div ~., data = microRFdata %>% select(., -c(Fungal.div, EMF, Nitrifiers)) %>% na.omit,
                        ntree = 5000, mtry = 5, nodesize = 100, maxnodes = 10)



# PARTIAL DEPENDENCE PLOTS ------------------------------------------------

partialPlot(fungDivRF, uniqueRFdata, PO4)
partialPlot(emfRF, uniqueRFdata, NO3)


# MODEL PERFORMANCE -------------------------------------------------------
cor(uniqueRFdata$Nitrifiers, predict(nitriRF, uniqueRFdata))^2
cor(uniqueRFdata$EMF, predict(emfRF, uniqueRFdata))^2
cor(uniqueRFdata$Fungal.div, predict(fungDivRF, uniqueRFdata))^2
cor(uniqueRFdata$Bacterial.div, predict(bactDivRF, uniqueRFdata))^2



emfVarImp <- emfRF$importance %>%
  data.frame() %>%
  rownames_to_column() %>%
  subset(., rowname != 'species_full') %>%
  mutate(RelImportance = IncNodePurity/max(IncNodePurity),
         DepVar = 'EMF',
         RelImportance = round(RelImportance, 2)) %>%
  arrange(RelImportance) %>%
  tail(n=10)

nitVarImp <- nitriRF$importance %>%
  data.frame() %>%
  rownames_to_column() %>%
  subset(., rowname != 'species_full') %>%
  mutate(RelImportance = IncNodePurity/max(IncNodePurity),
         DepVar = 'Nitrifiers',
         RelImportance = round(RelImportance, 2)) %>%
  arrange(RelImportance) %>%
  tail(n=10)

bacDivVarImp <- bactDivRF$importance %>%
  data.frame() %>%
  rownames_to_column() %>%
  subset(., rowname != 'species_full') %>%
  mutate(RelImportance = IncNodePurity/max(IncNodePurity),
         DepVar = 'BA.div',
         RelImportance = round(RelImportance, 2)) %>%
  arrange(RelImportance) %>%
  tail(n=10)

fungDivVarImp <- fungDivRF$importance %>%
  data.frame() %>%
  rownames_to_column() %>%
  subset(., rowname != 'species_full') %>%
  mutate(RelImportance = IncNodePurity/max(IncNodePurity),
         DepVar = 'Fungal.div',
         RelImportance = round(RelImportance, 2)) %>%
  arrange(RelImportance) %>%
  tail(n=10)


plotData <- rbind(emfVarImp, nitVarImp, bacDivVarImp, fungDivVarImp) %>%
  mutate(DepVar = factor(DepVar, levels = c('EMF', 'Nitrifiers', 'BA.div', 'Fungal.div'), ordered=T))

ggplot(plotData, aes(x=DepVar, y=reorder_within(rowname, RelImportance, DepVar), fill = RelImportance, label = RelImportance)) +
  geom_tile() + theme_bw(base_size = 20) +
  geom_text(color = 'white',fontface='bold', size = 5) +
  paletteer::scale_fill_paletteer_c("ggthemes::Classic Orange-Blue") +
  facet_wrap(~DepVar, ncol = 2, scales = 'free') +
  scale_y_reordered() +
  labs(fill = 'Standardized\nvariable importance') +
  theme(legend.position = 'bottom',
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.width=unit(2.5,"cm"))

ggsave(here('Figures', 'Fig5_Microbe RFs.png'), height = 25, width = 25, dpi = 300, units = 'cm')


# emfPlot <- ggplot(emfVarImp, aes(x=DepVar, y=reorder(rowname, RelImportance), fill = RelImportance, label = RelImportance)) +
#   geom_tile() + theme_bw(base_size = 20) +
#   geom_text(color = 'white',fontface='bold', size = 5) +
#   paletteer::scale_fill_paletteer_c("ggthemes::Classic Orange-Blue") +
#   facet_wrap(~DepVar) +
#   theme(legend.position = 'none',
#         axis.title = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())
#   
# nitPlot <- ggplot(nitVarImp, aes(x=DepVar, y=reorder(rowname, RelImportance), fill = RelImportance, label = RelImportance)) +
#   geom_tile() + theme_bw(base_size = 20) +
#   geom_text(color = 'white',fontface='bold', size = 5) +
#   paletteer::scale_fill_paletteer_c("ggthemes::Classic Orange-Blue") +
#   scale_y_discrete(position = 'right') +
#   facet_wrap(~DepVar) +
#   theme(legend.position = 'none',
#         axis.title = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())
# 
# 
# fungDivPlot <- ggplot(fungDivVarImp, aes(x=DepVar, y=reorder(rowname, RelImportance), fill = RelImportance, label = RelImportance)) +
#   geom_tile() + theme_bw(base_size = 20) +
#   geom_text(color = 'white',fontface='bold', size = 5) +
#   paletteer::scale_fill_paletteer_c("ggthemes::Classic Orange-Blue") +
#   facet_wrap(~DepVar) +
#   theme(legend.position = 'none',
#         axis.title = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())
# 
# bactDivPlot <- ggplot(bacDivVarImp, aes(x=DepVar, y=reorder(rowname, RelImportance), fill = RelImportance, label = RelImportance)) +
#   geom_tile() + theme_bw(base_size = 20) +
#   geom_text(color = 'white',fontface='bold', size = 5) +
#   paletteer::scale_fill_paletteer_c("ggthemes::Classic Orange-Blue") +
#   scale_y_discrete(position = 'right') +
#   facet_wrap(~DepVar) +
#   theme(legend.position = 'none',
#         axis.title = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())

#plot_grid(emfPlot, nitPlot, fungDivPlot, bactDivPlot, ncol = 2, align = 'vh', axis = 'r', rel_widths = c(0.55, 0.45))




# ANOVA of EMF and nitrifiers ---------------------------------------------

microAnova <- iData %>%
  select(., c(treatment, EMF_0_5, EMF_5_15, nitrifiers)) %>%
  mutate(EMF = (EMF_0_5 + EMF_5_15)/2)

micAOVunique <- microAnova[-which(duplicated(microAnova$EMF)),]

subMicro <- microAnova[-which(duplicated(microAnova$EMF)),] %>%
  select(., -c(EMF_0_5, EMF_5_15)) %>%
  melt(., id = 'treatment')

ggplot(subMicro, aes(x=treatment, y=value, fill = treatment)) + theme_bw(base_size = 20) +
  geom_boxplot(outliers = F) +
  scale_fill_brewer(palette = 'Set1') +
  labs(y = 'Relative abundance') +
  facet_wrap(~variable, scales = 'free_y') +
  theme(axis.title.x = element_blank())
  

emfAOV <- aov(EMF ~ treatment, data = micAOVunique)
summary(emfAOV)

nitAOV <- aov(nitrifiers ~ treatment, data = micAOVunique)
summary(nitAOV)

TukeyHSD(nitAOV)
