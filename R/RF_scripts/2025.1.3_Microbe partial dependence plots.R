library(pacman)
p_load(tidyverse, here, reshape2, randomForest, paletteer, iml, patchwork, ggh4x, tidytext, cowplot, pdp, grid, gridExtra)

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
            #SLA = sla,
            #LDMC = ldmc,
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




# RUN MODELS --------------------------------------------------------------
emfRF <- randomForest(EMF ~., data = uniqueRFdata %>% select(., -c(BA.div, Fungal.div, Nitrifiers)) %>% na.omit,
                      ntree = 5000, mtry = 5, nodesize = 100, maxnodes = 10)

nitriRF <- randomForest(Nitrifiers ~., data = uniqueRFdata %>% select(., -c(BA.div, Fungal.div, EMF)) %>% na.omit,
                        ntree = 5000, mtry = 5, nodesize = 100, maxnodes = 10)

fungDivRF <- randomForest(Fungal.div ~., data = uniqueRFdata %>% select(., -c(BA.div, EMF, Nitrifiers)) %>% na.omit,
                        ntree = 5000, mtry = 5, nodesize = 100, maxnodes = 10)

bactDivRF <- randomForest(BA.div ~., data = uniqueRFdata %>% select(., -c(Fungal.div, EMF, Nitrifiers)) %>% na.omit,
                        ntree = 5000, mtry = 5, nodesize = 100, maxnodes = 10)

# rank variable importance ------------------------------------------------

emfVarImp <- emfRF$importance %>%
  data.frame() %>%
  rownames_to_column() %>%
  subset(., rowname != 'species_full') %>%
  mutate(RelImportance = IncNodePurity/max(IncNodePurity),
         DepVar = 'EMF',
         RelImportance = round(RelImportance, 2)) %>%
  arrange(RelImportance) %>%
  tail(n=5)

nitVarImp <- nitriRF$importance %>%
  data.frame() %>%
  rownames_to_column() %>%
  subset(., rowname != 'species_full') %>%
  mutate(RelImportance = IncNodePurity/max(IncNodePurity),
         DepVar = 'Nitrifiers',
         RelImportance = round(RelImportance, 2)) %>%
  arrange(RelImportance) %>%
  tail(n=5)

bacDivVarImp <- bactDivRF$importance %>%
  data.frame() %>%
  rownames_to_column() %>%
  subset(., rowname != 'species_full') %>%
  mutate(RelImportance = IncNodePurity/max(IncNodePurity),
         DepVar = 'BA.div',
         RelImportance = round(RelImportance, 2)) %>%
  arrange(RelImportance) %>%
  tail(n=5)

fungDivVarImp <- fungDivRF$importance %>%
  data.frame() %>%
  rownames_to_column() %>%
  subset(., rowname != 'species_full') %>%
  mutate(RelImportance = IncNodePurity/max(IncNodePurity),
         DepVar = 'Fungal.div',
         RelImportance = round(RelImportance, 2)) %>%
  arrange(RelImportance) %>%
  tail(n=5)

allVarImp <- rbind(
  emfVarImp %>% select(rowname, DepVar, RelImportance),
  nitVarImp %>% select(rowname, DepVar, RelImportance),
  bacDivVarImp %>% select(rowname, DepVar, RelImportance),
  fungDivVarImp %>% select(rowname, DepVar, RelImportance)
) %>%
  setNames(., nm = c('Predictor', 'microbeType', 'relImportance'))



# extract pdp values ------------------------------------------------------

pdpFunc <- function(model, var, microbe){
  #model = emfRF
  #var = 'NO3'
  #microbe = 'EMF'
  partial(model, var, grid.resolution = 50) %>%
    setNames(., nm = c('x', 'y')) %>%
    mutate(normY = ((y-mean(y))/mean(y))*100,
           Predictor = var,
           microbeType = microbe)
  
}

allPDP <- rbind(
map_dfr(emfVarImp$rowname, ~pdpFunc(emfRF, var = .x, microbe = 'EMF')),
map_dfr(nitVarImp$rowname, ~pdpFunc(nitriRF, var = .x, microbe = 'Nitrifiers')),
map_dfr(bacDivVarImp$rowname, ~pdpFunc(fungDivRF, var = .x, microbe = 'BA.div')),
map_dfr(fungDivVarImp$rowname, ~pdpFunc(bactDivRF, var = .x, microbe = 'Fungal.div'))
)


joinData <- left_join(allPDP, allVarImp)

baseSize <-  27
axisFontSize <- 15


EMFplot <- joinData %>% subset(microbeType == 'EMF') %>%
  mutate(Predictor = factor(Predictor, levels = emfVarImp$rowname) %>% fct_rev) %>%
  ggplot(., aes(x=x, y=normY)) + theme_bw(base_size = baseSize) +
  geom_hline(yintercept = 0, color = 'black') +
  geom_line(size = 3, color = '#A16928FF') + 
  #scale_alpha_continuous(range = c(0.25, 1)) +
  #paletteer::scale_color_paletteer_c("ggthemes::Orange") +
  facet_grid2(microbeType~Predictor, scales = 'free_x', independent = 'x') +
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_text(size = axisFontSize))

nitriPlot <- joinData %>% subset(microbeType == 'Nitrifiers') %>%
  mutate(Predictor = factor(Predictor, levels = nitVarImp$rowname) %>% fct_rev) %>%
  ggplot(., aes(x=x, y=normY)) + theme_bw(base_size = baseSize) +
  geom_hline(yintercept = 0, color = 'black') +
  geom_line(size = 3, color = "#D6BD8DFF") + 
  #scale_alpha_continuous(range = c(0.25, 1)) +
  #paletteer::scale_color_paletteer_c("ggthemes::Blue") +
  facet_grid2(microbeType~Predictor, scales = 'free_x', independent = 'x') +
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_text(size = axisFontSize))

BAdivPlot <- joinData %>% subset(microbeType == 'BA.div') %>%
  mutate(Predictor = factor(Predictor, levels = bacDivVarImp$rowname) %>% fct_rev) %>%
  ggplot(., aes(x=x, y=normY)) + theme_bw(base_size = baseSize) +
  geom_hline(yintercept = 0, color = 'black') +
  geom_line(size = 3, color = "#B5C8B8FF") + 
  #scale_alpha_continuous(range = c(0.25, 1)) +
  #paletteer::scale_color_paletteer_c("ggthemes::Green") +
  facet_grid2(microbeType~Predictor, scales = 'free_x', independent = 'x') +
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_text(size = axisFontSize))

FungdivPlot <- joinData %>% subset(microbeType == 'Fungal.div') %>%
  mutate(Predictor = factor(Predictor, levels = fungDivVarImp$rowname) %>% fct_rev) %>%
  ggplot(., aes(x=x, y=normY)) + theme_bw(base_size = baseSize) +
  geom_hline(yintercept = 0, color = 'black') +
  geom_line(size = 3, color = "#2887A1FF") + 
  #scale_alpha_continuous(range = c(0.25, 1)) +
  #paletteer::scale_color_paletteer_c("ggthemes::Brown") +
  facet_grid2(microbeType~Predictor, scales = 'free_x', independent = 'x') +
  theme(axis.title = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = axisFontSize))


plot <- plot_grid(EMFplot, nitriPlot, BAdivPlot, FungdivPlot, ncol = 1, align = 'vh')
title <- ggdraw() + draw_label(expression(Decreasing.Variable.Importance * symbol('\256')), size = 25)


#plot <- plot_grid(title, plot, ncol = 1, rel_heights = c(0.05, 1), align = 'h')

y.grob <- textGrob('Standardized change (%)', gp=gpar(fontsize = 30), rot = 90)
x.grob <- textGrob(expression(Decreasing.Variable.Importance * symbol('\256')), gp=gpar(fontsize=25))
finalPlot <- grid.arrange(arrangeGrob(plot, left = y.grob), top = x.grob)

ggsave(here('Figures', 'Fig 5_microbe pdp.png'), finalPlot, height = 30, width = 40, units = 'cm', dpi = 300)

