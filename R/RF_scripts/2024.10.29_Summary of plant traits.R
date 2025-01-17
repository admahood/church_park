library(pacman)
p_load(tidyverse, here, reshape2, randomForest, tidytext, ggh4x, cowplot)

iData <- read.csv(here('Data', 'cp_trait_soil_data.csv')) %>%
  subset(., species_full != 'Carex sp')

# statistics for results
slaAOV <- aov(sla ~ species_full * treatment, data = iData)
summary(slaAOV)
TukeyHSD(slaAOV, which = 'treatment')

ldmcAOV <- aov(ldmc ~ species_full + treatment, data = iData)
summary(ldmcAOV)
TukeyHSD(htAOV, which = 'treatment')

htAOV <- aov(height_cm ~ species_full * treatment, data = iData)
summary(htAOV)
TukeyHSD(htAOV, which = 'treatment')

lm(sla ~ ldmc, data = iData %>% subset(species_full == 'Oreochrysum parryi')) %>% summary




# Generate figure ---------------------------------------------------------

reData <- iData %>%
  mutate(
    species_full = case_when(
    str_detect(species_full, regex('Vaccinium sp')) ~ 'V. scoparium',
    str_detect(species_full, regex('Oreochrysum parryi')) ~ 'O. parryi'),
    treatment = case_when(
    str_detect(treatment, regex('\\bm\\b')) ~ 'M',
    str_detect(treatment, regex('\\bb\\b')) ~ 'B',
    str_detect(treatment, regex('\\bbm\\b')) ~ 'BM',
    str_detect(treatment, regex('0')) ~ 'CTL'
   )) %>%
  mutate(treatment = factor(treatment, levels = c('CTL', 'B', 'BM', 'M'))) %>%
  rename(SLA = sla,
         LDMC = ldmc,
         Height = height_cm)



#### only boxplot and two significant relationships

boxplots_SLA <- reData %>%
  select(., species_full, treatment, SLA) %>%
  melt(id = c('species_full', 'treatment')) %>%
  ggplot(., aes(x=treatment, y=value, fill = treatment)) + theme_bw(base_size = 20) +
  geom_boxplot(outliers = F) +
  scale_fill_brewer(palette = 'Set1') +
  labs(y = 'Plant functional trait') +
  facet_grid2(variable ~ species_full) +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank())

boxplots_LDMC <- reData %>%
  select(., species_full, treatment, LDMC) %>%
  melt(id = c('species_full', 'treatment')) %>%
  ggplot(., aes(x=treatment, y=value, fill = treatment)) + theme_bw(base_size = 20) +
  geom_boxplot(outliers = F) +
  scale_fill_brewer(palette = 'Set1') +
  labs(y = 'Plant functional trait') +
  facet_grid2(variable ~ species_full) +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

boxplots_ht <- reData %>%
  select(., species_full, treatment, Height) %>%
  melt(id = c('species_full', 'treatment')) %>%
  ggplot(., aes(x=treatment, y=value, fill = treatment)) + theme_bw(base_size = 20) +
  geom_boxplot(outliers = F) +
  scale_fill_brewer(palette = 'Set1') +
  labs(y = 'Plant functional trait') +
  facet_grid2(variable ~ species_full) +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_blank())

boxplots_all <- plot_grid(boxplots_SLA, boxplots_LDMC, boxplots_ht, ncol = 1, align = 'v', rel_heights = c(0.4, 0.35, 0.405))

sla_v_ldmc_significant <- ggplot(reData, aes(x=SLA, y=LDMC)) + 
  theme_bw(base_size = 20) +
  geom_point(aes(fill = treatment), shape = 21, color = 'gray10', size = 3) + 
  geom_smooth(method = 'lm', se = F, color = 'gray90', size = 2) +
  geom_smooth(method = 'lm', se = F, color = 'black') +
  scale_fill_brewer(palette = 'Set1') +
  geom_rug(color = 'gray35') + 
  facet_wrap(~species_full)+
  labs(x='SLA', y='LDMC') + 
  theme(legend.position = 'none')


plot_grid(boxplots_all, sla_v_ldmc_significant, ncol = 1, align = 'h', axis = 'l', rel_heights = c(0.7, 0.3))

ggsave(here('Figures', 'Fig3_traits boxplots.png'), height = 35, width = 25, dpi = 300, units = 'cm')










