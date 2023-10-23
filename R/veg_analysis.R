# Church Park vegetation composition analysis

library(tidyverse)
library(vegan)

sites <- readxl::read_xlsx("data/church_park_cover.xlsx",sheet = "site_info") %>%
  transmute(plot = str_c("b", block, plot), treatment = treatment)

comm_long <- readxl::read_xlsx("data/church_park_cover.xlsx") %>%
  pivot_longer(cols = names(.)[2:length(.)],values_drop_na = T) %>%
  tidyr::separate(name, c("block", "treatment", "quadrat"), sep = "_");comm_long


comm_wide <- comm_long %>%
  group_by(block, treatment, species) %>%
  summarise(cover = sum(value)/2) %>% # just doing mean does not work
  ungroup() %>%
  mutate(species = str_remove_all(species,"\\.") %>% str_replace_all(" ", "_")) %>%
  pivot_wider(id_cols = c(block, treatment), 
              names_from = species, 
              values_fill = 0,
              values_from = cover); comm_wide 

comm <- comm_wide %>%
  mutate(row = str_c(block, treatment)) %>%
  dplyr::select(-block, -treatment) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("row"); comm 
 
nmds <- comm %>%
  wisconsin() %>%
  metaMDS()

stressplot(nmds)

site_scores <- as.data.frame(vegan::scores(nmds)$sites) %>%
  as_tibble(rownames = "plot") %>%
  left_join(sites); site_scores

ef <- envfit(nmds, comm, na.rm = T, permutations = 9999)
eft <- envfit(nmds, sites[,2], na.rm = T, permutations = 9999)
sp <-as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))
species <- as.data.frame(cbind(sp, p=ef$vectors$pvals)) %>%
  tibble::rownames_to_column("species")  %>%
  filter(p < 0.05)


p_ab <- ggplot(site_scores, aes(x=NMDS1, y=NMDS2)) +
  coord_fixed() +
  geom_point(size=2, aes(color = treatment)) +
  geom_segment(data = species,x=0,y=0, color = "grey",arrow = arrow(),
               aes(yend = NMDS2, xend = NMDS1), lwd=1)+
  geom_text(data = species,size=4, aes(label = species), color = "grey40") +
  theme_classic() +
  theme(panel.background = element_rect(fill="transparent", color = "black")) +
  ggtitle("Abundance-Based");p_ab

adonis2(comm ~ treatment, data = sites)



