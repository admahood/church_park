# Church Park vegetation composition analysis

library(tidyverse)
library(vegan)
library(brms)

sites <- readxl::read_xlsx("data/church_park_cover.xlsx",sheet = "site_info") %>%
  mutate(plot = str_c("b", block, plot), treatment = treatment %>% str_replace_all("c", "0"),
         block = as.factor(block)) 

comm_long <- readxl::read_xlsx("data/church_park_cover.xlsx") %>%
  pivot_longer(cols = names(.)[2:length(.)],values_drop_na = T) %>%
  tidyr::separate(name, c("block", "treatment", "quadrat"), sep = "_") %>%
  mutate(species = str_remove_all(species," \\(native\\)") %>% str_to_title());comm_long
unique(comm_long$species)

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

# sites <-
#   left_join(comm_wide %>% dplyr::select(block, treatment, quadrat) %>%
#               transmute(plot = str_c(block, treatment), quadrat=quadrat),
#             sites0 %>% mutate(block = as.character(block))) %>%
#   mutate(rowname = str_c(plot,quadrat))

# sites$rowname == rownames(comm)
# brms jsdm ====================================================================
bformula <- formula(str_c("mvbind(",
  paste(names(comm), collapse = ","),
  ") ~ treatment + (1|block)"
))

# bformula <- str_c("mvbind(",
#                           paste(names(comm)[1:10], collapse = ","),
#                           ") ~ treatment + (1|block)"
# )
# prior1 <- prior(normal(0,10), class = b)

# 44 minutes for 2000
to <- Sys.time()
zbmd <- brm(bformula, # + set_rescor(rescor = TRUE), 
           data = cbind(sites, comm/100),
           iter = 4000,
           family = "gaussian"
           )
t1 <- Sys.time()
print(t1-to)

summary(zbmd)

t2 <- Sys.time()
zbmd1 <- brm(bformula, # + set_rescor(rescor = TRUE), 
            data = cbind(sites, comm/100),
            iter = 6000,
            stan_model_args = list(max_treedepth = 20),
            family = "zero_inflated_beta")
t3 <- Sys.time()
print(t3-t2)

# with random effects
broom.mixed::tidy(zbmd) %>%
  filter(str_sub(term, 1,2) != "sd", term != "(Intercept)") %>% 
  mutate(term = str_remove_all(term, "treatment")) %>%
  ggplot(aes(x=estimate, y=response)) +
  geom_vline(xintercept = 0, lty=2) +
  geom_point() +
  facet_wrap(~term, scales="free_x", nrow=1) +
  geom_segment(aes(x = conf.low, xend=conf.high, y=response, yend=response), lwd=1)

