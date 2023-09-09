library(tidyverse)



lut_trt <- c("1A" = "m","1B" = "m","1C" = "bm","1D" = "b",
             "2A" = "bm","2B" = "b","2C" = "m","2D" = "c",
             "3A" = "m","3B" = "c","3C" = "bm","3D" = "b",
             "4A" = "b","4B" = "bm","4C" = "c","4D" = "m",
             "5A" = "b","5B" = "m","5C" = "c","5D" = "bm",
             "6A" = "c","6B" = "b","6C" = "m","6D" = "bm")

d <- read_csv("data/fraser_veg_trait_data.csv") %>%
  mutate(LDMC = dry_weight_g*1000/wet_weight_g,
         trt = lut_trt[plot],
         block = str_sub(plot, 1,1),
         SLA = (leaf_area * (real_val_area/val_area))/dry_weight_g)


# looking at LDMC======================
d %>%
  filter(trait != "SLA") %>%
  ggplot(aes(x=trt, y=LDMC, fill = species_code)) +
  geom_boxplot() +
  facet_wrap(~species_code, scales = "free")

library(lmerTest)
library(emmeans)
library(car)
mc<- d %>%
  filter(trait != "SLA",
         species_code == "CASP") %>%
  lmer((LDMC) ~ trt + (1|block/plot), data = .)
Anova(mc)
performance::check_model(mc)
emmeans::emmeans(mc, ~ trt) %>% pairs()

mv <- d %>%
  filter(trait != "SLA",
         species_code == "VASP") %>%
  lmer(LDMC ~ trt + (1|block/plot), data = .) 
Anova(mv)
summary(mv)
emmeans::emmeans(mv, ~ trt)
emmeans::emmeans(mv, ~ trt) %>% pairs()

mo <- d %>%
  filter(trait != "SLA",
         species_code == "SOSP") %>%
  lmer(LDMC ~ trt + (1|block/plot), data = .)
Anova(mo)
summary(mo)
emmeans::emmeans(mo, ~ trt) %>% pairs()

ma <- d %>%
  filter(trait != "SLA") %>%
  lmer(LDMC ~ trt+species_code + (1|plot), data = .) %>%
  summary()

# looking at SLA ===============

d %>%
  filter(trait != "LDMC",
         species_code != "CASP") %>%
  ggplot(aes(x=trt, y=SLA, fill = species_code)) +
  geom_boxplot() +
  facet_wrap(~species_code, scales = "free_y")

ggsave("out/sla_fig.png")

library(lmerTest)
library(emmeans)
library(car)

mv <- d %>%
  filter(species_code == "VASP") %>%
  lmer(SLA ~ trt + (1|block), data = .) 
Anova(mv)
summary(mv)
emmeans::emmeans(mv, ~ trt)
emmeans::emmeans(mv, ~ trt) %>% pairs()
emmeans::emmeans(mv, ~ trt) %>% contrast()

mo <- d %>%
  filter(species_code == "SOSP") %>%
  lmer(log(SLA) ~ trt + (1|block), data = .)
performance::check_model(mo)
Anova(mo, test.statistic = "F",type = "II")
anova(mo)
summary(mo)
emmeans::emmeans(mo, ~ trt)
emmeans::emmeans(mo, ~ trt) %>% pairs()
emmeans::emmeans(mo, ~ trt) %>% contrast()

agricolae::kruskal(d %>% filter(species_code == "SOSP") %>% pull(SLA), 
                   d %>% filter(species_code == "SOSP") %>% pull(trt)) -> xx

xx

agricolae::kruskal(d %>% filter(species_code == "VASP") %>% pull(SLA), 
                   d %>% filter(species_code == "VASP") %>% pull(trt)) -> xx

xx

# all spp =============
ma <- d %>%
  filter(species_code != "CASP") %>%
  lmer(sqrt(SLA) ~ trt + species_code + (1|block), data = .)
performance::check_model(ma)
Anova(ma, test.statistic = "F",type = "II")
anova(ma)
summary(ma)

ma <- d %>%
  filter(species_code != "CASP") %>%
  lmer(sqrt(SLA) ~ trt * species_code + (1|block), data = .)
performance::check_model(ma)
Anova(ma, test.statistic = "F",type = "II")
anova(ma)
summary(ma)

