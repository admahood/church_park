# sem for traits
source("R/a_church_data_prep.R")
source("R/a_itv_data_prep.R")
library(lavaan)
library(ggsem) # devtools::install_github("admahood/ggsem")
library(ggtext)

# data_prep

d <- d_church |>
  dplyr::select(-plot) |>
  left_join(sites_w_23_soil, by = c("treatment", "block")) |>
  left_join(gc_by_plot_23) |>
  dplyr::rename(vwc = mean_vwc, mulch = wood_mulch) |> 
  mutate_if(is.numeric, datawizard::standardise) |>
  mutate_if(is.character, as.factor) |>
  left_join(pcab) |>
  left_join(pcaf)
summary(d)

glimpse(d)

## plot of sla ~ ldmc after accounting for folded aspect
# 
# ggplot(d, aes(ldmc, sla, color = folded_aspect)) +
#   geom_text(aes(label = block)) +
#   scale_color_viridis_c()+
#   facet_wrap(~species_full, scales = "free")
# 
# lm(ldmc ~ sla + folded_aspect, data = d |> filter(species_full == "Carex sp")) |>
#   summary()
# lm(ldmc ~ sla + folded_aspect, data = d |> filter(species_full == "Oreochrysum parryi")) |>
#   summary()
# lm(ldmc ~ sla + folded_aspect, data = d |> filter(species_full == "Vaccinium sp")) |>
#   summary()


###### all spp together, grouped ===============================================

moda <- 'sla ~  doc_mg_g  + folded_aspect+mulch + din_mg_g + soil_p_h + height_cm +vwc + percent_sand_0to5
         ldmc ~  doc_mg_g  + folded_aspect+mulch + din_mg_g + soil_p_h + height_cm +vwc + percent_sand_0to5
        vwc ~ mulch  + folded_aspect + biochar
        tdn_mg_g ~ mulch + biochar + folded_aspect + din_mg_g + soil_p_h + ca_mg_g  + vwc
        ca_mg_g ~ mulch + biochar + vwc + doc_mg_g + folded_aspect
        soil_p_h ~ biochar + vwc + folded_aspect + ca_mg_g#+ tdn_mg_g mulch + 
        din_mg_g ~ mulch + biochar + vwc + folded_aspect + soil_p_h + ca_mg_g + pcab 
        doc_mg_g ~ mulch + biochar + vwc + folded_aspect + soil_p_h + tdn_mg_g+ pcaf# + shannon_bacteria_5
        height_cm ~ vwc +  doc_mg_g + tdn_mg_g + ca_mg_g +  soil_p_h + mulch + pcaf
        #shannon_bacteria_5 ~ vwc + mulch + biochar + soil_p_h
        tdn_mg_g ~~ doc_mg_g
        ' 
fit <- lavaan::sem(moda, data = d, group = 'species_full')
fit
summary(fit)
# ggsem(fit, title = paste("all, p =", round(fit@test$standard$pvalue,2)))

modificationindices(fit) |> arrange(desc(mi)) |> head(n=20)
resid(fit, "cor")#$cov
lavaan::fitMeasures(fit) %>%
  round(3) %>%
  as_tibble(rownames = "metric") %>%
  filter(metric %in% c("tli", "cfi", "rmsea", "srmr"))

# oreochrysum ==========

modo <-'ldmc ~  doc_mg_g  + folded_aspect + din_mg_g + soil_p_h + pcaf + vwc
        sla ~  doc_mg_g  + folded_aspect+ din_mg_g + soil_p_h+ tdn_mg_g
        vwc ~ mulch  + folded_aspect + biochar + percent_sand_0to5 + twi + ca_mg_g
         pcab ~ ca_mg_g + mulch + tdn_mg_g+ percent_sand_0to5 + folded_aspect + twi
         pcaf ~  mulch + din_mg_g + ca_mg_g + tdn_mg_g + soil_p_h
         tdn_mg_g ~ twi + mulch + biochar + folded_aspect + soil_p_h + ca_mg_g  + vwc + pcab + pcaf + doc_mg_g+ percent_sand_0to5 + din_mg_g
         ca_mg_g ~ biochar + vwc + twi  + pcab
         soil_p_h ~ twi + biochar + vwc + folded_aspect + ca_mg_g+doc_mg_g#+ tdn_mg_g mulch + 
         din_mg_g ~ biochar + vwc + folded_aspect + soil_p_h + ca_mg_g + tdn_mg_g + pcab + pcaf + percent_sand_0to5 + twi
         doc_mg_g ~ mulch + biochar  + folded_aspect  + tdn_mg_g+ percent_sand_0to5 + twi + ca_mg_g
         sla ~~ vwc
         ' 
fito <- lavaan::sem(modo, data = d |> filter(species_full == "Oreochrysum parryi"))
fito

modificationindices(fito) |> arrange(desc(mi)) |> head()
summary(fito, rsquare=T, fit.measures =T)
resid(fito, "cor")$cov
lavaan::fitMeasures(fito) %>%
  round(3) %>%
  as_tibble(rownames = "metric") %>%
  filter(metric %in% c("tli", "cfi", "rmsea", "srmr"))

layout_df <- ggsem::random_layout(fito) |>
  mutate(x = case_match(metric,
                        'ldmc' ~ 1,
                        "sla" ~ 1,
                        'pcaf' ~.45,
                        'tdn_mg_g' ~ .25,
                        'vwc' ~ .5,
                        'pcab' ~ .45,
                        'ca_mg_g' ~ .75, 
                        'soil_p_h' ~ .75,
                        "din_mg_g" ~ .75,
                        'doc_mg_g' ~ .5,
                        'percent_sand_0to5' ~ 0,
                        'folded_aspect' ~ 0,
                        'mulch' ~ 0,
                        'twi' ~ 0,
                        'biochar'~ 0),
         y = case_match(metric,
                        'ldmc' ~ .75,
                        "sla" ~ 0.25,
                        'vwc' ~ 1,
                        'pcaf' ~.7,
                        'tdn_mg_g' ~ .5,
                        'pcab' ~ .25,
                        'ca_mg_g' ~ .8, 
                        'soil_p_h' ~ .25,
                        "din_mg_g" ~ .5,
                        'doc_mg_g' ~ 0,
                        'folded_aspect' ~ 0,
                        'mulch' ~ .5,
                        'twi' ~ .25,
                        'percent_sand_0to5' ~ .75,
                        'biochar'~ 1))

lut_names <- c("LDMC", "SLA", "VWC", "Bacteria", "Fungi", "TDN", "Ca",
                      "Soil pH", "DIN", "DOC", "Aspect", "%Sand", "Mulch", "Biochar", "TWI" )
names(lut_names) <- layout_df$metric
layout_df <- layout_df |> mutate(new_node_names = lut_names[metric])
po <- ggsem(fito, layout_df = layout_df, layout = 'manual', rename_nodes = T,
            new_node_names = lut_names, show_legend = TRUE,
            labels = F, cols = c("transparent", "#E41A1C", "#377EB8"),
            title = paste("<i>Oreochrysum parryi</i>, p =", 
                          round(fito@test$standard$pvalue,2))) +
  theme(plot.title = element_markdown(),
        legend.position = c(1,1),
        legend.justification = c(1,1));po

po1 <- ggsem(fito, layout_df = layout_df, layout = 'manual', rename_nodes = T,
            new_node_names = lut_names, show_legend = TRUE, 
            exclude = c("pcab", "vwc", "biochar", "tdn_mg_g", "ca_mg_g"),
            labels = F, cols = c("transparent", "#E41A1C", "#377EB8"),
            title = paste("<i>Oreochrysum parryi</i>, p =", 
                          round(fito@test$standard$pvalue,2))) +
  theme(plot.title = element_markdown(),
        legend.position = c(1,1),
        legend.justification = c(1,1));po1

ggsave("out/oreo_sem.png", bg = "white")


# vaccinium ===============

modv<-'ldmc ~  doc_mg_g  + folded_aspect + din_mg_g + soil_p_h + pcaf + vwc + twi + ca_mg_g 
        sla ~  vwc + doc_mg_g  + folded_aspect+ din_mg_g + soil_p_h+ tdn_mg_g + percent_sand_0to5 + pcaf + ca_mg_g 
        vwc ~ mulch  + folded_aspect + biochar + percent_sand_0to5 + twi + ca_mg_g
         pcab ~ ca_mg_g + mulch + tdn_mg_g+ percent_sand_0to5 + folded_aspect 
         pcaf ~  mulch + ca_mg_g + tdn_mg_g + soil_p_h + biochar + folded_aspect
         tdn_mg_g ~ twi + mulch + biochar + folded_aspect + soil_p_h + ca_mg_g  + vwc + pcab + pcaf + doc_mg_g+ percent_sand_0to5 + din_mg_g
         ca_mg_g ~ biochar + vwc  + pcab + folded_aspect + percent_sand_0to5
         soil_p_h ~ twi + biochar + vwc + folded_aspect + ca_mg_g+doc_mg_g#+ tdn_mg_g mulch + 
         din_mg_g ~ biochar + vwc + folded_aspect + soil_p_h + ca_mg_g + tdn_mg_g + pcab + pcaf + percent_sand_0to5 + twi
         doc_mg_g ~ mulch + biochar  + folded_aspect  + tdn_mg_g+ percent_sand_0to5 + twi + ca_mg_g
         ' 

fitv <- lavaan::sem(modv, data = d |> filter(species_full == "Vaccinium sp"))
fitv
modificationindices(fitv) |> arrange(desc(mi)) 
summary(fitv, rsquare=T, fit.measures =T)
resid(fitv, "cor")$cov

pv<- ggsem(fitv, layout_df = layout_df, layout = 'manual',rename_nodes = T,
           labels = FALSE,
           cols = c("transparent", "#E41A1C", "#377EB8"), new_node_names = lut_names,
           title = paste("<i>Vaccinium</i> sp., p =", round(fitv@test$standard$pvalue,2))) +
  theme(plot.title = element_markdown());pv
ggsave("out/vaccinium_sem.png", bg = "white")

pv1<- ggsem(fitv, layout_df = layout_df, layout = 'manual',rename_nodes = T,
           labels = FALSE,
           exclude = c("din_mg_g", "mulch", "tdn_mg_g", "pcab"),
           cols = c("transparent", "#E41A1C", "#377EB8"), new_node_names = lut_names,
           title = paste("<i>Vaccinium</i> sp., p =", round(fitv@test$standard$pvalue,2))) +
  theme(plot.title = element_markdown());pv1



# carex =================================

modc <- 'ldmc ~  doc_mg_g  + folded_aspect + din_mg_g + soil_p_h + pcaf + vwc + twi + ca_mg_g + na_mg_g + simpson_fungi_5
        sla ~  simpson_fungi_5 + vwc + doc_mg_g  + folded_aspect+ din_mg_g + soil_p_h+ tdn_mg_g + percent_sand_0to5 + pcaf + ca_mg_g 
        vwc ~ mulch  + folded_aspect + biochar + percent_sand_0to5 + twi + ca_mg_g
         na_mg_g ~ twi + ca_mg_g + tdn_mg_g + folded_aspect
         pcab ~ ca_mg_g  + tdn_mg_g+ percent_sand_0to5 + folded_aspect 
         pcaf ~  mulch + ca_mg_g + tdn_mg_g + soil_p_h + biochar + folded_aspect
         tdn_mg_g ~ simpson_fungi_5 + na_mg_g + twi + mulch + biochar + folded_aspect + soil_p_h + ca_mg_g  + vwc + pcab + pcaf + doc_mg_g+ percent_sand_0to5 + din_mg_g
         ca_mg_g ~ doc_mg_g + biochar + vwc  + pcab + folded_aspect + percent_sand_0to5
         soil_p_h ~ percent_sand_0to5 + na_mg_g + twi + biochar + vwc + folded_aspect + ca_mg_g+doc_mg_g#+ tdn_mg_g mulch + 
         din_mg_g ~ biochar + vwc + folded_aspect + soil_p_h + ca_mg_g + tdn_mg_g + pcab + pcaf + percent_sand_0to5 + twi
         doc_mg_g ~ din_mg_g + simpson_fungi_5 + mulch + biochar  + folded_aspect  + tdn_mg_g+ percent_sand_0to5 + twi + ca_mg_g
         ' 

fitc <- lavaan::sem(modc, data = d |> filter(species_full == "Carex sp"))
fitc

modificationindices(fitc) |> arrange(desc(mi)) |> head()
summary(fitc, rsquare=T, fit.measures =T)
resid(fitc, "cor")$cov
lavaan::fitMeasures(fitc) %>%
  round(3) %>%
  as_tibble(rownames = "metric") %>%
  filter(metric %in% c("tli", "cfi", "rmsea", "srmr"))

pc <- ggsem(fitc, layout_df = layout_df, layout = 'auto',#rename_nodes = T,
            labels = FALSE,
            cols = c("transparent", "#E41A1C", "#377EB8"), new_node_names = lut_names,
            title = paste("<i>Carex</i> sp., p =", round(fitc@test$standard$pvalue,2))) +
  theme(plot.title = element_markdown());pc
ggsave("out/carex_sem.png", bg = "white")

pc1 <- ggsem(fitc, layout_df = layout_df, layout = 'manual',rename_nodes = T,
            labels = FALSE, exclude =c(  "folded_aspect" , "vwc","pcab","pcaf","tdn_mg_g","ca_mg_g","soil_p_h",     
                                         "din_mg_g" ,           "biochar"   ),
            cols = c("transparent", "#E41A1C", "#377EB8"), new_node_names = lut_names,
            title = paste("<i>Carex</i> sp., p =", round(fitc@test$standard$pvalue,2))) +
  theme(plot.title = element_markdown());pc1



# no traits =================
mod0 <- 'vwc ~ mulch  + folded_aspect + biochar + percent_sand_0to5 + twi
         pcab ~ vwc + ca_mg_g + mulch + tdn_mg_g+ percent_sand_0to5 + folded_aspect + pcaf + twi
         pcaf ~ vwc + ca_mg_g + pcab + tdn_mg_g+ percent_sand_0to5 + twi
         tdn_mg_g ~ twi + mulch + biochar + folded_aspect + soil_p_h + ca_mg_g  + vwc + pcab + pcaf + doc_mg_g+ percent_sand_0to5 + din_mg_g
         ca_mg_g ~ biochar + vwc + doc_mg_g + folded_aspect + twi + mulch
         soil_p_h ~ twi + biochar + vwc + folded_aspect + ca_mg_g+doc_mg_g#+ tdn_mg_g mulch + 
         din_mg_g ~ biochar + vwc + folded_aspect + soil_p_h + ca_mg_g + tdn_mg_g + pcab + pcaf + percent_sand_0to5 + twi
         doc_mg_g ~ mulch + biochar + vwc + folded_aspect + soil_p_h + tdn_mg_g+ percent_sand_0to5 + twi + ca_mg_g
         ' 

fit0 <- lavaan::sem(mod0, data = d)
fit0


p0 <- ggsem(fit0, layout_df = layout_df, layout = 'manual',rename_nodes = T,
            labels = FALSE,
            cols = c("transparent", "#377EB8", "#E41A1C"),
            new_node_names = lut_names,
            title = paste("Soil Only, p =", round(fit0@test$standard$pvalue,2))) +
  theme(plot.title = element_markdown());p0

ggsave(plot = p0, filename = "out/soil_only_sem.png", width =9, height = 7, bg="white")

modificationindices(fit0) |> arrange((mi)) 
summary(fit0)

# multipanel===========
ggpubr::ggarrange(po, pv, pc, p0, nrow=2, ncol=2)
ggsave(filename = "out/multipanel_sem.png", bg="white", height = 15, width = 17)


ggpubr::ggarrange(po1, pv1, pc1, p0, nrow=2, ncol=2, labels = "auto")
ggsave(filename = "out/multipanel_sem1.png", bg="white", height = 15, width = 14)


# random forest

library(randomForest)
rfd <- d |> dplyr::select(-plot, -sample_id, -lma, -individual, -ldmc, -latitude, -height_cm,
                          -longitude, -PC2)

randomForest(sla ~ ., rfd |> na.omit() |> filter(species_full == "Carex sp"), ntree = 2000) -> sla_rf
sla_rf; plot(sla_rf)
varImpPlot(sla_rf)
randomForest(sla ~ ., rfd |> na.omit() |> filter(species_full == "Oreochrysum parryi"), ntree = 2000) -> sla_rf
sla_rf; plot(sla_rf)
varImpPlot(sla_rf)
randomForest(sla ~ ., rfd |> na.omit() |> filter(species_full == "Vaccinium sp"), ntree = 2000) -> sla_rf
sla_rf; plot(sla_rf)
varImpPlot(sla_rf)

rfd <- d |> dplyr::select(-plot, -sample_id, -lma, -individual, -sla, -latitude,
                          -longitude, -PC2, -height_cm)

randomForest(ldmc ~ ., rfd |> na.omit() |> filter(species_full == "Carex sp"), ntree = 2000) -> sla_rf
sla_rf
varImpPlot(sla_rf)
randomForest(ldmc ~ ., rfd |> na.omit() |> filter(species_full == "Oreochrysum parryi"), ntree = 2000) -> sla_rf
sla_rf
varImpPlot(sla_rf)
randomForest(ldmc ~ ., rfd |> na.omit() |> filter(species_full == "Vaccinium sp"), ntree = 2000) -> sla_rf
sla_rf
varImpPlot(sla_rf)
plot(sla_rf)


# (g)lmms =============
library(lmerTest)
library(performance)
library(broom.mixed)

# ldmc
lme_o <- lmer(ldmc ~  doc_mg_g  + folded_aspect + din_mg_g + soil_p_h + pcaf + 
                vwc + (1|block/plot),
             data = d |> filter(species_full == "Oreochrysum parryi"))
summary(lme_o); check_model(lme_o); car::Anova(lme_o)

lme_v <- lmer(ldmc ~  folded_aspect + din_mg_g + 
                soil_p_h + pcaf + twi +(1|block/plot),
              data = d |> filter(species_full == "Vaccinium sp"))
summary(lme_v); check_model(lme_v); car::Anova(lme_v)

# Carex ldmc
lme_c <- lmer(ldmc ~   doc_mg_g  + folded_aspect + din_mg_g + soil_p_h + pcaf + vwc + twi + ca_mg_g + na_mg_g + simpson_fungi_5
              + (1|block/plot),
              data = d |> filter(species_full == "Carex sp"))
summary(lme_c); performance::check_model(lme_c); car::Anova(lme_c)

# plotting lmers

tidy(lme_o) |>
  mutate(species = "Oreochrysum") |>
  bind_rows(tidy(lme_c1) |> mutate(species = "Carex")) |>
  bind_rows(tidy(lme_v) |> mutate(species = "Vaccinium")) |>
  dplyr::filter(effect == "fixed", term !='(Intercept)') |>
  ggplot(aes(x=estimate, y = term)) +
  geom_point() +
  geom_segment(aes(x = estimate-std.error, xend = estimate+std.error)) +
  geom_vline(xintercept = 0, lty=2) +
  facet_wrap(~species) +
  ggtitle("LDMC")
ggsave("out/lmers_ldmc.png")

# sla

lme_o <- lmer(sla ~ doc_mg_g  + folded_aspect + din_mg_g + soil_p_h+ tdn_mg_g + (1|block/plot),
              data = d |> filter(species_full == "Oreochrysum parryi"))
summary(lme_o); check_model(lme_o); car::Anova(lme_o)

lme_v <- lmer(sla ~ folded_aspect + 
                din_mg_g + soil_p_h+ tdn_mg_g + percent_sand_0to5 + pcaf 
              + (1|block/plot),
              data = d |> filter(species_full == "Vaccinium sp"))
summary(lme_v); check_model(lme_v); car::Anova(lme_v)

# Carex 
lme_c <- lmer(sla ~  simpson_fungi_5 + vwc + soil_p_h  + 
                 + simpson_bacteria_5 + biochar  + (1|block/plot),
              data = d |> filter(species_full == "Carex sp"))

summary(lme_c); performance::check_model(lme_c); car::Anova(lme_c)

AIC(lme_c, lme_v)

# plotting lmers

tidy(lme_o) |>
  mutate(species = "Oreochrysum") |>
  bind_rows(tidy(lme_c) |> mutate(species = "Carex")) |>
  bind_rows(tidy(lme_v) |> mutate(species = "Vaccinium")) |>
  dplyr::filter(effect == "fixed", term !='(Intercept)') |>
  ggplot(aes(x=estimate, y = term)) +
  geom_point() +
  geom_segment(aes(x = estimate-std.error, xend = estimate+std.error)) +
  geom_vline(xintercept = 0, lty=2) +
  facet_wrap(~species) +
  ggtitle("SLA")
ggsave("out/lmers_sla.png")


