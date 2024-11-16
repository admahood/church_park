# piecewise sem
library(lmerTest)
library(piecewiseSEM)
library(ggsem)
library(ggeffects)
library(ggpubr)
library(ggdist)
# library(semEff)
library(MuMIn)
source("R/a_church_data_prep.R")
source("R/a_itv_data_prep.R")
library(lavaan)
library(ggsem) # devtools::install_github("admahood/ggsem")
library(ggtext)
library(semEff)
# try it with simulated values of mulch and biochar
# data_prep

d <- d_church |>
  dplyr::select(-plot) |>
  left_join(sites_w_23_soil, by = c("treatment", "block")) |>
  left_join(gc_by_plot_23) |>
  dplyr::rename(vwc = mean_vwc, mulch = wood_mulch) |> 
  # mutate_if(is.numeric, datawizard::standardise) |>
  mutate_if(is.character, as.factor) |>
  left_join(nitrifiers) |>
  left_join(discriminant_taxa) |>
  left_join(soil_23_tf_both_depths) |>
  mutate(DOCtoTDN_0_5 = dissolved_organic_carbon_0_5/dissolved_organic_carbon_0_5,
         DOCtoTDN_5_15 = dissolved_organic_carbon_5_15/dissolved_organic_carbon_5_15,
         CN = total_c_0_5/total_n_0_5,
         CN_5_15 = total_n_5_15/total_n_5_15,
         din_tf = ammonium_0_5 + nitrate_0_5) |>
  dplyr::rename(DOC = dissolved_organic_carbon_0_5,
                DOC_5_15 = dissolved_organic_carbon_5_15,
                PPF = plant_pathogen_5,
                Aspect = folded_aspect,
                DIN = din_tf,
                TDN = dissolved_nitrogen_0_5,
                TDN_5_15 = dissolved_nitrogen_5_15,
                EMF = ectomycorrhizae_5,
                EMF_5_15 = ectomycorrhizae_15,
                pH = soil_p_h,
                Cations = cation_sum_0_5,
                cations_5_15 = cation_sum_5_15,
                Sand = percent_sand_0to5) |>
  mutate(DON = TDN - DIN,
         Clay = percent_clay_0to5,
         NH4_sum = ammonium_0_5 + ammonium_5_15,
         ammonium_mean = (ammonium_0_5 + ammonium_5_15)/2,
         NO3_sum = nitrate_0_5 + nitrate_5_15,
         a_n_ratio = NH4_sum/NO3_sum,
         an5 = ammonium_0_5/nitrate_0_5,
         an15 = ammonium_5_15/nitrate_5_15) |>
  mutate_if(is.character, as.factor) |>
  dplyr::select(-ends_with('mg_g'), -c, -n, -lma, -sample_id) |>
  left_join(tvc) |>
  dplyr::rename(NH4 = ammonium_0_5,
                NH4_15 = ammonium_5_15,
                NO3 = nitrate_0_5,
                NO3_15 = nitrate_5_15,
                Mulch = mulch,
                Biochar = biochar,
                VWC = vwc,
                Nitrifiers = nitrifiers,
                Bare = bare,
                TWI = twi,
                Fungal_div = simpson_fungi_15) #|>
  #datawizard::standardise()
# write_csv(d, "data/cp_trait_soil_data.csv") 
oreo <- filter(d, species_full == "Oreochrysum parryi"); summary(oreo)
vac <- filter(d, species_full == "Vaccinium sp"); summary(vac)
car <- filter(d, species_full == "Carex sp") |>
  filter(!is.na(sla)); summary(car)

# lmers for the traits

# oreochrysum parryi ===========================================================
# LDMC
mol <- lmer(ldmc ~ NO3 + Aspect + CN + # Sand +#TWI +# EMF +#pH + EMF + DOC + PPF + TDN + 
               (1|block),
            data = oreo, na.action = na.fail); summary(mol)
mnitrate <- lmer(NO3 ~ pH + Nitrifiers + EMF + Cations + Mulch + (1|block),
              data = oreo, na.action = na.fail)
mph <- lmer(pH ~  Nitrifiers + Cations + Mulch + (1|block), 
            data = oreo, na.action = na.fail)
mcn <- lmer(CN ~ Biochar + Aspect + (1|block), data = oreo)
memf <- lmer(EMF ~ pH + Nitrifiers + Mulch + (1|block), 
             data = oreo, na.action = na.fail)
# mppf <- lmer(PPF ~ Biochar + DIN + pH  + TDN + (1|block), 
#              data = oreo, na.action = na.fail)
mtdn <- lmer(TDN ~ Mulch + Biochar + EMF + VWC + pH + DIN + Nitrifiers +(1|block),
             data = oreo, na.action = na.fail)
mnit <- lmer(Nitrifiers ~ Mulch + Biochar + VWC + Clay + (1|block),
             data = oreo, na.action = na.fail)
mcat <- lmer(Cations ~ Biochar + VWC + Mulch + Clay + (1|block), data = oreo)
mvc <- lmer(VWC ~ Mulch + Clay + CN + (1|block), data = oreo)

piecewiseSEM::psem(mol, mnitrate, mph, memf, mnit, mcat, mvc, mcn) -> psem
sol <- summary(psem)
sol$Cstat
sol$R2$Marginal[1]
# ggsem::random_layout(psem)
pol <- ggsem(psem, cols =RColorBrewer::brewer.pal(2, "Set1"), alpha = 0.05, labels = F,
      title = c("LDMC: Oreochrysum parryi", 
                paste('LDMC R2 =' ,sol$R2$Marginal[1], 
                ', p = ',sol$Cstat[1,3], ', C =', sol$Cstat[1,1])));pol
ggsave('out/piecewise_oreo_ldmc.png', bg = 'white')


sem.boot_ol <- bootEff(psem, R = 500, seed = 13, parallel = "no", ran.eff = "block")
sem_eff_ol <- semEff(sem.boot_ol)
summary(sem_eff_ol)
summary(
  semEff(sem.boot_ol, predictor = c("Mulch")),
  response = "ldmc"
)
# oreo sla ==========================================================================

mos <- lmer(sla ~  VWC  + Aspect + Cations +# Nitrifiers + #NO3_15 +# NO3_sum +# an5 +# CN + #ammonium + Sand +#+ nitrate_5_15 +
              (1|block),
            data = oreo); summary(mos); car::vif(mos)
mnit15 <- lmer(NO3_15 ~ Mulch + Biochar + Nitrifiers + Aspect + Cations + (1|block), data = oreo)
mamon <- lmer(NH4 ~ EMF + Cations + Nitrifiers + NO3_15 +(1|block),
             data = oreo, na.action = na.fail)
mph <- lmer(pH ~  Biochar + (1|block), 
            data = oreo, na.action = na.fail)
memf <- lmer(EMF ~ pH + Nitrifiers + Mulch + (1|block), 
             data = oreo, na.action = na.fail)
# mppf <- lmer(PPF ~ Biochar + DIN + pH  + TDN + (1|block), 
#              data = oreo, na.action = na.fail)
mtdn <- lmer(TDN ~ Mulch + Biochar + EMF + VWC + pH + DIN + Nitrifiers +(1|block),
             data = oreo, na.action = na.fail)
mnit <- lmer(Nitrifiers ~ Mulch + Biochar + VWC + pH + Clay + (1|block),
             data = oreo, na.action = na.fail)
mcat <- lmer(Cations ~ Biochar + VWC + Mulch + Clay + pH + (1|block), data = oreo)
mvc <- lmer(VWC ~ Mulch + Biochar + Clay  + (1|block), data = oreo)
mph <- lmer(pH ~ Biochar  + (1|block), data = oreo)


piecewiseSEM::psem(mos, mvc, mcat, mph#,# mnit
                   ) -> psem_os
sos <- summary(psem_os);sos
pos <- ggsem(psem_os, #cols = RColorBrewer::brewer.pal(3, "Set1"), 
             alpha = 0.05, labels = F,
      title = c("SLA: Oreochrysum parryi", paste('SLA R2 =' ,sos$R2$Marginal[1], 
                                                 ', p = ',sos$Cstat[1,3], ', C =', sos$Cstat[1,1])));pos
ggsave('out/piecewise_oreo_sla.png', bg = 'white')

sem.boot_os <- bootEff(psem_os, R = 500, seed = 13, parallel = "no", ran.eff = "block")
sem_eff_os <- semEff(sem.boot_os)
summary(sem_eff_os)
summary(
  semEff(sem.boot_os, predictor = c('Biochar', 'Mulch')),
  response = "sla"
)

# oreo height ==================================================================
moh <- lmer(height_cm ~  VWC + Cations + EMF + NO3+ 
              pH + 
              (1|block),
            data = oreo, na.action = na.fail)
mamon <- lmer(NO3 ~ EMF + pH + Cations + #Nitrifiers+ 
                (1|block),
              data = oreo, na.action = na.fail)
mph <- lmer(pH ~  Biochar + (1|block), 
            data = oreo, na.action = na.fail)
memf <- lmer(EMF ~ pH + Nitrifiers + Mulch + (1|block), 
             data = oreo, na.action = na.fail)
mppf <- lmer(PPF ~ Biochar + DIN + pH  + TDN + (1|block), 
             data = oreo, na.action = na.fail)
mtdn <- lmer(TDN ~ Mulch + Biochar + EMF + VWC + pH + DIN + Nitrifiers +(1|block),
             data = oreo, na.action = na.fail)
mnit <- lmer(Nitrifiers ~ Mulch + Biochar + VWC + pH + (1|block),
             data = oreo, na.action = na.fail)
mcat <- lmer(Cations ~ Biochar + VWC + Mulch + Clay + pH + (1|block), data = oreo)
mvc <- lmer(VWC ~ Mulch + Biochar + Clay  + (1|block), data = oreo)

piecewiseSEM::psem(moh, mcat, mvc, mnit,
                   memf, mamon, mph) -> psem_oh
soh <- summary(psem_oh); soh
poh <- ggsem(psem_oh, cols = RColorBrewer::brewer.pal(2, "Set1"), alpha = 0.05, labels = F, show_legend = T,
      title = c("Height: Oreochrysum parryi",paste('Height R2 =' ,soh$R2$Marginal[1], 
                                                   ', p = ',soh$Cstat[1,3], ', C =', soh$Cstat[1,1])));poh
ggsave('out/piecewise_oreo_height.png', bg = 'white')

sem.boot_oh <- bootEff(psem_oh, R = 500, seed = 13, parallel = "no", ran.eff = "block")
sem_eff_oh <- semEff(sem.boot_oh)
summary(sem_eff_oh)
summary(
  semEff(sem.boot_oh, predictor = c("Mulch")),
  response = "height.cm"
)

# halfeye plot =================================================================
# maybe bootsratp to 1000 and see how things look
efol <- getEff(sem_eff_ol, responses = 'ldmc', type = 'boot')$ldmc$Total |> as.data.frame() |>
  dplyr::select(Mulch, Biochar) %>%
  pivot_longer(cols = names(.)) |>
  mutate(response = "ldmc")
efos <- getEff(sem_eff_os, responses = 'sla', type = 'boot')$sla$Total |> as.data.frame() |>
  dplyr::select(Mulch, Biochar) %>%
  pivot_longer(cols = names(.)) |>
  mutate(response = 'sla')
efoh <- getEff(sem_eff_oh, responses = 'height.cm', type = 'boot')$height.cm$Total |> as.data.frame() |>
  dplyr::select(Mulch, Biochar) %>%
  pivot_longer(cols = names(.)) |>
  mutate(response = 'height')

bind_rows(efol, efos, efoh) |>
  ggplot() +
  ggdist::stat_dotsinterval(aes(x=value, y=name)) +
  facet_wrap(~response, ncol=1) +
  geom_vline(xintercept =0, linetype =2) +
  ggthemes::theme_clean()

efl <- getEffTable(sem_eff_ol, responses = 'ldmc') |> dplyr::filter(effect_type == 'total', predictor %in% c('Mulch', 'Biochar'))
efs <- getEffTable(sem_eff_os, responses = 'sla') |> dplyr::filter(effect_type == 'total', predictor %in% c('Mulch', 'Biochar'))
efh <- getEffTable(sem_eff_oh, responses = 'height.cm') |> dplyr::filter(effect_type == 'total', predictor %in% c('Mulch', 'Biochar'))

poe <- bind_rows(efl, efs, efh) |>
  mutate(sig = ifelse(lower_ci * upper_ci >0, "*", ""),
         response = str_to_upper(response),
         predictor = str_to_title(predictor),
         response = ifelse(response== "HEIGHT.CM", 'Height', response)) |>
  ggplot(#aes(color = sig)
         ) +
  geom_point(aes(x=effect, y=predictor)) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y= predictor), linewidth = 1, height = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(~response, ncol=1) +
  xlab("Standardized Total Effect") +
  ggtitle('Total Effects of Mulch and Biochar on Traits', 'After Accounting for the Effects of Mediators') +
  scale_color_manual(values = c('grey50', 'firebrick')) +
  ggthemes::theme_clean() +
  theme(axis.title.y = element_blank()); poe

# multipanel oreochrysum =======================================================
ggarrange(ncol=1, nrow=2,
ggarrange(pol, pos, nrow =1, ncol=2),
ggarrange(poh, poe, nrow=1, ncol=2, widths = c(1.4,1)))
ggsave(filename = 'out/piecewise_oreo_multi.png', width = 13, height = 13, bg='white')

# vaccinium sp ldmc =============

mvl <- lmer(ldmc ~ pH + Fungal_div + VWC +
              NO3 + (1|block), 
            data = vac, na.action = na.fail); summary(mvl); car::Anova(mvl)
msf <- lmer(Fungal_div ~ TWI + Mulch + Nitrifiers + Biochar + (1|block), data = vac)
vnitra <- lmer(NO3 ~ pH + EMF + Cations + Bare + Mulch + TWI+ Nitrifiers +(1|block),
             data = vac)

vph <- lmer(pH ~  Mulch  + Nitrifiers + TWI + EMF + Cations +(1|block), 
            data = vac, na.action = na.fail)

vcat <- lmer(Cations ~ Biochar + Bare + Mulch + EMF +(1|block), data = vac)

vemf <- lmer(EMF ~  Mulch + Bare + (1|block), 
             data = vac, na.action = na.fail)
vvc <- lmer(VWC ~ Mulch  + Cations + Bare + (1|block), data = vac)

vtdn <- lmer(TDN ~ Mulch + Biochar + Bare + (1|block),
             data = vac, na.action = na.fail)
vnit <- lmer(Nitrifiers ~ Biochar + VWC + EMF + #TDN + ammonium_5_15 + 
               Bare + (1|block),
             data = vac, na.action = na.fail)

piecewiseSEM::psem(mvl, vnitra, msf, vph, vcat, vvc, vemf, vnit#,vamm15, vnit, vph, vemf, vvc, vcat, vtdn, vnit15
                   ) -> psemv

svl <- summary(psemv);svl
pvl <- ggsem(psemv, cols = RColorBrewer::brewer.pal(2, "Set1"), labels = F,
      alpha = 0.05, title = c("LDMC: Vaccinium scoparium", 
                              paste('LDMC R2 =' ,svl$R2$Marginal[1], 
                                    ', p = ',svl$Cstat[1,3], ', C =', svl$Cstat[1,1])));pvl
ggsave('out/piecewise_vacc.png', bg = 'white')

sem.boot_vl <- bootEff(psemv, R = 500, seed = 13, parallel = "no", ran.eff = "block")
sem_eff_vl <- semEff(sem.boot_vl)
summary(sem_eff_vl)
summary(
  semEff(sem.boot_vl, predictor = c("Mulch")),
  response = "ldmc"
)

# vaccinium sp sla =============================================================

# ggplot(vac, aes(x=nitrate, y=nitrate_5_15, size=nitrate_5_15, color = nitrate)) +geom_point()

mvs <- lmer(sla ~ DON + EMF + NH4_15 + NH4 + #Clay +
              (1|block),
            data = vac); summary(mvs); car::vif(mvs); car::Anova(mvs); AIC(mvs)
# ggeffects::ggpredict(mvs, terms = 'ammonium_mean') |> plot()
mdm <- lmer(DON ~ Mulch + Biochar + VWC + #ammonium_5_15 +# Bare +#ammonium + 
              (1|block),
            data = vac)

mdc <- lmer(DOC ~ Mulch + Biochar + VWC + #ammonium_5_15 +# Bare +#ammonium + 
              (1|block),
            data = vac)

mtdn <- lmer(TDN ~ Mulch + Biochar +  EMF + DON   +#VWC +
              (1|block),
            data = vac)
mnit15 <- lmer(NO3_15 ~ TWI + Mulch + NH4 + NH4_15 +#Biochar +
              (1|block),
            data = vac)
mVWC <- lmer(VWC ~ Mulch + Biochar + Bare +# TWI + 
               Clay +(1|block), data = vac)
vammon <- lmer(NH4 ~ #TDN +
                 Bare + #Mulch+ 
                 Biochar + NH4_15 + DON + EMF +
                 (1|block),
               data = vac)
vamm15 <- lmer(NH4_15 ~ Bare + VWC +Mulch + DON +
                 (1|block),
               data = vac)
# vammon <- lmer(ammonium_mean ~ nitrate_5_15 + TDN + Bare + Mulch+ Biochar +# ammonium_5_15 + 
#                  (1|block),
#                data = vac)
vemf <- lmer(EMF ~  Biochar + VWC + Bare + #nitrate_5_15 + #Sand + 
               (1|block), 
             data = vac, na.action = na.fail)

piecewiseSEM::psem(mvs, vammon, vamm15, 
                   vemf, mdm, #mtdn, #mnit15,
                   mVWC) -> psemvs
svs <- summary(psemvs);svs
pvs <- ggsem(psemvs, cols = RColorBrewer::brewer.pal(2, "Set1"), labels = F,
      alpha = 0.05, title = c("SLA: Vaccinium scoparium", 
                              paste('SLA R2 =' ,svs$R2$Marginal[1], 
                                    ', p = ',svs$Cstat[1,3], ', C =', svs$Cstat[1,1])));pvs

ggsave('out/piecewise_vacc.png', bg = 'white')

sem.boot_vs <- bootEff(psemvs, R = 500, seed = 13, parallel = "no", ran.eff = "block")
sem_eff_vs <- semEff(sem.boot_vs)
summary(sem_eff_vs)
summary(
  semEff(sem.boot_vs, predictor = c("Mulch")),
  response = "sla"
)

# vaccinium height =============================================================
mvh <- lmer(height_cm ~ NH4 + NH4_15 + Nitrifiers + NO3_15+  #EMF +# DON + # Cations+
              (1|block),
            data = vac); summary(mvh); car::vif(mvh)


mniters <- lmer(Nitrifiers ~ Mulch + Biochar + VWC  + NH4_15+ #Bare +
                 (1|block),
               data = vac)

mdm <- lmer(DON ~ Mulch + Biochar + VWC + #ammonium_5_15 +# Bare +#ammonium + 
              (1|block),
            data = vac)

mdc <- lmer(DOC ~ Mulch + Biochar + VWC + #ammonium_5_15 +# Bare +#ammonium + 
              (1|block),
            data = vac)

mtdn <- lmer(TDN ~ Mulch + Biochar +  EMF + DON   +#VWC +
               (1|block),
             data = vac)
mnit15 <- lmer(NO3_15 ~ TWI + Mulch + NH4 + Bare + Nitrifiers +#NH4_15 +#Biochar +
                 (1|block),
               data = vac)
mVWC <- lmer(VWC ~ Mulch + Biochar + Bare + TWI + 
               Clay +(1|block), data = vac)
vammon <- lmer(NH4 ~ #TDN +
                 Bare + #Mulch+ 
                 #Biochar +
                 NH4_15 + DON +# EMF +
                 (1|block),
               data = vac)
vamm15 <- lmer(NH4_15 ~ Bare + VWC +Mulch + DON +
                 (1|block),
               data = vac)
vemf <- lmer(EMF ~  VWC + Bare + Nitrifiers + #nitrate_5_15 + #Sand + 
               (1|block), 
             data = vac, na.action = na.fail)

piecewiseSEM::psem(mvh, vammon, vamm15, #vemf, #mtdn, 
                   mnit15, mVWC, mniters, mdm) -> psemvh
svh <- summary(psemvh);svh
pvh <- ggsem(psemvh, cols = RColorBrewer::brewer.pal(2, "Set1"), labels = F,show_legend = T,
             alpha = 0.05, title = c("Height: Vaccinium scoparium", 
                                     paste('Height R2 =' ,svh$R2$Marginal[1], 
                                           ', p = ',svh$Cstat[1,3], ', C =', svh$Cstat[1,1])));pvh

sem.boot_vh <- bootEff(psemvh, R = 500, seed = 13, parallel = "no", ran.eff = "block")
sem_eff_vh <- semEff(sem.boot_vh)
summary(sem_eff_vh)
summary(
  semEff(sem.boot_vh, predictor = c("Mulch")),
  response = "height.cm"
)

eflv <- getEffTable(sem_eff_vl, responses = 'ldmc') |> dplyr::filter(effect_type == 'total', predictor %in% c('Mulch', 'Biochar'))
efsv <- getEffTable(sem_eff_vs, responses = 'sla') |> dplyr::filter(effect_type == 'total', predictor %in% c('Mulch', 'Biochar'))
efhv <- getEffTable(sem_eff_vh, responses = 'height.cm') |> dplyr::filter(effect_type == 'total', predictor %in% c('Mulch', 'Biochar'))

pve <- bind_rows(eflv, efsv, efhv) |>
  mutate(sig = ifelse(lower_ci * upper_ci >0, "*", ""),
         response = str_to_upper(response),
         predictor = str_to_title(predictor),
         response = ifelse(response== "HEIGHT.CM", 'Height', response)) |>
  ggplot(#aes(color = sig)
  ) +
  geom_point(aes(x=effect, y=predictor)) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y= predictor), linewidth = 1, height = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(~response, ncol=1) +
  xlab("Standardized Total Effect") +
  ggtitle('Total Effects of Mulch and Biochar on Traits', 'After Accounting for the Effects of Mediators') +
  scale_color_manual(values = c('grey50', 'firebrick')) +
  ggthemes::theme_clean() +
  theme(axis.title.y = element_blank())


# multipanel vaccinium =======================================================
ggarrange(nrow=2, ncol=1,
ggarrange(pvl, pvs, nrow=1, ncol=2),
ggarrange(pvh, pve, nrow=1, ncol=2, widths = c(1.4, 1))) 
ggsave(filename = 'out/piecewise_vac_multi.png', width = 13, height = 13, bg='white')


# make an R2 table

rv <- piecewiseSEM::rsquared(psemv) |> mutate(species = 'Vaccinium')
ro <- piecewiseSEM::rsquared(psem) |> mutate(species = "Oreochrysum")

bind_rows(rv, ro) |>
  dplyr::select(1,5,6,7) |>
  mutate_if(is.numeric, round, 2) |>
  write_csv("out/psem_r2_table.csv")



# emf sem ======================================================================
# need to do 5-15
# ggplot(d, aes(x=log(EMF +1), y=log(EMF_5_15+1))) +
#   geom_point() 
# 

mniters <- lmer(Nitrifiers ~ Mulch + Biochar + VWC + pH + TWI +
                  (1|block),
                data = d |> filter(!is.na(sla)))
mVWC <- lmer(VWC ~ Mulch + Biochar + Bare + TWI + pH + (1|block), data = d |> filter(!is.na(sla)))
mnitr <- lmer(NO3 ~ VWC + Bare  + TWI + pH + Nitrifiers + (1|block), data = d |> filter(!is.na(sla)))

memf <- lmer(EMF ~  Bare +Nitrifiers + NO3 + Mulch + pH +(1|block),
             data = d |> filter(!is.na(sla)), na.action = na.fail); summary(memf)
mph <- lmer(pH ~  Biochar + TWI + Bare + (1|block), 
            data = d |> filter(!is.na(sla)))
piecewiseSEM::psem(memf, mVWC, mniters, mnitr#, mph
                   ) -> ps_emf
summary(ps_emf)

pemf <- ggsem(ps_emf, cols = RColorBrewer::brewer.pal(2, "Set1"), labels = F,
             alpha = 0.05, title = c("EMF", "Marginal R2 = 0.56; C = 25, p = 0.197"));pemf
ggsave('out/emf_sem.png', width =8, height = 8, bg = 'white')

sem.boot_emf <- bootEff(ps_emf, R = 500, seed = 13, parallel = "no", ran.eff = "block")
sem_eff_emf <- semEff(sem.boot_emf)
summary(sem_eff_emf)
summary(
  semEff(sem.boot_emf, predictor = c("Biochar", 'Mulch')),
  response = "EMF"
)

pinset <- getEffTable(sem_eff_emf) |>
  filter(effect_type == 'total', predictor %in% c("Biochar", 'Mulch'), response == 'EMF') |>
  mutate(sig = ifelse(lower_ci * upper_ci >0, "*", ""),
         response = str_to_upper(response),
         predictor = str_to_title(predictor),
         response = ifelse(response== "HEIGHT.CM", 'Height', response)) |>
  ggplot(#aes(color = sig)
         ) +
  geom_point(aes(x=effect, y= predictor)) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y= predictor), linewidth = 1, height = 0) +
  geom_vline(xintercept = 0) +
  xlab("Total Effect") +
  scale_color_manual(values = c('grey50', 'firebrick')) +
  ggthemes::theme_clean() +
  theme(axis.title.y = element_blank(), 
        legend.position = 'none')

library(cowplot)

ggdraw(pemf) +
  draw_plot(pinset, x=.6, y=0, width = .4, height = .2)
ggsave('out/emf_sem.png', width = 6, height = 6, bg='white')
