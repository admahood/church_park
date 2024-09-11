# piecewise sem
library(lmerTest)
library(piecewiseSEM)
library(ggsem)
library(ggeffects)
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
  # mutate_if(is.numeric, datawizard::standardise) |>
  mutate_if(is.character, as.factor) |>
  left_join(nitrifiers) |>
  left_join(pcaf) |>
  left_join(nitrifiers) |>
  left_join(discriminant_taxa) |>
  left_join(soil_23_tf) |>
  mutate(din_tf = ammonium + nitrate) |>
  dplyr::rename(DOC = dissolved_organic_carbon,
                PPF = plant_pathogen_5,
                Aspect = folded_aspect,
                DIN = din_tf,
                TDN = dissolved_nitrogen,
                EMF = ectomycorrhizae_5,
                pH = soil_p_h,
                cations = cation_sum,
                sand = percent_sand_0to5)

oreo <- filter(d, species_full == "Oreochrysum parryi"); summary(oreo)
vac <- filter(d, species_full == "Vaccinium sp"); summary(vac)
car <- filter(d, species_full == "Carex sp") |>
  filter(!is.na(sla)); summary(car)

# lmers for the traits

# oreochrysum parryi ========================
library(MuMIn)
mol <- lmer(ldmc ~ Aspect + DIN +# EMF +#pH + EMF + DOC + PPF + TDN + 
               (1|block),
            data = oreo, na.action = na.fail)
mos <- lmer(sla ~ Aspect + ldmc+ cations + vwc +
              (1|block),
            data = oreo, na.action = na.fail)

mdin <- lmer(DIN ~ vwc + pH + nitrifiers + EMF + (1|block),
             data = oreo, na.action = na.fail)

mph <- lmer(pH ~  nitrifiers + cations + mulch +(1|block), 
            data = oreo, na.action = na.fail)
memf <- lmer(EMF ~ pH + nitrifiers + mulch + (1|block), 
             data = oreo, na.action = na.fail)
mppf <- lmer(PPF ~ biochar + DIN + pH  + TDN + (1|block), 
             data = oreo, na.action = na.fail)
mtdn <- lmer(TDN ~ mulch + biochar + EMF + vwc + pH + DIN + nitrifiers +(1|block),
             data = oreo, na.action = na.fail)
mnit <- lmer(nitrifiers ~ mulch + biochar + vwc +(1|block),
             data = oreo, na.action = na.fail)
mcat <- lmer(cations ~ biochar + vwc + mulch + (1|block), data = oreo)
mvc <- lmer(vwc ~ mulch + biochar + (1|block), data = oreo)

piecewiseSEM::psem(mol, mdin, mph, memf, mnit, mos, mcat, mvc) -> psem
summary(psem)
ggsem(psem, cols = c("black", "gold"), alpha = 0.05, 
      title = c("LDMC: Oreochrysum parryi", 
                'R2 = 0.28'))
ggsave('out/piecewise_oreo_ldmc.png', bg = 'white')

# vaccinium sp =============

mvl <- lmer(ldmc ~ mulch  + twi + pH + sand + DIN +(1|block), 
            data = vac, na.action = na.fail)
summary(mvl)

mvs <- lmer(sla ~ ldmc+ 
              (1|block),
            data = vac)

vdin <- lmer(DIN ~ pH + EMF + cations + (1|block),
             data = vac)
summary(vdin)

vph <- lmer(pH ~  nitrifiers + vwc  + nitrifiers + twi + cations + (1|block), 
            data = vac, na.action = na.fail)


vemf <- lmer(EMF ~ pH + nitrifiers + mulch + (1|block), 
             data = vac, na.action = na.fail)

vvc <- lmer(vwc ~ mulch + sand + cations + (1|block), data = vac)

vppf <- lmer(PPF ~ biochar + vwc + EMF +(1|block), 
             data = vac, na.action = na.fail)

vtdn <- lmer(TDN ~ mulch + biochar + EMF + vwc + pH + DIN + nitrifiers +(1|block),
             data = vac, na.action = na.fail)

vcat <- lmer(cations ~ biochar + (1|block), data = vac)

vnit <- lmer(nitrifiers ~ mulch + biochar + vwc + (1|block),
             data = vac, na.action = na.fail)

piecewiseSEM::psem(mvl, vdin, vnit, vph, mvs, vemf, vvc, vcat) -> psemv
summary(psemv)
ggsem(psemv, cols = c("black", "gold"), alpha = 0.05, title = c("LDMC: Vaccinium sp", 'R2 = 0.16'))
ggsave('out/piecewise_vac_ldmc.png', bg = 'white')

# Carex sp =============

mcl <- lmer(ldmc ~  vwc + total_n + DOC + TDN + (1|block), 
            data = car, na.action = na.fail)
summary(mcl)

mcs <- lmer(sla ~ DIN + sand + cations + total_c +
              (1|block),
            data = car)

cdin <- lmer(DIN ~  pH + EMF + cations + vwc + nitrifiers + (1|block),
             data = car)
summary(cdin)

cph <- lmer(pH ~  nitrifiers + vwc  + nitrifiers + twi + (1|block), 
            data = car, na.action = na.fail)


cemf <- lmer(EMF ~ pH + nitrifiers + mulch + (1|block), 
             data = car, na.action = na.fail)


cppf <- lmer(PPF ~ biochar + vwc + EMF +(1|block), 
             data = car, na.action = na.fail)

ctdn <- lmer(TDN ~ mulch + biochar + EMF + vwc + pH + DIN + nitrifiers +(1|block),
             data = car, na.action = na.fail)

cvc <- lmer(vwc ~ mulch + biochar + sand + cations + (1|block), data = car)

cnit <- lmer(nitrifiers ~ mulch + biochar + vwc  +(1|block),
             data = car, na.action = na.fail)

piecewiseSEM::psem(mcl, cdin, cnit, cph, mcs, cemf, cvc) -> psemc
summary(psemc)
ggsem(psemc, cols = c("black", "gold"), alpha = 0.05, title = c("LDMC: Carex sp", 'R2 = 0.05'))
ggsave('out/piecewise_car_ldmc.png', bg = 'white')
