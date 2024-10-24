# piecewise sem
library(lmerTest)
library(piecewiseSEM)
library(ggsem)
library(ggeffects)
# library(semEff)
library(MuMIn)
source("R/a_church_data_prep.R")
source("R/a_itv_data_prep.R")
library(lavaan)
library(ggsem) # devtools::install_github("admahood/ggsem")
library(ggtext)
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
                sand = percent_sand_0to5) |>
  mutate_if(is.character, as.factor)

oreo <- filter(d, species_full == "Oreochrysum parryi"); summary(oreo)
vac <- filter(d, species_full == "Vaccinium sp"); summary(vac)
car <- filter(d, species_full == "Carex sp") |>
  filter(!is.na(sla)); summary(car)

# lmers for the traits

# oreochrysum parryi ========================
mol <- lmer(ldmc ~ DIN + Aspect + # sand +#twi +# EMF +#pH + EMF + DOC + PPF + TDN + 
               (1|block),
            data = oreo, na.action = na.fail); summary(mol)
# mol <- lm(ldmc ~ DIN + Aspect, data = oreo)
mos <- lmer(sla ~  ldmc + vwc + cations + Aspect +
              (1|block),
            data = oreo, na.action = na.fail)
# mos <- lm(sla ~  ldmc + vwc + cations + Aspect, data = oreo)
mdin <- lmer(DIN ~ vwc + pH + nitrifiers + EMF + (1|block),
             data = oreo, na.action = na.fail)
mph <- lmer(pH ~  nitrifiers + cations + mulch + (1|block), 
            data = oreo, na.action = na.fail)
memf <- lmer(EMF ~ pH + nitrifiers + mulch + (1|block), 
             data = oreo, na.action = na.fail)
mppf <- lmer(PPF ~ biochar + DIN + pH  + TDN + (1|block), 
             data = oreo, na.action = na.fail)
mtdn <- lmer(TDN ~ mulch + biochar + EMF + vwc + pH + DIN + nitrifiers +(1|block),
             data = oreo, na.action = na.fail)
mnit <- lmer(nitrifiers ~ mulch + biochar + vwc +(1|block),
             data = oreo, na.action = na.fail)
mcat <- lmer(cations ~ biochar + vwc + mulch + sand + (1|block), data = oreo)
mvc <- lmer(vwc ~ mulch + biochar + sand + (1|block), data = oreo)

lapply(c(mdin, mph, memf, mppf, mtdn, mnit, mcat), performance::r2)
piecewiseSEM::psem(mol, mdin, mph, memf, mnit, mos, mcat, mvc) -> psem

summary(psem)
RColorBrewer::display.brewer.all()
ggsem(psem, cols =RColorBrewer::brewer.pal(2, "Set1"), alpha = 0.05, 
      title = c("LDMC: Oreochrysum parryi", 
                'LDMC R2 = 0.28, SLA R2 = 0.55'))
ggsave('out/piecewise_oreo_ldmc.png', bg = 'white')

# vaccinium sp =============

mvl <- lmer(ldmc ~ mulch  + pH + DIN +(1|block), 
            data = vac, na.action = na.fail); summary(mvl); car::Anova(mvl)
mvs <- lmer(sla ~ ldmc+ 
              (1|block),
            data = vac)
vdin <- lmer(DIN ~ pH + EMF + cations + (1|block),
             data = vac)
vph <- lmer(pH ~  mulch  + nitrifiers + twi + cations + (1|block), 
            data = vac, na.action = na.fail)
vemf <- lmer(EMF ~ pH + nitrifiers + mulch + TDN +(1|block), 
             data = vac, na.action = na.fail)
vvc <- lmer(vwc ~ mulch + sand + cations + (1|block), data = vac)
vppf <- lmer(PPF ~ biochar + vwc + EMF +(1|block), 
             data = vac, na.action = na.fail)
vtdn <- lmer(TDN ~ mulch + biochar + EMF + vwc + pH + DIN + nitrifiers +(1|block),
             data = vac, na.action = na.fail)
vcat <- lmer(cations ~ biochar + TDN + (1|block), data = vac)
vnit <- lmer(nitrifiers ~ biochar + vwc + TDN + (1|block),
             data = vac, na.action = na.fail)

piecewiseSEM::psem(mvl, vdin, vnit, vph, vemf, vvc, vcat) -> psemv
summary(psemv)
ggsem(psemv, cols = c("black", "gold"), alpha = 0.05, title = c("LDMC: Vaccinium sp", 'LDMC R2 = 0.09'))
ggsave('out/piecewise_vacc.png', bg = 'white')

# make an R2 table

rv <- piecewiseSEM::rsquared(psemv) |> mutate(species = 'Vaccinium')
ro <- piecewiseSEM::rsquared(psem) |> mutate(species = "Oreochrysum")

bind_rows(rv, ro) |>
  dplyr::select(1,5,6,7) |>
  mutate_if(is.numeric, round, 2) |>
  write_csv("out/psem_r2_table.csv")

# Carex sp =============

mcl <- lmer(ldmc ~  vwc + total_n + DOC + TDN + (1|block), 
            data = car, na.action = na.fail)
summary(mcl)

# try a glmer with a gamma distribution
mcs <- lmer(sla ~ DIN + 
              (1|block),
            data = car); summary(mcs)

library(brms)
bcs <- brm(sla ~ DIN + 
              (1|block), family = 'gamma', 
            data = car); summary(bcs)
conditional_effects(bcs)
performance::r2_bayes(bcs)
cdin <- lmer(DIN ~  pH + EMF + cations + vwc + (1|block),
             data = car)
summary(cdin)

cph <- lmer(pH ~  nitrifiers + vwc  + nitrifiers + twi + (1|block), 
            data = car, na.action = na.fail)


cemf <- lmer(EMF ~ pH + nitrifiers + mulch + (1|block), 
             data = car, na.action = na.fail)


cppf <- lmer(PPF ~ biochar + vwc + EMF +(1|block), 
             data = car, na.action = na.fail)

ctdn <- lmer(TDN ~ mulch + biochar + EMF + pH + nitrifiers + cations + (1|block),
             data = car, na.action = na.fail)

cvc <- lmer(vwc ~ mulch + biochar + sand + cations + (1|block), data = car)

cnit <- lmer(nitrifiers ~ mulch + biochar + vwc +(1|block),
             data = car, na.action = na.fail)
ccat <- lmer(cations ~ mulch + biochar + (1|block), data = car)


piecewiseSEM::psem(mcs, cdin, cnit, cph,  cemf, cvc) -> psemc
summary(psemc)
ggsem(psemc, cols = c("black", "gold"), 
      alpha = 0.05, title = c("LDMC: Carex sp", 'R2 = 0.05'))
ggsave('out/piecewise_car_ldmc.png', bg = 'white')
