# gllvm

source("R/a_church_data_prep.R")
library(gllvm)
y=as.matrix(comm/100)

mod0 <- gllvm(y = y, 
              family = "beta", method = "EVA")

mod <- gllvm(y = y, X = sites_w_23_soil, 
             formula = ~ treatment,
             family = "beta", method = "EVA")

mod1 <- gllvm(y = y, X = sites_w_23_soil, 
             formula = ~ treatment + percent_sand_0to5 + n + c+ mean_vwc + ectomycorrhizae_5,
             family = "beta", method = "EVA")

mod2 <- gllvm(y = y, X = sites_w_23_soil, 
              formula = ~  percent_sand_0to5 + n + c + mean_vwc,
              family = "beta", method = "EVA")

modn <- gllvm(y = nitrifier_binary_matrix_5, X = sites_w_23_soil, 
              formula = ~  mean_vwc + treatment + c + n + shannon_fungi_5,
              family = "binomial")
summary(modn)
MuMIn::sw(list(mod, mod1, mod2))
getResidualCov(mod, adjust = 0)$trace
getResidualCov(mod0, adjust = 0)$trace
getResidualCov(mod1, adjust = 0)$trace
getResidualCov(mod2, adjust = 0)$trace

AIC(mod, mod0, mod1, mod2)

summary(mod1)
plot(mod)
confint(mod)
ordiplot(mod)
# coefplot(mod)

MuMIn::dredge(mod)

means <- coefficients(mod)
cis <- confint(mod)


# ggplot(means$Species.scores |> as_tibble(rownames = 'name'), aes(LV1, LV2)) +
#   geom_text(aes(label=name))

sigs <- cis |>
  as_tibble(rownames="name") |>
  filter(str_sub(name, 1,5) == "Xcoef") |>
  mutate(sig = ifelse(`2.5 %` * `97.5 %` > 0, "yes", "no"),
         name = str_remove_all(name, "Xcoef.")) |>
  tidyr::separate(name, into = c("name", "species"),sep = "\\:") 

means$Xcoef |> 
  as_tibble(rownames = 'species') |> 
  pivot_longer(-species) |>
  left_join(sigs) |>
  # group_by(name) |>
  # mutate(value = scale(value, center = F)) |>
  # ungroup() |>
  ggplot(aes(x=name, y = species)) +
  geom_tile(aes(fill = value, alpha = sig)) +
  geom_tile(aes(color = sig), fill = "transparent") +
  scale_fill_gradient2() +
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_color_manual(values = c("transparent", "black")) +
  theme(axis.text.x = element_text(angle =45, hjust=1))



gllvm::getResidualCor(mod) |> 
  ggcorrplot::ggcorrplot(type = "lower", hc.order = T)

gllvm::getResidualCov(mod)$cov |> 
  ggcorrplot::ggcorrplot(type = "lower", hc.order = T)
