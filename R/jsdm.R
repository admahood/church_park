# Church Park vegetation composition analysis
source("R/a_church_data_prep.R")
library(brms)

# brms jsdm ====================================================================
bformula <- formula(str_c("mvbind(Carex_sp, Oreochrysum_parryi,",
                          " Vaccinium_scoparium)",
                          " ~ treatment + (1|block)"
))

# gaussian (with rescor)
to <- Sys.time()
gmd <- brm(bformula, # + set_rescor(rescor = TRUE), 
           data = cbind(sites_w_23_soil, comm/100),
           iter = 4000,
           family = "gaussian"
           )
t1 <- Sys.time()
print(t1-to)

summary(gmd)
# with random effects
broom.mixed::tidy(gmd) %>%
  filter(str_sub(term, 1,2) != "sd", term != "(Intercept)") %>% 
  mutate(term = str_remove_all(term, "treatment")) %>%
  ggplot(aes(x=estimate, y=response)) +
  geom_vline(xintercept = 0, lty=2) +
  geom_point() +
  facet_wrap(~term, scales="free_x", nrow=1) +
  geom_segment(aes(x = conf.low, xend=conf.high, y=response, yend=response), lwd=1) +
  ggtitle("gaussian model")

summary(gmd)$rescor_pars %>%
  tibble::rownames_to_column("name") |>
  janitor::clean_names() |>
  ggplot(aes(x=estimate, y=name)) +
  geom_vline(xintercept = 0, lty=2) +
  geom_point() +
  geom_segment(aes(x = l_95_percent_ci , xend=u_95_percent_ci , y=name, yend=name), lwd=1) +
  ggtitle("gaussian model")


# zero inflated beta
t2 <- Sys.time()
zbmd1 <- brm(bformula, # + set_rescor(rescor = TRUE), 
            data = cbind(sites_w_23_soil, comm/100),
            iter = 4000,
            family = "zero_inflated_beta")
t3 <- Sys.time()
print(t3-t2)

summary(zbmd1)

broom.mixed::tidy(zbmd1) %>%
  filter(str_sub(term, 1,2) != "sd", term != "(Intercept)") %>%
  mutate(term = str_remove_all(term, "treatment")) %>%
  ggplot(aes(x=estimate, y=response)) +
  geom_vline(xintercept = 0, lty=2) +
  geom_point() +
  facet_wrap(~term, scales="free_x", nrow=1) +
  geom_segment(aes(x = conf.low, xend=conf.high, y=response, yend=response), lwd=1) +
  ggtitle("zero-inflated beta model")




library(gllvm)

mod <- gllvm(y = comm/100, 
      family = "beta",
      method = "LA")

plot(mod)
ordiplot(mod, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), 
         ylim = c(-3, 3), 
         main = "No predictors")
# with preds
mod1 <- gllvm(y = comm/100, 
             X = sites_w_23_soil,
             formula = ~ treatment,
             row.eff = ~ (1|block),
             family = "beta",
             method = "LA")

plot(mod1)
ordiplot(mod1, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), 
         ylim = c(-3, 3), 
         main = "w treatment")
coefplot(mod1)
AIC(mod, mod1)
# sjsdm? ===================