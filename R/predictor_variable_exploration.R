# explore predictor variables

source("R/a_church_data_prep.R")

# treatment effects
sites_w_23_soil |>
  dplyr::select(-aspect, -elevation, -flowdir, -folded_aspect,
                -latitude, -longitude, -twi, -TPI, -TRI,
                -slope,
                -starts_with("percent")) |>
  pivot_longer(-c(block, plot, treatment, sample_id)) |>
  ggplot(aes(x=treatment, y = value, fill=treatment)) +
  geom_boxplot(outliers = F) +
  facet_wrap(~name, scales = "free_y")  +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())

ggsave("out/predictor_variables.png", width=14, height=10, bg="white")

# block effects?
sites_w_23_soil |>
  # dplyr::select(-aspect, -elevation, -flowdir, -folded_aspect,
  #               -latitude, -longitude, -twi, -TPI, -TRI,
  #               -slope,
  #               -starts_with("percent")) |>
  pivot_longer(-c(block, plot, treatment, sample_id)) |>
  ggplot(aes(x=block, y = value)) +
  geom_boxplot() +
  facet_wrap(~name, scales = "free_y") #+
  # theme(axis.text.y = element_blank(), 
  #       axis.ticks.y = element_blank()) +
  # scale_y_log10()


funguild |>
  pivot_longer(-c(block, treatment, depth)) |>
  filter(!name %in% c("unassigned", "amf")) |>
  ggplot(aes(x=treatment, y=value)) +
  geom_boxplot(outliers=F) +
  # geom_jitter() +
  facet_grid(name~ depth, scales = "free")
