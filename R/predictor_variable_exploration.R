# explore predictor variables

source("R/a_church_data_prep.R")

# treatment effects
sites_w_23_soil |>
  dplyr::select(-aspect, -elevation, -flowdir, -folded_aspect,
                -latitude, -longitude, -twi, -TPI, -TRI,
                -slope, -starts_with(c("amf", "shannon", "simpson", "species", "unassigned")),
                -starts_with(c("percent", "other", "plant_pathogen", "ecto"))) |>
  pivot_longer(-c(block, plot, treatment, sample_id)) |>
  ggplot(aes(x=treatment, y = value, fill=treatment)) +
  geom_boxplot(outliers = F) +
  facet_wrap(~name, scales = "free_y")  +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())

ggsave("out/soil_variables.png", width=9, height=6, bg="white")

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


# common plants =======

comm_both |>
  as_tibble(rownames = "row") |>
  mutate(block = str_sub(row, 1,2),
         year = str_extract(row, "\\d{4}"),
         treatment = ifelse(nchar(row) == 7, str_sub(row, 3,3), str_sub(row, 3,4))) |>
  dplyr::select(block, year, treatment, Oreochrysum_parryi, Vaccinium_scoparium, Carex_sp) |>
  pivot_longer(-c(block, year, treatment)) |>
  mutate(name = str_replace_all(name, "_", " ")) |>
  ggplot(aes(x=treatment, y=value, fill = year)) +
  geom_boxplot() +
  facet_wrap(~name) +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(1,0.95),
        legend.justification = c(1,1),
        legend.background = element_rect(color = "black", fill=NA))

ggsave("out/common_plants.png", width = 7, height = 3.5, bg="white")


common_plant_differences <- comm_both |>
  as_tibble(rownames = "row") |>
  mutate(block = str_sub(row, 1,2),
         year = str_extract(row, "\\d{4}"),
         treatment = ifelse(nchar(row) == 7, str_sub(row, 3,3), str_sub(row, 3,4))) |>
  dplyr::select(block, year, treatment, Oreochrysum_parryi, Vaccinium_scoparium, Carex_sp) |>
  pivot_longer(-c(block, year, treatment)) |>
  pivot_wider(names_from = year, values_from = value) |>
  mutate(diff = `2023` - `2016`,
         treatment = str_replace_all(treatment, "c", "0"))


ts <- data.frame(treatment = NA, name = NA, p = NA, ypos = NA)
counter <- 1
for(i in unique(common_plant_differences$name)){
  for(t in unique(common_plant_differences$treatment)){
    p <- t.test(x = common_plant_differences |> filter(name == i, treatment == t) |> pull(diff))$p.value
    ts[counter,1] <- t
    ts[counter,2] <- i
    ts[counter,3] <- p
    ts[counter,4] <- common_plant_differences |> filter(name == i, treatment == t) |> pull(diff) |> max()
    
    counter <- counter + 1
    
  }
}
sigs <- filter(ts, p < 0.05)
ggplot(common_plant_differences, aes(y=diff, fill=treatment, x=name)) +
  geom_boxplot(outliers = F) +
  geom_hline(yintercept = 0) +
  geom_point(data = sigs, y=c(4, 9, 7), x = c(1.9, 1.72, 2.28), shape=8,
             size=3, stroke = 1, show.legend = F) +
  ylab("Change in Percent Cover") +
  theme_bw() +
  theme(axis.title.x= element_blank()) 

ggsave("out/common_plants.png", width = 7, height = 3.5, bg="white")
