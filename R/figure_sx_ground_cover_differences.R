#figure sx ground cover differences
library(tidyverse)
source("R/a_church_data_prep.R")
# ground cover =================================================================

gc23 <- readxl::read_xlsx("data/church_park_cover.xlsx", sheet = "ground_cover") |>
  dplyr::filter(category != "Basal Total") |>
  dplyr::mutate(category = str_to_lower(category) |> 
                  str_replace_all(" ", "_") |> 
                  str_replace_all("/", "_") |> 
                  str_remove_all("\\(>") |>
                  str_remove_all("\\)") |>
                  str_remove_all("&") |>
                  str_remove_all("_$")) |>
  pivot_longer(-category) |>
  filter(!is.na(value)) |>
  tidyr::separate(name, into = c("block", "plot", "quadrat")) |>
  mutate(plot = str_c(block, plot)) |>
  left_join(sites) |>
  pivot_wider(names_from = category, values_from = value, values_fill = 0) |>
  dplyr::mutate(biochar = biochar + (biochar__wood_mulch/2),
                wood_mulch = wood_mulch + (biochar__wood_mulch/2)) |>
  dplyr::select(-biochar__wood_mulch, -plot) |>
  dplyr::rename(cwd_fine = `cwd_1,10,100_hr`,
                cwd_coarse =  cwd_1000_hr) |>
  dplyr::mutate(sample_year = 2023); print(gc23)


gc16 <- readxl::read_xlsx("data/Church Park Botany 2016_ccr.xlsx", sheet = "ground_cover_2016") |>
  janitor::clean_names() |>
  dplyr::select(-date) |>
  dplyr::mutate(treatment = case_match(plot,
                                       "Biochar" ~ "b",
                                       "Biochar + Mulch" ~ "bm",
                                       "Control" ~ "c",
                                       "Mulch" ~ "m"), 
                quadrat = str_c("q", quadrat),
                block = str_c("b", block)) |>
  dplyr::select(-plot, -moss_lichen) |>
  dplyr::rename(root_stump = other_root_stump, cwd_coarse = cwd_course) |>
  dplyr::mutate(sample_year = 2016);gc16

gc_diff <- dplyr::bind_rows(gc16, gc23) |>
  pivot_longer(-c(block, quadrat, sample_year, treatment)) |>
  filter(!is.na(value)) |>
  group_by(block, treatment, sample_year, name) |>
  summarise(value = sum(value)/2) |>
  ungroup() |>
  pivot_wider(names_from = sample_year, values_from = value, values_fill = 0) |>
  mutate(diff = `2023` - `2016`) 

pgcd <- gc_diff |>
  filter(!name %in% c("root_stump", "cwd_coarse", "cwd_fine", "basal_live_plant", 'Tree')) |>
  pivot_longer(-c(block, treatment, name), names_to = 'year') |>
  filter(year != "diff") |>
  mutate(treatment = str_replace_all(treatment, "c", "0"),
         name = ifelse(name == "bare_soil_gravel", "Bare Ground", name)) |>
  mutate(treatment = str_to_upper(treatment) |> str_replace_all("0", "CTL"),
         treatment = fct_relevel(treatment, 'CTL', 'B', "BM", "M"),
         name = str_replace_all(name, "_", " ") |> 
           str_to_title() |> 
           str_remove_all(c(" 1cm")) |>
           str_remove_all(" Duff") |>
           str_remove_all("Wood "))|>
  ggplot(aes(y=value, x=year, fill = treatment)) +
  geom_boxplot() +
  ylab('Percent Cover') +
  facet_grid(name~treatment, scales = 'free') +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position ='none',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.justification = c(1,1),
        legend.background = element_rect(color = "black"))

pvd <- veg_diff |>
  bind_rows(oreo_vac |> dplyr::rename(diff = cover_change)) |>
  filter(!name %in% c("root_stump", "cwd_coarse", "cwd_fine", "basal_live_plant", 'Tree')) |>
  pivot_longer(-c(block, treatment, name), names_to = 'year') |>
  filter(year != "diff") |>
  mutate(treatment = str_replace_all(treatment, "c", "0"))|>
  mutate(treatment = str_to_upper(treatment) |> str_replace_all("0", "CTL"),
         treatment = fct_relevel(treatment, 'CTL', 'B', "BM", "M"),
         name = fct_relevel(name, 'Graminoid', 'Forb', "Shrub", "Invasive", "Oreochrysum", "Vaccinium"))|>
  ggplot(aes(y=value, x=year, fill = treatment)) +
  geom_boxplot() +
  ylab('Percent Cover') +
  facet_grid(name~treatment, scales = 'free') +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position ='none',
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.justification = c(1,1),
        legend.background = element_rect(color = "black"))

ggarrange(pgcd, pvd, nrow =2, heights = c(5,6))

ggsave("out/figure_s1_ground_cover_differences.png", bg='white', width=6, height=11)
