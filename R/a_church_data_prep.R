# church park data prep
library(tidyverse)
library(janitor)

soil_23 <- read_csv('data/cp_soil_nutrients_almModified.xlsx - Sheet1.csv') |>
  janitor::clean_names() |>
  dplyr::select(-ends_with('mg_l'), -sample_num, -plot_id) |>
  mutate(sample_id = str_replace_all(sample_id, "B_M", "BM"),
         depth = ifelse(depth == "0-5", "5", "15"),
         treatment = str_to_lower(trt) |> as.factor(),
         block = str_c("b", blk),
         depth = as.factor(depth)) |>
  dplyr::select(-blk, -trt)  ; glimpse(soil_23)

vwc_23 <- read_csv("data/VWC_22Jun.csv") |>
  dplyr::mutate(treatment = case_match(Trt,
                                       "Biochar" ~ "b",
                                       "B+M" ~ "bm",
                                       "Ctl" ~ "c",
                                       "Mulch" ~ "m"    
  ),block = str_c("b", Block)) |>
  dplyr::select(block, treatment, mean_vwc = Mean); vwc_23

soil_texture_23 <- read_csv("data/2023_ChurchPark_SoilTexture.csv") |>
  janitor::clean_names() |>
  dplyr::select(block = blk, starts_with("percent"), depth) |>
  dplyr::mutate(block = str_c("b", block), 
                depth = str_replace_all(depth,"-", "to")) |>
  pivot_wider(id_cols = block, names_from = depth,
              values_from = c(percent_clay, 
                              percent_silt_plus_clay, 
                              percent_silt, 
                              percent_sand)); soil_texture_23
  
soil_0_5in <- soil_23 |> filter(depth == "5") |> dplyr::select(-depth)

# microbial diversity 

fungi_div <- readxl::read_xlsx("data/Church_Park_Microbial_Data.xlsx", sheet = 3) |>
  janitor::clean_names() |>
  mutate(sample_id = str_replace_all(sample_id, "B_M", "BM")) |>
  tidyr::separate(sample_id,sep = "_", into = c("study", "block", "treatment", "plot", "depth")) |>
  mutate(treatment = str_to_lower(treatment) |> str_replace_all("ctl", "c"), 
         block = str_to_lower(block)) |>
  dplyr::select(-study, -simpson_5) |>
  dplyr::rename(simpson = simpson_4) |>
  filter(!is.na(treatment), depth != '15.2') |>
  pivot_wider(names_from = depth, 
              values_from =c("shannon", "species_richness", "simpson"),
              names_prefix = "fungi_")|>
  mutate(bt = str_c(block, treatment)) |>
  dplyr::select(-block, -treatment, -plot); fungi_div

bacteria_div <- readxl::read_xlsx("data/Church_Park_Microbial_Data.xlsx", sheet = 5) |>
  janitor::clean_names() |>
  mutate(sample_id = str_replace_all(sample_id, "B_M", "BM")) |>
  tidyr::separate(sample_id,sep = "_", into = c("study", "block", "treatment", "plot", "depth")) |>
  mutate(treatment = str_to_lower(treatment) |> str_replace_all("ctl", "c"), 
         block = str_to_lower(block)) |>
  dplyr::select(-study, -simpson_5) |>
  dplyr::rename(simpson = simpson_4) |>
  filter(!is.na(treatment), depth != '15.2') |>
  pivot_wider(names_from = depth, 
              values_from =c("shannon", "species_richness", "simpson"),
              names_prefix = "bacteria_")|>
  mutate(bt = str_c(block, treatment))|>
  mutate(bt = str_c(block, treatment)) |>
  dplyr::select(-block, -treatment, -plot); bacteria_div

sites <- readxl::read_xlsx("data/church_park_cover.xlsx",sheet = "site_info") %>%
  mutate(plot = str_c("b", block, plot), treatment = treatment,
         block = str_c("b", block)) 

sites_w_23_soil <- 
  sites |>
  left_join(vwc_23) |>
  left_join(soil_0_5in) |>
  left_join(soil_texture_23) |>
  mutate(bt = str_c(block, treatment)) |>
  left_join(bacteria_div) |>
  left_join(fungi_div) |>
  tibble::column_to_rownames("bt") |>
  mutate(treatment = treatment |> str_replace_all("c", "0"))

# veg community 
comm_16_long <- readxl::read_xlsx("data/Church Park Botany 2016_ccr.xlsx", sheet = "botany2016") |>
  dplyr::mutate(treatment = case_match(Plot,
    "Biochar" ~ "b",
    "Biochar + Mulch" ~ "bm",
    "Control" ~ "c",
    "Mulch" ~ "m"    
  ), Quadrat = str_c("q", Quadrat),
  Block = str_c("b", Block)) |>
  dplyr::select(-Plot)  %>%
  pivot_longer(cols = names(.)[4:ncol(.)-1],values_drop_na = T, names_to = "species") |>
  dplyr::rename(block=Block, quadrat=Quadrat);glimpse(comm_16_long)

comm_16_wide <- comm_16_long %>%
  group_by(block, treatment, species) %>%
  summarise(cover = sum(value)/2) %>% # just doing mean does not work
  ungroup() %>%
  mutate(species = str_remove_all(species,"\\.") %>% str_replace_all(" ", "_")) %>%
  pivot_wider(id_cols = c(block, treatment), 
              names_from = species, 
              values_fill = 0,
              values_from = cover); comm_16_wide 

comm16 <- comm_16_wide %>%
  mutate(row = str_c(block, treatment)) %>%
  dplyr::select(-block, -treatment) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("row"); comm16 

comm_long <- readxl::read_xlsx("data/church_park_cover.xlsx") %>%
  pivot_longer(cols = names(.)[2:length(.)],values_drop_na = T) %>%
  tidyr::separate(name, c("block", "plot", "quadrat"), sep = "_") %>%
  dplyr::mutate(plot = str_c(block, plot) |> as.factor()) |>
  left_join(sites) |>
  dplyr::select(-plot) |>
  mutate(species = str_remove_all(species,"\\.") %>% str_replace_all(" ", "_")) |>
  mutate(species = ifelse(species == "Achillea_millifolium", "Achillea_millefolium", species),
         species = ifelse(species == "Lodgepole", "Pinus_contorta", species),
         species = ifelse(str_sub(species, 1,10)== "Cirsium_sp", "Cirsium_sp", species));comm_long

comm_wide <- comm_long %>%
  group_by(block, treatment, species) %>%
  summarise(cover = sum(value)/2) %>% # just doing mean does not work
  ungroup() %>%
  pivot_wider(id_cols = c(block, treatment), 
              names_from = species, 
              values_fill = 0,
              values_from = cover); comm_wide 

comm <- comm_wide %>%
  mutate(row = str_c(block, treatment)) %>%
  dplyr::select(-block, -treatment) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("row"); comm 

# both years joined
comm_both <- comm_long %>%
  mutate(year = "2023") |>
  bind_rows(comm_16_long %>%
              mutate(year = "2016")) |>
  mutate(species = ifelse(str_sub(species,1,5) == "Collo", "Other_Forbs", species),
         species = ifelse(str_sub(species,1,5) == "Epilo", "Other_Forbs", species),
         species = ifelse(str_sub(species,1,5) == "Gayop", "Other_Forbs", species),
         species = ifelse(str_sub(species,1,5) == "Nasel", "Other_Graminoids", species),
         species = ifelse(str_sub(species,1,5) == "Elymu", "Elymus_spp", species),
         species = ifelse(str_sub(species,1,5) == "Pachy", "Other_Shrubs", species),
         species = ifelse(str_sub(species,1,5) == "Cheno", "Other_Forbs", species),
         species = ifelse(str_sub(species,1,5) == "Madia", "Other_Forbs", species),
         species = ifelse(str_sub(species,1,5) == "Anaph", "Other_Forbs", species),
         species = ifelse(str_sub(species,1,5) == "Hiera", "Other_Forbs", species),
         species = ifelse(str_sub(species,1,5) == "Hacke", "Other_Forbs", species),
         species = ifelse(str_sub(species,1,5) == "Polyg", "Other_Forbs", species),
         species = ifelse(str_sub(species,1,5) == "Boech", "Other_Forbs", species),
         species = ifelse(str_sub(species,1,5) == "Penst", "Other_Forbs", species),
         species = ifelse(str_sub(species,1,5) == "Sperg", "Other_Forbs", species),
         species = ifelse(str_sub(species,1,5) == "Agose", "Other_Forbs", species),
         species = ifelse(str_sub(species,1,5) == "Senec", "Senecio_spp", species),
         species = ifelse(str_sub(species,1,5) == "Arnic", "Arnica_spp", species),
         species = ifelse(str_sub(species,1,5) == "Bromu", "Other_Graminoids", species),
         species = ifelse(species == "FB Tooth Top", "Other_Forbs", species),
         species = ifelse(species == "Unk_grass", "Other_Graminoids", species),
         species = ifelse(species == "Helianthella_sp", "Other_Forbs", species),
         species = ifelse(species == "Nasella viridulis", "Other_Graminoids", species),
  ) |>
  group_by(block, treatment, species, year) %>%
  summarise(cover = sum(value)/2) %>% # just doing mean does not work
  ungroup() %>%
  pivot_wider(id_cols = c(block, treatment, year), 
              names_from = species, 
              values_fill = 0,
              values_from = cover)%>%
  mutate(row = str_c(block, treatment, year)) %>%
  dplyr::select(-block, -treatment, -year) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("row")
