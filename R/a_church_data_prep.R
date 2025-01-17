# church park data prep
library(tidyverse)
library(janitor)
library(sf)
library(terra)
library(topomicro)
library(solrad)
# getting topography stuff

ext <- c(-106, -105.8, 39.9, 40)
dem <- terra::rast("data/dem/n39_w106_1arc_v3.bil") |>
  terra::crop(ext)
names(dem) <- "elevation"

terrain <- terra::terrain(dem, v = c("slope", "aspect", "TPI", "TRI", "flowdir"))
twi <- topomicro::get_twi(dem, resolution = 30)

topo <- read_csv("data/fraser - locations.csv") |>
  st_as_sf(crs = 4326, coords = c("longitude", "latitude"), remove=F) %>%
  mutate(terra::extract(terrain, .),
         terra::extract(twi$twi, .),
         terra::extract(dem, .),
         folded_aspect = topomicro::get_folded_aspect(aspect),
         block = str_c("b", plot)) |>
  st_set_geometry(NULL) |>
  dplyr::select(-ID, -plot)


soil_23 <- read_csv('data/cp_soil_nutrients_almModified.xlsx - Sheet1.csv') |>
  janitor::clean_names() |>
  dplyr::select(-ends_with('mg_l'), -sample_num, -plot_id) |>
  mutate(sample_id = str_replace_all(sample_id, "B_M", "BM"),
         depth = ifelse(depth == "0-5", "5", "15"),
         treatment = str_to_lower(trt) |> as.factor(),
         block = str_c("b", blk),
         depth = as.factor(depth)) |>
  dplyr::select(-blk, -trt)  ; glimpse(soil_23)
glimpse(soil_23)

soil_23_tf <- readxl::read_xlsx('data/ConsolidatedChurchPark_for_Adam_070824.xlsx') |>
  dplyr::select(-1) |>
  dplyr::rename(block = BLK, plot = PlotID, treatment = Trt, depth = Depth) |>
  janitor::clean_names() |>
  mutate_if(is.character, str_to_lower) |>
  dplyr::filter(depth == '0-5') |>
  dplyr::select(-depth) |>
  dplyr::mutate(block = str_c('b', block),
                plot = str_c(block, plot),
                treatment = str_replace_all(treatment, 'c', '0'))
glimpse(soil_23_tf)

soil_23_tf_both_depths <- readxl::read_xlsx('data/ConsolidatedChurchPark_for_Adam_070824.xlsx') |>
  dplyr::select(-1) |>
  dplyr::rename(block = BLK, plot = PlotID, treatment = Trt, depth = Depth) |>
  janitor::clean_names() |>
  mutate_if(is.character, str_to_lower) %>%
  tidyr::pivot_wider(id_cols = c('block', 'plot', 'treatment'), names_from = depth, values_from = names(.)[5:length(.)]) |>
  janitor::clean_names() |>
  dplyr::mutate(block = str_c('b', block),
                plot = str_c(block, plot),
                treatment = str_replace_all(treatment, 'c', '0'))
glimpse(soil_23_tf_both_depths)
# ================================
# soil_pca <- prcomp(soil_23[,3:20], scale = T, center = T)
# round((soil_pca$sdev^2 / sum(soil_pca$sdev^2))*100, 2)
# biplot(soil_pca)
# 
# rotation = as_tibble(soil_pca$rotation, rownames = 'row') 
# 
# soil_23$PC1 <- soil_pca$x[,1]
# soil_23$PC2 <- soil_pca$x[,2]
# 
# ggplot(soil_23, aes(color=depth, x=PC1, y=PC2)) +
#   geom_point() 
# library(ggrepel)
# ggplot(soil_23, aes(x=PC1, y=PC2)) +
#   geom_point(size=3, aes(color=treatment)) +
#   facet_wrap(~depth) +
#   stat_ellipse(aes(color = treatment)) +
#   geom_text_repel(data = rotation, aes(x=PC1*10, y=PC2*10, label = row))
# 
# ##
# 
# soil5 <- soil_23 |> filter(depth == "5") |> mutate(soil_p_h = soil_p_h * -1)
# 
# soil_pca <- prcomp(soil5[,3:20], scale = T, center = T)
# round((soil_pca$sdev^2 / sum(soil_pca$sdev^2))*100, 2)
# biplot(soil_pca)
# 
# rotation = as_tibble(soil_pca$rotation, rownames = 'row') 
# 
# soil5$PC1 <- soil_pca$x[,1]
# soil5$PC2 <- soil_pca$x[,2]
# 
# library(ggrepel)
# ggplot(soil5, aes(x=PC1, y=PC2)) +
#   geom_point(size=3, aes(color=treatment)) +
#   facet_wrap(~depth) +
#   stat_ellipse(aes(color = treatment)) +
#   geom_text_repel(data = rotation, aes(x=PC1*10, y=PC2*10, label = row)) +
#   theme_bw()
# ggsave("out/soil_pca_0-5.png", width = 7, height=7, bg='white')
# 
# ggplot(sites_w_23_soil, aes(x=tdn_mg_g, y=shannon_fungi_5)) +
#   geom_point()


# =============================
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

# funguild

funguild <- readxl::read_xlsx("data/simplified_FUNguild_ks_2024.xlsx") |>
  janitor::clean_names() |>
  pivot_longer(-guild) |>
  mutate(name = str_replace_all(name, "b_m", "bm") |>
           str_replace_all("ctl", "0")) |>
  tidyr::separate(name, into = c("study", "block", "treatment","plot", "depth")) |>
  dplyr::select(-study, -plot) |>
  mutate(guild = str_to_lower(guild) |>
           str_replace_all(" ", "_")) |>
  pivot_wider(names_from = guild, values_from = value) |>
  dplyr::select(-sums)

funguild_wider <-  readxl::read_xlsx("data/simplified_FUNguild_ks_2024.xlsx") |>
  janitor::clean_names() |>
  pivot_longer(-guild) |>
  mutate(name = str_replace_all(name, "b_m", "bm") |>
           str_replace_all("ctl", "0")) |>
  tidyr::separate(name, into = c("study", "block", "treatment","plot", "depth")) |>
  dplyr::select(-study, -plot) |>
  mutate(guild = str_to_lower(guild) |>
           str_replace_all(" ", "_")) |>
  pivot_wider(names_from = c(guild, depth), values_from = value) |>
  dplyr::select(-starts_with('sums'))

# nitrifier data

# nitrifier_input <- readxl::read_xlsx("data/Nitrifier_Feature_table.xlsx") |>
#   janitor::clean_names() |>
#   pivot_longer(-c(full_taxonomy_string, family)) |>
#   mutate(name = str_replace_all(name, "b_m", "bm") |>
#            str_replace_all("ctl", "0")) |>
#   tidyr::separate(name, into = c("study", "block", "treatment","plot", "depth")) |>
#   dplyr::select(-study, -plot) %>%
#   mutate(dummy_name = 1:nrow(.) |> as.character(),
#          dummy_name = str_c("sp", dummy_name)) |>
#   filter(value > 0)
#   
# nitrifier_binary_matrix_5 <- 
#   nitrifier_input |>
#   filter(depth == 5) |>
#   dplyr::select(-depth) |>
#   arrange(block, treatment) |>
#   mutate(row = str_c(block, "_", treatment),
#          value = ifelse(value > 0, 1, 0)) |>
#   dplyr::select(-full_taxonomy_string, -family, -block, -treatment) |>
#   pivot_wider(names_from = dummy_name, values_from = value, values_fill = 0) |>
#   tibble::column_to_rownames("row")
# 
# nitrifier_traits <-
#   nitrifier_input |> dplyr::select(dummy_name, full_taxonomy_string, family)
# 
# 
# summary(nitrifier_binary_matrix)
# colSums(nitrifier_binary_matrix)
# rowSums(nitrifier_binary_matrix)

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
  mutate(treatment = treatment |> str_replace_all("c", "0")) |>
  left_join(topo) |>
  left_join(funguild_wider) |>
  arrange(block, treatment)

# write_csv(sites_w_23_soil, "data/sites_23_w_soil.csv")

# veg community 2016 =============================================
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

# veg comm 2023 =========================================
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

# both years joined ============================================================
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
  mutate(row = str_c(block, "_", treatment, "_", year)) %>%
  dplyr::select(-block, -treatment, -year) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("row")

# veg differences year to year

veg_compare <- comm_long %>%
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
  ungroup() |>
  pivot_wider(names_from = year, values_from = cover, values_fill = 0) |>
  mutate(cover_change = `2023` - `2016`,
         treatment = ifelse(treatment == 'c', '0', treatment))
# library(lmerTest)
# library(emmeans)
# mv <- lmer(cover_change ~ treatment + (1|block), data = veg_compare |> filter(species =="Vaccinium_scoparium")) 
# mc <- lmer(cover_change ~ treatment + (1|block), data = veg_compare |> filter(species =="Carex_sp")) 
# mo <- lmer(cover_change ~ treatment + (1|block), data = veg_compare |> filter(species =="Oreochrysum_parryi"))
# 
# lapply(list(mv, mo, mc), summary)
# 
# emmeans(mv, ~treatment) |> pairs()
# emmeans(mc, ~treatment) |> pairs()
# emmeans(mo, ~treatment) |> pairs()

unique(veg_compare$species) |> paste(collapse = "', '")

lut_fg <- c('Abies_lasiocarpa' = 'Tree', 
  'Achillea_millefolium' = "Forb", 
  'Arctostaphylus_uvaursi' = 'Shrub', 
  'Arnica_spp' = 'Forb',
  'Carex_sp' = 'Graminoid', 
  'Ceanothus_velutinus' = 'Shrub', 
  'Chamerion_angustifolium' = 'Forb', 
  'Cirsium_arvense' = 'Invasive', 
  'Elymus_spp' = 'Graminoid', 
  'Erigeron_speciosus' = 'Forb', 
  'Lactuca_serriola' = "Invasive",
  'Lathyrus_laetivirens' = 'Forb', 
  'Lupinus_sp' = 'Shrub', 
  'Mahonia_repens' = 'Shrub',
  'Oreochrysum_parryi' = 'Forb',
  'Other_Forbs' = 'Forb',
  'Other_Graminoids' = 'Graminoid', 
  'Other_Shrubs' = "Shrub",
  'Phacelia_heterophylla' = 'Forb', 
  'Phleum_pratense' = "Invasive", 
  'Pinus_contorta' = "Tree", 
  'Populus_tremuloides' = "Tree", 
  'Rosa_sp' = "Shrub", 
  'Senecio_spp' = 'Forb', 
  'Solidago_multiradiata' = 'Forb', 
  'Taraxacum_officinale' = "Invasive",
  'Tragopogon_dubius' = 'Invasive',
  'Vaccinium_scoparium' = 'Shrub',
  'Verbascum_thapsus' = 'Invasive',
  'Cirsium_sp' = "Invasive",
  'Shepherida_canadensis' = 'Shrub',
  'Lupinus_argenteus' = "Shrub")

oreo_vac <- veg_compare |>
  dplyr::filter(species %in% c("Oreochrysum_parryi", "Vaccinium_scoparium"))|>
  dplyr::mutate(species = ifelse(species == 'Oreochrysum_parryi', "Oreochrysum", "Vaccinium")) |>
  dplyr::rename(name = species) |>
  group_by(block, treatment, name) |>
  summarise_all(mean)

veg_diff <- veg_compare |>
  mutate(name = lut_fg[species]) |>
  dplyr::select(-species) |>
  dplyr::rename(diff = cover_change) |>
  group_by(block, treatment, name) |>
  summarise_all(mean)

# veg 2023 quadrat level =======================================================

comm_wide_q <- comm_long %>%
  group_by(block, treatment, quadrat, species) %>%
  summarise(cover = sum(value)/2) %>% # just doing mean does not work
  ungroup() %>%
  pivot_wider(id_cols = c(block, treatment, quadrat), 
              names_from = species, 
              values_fill = 0,
              values_from = cover); comm_wide_q

comm_q <- comm_wide_q %>%
  mutate(row = str_c(block, treatment, quadrat)) %>%
  dplyr::select(-block, -treatment, -quadrat) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("row"); comm_q



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
  ggplot(aes(y=value, fill = year)) +
  geom_boxplot() +
  ylab('Percent Cover') +
  facet_grid(treatment~name, scales = 'free') +
  theme_bw() +
  theme(legend.position ='none',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.justification = c(1,1),
        legend.background = element_rect(color = "black"))
pvd <- veg_diff |>
  filter(!name %in% c("root_stump", "cwd_coarse", "cwd_fine", "basal_live_plant", 'Tree')) |>
  pivot_longer(-c(block, treatment, name), names_to = 'year') |>
  filter(year != "diff") |>
  mutate(treatment = str_replace_all(treatment, "c", "0"),
         name = ifelse(name == "bare_soil_gravel", "Bare Ground", name)) |>
  ggplot(aes(y=value, fill = year)) +
  geom_boxplot() +
  ylab('Percent Cover') +
  facet_grid(treatment~name, scales = 'free') +
  theme_bw() +
  theme(legend.position = c(1,1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.justification = c(1,1),
        legend.background = element_rect(color = "black"))

ggarrange(pgcd, pvd, nrow =1, widths = c(5,4.5))

ggsave("out/ground_cover_differences.png", bg='white', width=9, height=6)

sort(names(gc16));sort(names(gc23))

ggplot(gc23, aes(x=treatment, y = value)) +
  geom_boxplot() +
  facet_wrap(~category, scales = "free")

ggplot(gc16, aes(x=treatment, y = value)) +
  geom_boxplot() +
  facet_wrap(~category, scales = "free")

gc_by_plot_23 <- gc23 |>
  group_by(block, treatment) |>
  summarise(litter_duff = mean(litter_duff),
            rock = mean(rock_1cm),
            bare = mean(bare_soil_gravel),
            wood_mulch = mean(wood_mulch),
            biochar = mean(biochar)) |>
  ungroup() |>
  mutate(treatment = ifelse(treatment == "c", "0", treatment))

gc_by_plot_16 <- gc16 |>
  group_by(block, treatment) |>
  summarise(litter_duff = mean(litter_duff),
            rock = mean(rock_1cm),
            bare = mean(bare_soil_gravel),
            wood_mulch = mean(wood_mulch),
            biochar = mean(biochar))|>
  ungroup() |>
  mutate(treatment = ifelse(treatment == "c", "0", treatment))

# fungal and bacterial communities

fungi <- readxl::read_xlsx("data/Church_Park_Microbial_Data.xlsx", sheet =2)%>%
  pivot_longer(cols = names(.)[2:ncol(.)]) %>%
  pivot_wider(id_cols = name, names_from = taxonomy, values_from = value, values_fn = sum)|>
  filter(str_sub(name,1,3) != "C_K", !str_detect(name,"\\.2"), !str_detect(name, "_15$")) |>
  mutate(name = str_replace_all(name,"B_M", "BM"))|>
  tibble::column_to_rownames("name") |>
  vegan::decostand("hellinger")

bacteria <- readxl::read_xlsx("data/Church_Park_Microbial_Data.xlsx", sheet =4)%>%
  pivot_longer(cols = names(.)[2:ncol(.)]) %>%
  pivot_wider(id_cols = name, names_from = taxonomy, values_from = value, values_fn = sum)|>
  filter(str_sub(name,1,3) != "C_K", !str_detect(name,"\\.2")
         ,!str_detect(name, "_15")
  ) |>
  mutate(name = str_replace_all(name,"B_M", "BM"))|>
  tibble::column_to_rownames("name") |>
  vegan::decostand("hellinger")

pcab <- vegan::rda(bacteria) |>
  vegan::scores() |>
  purrr::pluck(2) |>
  tibble::as_tibble(rownames = "name") |>
  dplyr::mutate(name = str_to_lower(name)) |>
  tidyr::separate(name, into = c("study", "block", "treatment", "plot", 'depth')) |>
  dplyr::select(pcab=PC1, block, treatment, PC2)|>
  ungroup() |>
  mutate(treatment = ifelse(treatment == "ctl", "0", treatment))
pcaf <- vegan::rda(fungi)|>
  vegan::scores() |>
  purrr::pluck(2) |>
  tibble::as_tibble(rownames = "name") |>
  dplyr::mutate(name = str_to_lower(name)) |>
  tidyr::separate(name, into = c("study", "block", "treatment", "plot", 'depth')) |>
  dplyr::select(pcaf = PC1, block, treatment)|>
  ungroup() |>
  mutate(treatment = ifelse(treatment == "ctl", "0", treatment))

nitrifiers <- readxl::read_xlsx('data/NItrifierFeatureTable.xlsx') |>
  filter(Full_Taxonomy_String == "summed relabund") |>
  dplyr::select(-c(1:2)) %>%
  pivot_longer(cols = names(.), values_to = "nitrifiers") |>
  mutate(name = str_replace_all(name, '_B_M_', "_BM_") |> str_to_lower() |> str_replace_all('ctl', '0')) |>
  tidyr::separate(name, into = c('project', 'block', 'treatment', 'plot', 'depth')) |>
  mutate(plot = str_c(block, plot)) |>
  filter(depth == '5') |>
  dplyr::select(-project, -depth)

discriminant_taxa <-
  readxl::read_xlsx('data/DiscrimnantTaxa_ChurchPark_071124.xlsx') |>
  tidyr::separate(taxonomy, into= c('d', 'phyllum', 'c', 'o', 'f', 'g', 's'), sep = "; ") |>
  dplyr::select(-c(d, c, o, f, g, s)) |>
  mutate(phyllum = str_remove_all(phyllum, 'p__')) %>%
  pivot_longer(cols = names(.)[2:ncol(.)]) |>
  mutate(name = str_replace_all(name, '_B_M_', "_BM_") |> str_to_lower() |> str_replace_all('ctl', '0')) |>
  tidyr::separate(name, into = c('project', 'block', 'treatment', 'plot', 'depth')) |>
  mutate(plot = str_c(block, plot)) |>
  filter(depth == '5') |>
  dplyr::select(-project, -depth) |>
  pivot_wider(names_from = phyllum, values_from = value)

# apply(nitrifiers, 2, unique)
# 
# 
# ggplot(pcab, aes(pcab, PC2,color = treatment)) +
#   geom_point(aes(), size=3) +
#   stat_ellipse()

# total veg cover =====================

tvc <- comm_long |>
  mutate(treatment = ifelse(treatment == 'c', '0', treatment)) |>
  group_by(block, treatment) |>
  summarise(total_veg_cover = sum(value)/2) |>
  ungroup()

