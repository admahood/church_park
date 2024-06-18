# data prep

library(tidyverse)
library(janitor)
library(lmerTest)

# church park ingest SLA, Height=============
if(!file.exists("data/church_traits.csv")){
  ijc <- list.files("data/IJ_church_park/", full.names=T) %>%
  lapply(FUN = read_csv)

  lut_trt <- c("1A" = "m","1B" = "m","1C" = "bm","1D" = "b",
               "2A" = "bm","2B" = "b","2C" = "m","2D" = "0",
               "3A" = "m","3B" = "0","3C" = "bm","3D" = "b",
               "4A" = "b","4B" = "bm","4C" = "0","4D" = "m",
               "5A" = "b","5B" = "m","5C" = "0","5D" = "bm",
               "6A" = "0","6B" = "b","6C" = "m","6D" = "bm")
  
  la_church <- ijc %>%
    bind_rows() %>%
    janitor::clean_names() %>%
    tidyr::separate(label, into = c("plot", "species", "individual","ext", "end1", "end2")) %>%
    dplyr::group_by(plot, species, individual) %>%
    summarise(leaf_area_scanner = sum(area*percent_area/100)) %>%
    ungroup() %>%
    mutate(species_full = case_match(species, 
                                     "CASP" ~ "Carex sp")) %>%
    dplyr::select(-species)
  
  d_church_carex_sla <- readxl::read_xlsx("data/church_veg_trait_data.xlsx") |>
    filter(trait == "SLA") %>%
    mutate(species_full = case_match(species_code, 
                                     "CASP" ~ "Carex sp",
                                     "VASP" ~ "Vaccinium sp",
                                     "SOSP" ~ "Oreochrysum parryi"), 
           individual = as.character(individual)) %>%
    dplyr::filter(species_full == "Carex sp") %>%
    left_join(la_church, by = c("plot", "species_full", "individual")) %>%
    dplyr::select(-date_yyyymmdd, -observed_by, -location, -leaf_area, -val_area,
                  -real_val_area, -species_code, -trait, -wet_weight_g, -length_cm) %>%
    mutate(lma = dry_weight_g/leaf_area_scanner,
           sla = (leaf_area_scanner/dry_weight_g)/1000) |>
    dplyr::select(plot, individual, lma, sla)
  
  d_church_carex <- readxl::read_xlsx("data/church_veg_trait_data.xlsx") |>
    filter(trait == "LDMC") %>%
    mutate(species_full = case_match(species_code, 
                                     "CASP" ~ "Carex sp",
                                     "VASP" ~ "Vaccinium sp",
                                     "SOSP" ~ "Oreochrysum parryi"), 
           individual = as.character(individual)) %>%
    dplyr::filter(species_full == "Carex sp") %>%
    dplyr::select(-date_yyyymmdd, -observed_by, -location, -leaf_area, -val_area,
                  -real_val_area, -species_code, -trait, -entered_by) %>%
    dplyr::rename(height_cm = length_cm) |>
    mutate(height_cm = as.numeric(height_cm),
           ldmc = (dry_weight_g*1000)/as.numeric(wet_weight_g)) %>%
    dplyr::select(-wet_weight_g, -dry_weight_g) |>
    left_join(d_church_carex_sla)
  
  d_church <- readxl::read_xlsx("data/church_veg_trait_data.xlsx",
                                    col_types = c("date", "text", "text", 
                                                  "text", "text", "text",
                                                  "text", "numeric", "numeric",
                                                  "numeric", "numeric", "numeric", "numeric",
                                                  "text")) |> 
    filter(trait == "all") %>%
    mutate(species_full = case_match(species_code,
                                     "VASP" ~ "Vaccinium sp",
                                     "SOSP" ~ "Oreochrysum parryi"), 
           individual = as.character(individual)) %>%
    dplyr::select(-date_yyyymmdd, -observed_by, -location, -species_code, -trait, -entered_by) %>%
    mutate(adjusted_la = leaf_area * lm(real_val_area ~ val_area, data = .)$coefficients[2] +
             rnorm(sd = 0.0007076 * sqrt(nrow(.)), n = nrow(.))) %>% # drawn from the prior lm
    dplyr::rename(height_cm = length_cm) %>%
    mutate(lma = dry_weight_g/adjusted_la,
           height_cm = height_cm |> as.numeric(),
           sla = (adjusted_la/dry_weight_g)/1000, # something's off with the units... doesn't actually matter
           ldmc = (dry_weight_g*1000)/wet_weight_g) %>%
    dplyr::select(plot, individual, height_cm, sla, ldmc, species_full, lma) %>%
    bind_rows(d_church_carex) %>%
    mutate(treatment = lut_trt[plot],
           block = str_c("b", str_sub(plot, 1,1)))
  
  write_csv(d_church, "data/church_traits.csv")}else{
  d_church <- read_csv("data/church_traits.csv")
}

# church park visualization ====================================================
# 
# p1 <- ggplot(d_church, aes(x = sla , fill = species_full)) +
#   geom_density(alpha=0.5) +
#   geom_rug(aes(color = species_full)) +
#   guides(color = "none")
# 
# p2 <- ggplot(d_church, aes(x = ldmc , fill = species_full)) +
#   geom_density(alpha=0.5) +
#   geom_rug(aes(color = species_full)) +
#   guides(color = "none")
# 
# p3 <- ggplot(d_church, aes(x = height_cm , fill = species_full)) +
#   geom_density(alpha=0.5) +
#   geom_rug(aes(color = species_full)) +
#   guides(color = "none")
# 
# ggpubr::ggarrange(p1, p2, p3, nrow = 1, ncol=3, common.legend = TRUE) %>%
#   ggsave(plot = ., filename = "out/density_plots.png", width =9, height =5)
# 
# ggplot(d_church, aes(x = height_cm, y=ldmc, color = species_full)) +
#   geom_point()
# ggplot(d_church, aes(x = height_cm, y=sla, color = species_full)) +
#   geom_point()
# ggplot(d_church, aes(x = sla, y=ldmc, color = species_full)) +
#   geom_point() +
#   ggtitle("LDMC tends to be inversely related to SLA and Leaf thickness (P-H 2013).",
#           "species with low LDMC tend to be associated with productive, often highly disturbed environments.")
# 
# p1 <- ggplot(d_church, aes(x = trt, y=sla, fill = trt)) +
#   geom_boxplot() +
#   facet_wrap(~species_full, scales="free") +
#   theme(legend.position = "none");p1
# p2 <- ggplot(d_church, aes(x = trt, y=ldmc, fill = trt)) +
#   geom_boxplot() +
#   facet_wrap(~species_full, scales="free") +
#   theme(legend.position = "none")
# p3 <- ggplot(d_church, aes(x = trt, y=height_cm, fill = trt)) +
#   geom_boxplot() +
#   facet_wrap(~species_full, scales="free") +
#   theme(legend.position = "none")
# 
# ggpubr::ggarrange(p1, p2, p3, nrow = 3, ncol=1, common.legend = TRUE) %>%
#   ggsave(plot = ., filename = "out/trait_x_treatment.png", width =7, height =10)
