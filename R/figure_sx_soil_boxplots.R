# setup
library(tidyverse)
library(cowplot)

# look up tables ===============================================================

lut_varcats <- c("CN" = "Soil Characteristics",
                 "TDN" = "Soil Characteristics",
                 "NO3" = "Soil Characteristics",
                 "DON" = "Soil Characteristics",
                 "DOC" = "Soil Characteristics",
                 "DIN" = "Soil Characteristics",
                 "NH4" = "Soil Characteristics",
                 "Mulch" = "Cover",
                 "Biochar" = "Cover",
                 "Vegetation" = "Cover",
                 "Bare_Ground" = "Cover",
                 "Fungal_Div" = "Mic.",
                 "BA_Div" = "Mic.",
                 "Nitrifiers" = "Mic.",
                 "EMF" = "Mic.",
                 "VWC" = "Soil Characteristics",
                 "TWI" = "Moisture",
                 "Sand" = "Moisture",
                 "Clay" = "Moisture",
                 "Aspect" = "Moisture",
                 "K" = "Other Chem",
                 "PO4" = "Soil Characteristics",
                 "pH" = "Soil Characteristics",
                 "Cations" = "Soil Characteristics")
# figure sx soil boxplots


ddd <-  read_csv("data/cp_trait_soil_data.csv") |>
  dplyr::select(
    treatment, VWC = vwc, 
    BA_Div = simpson_bacteria_5, EMF = EMF_0_5,
    Fungal_Div = simpson_fungi_5, pH,
    Bare_Ground = bare, Mulch = mulch, Biochar = biochar, Nitrifiers = nitrifiers,
    # TN = total_n_0_5, TC = total_c_0_5, 
    DOC = DOC_0_5, TDN = TDN_0_5, 
    NH4 = ammonium_0_5, NO3 = nitrate_0_5, PO4 = phosphate_0_5, Cations = cations_0_5,
    # K = potassium_0_5, 
    DIN = DIN_0_5, DON = DON_0_5, Vegetation = total_veg_cover
  ) |>
  mutate(treatment = str_to_upper(treatment) |> str_replace_all("0", "CTL"),
         treatment = fct_relevel(treatment, 'CTL', 'B', "BM", "M")) |>
  unique() %>%
  pivot_longer(cols = names(.)[2:ncol(.)]) |>
  mutate(cat = lut_varcats[name]) 

row1 <- ddd |>
  filter(cat == 'Soil Characteristics') |>
  ggplot() +
  geom_boxplot(aes(x=treatment, y=value, fill = treatment), outliers = F) +
  facet_wrap(~name, scales = 'free_y', nrow = 2) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  ggtitle('Soil Characteristics')+
  theme(axis.title = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none'); row1

row2 <- ddd |>
  filter(cat == 'Cover') |>
  mutate(name = str_replace_all(name, "_", ' ')) |>
  ggplot() +
  geom_boxplot(aes(x=treatment, y=value, fill = treatment), outliers = F) +
  facet_wrap(~name, scales = 'free_y', nrow = 1) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  ggtitle('Ground Cover')+
  theme(axis.title = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none'); row2

row3 <- ddd |>
  filter(cat == 'Mic.') |>
  mutate(name = str_replace_all(name, "_", ' ')) |>
  ggplot() +
  geom_boxplot(aes(x=treatment, y=value, fill = treatment), outliers = F) +
  facet_wrap(~name, scales = 'free_y', nrow = 1) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  ggtitle('Microbes')+
  theme(axis.title = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none'); row3


full <- cowplot::ggdraw(xlim =c(0,5), ylim =c(0,3.8)) +
  draw_plot(row1, 0, 2, 5, 1.8) +  
  draw_plot(row2, 0, 1, 4, 1) +  
  draw_plot(row3, 0, 0, 4, 1) 
ggsave(plot = full, bg = 'white', filename = "out/figure_S3_soil_boxes.png", width =10, height =9)
