# soil pca

source("R/a_church_data_prep.R")
source("R/a_itv_data_prep.R")
library(ggtext)
library(ggrepel)
# data_prep

d <- sites_w_23_soil |>
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


glimpse(d)

soil_pca <- d |>
  dplyr::select(potassium, calcium, phosphate, sulfate, cations, DIN,total_c, 
                total_n, DOC, TDN, ammonium, nitrate) |>
  datawizard::standardise() |>
  prcomp()


rotation = as_tibble(soil_pca$rotation, rownames = 'row') |>
  mutate_if(is.numeric, function(x){x*10})


soil_pca$x |>
  as.data.frame() |>
  cbind(d$treatment) |>
  janitor::clean_names() |>
  ggplot(aes(x=pc1, y=pc2)) +
  geom_point(aes(color = d_treatment)) +
  stat_ellipse(aes(color = d_treatment)) +
  geom_text_repel(data = rotation, aes(x=PC1, y=PC2, label  = row)) +
  ggtitle("PCA on Soil Variables") +
  xlab('Component 1') +
  ylab('Component 2') +
  coord_equal() +
  theme_classic() +
  theme(panel.background = element_rect(color = 'black'),
        legend.title = element_blank(),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_rect(fill = NA))

ggsave('out/soil_pca_tf.png', width = 5, height = 5.5, bg = 'white')
