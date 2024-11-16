# ordinations 
library(tidyverse)
library(vegan)
library(ggpubr)
library(ggrepel)
source("R/a_church_data_prep.R")


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
  mutate_if(is.character, as.factor) |>
  mutate(treatment = str_replace_all(treatment, "0", "ctl") |> str_to_upper())

# plant community ==============================================================
nms_a <- metaMDS(comm_both |> decostand("pa"), trymax = 100)

stressplot(nms_a)

site_scores <- as.data.frame(vegan::scores(nms_a)$sites) %>%
  as_tibble(rownames = "plot") %>%
  mutate(year = str_extract(plot, "\\d{4}"),
         plot = str_remove_all(plot,"\\d{4}"),
         treatment = str_remove_all(plot, "b\\d{1}") |> str_remove_all("_"),
         treatment = str_replace_all(treatment, "c", "CTL") |> str_to_upper()) |>
  mutate(treatment = fct_relevel(treatment, 'CTL', 'B', "BM", "M")) 
           
           

ef <- envfit(nms_a, comm_both, na.rm = T, permutations = 9999)

sp <-as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))

species <- as.data.frame(cbind(sp, p=ef$vectors$pvals)) %>%
  tibble::rownames_to_column("species")  %>%
  filter(p < 0.005)

adoc <- adonis2(comm_both |> decostand("pa") ~ treatment + year, data = site_scores)

p_oc <- ggplot(site_scores, aes(x=NMDS1, y=NMDS2)) +
  # coord_fixed() +
  geom_point(size=2, 
             aes(color = treatment, shape = year)) +
  # geom_segment(data = species,x=0,y=0, color = "grey",arrow = arrow(),
  #              aes(yend = NMDS2, xend = NMDS1), lwd=1) +
  # ggrepel::geom_text_repel(data = species,size=4, aes(label = species), color = "grey40") +
  theme_classic() +
  stat_ellipse(aes(group = year)) +
  scale_shape_manual(values = c(17,19))+
  scale_color_brewer(palette = "Set1") +
  theme(panel.background = element_rect(fill="transparent", color = "black"),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        # legend.background = element_rect(fill = 'transparent', color = 'black')
        ) +
  ggtitle("a. Plant community");p_oc  

# microbial community ==========================================================
key <- readxl::read_xlsx("data/Church_Park_Microbial_Data.xlsx", sheet = 1) 
ra_ITS <- readxl::read_xlsx("data/Church_Park_Microbial_Data.xlsx", sheet = 2) |>
  mutate(taxonomy = str_replace_all(taxonomy, "; ", "_"))
ra_16s <- readxl::read_xlsx("data/Church_Park_Microbial_Data.xlsx", sheet = 4) |>
  mutate(taxonomy = str_replace_all(taxonomy, "; ", "_"))

commf <- ra_ITS %>%
  pivot_longer(cols = names(.)[2:ncol(.)]) %>%
  pivot_wider(id_cols = name, names_from = taxonomy, values_from = value, values_fn = sum)|>
  filter(str_sub(name,1,3) != "C_K") |>
  mutate(name = str_replace_all(name,"B_M", "BM"))|>
  filter(!str_detect(name, '_15')) |>
  tibble::column_to_rownames("name")

commb <- ra_16s %>%
  pivot_longer(cols = names(.)[2:ncol(.)]) %>%
  pivot_wider(id_cols = name, names_from = taxonomy, values_from = value, values_fn = sum)|>
  filter(str_sub(name,1,3) != "C_K") |>
  mutate(name = str_replace_all(name,"B_M", "BM"))|>
  filter(!str_detect(name, '_15')) |>
  tibble::column_to_rownames("name")

nmds <- vegan::metaMDS(commf, trymax=100)
nmdsb <- vegan::metaMDS(commb, trymax = 100)

envfit(nmds, d |> dplyr::select(nitrifiers, EMF), 9999)
envfit(nmdsb, d |> dplyr::select(nitrifiers, EMF), 9999)


paf <- nmds$points |>
  as.data.frame() |>
  tibble::rownames_to_column("id") |>
  tidyr::separate(id,into = c("project", "block", "treatment", "plot", "depth"))  |>
  mutate(treatment = fct_relevel(treatment, 'CTL', 'B', "BM", "M")) |>
  ggplot(aes(x=MDS1, y=MDS2)) +
  geom_point(size=2, aes(color = treatment)) +
  stat_ellipse(aes(color = treatment)) +
  ggtitle("c. Fungal community (ITS)") +
  theme_classic() +
  xlab("NMDS1") +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(color="black", fill=NA));paf

pof <- nmdsb$points |>
  as.data.frame() |>
  tibble::rownames_to_column("id") |>
  tidyr::separate(id,into = c("project", "block", "treatment", "plot", "depth")) |>
  mutate(treatment = fct_relevel(treatment, 'CTL', 'B', "BM", "M")) |>
  ggplot(aes(x=MDS1, y=MDS2)) +
  geom_point(size=2, aes(color = treatment)) +
  stat_ellipse(aes(color = treatment)) +
  scale_color_brewer(palette = "Set1") +
  theme_classic() +
  xlab("NMDS1") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(color="black", fill=NA)) +
  ggtitle("b. Bacterial community (16s)");pof

top <- ggarrange(p_oc, pof, paf, nrow = 1, ncol =3, common.legend = TRUE, legend = 'bottom')
ggsave('out/top_panel.png', width = 7.5, height = 4, bg='white')

# bottom panel: soil ===========================================================




glimpse(d)

soil_pca <- d |>
  mutate(DON = TDN - DIN,
         `C:N` = total_c/total_n) |>
  dplyr::select(K = potassium, Ca = calcium, PO4 = phosphate, SO4 = sulfate, Na = na_mg_g,
                Cations = cations, DON, `C:N`,# DIN, 
                TC = total_c,
                TN = total_n, DOC, TDN, NH4 = ammonium, NO3 = nitrate) |>
  datawizard::standardise() |>
  prcomp()


rotation = as_tibble(soil_pca$rotation, rownames = 'row') |>
  mutate_if(is.numeric, function(x){x*10})


spca <- soil_pca$x |>
  as.data.frame() |>
  cbind(d$treatment) |>
  janitor::clean_names() |>
  mutate(d_treatment = fct_relevel(d_treatment, 'CTL', 'B', "BM", "M")) |>
  ggplot(aes(x=pc1, y=pc2)) +
  geom_point(aes(color = d_treatment), size=2) +
  stat_ellipse(aes(color = d_treatment)) +
  geom_text_repel(data = rotation, aes(x=PC1, y=PC2, label  = row), fontface = 'bold') +
  ggtitle("d. PCA on Soil Variables") +
  xlab('Component 1') +
  ylab('Component 2') +
  scale_color_brewer(palette = "Set1") +
  coord_equal() +
  theme_classic() +
  theme(panel.background = element_rect(color = 'black'),
        legend.title = element_blank(),
        legend.position = 'none',
        legend.justification = c(1,1),
        legend.background = element_rect(fill = NA));spca

soil_boxes <- d |>
  mutate(CN = total_c/total_n) |>
  dplyr::select(pH, TN = total_n, VWC = vwc, TC = total_c, TDN, Cations = cations, CN, DOC, TDN, NH4 = ammonium, NO3 = nitrate, treatment, EMF, Nitrifiers = nitrifiers) |>
  pivot_longer(cols = c(pH, TN, TC, TDN, Cations, EMF, CN, DOC, TDN, NH4, NO3, Nitrifiers, VWC)) |>
  mutate(treatment = fct_relevel(treatment, 'CTL', 'B', "BM", "M")) |>
  ggplot(aes(x=treatment, y = value, fill = treatment)) +
    scale_fill_brewer(palette = "Set1") +
    geom_boxplot(outliers= F) +
    facet_wrap(~name, scales = 'free_y') +
  theme_bw() +
  ggtitle('e. Soil variables by treatment')+
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none'); soil_boxes

bottom <- ggarrange(spca, soil_boxes, nrow = 1, ncol =2);bottom

ggarrange(top, bottom, nrow =2, ncol =1) -> full

ggsave(plot = full, bg = 'white', filename = "out/figure_2_multipanel.png", width =8, height = 8)


