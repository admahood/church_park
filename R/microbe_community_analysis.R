# microbe data
library(tidyverse)
library(vegan)
library(topomicro)
library(ggrepel)
library(janitor)
key <- readxl::read_xlsx("data/Church_Park_Microbial_Data.xlsx", sheet = 1) 
ra_ITS <- readxl::read_xlsx("data/Church_Park_Microbial_Data.xlsx", sheet = 2) |>
  mutate(taxonomy = str_replace_all(taxonomy, "; ", "_"))
ra_16s <- readxl::read_xlsx("data/Church_Park_Microbial_Data.xlsx", sheet = 4) |>
  mutate(taxonomy = str_replace_all(taxonomy, "; ", "_"))

soil <- read_csv('data/cp_soil_nutrients_almModified.xlsx - Sheet1.csv') |>
  janitor::clean_names() |>
  dplyr::select(-ends_with('mg_l'), -sample_num, -plot_id) |>
  mutate(sample_id = str_replace_all(sample_id, "B_M", "BM"),
         depth = ifelse(depth == "0-5", "5", "15"),
         trt = as.factor(trt),
         blk = as.factor(blk),
         depth = as.factor(depth)) |>
  tibble::column_to_rownames("sample_id"); glimpse(soil)

# ITS ==========================================================================
commf <- ra_ITS %>%
  pivot_longer(cols = names(.)[2:ncol(.)]) %>%
  pivot_wider(id_cols = name, names_from = taxonomy, values_from = value, values_fn = sum)|>
  filter(str_sub(name,1,3) != "C_K") |>
  mutate(name = str_replace_all(name,"B_M", "BM"))|>
  filter(!str_detect(name, '_15')) |>
  tibble::column_to_rownames("name")

nmds <- vegan::metaMDS(commf, trymax=100)
nmdspa <- vegan::metaMDS(decostand(commf, 'pa'), trymax = 100)

fitf<-envfit(nmds, 
            soil |> dplyr::select(-blk) |> filter(depth != "15"), 
            permutations = 9999, 
            strata=soil |> filter(depth != "15") |> pull(blk)) |>
  topomicro::tidy_envfit()|>
  filter(p < 0.05) |>
  mutate(NMDS1 = (NMDS1)* .75,
         NMDS2 = (NMDS2)* .75);fitf
fitpaf<-envfit(nmdspa, 
              soil |> dplyr::select(-blk) |> filter(depth != "15"), 
              permutations = 9999, 
              strata=soil |> filter(depth != "15") |> pull(blk)) |>
  topomicro::tidy_envfit()|>
  filter(p < 0.05) |>
  mutate(NMDS1 = (NMDS1)* .5,
         NMDS2 = (NMDS2)* .5);fitpaf



paf <- nmds$points |>
  as.data.frame() |>
  tibble::rownames_to_column("id") |>
  tidyr::separate(id,into = c("project", "block", "treatment", "plot", "depth")) |>
  ggplot(aes(x=MDS1, y=MDS2)) +
  geom_point(size=2, aes(color = treatment)) +
  facet_wrap(~depth, scales="free") +
  stat_ellipse(aes(color = treatment)) +
  geom_segment(data = fitf, aes(x=0, y=0, xend = NMDS1, yend=NMDS2), arrow = arrow()) +
  geom_text_repel(data = fitf, aes(x=NMDS1, y=NMDS2, label = var)) +
  ggtitle("Abundance, ITS (Fungi)") +
  theme_classic() +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(color="black", fill=NA));paf

pof <- nmdspa$points |>
  as.data.frame() |>
  tibble::rownames_to_column("id") |>
  tidyr::separate(id,into = c("project", "block", "treatment", "plot", "depth")) |>
  ggplot(aes(x=MDS1, y=MDS2)) +
  geom_point(size=2, aes(color = treatment)) +
  facet_wrap(~depth, scales="free") +
  stat_ellipse(aes(color = treatment)) +
  geom_segment(data = fitpaf, aes(x=0, y=0, xend = NMDS1, yend=NMDS2), arrow = arrow()) +
  geom_text_repel(data = fitpaf, aes(x=NMDS1, y=NMDS2, label = var), nudge_x = .1) +
  scale_color_brewer(palette = "Set1") +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(color="black", fill=NA)) +
  ggtitle("Occurrence, ITS (Fungi)");pof

adonis2(commf ~ treatment+ block + soil_p_h + na_mg_g + po4_mg_g , 
        data = sites_w_23_soil)

# ggpubr::ggarrange(paf, pof, nrow=2) |>
#   ggsave(filename = "out/nmds_microbes_fungi_d5.png", width=7, height=7)



# 16s ==========================================================================
commb <- ra_16s  %>%
  pivot_longer(cols = names(.)[2:ncol(.)]) %>%
  pivot_wider(id_cols = name, names_from = taxonomy, values_from = value, values_fn = sum)|>
  filter(str_sub(name,1,3) != "C_K") |>
  mutate(name = str_replace_all(name,"B_M", "BM"))|>
  filter(!str_detect(name, '_15')) |>
  tibble::column_to_rownames("name");commb

nmds <- vegan::metaMDS(commb, trymax=100)
nmdspa <- vegan::metaMDS(decostand(commb, 'pa'), trymax = 100)

fit<-envfit(nmds, 
            soil |> dplyr::select(-blk) |> filter(depth != "15"), 
            permutations = 9999, 
            strata=soil |> filter(depth != "15") |> pull(blk)) |>
  topomicro::tidy_envfit()|>
  filter(p < 0.05) |>
  mutate(NMDS1 = (NMDS1)* .75,
         NMDS2 = (NMDS2)* .75);fit
fitpa<-envfit(nmdspa, 
              soil |> dplyr::select(-blk) |> filter(depth != "15"), 
              permutations = 9999, 
              strata=soil |> filter(depth != "15") |> pull(blk)) |>
  topomicro::tidy_envfit()|>
  filter(p < 0.05) |>
  mutate(NMDS1 = (NMDS1)* .5,
         NMDS2 = (NMDS2)* .5);fitpa



pa <- nmds$points |>
  as.data.frame() |>
  tibble::rownames_to_column("id") |>
  tidyr::separate(id,into = c("project", "block", "treatment", "plot", "depth")) |>
  ggplot(aes(x=MDS1, y=MDS2)) +
  geom_point(size=2, aes(color = treatment)) +
  facet_wrap(~depth, scales="free") +
  stat_ellipse(aes(color = treatment)) +
  geom_segment(data = fit, aes(x=0, y=0, xend = NMDS1, yend=NMDS2), arrow = arrow()) +
  geom_text_repel(data = fit, aes(x=NMDS1, y=NMDS2, label = var)) +
  ggtitle("Abundance, 16S (Bacteria)") +
  theme_classic() +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(color="black", fill=NA));pa

po <- nmdspa$points |>
  as.data.frame() |>
  tibble::rownames_to_column("id") |>
  tidyr::separate(id,into = c("project", "block", "treatment", "plot", "depth")) |>
  ggplot(aes(x=MDS1, y=MDS2)) +
  geom_point(size=2, aes(color = treatment)) +
  facet_wrap(~depth) +
  stat_ellipse(aes(color = treatment))+
  geom_segment(data = fitpa, aes(x=0, y=0, xend = NMDS1, yend=NMDS2), arrow = arrow()) +
  geom_text_repel(data = fitpa, aes(x=NMDS1, y=NMDS2, label = var)) +
  scale_color_brewer(palette = "Set1") +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(color="black", fill=NA)) +
  ggtitle("Occurrence, 16S (Bacteria)");po

ggpubr::ggarrange(pa, po, paf, pof, ncol=2, nrow=2, common.legend = T, legend="bottom") |>
  ggsave(filename = "out/nmds_microbes_5in.png", width=7, height=7, bg="white")

adonis2(commb ~ treatment+ block + soil_p_h + na_mg_g + po4_mg_g, 
        data = sites_w_23_soil)
