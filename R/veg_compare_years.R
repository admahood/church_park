source("R/a_church_data_prep.R")

# prevalence by species=========================================================
# idea, one of those flow plots
prevalences <- bind_rows(
  comm_long |>
    group_by(species) |>
    reframe(prevalence_plot = length(unique(str_c(block,treatment)))) |>
    arrange(desc(prevalence_plot)) |> mutate(year = "2023")
  ,
  comm_16_long |> 
    filter(value>0) |>
    group_by(species) |>
    reframe(prevalence_plot = length(unique(str_c(block,treatment)))) |>
    arrange(desc(prevalence_plot))  |> mutate(year = "2016")
)

prevalences |>
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
ggplot(aes(x=year, y=species, fill=prevalence_plot)) +
  geom_tile(color = "black") +
  scale_fill_viridis_c()
ggsave("out/species_x_year.png")


# abundances====================================================================
ab <- comm_long |>
  mutate(year = "2023") |>
  group_by(treatment, species, year) |>
  summarise(cover = mean(value)) |>
  bind_rows(comm_16_long |>
              mutate(year = "2016") |>
              group_by(treatment, species, year) |>
              summarise(cover = mean(value)))

ggplot(ab, aes(x=year, y=cover, fill=species)) +
  geom_bar(stat="identity", color = "black") +
  facet_wrap(~treatment)

longs <- comm_long |>
  mutate(year = "2023") |>
  bind_rows(comm_16_long |>
              mutate(year = "2016")) |>
  group_by(year, block, treatment) |>
  reframe(cover = sum(value)/2)# |>
  # group_by(year, treatment) |>
  # reframe(cover=mean(cover))

ggplot(longs, aes(x=treatment, y=cover, fill=year)) +
  geom_boxplot() +
  # geom_bar(stat="identity", position = "dodge", color = "black") +
  ggtitle("total veg cover")
ggsave(filename = "out/total_veg_cover_comparison.png")

# nmds =========================================================================
library(vegan)

nms_a <- metaMDS(comm_both, trymax = 100)

stressplot(nms_a)

site_scores <- as.data.frame(vegan::scores(nms_a)$sites) %>%
  as_tibble(rownames = "plot") %>%
  mutate(year = str_extract(plot, "\\d{4}"),
         plot = str_remove_all(plot,"\\d{4}"),
         treatment = str_remove_all(plot, "b\\d{1}")) 

ef <- envfit(nms_a, comm_both, na.rm = T, permutations = 9999)

sp <-as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))

species <- as.data.frame(cbind(sp, p=ef$vectors$pvals)) %>%
  tibble::rownames_to_column("species")  %>%
  filter(p < 0.005)

adoc <- adonis2(comm_both ~ treatment + year, data = site_scores)

p_ab <- ggplot(site_scores, aes(x=NMDS1, y=NMDS2)) +
  coord_fixed() +
  geom_point(size=4, aes(color = treatment, shape = year)) +
  geom_segment(data = species,x=0,y=0, color = "grey",arrow = arrow(),
               aes(yend = NMDS2, xend = NMDS1), lwd=1)+
  ggrepel::geom_text_repel(data = species,size=4, aes(label = species), color = "grey40") +
  theme_classic() +
  stat_ellipse(aes(group = year)) +
  theme(panel.background = element_rect(fill="transparent", color = "black")) +
  ggtitle("Abundance-Based", paste("Treatment explains", 
                                    signif(adoc$R2[1], 3) * 100, 
                                    "% of the variation, year explains", 
                                   signif(adoc$R2[2], 3) * 100,"%"));p_ab  
ggsave("out/nmds_veg_year_comparison.png")

# occurrence 
nms_a <- metaMDS(comm_both |> decostand("pa"), trymax = 100)

stressplot(nms_a)

site_scores <- as.data.frame(vegan::scores(nms_a)$sites) %>%
  as_tibble(rownames = "plot") %>%
  mutate(year = str_extract(plot, "\\d{4}"),
         plot = str_remove_all(plot,"\\d{4}"),
         treatment = str_remove_all(plot, "b\\d{1}")) 

ef <- envfit(nms_a, comm_both, na.rm = T, permutations = 9999)

sp <-as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))

species <- as.data.frame(cbind(sp, p=ef$vectors$pvals)) %>%
  tibble::rownames_to_column("species")  %>%
  filter(p < 0.005)

adoc <- adonis2(comm_both |> decostand("pa") ~ treatment + year, data = site_scores)

p_oc <- ggplot(site_scores, aes(x=NMDS1, y=NMDS2)) +
  coord_fixed() +
  geom_point(size=4, aes(color = treatment, shape = year)) +
  geom_segment(data = species,x=0,y=0, color = "grey",arrow = arrow(),
               aes(yend = NMDS2, xend = NMDS1), lwd=1) +
  ggrepel::geom_text_repel(data = species,size=4, aes(label = species), color = "grey40") +
  theme_classic() +
  stat_ellipse(aes(group = year)) +
  theme(panel.background = element_rect(fill="transparent", color = "black")) +
  ggtitle("Occurrence-Based", paste("Treatment explains", 
                                   signif(adoc$R2[1], 3) * 100, 
                                   "% of the variation, year explains", 
                                   signif(adoc$R2[2], 3) * 100,"%"));p_oc  
ggsave("out/nmds_veg_year_comparison_oc.png")
