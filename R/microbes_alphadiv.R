# microbial diversity
source("R/a_church_data_prep.R")
fungi <- readxl::read_xlsx("data/Church_Park_Microbial_Data.xlsx", sheet =2)%>%
  pivot_longer(cols = names(.)[2:ncol(.)]) %>%
  pivot_wider(id_cols = name, names_from = taxonomy, values_from = value, values_fn = sum)|>
  filter(str_sub(name,1,3) != "C_K", !str_detect(name,".2")) |>
  mutate(name = str_replace_all(name,"B_M", "BM"))|>
  tibble::column_to_rownames("name")

fundiv <- data.frame(
shannon = fungi |> vegan::diversity(index = "shannon"),
simpson = fungi |> vegan::diversity(index = "simpson"),
invsimpson = fungi |> vegan::diversity(index = "invsimpson"),
nspp = fungi |> vegan::specnumber()
) |>
  mutate(pielou = shannon/max(shannon))

bacteria <- readxl::read_xlsx("data/Church_Park_Microbial_Data.xlsx", sheet =4)%>%
  pivot_longer(cols = names(.)[2:ncol(.)]) %>%
  pivot_wider(id_cols = name, names_from = taxonomy, values_from = value, values_fn = sum)|>
  filter(str_sub(name,1,3) != "C_K", !str_detect(name,".2")
         #,!str_detect(name, "_15")
         ) |>
  mutate(name = str_replace_all(name,"B_M", "BM"))|>
  tibble::column_to_rownames("name")

bacdiv <- data.frame(
  shannon = bacteria |> vegan::diversity(index = "shannon"),
  simpson = bacteria |> vegan::diversity(index = "simpson"),
  invsimpson = bacteria |> vegan::diversity(index = "invsimpson"),
  nspp = bacteria |> vegan::specnumber()
) |>
  mutate(pielou = shannon/max(shannon))


library(vegan)

# bacteria
sitez <- rownames(bacteria) |> str_to_lower() |> as_tibble() |> 
  tidyr::separate(value, into = c("study", "block", "treatment", "plot", "depth")) |>
  mutate(treatment = ifelse(treatment == "ctl", "c", treatment)) |>
  left_join(soil_23) |>
  left_join(gc_by_plot_23)
bacteriah <- decostand(bacteria, "hellinger")
vegan::rda(bacteriah ~ treatment + depth, data = sitez, scale = T) -> rda
pct_exp <- round(summary(rda)$constr.chi / summary(rda)$tot.chi *100, 2)
RsquareAdj(rda)
anova.cca(rda, step = 1000)
anova.cca(rda, step = 1000, by = "term")
anova.cca(rda, step = 1000, by = "axis")

sitez <- sitez |>
  cbind(scores(rda)$sites)
ggplot(sitez, aes(x=RDA1, y=RDA2,color = treatment)) +
  geom_point(aes()) +
  stat_ellipse() +
  facet_wrap(~depth)

# fungi
sitez <- rownames(fungi) |> str_to_lower() |> as_tibble() |> 
  tidyr::separate(value, into = c("study", "block", "treatment", "plot", "depth")) |>
  mutate(treatment = ifelse(treatment == "ctl", "c", treatment)) |>
  left_join(soil_23)|>
  left_join(gc_by_plot_23)
fungih <- decostand(fungi, "hellinger")

vegan::rda(fungih ~ treatment + depth, data = sitez, scale = T) -> rda
pct_exp <- round(summary(rda)$constr.chi / summary(rda)$tot.chi *100, 2)
RsquareAdj(rda)
anova.cca(rda, step = 1000)
anova.cca(rda, step = 1000, by = "term")
anova.cca(rda, step = 1000, by = "axis")

sitez <- sitez |>
  cbind(scores(rda)$sites)
ggplot(sitez, aes(x=RDA1, y=RDA2,color = treatment)) +
  geom_point(aes()) +
  stat_ellipse() +
  facet_wrap(~depth)

vegan::rda(fungih, data = sitez, scale =T) -> rda
sitez <- sitez |>
  cbind(scores(rda)$sites)
ggplot(sitez, aes(x=PC1, y=PC2,color = treatment)) +
  geom_point(aes()) +
  stat_ellipse() +
  facet_wrap(~depth)


