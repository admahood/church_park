# microbial diversity
library(tidyverse)
fungi <- readxl::read_xlsx("data/Church_Park_Microbial_Data.xlsx", sheet =2)%>%
  pivot_longer(cols = names(.)[2:ncol(.)]) %>%
  pivot_wider(id_cols = name, names_from = taxonomy, values_from = value, values_fn = sum)|>
  filter(str_sub(name,1,3) != "C_K") |>
  mutate(name = str_replace_all(name,"B_M", "BM"))|>
  tibble::column_to_rownames("name")

data.frame(
shannon = fungi |> vegan::diversity(index = "shannon"),
simpson = fungi |> vegan::diversity(index = "simpson"),
invsimpson = fungi |> vegan::diversity(index = "invsimpson"),
nspp = fungi |> vegan::specnumber()
) |>
  mutate(pielou = shannon/max(shannon))

