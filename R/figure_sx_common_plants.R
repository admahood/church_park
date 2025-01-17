library(tidyverse)
source('R/a_church_data_prep.R')
common_plant_differences <- comm_both |>
  as_tibble(rownames = "row") |> 
  tidyr::separate(row, into = c('block', "treatment", 'year'), sep = "_") |>
  dplyr::select(block, year, treatment, Oreochrysum_parryi, Vaccinium_scoparium) |>
  mutate(treatment = str_remove_all(treatment,"_")) |>
  pivot_longer(-c(block, year, treatment)) |>
  pivot_wider(names_from = year, values_from = value, values_fn = mean) |>
  mutate(diff = `2023` - `2016`,
         treatment = str_replace_all(treatment, "c", "0"),
         name = str_replace_all(name, "_", " "))|>
  mutate(treatment = str_to_upper(treatment) |> str_replace_all("0", "CTL"),
         treatment = fct_relevel(treatment, 'CTL', 'B', "BM", "M"))


ts <- data.frame(treatment = NA, name = NA, p = NA, ypos = NA)
counter <- 1
for(i in unique(common_plant_differences$name)){
  for(t in unique(common_plant_differences$treatment)){
    p <- t.test(x = common_plant_differences |> filter(name == i, treatment == t) |> pull(diff))$p.value
    ts[counter,1] <- t
    ts[counter,2] <- i
    ts[counter,3] <- p
    ts[counter,4] <- common_plant_differences |> filter(name == i, treatment == t) |> pull(diff) |> max()
    
    counter <- counter + 1
    
  }
}
sigs <- filter(ts, p < 0.05)
ggplot(common_plant_differences |>
         mutate(treatment = str_to_upper(treatment) |> str_replace_all("0", "CTL"),
                treatment = fct_relevel(treatment, 'CTL', 'B', "BM", "M")) ,
       aes(y=diff, fill=treatment, x=name)) +
  geom_boxplot(outliers = F) +
  geom_hline(yintercept = 0) +
  geom_point(data = sigs, y=c(4, 9, 7), x = c(.91, .72, 1.28), shape=8,
             size=3, stroke = 1, show.legend = F) +
  ylab("Change in Percent Cover") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.title.x= element_blank()) 

ggsave("out/figure_s2_common_plants.png", width = 5, height = 3.5, bg="white")
