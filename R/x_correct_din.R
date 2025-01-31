# getting the properly calculated DIN
library(tidyverse)

d_old <- read_csv('data/cp_trait_soil_data.csv') |>
  dplyr::select(-starts_with('DIN'), -starts_with('DON')) # one column
d_tf <- readxl::read_xlsx('data/ConsolidatedChurchPark_for_Adam_070824_Updated.xlsx') |>
  dplyr::select(treatment = Trt, DIN, Depth, block = BLK) |>
  dplyr::mutate(block = str_c('b', block),
                treatment = str_to_lower(treatment) |> str_replace_all("c", '0'),
                Depth = str_replace_all(Depth,"-", "_")) |>
  pivot_wider(names_from = Depth, values_from = DIN, names_glue = "DIN_{Depth}")
d_new <- left_join(d_old, d_tf) |>
  mutate(DON_0_5 = TDN_0_5 - DIN_0_5,
         DON_5_15 = TDN_5_15 - DIN_5_15)

write_csv(d_new, 'data/cp_trait_soil_data_20250128_update.csv')

din_compare <- read_csv('data/cp_trait_soil_data.csv') |>
  dplyr::rename(old_din = DIN_0_5, old_don = DON_0_5) |>
  left_join(d_new)


ggplot(din_compare, aes(x=DIN_0_5, y=old_din)) +
  geom_point() +
  xlab('new_din') +
  geom_smooth()

ggplot(din_compare, aes(x=DON_0_5, y=old_don)) +
  geom_point() +
  xlab('new_din') +
  geom_smooth()

ggsave('out/din_compare.png', bg='white', width = 3.5, height =3.5)
