library(tidyverse)
library(sf)
library(terra)
library(ggnewscale)

# data in ======================================================================
plots <- read_csv('data/fraser - locations.csv') |>
  dplyr::rename(x=longitude, y=latitude)
  
dem <- terra::rast('data/dem/n39_w106_1arc_v3.bil') |>
  as.data.frame(xy=T) |>
  filter(x < -105.935,
         x > -105.96,
         y > 39.93,
         y < 39.955) |>
  dplyr::rename(`Elevation (m)` = 3)

slope <- terra::rast('data/dem/n39_w106_1arc_v3.bil') |> terra::terrain(v='slope', unit = 'radians')
aspect <- terra::rast('data/dem/n39_w106_1arc_v3.bil')|> terra::terrain(v='aspect', unit = 'radians')

hs <- terra::shade(slope = slope, aspect = aspect) |>
  as.data.frame(xy=T) |>
  filter(x < -105.935,
         x > -105.96,
         y > 39.93,
         y < 39.955) |>
  dplyr::rename(hs = 3)
  

ggplot(dem, aes(x=x,y=y)) +
  geom_raster(data = hs, aes(fill=hs)) +
  scale_fill_gradient(low = 'grey10', high = "grey90", guide = F) +
  ggnewscale::new_scale_fill() +
  geom_raster(aes(fill = `Elevation (m)`), alpha = 0.5) +
  scale_fill_binned(type = 'viridis') +
  coord_equal(expand = F) +
  geom_point(data = plots) +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_rect(color = 'black'),
        axis.title = element_blank())

ggsave('out/map_figure.png', width = 6, height = 6, bg = 'white')

plot(dem); plot(plots, add=T)
