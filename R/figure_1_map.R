library(tidyverse)
library(sf)
library(terra)
library(ggnewscale)
library(elevatr)

# data in ======================================================================
plots <- read_csv('data/fraser - locations.csv') |>
  dplyr::rename(x=longitude, y=latitude)

roads <- st_read('data/grand_county_roads/')

psf <- plots |>
  st_as_sf(coords = c('x', 'y'), crs = 4326) |>
  st_buffer(500)

dem <- elevatr::get_elev_raster(psf, 14)
plot(dem); plot(psf, add=T)

dem_srtm <- terra::rast('data/dem/n39_w106_1arc_v3.bil')
  

dem_df <- dem |> as.data.frame(xy=T) |>
  filter(x < -105.935,
         x > -105.96,
         y > 39.93,
         y < 39.955) |>
  dplyr::rename(`Elevation (m)` = 3)

slope <- dem |> raster::terrain(v='slope', unit = 'radians')
aspect <- raster::terrain(dem, 'aspect' , unit = 'radians')
plot(aspect)
hs <- raster::hillShade(slope = slope, aspect = aspect) |>
  as.data.frame(xy=T) |>
  filter(x < -105.935,
         x > -105.96,
         y > 39.93,
         y < 39.955) |>
  dplyr::rename(hs = 3)
  

ggplot(dem_df) +
  geom_raster(data = hs, aes(fill=hs,x=x,y=y)) +
  scale_fill_gradient(low = 'grey10', high = "grey90", guide = 'none') +
  ggnewscale::new_scale_fill() +
  geom_raster(aes(fill = `Elevation (m)`, x=x,y=y), alpha = 0.5) +
  scale_fill_gradientn(colours = terrain.colors(10)) +
  geom_sf(data = roads, color = 'black', lwd=1) +
  geom_contour(aes(x=x,y=y,z=`Elevation (m)`), color = 'grey30', binwidth = 50, lwd = .25) +
  coord_sf(expand = F, xlim = c(-105.935, -105.96), ylim = c(39.93, 39.952)) +
  geom_point(data = plots, aes(x=x,y=y)) +
  # ggthemes::theme_clean() +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_rect(color = 'black'),
        axis.title = element_blank())

ggsave('out/map_figure.png', width = 6, height = 6, bg = 'white')

