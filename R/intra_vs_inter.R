# intra vs inter specific trait variation
library(tidyverse)
library(cati)

d <- read_csv("data/cp_trait_soil_data.csv")
  
oreo <- filter(d, species_full == "Oreochrysum parryi"); summary(oreo)
vac <- filter(d, species_full == "Vaccinium sp"); summary(vac)
car <- filter(d, species_full == "Carex sp") |>
  filter(!is.na(sla)); summary(car)

cati::RaoRel()

?cati::partvar

traits.cp <- dplyr::select(d, sla, ldmc, height = height_cm)
factors.cp <- dplyr::select(d, plot, block, species_full) |>
  mutate_all(as.factor)

pv <- partvar(traits.cp, factors.cp)
barPartvar(pv)

df <- data.frame(name = c('plot', 'block', 'species', "within"), sla = pv[,1], ldmc = pv[,2], height = pv[,3])

df |>
  pivot_longer(cols = c('sla', 'ldmc', 'height'), names_to = 'trait') |>
  dplyr::filter(name %in% c('species', "within")) |>
  dplyr::mutate(name = str_replace_all(name, "species", 'between species'),
                name = str_replace_all(name, "within", 'within species')) |>
  ggplot(aes(x=trait, y=value, fill = name)) +
  geom_bar(stat = 'identity', color = 'black') +
  ylab("Proportion of Variance") +
  theme_bw() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.title = element_blank(),
        legend.background = element_rect(color = 'black'))
ggsave('out/trait_variance_partitioning.png', width = 4, height = 4, bg = 'white')
# example from docs
data(finch.ind)
## Not run: 
cond<-seq(1,length(sp.finch)*2, by = 2)
genus <- as.vector(unlist(strsplit(as.vector(sp.finch),"_"))[cond])

res.partvar.finch <- partvar(traits = traits.finch, 
                             factors = cbind(sites = as.factor(as.vector(ind.plot.finch)), 
                                             species = as.factor(as.vector(sp.finch)), genus = as.factor(genus)))

res.partvar.finch

oldpar<-par()
par(mfrow = c(2,2), mai = c(0.2,0.2,0.2,0.2))
piePartvar(res.partvar.finch)
par(oldpar)



# just do sds

intra <- d |>
  filter(species_full != 'Carex sp') |>
  group_by(species_full) |>
  summarise(sla_intra = sd(sla, na.rm=T)/mean(sla, na.rm=T),
            ldmc_intra = sd(ldmc)/mean(ldmc),
            height_intra = sd(height_cm)/mean(height_cm)) |>
  ungroup() |>
  mutate_if(is.numeric, round, 3)


inter <- d |>
  filter(species_full != 'Carex sp') |>
  summarise(sla_inter = sd(sla, na.rm=T)/mean(sla, na.rm=T),
                       ldmc_inter = sd(ldmc)/mean(ldmc),
                       height_inter = sd(height_cm)/mean(height_cm))
  
intra |>
  write_csv('out/coeficient_of_variations.csv')


