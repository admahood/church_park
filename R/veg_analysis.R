# Church Park vegetation comp
source("R/a_church_data_prep.R")

# nmds ===========
nmds <- comm %>%
  wisconsin() %>%
  metaMDS()

stressplot(nmds)

site_scores <- as.data.frame(vegan::scores(nmds)$sites) %>%
  as_tibble(rownames = "plot") %>%
  left_join(sites |> mutate(plot = str_c(block, treatment))); site_scores

ef <- envfit(nmds, comm, na.rm = T, permutations = 9999)
sp <-as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))
species <- as.data.frame(cbind(sp, p=ef$vectors$pvals)) %>%
  tibble::rownames_to_column("species")  %>%
  filter(p < 0.05)

adab <- adonis2(comm ~ mean_vwc + treatment + no3_mg_g + soil_p_h + c + n + 
                  doc_mg_g +tdn_mg_g +din_mg_g +don_mg_g +nh4_mg_g+ nh4_n_mg_g, data = sites_w_23_soil);adab

p_ab <- ggplot(site_scores, aes(x=NMDS1, y=NMDS2)) +
  coord_fixed() +
  geom_point(size=4, aes(color = treatment)) +
  geom_segment(data = species,x=0,y=0, color = "grey",arrow = arrow(),
               aes(yend = NMDS2, xend = NMDS1), lwd=1)+
  ggrepel::geom_text_repel(data = species,size=4, aes(label = species), color = "grey40") +
  theme_classic() +
  theme(panel.background = element_rect(fill="transparent", color = "black")) +
  ggtitle("Abundance-Based", paste("Treatement explaines", 
                                   signif(adab$R2[1], 3) * 100, 
                                   "% of the variation"));p_ab



# occurrence based ====================

nmds <- comm %>%
  decostand(method = "pa") %>%
  metaMDS()

stressplot(nmds)

site_scores <- as.data.frame(vegan::scores(nmds)$sites) %>%
  as_tibble(rownames = "plot") %>%
  left_join(sites); site_scores

ef <- envfit(nmds, comm %>%
               decostand(method = "pa"), na.rm = T, permutations = 9999)
eft <- envfit(nmds, sites[,2], na.rm = T, permutations = 9999)
sp <-as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))
species <- as.data.frame(cbind(sp, p=ef$vectors$pvals)) %>%
  tibble::rownames_to_column("species")  %>%
  filter(p < 0.01)

adoc <- adonis2(comm %>% decostand("pa") ~ treatment, data = sites)

p_oc <- ggplot(site_scores, aes(x=NMDS1, y=NMDS2)) +
  coord_fixed() +
  geom_point(size=4, aes(color = treatment)) +
  geom_segment(data = species,x=0,y=0, color = "grey",arrow = arrow(),
               aes(yend = NMDS2, xend = NMDS1), lwd=1)+
  ggrepel::geom_text_repel(data = species,size=4, aes(label = species), color = "grey40") +
  theme_classic() +
  theme(panel.background = element_rect(fill="transparent", color = "black")) +
  ggtitle("Occurrence-Based", paste("Treatement explaines", 
                                    signif(adoc$R2[1], 3) * 100, 
                                    "% of the variation"));p_oc


ggpubr::ggarrange(p_ab, p_oc, common.legend = TRUE, widths = c(1,1.15)) %>%
  ggsave(width =10, height =6, filename = "out/nmds.png", plot = ., bg="white")


# biodiversity indexes

diversity(comm, groups = sites$treatment)
df <- data.frame(turnover = NA, nestedness = NA, jaccard = NA)
df[1,] <- nestedbetajac(comm[sites$treatment == "b"])
df[2,] <- nestedbetajac(comm[sites$treatment == "bm"])
df[3,] <- nestedbetajac(comm[sites$treatment == "c"])
df[4,] <- nestedbetajac(comm[sites$treatment == "m"])
df$trt <- c("b", "bm", "c", "m")

ggplot(df, aes(x = turnover, y=trt)) +
  geom_bar(stat = "identity")
ggplot(df, aes(x = nestedness, y=trt)) +
  geom_bar(stat = "identity")
ggplot(df, aes(x = jaccard, y=trt)) +
  geom_bar(stat = "identity")

pdiv <- sites %>%
 mutate(shannon = diversity(comm, index = "shannon"),
        simpson = diversity(comm, index = "simpson"),
        invsimpson = diversity(comm, index = "invsimpson"),
        richness = vegan::specnumber(comm)) %>%
  pivot_longer(cols = c("shannon", "simpson", "invsimpson", "richness")) %>%
  ggplot(aes(x=treatment, y=value, fill = treatment)) +
  geom_boxplot() +
  facet_wrap(~name, scales = "free_y") ;pdiv

ggsave(plot = pdiv, filename = "out/alphadiversity_x_treatment.png", width =5, height = 5)
