source("R/a_church_data_prep.R")

library(Hmsc)

# load wide community data
Y <- comm
# load ancillary data

XData <- sites_w_23_soil |>
  mutate_if(is.character, as.factor)

# model formula

XFormula <- ~ 
  treatment +
  mean_vwc +
  c +
  n +
  soil_p_h +
  doc_mg_g +
  percent_sand_0to5 +
  nh4_n_mg_g +
  no3_n_mg_g +
  na_mg_g +
  k_mg_g +
  mg_mg_g +
  ca_mg_g +
  cl_mg_g +
  po4_mg_g +
  so4_mg_g +
  shannon_bacteria_15 + 
  shannon_bacteria_5 +
  shannon_fungi_5 +
  shannon_fungi_15

# random level setup

studyDesign <- data.frame(block = (XData$block),
                          plot = (XData$plot))
rLb <- HmscRandomLevel(units = studyDesign$block)
rLp <- HmscRandomLevel(units = studyDesign$plot)


mod = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr="probit",
           #TrData = traits,
           #TrFormula = tr_form,
           studyDesign = studyDesign,
           ranLevels = list("block" = rLb, "plot" = rLp))

# mcmc setup ===================================================================

nChains = 4
test.run = TRUE 

if (test.run){
  #with this option, the vignette evaluates in ca. 1 minute in adam's laptop
  thin = 10
  samples = 100
  transient = ceiling(thin*samples*.5)
}else{
  # with a spatial random effect, evaluates in 20 minutes on adam's laptop
  thin = 10
  samples = 1000
  transient = ceiling(thin*samples*.5)
}

hmsc_file <- str_c("data/hmsc/hmsc_probit", 
                   str_c(
                     "_",thin * samples + transient, "_iterations"
                   )
                   ,".rda") |>
  str_replace_all("000000", "M") |> 
  str_replace_all("000", "K")

# run hmsc =====================================================================

# Sodium by treatment?

t0 <- Sys.time()
print(t0)
if(!dir.exists("data/hmsc")) dir.create("data/hmsc")
if(!file.exists(hmsc_file)){
  m = sampleMcmc(mod, thin = thin, 
                 samples = samples, 
                 transient = transient,
                 nChains = nChains, 
                 nParallel = nChains)
  print(Sys.time()-t0)
  save(m, file=hmsc_file)
}else{load(hmsc_file)}

# visualize ===========
library(gghmsc)
vp_cols <- c(
  RColorBrewer::brewer.pal(9, "YlOrRd"),
  RColorBrewer::brewer.pal(8, "Blues")[c(1:4)],
  RColorBrewer::brewer.pal(5, "Greys")[c(3,4)],
  RColorBrewer::brewer.pal(8, "Blues")[c(5:8)],
  RColorBrewer::brewer.pal(3, "Purples")
)

gghmsc::gghmsc_convergence(m)
gghmsc::gghmsc_vp(m, cols = vp_cols)
gghmsc::gghmsc_beta(m)
gghmsc::gghmsc_omega(m)
