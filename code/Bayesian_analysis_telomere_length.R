#===== Information =====#

### Author: Agnes Szwarczynska & Yolanda Yolanda Qian
### Date: 29/06/2026
### Aim: Run Bayesian models using the telomere length dataset and extract data neccessary to build a phylogenetic tree

#===== Required packages =====#

library(readr)
library(readxl)   # needed for read_excel
library(brms)
library(rotl)
library(ape)
library(dplyr)

#===== Loading data =====#

data_tl <- read_excel("../data/telomere_length.xlsx")

#===== Preparing data =====#

data_tl <- as.data.frame(data_tl)
data_tl$Study_ID <- as.factor(data_tl$Study_ID)

#===== Priors =====#

# Half-Cauchy achieved by setting lb = 0
priors_main <- c(
  set_prior("normal(0, 0.6)", class = "Intercept", lb = -1, ub = 1),
  set_prior("cauchy(0, 0.5)", class = "sd", lb = 0)   # Half-Cauchy for tau
)

priors <- c(
  set_prior("normal(0, 0.6)", class = "Intercept", lb = -1, ub = 1),
  set_prior("normal(0, 0.6)", class = "b"),
  set_prior("cauchy(0, 0.5)", class = "sd", lb = 0)   # Half-Cauchy for tau
)

#===== Bayesian models =====#

### Intercept-only model

m.brm_intercept <- brm(
  Zr | se(SE_Zr) ~ 1 + (1 | Study_ID),
  family  = gaussian(link = "identity"),
  data    = data_tl,           
  prior   = priors_main,
  init    = 0,
  iter    = 40000,
  warmup  = 1000,
  chains  = 4,                 
  thin    = 1,                 
  control = list(
    adapt_delta   = 0.95,      
    max_treedepth = 15         
  ),
  seed = 14
)

saveRDS(m.brm_intercept, file = "../rds_objects/m.brm_intercept.rds")
sum_m.brm_intercept <- summary(m.brm_intercept)   
saveRDS(sum_m.brm_intercept, file = "../rds_objects/sum_m.brm_intercept.rds")

### Model with all moderators

m.brm_moderators <- brm(
  Zr | se(SE_Zr) ~ 1 + Sex + Age_accounted + Method + Long_Cross + Stage + Proxy + (1 | Study_ID),
  family  = gaussian(link = "identity"),
  data    = data_tl,           
  prior   = priors,
  init    = 0,
  iter    = 40000,
  warmup  = 1000,
  chains  = 4,                 
  thin    = 1,                 
  control = list(
    adapt_delta   = 0.95,      
    max_treedepth = 15         
  ),
  seed = 14
)

saveRDS(m.brm_moderators, file = "../rds_objects/m.brm_moderators.rds")
sum_m.brm_moderators <- summary(m.brm_moderators)   
saveRDS(sum_m.brm_moderators, file = "../rds_objects/sum_m.brm_moderators.rds")

### Model with "Sex" and "Age" moderator

m.brm_sig_moderators <- brm(
  Zr | se(SE_Zr) ~ 1 + Sex + Age_accounted + (1 | Study_ID),
  family  = gaussian(link = "identity"),
  data    = data_tl,           
  prior   = priors,
  init    = 0,
  iter    = 40000,
  warmup  = 1000,
  chains  = 4,                 
  thin    = 1,                 
  control = list(
    adapt_delta   = 0.95,      
    max_treedepth = 15         
  ),
  seed = 14
)

saveRDS(m.brm_sig_moderators, file = "../rds_objects/m.brm_sig_moderators.rds")
sum_m.brm_sig_moderators <- summary(m.brm_sig_moderators)   
saveRDS(sum_m.brm_sig_moderators, file = "../rds_objects/sum_m.brm_sig_moderators.rds")

### Model with "Sex" and "Age" moderator: sensitivity analysis (removing studies that use hormones as proxies)

data_tl_sens <- data_tl %>% filter(!Study_ID %in% c(18,20))

m.brm_sig_moderators_sens <- brm(
  Zr | se(SE_Zr) ~ 1 + Sex + Age_accounted + (1 | Study_ID),
  family  = gaussian(link = "identity"),
  data    = data_tl_sens,           
  prior   = priors,
  init    = 0,
  iter    = 40000,
  warmup  = 1000,
  chains  = 4,                 
  thin    = 1,                 
  control = list(
    adapt_delta   = 0.95,      
    max_treedepth = 15         
  ),
  seed = 14
)

saveRDS(m.brm_sig_moderators_sens, file = "../rds_objects/m.brm_sig_moderator_sens.rds")
sum_m.brm_sig_moderators_sens <- summary(m.brm_sig_moderators_sens)   
saveRDS(sum_m.brm_sig_moderators_sens, file = "../rds_objects/sum_m.brm_sig_moderators_sens.rds")

#===== Bayesian phylogenetic models =====#

### Building a phylogenetic tree

names_species <- tnrs_match_names(data_tl$Species)  # extracting species names
phylo_tree <- tol_induced_subtree(ott_ids = names_species$ott_id, label_format = "name")
plot(phylo_tree, no.margin = TRUE)

# svg("Phylo_tree.svg", width = 30, height = 24)
# plot(phylo_tree, no.margin = TRUE)
# dev.off()

### Model with a phylogenetic tree as a moderator

phylo_tree_brl <- compute.brlen(phylo_tree)           # computing branch lengths of a tree
phylo_matrix  <- vcv(phylo_tree_brl, corr = TRUE)     # building distance matrix
data_tl$obs   <- 1:nrow(data_tl)  

m.brm_phylogenetic <- brm(
  Zr | se(SE_Zr) ~ 1 + (1 | Study_ID) + (1 | gr(Species2, cov = phylo_matrix)) + (1 | obs),
  family   = gaussian(link = "identity"),
  data     = data_tl,                                  
  data2    = list(phylo_matrix = phylo_matrix),
  prior    = priors_main, #using default priors doesn't help with divergent transitions
  control  = list(adapt_delta = 0.99,
                  max_treedepth = 15),
  init     = 0,
  iter     = 40000,
  warmup   = 1000,
  seed     = 14,
  thin     = 1
)

saveRDS(m.brm_phylogenetic, file = "../rds_objects/m.brm_phylogenetic.rds")
sum_m.brm_phylogenetic <- summary(m.brm_phylogenetic)   
saveRDS(sum_m.brm_phylogenetic, file = "../rds_objects/sum_m.brm_phylogenetic.rds")

### Extracting values for the phylogenetic tree figure
m.brm_zr_values <- brm(Zr|se(SE_Zr) ~ 1 + Species + (1|Study_ID),
                  family = gaussian(link = "identity"),
                  data = data_tl,
                  prior = priors,
                  control = list(adapt_delta = 0.99),
                  init = 0,
                  iter = 40000,
                  warmup = 1000,
                  seed=14,
                  thin=1)

saveRDS(m.brm_zr_values, file = "../rds_objects/m.brm_zr_values.rds")
sum_m.brm_phylogenetic_values <- summary(m.brm_zr_values)
saveRDS(sum_m.brm_phylogenetic_values, file = "../rds_objects/sum_m.brm_phylogenetic_values.rds")
