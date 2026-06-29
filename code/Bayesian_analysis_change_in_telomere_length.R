#===== Information =====#

### Author: Agnes Szwarczynska & Yolanda Yolanda Qian
### Date: 29/06/2026
### Aim: Run Bayesian models using the change in telomere length dataset and extract data neccessary to build a phylogenetic tree


#===== Required packages =====#

library(readr)
library(readxl)   # needed for read_excel
library(brms)
library(rotl)
library(ape)
library(dplyr)

#===== Loading data =====#

data_ta <- read_excel("../data/change_in_telomere_length.xlsx")

#===== Preparing data =====#

data_ta <- as.data.frame(data_ta)
data_ta$Study_ID <- as.factor(data_ta$Study_ID)

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

### Intercept-only model

m.brm_intercept_ta <- brm(
  Zr | se(SE_Zr) ~ 1 + (1 | Study_ID),
  family  = gaussian(link = "identity"),
  data    = data_ta,           
  prior   = priors_main,
  init    = 0,
  iter    = 40000,
  warmup  = 1000,
  chains  = 4,                 
  thin    = 1,                 
  control = list(
    adapt_delta   = 0.99,      
    max_treedepth = 15         
  ),
  seed = 14
)

saveRDS(m.brm_intercept_ta, file = "../rds_objects/m.brm_intercept_ta.rds")
sum_m.brm_intercept_ta <- summary(m.brm_intercept_ta)   
saveRDS(sum_m.brm_intercept_ta, file = "../rds_objects/sum_m.brm_intercept_ta.rds")

### Model with all moderators

m.brm_moderators_ta <- brm(
  Zr | se(SE_Zr) ~ 1 + Sex + Age_accounted + Method + Proxy + (1 | Study_ID),
  family  = gaussian(link = "identity"),
  data    = data_ta,           
  prior   = priors,
  init    = 0,
  iter    = 40000,
  warmup  = 1000,
  chains  = 4,                 
  thin    = 1,                 
  control = list(
    adapt_delta   = 0.99,      
    max_treedepth = 15         
  ),
  seed = 14
)

saveRDS(m.brm_moderators_ta, file = "../rds_objects/m.brm_moderators_ta.rds")
sum_m.brm_moderators_ta <- summary(m.brm_moderators_ta)   
saveRDS(sum_m.brm_moderators_ta, file = "../rds_objects/sum_m.brm_moderators_ta.rds")

### Model with all moderators : sensitivity analysis

data_ta_sens <- data_ta %>% filter(!Study_ID %in% c(1,10)) #conducting sensitivity analysis is impossible in this case, because we are left with only 1 study in the "indirect" proxy group

# m.brm_moderators_ta_sens <- brm(
#   Zr | se(SE_Zr) ~ 1 + Sex + Age_accounted + Method + Proxy + (1 | Study_ID),
#   family  = gaussian(link = "identity"),
#   data    = data_ta_sens,           
#   prior   = priors,
#   init    = 0,
#   iter    = 40000,
#   warmup  = 1000,
#   chains  = 4,                 
#   thin    = 1,                 
#   control = list(
#     adapt_delta   = 0.99,      
#     max_treedepth = 15         
#   ),
#   seed = 14
# )
# 
# saveRDS(m.brm_moderators_ta_sens, file = "../rds_objects/m.brm_moderators_ta_sens.rds")
# sum_m.brm_moderators_ta_sens <- summary(m.brm_moderators_ta_sens)   
# saveRDS(sum_m.brm_moderators_ta_sens, file = "../rds_objects/sum_m.brm_moderators_ta_sens.rds")

#===== Bayesian phylogenetic models =====#

### Building a phylogenetic tree

names_species <- tnrs_match_names(data_ta$Species)  ### extracting species names
phylo_tree <- tol_induced_subtree(ott_ids = names_species$ott_id, label_format = "name")
plot(phylo_tree, no.margin = TRUE)

# svg("Phylo_tree_ta.svg", width = 30, height = 24)
# plot(phylo_tree, no.margin = TRUE)
# dev.off()

### Model with a phylogenetic tree as a moderator

phylo_tree_brl <- compute.brlen(phylo_tree)           # computing branch lengths of a tree
phylo_matrix  <- vcv(phylo_tree_brl, corr = TRUE)     # building distance matrix
data_ta$obs   <- 1:nrow(data_ta)  

m.brm_phylogenetic_ta <- brm(
  Zr | se(SE_Zr) ~ 1 + (1 | Study_ID) + (1 | gr(Species2, cov = phylo_matrix)) + (1 | obs),
  family   = gaussian(link = "identity"),
  data     = data_ta,                                  
  data2    = list(phylo_matrix = phylo_matrix),
  prior    = priors_main,
  control  = list(adapt_delta = 0.99,
                  max_treedepth = 15),
  init     = 0,
  iter     = 40000,
  warmup   = 1000,
  seed     = 14,
  thin     = 1
)

saveRDS(m.brm_phylogenetic_ta, file = "../rds_objects/m.brm_phylogenetic_ta.rds")
sum_m.brm_phylogenetic_ta <- summary(m.brm_phylogenetic_ta)   
saveRDS(sum_m.brm_phylogenetic_ta, file = "../rds_objects/sum_m.brm_phylogenetic_ta.rds")

### Extracting values for the phylogenetic tree figure

m.brm_zr_values_ta <- brm(Zr|se(SE_Zr) ~ 1 + Species + (1|Study_ID),
                       family = gaussian(link = "identity"),
                       data = data_ta,
                       prior = priors,
                       control = list(adapt_delta = 0.99),
                       init = 0,
                       iter = 40000,
                       warmup = 1000,
                       seed=14,
                       thin=1)

saveRDS(m.brm_zr_values_ta, file = "../rds_objects/m.brm_zr_values_ta.rds")
sum_m.brm_zr_values_ta <- summary(m.brm_zr_values_ta)
saveRDS(sum_m.brm_zr_values_ta, file = "../rds_objects/sum_m.brm_zr_values_ta.rds")
