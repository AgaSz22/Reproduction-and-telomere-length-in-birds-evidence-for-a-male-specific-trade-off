#===== Information =====#

### Author: Agnes Szwarczynska & Yolanda Yolanda Qian
### Date: 29/06/2026
### Aim: Build and visualise a phylogenetic tree for the change in telomere length dataset


#===== Required packages =====#

library(readr)
library(readxl)
library(brms)
library(rotl)
library(ape)
library(dplyr)
library(ggtree)
library(ggplot2)

#===== Loading data =====#


data_ta <- read_excel("../data/change_in_telomere_length.xlsx")

#===== Preparing data =====#

data_ta           <- as.data.frame(data_ta)
data_ta$Study_ID  <- as.factor(data_ta$Study_ID)

#===== Phylogenetic tree =====#

names_species    <- tnrs_match_names(data_ta$Species)
phylo_tree       <- tol_induced_subtree(ott_ids = names_species$ott_id, label_format = "name")
phylo_tree_brl   <- compute.brlen(phylo_tree)
phylo_matrix     <- vcv(phylo_tree_brl, corr = TRUE)

#===== Load model =====#

m.brm_zr_values_ta <- readRDS("../rds_objects/m.brm_zr_values_ta.rds")

# ── Extract fixed effects & CrIs ────────────────────────────────────────────

fixed_summary <- as.data.frame(summary(m.brm_zr_values_ta)$fixed)
n_species     <- nrow(fixed_summary)          

telomeres     <- fixef(m.brm_zr_values_ta)
telomeres_est <- telomeres[1:n_species, "Estimate"]   
ci_lower      <- fixed_summary[, "l-95% CI"]
ci_upper      <- fixed_summary[, "u-95% CI"]

# ── Build common_names lookup ──────────────────────────────────────────────

common_names <- unique(dplyr::select(data_ta, Species, Common_name))

# Create the new tip labels: "Common name (Latin name)"
common_names$new_tip_labels <- paste0(
  common_names$Common_name, " (", common_names$Species, ")"
)

# ── Relabel phylogenetic tree ──────────────────────────────────────────────

phylo_for_gradient           <- phylo_tree_brl
tip_species                  <- gsub("_", " ", phylo_tree_brl$tip.label)

phylo_for_gradient$tip.label <- common_names$new_tip_labels[
  match(tip_species, common_names$Species)
]

# ── Build species_est data frame ──────────────────────────────────────────

common_names <- common_names[order(common_names$Species), ]

species_est <- data.frame(
  label    = common_names$new_tip_labels,   
  Zr_est   = telomeres_est,
  ci_lower = ci_lower,
  ci_upper = ci_upper,
  row.names = NULL
)

# ── Plot ──────────────────────────────────────────────────────────────────────
tree <- ggtree(
  phylo_for_gradient,
  ladderize = FALSE,
  linewidth = 1.5
) %<+% species_est +                      
  aes(color = Zr_est) +                     
  scale_color_continuous(
    low  =  "#9ebcda",
    high = "#5f1c8c",  
    na.value = "grey50"                     
  ) +
  geom_tiplab(as_ylab = TRUE) +
  theme(
    legend.position = "bottom",
    axis.text.y     = element_text(size = 12, colour = "black")
  ) +
  labs(color = "Zr\n")

tree

# ── Save figure ───────────────────────────────────────────────────────────────

ggsave(
  filename = "../output_figures/telomere_tree_ta.svg",
  plot     = tree,
  device   = "svg",
  width    = 190,       
  height   = 148,        
  units    = "mm",
  dpi      = 300
)
