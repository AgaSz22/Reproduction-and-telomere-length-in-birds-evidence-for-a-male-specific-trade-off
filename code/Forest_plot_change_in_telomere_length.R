#===== Information =====#

### Author: Agnes Szwarczynska & Yolanda Qian
### Date: 29/06/2026
### Aim: Build and visualise forest plot for studies on change in telomere length (Fig.5)

#===== Libraries =====#

library(readr)
library(brms)
library(ggplot2)
library(tidybayes)
library(dplyr)
library(glue)
library(stringr)
library(forcats)
library(metafor)
library(orchaRd)

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

#===== Model =====#

### Intercept-only model

# m.brm_author_ta <- brm(
#   Zr | se(SE_Zr) ~ 1 + (1 | Author),
#   family  = gaussian(link = "identity"),
#   data    = data_ta,
#   prior   = priors_main,
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
# saveRDS(m.brm_author_ta, file = "../rds_objects/m.brm_author_ta.rds")
# sum_m.brm_author_ta <- summary(m.brm_author_ta)
# saveRDS(sum_m.brm_author_ta, file = "../rds_objects/sum_m.brm_author_ta.rds")

m.brm_author_ta <- readRDS(file = "../rds_objects/m.brm_author_ta.rds")

# --- Pre-compute releveled factors ---

study.draws <- spread_draws(m.brm_author_ta, r_Author[Author,], b_Intercept) %>% 
  mutate(b_Intercept = r_Author + b_Intercept)

pooled.effect.draws <- spread_draws(m.brm_author_ta, b_Intercept) %>% 
  mutate(Author = "Pooled Effect")

forest.data <- bind_rows(study.draws, 
                         pooled.effect.draws) %>% 
  ungroup() %>%
  mutate(Author = str_replace_all(Author, "[.]", " ")) %>% 
  mutate(Author = reorder(Author, b_Intercept)) %>%
  mutate(Author = relevel(factor(Author), "Pooled Effect", after = Inf))

forest.data.summary <- group_by(forest.data, Author) %>% 
  mean_qi(b_Intercept)

# --- Plot ---
ggplot(aes(b_Intercept, Author),
       data = forest.data) +
  
  geom_vline(xintercept = fixef(m.brm_author_ta)[1, 1],
             color = "#9ebcda", linewidth = 1) +
  geom_vline(xintercept = fixef(m.brm_author_ta)[1, 3:4],
             color = "black", linetype = 2) +
  geom_vline(xintercept = 0, color = "black",
             linewidth = 1) +
  
  geom_pointinterval(data = forest.data.summary,
                     linewidth = 1,
                     aes(xmin =.lower, xmax =.upper)) +
  
  geom_text(data = mutate(forest.data.summary,
                          across(where(is.numeric), \(x) round(x, 2)),
                          label_text = glue("{b_Intercept} [{.lower}, {.upper}]")),
            aes(label = label_text, x = Inf),
            hjust = "inward") +
  
  labs(x = "Zr value",
       y = NULL) +
  
  theme_minimal() +
  theme(axis.text.y  = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 16, color = "black"),
        axis.text.x  = element_text(size = 14, color = "black")) +
  
  coord_cartesian(xlim = c(-0.8, 0.8))

# --- Save as SVG at 1/3 A4 ---
ggsave("../output_figures/forest_plot_change_in_telomere_length.svg",
       plot = last_plot(),
       device = "svg",
       width = 10.7,
       height = 5.9,
       units = "in",
       dpi = 300)

#===== Posterior predictive checks =====#  

pp_check(m.brm_author_ta, ndraws = 100)
plot(m.brm_author_ta)

#===== Calculate the I^2 statistics =====#  

m.rma <- rma.mv(yi     = Zr,
                V      = SE_Zr^2,       # variance = SE²
                random = ~ 1 | Author,
                method = "REML",
                data   = data_ta)

i2_ml(m.rma)
