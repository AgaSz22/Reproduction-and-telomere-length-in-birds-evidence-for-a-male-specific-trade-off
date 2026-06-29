#===== Information =====#

### Author: Agnes Szwarczynska & Yolanda Yolanda Qian
### Date: 29/06/2026
### Aim: Build and visualise forest plot for studies on telomere length (Fig.3)

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

#===== Model =====#

# m.brm_author <- brm(
#   Zr | se(SE_Zr) ~ 1 + (1 | Author),
#   family  = gaussian(link = "identity"),
#   data    = data_tl,           
#   prior   = priors_main,
#   init    = 0,
#   iter    = 40000,
#   warmup  = 1000,
#   chains  = 4,                 
#   thin    = 1,                 
#   control = list(
#     adapt_delta   = 0.95,      
#     max_treedepth = 15         
#   ),
#   seed = 14
# )

# saveRDS(m.brm_author, file = "../rds_objects/m.brm_author.rds")
# sum_m.brm_author <- summary(m.brm_author)   
# saveRDS(sum_m.brm_author, file = "../rds_objects/sum_m.brm_author.rds")

m.brm_author <- readRDS(file = "../rds_objects/m.brm_author.rds")

# --- Pre-compute releveled factors ---

study.draws <- spread_draws(m.brm_author, r_Author[Author,], b_Intercept) %>% 
  mutate(b_Intercept = r_Author + b_Intercept)

pooled.effect.draws <- spread_draws(m.brm_author, b_Intercept) %>% 
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
  
  geom_vline(xintercept = fixef(m.brm_author)[1, 1],
             color = "#9ebcda", linewidth = 1) +
  geom_vline(xintercept = fixef(m.brm_author)[1, 3:4],
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

#===== Save as SVG at 1/3 A4 =====#  

ggsave("../output_figures/forest_plot_telomere_length.svg",
       plot = last_plot(),
       device = "svg",
       width = 8.3,
       height = 5.9,
       units = "in",
       dpi = 300)

#===== Posterior predictive checks =====#  

pp_check(m.brm_author, ndraws = 100)
plot(m.brm_author)

#===== Calculate the I^2 statistics =====#  

m.rma <- rma.mv(yi     = Zr,
                V      = SE_Zr^2,       # variance = SE²
                random = ~ 1 | Author,
                method = "REML",
                data   = data_tl)

i2_ml(m.rma)

