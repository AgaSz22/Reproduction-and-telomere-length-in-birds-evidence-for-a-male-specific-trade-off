#===== Information =====#

### Author: Agnes Szwarczynska & Yolanda Yolanda Qian
### Date: 29/06/2026
### Aim: Producing a publication bias figure

#===== Required packages =====#

library(readxl)
library(ggplot2)
library(sjPlot)
library(patchwork)
library(dplyr)
library(lme4)

#===== Loading data =====#

data_tl <- read_excel("../data/telomere_length.xlsx")
data_ta <- read_excel("../data/change_in_telomere_length.xlsx")

#===== Add Zr_type label to each dataset =====#

data_tl <- data_tl %>% mutate(Zr_type = "Telomere length")
data_ta <- data_ta %>% mutate(Zr_type = "Change in telomere length")

#===== Combine into one dataframe =====#

data_combined <- bind_rows(data_tl, data_ta)

#===== Colour & shape palette for Zr type =====#

zr_colours <- c("Telomere length"    = "#9ebcda",   
                "Change in telomere length" = "#5f1c8c")   

zr_shapes  <- c("Telomere length"    = 16,           # filled circle
                "Change in telomere length" = 17)           # filled triangle

#===== Shared theme =====#

shared_theme <- theme_bw() +
  theme(
    axis.text        = element_text(size = 21),
    axis.title.x     = element_text(size = 22,
                                    margin = margin(t = 10)),  # pushes x title down away from labels
    axis.title.y     = element_text(size = 22,
                                    margin = margin(r = 10)),  # pushes y title left away from labels
    legend.position  = "bottom",
    legend.title     = element_text(size = 25, face = "bold"),
    legend.text      = element_text(size = 24),
    plot.tag         = element_text(size = 24, face = "bold"),
    panel.grid.minor = element_blank()
  )

#===== Panel a) — Zr vs Number of Individuals =====#

p_a <- ggplot(data_combined,
              aes(x      = Individuals,
                  y      = Zr,
                  colour = Zr_type,
                  shape  = Zr_type)) +
  geom_hline(yintercept = 0,
             linetype   = "dashed",
             colour     = "grey50") +
  geom_point(size  = 2.5,
             alpha = 0.85) +
  scale_colour_manual(values = zr_colours, name = "Zr type") +
  scale_shape_manual(values  = zr_shapes,  name = "Zr type") +
  xlim(0, 1300) +
  ylim(-1, 1) +
  labs(x   = "Number of Individuals",
       y   = "Zr value",
       tag = "a)") +
  shared_theme

#===== Panel b) — Zr vs Year of Publication =====#

p_b <- ggplot(data_combined,
              aes(x      = Year,
                  y      = Zr,
                  colour = Zr_type,
                  shape  = Zr_type)) +
  geom_hline(yintercept = 0,
             linetype   = "dashed",
             colour     = "grey50") +
  geom_point(size  = 2.5,
             alpha = 0.85) +
  scale_colour_manual(values = zr_colours, name = "Zr type") +
  scale_shape_manual(values  = zr_shapes,  name = "Zr type") +
  scale_x_continuous(breaks = seq(2013, 2026, by = 2),
                     limits = c(2012, 2026)) +
  ylim(-1, 1) +
  labs(x   = "Year of Publication",
       y   = "Zr value",
       tag = "b)") +
  shared_theme

#===== Panel c) — Funnel Plot (SE vs Zr) =====#

se_max    <- max(data_combined$SE_Zr, na.rm = TRUE) * 1.1
se_seq    <- seq(0, se_max, length.out = 200)

funnel_df <- data.frame(
  se    = se_seq,
  lower = -1.96 * se_seq,
  upper =  1.96 * se_seq
)

p_c <- ggplot() +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = Inf,
           fill = "#F2F0EF") +
  geom_ribbon(data = funnel_df,
              aes(x    = se,
                  ymin = lower,
                  ymax = upper),
              fill  = "white",
              alpha = 1) +
  geom_line(data     = funnel_df,
            aes(x = se, y = lower),
            linetype  = "dashed",
            colour    = "black",
            linewidth = 0.5) +
  geom_line(data     = funnel_df,
            aes(x = se, y = upper),
            linetype  = "dashed",
            colour    = "black",
            linewidth = 0.5) +
  geom_point(data  = data_combined,
             aes(x      = SE_Zr,
                 y      = Zr,
                 colour = Zr_type,
                 shape  = Zr_type),
             size  = 3,
             alpha = 0.9) +
  scale_colour_manual(values = zr_colours, name = "Zr type") +
  scale_shape_manual(values  = zr_shapes,  name = "Zr type") +
  scale_x_reverse(breaks = seq(0, 0.6, by = 0.1)) +
  coord_flip() +
  labs(x   = "Standard error",
       y   = "Zr value",
       tag = "c)") +
  shared_theme + 
  theme(legend.position = "bottom")

#===== Combine all panels =====#

p_a <- p_a + theme(legend.position = "none")   # no legend on panel a
p_b <- p_b + theme(legend.position = "none")   # no legend on panel b
p_c <- p_c + theme(legend.position = "none")

combined_plot <- (p_a + p_b) /
  p_c +
  theme(legend.position = "bottom")

#===== Save =====#

ggsave("../output_figures/publication_bias_figure.svg",
       plot   = combined_plot,
       width  = 14,
       height = 16,
       units  = "in",
       dpi    = 300)

#===== Test the relationship between year of publication and Zr value =====#

pub_year <- lmer(Zr ~ Year + (1| Study_ID ), data = data_combined)
summary(pub_year)
hist(resid(pub_year))
