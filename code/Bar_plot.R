# ==========================================
# Information
# ==========================================

### Author: Agnes Szwarczynska & Yolanda Qian
### Date: 29/06/2026
### Aim: Visualise descriptive statistics (Fig.2)

# ==========================================
# LOAD libraries
# ==========================================

library(readxl)

# ==========================================
# LOAD DATASETS
# ==========================================

data_tl <- read_excel("../data/telomere_length.xlsx")
data_ta <- read_excel("../data/change_in_telomere_length.xlsx")

# ==========================================
# PLOTTING FUNCTION
# ==========================================

plot_custom_bar <- function(raw_column, x_label, tag, y_max = NULL) {
  
  class_counts <- table(raw_column)
  
  if (is.null(y_max)) y_max <- max(class_counts) * 1.4
  
  par(mar = c(5, 5, 3, 2))
  
  opar <- par(lwd = 2.5)
  
  class_barplot <- barplot(class_counts,
                           ylab = "Frequency",
                           xlab = x_label,
                           cex.names = 1.2,
                           cex.axis = 1.2,
                           cex.lab = 1.4,
                           col = "black",
                           ylim = c(0, y_max),
                           space = 0.4,
                           main = ""
  )
  
  text(x = class_barplot, y = class_counts + y_max * 0.08, labels = class_counts,
       cex = 1.5, col = "black", font = 2)
  
  mtext(tag, side = 3, adj = -0.1, line = 1.2, cex = 1.2, font = 2)
  
  par(opar)
}

# ==========================================
# LAYOUT & EXECUTION (A4 Setup)
# ==========================================

# --- Reusable plotting function ---
draw_figure <- function() {
  
  layout_matrix <- matrix(c(
    1, 1, 1, 1, 1, 1,
    2, 2, 3, 3, 4, 4,
    5, 5, 6, 6, 7, 7,
    8, 8, 8, 8, 8, 8,
    9, 9, 9, 10,10,10,
    11,11,11, 12,12,12
  ), nrow = 6, byrow = TRUE)
  
  layout(layout_matrix, heights = c(0.15, 1, 1, 0.15, 1, 1))
  
  # --- Plot 1: Title I ---
  par(mar = c(0, 4, 0, 0))
  plot.new()
  text(x = 0, y = 0.5, "I)", cex = 2, font = 2, adj = 0)
  
  # --- Plots 2-7: Section I (shared y_max = 80) ---
  plot_custom_bar(data_tl$Proxy,        "Proxy",                   "(a)", y_max = 80)
  plot_custom_bar(data_tl$Sex,          "Sex",                     "(b)", y_max = 80)
  plot_custom_bar(data_tl$Stage,        "Age at measurement",      "(c)", y_max = 80)
  plot_custom_bar(data_tl$Method,       "Laboratory method",       "(d)", y_max = 80)
  plot_custom_bar(data_tl$Age_accounted,"Accounted for age",       "(e)", y_max = 80)
  plot_custom_bar(data_tl$Long_Cross,   "Sample collection design","(f)", y_max = 80)
  
  # --- Plot 8: Title II ---
  par(mar = c(0, 4, 0, 0))
  plot.new()
  text(x = 0, y = 0.5, "II)", cex = 2, font = 2, adj = 0)
  
  # --- Plots 9-12: Section II (shared y_max = 25) ---
  plot_custom_bar(data_ta$Proxy,        "Proxy",             "(a)", y_max = 25)
  plot_custom_bar(data_ta$Sex,          "Sex",               "(b)", y_max = 25)
  plot_custom_bar(data_ta$Method,       "Laboratory method", "(c)", y_max = 25)
  plot_custom_bar(data_ta$Age_accounted,"Accounted for age", "(d)", y_max = 25)
  
}

# --- Save as SVG---
svg("../output_figures/Bar_plot_descriptive.svg", width = 8.27, height = 11.69)
draw_figure()
dev.off()

