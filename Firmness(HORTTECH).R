# ============================================================
# Blackberry firmness & timing — FINAL (pub style, overwrites)
# - Boxplots by genotype (FF / TA.XT / SUBJ) with Tukey letters
# - Correlation "cowplots" (FF, SUBJ, TA.XT) — no title, smaller method labels
# - Efficiency bars (time per berry) with mean labels
# - NEW: Time ANOVA (time ~ type) + Tukey letters + time boxplots by method
# - High-res PNGs saved to Desktop/Blackberry_Firmness_Exports
# ============================================================

rm(list = ls()); options(stringsAsFactors = FALSE)

# ---- Packages ----
pkgs <- c(
  "dplyr","tidyr","ggplot2","cowplot","agricolae","stringr","readr",
  "tibble","grid","RColorBrewer","scales","broom"
)
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- File paths (EDIT if needed) ----
file_2023 <- "C:/Users/RhysB/OneDrive/Desktop/2023FirmnessData.csv"
file_2024 <- "C:/Users/RhysB/OneDrive/Desktop/2024FirmnessData.csv"
"C:\Users\RhysB\OneDrive\Desktop\Y1_GWAS_Collection_List.xlsx"
# ---- Desktop export folder (robust to OneDrive) ----
desktop_candidates <- c(
  file.path(Sys.getenv("USERPROFILE"), "OneDrive", "Desktop"),
  file.path(Sys.getenv("USERPROFILE"), "Desktop")
)
desktop_path <- desktop_candidates[dir.exists(desktop_candidates)][1]
out_dir <- file.path(desktop_path, "Blackberry_Firmness_Exports")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Helpers ----
`%||%` <- function(a,b) if (is.null(a)) b else a

read_firm <- function(path, year){
  df <- read.csv(path, check.names = FALSE)
  keep <- intersect(names(df), c("seln","harv_no","rep","hdate","evaldate","type","time","force","meansep"))
  df <- df[, keep, drop = FALSE]
  df %>%
    mutate(
      seln     = factor(seln),
      harv_no  = suppressWarnings(factor(harv_no)),
      rep      = suppressWarnings(factor(rep)),
      type     = factor(tolower(type), levels = c("ff","taxt","subj")),
      time     = suppressWarnings(as.numeric(time)),
      force    = suppressWarnings(as.numeric(force)),
      year     = year
    )
}

dat23 <- read_firm(file_2023, 2023)
dat24 <- read_firm(file_2024, 2024)

# ---- Robust palette for 10+ genotypes ----
genotype_palette <- function(levels_vec){
  n <- length(levels_vec)
  if (n <= 12) {
    RColorBrewer::brewer.pal(max(3, n), "Set3")
  } else {
    scales::hue_pal(l = 65, c = 100)(n)
  }
}

# ---- HSD letters for force ~ seln (per method/year) ----
hsd_letters_force <- function(d) {
  if (nrow(d) == 0) return(NULL)
  fit <- aov(force ~ seln, data = d)
  out <- agricolae::HSD.test(fit, "seln", group = TRUE)$groups
  tibble(seln = rownames(out),
         mean = out$means %||% out$`means`,
         letters = out$groups) %>% arrange(seln)
}

# Y position for force letters
box_y_positions <- function(d) {
  d %>%
    group_by(seln) %>%
    summarise(y_max = max(force, na.rm = TRUE), .groups = "drop") %>%
    mutate(y_lab = y_max + 0.12 * (max(y_max, na.rm = TRUE)))
}

# ---- Themes ----
theme_gray_pub <- function(base_size = 14){
  theme_gray(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major = element_line(size = 0.4, colour = "white"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold")
    )
}
theme_classic_pub <- function(base_size = 14){
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold")
    )
}

# ---- Box-and-whisker (colored; tight margins) ----
plot_box_with_letters <- function(d, title, ylab,
                                  show_x_title = FALSE,
                                  show_legend  = FALSE){
  if (nrow(d) == 0) return(NULL)
  letters <- hsd_letters_force(d)
  ypos    <- box_y_positions(d)
  lab_df  <- dplyr::left_join(letters, ypos, by = "seln")
  
  pal_map <- setNames(genotype_palette(levels(d$seln)), levels(d$seln))
  
  y_max  <- max(d$force, na.rm = TRUE)
  offset <- 0.15 * y_max
  
  ggplot(d, aes(x = seln, y = force, fill = seln)) +
    geom_boxplot(outlier.shape = 16, width = 0.72, size = 0.6, fatten = 2) +
    geom_text(data = lab_df,
              aes(x = seln, y = y_lab, label = letters),
              vjust = -0.15, size = 4.8, inherit.aes = FALSE) +
    scale_fill_manual(values = pal_map, name = "Genotype") +
    labs(x = if (show_x_title) "Genotype" else NULL,
         y = ylab, title = title) +
    theme_gray_pub(base_size = 14) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      plot.margin     = margin(t = 6, r = 4, b = if (show_x_title) 4 else 0, l = 6),
      legend.position = if (show_legend) "right" else "none",
      legend.key.size = grid::unit(0.70, "cm"),
      legend.text     = element_text(size = 13),
      legend.title    = element_text(size = 14, face = "bold")
    ) +
    expand_limits(y = y_max + offset) +
    coord_cartesian(ylim = c(0, y_max + offset * 1.6), clip = "off")
}

# ---- Efficiency bars (time) with mean labels (bigger fonts) ----
plot_time_bar <- function(df_year, year){
  stats <- df_year %>%
    group_by(type) %>%
    summarise(
      n        = n(),
      avg_time = mean(time, na.rm = TRUE),
      sd_time  = sd(time,   na.rm = TRUE),
      se_time  = sd_time / sqrt(n),
      .groups = "drop"
    ) %>%
    mutate(
      type_nice = factor(recode(tolower(type),
                                "ff"="FF","taxt"="TA.XT","subj"="SUBJ"),
                         levels = c("FF","TA.XT","SUBJ"))
    )
  
  y_max <- max(stats$avg_time + stats$se_time, na.rm = TRUE) * 1.05
  
  ggplot(stats, aes(x = type_nice, y = avg_time)) +
    geom_col(width = 0.6, fill = "#BFD7EA", colour = "black", linewidth = 0.5) +
    geom_errorbar(aes(ymin = avg_time - se_time, ymax = avg_time + se_time),
                  width = 0.28, linewidth = 1.2, colour = "black") +
    geom_text(aes(label = sprintf("%.1f s", avg_time), y = avg_time + se_time + 0.5),
              size = 5.2, fontface = "bold", hjust = 0) +
    coord_flip() +
    labs(x = "Method", y = "Seconds per berry (s)",
         title = paste0("Efficiency by Method (", year, ")")) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 18),
      axis.text  = element_text(size = 14),
      plot.margin = margin(8, 16, 8, 8)
    ) +
    expand_limits(y = y_max)
}

# ---- Correlation panels (no main title; smaller method labels) ----
plot_corr_pairs <- function(df_year, year, outdir, save_outputs,
                            method_label_cex = 1.4, r_cex = 1.6){
  # genotype means wide
  M <- df_year %>%
    group_by(seln, type) %>%
    summarise(force = mean(force, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = type, values_from = force) %>%
    dplyr::select(ff, subj, taxt)
  
  if (!nrow(M) || ncol(M) < 2) return(invisible(NULL))
  
  panel_cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0,1,0,1))
    r <- suppressWarnings(cor(x, y, use = "pairwise.complete.obs"))
    text(0.5, 0.5, paste0("R = ", sprintf("%.2f", r)),
         cex = r_cex, font = 2)
  }
  panel_pts <- function(x, y) points(x, y, pch = 19, cex = 1.2)
  
  # Save PNG (no main title). Method labels trimmed font size.
  fig_file <- file.path(outdir, paste0("Fig2_Correlation_", year, ".png"))
  if (save_outputs) {
    if (file.exists(fig_file)) file.remove(fig_file)
    png(fig_file, width = 1800, height = 1800, res = 200)
    pairs(M,
          labels = c("FruitFirm 1000","Subjective","TA.XT"),
          lower.panel = panel_cor, upper.panel = panel_pts,
          cex.labels = method_label_cex, font.labels = 2)
    dev.off()
    message("Saved: ", normalizePath(fig_file))
  } else {
    pairs(M,
          labels = c("FruitFirm 1000","Subjective","TA.XT"),
          lower.panel = panel_cor, upper.panel = panel_pts,
          cex.labels = method_label_cex, font.labels = 2)
  }
}

# ---- Yearly boxplot panel with legend column (centered near TA.XT) ----
assemble_year_panel <- function(df_year, year){
  ff   <- df_year %>% filter(type == "ff")
  taxt <- df_year %>% filter(type == "taxt")
  subj <- df_year %>% filter(type == "subj")
  
  p_ff   <- plot_box_with_letters(ff,   paste0("FruitFirm 1000 (", year, ")"), "Force (g/mm)",
                                  show_x_title = FALSE, show_legend = FALSE)
  p_taxt <- plot_box_with_letters(taxt, paste0("TA.XT (",          year, ")"), "Force (g)",
                                  show_x_title = FALSE, show_legend = FALSE)
  p_subj <- plot_box_with_letters(subj, paste0("Subjective (",     year, ")"), "Firmness (1–5)",
                                  show_x_title = TRUE,  show_legend = TRUE)
  
  legend_right <- cowplot::get_legend(
    p_subj + theme(
      legend.position   = "right",
      legend.box.margin = margin(2, 2, 2, 2),
      legend.key.size   = grid::unit(0.70, "cm"),
      legend.text       = element_text(size = 13),
      legend.title      = element_text(size = 14, face = "bold")
    )
  )
  
  p_subj_no_legend <- p_subj + theme(legend.position = "none")
  
  left_stack <- plot_grid(
    p_ff, p_taxt, p_subj_no_legend,
    labels = c("a.", "b.", "c."),
    label_size = 14, label_fontface = "bold",
    ncol = 1, align = "v", axis = "l",
    label_x = 0.07, label_y = 0.98, label_hjust = 0,
    rel_heights = c(1.35, 1.35, 1.45)
  )
  
  # Legend vertical nudge (you preferred 0.85)
  legend_shift <- 0.85
  right_column <- plot_grid(
    NULL, legend_right, NULL,
    ncol = 1,
    rel_heights = c(1.35 - legend_shift, 1.35 + legend_shift, 1.45)
  )
  
  final <- plot_grid(
    left_stack, right_column,
    ncol = 2, rel_widths = c(6.1, 0.9), align = "h"
  )
  
  ggdraw(final) + theme(plot.margin = margin(2, 2, 2, 2))
}

# ---- TIME ANOVA + Tukey letters + time boxplot (NEW) ----
time_letters_positions <- function(df_year){
  df_year %>%
    mutate(type_nice = factor(recode(tolower(type),
                                     "ff"="FF","taxt"="TA.XT","subj"="SUBJ"),
                              levels = c("FF","TA.XT","SUBJ"))) %>%
    group_by(type_nice) %>%
    summarise(y_max = max(time, na.rm = TRUE), .groups="drop") %>%
    mutate(y_lab = y_max + 0.06 * max(y_max, na.rm = TRUE))
}

plot_time_box_letters <- function(df_year, year){
  if (!nrow(df_year)) return(NULL)
  
  fit  <- aov(time ~ type, data = df_year)
  grp  <- agricolae::HSD.test(fit, "type", group = TRUE)$groups
  labs <- tibble(method = rownames(grp),
                 mean   = grp$means,
                 letters= grp$groups) %>%
    mutate(type_nice = recode(tolower(method),
                              "ff"="FF","taxt"="TA.XT","subj"="SUBJ"))
  
  ypos <- time_letters_positions(df_year)
  
  dfp <- df_year %>%
    mutate(type_nice = factor(recode(tolower(type),
                                     "ff"="FF","taxt"="TA.XT","subj"="SUBJ"),
                              levels = c("FF","TA.XT","SUBJ")))
  
  y_max  <- max(dfp$time, na.rm = TRUE)
  offset <- 0.15 * y_max
  
  ggplot(dfp, aes(x = type_nice, y = time)) +
    geom_boxplot(outlier.shape = 16, width = 0.65, size = 0.6, fatten = 2, fill = "grey90") +
    geom_text(data = dplyr::left_join(labs, ypos, by = "type_nice"),
              aes(x = type_nice, y = y_lab, label = letters),
              vjust = -0.2, size = 5, fontface = "bold") +
    labs(x = "Method", y = "Seconds per berry (s)",
         title = paste0("Measurement Time by Method (", year, ")")) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(face = "bold", size = 16),
      axis.text  = element_text(size = 14),
      plot.margin = margin(6, 8, 6, 6)
    ) +
    expand_limits(y = y_max + offset)
}

analyze_time_anova <- function(df_year, year, outdir, save_outputs = TRUE){
  if (save_outputs && !dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  fit <- aov(time ~ type, data = df_year)
  aov_tidy <- broom::tidy(fit)
  
  aov_file <- file.path(outdir, paste0("ANOVA_Time_", year, ".csv"))
  if (file.exists(aov_file)) file.remove(aov_file)
  readr::write_csv(aov_tidy, aov_file)
  
  hsd <- agricolae::HSD.test(fit, "type", group = TRUE)
  hsd_file <- file.path(outdir, paste0("ANOVA_Time_TukeyGroups_", year, ".csv"))
  if (file.exists(hsd_file)) file.remove(hsd_file)
  readr::write_csv(
    tibble::tibble(method = rownames(hsd$groups),
                   mean   = hsd$groups$means,
                   group  = hsd$groups$groups),
    hsd_file
  )
  
  p_time_box <- plot_time_box_letters(df_year, year)
  fig_file <- file.path(outdir, paste0("Fig4_TimeBox_", year, ".png"))
  if (file.exists(fig_file)) file.remove(fig_file)
  ggsave(fig_file, p_time_box, width = 7.5, height = 6.2, dpi = 600)
  
  message("Saved: ", normalizePath(aov_file))
  message("Saved: ", normalizePath(hsd_file))
  message("Saved: ", normalizePath(fig_file))
  
  print(p_time_box)
}

# ---- Year driver (saves boxpanel, time bar, correlation) ----
analyze_year <- function(df_year, year, outdir, save_outputs = TRUE){
  if (save_outputs && !dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  panel   <- assemble_year_panel(df_year, year)
  p_time  <- plot_time_bar(df_year, year)
  
  f1 <- file.path(outdir, paste0("Fig1_Boxplots_", year, ".png"))
  f3 <- file.path(outdir, paste0("Fig3_TimeBar_", year, ".png"))
  if (file.exists(f1)) file.remove(f1)
  if (file.exists(f3)) file.remove(f3)
  ggsave(filename = f1, plot = panel, width = 12.5, height = 8.8, dpi = 600)
  ggsave(filename = f3, plot = p_time, width = 10, height = 7, dpi = 600)
  message("Saved: ", normalizePath(f1))
  message("Saved: ", normalizePath(f3))
  
  # correlation panels (no title; smaller method labels)
  plot_corr_pairs(df_year, year, outdir, save_outputs,
                  method_label_cex = 1.4, r_cex = 1.6)
  
  print(panel); print(p_time)
}

# ---- RUN (both years) ----
analyze_year(dat23, 2023, outdir = out_dir, save_outputs = TRUE)
analyze_year(dat24, 2024, outdir = out_dir, save_outputs = TRUE)

# NEW: time ANOVA + Tukey + time boxplots per year
analyze_time_anova(dat23, 2023, outdir = out_dir, save_outputs = TRUE)
analyze_time_anova(dat24, 2024, outdir = out_dir, save_outputs = TRUE)

cat("\nAll files saved to:\n", normalizePath(out_dir), "\n")
print(list.files(out_dir, full.names = TRUE))
