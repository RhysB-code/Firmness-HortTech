# ============================================================
# Blackberry firmness & timing — FINAL (pub style, robust)
# - Boxplots by genotype (FF / TA.XT / SUBJ) with Tukey letters (force)
# - Time efficiency bars (seconds/berry) with SE + Tukey letters (robust)
# - One-way TIME ANOVA (time ~ method) per year + assumptions + eta^2
# - High-res PNGs saved to Desktop/Blackberry_Firmness_Exports
# ============================================================

rm(list = ls()); options(stringsAsFactors = FALSE)

# ---- Packages ----
pkgs <- c(
  "dplyr","tidyr","ggplot2","cowplot","agricolae","stringr","readr",
  "tibble","grid","RColorBrewer","scales","broom","car","effectsize","rlang"
)
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# Prefer dplyr::recode explicitly (avoid car::recode collisions)
recode <- dplyr::recode

# ---- File paths (EDIT if needed) ----
file_2023 <- "C:/Users/RhysB/OneDrive/Desktop/2023FirmnessData.csv"
file_2024 <- "C:/Users/RhysB/OneDrive/Desktop/2024FirmnessData.csv"

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

# Safe printer to avoid "numeric -> grob" warnings
safe_print <- function(p){
  if (inherits(p, c("gg","ggplot"))) print(p)
  invisible(NULL)
}

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
  tibble(
    seln = rownames(out),
    # agricolae sometimes stores 'means' or `means`; either way, %||% covers it
    mean = out$means %||% out$`means`,
    letters = out$groups
  ) %>% arrange(seln)
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

# ---- Genotype boxplots (force) with letters ----
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

# ---- Build Tukey labels for time (robust to agricolae variants) ----
build_time_tukey_labels <- function(df_year){
  df_year <- df_year %>% mutate(type = factor(tolower(type), levels = c("ff","taxt","subj")))
  
  # Always compute means as fallback
  stats_means <- df_year %>%
    group_by(type) %>%
    summarise(avg_time = mean(time, na.rm = TRUE), .groups = "drop")
  
  tukey_groups <- tryCatch({
    fit <- aov(time ~ type, data = df_year)
    g   <- agricolae::HSD.test(fit, "type", group = TRUE)$groups
    if (is.null(g) || nrow(g) == 0) stop("Empty HSD groups")
    
    # Agricolae may name the mean column after the response (e.g., 'time') or as 'means'/'mean'
    mean_col <- if ("time"  %in% names(g)) "time"  else
      if ("means" %in% names(g)) "means" else
        if ("mean"  %in% names(g)) "mean"  else NA_character_
    if (is.na(mean_col)) stop("No mean column found in HSD groups")
    
    g_df <- tibble::tibble(
      type_raw     = rownames(g),
      avg_from_hsd = as.numeric(g[[mean_col]]),
      letters      = g[["groups"]]
    )
    
    stats_means %>%
      rename(type_raw = type) %>%
      left_join(g_df, by = "type_raw") %>%
      transmute(
        type      = type_raw,
        avg_time  = ifelse(is.na(avg_from_hsd), avg_time, avg_from_hsd),
        letters   = letters %||% ""   # empty if missing
      )
  }, error = function(e){
    # Fallback: no letters, keep means
    stats_means %>% transmute(type = type, avg_time = avg_time, letters = "")
  })
  
  tukey_groups %>%
    mutate(type_nice = factor(dplyr::recode(as.character(type),
                                            "ff"="FF","taxt"="TA.XT","subj"="SUBJ"),
                              levels = c("FF","TA.XT","SUBJ")))
}

# ---- Efficiency bars (time) with mean labels + Tukey letters (robust) ----
plot_time_bar_with_letters <- function(df_year, year){
  stats <- df_year %>%
    mutate(type = factor(tolower(type), levels = c("ff","taxt","subj"))) %>%
    group_by(type) %>%
    summarise(
      n        = n(),
      avg_time = mean(time, na.rm = TRUE),
      sd_time  = sd(time,   na.rm = TRUE),
      se_time  = ifelse(is.na(sd_time) | n <= 1, NA_real_, sd_time / sqrt(n)),
      .groups  = "drop"
    ) %>%
    mutate(type_nice = factor(dplyr::recode(as.character(type),
                                            "ff"="FF","taxt"="TA.XT","subj"="SUBJ"),
                              levels = c("FF","TA.XT","SUBJ")))
  
  labs  <- build_time_tukey_labels(df_year)
  
  y_max <- max(stats$avg_time + ifelse(is.na(stats$se_time), 0, stats$se_time), na.rm = TRUE)
  if (!is.finite(y_max)) y_max <- max(stats$avg_time, na.rm = TRUE)
  y_max <- y_max * 1.15
  
  p <- ggplot(stats, aes(x = type_nice, y = avg_time)) +
    geom_col(width = 0.6, fill = "#BFD7EA", colour = "black", linewidth = 0.5) +
    geom_errorbar(aes(ymin = avg_time - se_time, ymax = avg_time + se_time),
                  width = 0.28, linewidth = 1.2, colour = "black") +
    geom_text(aes(label = sprintf("%.1f s", avg_time),
                  y = avg_time + ifelse(is.na(se_time), 0, se_time) + 0.5),
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
  
  # Add letters only if present and non-empty
  if (!is.null(labs) && nrow(labs)){
    labs_nonempty <- labs %>%
      filter(!is.na(letters) & letters != "") %>%
      select(type_nice, letters) %>%
      distinct()
    if (nrow(labs_nonempty)){
      p <- p + geom_text(data = labs_nonempty,
                         aes(x = type_nice, y = y_max*0.98, label = letters),
                         size = 5.2, fontface = "bold")
    }
  }
  p
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
  
  left_stack <- cowplot::plot_grid(
    p_ff, p_taxt, p_subj_no_legend,
    labels = c("a.", "b.", "c."),
    label_size = 14, label_fontface = "bold",
    ncol = 1, align = "v", axis = "l",
    label_x = 0.07, label_y = 0.98, label_hjust = 0,
    rel_heights = c(1.35, 1.35, 1.45)
  )
  
  # Legend vertical nudge
  legend_shift <- 0.85
  right_column <- cowplot::plot_grid(
    NULL, legend_right, NULL,
    ncol = 1,
    rel_heights = c(1.35 - legend_shift, 1.35 + legend_shift, 1.45)
  )
  
  final <- cowplot::plot_grid(
    left_stack, right_column,
    ncol = 2, rel_widths = c(6.1, 0.9), align = "h"
  )
  
  cowplot::ggdraw(final) + theme(plot.margin = margin(2, 2, 2, 2))
}

# ---- TIME ANOVA (simple) + assumptions + effect size + Tukey letters (robust export) ----
time_anova_simple <- function(df_year, year, outdir, save_outputs = TRUE){
  if (!nrow(df_year)) return(invisible(NULL))
  if (save_outputs && !dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  df_year <- df_year |> mutate(type = factor(tolower(type), levels = c("ff","taxt","subj")))
  fit <- aov(time ~ type, data = df_year)
  
  # Tidy ANOVA
  aov_tidy <- broom::tidy(fit)
  
  # Assumption checks
  res <- residuals(fit)
  shapiro_p <- tryCatch(shapiro.test(res)$p.value, error = function(e) NA_real_)
  lev_p     <- tryCatch(car::leveneTest(time ~ type, data = df_year)[["Pr(>F)"]][1], error = function(e) NA_real_)
  
  # Effect size (eta^2)
  eta <- tryCatch(effectsize::eta_squared(fit, partial = FALSE), error = function(e) NULL)
  
  # Tukey group letters (robust to column naming)
  tukey_tbl <- tryCatch({
    g <- agricolae::HSD.test(fit, "type", group = TRUE)$groups
    if (is.null(g) || nrow(g) == 0) stop("Empty HSD groups")
    mean_col <- if ("time"  %in% names(g)) "time"  else
      if ("means" %in% names(g)) "means" else
        if ("mean"  %in% names(g)) "mean"  else NA_character_
    if (is.na(mean_col)) stop("No mean column in Tukey groups")
    tibble::tibble(
      method = rownames(g),
      mean   = as.numeric(g[[mean_col]]),
      group  = g[["groups"]]
    ) %>%
      mutate(method = dplyr::recode(tolower(method), ff = "FF", taxt = "TA.XT", subj = "SUBJ"))
  }, error = function(e){
    # Fallback: means by method without Tukey groups
    df_year %>%
      group_by(type) %>%
      summarise(mean = mean(time, na.rm = TRUE), .groups = "drop") %>%
      transmute(method = dplyr::recode(as.character(type),
                                       "ff"="FF","taxt"="TA.XT","subj"="SUBJ"),
                mean = mean,
                group = NA_character_)
  })
  
  # Save CSVs
  if (save_outputs){
    readr::write_csv(aov_tidy,   file.path(outdir, paste0("ANOVA_Time_", year, ".csv")))
    readr::write_csv(tukey_tbl,  file.path(outdir, paste0("ANOVA_Time_TukeyGroups_", year, ".csv")))
    readr::write_csv(
      tibble::tibble(
        shapiro_p = shapiro_p,
        levene_p  = lev_p,
        eta2      = if (!is.null(eta)) eta$Eta2[eta$Parameter == "type"] else NA_real_
      ),
      file.path(outdir, paste0("ANOVA_Time_Assumptions_", year, ".csv"))
    )
  }
  
  # Console summary (copy-paste friendly)
  cat("\n--- Time ANOVA (", year, ") ---\n", sep = "")
  print(aov_tidy)
  cat("\nAssumptions: Shapiro p =", sprintf("%.4f", shapiro_p),
      " | Levene p =", sprintf("%.4f", lev_p), "\n")
  if (!is.null(eta)){
    cat("Effect size (eta^2 for type):", round(eta$Eta2[eta$Parameter == "type"], 3), "\n")
  }
  cat("\nTukey (method means with groups):\n")
  print(tukey_tbl)
}

# ---- Year driver (saves boxpanel, time bar) ----
analyze_year <- function(df_year, year, outdir, save_outputs = TRUE){
  if (save_outputs && !dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  panel   <- assemble_year_panel(df_year, year)
  p_time  <- plot_time_bar_with_letters(df_year, year)
  
  f1 <- file.path(outdir, paste0("Fig1_Boxplots_", year, ".png"))
  f3 <- file.path(outdir, paste0("Fig3_TimeBar_", year, ".png"))
  if (file.exists(f1)) file.remove(f1)
  if (file.exists(f3)) file.remove(f3)
  ggsave(filename = f1, plot = panel, width = 12.5, height = 8.8, dpi = 600)
  ggsave(filename = f3, plot = p_time, width = 10, height = 7, dpi = 600)
  message("Saved: ", normalizePath(f1))
  message("Saved: ", normalizePath(f3))
  
  safe_print(panel)
  safe_print(p_time)
}

# ---- RUN (both years) ----
analyze_year(dat23, 2023, outdir = out_dir, save_outputs = TRUE)
analyze_year(dat24, 2024, outdir = out_dir, save_outputs = TRUE)

# SIMPLE time ANOVA + Tukey + assumptions per year
time_anova_simple(dat23, 2023, outdir = out_dir, save_outputs = TRUE)
time_anova_simple(dat24, 2024, outdir = out_dir, save_outputs = TRUE)

cat("\nAll files saved to:\n", normalizePath(out_dir), "\n")
print(list.files(out_dir, full.names = TRUE))
