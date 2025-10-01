# ============================================================
# Blackberry firmness — two-year pipeline (publication style)
# Shows plots in RStudio "Plots" pane by default.
# Set save_outputs <- TRUE to also write PNG/CSV files.
# ============================================================

rm(list = ls())

# --- Packages ---
pkgs <- c("dplyr","tidyr","ggplot2","cowplot","agricolae","stringr","readr","tibble")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# --- Toggle: show vs save ---
save_outputs <- FALSE  # TRUE to also save PNG/CSV outputs

# --- File paths (EDIT if needed) ---
file_2023 <- "C:/Users/RhysB/OneDrive/Desktop/2023FirmnessData.csv"
file_2024 <- "C:/Users/RhysB/OneDrive/Desktop/2024FirmnessData.csv"

# --- Helpers ---
`%||%` <- function(a,b) if (is.null(a)) b else a

read_firm <- function(path, year){
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  keep <- intersect(names(df), c("seln","harv_no","rep","hdate","evaldate","type","time","force","meansep"))
  df <- df[, keep, drop = FALSE]
  df %>%
    mutate(
      seln     = as.factor(seln),
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

# HSD letters for force ~ seln within a method
hsd_letters_force <- function(d) {
  if (nrow(d) == 0) return(NULL)
  fit <- aov(force ~ seln, data = d)
  out <- agricolae::HSD.test(fit, "seln", group = TRUE)$groups
  tibble(seln = rownames(out),
         mean = out$means %||% out$`means`,
         letters = out$groups) %>% arrange(seln)
}

# Y positions for letters above boxes
box_y_positions <- function(d) {
  d %>%
    group_by(seln) %>%
    summarise(y_max = max(force, na.rm = TRUE), .groups = "drop") %>%
    mutate(y_lab = y_max + 0.06 * (max(y_max, na.rm = TRUE)))
}

# Mean-by-(seln x type) wide
mean_wide <- function(df) {
  df %>%
    group_by(seln, type) %>%
    summarise(force = mean(force, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = type, values_from = force, names_prefix = "force.") %>%
    arrange(seln)
}

# --- Publication themes ---
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

# --- CORRELATION PANELS ---
plot_corr_pairs <- function(df_year, year, outdir, save_outputs){
  mw <- mean_wide(df_year)
  cols <- c("force.ff","force.subj","force.taxt")
  cols <- cols[cols %in% names(mw)]
  if (length(cols) < 2) return(invisible(NULL))
  M <- dplyr::select(mw, all_of(cols))
  
  diag_labels <- c("FruitFirm 1000", "Subjective", "TA.XT")
  
  cor_mat <- suppressWarnings(cor(M, use = "pairwise.complete.obs", method = "pearson"))
  if (save_outputs) {
    write.csv(cor_mat,
              file = file.path(outdir, paste0("Fig2_Correlation_", year, "_correlations.csv")),
              row.names = TRUE)
  }
  
  panel_cor <- function(x, y, digits = 2){
    usr <- par("usr"); on.exit(par(usr = usr))
    par(usr = c(0,1,0,1))
    r <- suppressWarnings(cor(x, y, use = "pairwise.complete.obs"))
    text(0.5, 0.5, paste0("r = ", format(round(r, digits), nsmall = digits)),
         cex = 1.9, font = 2)
  }
  panel_pts <- function(x, y) points(x, y, pch = 19, cex = 1)
  
  main_title <- paste0("Correlation of Mean Firmness by Methods (", year, ")")
  
  if (save_outputs) {
    png(file.path(outdir, paste0("Fig2_Correlation_", year, ".png")),
        width = 1600, height = 1600, res = 200)
    pairs(M, labels = diag_labels, lower.panel = panel_cor, upper.panel = panel_pts,
          main = main_title, cex.labels = 1.9, font.labels = 2)
    dev.off()
  }
  pairs(M, labels = diag_labels, lower.panel = panel_cor, upper.panel = panel_pts,
        main = main_title, cex.labels = 1.9, font.labels = 2)
}

# --- BOX-AND-WHISKER (now with aspect.ratio to make them taller) ---
plot_box_with_letters <- function(d, title, ylab){
  if (nrow(d) == 0) return(NULL)
  
  letters <- hsd_letters_force(d)
  ypos    <- box_y_positions(d)
  lab_df  <- dplyr::left_join(letters, ypos, by = "seln")
  
  # Headroom so mean-sep letters never clip
  y_max  <- max(d$force, na.rm = TRUE)
  offset <- 0.15 * y_max
  
  ggplot(d, aes(x = seln, y = force)) +
    geom_boxplot(outlier.shape = 16, width = 0.72, size = 0.6, fatten = 2) +
    geom_text(data = lab_df,
              aes(x = seln, y = y_lab, label = letters),
              vjust = -0.15, size = 4.5) +
    labs(x = "Genotype", y = ylab, title = title) +
    theme_gray_pub(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    expand_limits(y = y_max + offset)
}


# --- TIME BAR CHART ---
plot_time_bar <- function(df_year, year){
  stats <- df_year %>%
    dplyr::group_by(type) %>%
    dplyr::summarise(
      n        = dplyr::n(),
      avg_time = mean(time, na.rm = TRUE),
      sd_time  = sd(time,   na.rm = TRUE),
      se_time  = sd_time / sqrt(n),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      type_nice = factor(dplyr::recode(tolower(type),
                                       "ff"   = "FF",
                                       "taxt" = "TA.XT",
                                       "subj" = "SUBJ"),
                         levels = c("FF","TA.XT","SUBJ"))
    )
  
  y_max  <- max(stats$avg_time + stats$se_time, na.rm = TRUE)
  pad    <- 0.15 * y_max
  
  ggplot(stats, aes(x = type_nice, y = avg_time)) +
    geom_col(width = 0.6, fill = "#A6CEE3", colour = "black", linewidth = 0.4) +
    geom_errorbar(aes(ymin = avg_time - se_time, ymax = avg_time + se_time),
                  width = 0.3, linewidth = 1.4, colour = "black") +
    coord_flip() +
    labs(x = "Method", y = "Seconds per berry (s)",
         title = paste0("Efficiency by Method — ", year)) +
    theme_classic_pub(base_size = 14) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title   = element_text(face = "bold", size = 14),
      axis.text    = element_text(size = 12),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA)
    ) +
    expand_limits(y = y_max + pad)
}

# --- Yearly driver ---
analyze_year <- function(df_year, year, outdir = NULL, save_outputs = FALSE){
  outdir <- outdir %||% file.path(getwd(), paste0("firmness_outputs_", year))
  if (save_outputs && !dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  ff   <- df_year %>% filter(type == "ff")
  taxt <- df_year %>% filter(type == "taxt")
  subj <- df_year %>% filter(type == "subj")
  
  if (save_outputs) {
    ff_letters   <- hsd_letters_force(ff);   if(!is.null(ff_letters))   write.csv(ff_letters,   file.path(outdir, paste0("HSD_letters_FF_",   year, ".csv")), row.names = FALSE)
    taxt_letters <- hsd_letters_force(taxt); if(!is.null(taxt_letters)) write.csv(taxt_letters, file.path(outdir, paste0("HSD_letters_TAXT_", year, ".csv")), row.names = FALSE)
    subj_letters <- hsd_letters_force(subj); if(!is.null(subj_letters)) write.csv(subj_letters, file.path(outdir, paste0("HSD_letters_SUBJ_", year, ".csv")), row.names = FALSE)
  } else {
    invisible(hsd_letters_force(ff)); invisible(hsd_letters_force(taxt)); invisible(hsd_letters_force(subj))
  }
  
  p_ff   <- plot_box_with_letters(ff,   paste0("FruitFirm 1000 — ", year), "Force (g/mm)")
  p_taxt <- plot_box_with_letters(taxt, paste0("TA.XT — ", year),          "Force (g)")
  p_subj <- plot_box_with_letters(subj, paste0("Subjective — ", year),     "Firmness (1-5)")
  
  panel <- plot_grid(
    p_ff, p_taxt, p_subj,
    labels = c("a.", "b.", "c."),
    label_size = 14,
    label_fontface = "bold",
    ncol = 1,
    align = "v",
    label_x = 0.07,
    label_y = 0.98,
    label_hjust = 0
  )
  
  if (save_outputs) {
    ggsave(filename = file.path(outdir, paste0("Fig1_Boxplots_", year, ".png")),
           plot = panel, width = 10, height = 12, dpi = 300)
  }
  print(panel)
  
  p_time <- plot_time_bar(df_year, year)
  if (save_outputs) {
    ggsave(filename = file.path(outdir, paste0("Fig3_TimeBar_", year, ".png")),
           plot = p_time, width = 7, height = 6, dpi = 300)
  }
  print(p_time)
  
  plot_corr_pairs(df_year, year, outdir, save_outputs)
}

# --- Run ---
analyze_year(dat23, 2023, save_outputs = save_outputs)
analyze_year(dat24, 2024, save_outputs = save_outputs)
