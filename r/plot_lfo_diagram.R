# LFO Cross-Validation Diagram for Survival Data
# Visualizes LFO vs LOO CV with staggered patient entry

library(tidyverse)
library(patchwork)

# AZ Color Palette (from .Rprofile)
AZ_plum <- "#830051"
AZ_gold <- "#F0AB00"
AZ_turquoise <- "#68D2DF"
AZ_darkpurple <- "#3C1053"
AZ_green <- "#C4D600"
AZ_navy <- "#003865"
AZ_platinum <- "#9DB0AC"
AZ_darkgrey <- "#3F4444"
AZ_pink <- "#D0006F"
AZ_lightpurple <- "#AE53DE"

# Create synthetic patient data for demonstration
set.seed(123)
n_patients <- 8
study_start <- 0
study_end <- 365

# Generate patient entry times and observation periods
patients <- tibble(
  patient_id = 1:n_patients,
  entry_day = sort(runif(n_patients, 0, 200)),
  observation_end = entry_day + runif(n_patients, 100, 300)
) %>%
  mutate(observation_end = pmin(observation_end, study_end))

# Generate visit times (roughly every 4 weeks with some variation)
generate_visits <- function(entry, end, mean_interval = 28) {
  visits <- entry
  current_day <- entry
  while (current_day < end) {
    # Add some randomness: +/- 7 days around the mean interval
    next_interval <- rnorm(1, mean = mean_interval, sd = 7)
    current_day <- current_day + next_interval
    if (current_day <= end) {
      visits <- c(visits, current_day)
    }
  }
  return(visits)
}

# Create visit data for all patients
visit_data <- map_dfr(1:n_patients, function(i) {
  visits <- generate_visits(patients$entry_day[i], patients$observation_end[i])
  tibble(
    patient_id = i,
    visit_day = visits
  )
})

# Define LFO cutoffs
lfo_cutoffs <- seq(100, 300, by = 50)
lfo_window <- 60 # Days after cutoff to evaluate

# ===== PLOT 1: LFO Cross-Validation =====
create_lfo_plot <- function(cutoff_idx = 1) {
  cutoff <- lfo_cutoffs[cutoff_idx]

  # Prepare data for plotting
  plot_data <- patients %>%
    mutate(
      train_end = pmin(observation_end, cutoff),
      test_start = cutoff,
      test_end = pmin(observation_end, cutoff + lfo_window),
      has_train = entry_day < cutoff,
      has_test = observation_end > cutoff,
      train_length = ifelse(has_train, train_end - entry_day, 0),
      test_length = ifelse(has_test, test_end - test_start, 0)
    )

  p <- ggplot(plot_data) +
    # Training data (before cutoff)
    geom_segment(
      aes(x = entry_day, xend = train_end, y = patient_id, yend = patient_id),
      data = . %>% filter(has_train),
      color = AZ_navy,
      linewidth = 3,
      alpha = 0.8
    ) +
    # Test data (after cutoff, within window)
    geom_segment(
      aes(x = test_start, xend = test_end, y = patient_id, yend = patient_id),
      data = . %>% filter(has_test, test_length > 0),
      color = AZ_gold,
      linewidth = 3,
      alpha = 0.8
    ) +
    # Visit points - training period
    geom_point(
      data = visit_data %>%
        filter(visit_day < cutoff) %>%
        mutate(
          y_offset = patient_id + ifelse(row_number() %% 2 == 0, 0.15, -0.15)
        ),
      aes(x = visit_day, y = y_offset),
      color = AZ_navy,
      size = 2.5,
      alpha = 0.9,
      shape = 23,
      fill = AZ_navy
    ) +
    # Visit points - test period
    geom_point(
      data = visit_data %>%
        filter(visit_day >= cutoff, visit_day <= cutoff + lfo_window) %>%
        mutate(
          y_offset = patient_id + ifelse(row_number() %% 2 == 0, 0.15, -0.15)
        ),
      aes(x = visit_day, y = y_offset),
      color = AZ_gold,
      size = 2.5,
      alpha = 1,
      shape = 23,
      fill = AZ_gold
    ) +
    # Patient entry points
    geom_point(
      aes(x = entry_day, y = patient_id),
      color = AZ_plum,
      size = 3,
      shape = 16
    ) +
    # Cutoff line
    geom_vline(
      xintercept = cutoff,
      linetype = "dashed",
      color = AZ_darkgrey,
      linewidth = 1
    ) +
    # Test window end
    geom_vline(
      xintercept = cutoff + lfo_window,
      linetype = "dotted",
      color = AZ_darkgrey,
      linewidth = 0.8
    ) +
    # Annotations
    annotate(
      "text",
      x = cutoff,
      y = n_patients + 0.5,
      label = "Cutoff",
      color = AZ_darkgrey,
      hjust = -0.1,
      fontface = "bold"
    ) +
    annotate(
      "rect",
      xmin = cutoff,
      xmax = cutoff + lfo_window,
      ymin = 0,
      ymax = n_patients + 1,
      fill = AZ_gold,
      alpha = 0.1
    ) +
    annotate(
      "text",
      x = cutoff + lfo_window / 2,
      y = 0.3,
      label = "Test Window",
      color = AZ_gold,
      fontface = "bold"
    ) +
    scale_y_continuous(
      breaks = 1:n_patients,
      labels = paste("Patient", 1:n_patients)
    ) +
    scale_x_continuous(breaks = seq(0, study_end, by = 50)) +
    labs(
      title = sprintf("LFO Cross-Validation: Cutoff %d", cutoff_idx),
      subtitle = sprintf(
        "Train on data before day %d | Test on next %d days",
        cutoff,
        lfo_window
      ),
      x = "Study Day",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", color = AZ_navy),
      plot.subtitle = element_text(color = AZ_darkgrey),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(color = AZ_darkgrey, hjust = 1)
    )

  return(p)
}

# ===== PLOT 2: All LFO Cutoffs Overview =====
create_lfo_overview <- function() {
  # Create data for all cutoffs
  all_cutoffs <- map_dfr(seq_along(lfo_cutoffs), function(i) {
    cutoff <- lfo_cutoffs[i]
    patients %>%
      mutate(
        cutoff_id = i,
        cutoff_day = cutoff,
        train_end = pmin(observation_end, cutoff),
        test_start = cutoff,
        test_end = pmin(observation_end, cutoff + lfo_window),
        has_train = entry_day < cutoff,
        has_test = observation_end > cutoff & observation_end > cutoff,
        phase = "both"
      )
  })

  p <- ggplot() +
    # Patient timelines
    geom_segment(
      data = patients,
      aes(
        x = entry_day,
        xend = observation_end,
        y = patient_id,
        yend = patient_id
      ),
      color = AZ_platinum,
      linewidth = 2,
      alpha = 0.6
    ) +
    # Visit points
    geom_point(
      data = visit_data %>%
        group_by(patient_id) %>%
        mutate(
          y_offset = patient_id + ifelse(row_number() %% 2 == 0, 0.15, -0.15)
        ) %>%
        ungroup(),
      aes(x = visit_day, y = y_offset),
      color = AZ_darkgrey,
      size = 2,
      alpha = 0.8,
      shape = 23,
      fill = AZ_darkgrey
    ) +
    # Cutoff lines
    geom_vline(
      xintercept = lfo_cutoffs,
      linetype = "dashed",
      color = AZ_navy,
      linewidth = 0.7,
      alpha = 0.7
    ) +
    # Test windows - use single color for clarity
    geom_rect(
      data = tibble(cutoff = lfo_cutoffs, cutoff_id = seq_along(lfo_cutoffs)),
      aes(
        xmin = cutoff,
        xmax = cutoff + lfo_window,
        ymin = 0.3,
        ymax = n_patients + 0.7
      ),
      fill = AZ_gold,
      alpha = 0.12
    ) +
    # Patient entry points
    geom_point(
      data = patients,
      aes(x = entry_day, y = patient_id),
      color = AZ_plum,
      size = 2.5,
      shape = 16
    ) +
    scale_y_continuous(
      breaks = 1:n_patients,
      labels = paste("P", 1:n_patients)
    ) +
    scale_x_continuous(breaks = seq(0, study_end, by = 50)) +
    labs(
      title = "LFO Cross-Validation: Sequential Forward Validation",
      subtitle = "Multiple cutoffs with fixed-width test windows",
      x = "Study Day",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", color = AZ_navy),
      plot.subtitle = element_text(color = AZ_darkgrey),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(color = AZ_darkgrey),
      legend.position = "none"
    )

  return(p)
}

# ===== PLOT 3: LOO Cross-Validation =====
create_loo_plot <- function(held_out_patient = 1) {
  plot_data <- patients %>%
    mutate(
      is_held_out = patient_id == held_out_patient,
      color_group = ifelse(is_held_out, "held_out", "training")
    )

  p <- ggplot(plot_data) +
    # Training patients
    geom_segment(
      aes(
        x = entry_day,
        xend = observation_end,
        y = patient_id,
        yend = patient_id,
        color = color_group,
        linewidth = is_held_out
      ),
      alpha = 0.85
    ) +
    # Visit points
    geom_point(
      data = visit_data %>%
        left_join(
          plot_data %>% select(patient_id, color_group),
          by = "patient_id"
        ) %>%
        group_by(patient_id) %>%
        mutate(
          y_offset = patient_id + ifelse(row_number() %% 2 == 0, 0.15, -0.15)
        ) %>%
        ungroup(),
      aes(x = visit_day, y = y_offset, color = color_group, fill = color_group),
      size = 2,
      alpha = 0.9,
      shape = 23
    ) +
    # Patient entry points
    geom_point(
      aes(
        x = entry_day,
        y = patient_id,
        color = color_group,
        size = is_held_out
      ),
      shape = 16
    ) +
    scale_color_manual(values = c("training" = AZ_navy, "held_out" = AZ_pink)) +
    scale_linewidth_manual(values = c("FALSE" = 2.5, "TRUE" = 3.5)) +
    scale_size_manual(values = c("FALSE" = 2.5, "TRUE" = 4)) +
    scale_y_continuous(
      breaks = 1:n_patients,
      labels = paste("Patient", 1:n_patients)
    ) +
    scale_x_continuous(breaks = seq(0, study_end, by = 50)) +
    labs(
      title = sprintf(
        "LOO Cross-Validation: Hold Out Patient %d",
        held_out_patient
      ),
      subtitle = "Train on all other patients | Test on entire held-out patient timeline",
      x = "Study Day",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", color = AZ_navy),
      plot.subtitle = element_text(color = AZ_darkgrey),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(color = AZ_darkgrey, hjust = 1),
      legend.position = "none"
    )

  return(p)
}

# ===== PLOT 4: Side-by-side Comparison =====
create_comparison_plot <- function() {
  # LFO panel
  lfo_data <- map_dfr(seq_along(lfo_cutoffs), function(i) {
    cutoff <- lfo_cutoffs[i]
    patients %>%
      mutate(
        cutoff_id = i,
        cutoff_day = cutoff,
        method = "LFO",
        fold = i
      )
  })

  # LOO panel
  loo_data <- map_dfr(1:n_patients, function(i) {
    patients %>%
      mutate(
        held_out = patient_id == i,
        method = "LOO",
        fold = i
      )
  })

  # LFO plot
  p_lfo <- ggplot() +
    geom_segment(
      data = patients,
      aes(
        x = entry_day,
        xend = observation_end,
        y = patient_id,
        yend = patient_id
      ),
      color = AZ_platinum,
      linewidth = 2,
      alpha = 0.4
    ) +
    # Visit points
    geom_point(
      data = visit_data %>%
        group_by(patient_id) %>%
        mutate(
          y_offset = patient_id + ifelse(row_number() %% 2 == 0, 0.15, -0.15)
        ) %>%
        ungroup(),
      aes(x = visit_day, y = y_offset),
      color = AZ_darkgrey,
      size = 1.5,
      alpha = 0.7,
      shape = 23,
      fill = AZ_darkgrey
    ) +
    geom_vline(
      data = tibble(cutoff = lfo_cutoffs),
      aes(xintercept = cutoff),
      linetype = "dashed",
      color = AZ_navy,
      linewidth = 0.6,
      alpha = 0.7
    ) +
    geom_rect(
      data = tibble(cutoff = lfo_cutoffs, id = seq_along(lfo_cutoffs)),
      aes(
        xmin = cutoff,
        xmax = cutoff + lfo_window,
        ymin = 0.3,
        ymax = n_patients + 0.7
      ),
      fill = AZ_gold,
      alpha = 0.15
    ) +
    geom_point(
      data = patients,
      aes(x = entry_day, y = patient_id),
      color = AZ_plum,
      size = 2
    ) +
    scale_y_continuous(
      breaks = 1:n_patients,
      labels = paste("P", 1:n_patients)
    ) +
    labs(
      title = "LFO: Time-Based Folds",
      subtitle = "Train → Cutoff → Test Window",
      x = "Study Day",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", color = AZ_navy),
      plot.subtitle = element_text(color = AZ_darkgrey),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    )

  # LOO plot - show concept with rectangles for each fold
  loo_concept <- tibble(
    fold = 1:n_patients,
    patient_id = 1:n_patients
  )

  p_loo <- ggplot(patients) +
    geom_segment(
      aes(
        x = entry_day,
        xend = observation_end,
        y = patient_id,
        yend = patient_id
      ),
      color = AZ_navy,
      linewidth = 2,
      alpha = 0.4
    ) +
    # Visit points
    geom_point(
      data = visit_data %>%
        group_by(patient_id) %>%
        mutate(
          y_offset = patient_id + ifelse(row_number() %% 2 == 0, 0.15, -0.15)
        ) %>%
        ungroup(),
      aes(x = visit_day, y = y_offset),
      color = AZ_navy,
      size = 1.5,
      alpha = 0.7,
      shape = 23,
      fill = AZ_navy
    ) +
    geom_segment(
      data = \(d) {
        d %>% mutate(highlight = patient_id %% 2 == 1) %>% filter(highlight)
      },
      aes(
        x = entry_day,
        xend = observation_end,
        y = patient_id,
        yend = patient_id
      ),
      color = AZ_pink,
      linewidth = 2,
      alpha = 0.8
    ) +
    # Visit points for highlighted patients
    geom_point(
      data = visit_data %>%
        filter(patient_id %% 2 == 1) %>%
        group_by(patient_id) %>%
        mutate(
          y_offset = patient_id + ifelse(row_number() %% 2 == 0, 0.15, -0.15)
        ) %>%
        ungroup(),
      aes(x = visit_day, y = y_offset),
      color = AZ_pink,
      size = 1.5,
      alpha = 0.9,
      shape = 23,
      fill = AZ_pink
    ) +
    geom_point(aes(x = entry_day, y = patient_id), color = AZ_plum, size = 2) +
    annotate(
      "text",
      x = study_end - 50,
      y = n_patients / 2,
      label = "Each patient\nheld out once",
      color = AZ_pink,
      fontface = "bold",
      hjust = 0.5
    ) +
    scale_y_continuous(
      breaks = 1:n_patients,
      labels = paste("P", 1:n_patients)
    ) +
    labs(
      title = "LOO: Patient-Based Folds",
      subtitle = "N folds, 1 patient held out each time",
      x = "Study Day",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", color = AZ_navy),
      plot.subtitle = element_text(color = AZ_darkgrey),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    )

  # Combine
  p_combined <- p_lfo +
    p_loo +
    plot_annotation(
      title = "Cross-Validation Strategies for Survival Data",
      subtitle = "LFO evaluates temporal predictive performance | LOO evaluates patient-level generalization",
      theme = theme(
        plot.title = element_text(face = "bold", color = AZ_navy),
        plot.subtitle = element_text(color = AZ_darkgrey)
      )
    )

  return(p_combined)
}

# ===== Generate Example Plots =====

# Individual LFO cutoff examples
plot_lfo_cutoff_1 <- create_lfo_plot(1)
plot_lfo_cutoff_2 <- create_lfo_plot(2)
plot_lfo_cutoff_3 <- create_lfo_plot(3)

# LFO overview with all cutoffs
plot_lfo_overview <- create_lfo_overview()

# LOO examples
plot_loo_patient_1 <- create_loo_plot(1)
plot_loo_patient_3 <- create_loo_plot(3)

# Comparison plot
plot_comparison <- create_comparison_plot()

# Display main comparison
print(plot_comparison)

# Optional: Save plots
# ggsave("lfo_vs_loo_comparison.png", plot_comparison,
#        width = 12, height = 6, dpi = 300, bg = "white")
# ggsave("lfo_overview.png", plot_lfo_overview,
#        width = 10, height = 6, dpi = 300, bg = "white")
