####Function to create a rainbow plot for the proportion of different level of CD4 count at diagnosis####
create_cd4_rainbow_plot <- function(data, 
                                    start_year = NULL,
                                    end_year = NULL, 
                                    show_n = TRUE,          # Show sample size on x-axis
                                    # Updated filtering options
                                    previ_diag_overseas = NULL,
                                    exposure_filter = NULL,       # hiv_exposure_risk
                                    gender_filter = NULL,         # gender
                                    region_birth_filter = NULL,   # region_birth
                                    lang_home_filter = NULL,      # lang_home
                                    place_acq_filter = NULL,      # place_acq
                                    state_res_filter = NULL,      # state_res
                                    phn_res_filter = NULL,        # phn_res
                                    cd4_var = "cd4_count_diag",
                                    year_var = "year_diag",
                                    title = "CD4 Count at HIV Diagnosis by Year",
                                    subtitle = NULL,
                                    # Updated colors to accommodate the extra CD4 category
                                    colors = c("#CC79A7","#0072B2","#009E73","#F0E442","#E69F00","#D55E00","#D55E00"),
                                    legend_position = "bottom",
                                    dash_line = FALSE) {          # NEW: Option to add median CD4 line
  
  # Load required libraries
  require(tidyverse)
  require(lubridate)
  
  # Clone the data to avoid modifying the original
  plot_data <- data %>% as.data.frame()
  
  plot_data <- plot_data %>%
    filter(
      !is.na(!!sym(cd4_var)),  # Ensure CD4 count is not NA
      !is.na(!!sym(year_var))  # Ensure year is not NA
    )
  
  # Filter by year if specified
  if (!is.null(start_year)) {
    plot_data <- plot_data %>% 
      filter(!!sym(year_var) >= start_year)
  }
  
  if (!is.null(end_year)) {
    plot_data <- plot_data %>% 
      filter(!!sym(year_var) <= end_year)
  }
  # Apply filter for previously diagnosed overseas if specified
  if (!is.null(previ_diag_overseas)) {
    if (previ_diag_overseas == TRUE) {
      plot_data <- plot_data %>%
        filter(previ_diag_overseas == 1)
    } else if (previ_diag_overseas == FALSE) {
      plot_data <- plot_data %>%
        filter(previ_diag_overseas == 0)
    }
  }
  
  # Apply filters based on various demographic variables
  # HIV exposure risk
  if (!is.null(exposure_filter)) {
    if ("hiv_exposure_risk" %in% names(plot_data)) {
      plot_data <- plot_data %>% 
        filter(hiv_exposure_risk == exposure_filter)
    } else {
      warning("hiv_exposure_risk variable not found in dataset")
    }
  }
  
  # Gender
  if (!is.null(gender_filter)) {
    if ("gender" %in% names(plot_data)) {
      plot_data <- plot_data %>% 
        filter(gender == gender_filter)
    } else {
      warning("gender variable not found in dataset")
    }
  }
  
  # Region of birth
  if (!is.null(region_birth_filter)) {
    if ("region_birth" %in% names(plot_data)) {
      plot_data <- plot_data %>% 
        filter(region_birth == region_birth_filter)
    } else {
      warning("region_birth variable not found in dataset")
    }
  }
  
  # Language spoken at home
  if (!is.null(lang_home_filter)) {
    if ("lang_home" %in% names(plot_data)) {
      plot_data <- plot_data %>% 
        filter(lang_home == lang_home_filter)
    } else {
      warning("lang_home variable not found in dataset")
    }
  }
  
  # Place of acquisition
  if (!is.null(place_acq_filter)) {
    if ("place_acq" %in% names(plot_data)) {
      plot_data <- plot_data %>% 
        filter(place_acq == place_acq_filter)
    } else {
      warning("place_acq variable not found in dataset")
    }
  }
  
  # State of residence
  if (!is.null(state_res_filter)) {
    if ("state_res" %in% names(plot_data)) {
      plot_data <- plot_data %>% 
        filter(state_res == state_res_filter)
    } else {
      warning("state_res variable not found in dataset")
    }
  }
  
  # PHN of residence
  if (!is.null(phn_res_filter)) {
    if ("phn_res" %in% names(plot_data)) {
      plot_data <- plot_data %>% 
        filter(phn_res == phn_res_filter)
    } else {
      warning("phn_res variable not found in dataset")
    }
  }
  
  # Check if we have data left after all the filtering
  if (nrow(plot_data) == 0) {
    stop("No data remaining after applying all filters. Please adjust your filter criteria.")
  }
  
  # Calculate median CD4 count by year (if dash_line is TRUE)
  if (dash_line) {
    median_cd4_by_year <- plot_data %>%
      group_by(!!sym(year_var)) %>%
      summarise(
        median_cd4 = median(!!sym(cd4_var), na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  # Create CD4 categories - MODIFIED to split <200 into <50 and 50-199
  plot_data <- plot_data %>%
    mutate(cd4_category = case_when(
      newly_acq_hiv_flag == 1 & !!sym(cd4_var) < 350 ~ "Seroconversion",
      !!sym(cd4_var) < 50 ~ "< 50",           # New category
      !!sym(cd4_var) >= 50 & !!sym(cd4_var) < 200 ~ "50-199", # New category
      !!sym(cd4_var) >= 200 & !!sym(cd4_var) < 350 ~ "200-349",
      !!sym(cd4_var) >= 350 & !!sym(cd4_var) < 500 ~ "350-499",
      !!sym(cd4_var) >= 500 ~ "≥ 500",
      TRUE ~ "Unknown"
    ))
  
  # Set the order of CD4 categories - UPDATED for new categories
  plot_data$cd4_category <- factor(
    plot_data$cd4_category,
    levels = c("Unknown","Seroconversion", "≥ 500", "350-499","200-349", "50-199","< 50")
  )
  
  # Create summary data for plotting
  summary_data <- plot_data %>%
    group_by(!!sym(year_var), cd4_category) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(!!sym(year_var)) %>%
    mutate(total = sum(count),
           proportion = count / total * 100) %>%
    ungroup()
  
  # Create a complete grid of all year x category combinations to ensure no gaps
  years <- sort(unique(summary_data[[year_var]]))
  category_levels <- levels(plot_data$cd4_category)
  
  complete_grid <- expand.grid(
    year_diag = years,
    cd4_category = category_levels,
    stringsAsFactors = FALSE
  )
  
  # Make cd4_category a factor with the correct levels
  complete_grid$cd4_category <- factor(complete_grid$cd4_category, levels = category_levels)
  
  # Join with the summarized data to fill in any missing combinations
  summary_data_complete <- complete_grid %>%
    left_join(summary_data, by = c("year_diag" = year_var, "cd4_category")) %>%
    # Replace NA values with zeros
    mutate(
      count = ifelse(is.na(count), 0, count)
    ) %>%
    group_by(year_diag) %>%
    # Recalculate total and proportion
    mutate(
      total = sum(count, na.rm = TRUE),
      proportion = ifelse(total > 0, count / total * 100, 0)
    ) %>%
    ungroup()
  
  # Check if we have any "Unknown" category data
  has_unknown <- any(summary_data_complete$cd4_category == "Unknown" & summary_data_complete$count > 0)
  
  # Filter out the "Unknown" category if it has no data
  if (!has_unknown) {
    summary_data_complete <- summary_data_complete %>%
      filter(cd4_category != "Unknown")
  }
  
  # Adjust colors based on whether we have "Unknown" data
  # Now we have 5 CD4 categories plus potentially Unknown
  plot_colors <- if (has_unknown) colors[1:7] else colors[1:6]  # Just the 5 CD4 categories
  
  
  # Safety check - make sure we have enough colors
  num_categories_needed <- length(unique(summary_data_complete$cd4_category))
  if (length(plot_colors) < num_categories_needed) {
    # If not enough colors provided, extend with default colors
    missing_colors <- num_categories_needed - length(plot_colors)
    default_colors <- c("#999999", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")
    plot_colors <- c(plot_colors, default_colors[1:missing_colors])
  }
  
  # Get unique years for x-axis
  years <- sort(unique(summary_data_complete$year_diag))
  
  # Calculate counts for each year reliably
  year_counts <- numeric(length(years))
  
  
  for (i in 1:length(years)) {
    year_counts[i] <- sum(summary_data$count[summary_data[[year_var]] == years[i]])
  }
  
  if (show_n) {
    x_labels <- paste0(years, "\nn=", year_counts)
  } else {
    x_labels <- as.character(years)
  }
  
  # Create the base plot
  p <- ggplot(summary_data_complete, aes(x = year_diag, y = proportion, fill = cd4_category)) +
    # Use geom_area with group parameter to ensure continuous areas
    geom_area(position = "stack", alpha = 1, aes(group = cd4_category)) +
    scale_fill_manual(values = plot_colors) +
    # Add exact x and y limits
    scale_x_continuous(
      breaks = years,
      labels = x_labels,
      expand = c(0, 0),
      limits = c(min(years) - 0.0005, max(years) + 0.0005)
    ) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      breaks = seq(0, 100, 20),
      expand = c(0, 0),
      limits = c(0, 100.1),
      name = "Proportion (%)"
    )
  
  # Add median CD4 line and secondary y-axis if requested
  if (dash_line) {
    # Calculate reasonable y-axis limits for CD4 counts
    y2_max <- max(median_cd4_by_year$median_cd4, na.rm = TRUE) * 1.2
    y2_max <- ceiling(y2_max / 100) * 100  # Round up to nearest 100
    
    # Identify first and last years for special text positioning
    first_year <- min(median_cd4_by_year[[year_var]])
    last_year <- max(median_cd4_by_year[[year_var]])
    
    # Prepare data for label positioning
    median_cd4_labels <- median_cd4_by_year %>%
      mutate(
        # Adjust horizontal position for first and last year
        label_x = case_when(
          !!sym(year_var) == first_year ~ first_year + 0.01,  # Move first year label right
          !!sym(year_var) == last_year ~ last_year - 0.01,    # Move last year label left
          TRUE ~ !!sym(year_var)                             # Keep others centered
        ),
        # Adjust label position for vertical nudge
        label_hjust = case_when(
          !!sym(year_var) == first_year ~ 0,  # Left-align for first year label
          !!sym(year_var) == last_year ~ 1,   # Right-align for last year label
          TRUE ~ 0.5                         # Center-align for others
        )
      )
    
    # Add secondary y-axis and dash line
    p <- p + 
      geom_line(
        data = median_cd4_by_year, 
        aes(x = !!sym(year_var), y = median_cd4 / y2_max * 100),
        linetype = "dashed", 
        color = "black", 
        size = 1.2,
        inherit.aes = FALSE  # Important: Don't inherit the fill aesthetic
      ) +
      geom_point(
        data = median_cd4_by_year, 
        aes(x = !!sym(year_var), y = median_cd4 / y2_max * 100),
        color = "black", 
        size = 3,
        inherit.aes = FALSE  # Don't inherit aesthetics
      ) +
      # Add point labels with CD4 values - adjusted positioning
      geom_text(
        data = median_cd4_labels, 
        aes(
          x = label_x, 
          y = median_cd4 / y2_max * 100, 
          label = round(median_cd4),
          hjust = label_hjust
        ),
        nudge_y = 5,
        size = 3.5,
        inherit.aes = FALSE  # Don't inherit aesthetics
      ) +
      # Add secondary y-axis for CD4 values
      scale_y_continuous(
        labels = function(x) paste0(x, "%"),
        breaks = seq(0, 100, 20),
        expand = c(0, 0),
        limits = c(-0.1, 100.1),
        name = "Proportion (%)",
        sec.axis = sec_axis(
          ~ . * y2_max / 100, 
          name = "Median CD4 Count (cells/μL)",
          breaks = seq(0, y2_max, 100)
        )
      )
  }
  
  # Complete the plot with labels and theme
  p <- p + 
    labs(
      title = title,
      subtitle = subtitle,
      x = "Year",
      fill = "CD4 Count"
    ) +
    # Clean theme with minimal gridlines
    theme_minimal() +
    theme(
      legend.position = legend_position,
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.border = element_rect(fill = NA, color = "gray80"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      # Add proper margin for the x-axis labels with counts
      axis.text.x.bottom = element_text(margin = margin(t = 5, b = 5)),
      # Clean white background
      plot.background = element_rect(fill = "white", color = NA),
      # Reasonable margins - increased right margin when dash_line is TRUE
      plot.margin = margin(t = 20, r = if(dash_line) 50 else 30, b = 20, l = 20, unit = "pt")
    ) +
    # Use regular coordinate system with no clipping
    coord_cartesian(clip = "off")
  
  return(list(p=p,
              summary_data_complete = summary_data_complete))
}

####Function to create a rainbow plot for absolute counts of CD4 count at diagnosis####
create_cd4_absolute_rainbow_plot <- function(data, 
                                             start_year = NULL,
                                             end_year = NULL, 
                                             show_n = TRUE,          # Show sample size on x-axis
                                             # Filtering options
                                             previ_diag_overseas = NULL,
                                             exposure_filter = NULL,       # hiv_exposure_risk
                                             gender_filter = NULL,         # gender
                                             region_birth_filter = NULL,   # region_birth
                                             lang_home_filter = NULL,      # lang_home
                                             place_acq_filter = NULL,      # place_acq
                                             state_res_filter = NULL,      # state_res
                                             phn_res_filter = NULL,        # phn_res
                                             cd4_var = "cd4_count_diag",
                                             year_var = "year_diag",
                                             title = "Absolute CD4 Count at HIV Diagnosis by Year",
                                             subtitle = NULL,
                                             # Colors for the CD4 categories from bottom to top
                                             colors = c("#999999","#CC79A7","#0072B2","#009E73","#F0E442","#E69F00","#D55E00"),
                                             legend_position = "bottom",
                                             dash_line = FALSE) {          # Option to add median CD4 line
  
  # Load required libraries
  require(tidyverse)
  require(lubridate)
  
  # Clone the data to avoid modifying the original
  plot_data <- data %>% as.data.frame()
  
  plot_data <- plot_data %>%
    filter(
      #      !is.na(!!sym(cd4_var)),  # Ensure CD4 count is not NA
      !is.na(!!sym(year_var))  # Ensure year is not NA
    )
  
  # Filter by year if specified
  if (!is.null(start_year)) {
    plot_data <- plot_data %>% 
      filter(!!sym(year_var) >= start_year)
  }
  
  if (!is.null(end_year)) {
    plot_data <- plot_data %>% 
      filter(!!sym(year_var) <= end_year)
  }
  # Apply filter for previously diagnosed overseas if specified
  if (!is.null(previ_diag_overseas)) {
    if (previ_diag_overseas == TRUE) {
      plot_data <- plot_data %>%
        filter(previ_diag_overseas == 1)
    } else if (previ_diag_overseas == FALSE) {
      plot_data <- plot_data %>%
        filter(previ_diag_overseas == 0)
    }
  }
  
  # Apply filters based on various demographic variables
  # HIV exposure risk
  if (!is.null(exposure_filter)) {
    if ("hiv_exposure_risk" %in% names(plot_data)) {
      plot_data <- plot_data %>% 
        filter(hiv_exposure_risk == exposure_filter)
    } else {
      warning("hiv_exposure_risk variable not found in dataset")
    }
  }
  
  # Gender
  if (!is.null(gender_filter)) {
    if ("gender" %in% names(plot_data)) {
      plot_data <- plot_data %>% 
        filter(gender == gender_filter)
    } else {
      warning("gender variable not found in dataset")
    }
  }
  
  # Region of birth
  if (!is.null(region_birth_filter)) {
    if ("region_birth" %in% names(plot_data)) {
      plot_data <- plot_data %>% 
        filter(region_birth == region_birth_filter)
    } else {
      warning("region_birth variable not found in dataset")
    }
  }
  
  # Language spoken at home
  if (!is.null(lang_home_filter)) {
    if ("lang_home" %in% names(plot_data)) {
      plot_data <- plot_data %>% 
        filter(lang_home == lang_home_filter)
    } else {
      warning("lang_home variable not found in dataset")
    }
  }
  
  # Place of acquisition
  if (!is.null(place_acq_filter)) {
    if ("place_acq" %in% names(plot_data)) {
      plot_data <- plot_data %>% 
        filter(place_acq == place_acq_filter)
    } else {
      warning("place_acq variable not found in dataset")
    }
  }
  
  # State of residence
  if (!is.null(state_res_filter)) {
    if ("state_res" %in% names(plot_data)) {
      plot_data <- plot_data %>% 
        filter(state_res == state_res_filter)
    } else {
      warning("state_res variable not found in dataset")
    }
  }
  
  # PHN of residence
  if (!is.null(phn_res_filter)) {
    if ("phn_res" %in% names(plot_data)) {
      plot_data <- plot_data %>% 
        filter(phn_res == phn_res_filter)
    } else {
      warning("phn_res variable not found in dataset")
    }
  }
  
  # Check if we have data left after all the filtering
  if (nrow(plot_data) == 0) {
    stop("No data remaining after applying all filters. Please adjust your filter criteria.")
  }
  
  # Calculate median CD4 count by year (if dash_line is TRUE)
  if (dash_line) {
    median_cd4_by_year <- plot_data %>%
      group_by(!!sym(year_var)) %>%
      summarise(
        median_cd4 = median(!!sym(cd4_var), na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  # Create CD4 categories with the same categories as the original function
  # Create CD4 categories - MODIFIED to split <200 into <50 and 50-199
  plot_data <- plot_data %>%
    mutate(cd4_category = case_when(
      newly_acq_hiv_flag == 1 & !!sym(cd4_var) < 350 ~ "Seroconversion",
      !!sym(cd4_var) < 50 ~ "< 50",           # New category
      !!sym(cd4_var) >= 50 & !!sym(cd4_var) < 200 ~ "50-199", # New category
      !!sym(cd4_var) >= 200 & !!sym(cd4_var) < 350 ~ "200-349",
      !!sym(cd4_var) >= 350 & !!sym(cd4_var) < 500 ~ "350-499",
      !!sym(cd4_var) >= 500 ~ "≥ 500",
      TRUE ~ "Unknown"
    ))
  
  # Set the order of CD4 categories - UPDATED for new categories
  plot_data$cd4_category <- factor(
    plot_data$cd4_category,
    levels = c("Unknown","Seroconversion", "≥ 500", "350-499","200-349", "50-199","< 50")
  )
  
  # Create summary data for plotting - now using absolute counts
  summary_data <- plot_data %>%
    group_by(!!sym(year_var), cd4_category) %>%
    summarise(count = n(), .groups = "drop") %>%
    ungroup()
  
  # Create a complete grid of all year x category combinations to ensure no gaps
  years <- sort(unique(plot_data[[year_var]]))
  category_levels <- levels(plot_data$cd4_category)
  
  complete_grid <- expand.grid(
    year_diag = years,
    cd4_category = category_levels,
    stringsAsFactors = FALSE
  )
  
  # Make cd4_category a factor with the correct levels
  complete_grid$cd4_category <- factor(complete_grid$cd4_category, levels = category_levels)
  
  # Join with the summarized data to fill in any missing combinations
  summary_data_complete <- complete_grid %>%
    left_join(summary_data, by = c("year_diag" = year_var, "cd4_category")) %>%
    # Replace NA values with zeros
    mutate(
      count = ifelse(is.na(count), 0, count)
    )
  
  # Check if we have any "Unknown" category data
  has_unknown <- any(summary_data_complete$cd4_category == "Unknown" & summary_data_complete$count > 0)
  
  # Filter out the "Unknown" category if it has no data
  if (!has_unknown) {
    summary_data_complete <- summary_data_complete %>%
      filter(cd4_category != "Unknown")
  }
  
  # Adjust colors based on whether we have "Unknown" data
  plot_colors <- if (has_unknown) colors[1:7] else colors[1:6]
  
  # Safety check - make sure we have enough colors
  num_categories_needed <- length(unique(summary_data_complete$cd4_category))
  if (length(plot_colors) < num_categories_needed) {
    # If not enough colors provided, extend with default colors
    missing_colors <- num_categories_needed - length(plot_colors)
    default_colors <- c("#999999", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")
    plot_colors <- c(plot_colors, default_colors[1:missing_colors])
  }
  
  # Get unique years for x-axis
  years <- sort(unique(summary_data_complete$year_diag))
  
  # Calculate total counts for each year to display in labels
  year_totals <- summary_data_complete %>%
    group_by(year_diag) %>%
    summarise(total = sum(count, na.rm = TRUE), .groups = "drop")
  
  # Prepare x-axis labels with sample size
  if (show_n) {
    x_labels <- paste0(years, "\nn=", year_counts)
  } else {
    x_labels <- as.character(years)
  }
  
  # Calculate the maximum count sum to set appropriate y-axis limits
  max_year_count <- max(year_totals$total)
  max_y_value <- ceiling(max_year_count * 1.1 / 50) * 50  # Round up to nearest 50 with 10% buffer
  
  # Create the base plot - now using absolute counts
  p <- ggplot(summary_data_complete, aes(x = year_diag, y = count, fill = cd4_category)) +
    # Use geom_area with group parameter to ensure continuous areas
    geom_area(position = "stack", alpha = 1, aes(group = cd4_category)) +
    scale_fill_manual(values = plot_colors) +
    # Add exact x and y limits
    scale_x_continuous(
      breaks = years,
      labels = x_labels,
      expand = c(0, 0),
      limits = c(min(years) - 0.0005, max(years) + 0.0005)
    ) +
    scale_y_continuous(
      breaks = seq(0, max_y_value, by = ifelse(max_y_value > 500, 100, 50)),
      expand = c(0, 0),
      limits = c(0, max_y_value),
      name = "Number of Diagnoses"
    )
  
  # Add median CD4 line and secondary y-axis if requested
  if (dash_line) {
    # Calculate reasonable y-axis limits for CD4 counts
    y2_max <- max(median_cd4_by_year$median_cd4, na.rm = TRUE) * 1.2
    y2_max <- ceiling(y2_max / 100) * 100  # Round up to nearest 100
    
    # Identify first and last years for special text positioning
    first_year <- min(median_cd4_by_year[[year_var]])
    last_year <- max(median_cd4_by_year[[year_var]])
    
    # Prepare data for label positioning
    median_cd4_labels <- median_cd4_by_year %>%
      mutate(
        # Adjust horizontal position for first and last year
        label_x = case_when(
          !!sym(year_var) == first_year ~ first_year + 0.01,  # Move first year label right
          !!sym(year_var) == last_year ~ last_year - 0.01,    # Move last year label left
          TRUE ~ !!sym(year_var)                             # Keep others centered
        ),
        # Adjust label position for vertical nudge
        label_hjust = case_when(
          !!sym(year_var) == first_year ~ 0,  # Left-align for first year label
          !!sym(year_var) == last_year ~ 1,   # Right-align for last year label
          TRUE ~ 0.5                         # Center-align for others
        )
      )
    
    # Add secondary y-axis and dash line
    p <- p + 
      geom_line(
        data = median_cd4_by_year, 
        aes(x = !!sym(year_var), y = median_cd4 / y2_max * max_y_value),
        linetype = "dashed", 
        color = "black", 
        size = 1.2,
        inherit.aes = FALSE  # Important: Don't inherit the fill aesthetic
      ) +
      geom_point(
        data = median_cd4_by_year, 
        aes(x = !!sym(year_var), y = median_cd4 / y2_max * max_y_value),
        color = "black", 
        size = 3,
        inherit.aes = FALSE  # Don't inherit aesthetics
      ) +
      # Add point labels with CD4 values - adjusted positioning
      geom_text(
        data = median_cd4_labels, 
        aes(
          x = label_x, 
          y = median_cd4 / y2_max * max_y_value, 
          label = round(median_cd4),
          hjust = label_hjust
        ),
        nudge_y = max_y_value * 0.05,  # Nudge labels up by 5% of the y-axis range
        size = 3.5,
        inherit.aes = FALSE  # Don't inherit aesthetics
      ) +
      # Add secondary y-axis for CD4 values
      scale_y_continuous(
        breaks = seq(0, max_y_value, by = ifelse(max_y_value > 500, 100, 50)),
        expand = c(0, 0),
        limits = c(0, max_y_value),
        name = "Number of Diagnoses",
        sec.axis = sec_axis(
          ~ . * y2_max / max_y_value, 
          name = "Median CD4 Count (cells/μL)",
          breaks = seq(0, y2_max, 100)
        )
      )
  }
  
  # Complete the plot with labels and theme
  p <- p + 
    labs(
      title = title,
      subtitle = subtitle,
      x = "Year",
      fill = "CD4 Count"
    ) +
    # Clean theme with minimal gridlines
    theme_minimal() +
    theme(
      legend.position = legend_position,
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.border = element_rect(fill = NA, color = "gray80"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      # Add proper margin for the x-axis labels with counts
      axis.text.x.bottom = element_text(margin = margin(t = 5, b = 5)),
      # Clean white background
      plot.background = element_rect(fill = "white", color = NA),
      # Reasonable margins - increased right margin when dash_line is TRUE
      plot.margin = margin(t = 20, r = if(dash_line) 50 else 30, b = 20, l = 20, unit = "pt")
    ) +
    # Use regular coordinate system with no clipping
    coord_cartesian(clip = "off")
  
  return(p)
}
