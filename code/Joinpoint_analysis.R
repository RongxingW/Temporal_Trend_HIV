# Helper function to select every nth element (for plot labels)
# Keeps only every n-th label (or hides every n-th) so time-axis labels stay readable in plots.
EveryNth <- function(x, n, inverse = FALSE) {
  if (inverse) {
    if (n == 0) return(x)
    # Return the elements that are NOT every nth
    result <- x
    indices <- seq(from = n, to = length(x), by = n)
    result[indices] <- ""
    return(result)
  } else {
    # Return only every nth element
    result <- rep("", length(x))
    indices <- seq(from = n, to = length(x), by = n)
    result[indices] <- x[indices]
    return(result)
  }
}

# Function to get data in correct format for trend analysis
# Updated to handle YYYY-MM format properly
get_trend_data <- function(data, 
                           start_year = NULL,
                           end_year = NULL,
                           # Filtering options
                           exposure_filter = NULL,       # hiv_exposure_risk
                           gender_filter = NULL,         # gender
                           region_birth_filter = NULL,   # region_birth
                           au_birth_filter = NULL,       # au_birth
                           lang_home_filter = NULL,      # lang_home
                           place_acq_filter = NULL,      # place_acq
                           state_res_filter = NULL,      # state_res
                           phn_res_filter = NULL,        # phn_res
                           cd4_count_filter = NULL,     # cd4_count_diag
                           count_period = "monthly",     # "monthly", "quarterly", or "annually"
                           late_diag_filter = NULL,      # late_hiv_diag_flag\
                           advanced_diag_filter = NULL,       # advanced_hiv_diag_flag
                           include_overseas = FALSE,
                           CD4 = FALSE,                  # CD4 count at diagnosis available
                           cd4_var = "cd4_count_diag") {    # Include previously diagnosed overseas
  
  # Clone the data to avoid modifying the original
  plot_data <- data %>% as.data.frame()
  
  # Check if year_diag already exists, if not create it
  if (!"year_diag" %in% names(plot_data)) {
    plot_data <- plot_data %>%
      mutate(year_diag = as.numeric(substr(month_hiv_diag, 1, 4)))
  }
  
  # Filter by year if specified
  if (!is.null(start_year)) {
    plot_data <- plot_data %>% 
      filter(year_diag >= start_year)
  }
  
  if (!is.null(end_year)) {
    plot_data <- plot_data %>% 
      filter(year_diag <= end_year)
  }
  
  # Apply filters based on various demographic variables
  # HIV exposure risk
  if (!is.null(exposure_filter)) {
    plot_data <- plot_data %>% 
      filter(hiv_exposure_risk %in% exposure_filter)
  }
  
  # Gender
  if (!is.null(gender_filter)) {
    plot_data <- plot_data %>% 
      filter(sex == gender_filter)
  }
  
  # Region of birth
  if (!is.null(region_birth_filter)) {
    plot_data <- plot_data %>% 
      filter(region_birth == region_birth_filter)
  }
  
  # Region of birth
  if (!is.null(au_birth_filter)) {
    plot_data <- plot_data %>% 
      filter(au_birth == au_birth_filter)
  }
  
  # Language spoken at home
  if (!is.null(lang_home_filter)) {
    plot_data <- plot_data %>% 
      filter(lang_home == lang_home_filter)
  }
  
  # Place of acquisition
  if (!is.null(place_acq_filter)) {
    plot_data <- plot_data %>% 
      filter(place_acq == place_acq_filter)
  }
  
  # State of residence
  if (!is.null(state_res_filter)) {
    plot_data <- plot_data %>% 
      filter(state_res == state_res_filter)
  }
  
  # PHN of residence
  if (!is.null(phn_res_filter)) {
    plot_data <- plot_data %>% 
      filter(phn_res == phn_res_filter)
  }
  
  # cd4 count at diagnosis
  if (!is.null(cd4_count_filter)) {
    plot_data <- plot_data %>% 
      filter(cd4_count_diag %in% cd4_count_filter)
  }
  
  # late HIV diagnosis
  if (!is.null(late_diag_filter)) {
    plot_data <- plot_data %>% 
      filter(late_hiv_diag_flag == late_diag_filter)
  }
  
  # advanced HIV diagnosis
  if (!is.null(advanced_diag_filter)) {
    plot_data <- plot_data %>% 
      filter(advanced_hiv_diag_flag == advanced_diag_filter)
  }
  
  # Filter by overseas diagnosis if specified
  if (!is.null(include_overseas)) {
    if (include_overseas == TRUE) {
      plot_data <- plot_data %>%
        filter(previ_diag_overseas == 1)
    } else if (include_overseas == FALSE) {
      plot_data <- plot_data %>%
        filter(previ_diag_overseas == 0)
    }
  }
  
  # If CD4 mode, keep only rows with CD4 available
  if (CD4) {
    if (!cd4_var %in% names(plot_data)) {
      stop(sprintf("Variable '%s' not found in data.", cd4_var))
    }
    plot_data <- plot_data %>% filter(!is.na(.data[[cd4_var]]))
  }
  
  # Check if we have data left after all the filtering
  if (nrow(plot_data) == 0) {
    stop("No data remaining after applying all filters. Please adjust your filter criteria.")
  }
  
  # Convert month_hiv_diag (YYYY-MM format) to Date for time series analysis
  plot_data <- plot_data %>%
    mutate(
      # Convert YYYY-MM to Date by adding "-01" for the day
      month_date = as.Date(paste0(month_hiv_diag, "-01")),
      # Convert to decimal date for time series
      decimal_month = lubridate::decimal_date(month_date),
      # Extract quarter
      quarter = paste0(year_diag, "-Q", lubridate::quarter(month_date))
    )
  
  # Helper: summariser depends on CD4 vs count
  summarise_fun <- if (CD4) {
    # use mean CD4 (keep output column named 'notifications' to avoid changing downstream code)
    function(df) dplyr::summarise(df, notifications = mean(.data[[cd4_var]], na.rm = TRUE))
  } else {
    function(df) dplyr::summarise(df, notifications = dplyr::n())
  }
  
  if (count_period == "monthly") {
    trend_data <- plot_data %>%
      group_by(month_hiv_diag) %>%
      summarise_fun() %>%
      mutate(
        month_date = as.Date(paste0(month_hiv_diag, "-01")),
        time = lubridate::decimal_date(month_date)
      ) %>%
      dplyr::select(time, notifications) %>%
      mutate(period = "month") %>%
      arrange(time) %>%
      ungroup()
  } else if (count_period == "quarterly") {
    trend_data <- plot_data %>%
      group_by(quarter) %>%
      summarise_fun() %>%
      mutate(
        year = as.numeric(substr(quarter, 1, 4)),
        q_num = as.numeric(substr(quarter, 7, 7)),
        time = year + (q_num - 1) * 0.25 + 0.125
      ) %>%
      dplyr::select(time, notifications) %>%
      mutate(period = "quarter") %>%
      arrange(time) %>%
      ungroup()
  } else if (count_period == "annually") {
    trend_data <- plot_data %>%
      group_by(year_diag) %>%
      summarise_fun() %>%
      rename(time = year_diag) %>%
      mutate(period = "annual") %>%
      arrange(time) %>%
      ungroup()
  } else {
    stop("Unknown count period. Use 'monthly', 'quarterly', or 'annually'")
  }
  
  # Attach a small attribute so downstream knows we're in CD4 mode
  attr(trend_data, "cd4_mode") <- CD4
  return(trend_data)
}


# Build period-level CD4 quantiles with the SAME filters
get_cd4_quantile_trend_data <- function(
    data,
    start_year = NULL,
    end_year = NULL,
    exposure_filter = NULL,      # hiv_exposure_risk
    gender_filter = NULL,        # gender (sex)
    region_birth_filter = NULL,  # region_birth
    au_birth_filter = NULL,      # au_birth
    lang_home_filter = NULL,     # lang_home
    place_acq_filter = NULL,     # place_acq
    state_res_filter = NULL,     # state_res
    phn_res_filter = NULL,       # phn_res
    include_overseas = FALSE,    # previously diagnosed overseas
    count_period = "monthly",    # "monthly", "quarterly", "annually"
    cd4_var = "cd4_count_diag",
    probs = c(0.05, 0.50, 0.95)
) {
  require(dplyr); require(lubridate); require(rlang)
  
  plot_data <- data %>% as.data.frame()
  
  # Derive year if missing
  if (!"year_diag" %in% names(plot_data)) {
    plot_data <- plot_data %>%
      mutate(year_diag = as.numeric(substr(month_hiv_diag, 1, 4)))
  }
  
  # Year window
  if (!is.null(start_year)) plot_data <- plot_data %>% filter(year_diag >= start_year)
  if (!is.null(end_year))   plot_data <- plot_data %>% filter(year_diag <= end_year)
  
  # Filters
  if (!is.null(exposure_filter))    plot_data <- plot_data %>% filter(hiv_exposure_risk %in% exposure_filter)
  if (!is.null(gender_filter))      plot_data <- plot_data %>% filter(sex == gender_filter)
  if (!is.null(region_birth_filter))plot_data <- plot_data %>% filter(region_birth == region_birth_filter)
  if (!is.null(au_birth_filter))    plot_data <- plot_data %>% filter(au_birth == au_birth_filter)
  if (!is.null(lang_home_filter))   plot_data <- plot_data %>% filter(lang_home == lang_home_filter)
  if (!is.null(place_acq_filter))   plot_data <- plot_data %>% filter(place_acq == place_acq_filter)
  if (!is.null(state_res_filter))   plot_data <- plot_data %>% filter(state_res == state_res_filter)
  if (!is.null(phn_res_filter))     plot_data <- plot_data %>% filter(phn_res == phn_res_filter)
  
  # Overseas switch
  if (!is.null(include_overseas)) {
    if (include_overseas == TRUE)  plot_data <- plot_data %>% filter(previ_diag_overseas == 1)
    if (include_overseas == FALSE) plot_data <- plot_data %>% filter(previ_diag_overseas == 0)
  }
  
  # Keep CD4 rows only
  if (!cd4_var %in% names(plot_data)) stop(sprintf("Variable '%s' not found.", cd4_var))
  plot_data <- plot_data %>% filter(!is.na(.data[[cd4_var]]))
  
  if (nrow(plot_data) == 0) stop("No CD4 data remaining after filters.")
  
  # Dates & quarter label
  plot_data <- plot_data %>%
    mutate(
      month_date = as.Date(paste0(month_hiv_diag, "-01")),
      quarter    = paste0(year_diag, "-Q", lubridate::quarter(month_date))
    )
  
  # Summariser
  qfun <- function(x) as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE))
  
  if (count_period == "monthly") {
    out <- plot_data %>%
      group_by(month_hiv_diag) %>%
      summarise(q = list(qfun(.data[[cd4_var]])), .groups = "drop") %>%
      mutate(
        month_date = as.Date(paste0(month_hiv_diag, "-01")),
        time = lubridate::decimal_date(month_date)
      )
  } else if (count_period == "quarterly") {
    out <- plot_data %>%
      group_by(quarter) %>%
      summarise(q = list(qfun(.data[[cd4_var]])), .groups = "drop") %>%
      mutate(
        year  = as.numeric(substr(quarter, 1, 4)),
        q_num = as.numeric(substr(quarter, 7, 7)),
        time = year + (q_num - 1) * 0.25 + 0.125
      )
  } else if (count_period == "annually") {
    out <- plot_data %>%
      group_by(year_diag) %>%
      summarise(q = list(qfun(.data[[cd4_var]])), .groups = "drop") %>%
      rename(time = year_diag)
  } else {
    stop("count_period must be 'monthly', 'quarterly', or 'annually'.")
  }
  
  # Unnest to wide (p05, p50, p95)
  out <- out %>%
    dplyr::mutate(
      p05 = vapply(q, `[`, 1, FUN.VALUE = numeric(1)),
      p50 = vapply(q, `[`, 2, FUN.VALUE = numeric(1)),
      p95 = vapply(q, `[`, 3, FUN.VALUE = numeric(1))
    ) %>%
    dplyr::select(time, p05, p50, p95) %>%
    dplyr::arrange(time)
  
  attr(out, "cd4_mode") <- TRUE
  return(out)
}


# Fits the overall time trend using Poisson or Negative Binomial 
# and reports the estimated change, confidence interval, and p-value.
get_trend <- function(trend_data, family_type = "poisson") {
  
  # If trend_data came from CD4 mode, prefer gaussian
  cd4_mode <- isTRUE(attr(trend_data, "cd4_mode"))
  if (cd4_mode && tolower(family_type) %in% c("poisson","nb")) {
    family_type <- "gaussian"
  }
  
  # --- Gaussian path (for CD4 averages) ---
  if (tolower(family_type) %in% c("gaussian","linear","gauss")) {
    gaus_model <- glm(
      notifications ~ time,
      data = trend_data,
      family = gaussian(link = "identity"),
      na.action = na.exclude,
      model = TRUE, x = TRUE, y = TRUE
    )
    beta <- coef(gaus_model)["time"]
    se   <- summary(gaus_model)$coefficients["time", "Std. Error"]
    zcrit <- qnorm(0.975)
    ci_lin <- c(beta - zcrit * se, beta + zcrit * se)
    p_value <- 2 * pnorm(-abs(beta / se))
    
    # For Gaussian, 'coeff' is the slope (Δ CD4 per year)
    return(list(
      model = gaus_model,
      coeff = unname(beta),
      interval = unname(ci_lin),
      p_value = p_value,
      type_used = "gaussian",
      dispersion = NA_real_
    ))
  }
  if (family_type %in% c("poisson", "POISSON")) {
    pois_model <- glm(
      notifications ~ time,
      data = trend_data,
      family = poisson(),
      na.action = na.exclude,
      model = TRUE, x = TRUE, y = TRUE
    )
    dispersion <- sum(residuals(pois_model, type = "pearson")^2) / pois_model$df.residual
    if (is.finite(dispersion) && dispersion > 1.0) {
      stop(
        sprintf(
          "Overdispersion detected (phi = %.2f) under Poisson. Recommended to use family_type = 'NB' (negative binomial).",
          dispersion
        )
      )
    }
    beta <- coef(pois_model)["time"]
    se   <- summary(pois_model)$coefficients["time", "Std. Error"]
    zcrit <- qnorm(0.975)
    ci_lin <- c(beta - zcrit * se, beta + zcrit * se)
    p_value <- 2 * pnorm(-abs(beta / se))
    coeff   <- exp(beta)
    interval <- exp(ci_lin)
    
    return(list(
      model = pois_model,
      coeff = coeff,
      interval = interval,
      p_value = p_value,
      type_used = "poisson",
      dispersion = dispersion
    ))
  }
  
  if (tolower(family_type) == "nb") {
    nb_model <- MASS::glm.nb(
      notifications ~ time,
      data = trend_data,
      link = log,                  # explicit (default is log, but be clear)
      na.action = na.exclude,
      model = TRUE, x = TRUE, y = TRUE,
      control = glm.control(trace = FALSE)
    )
    beta <- coef(nb_model)["time"]
    se   <- summary(nb_model)$coefficients["time", "Std. Error"]
    zcrit <- qnorm(0.975)
    ci_lin <- c(beta - zcrit * se, beta + zcrit * se)
    p_value <- 2 * pnorm(-abs(beta / se))
    coeff   <- exp(beta)
    interval <- exp(ci_lin)
    
    return(list(
      model = nb_model,
      coeff = coeff,
      interval = interval,
      p_value = p_value,
      type_used = "NB",
      dispersion = NA_real_
    ))
  }
  
  stop("Unknown 'family_type'. Use 'poisson' or 'NB'.")
}


# Function for generating overall trend plot
overall_trend_plot <- function(trend_data, trend_model, count_period = "monthly") {
  cd4_mode <- isTRUE(attr(trend_data, "cd4_mode"))
  
  # x-axis label and breaks
  x_label_vec <- c("monthly" = "Year", "quarterly" = "Year", "annually" = "Year")
  x_label <- x_label_vec[count_period]
  
  if (count_period %in% c("monthly", "quarterly")) {
    if (count_period == "monthly") {
      years <- unique(floor(trend_data$time))
      xbreaks <- years + 0.0
      xlabels <- as.character(years)
    } else {
      xbreaks <- floor(min(trend_data$time)):ceiling(max(trend_data$time))
      xlabels <- xbreaks
    }
  } else {
    xbreaks <- unique(trend_data$time)
    xlabels <- xbreaks
  }
  
  # y label
  y_lab <- if (cd4_mode)
    paste0(stringr::str_to_title(count_period), " average CD4 (cells/μL)")
  else
    paste0(stringr::str_to_title(count_period), " notifications")
  
  # y-axis range
  max_y <- max(trend_data$notifications, na.rm = TRUE)
  auto_upper <- if (is.finite(max_y) && max_y > 0) {
    p <- floor(log10(max_y)); ceiling(max_y / 10^p) * 10^p
  } else 1
  
  # --- Plot ---
  ggplot(trend_data, aes(x = time, y = notifications)) +
    geom_point(aes(color = "data")) +
    geom_line(
      aes(y = predict(trend_model, type = "response"), color = "overall"),
      size = 1
    ) +
    scale_color_manual(
      "",
      values = c(data = "#82afda", overall = "#EE7C7A"),
      labels = c(data = if (cd4_mode) "Observed CD4" else "HIV data",
                 overall = if (cd4_mode) "Overall mean CD4 trend" else "Fitted trend")
    ) +
    scale_y_continuous(
      limits = c(0, auto_upper),
      labels = scales::comma,
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_x_continuous(breaks = xbreaks, labels = xlabels) +
    ylab(y_lab) + xlab(x_label) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
      legend.position = "bottom"
    )
}



cd4_quantile_plot <- function(
    cd4_quant_data,              # from get_cd4_quantile_trend_data()
    count_period = "monthly",
    segment_model = NULL,        # optional: add dashed vlines at joinpoints
    y_upper = NULL,
    y_break = NULL
) {
  require(ggplot2); require(scales)
  
  # X breaks & labels
  x_label_vec <- c("monthly" = "Year", "quarterly" = "Year", "annually" = "Year")
  x_label <- x_label_vec[count_period]
  
  if (count_period %in% c("monthly","quarterly")) {
    if (count_period == "monthly") {
      years <- unique(floor(cd4_quant_data$time))
      xbreaks <- years + 0.0
      xlabels <- as.character(years)
    } else {
      xbreaks <- floor(min(cd4_quant_data$time)):ceiling(max(cd4_quant_data$time))
      xlabels <- xbreaks
    }
  } else {
    xbreaks <- unique(cd4_quant_data$time)
    xlabels <- xbreaks
  }
  
  # y limits
  max_y <- max(cd4_quant_data$p95, na.rm = TRUE)
  auto_upper <- 1500
  y_upper_final <- if (is.null(y_upper)) auto_upper else y_upper
  y_breaks <- if (!is.null(y_break) && is.finite(y_break) && y_break > 0) seq(0, y_upper_final, by = y_break) else waiver()
  
  # Long format for one legend
  df_long <- tidyr::pivot_longer(cd4_quant_data, cols = c(p05, p50, p95),
                                 names_to = "series", values_to = "value")
  series_labs <- c(p05 = "5th", p50 = "Median", p95 = "95th")
  series_cols <- c(p05 = "#9bbf8a", p50 = "#F2AB3F", p95 = "#B0AEC6")  # soft green, salmon, sky-blue
  
  p <- ggplot(df_long, aes(x = time, y = value, color = series)) +
    geom_line(size = 1) +
    scale_color_manual("", values = series_cols, labels = series_labs) +
    ylab(paste0(stringr::str_to_title(count_period), " CD4 (cells/μL)")) +
    xlab(x_label) +
    scale_y_continuous(limits = c(0, y_upper_final), breaks = y_breaks, labels = scales::comma,
                       expand = expansion(mult = c(0, 0.02))) +
    scale_x_continuous(breaks = xbreaks, labels = xlabels) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
      legend.position = "bottom"
    )
  
  # Optional: add dashed lines at joinpoints (same visual cue as segmented plot)
  if (!is.null(segment_model) && !is.null(segment_model$psi) && nrow(segment_model$psi) > 0) {
    for (i in 1:nrow(segment_model$psi)) {
      p <- p + geom_vline(xintercept = segment_model$psi[i, "Est."],
                          linetype = "dashed", colour = "#EE7C7A")
    }
  }
  
  return(p)
}


# Runs the complete “overall trend” workflow: prepares data, 
# fits the model, makes the plot, and returns a tidy summary.
analyze_overall_trend <- function(data, 
                                  start_year = NULL,
                                  end_year = NULL,
                                  # Filtering options
                                  exposure_filter = NULL,      # hiv_exposure_risk
                                  gender_filter = NULL,        # gender
                                  region_birth_filter = NULL,  # region_birth
                                  au_birth_filter = NULL,      # au_birth
                                  lang_home_filter = NULL,     # lang_home
                                  place_acq_filter = NULL,     # place_acq
                                  state_res_filter = NULL,     # state_res
                                  phn_res_filter = NULL,       # phn_res
                                  cd4_count_filter = NULL,     # cd4_count_diag
                                  late_diag_filter = NULL,      # late_hiv_diag_flag
                                  advanced_diag_filter = NULL,       # advanced_hiv_diag_flag
                                  count_period = "monthly",    # "monthly", "quarterly", or "annually"
                                  include_overseas = FALSE,     # Include previously diagnosed overseas
                                  family_type = "poisson",            # "poisson" or "linear"
                                  save_plots = FALSE,          # Whether to save plots
                                  output_dir = NULL,           # Directory to save plots
                                  plot_name_prefix = "hiv_trend",
                                  CD4 = FALSE,                 # cd4 count at diagnosis available
                                  cd4_var = "cd4_count_diag") {
  
  # Check parameters
  if (save_plots && is.null(output_dir)) {
    stop("If save_plots = TRUE, output_dir must be specified.")
  }
  
  # Get trend data with all filters applied
  trend_data <- get_trend_data(
    data = data,
    start_year = start_year,
    end_year = end_year,
    exposure_filter = exposure_filter,
    gender_filter = gender_filter,
    region_birth_filter = region_birth_filter,
    au_birth_filter = au_birth_filter,
    lang_home_filter = lang_home_filter,
    place_acq_filter = place_acq_filter,
    state_res_filter = state_res_filter,
    phn_res_filter = phn_res_filter,
    cd4_count_filter = cd4_count_filter,
    count_period = count_period,
    include_overseas = include_overseas,
    late_diag_filter = late_diag_filter,
    advanced_diag_filter = advanced_diag_filter,
    CD4 = CD4,
    cd4_var = cd4_var
  )
  
  # Create title based on filters
  title_parts <- c()
  if (!is.null(exposure_filter)) title_parts <- c(title_parts, paste0("Exposure: ", exposure_filter))
  if (!is.null(gender_filter)) title_parts <- c(title_parts, paste0("Gender: ", gender_filter))
  if (!is.null(region_birth_filter)) title_parts <- c(title_parts, paste0("Birth region: ", region_birth_filter))
  if (!is.null(au_birth_filter)) title_parts <- c(title_parts, paste0("Birth place: ", au_birth_filter))
  if (!is.null(state_res_filter)) title_parts <- c(title_parts, paste0("State: ", state_res_filter))
  if (!is.null(phn_res_filter)) title_parts <- c(title_parts, paste0("PHN: ", phn_res_filter))
  if(!is.null(late_diag_filter)) title_parts <- c(title_parts, paste0("Late HIV Diagnosis"))
  if(!is.null(advanced_diag_filter)) title_parts <- c(title_parts, paste0("Advanced HIV Diagnosis"))
  if (include_overseas == TRUE) title_parts <- c(title_parts, paste0("Previously diagnosed overseas"))
  if (CD4) title_parts <- c(title_parts, "Outcome: Average CD4")
  if (length(title_parts) > 0) {
    subtitle <- paste(title_parts, collapse = ", ")
  } else {
    subtitle <- "All notifications"
  }
  
  if (!include_overseas) {
    subtitle <- paste0(subtitle, " (Excluding previously diagnosed overseas)")
  }
  
  # Auto-switch family when CD4
  fam <- if (CD4 && tolower(family_type) %in% c("poisson","nb")) "gaussian" else family_type
  
  overall_trend <- get_trend(trend_data, family_type = fam)
  
  # Create overall trend plot
  overall_plot <- overall_trend_plot(trend_data, overall_trend$model, count_period = count_period) +
    ggtitle(paste0(if (CD4) "Average CD4 Trend (" else "HIV Notification Trend (", count_period, ")"),
            subtitle = subtitle)
  
  # Save overall trend plot if requested
  if (save_plots && !is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(
      filename = file.path(output_dir, paste0(plot_name_prefix, "_overall_", count_period, ".png")),
      plot = overall_plot,
      width = 10,
      height = 6,
      dpi = 300
    )
  }
  
  # Results table: interpret based on family
  if (tolower(overall_trend$type_used) == "gaussian") {
    overall_results <- data.frame(
      Model = "Overall Trend (Gaussian)",
      Slope_per_Year = overall_trend$coeff,       # Δ CD4 per year
      Lower_CI = overall_trend$interval[1],
      Upper_CI = overall_trend$interval[2],
      P_value = overall_trend$p_value,
      stringsAsFactors = FALSE
    )
  } else {
    overall_results <- data.frame(
      Model = "Overall Trend",
      Estimate = overall_trend$coeff,
      Lower_CI = overall_trend$interval[1],
      Upper_CI = overall_trend$interval[2],
      P_value = overall_trend$p_value,
      stringsAsFactors = FALSE
    )
  }
  
  # Return results
  return(list(
    trend_data = trend_data,
    overall_trend = overall_trend,
    overall_plot = overall_plot,
    overall_results = overall_results
  ))
}

# Function for manually calculating confidence intervals for segmented model
seg_CI <- function(est, std_err) {
  return(c(est - 1.96 * std_err, est + 1.96 * std_err))
}

# Fits a joinpoint (segmented) model on top of the base trend, 
# returning breakpoints with CIs, segment-specific slopes, and BIC.
get_segmented_trend <- function(overall_model, initial_changes, family_type = "poisson") {
  require(segmented)
  
  base_df <- as.data.frame(stats::model.frame(overall_model))
  if (!is.numeric(base_df$time)) base_df$time <- as.numeric(base_df$time)
  base_formula <- stats::formula(overall_model)
  
  # Decide family for refit
  fam_lower <- tolower(family_type)
  use_gaussian <- fam_lower %in% c("gaussian","linear","gauss")
  use_nb <- fam_lower == "nb"
  
  if (inherits(overall_model, "negbin") || use_nb) {
    nb0 <- tryCatch(
      MASS::glm.nb(base_formula, data = base_df, link = log,
                   na.action = na.exclude, model = TRUE, x = TRUE, y = TRUE,
                   control = glm.control(maxit = 100, trace = FALSE)),
      error = function(e) NULL
    )
    if (is.null(nb0)) stop("NB base model failed; cannot proceed to segmented.")
    theta_hat <- nb0$theta
    re_fit <- stats::glm(
      base_formula, data = base_df,
      family = MASS::negative.binomial(theta_hat),
      na.action = na.exclude, model = TRUE, x = TRUE, y = TRUE,
      control = glm.control(maxit = 100)
    )
  } else if (use_gaussian) {
    re_fit <- stats::glm(
      base_formula, data = base_df, family = gaussian(),
      na.action = na.exclude, model = TRUE, x = TRUE, y = TRUE,
      control = glm.control(maxit = 100)
    )
  } else {
    re_fit <- stats::glm(
      base_formula, data = base_df, family = poisson(),
      na.action = na.exclude, model = TRUE, x = TRUE, y = TRUE,
      control = glm.control(maxit = 100)
    )
  }
  
  segment_trend <- tryCatch({
    segmented.glm(
      re_fit,
      seg.Z = ~ time,
      psi   = list(time = initial_changes),
      control = seg.control(n.boot = 0)
    )
  }, error = function(e) {
    cat("Error in segmented.glm():", conditionMessage(e), "\n")
    return(NULL)
  })
  if (is.null(segment_trend)) return(NULL)
  
  glm_fit <- if (!is.null(segment_trend$glm.fit)) segment_trend$glm.fit else segment_trend
  
  seg_BIC  <- tryCatch(BIC(glm_fit), error = function(e) NA_real_)
  base_BIC <- tryCatch(BIC(re_fit),  error = function(e) NA_real_)
  
  psi <- segment_trend$psi
  psi_interval <- tryCatch({
    confint.segmented(segment_trend)
  }, error = function(e) {
    se_col <- tryCatch({ summary(segment_trend)$psi[, "Std. Error"] }, error = function(e2) rep(0.05, nrow(psi)))
    cbind("Est." = psi[, "Est."],
          "CI(95%).low" = psi[, "Est."] - 1.96 * se_col,
          "CI(95%).up"  = psi[, "Est."] + 1.96 * se_col)
  })
  
  # Slopes on link scale:
  # - Poisson/NB (log link): exponentiate → rate ratio per unit time
  # - Gaussian (identity): DO NOT exponentiate → change in CD4 per unit time
  slopes_list <- tryCatch({
    s <- slope(segment_trend)$time
    if (use_gaussian) {
      list(
        est = s[, "Est."],
        lo  = s[, "CI(95%).l"],
        up  = s[, "CI(95%).u"]
      )
    } else {
      list(
        est = exp(s[, "Est."]),
        lo  = exp(s[, "CI(95%).l"]),
        up  = exp(s[, "CI(95%).u"])
      )
    }
  }, error = function(e) {
    k <- nrow(psi) + 1
    if (use_gaussian) {
      list(est = rep(0, k), lo = rep(-10, k), up = rep(10, k))
    } else {
      list(est = rep(1, k), lo = rep(0.9, k), up = rep(1.1, k))
    }
  })
  
  change_date <- data.frame(
    Estimate = round(psi[, "Est."], 4),
    Lower_CI = round(psi_interval[, "CI(95%).low"], 4),
    Upper_CI = round(psi_interval[, "CI(95%).up"], 4)
  )
  
  wald_tab <- tryCatch({
    coefs <- coef(summary(segment_trend))
    rn <- rownames(coefs)
    idx <- which(grepl("^U\\d*\\.time$", rn))
    if (length(idx) == 0) idx <- which(grepl("^U.*\\.time$", rn))
    if (length(idx) == 0) {
      data.frame(Joinpoint = integer(0), Term = character(0),
                 Estimate = numeric(0), SE = numeric(0), z = numeric(0),
                 p_value = numeric(0), stringsAsFactors = FALSE)
    } else {
      out <- coefs[idx, , drop = FALSE]
      zvals <- out[, "Estimate"] / out[, "Std. Error"]
      pvals <- 2 * pnorm(-abs(zvals))
      jp_time <- if (!is.null(psi)) psi[, "Est."] else rep(NA_real_, length(idx))
      to_month <- function(x) format(lubridate::round_date(lubridate::date_decimal(x), unit = "month"), "%Y-%m")
      data.frame(
        Joinpoint      = seq_len(nrow(out)),
        Term           = rownames(out),
        Estimate       = out[, "Estimate"],
        SE             = out[, "Std. Error"],
        z              = zvals,
        p_value        = pvals,
        JP_Time        = jp_time,
        JP_Time_Month  = if (length(jp_time)) to_month(jp_time) else character(0),
        Exp_Estimate   = if (use_gaussian) NA_real_ else exp(out[, "Estimate"]),
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {
    data.frame(Joinpoint = integer(0), Term = character(0),
               Estimate = numeric(0), SE = numeric(0), z = numeric(0),
               p_value = numeric(0), JP_Time = numeric(0),
               JP_Time_Month = character(0), Exp_Estimate = numeric(0),
               stringsAsFactors = FALSE)
  })
  
  # AAPC equivalent:
  # - Gaussian: weighted average of slopes (absolute Δ CD4 / year)
  # - Poisson/NB: weighted average of log-slopes → exp(...)
  aapc <- tryCatch({
    if (!is.null(psi) && nrow(psi) > 0) {
      time_range <- range(base_df$time, na.rm = TRUE)
      total_time <- diff(time_range)
      breaks <- c(time_range[1], psi[, "Est."], time_range[2])
      seg_lengths <- diff(breaks)
      weights <- seg_lengths / total_time
      
      if (use_gaussian) {
        est <- sum(weights * slopes_list$est)
        lo  <- sum(weights * slopes_list$lo)
        up  <- sum(weights * slopes_list$up)
        list(est = est, lo = lo, up = up)  # units: CD4 per year
      } else {
        log_slopes <- log(slopes_list$est)
        aapc_log <- sum(weights * log_slopes)
        aapc_est <- exp(aapc_log)
        log_slopes_lo <- log(slopes_list$lo)
        log_slopes_up <- log(slopes_list$up)
        aapc_lo <- exp(sum(weights * log_slopes_lo))
        aapc_up <- exp(sum(weights * log_slopes_up))
        list(est = aapc_est, lo = aapc_lo, up = aapc_up)
      }
    } else {
      if (use_gaussian) {
        list(est = slopes_list$est[1], lo = slopes_list$lo[1], up = slopes_list$up[1])
      } else {
        list(est = slopes_list$est[1], lo = slopes_list$lo[1], up = slopes_list$up[1])
      }
    }
  }, error = function(e) list(est = NA, lo = NA, up = NA))
  
  list(
    segment_model = segment_trend,
    psi = psi,
    psi_interval = psi_interval,
    change_date = change_date,
    slopes = slopes_list$est,
    slopes_lower = slopes_list$lo,
    slopes_upper = slopes_list$up,
    bics = c("segment_BIC" = seg_BIC, "overall_BIC" = base_BIC),
    model_type = if (use_gaussian) "gaussian" else if (use_nb) "NB" else "poisson",
    wald_slope_change = wald_tab,
    aapc = aapc
  )
}




# (Helper) Fits a model with exactly k joinpoints and returns that model’s BIC for comparison.
.fit_k_joinpoints <- function(base_model, data_years, k, family_type = "poisson") {
  # k = 0 -> no joinpoint: use base_model
  if (k == 0) {
    return(list(
      k = 0,
      fit = list(
        segment_model = NULL,  # no segmented object
        psi = NULL,
        psi_interval = NULL,
        change_date = NULL,
        slopes = NULL,
        slopes_lower = NULL,
        slopes_upper = NULL,
        bics = c("segment_BIC" = BIC(base_model), "overall_BIC" = BIC(base_model)),
        model_type = ifelse(tolower(family_type) == "nb", "NB", "poisson"),
        # carry the base model for plotting/prediction when needed
        base_model = base_model
      ),
      bic = BIC(base_model)
    ))
  }
  
  # Generate evenly spaced initial psi across observed time support
  min_year <- min(data_years); max_year <- max(data_years)
  initial_points <- seq(min_year, max_year, length.out = k + 2)[2:(k + 1)]
  
  seg_fit <- get_segmented_trend(base_model, initial_points, family_type = family_type)
  if (is.null(seg_fit)) {
    return(list(k = k, fit = NULL, bic = Inf))
  }
  
  seg_bic <- as.numeric(seg_fit$bics["segment_BIC"])
  if (!is.finite(seg_bic)) seg_bic <- Inf
  
  list(k = k, fit = seg_fit, bic = seg_bic)
}

# Tries 0..K joinpoints, compares BICs, and picks the model with the best (lowest) BIC.
find_best_segment_model <- function(model, data_years, num_changes = 1, family_type = "poisson") {
  require(segmented)
  
  # Fit all candidates
  candidates <- lapply(0:num_changes, function(k) .fit_k_joinpoints(model, data_years, k, family_type = family_type))
  
  # Pick the lowest BIC
  bics <- sapply(candidates, function(x) x$bic)
  best_idx <- which.min(bics)
  best <- candidates[[best_idx]]
  
  if (is.null(best$fit)) {
    cat("All segmented attempts failed; returning base model.\n")
    return(candidates[[1]]$fit)  # k=0
  }
  
  # Attach selected_k and a tidy small table of all BICs
  best$fit$selected_k <- best$k
  best$fit$all_bic <- data.frame(
    Joinpoints = 0:num_changes,
    BIC = bics
  )
  return(best$fit)
}

# Plots the data and the chosen segmented fit, marking each estimated joinpoint with a dashed line.
segment_plot <- function(trend_data, trend_model, count_period = "monthly", 
                         overall_plot = FALSE, overall_model = NULL, 
                         segment_plot = TRUE, y_upper = NULL, y_break = NULL) {
  cd4_mode <- isTRUE(attr(trend_data, "cd4_mode"))
  
  x_label_vec <- c("monthly" = "Year", "quarterly" = "Year", "annually" = "Year")
  x_label <- x_label_vec[count_period]
  
  segment_values <- tryCatch({
    data.frame(time = trend_data$time, 
               segvalue = fitted(trend_model$segment_model))
  }, error = function(e) {
    data.frame(time = trend_data$time, 
               segvalue = predict(overall_model, type = "response"))
  })
  
  if (count_period %in% c("monthly","quarterly")) {
    if (count_period == "monthly") {
      years <- unique(floor(trend_data$time))
      xbreaks <- years + 0.0
      xlabels <- as.character(years)
    } else {
      xbreaks <- floor(min(trend_data$time)):ceiling(max(trend_data$time))
      xlabels <- EveryNth(xbreaks, 1, inverse = FALSE)
    }
  } else {
    xbreaks <- unique(trend_data$time)
    xlabels <- xbreaks
  }
  
  max_y <- max(trend_data$notifications, na.rm = TRUE)
  auto_upper <- if (is.finite(max_y) && max_y > 0) {
    p <- floor(log10(max_y)); ceiling(max_y / 10^p) * 10^p
  } else 1
  y_upper_final <- if (is.null(y_upper)) auto_upper else y_upper
  y_breaks <- if (!is.null(y_break) && is.finite(y_break) && y_break > 0) seq(0, y_upper_final, by = y_break) else waiver()
  
  y_lab <- if (cd4_mode) {
    paste0(stringr::str_to_title(count_period), " average CD4 (cells/μL)")
  } else {
    paste0(stringr::str_to_title(count_period), " notifications")
  }
  
  trend_plot <- ggplot(trend_data, aes(x = time, y = notifications)) + 
    geom_point(aes(color = "data")) +
    ylab(y_lab) + xlab(x_label) +
    scale_y_continuous(
      limits = c(0, y_upper_final),
      breaks = y_breaks,
      labels = scales::comma,
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_x_continuous(breaks = xbreaks, labels = xlabels) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
      legend.position = "bottom"
    )
  
  line_layers <- list()
  color_values <- c("data" = "#82afda")
  color_labels <- c("data" = if (cd4_mode) "Mean CD4" else "HIV data")
  
  if (segment_plot && !is.null(trend_model$segment_model)) {
    line_layers <- c(line_layers, geom_line(data = segment_values, aes(x = time, y = segvalue, color = "segment"), size = 1))
    color_values["segment"] <- "#EE7C7A"
    color_labels["segment"] <- "Fitted trends"
  }
  if (overall_plot && !is.null(overall_model)) {
    line_layers <- c(line_layers, geom_line(aes(y = predict(overall_model, type = "response"), color = "overall"), size = 1))
    color_values["overall"] <- "#82afda"
    color_labels["overall"] <- "Overall trend"
  }
  for (layer in line_layers) trend_plot <- trend_plot + layer
  trend_plot <- trend_plot + scale_colour_manual("", values = color_values, labels = color_labels)
  
  if (!is.null(trend_model$psi)) {
    for (i in 1:nrow(trend_model$psi)) {
      trend_plot <- trend_plot + geom_vline(xintercept = trend_model$psi[i, "Est."], linetype = "dashed", colour = "#EE7C7A")
    }
  }
  return(trend_plot)
}


# Time-weighted mean slope and 95% CI for a Gaussian segmented.glm
# Uses linear contrast on the coefficient covariance (statistically correct).
weighted_gaussian_slope_CI <- function(selected, trend_data) {
  stopifnot(!is.null(selected$segment_model))
  seg <- selected$segment_model
  
  tmin <- min(trend_data$time, na.rm = TRUE)
  tmax <- max(trend_data$time, na.rm = TRUE)
  jp   <- if (!is.null(selected$psi) && nrow(selected$psi) > 0) as.numeric(selected$psi[,"Est."]) else numeric(0)
  breaks <- c(tmin, jp, tmax)
  seg_lengths <- diff(breaks)
  w <- seg_lengths / sum(seg_lengths)      # time-weights per segment
  K <- length(w)
  
  theta <- stats::coef(seg)
  V     <- stats::vcov(seg)
  
  nm <- names(theta)
  idx_time <- which(nm == "time")
  if (length(idx_time) != 1) idx_time <- grep("^time$", nm)
  u_idx <- grep("^U\\d+\\.time$", nm)
  
  if (length(idx_time) != 1) stop("Could not find base time slope in segmented coefs.")
  
  # Contrast vector: mu = (sum w) * beta_time + (w2+...+wK)*U1 + (w3+...+wK)*U2 + ...
  cvec <- rep(0, length(theta)); names(cvec) <- nm
  cvec[idx_time] <- sum(w)  # = 1
  if (length(u_idx) > 0) {
    tail_sums <- rev(cumsum(rev(w)))[-1]  # (w2+...+wK, w3+...+wK, ..., wK)
    cvec[u_idx] <- tail_sums[seq_along(u_idx)]
  }
  
  mu <- sum(cvec * theta)
  se <- sqrt(as.numeric(t(cvec) %*% V %*% cvec))
  ci <- c(mu - 1.96*se, mu + 1.96*se)
  
  list(
    weighted_slope = mu,   # Δ CD4 per year (time-weighted mean)
    se = se,
    ci95 = ci,
    weights = w,
    segments = data.frame(
      start = breaks[-length(breaks)],
      end   = breaks[-1],
      weight = w
    )
  )
}


# End-to-end segmented analysis: prepares data, fits overall model, 
# selects joinpoints by BIC, converts joinpoints to months, and returns plots and results.
analyze_segmented_trend_complete <- function(data, 
                                             start_year = NULL,
                                             end_year = NULL,
                                             exposure_filter = NULL,
                                             gender_filter = NULL,
                                             region_birth_filter = NULL,
                                             au_birth_filter = NULL,
                                             lang_home_filter = NULL,
                                             place_acq_filter = NULL,
                                             state_res_filter = NULL,
                                             phn_res_filter = NULL,
                                             cd4_count_filter = NULL,     # cd4_count_diag
                                             late_diag_filter = NULL,
                                             advanced_diag_filter = NULL,
                                             count_period = "monthly",
                                             include_overseas = FALSE,
                                             family_type = "poisson",
                                             num_changes = 1,
                                             show_overall = TRUE,
                                             show_segment = TRUE,
                                             y_upper = NULL, 
                                             y_break = NULL,
                                             plot_name_prefix = NULL,
                                             main_title = NULL,
                                             sub_title = NULL,
                                             CD4 = FALSE,                 # <<< NEW
                                             cd4_var = "cd4_count_diag")  # <<< NEW
{
  trend_data <- get_trend_data(
    data = data,
    start_year = start_year,
    end_year = end_year,
    exposure_filter = exposure_filter,
    gender_filter = gender_filter,
    region_birth_filter = region_birth_filter,
    au_birth_filter = au_birth_filter,
    lang_home_filter = lang_home_filter,
    place_acq_filter = place_acq_filter,
    state_res_filter = state_res_filter,
    phn_res_filter = phn_res_filter,
    cd4_count_filter = cd4_count_filter,
    count_period = count_period,
    include_overseas = include_overseas,
    late_diag_filter = late_diag_filter,
    advanced_diag_filter = advanced_diag_filter,
    CD4 = CD4,
    cd4_var = cd4_var
  )
  
  fam <- if (CD4 && tolower(family_type) %in% c("poisson","nb")) "gaussian" else family_type
  overall_results <- get_trend(trend_data, family_type = fam)
  overall_trend_model <- overall_results$model
  data_years <- unique(floor(trend_data$time))
  
  selected <- find_best_segment_model(model = overall_trend_model,
                                      data_years = data_years,
                                      num_changes = num_changes,
                                      family_type = fam)
  selected_k <- if (!is.null(selected$selected_k)) selected$selected_k else 0
  
  if (selected_k == 0) {
    overall_plot <- overall_trend_plot(trend_data, overall_trend_model, count_period = count_period)
    segment_plot_result <- overall_plot + ggtitle(
      ifelse(is.null(main_title),
             paste0(if (CD4) "Average CD4 Trend (" else "HIV Notification Trend (", count_period, ")"),
             main_title),
      subtitle = ifelse(is.null(sub_title), "Selected model: 0 joinpoints (BIC)", sub_title)
    )
    
    if (tolower(overall_results$type_used) == "gaussian") {
      segment_results <- data.frame(
        Segment = "Overall",
        Slope_per_Year = overall_results$coeff,
        Lower_CI = overall_results$interval[1],
        Upper_CI = overall_results$interval[2],
        Interpretation = ifelse(overall_results$coeff > 0, "Increasing CD4", ifelse(overall_results$coeff < 0, "Decreasing CD4", "Stable")),
        Significant = ifelse(overall_results$interval[1] > 0 | overall_results$interval[2] < 0, "Yes", "No")
      )
    } else {
      segment_results <- data.frame(
        Segment = "Overall",
        Rate_Ratio = overall_results$coeff,
        Lower_CI = overall_results$interval[1],
        Upper_CI = overall_results$interval[2],
        Percent_Change_Per_Year = sprintf("%.1f%% (%.1f%% to %.1f%%)",
                                          (overall_results$coeff - 1) * 100,
                                          (overall_results$interval[1] - 1) * 100,
                                          (overall_results$interval[2] - 1) * 100),
        Interpretation = ifelse(overall_results$coeff > 1, "Increasing", ifelse(overall_results$coeff < 1, "Decreasing", "Stable")),
        Significant = ifelse(overall_results$interval[1] > 1 | overall_results$interval[2] < 1, "Yes", "No")
      )
    }
    
    change_points_table <- data.frame(Change_Point = integer(0),
                                      Estimate = numeric(0),
                                      Lower_CI = numeric(0),
                                      Upper_CI = numeric(0))
  } else {
    segment_plot_result <- segment_plot(
      trend_data = trend_data,
      trend_model = selected,                 
      count_period = count_period,
      overall_plot = show_overall,
      overall_model = overall_trend_model,
      segment_plot = show_segment,
      y_upper = y_upper,
      y_break = y_break
    ) + ggtitle(
      ifelse(is.null(main_title),
             paste0(if (CD4) "Average CD4 Segmented Trend (" else "HIV Notification Segmented Trend (", count_period, ")"),
             main_title),
      subtitle = ifelse(is.null(sub_title), paste0("Selected model: ", selected_k, " joinpoint(s) (BIC)"), sub_title)
    )
    
    # change points
    if (!is.null(selected$change_date) && nrow(selected$change_date) > 0) {
      to_month <- function(x) format(lubridate::round_date(lubridate::date_decimal(x), unit = "month"), "%Y-%m")
      change_points_table <- data.frame(
        Change_Point  = seq_len(nrow(selected$change_date)),
        Estimate      = selected$change_date$Estimate,
        Lower_CI      = selected$change_date$Lower_CI,
        Upper_CI      = selected$change_date$Upper_CI,
        Estimate_Month = to_month(selected$change_date$Estimate),
        Lower_CI_Month = to_month(selected$change_date$Lower_CI),
        Upper_CI_Month = to_month(selected$change_date$Upper_CI)
      )
    } else {
      change_points_table <- data.frame(
        Change_Point  = integer(0),
        Estimate      = numeric(0),
        Lower_CI      = numeric(0),
        Upper_CI      = numeric(0),
        Estimate_Month = character(0),
        Lower_CI_Month = character(0),
        Upper_CI_Month = character(0)
      )
    }
    
    # Segment results table: branch on Gaussian vs log-link families
    if (tolower(selected$model_type) == "gaussian") {
      segment_results <- data.frame(
        Segment = paste0("Segment ", seq_along(selected$slopes)),
        Slope_per_Year = selected$slopes,
        Lower_CI = selected$slopes_lower,
        Upper_CI = selected$slopes_upper
      )
      segment_results$Change_per_Year <- paste0(
        sprintf("%.1f", selected$slopes), " (",
        sprintf("%.1f", selected$slopes_lower), " to ",
        sprintf("%.1f", selected$slopes_upper), ")"
      )
      segment_results$Interpretation <- ifelse(
        selected$slopes > 0, "Increasing CD4",
        ifelse(selected$slopes < 0, "Decreasing CD4", "Stable")
      )
      segment_results$Significant <- ifelse(
        (selected$slopes_lower > 0) | (selected$slopes_upper < 0),
        "Yes", "No"
      )
    } else {
      segment_results <- data.frame(
        Segment = paste0("Segment ", seq_along(selected$slopes)),
        Rate_Ratio = selected$slopes,
        Lower_CI = selected$slopes_lower,
        Upper_CI = selected$slopes_upper
      )
      segment_results$Percent_Change_Per_Year <- paste0(
        sprintf("%.1f%%", (selected$slopes - 1) * 100),
        " (", sprintf("%.1f%%", (selected$slopes_lower - 1) * 100),
        " to ", sprintf("%.1f%%", (selected$slopes_upper - 1) * 100), ")"
      )
      segment_results$Interpretation <- ifelse(
        selected$slopes > 1, "Increasing",
        ifelse(selected$slopes < 1, "Decreasing", "Stable")
      )
      segment_results$Significant <- ifelse(
        (selected$slopes_lower > 1) | (selected$slopes_upper < 1),
        "Yes", "No"
      )
    }
    
    wald_results <- if (!is.null(selected$wald_slope_change)) selected$wald_slope_change else
      data.frame(Joinpoint = integer(0), Term = character(0),
                 Estimate = numeric(0), SE = numeric(0), z = numeric(0),
                 p_value = numeric(0), JP_Time = numeric(0),
                 JP_Time_Month = character(0), Exp_Estimate = numeric(0))
    wald_results$Significant <- ifelse(wald_results$p_value < 0.05, "Yes", "No")
  }
  
  bic_table <- if (!is.null(selected$all_bic)) selected$all_bic else data.frame(Joinpoints = 0, BIC = BIC(overall_trend_model))
  
  if (selected_k == 0) {
    # No joinpoints: just report overall slope/CI (Gaussian expected for CD4)
    if (tolower(overall_results$type_used) == "gaussian") {
      overall_summary <- data.frame(
        Population = "Overall",
        Average_Slope_CD4_per_Year = unname(overall_results$coeff),
        Lower_CI = unname(overall_results$interval[1]),
        Upper_CI = unname(overall_results$interval[2]),
        stringsAsFactors = FALSE
      )
    } else {
      # Log-link case (for counts): standard APC/AAPC style
      overall_summary <- data.frame(
        Population = "Overall",
        Average_APC = unname(overall_results$coeff),
        Lower_CI = unname(overall_results$interval[1]),
        Upper_CI = unname(overall_results$interval[2]),
        Average_Percent_Change = sprintf(
          "%.1f%% (%.1f%% to %.1f%%)",
          (overall_results$coeff - 1) * 100,
          (overall_results$interval[1] - 1) * 100,
          (overall_results$interval[2] - 1) * 100
        ),
        stringsAsFactors = FALSE
      )
    }
  } else {
    # Joinpoints selected
    if (tolower(selected$model_type) == "gaussian") {
      # >>> Correct CI for time-weighted mean slope (no heuristic bounds)
      base_df <- as.data.frame(stats::model.frame(overall_trend_model))
      wg <- weighted_gaussian_slope_CI(selected, base_df)
      
      overall_summary <- data.frame(
        Population = "Overall",
        Weighted_Slope_CD4_per_Year = wg$weighted_slope,
        Lower_CI = wg$ci95[1],
        Upper_CI = wg$ci95[2],
        stringsAsFactors = FALSE
      )
    } else {
      # Log-link case: AAPC via weighted average of log-slopes
      overall_summary <- data.frame(
        Population = "Overall",
        Average_APC = selected$aapc$est,
        Lower_CI = selected$aapc$lo,
        Upper_CI = selected$aapc$up,
        Average_Percent_Change = sprintf(
          "%.1f%% (%.1f%% to %.1f%%)",
          (selected$aapc$est - 1) * 100,
          (selected$aapc$lo  - 1) * 100,
          (selected$aapc$up  - 1) * 100
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # --- Add CD4 quantile overlay when CD4 mode is on ---
  cd4_quantile_data <- NULL
  cd4_quantile_plot_obj <- NULL
  
  if (isTRUE(CD4)) {
    cd4_quantile_data <- get_cd4_quantile_trend_data(
      data = data,
      start_year = start_year,
      end_year = end_year,
      exposure_filter = exposure_filter,
      gender_filter = gender_filter,
      region_birth_filter = region_birth_filter,
      au_birth_filter = au_birth_filter,
      lang_home_filter = lang_home_filter,
      place_acq_filter = place_acq_filter,
      state_res_filter = state_res_filter,
      phn_res_filter = phn_res_filter,
      include_overseas = include_overseas,
      count_period = count_period,
      cd4_var = cd4_var,
      probs = c(0.05, 0.50, 0.95)
    )
    
    # pass 'selected' so joinpoints appear as dashed lines if available
    cd4_quantile_plot_obj <- cd4_quantile_plot(
      cd4_quant_data = cd4_quantile_data,
      count_period   = count_period,
      segment_model  = if (exists("selected")) selected else NULL,
      y_upper        = y_upper,
      y_break        = y_break
    )
  }
  
  
  
  list(
    trend_data = trend_data,
    overall_results = overall_results,
    selected_k = selected_k,
    segmented_model = if (selected_k == 0) NULL else selected,
    segment_plot = segment_plot_result,
    change_points = change_points_table,
    segment_results = segment_results,
    overall_summary = overall_summary,
    bic_candidates = bic_table,
    wald_adjacent_slopes = if (selected_k == 0) NULL else wald_results,
    cd4_quantile_data = cd4_quantile_data,
    cd4_quantile_plot = cd4_quantile_plot_obj
  )
}
