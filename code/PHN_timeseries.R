# =============================================================================
# PHN_timeseries.R  —  sourced by 1_timeseriesAnalysis.Rmd
# Primary Health Network (PHN) time-series helpers: notification rates per
# 100,000 population, negative-binomial rate trends (annual percent change),
# state-coloured rate heatmaps, and geographic inequality (population-weighted
# Gini + Lorenz curves) across PHNs.
# =============================================================================

# trend for national for fed and pdo, using negative binomial regression with log(pop) offset
fit_national_rate_trend <- function(df) {
  # Negative binomial model on counts with log(pop) offset
  mod   <- MASS::glm.nb(
    hiv_cases ~ year + offset(log(total_pop)),
    data = df
  )
  mtype <- "negative binomial"
  
  co <- summary(mod)$coefficients
  beta <- co["year", "Estimate"]
  se   <- co["year", "Std. Error"]
  p    <- co["year", grep("^Pr", colnames(co))]
  
  # Convert to rate ratio and APC (% change per year)
  RR    <- exp(beta)
  RR_lo <- exp(beta - 1.96 * se)
  RR_hi <- exp(beta + 1.96 * se)
  
  APC    <- (RR    - 1) * 100
  APC_lo <- (RR_lo - 1) * 100
  APC_hi <- (RR_hi - 1) * 100
  
  list(
    model   = mod,
    summary = tibble::tibble(
      model_type            = mtype,
      beta                  = beta,
      se                    = se,
      p_value               = p,
      RR_per_year           = RR,
      RR_lower_95CI         = RR_lo,
      RR_upper_95CI         = RR_hi,
      APC_percent_per_year  = APC,
      APC_lower_95CI        = APC_lo,
      APC_upper_95CI        = APC_hi
    )
  )
}

# Combine results into a tidy table
tidy_nat_rate_trends <- function(...) {
  dots <- list(...)
  
  # names like "PDO", "FED", "LATE" from the arguments
  grp_names <- names(dots)
  if (is.null(grp_names) || any(grp_names == "")) {
    grp_names <- paste0("group_", seq_along(dots))
  }
  
  purrr::imap_dfr(dots, ~ {
    out <- .x$summary   # each is the $summary tibble from fit_national_rate_trend()
    
    out %>%
      dplyr::mutate(
        group = .y,
        APC_label = sprintf(
          "%.1f (%.1f to %.1f)",
          APC_percent_per_year,
          APC_lower_95CI,
          APC_upper_95CI
        ),
        p_value = round(p_value, 3)
      ) %>%
      dplyr::select(
        group,
        model_type,
        RR_per_year,
        RR_lower_95CI,
        RR_upper_95CI,
        APC_percent_per_year,
        APC_lower_95CI,
        APC_upper_95CI,
        APC_label,
        p_value
      )
  })
}

# PHN case rates: cases × PHN × year, joined to PHN population 
compute_phn_case_rates <- function(hiv_data, previ_flag, phn_pop_yearly,
                                   time_period = 12,
                                   start_year = 2001, end_year = 2024) {
  cases <- hiv_data %>%
    dplyr::filter(year_diag >= start_year, year_diag <= end_year,
                  !is.na(phn_res), !is.na(year_diag),
                  previ_diag_overseas == previ_flag) %>%
    dplyr::group_by(phn_res, phn_name, year_diag) %>%
    dplyr::summarise(hiv_cases = dplyr::n(), .groups = "drop") %>%
    dplyr::rename(PHN_CODE = phn_res, PHN_NAME = phn_name, year = year_diag)
  
  phn_pop_yearly %>%
    dplyr::left_join(cases, by = c("PHN_CODE", "PHN_NAME", "year")) %>%
    dplyr::mutate(
      hiv_cases = ifelse(is.na(hiv_cases), 0, hiv_cases),
      case_rate_per_100k = (hiv_cases / total_pop) * 100000 / time_period
    ) %>%
    dplyr::filter(total_pop > 0)
}

# PHN-level trend: baseline rate + NB rate trend + APC + category
compute_phn_trend <- function(hiv_case_rates, start_year = 2001, end_year = 2024) {
  
  phn_baseline_rate <- hiv_case_rates %>%
    dplyr::filter(year %in% start_year:end_year) %>%
    dplyr::group_by(PHN_CODE, PHN_NAME) %>%
    dplyr::summarise(
      baseline_year = min(year, na.rm = TRUE),
      baseline_rate = case_rate_per_100k[which.min(year)],
      .groups = "drop"
    )
  
  na_row <- function() tibble::tibble(
    model_type = NA_character_, beta = NA_real_, se = NA_real_,
    p_value = NA_real_, RR_per_year = NA_real_,
    APC = NA_real_, lower_APC = NA_real_, upper_APC = NA_real_
  )
  
  phn_case_trend <- hiv_case_rates %>%
    dplyr::filter(year %in% start_year:end_year) %>%
    dplyr::group_by(PHN_CODE, PHN_NAME) %>%
    dplyr::group_modify(~ {
      df <- .x %>%
        dplyr::filter(!is.na(hiv_cases), !is.na(total_pop), total_pop > 0)
      
      if (nrow(df) == 0 || all(df$hiv_cases == 0, na.rm = TRUE)) return(na_row())
      
      mod <- tryCatch(
        MASS::glm.nb(hiv_cases ~ year + offset(log(total_pop)), data = df),
        error = function(e) NULL
      )
      if (is.null(mod)) return(na_row())
      
      co   <- summary(mod)$coefficients
      beta <- co["year", "Estimate"]
      se   <- co["year", "Std. Error"]
      p    <- co["year", grep("^Pr", colnames(co))]
      
      RR    <- exp(beta); RR_lo <- exp(beta - 1.96*se); RR_hi <- exp(beta + 1.96*se)
      
      tibble::tibble(
        model_type = "negative binomial",
        beta = beta, se = se, p_value = p, RR_per_year = RR,
        APC = (RR - 1)*100, lower_APC = (RR_lo - 1)*100, upper_APC = (RR_hi - 1)*100
      )
    }) %>%
    dplyr::ungroup()
  
  phn_case_trend %>%
    dplyr::left_join(phn_baseline_rate, by = c("PHN_CODE", "PHN_NAME")) %>%
    dplyr::mutate(
      annual_percent_change = round(APC, 1),
      apc_lower = round(lower_APC, 1),
      apc_upper = round(upper_APC, 1),
      ci_95 = sprintf("%.1f to %.1f", apc_lower, apc_upper),
      apc_with_ci = sprintf("%.1f (%.1f to %.1f)", annual_percent_change, apc_lower, apc_upper),
      p_value = round(p_value, 3),
      trend_category = dplyr::case_when(
        !is.na(p_value) & p_value < 0.05 & APC > 0 ~ "Increasing",
        !is.na(p_value) & p_value < 0.05 & APC < 0 ~ "Decreasing",
        is.na(p_value)                             ~ "Insufficient data",
        TRUE                                       ~ "Stable"
      )
    ) %>%
    dplyr::transmute(
      PHN_CODE, PHN_NAME, baseline_year, baseline_rate,
      annual_percent_change, ci_95, apc_with_ci,
      p_value, trend_category, model_type
    )
}

# PHN heatmap 
make_phn_heatmap <- function(hiv_case_rates, sub_title) {
  
  df <- hiv_case_rates %>%
    dplyr::select(PHN_NAME, year, case_rate = case_rate_per_100k) %>%
    dplyr::left_join(
      phn_map %>% dplyr::select(PHN_NAME, PHN_CODE, STATE, code_num),
      by = "PHN_NAME"
    )
  
  state_order <- c("NSW", "VIC", "QLD", "SA", "WA", "TAS", "NT", "ACT")
  state_legend_df <- tibble::tibble(STATE = factor(state_order, levels = state_order))
  
  ggplot(df, aes(x = year,
                 y = factor(PHN_NAME, levels = rev(ordered_names)),
                 fill = case_rate)) +
    geom_tile(color = "white", size = 0.2) +
    scale_fill_gradientn(
      colors = common_colors, breaks = common_breaks, labels = common_labels,
      limits = c(min(common_breaks), max(common_breaks)),
      name = "Notifications\nper 100,000"
    ) +
    scale_y_discrete(limits = rev(ordered_names),
                     labels = y_labels[rev(ordered_names)]) +
    scale_x_continuous(breaks = seq(2001, 2024, by = 5)) +
    geom_point(data = state_legend_df,
               aes(x = 2001, y = ordered_names[1], color = STATE),
               inherit.aes = FALSE, alpha = 0, size = 6) +
    scale_color_manual(values = state_cols[state_order], name = "State") +
    guides(
      fill  = guide_colorbar(order = 1, barwidth = unit(1.2, "cm"),
                             barheight = unit(8, "cm"),
                             title.position = "top", title.hjust = 0),
      color = guide_legend(order = 2, override.aes = list(alpha = 1, size = 18),
                           direction = "vertical", title.position = "top")
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
      axis.text.y = ggtext::element_markdown(size = 12),
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 18),
      legend.box = "vertical", legend.box.just = "left",
      legend.spacing.y = unit(4, "pt"),
      legend.direction = "vertical", legend.position = "right",
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18),
      legend.justification = "top",
      axis.title = element_text(size = 18)
    ) +
    labs(title = "", subtitle = sub_title,
         x = "Year", y = "Primary Health Network")
}

# Weighted Gini computation
weighted_gini <- function(x, w) {
  # Remove NAs
  ok <- !is.na(x) & !is.na(w)
  x <- x[ok]
  w <- w[ok]
  
  # Order by x ascending (low to high)
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  
  # Cumulative weighted proportions - include (0,0) point
  cw <- c(0, cumsum(w) / sum(w))
  px <- c(0, cumsum(w * x) / sum(w * x))
  
  # Numerical integration (trapezoid rule)
  # Area under Lorenz curve
  area <- sum((px[-1] + px[-length(px)]) / 2 * diff(cw))
  
  # Gini coefficient
  G <- 1 - 2 * area
  
  return(G)
}

plot_phn_lorenz <- function(hiv_case_rates_fed, phn_map, year_selected = 2001,
                            label_left = c("Perth North", "South Western Sydney", "Adelaide"),
                            label_right = c("Central and Eastern Sydney"),label_lefttop = c(),sub_title = "") {
  
  # Prepare data
  data_year <- hiv_case_rates_fed %>%
    filter(year == year_selected, total_pop > 0) %>%
    left_join(phn_map %>% dplyr::select(PHN_CODE, STATE), by = "PHN_CODE") %>%
    arrange(case_rate_per_100k) %>%
    mutate(
      cum_pop = cumsum(total_pop),
      cum_cases = cumsum(hiv_cases),
      prop_pop = cum_pop / sum(total_pop, na.rm = TRUE),
      prop_cases = cum_cases / sum(hiv_cases, na.rm = TRUE)
    )
  
  # Compute Gini
  gini_val <- weighted_gini(data_year$case_rate_per_100k, data_year$total_pop)
  
  # Prepare labels
  lefttop_labels <- data_year %>% filter(PHN_NAME %in% label_lefttop) %>%
    mutate(
      x_text = prop_pop - 0.03,
      y_text = prop_cases + 0.03
    )
  left_labels <- data_year %>% filter(PHN_NAME %in% label_left) %>%
    mutate(
      x_text = prop_pop - 0.03,
      y_text = prop_cases 
    )  
  right_labels <- data_year %>% filter(PHN_NAME %in% label_right) %>%
    mutate(
      x_text = prop_pop + 0.03,
      y_text = prop_cases - 0.03
    )
  
  #
  p <- ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey70") +
    
    # Colored step segments by state
    geom_segment(data = data_year,
                 aes(x = lag(prop_pop, default = 0),
                     y = lag(prop_cases, default = 0),
                     xend = prop_pop,
                     yend = lag(prop_cases, default = 0),
                     color = STATE), linewidth = 0.5) +
    geom_segment(data = data_year,
                 aes(x = prop_pop,
                     y = lag(prop_cases, default = 0),
                     xend = prop_pop,
                     yend = prop_cases,
                     color = STATE), linewidth = 0.5) +
    
    scale_color_manual(values = state_cols) +
    scale_x_continuous(limits = c(0, 1.05),
                       breaks = seq(0, 1, 0.25),
                       labels = c("0·00", "0·25", "0·50",  "0·75", "1·00"),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.05),
                       breaks = seq(0, 1, 0.25),
                       labels = c("0·00", "0·25", "0·50",  "0·75", "1·00"),
                       expand = c(0, 0)) +
    labs(
      title = paste0(""),
      subtitle = paste0(sub_title),
      x = "Cumulative proportion of population",
      y = "Cumulative proportion of notifications",
      color = "State/Territory"
    ) +
    
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      plot.margin = margin(t = 10, r = 50, b = 10, l = 10),
      plot.title = element_text(face = "bold", hjust = 0, size = 12),
      plot.subtitle = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10)
    ) +
    coord_fixed(ratio = 1, clip = "off")
  
  # Add label connectors and text
  if (nrow(left_labels)) {
    p <- p +
      geom_segment(
        data = left_labels,
        aes(x = x_text, y = y_text, xend = prop_pop, yend = prop_cases),
        color = "black"
      ) +
      geom_label(
        data = left_labels,
        aes(x = x_text, y = y_text, label = PHN_NAME),
        hjust = 1, vjust = 0.5, size = 3,
        fill = "white",        # White background behind label text
        label.size = 0         # No border line
      )
  }
  
  if (nrow(lefttop_labels)) {
    p <- p +
      geom_segment(
        data = lefttop_labels,
        aes(x = x_text, y = y_text, xend = prop_pop, yend = prop_cases),
        color = "black"
      ) +
      geom_label(
        data = lefttop_labels,
        aes(x = x_text, y = y_text, label = PHN_NAME),
        hjust = 1, vjust = 0, size = 3,
        fill = "white",        # White background behind label text
        label.size = 0         # No border line
      )
  }  
  
  if (nrow(right_labels)) {
    p <- p +
      geom_segment(
        data = right_labels,
        aes(x = x_text, y = y_text, xend = prop_pop, yend = prop_cases),
        color = "black"
      ) +
      geom_label(
        data = right_labels,
        aes(x = x_text, y = y_text, label = PHN_NAME),
        hjust = 0, vjust = 1, size = 3,
        fill = "white",        # White background behind text
        label.size = 0         # Removes border box outline
      )
  }
  
  return(list(plot = p, gini = gini_val, data = data_year))
}


phn_to_state <- list(
  "PHN108" = "NSW", "PHN110" = "NSW", "PHN107" = "NSW", "PHN103" = "NSW", 
  "PHN109" = "NSW", "PHN105" = "NSW", "PHN104" = "NSW", "PHN106" = "NSW",
  "PHN101" = "NSW", "PHN102" = "NSW",
  "PHN202" = "VIC", "PHN204" = "VIC", "PHN203" = "VIC", "PHN201" = "VIC", 
  "PHN206" = "VIC", "PHN205" = "VIC",
  "PHN302" = "QLD", "PHN304" = "QLD", "PHN307" = "QLD", "PHN305" = "QLD", 
  "PHN306" = "QLD", "PHN303" = "QLD", "PHN301" = "QLD",
  "PHN402" = "SA", "PHN401" = "SA",
  "PHN503" = "WA", "PHN502" = "WA", "PHN501" = "WA",
  "PHN601" = "TAS",
  "PHN801" = "ACT",
  "PHN701" = "NT"
)

# Create PHN to name mapping
phn_to_name <- list(
  "PHN108" = "Hunter New England and Central Coast", 
  "PHN110" = "Murrumbidgee",
  "PHN107" = "Western NSW",
  "PHN103" = "Western Sydney",
  "PHN109" = "North Coast",
  "PHN105" = "South Western Sydney",
  "PHN104" = "Nepean Blue Mountains",
  "PHN106" = "South Eastern NSW", 
  "PHN101" = "Central and Eastern Sydney",
  "PHN102" = "Northern Sydney",
  "PHN202" = "Eastern Melbourne",
  "PHN204" = "Gippsland",
  "PHN203" = "South Eastern Melbourne",
  "PHN201" = "North Western Melbourne",
  "PHN206" = "Western Victoria",
  "PHN205" = "Murray",
  "PHN302" = "Brisbane South",
  "PHN304" = "Darling Downs and West Moreton",
  "PHN307" = "Northern Queensland",
  "PHN305" = "Western Queensland",
  "PHN306" = "Central Queensland, Wide Bay, Sunshine Coast",
  "PHN303" = "Gold Coast",
  "PHN301" = "Brisbane North",
  "PHN402" = "Country SA",
  "PHN401" = "Adelaide",
  "PHN503" = "Country WA",
  "PHN502" = "Perth South",
  "PHN501" = "Perth North",
  "PHN601" = "Tasmania",
  "PHN801" = "Australian Capital Territory",
  "PHN701" = "Northern Territory"
)

state_vec <- unlist(phn_to_state)

phn_map <- tibble::tibble(
  PHN_CODE = names(phn_to_name),
  PHN_NAME = unname(unlist(phn_to_name))
) |>
  dplyr::mutate(
    STATE    = state_vec[PHN_CODE],      
    code_num = as.integer(sub("^PHN", "", PHN_CODE))
  )

# Order PHNs by numeric code (101, 102, ...)
ordered_names <- phn_map |>
  dplyr::arrange(code_num) |>
  dplyr::pull(PHN_NAME)

#Colors to use for state-coded y-axis labels
state_cols <- c(
  "NSW" = "#1f77b4",
  "VIC" = "#ff7f0e",
  "QLD" = "#2ca02c",
  "SA"  = "#d62728",
  "WA"  = "#9467bd",
  "TAS" = "#8c564b",
  "ACT" = "#e377c2",
  "NT"  = "#7f7f7f"
)

# Build HTML-styled labels for ggtext::element_markdown()
label_df <- phn_map |>
  dplyr::mutate(
    label_html = sprintf(
      "<span style='color:%s'>%s</span>",
      state_cols[STATE], PHN_NAME
    )
  )

# Named vector: names = PHN_NAME, values = colored HTML labels
y_labels <- stats::setNames(label_df$label_html, label_df$PHN_NAME)

