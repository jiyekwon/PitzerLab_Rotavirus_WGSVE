
# Function:  -------------------------------------------------------------------- 
##  fit VE model and extract summary
get_ve_summary <- function(data, label, adjusted = FALSE) {
  
  if (adjusted) {
    mod <- glm(INFECTION ~ VaccineBinary + ageYear + careLevel + colYear2, family = binomial, data = data)
  } else {
    mod <- glm(INFECTION ~ VaccineBinary, family = binomial, data = data)
  }
  conf <- confint(mod)
  ve <- 1 - exp(coef(mod)["VaccineBinary"])
  ci <- 1 - exp(conf["VaccineBinary", ])
  tibble(
    Group = label,
    Model = ifelse(adjusted, "Adjusted", "Unadjusted"),
    VE = round(ve, 4)*100,
    CI_lower = round(ci[2], 4)*100,
    CI_upper = round(ci[1], 4)*100
  )
}

run_ve_by_group <- function(data, comparison_label, group_vals, outcome_var = "GxPx") {
  purrr::map_dfr(group_vals, function(g) {
    dat <- if (outcome_var == "GxPx") {
      data %>% filter(GxPx == g | is.na(GxPx))
    } else if (outcome_var == "protection") {
      data %>% filter(protection == g | is.na(protection))
    } else if (outcome_var == "strainSpecific") {
      data %>% filter(strainSpecific == g | strainSpecific == "Control")
    } else if (outcome_var == "clusterAssignment.rv5") {
      data %>% filter(clusterAssignment.rv5 == g | clusterAssignment.rv5 == "Control")
    } else {
      stop("Unknown outcome_var")
    }
    
    bind_rows(
      get_ve_summary(dat, g, adjusted = FALSE),
      get_ve_summary(dat, g, adjusted = TRUE)
    ) %>% mutate(Comparison = comparison_label)
  })
}


run_multinom_ve_analysis <- function(
    data,
    comparison_label, 
    outcome_type = c("INFECTION","strainSpecific", "GxPx", "protection", "clusterAssignment.rv5"),
    covariates = c("ageYear", "careLevel", "colYear2"),
    predictor = "VaccineBinary"
) {
  outcome_type <- match.arg(outcome_type)
  
  # Define case_group based on the outcome_type
  data <- data %>%
    mutate(case_group = dplyr::case_when(
      outcome_type == "strainSpecific" & strainSpecific == "Control" ~ "Control",
      outcome_type == "strainSpecific" & strainSpecific == "StrainI" ~ "j=0",
      outcome_type == "strainSpecific" & strainSpecific == "StrainS" ~ "j=1",
      
      outcome_type == "clusterAssignment.rv5" & clusterAssignment.rv5 == "Control" ~ "Control",
      outcome_type == "clusterAssignment.rv5" & clusterAssignment.rv5 == "Cluster1" ~ "j=0",
      outcome_type == "clusterAssignment.rv5" & clusterAssignment.rv5 == "Cluster2" ~ "j=1",
      
      outcome_type == "GxPx" & is.na(GxPx) ~ "Control",
      outcome_type == "GxPx" & GxPx == "G1P[8]" ~ "G1P[8]",
      outcome_type == "GxPx" & GxPx == "G9P[8]" ~ "G9P[8]",
      outcome_type == "GxPx" & GxPx == "G12P[8]" ~ "G12P[8]",
      outcome_type == "GxPx" & GxPx == "G2P[4]" ~ "G2P[4]",
      outcome_type == "GxPx" ~ "other",
      
      outcome_type == "protection" & is.na(protection) ~ "Control",
      outcome_type == "protection" & protection == "homotypic" ~ "homotypic",
      outcome_type == "protection" & protection == "partially heterotypic" ~ "partially heterotypic",
      outcome_type == "protection" & protection == "fully heterotypic" ~ "fully heterotypic"
    )) %>%
    mutate(case_group = factor(case_group, levels = c(
      "Control", "AnyInfection", "j=0", "j=1", "G1P[8]", "G9P[8]", "G12P[8]", "G2P[4]","other",
      "homotypic", "partially heterotypic", "fully heterotypic"
    )))
  
  # Fit multinomial model
  formula_str <- paste("case_group ~", paste(c(predictor, covariates), collapse = " + "))
  mod_multi <- multinom(as.formula(formula_str), data = data)
  
  # Extract VE estimates
  coef_matrix <- coef(mod_multi)
  vcov_matrix <- vcov(mod_multi)
  
  ve_summary <- map_dfr(rownames(coef_matrix), function(outcome) {
    term <- paste0(outcome, ":", predictor)
    beta <- coef_matrix[outcome, predictor]
    se <- sqrt(diag(vcov_matrix)[term])
    
    ve <- 1 - exp(beta)
    ci_lower <- 1 - exp(beta + 1.96 * se)
    ci_upper <- 1 - exp(beta - 1.96 * se)
    
    tibble(
      OutcomeType = outcome_type,
      Outcome = outcome,
      Comparison= comparison_label,
      VE = round(ve * 100, 2),
      CI_lower = round(ci_lower * 100, 2),
      CI_upper = round(ci_upper * 100, 2)
    )
  })
  
  # Add sample sizes
  n_table <- data %>%
    group_by(case_group, !!sym(predictor)) %>%
    tally() %>%
    pivot_wider(names_from = !!sym(predictor), values_from = n, values_fill = 0) %>%
    rename( Outcome = case_group, n_unvax = `0`, n_vax = `1`)
  
  ve_summary <- ve_summary %>%
    left_join(n_table, by = "Outcome")
  
  return(list(ve_summary, mod_multi))
}



##  gd by location --- 
get_city_summary <- function(data, distance_col, vaxx_col) {
  data %>%
    group_by(studySite) %>%
    summarise(
      coverage = mean({{vaxx_col}}, na.rm = TRUE),
      mean_distance = mean({{distance_col}}, na.rm = TRUE),
      sd_distance = sd({{distance_col}}, na.rm = TRUE),
      n_samples = n(),
      se_distance = sd_distance / sqrt(n_samples),
      .groups = "drop"
    ) %>%
    mutate(
      ci_lower = mean_distance - 1.96 * se_distance,
      ci_upper = mean_distance + 1.96 * se_distance
    )
}



## sieve analysis functions --- 
cutBin <- function(col){
  as.character(cut(col, breaks = breaks, include.lowest = TRUE))
}

cutBin.f <- function(col){
  cut(col, breaks = breaks, include.lowest = TRUE)
}


run_sieve_analysis <- function(df, perc_col, vax_to_exclude, breaks_vec, label) {
  # Filter and prepare
  df <- df %>% filter(VaxxStatEver_4 != vax_to_exclude) %>%
    mutate(bin = cut(.data[[perc_col]], breaks = breaks_vec, include.lowest = TRUE))
  
  # Define unique bin levels
  df$bin <- factor(df$bin, levels = levels(df$bin))
  
  # Print contingency table
  cat(paste0("\n## ", label, " Bin Distribution Table ##\n"))
  print(table(df$VaxxStatEver_4, df$bin))
  
  # Run unadjusted model
  model_unadj <- multinom(bin ~ VaxxStatEver_4, data = df)
  
  # Run adjusted models
  model_c     <- multinom(bin ~ VaxxStatEver_4 + ageYear + careLevel + colYear2, data = df)
  model_2     <- multinom(bin ~ VaxxStatEver_4 + colYear2, data = df)
  model_3     <- multinom(bin ~ VaxxStatEver_4 + ageYear + colYear2, data = df)
  model_4     <- multinom(bin ~ VaxxStatEver_4 + careLevel + colYear2, data = df)
  
  # AIC comparison
  model_list <- list(unadjusted = model_unadj, adj_full = model_c, adj_colYear2 = model_2,
                     adj_age_colYear2 = model_3, adj_care_colYear2 = model_4)
  aic_values <- sapply(model_list, AIC)
  
  # Summary table of odds ratios
  summary_table <- function(model) {
    ors <- odds.ratio(model)
    data.frame(OR = round(ors[,1], 2),
               Lower = round(ors[,2], 2),
               Upper = round(ors[,3], 2))
  }
  
  or_unadj <- summary_table(model_unadj)
  or_adj2  <- summary_table(model_2)
  or_adjC  <- summary_table(model_c)
  
  list(
    aic = aic_values,
    OR_unadjusted = or_unadj,
    OR_adj_colYear2 = or_adj2,
    OR_adj_full = or_adjC
  )
}


