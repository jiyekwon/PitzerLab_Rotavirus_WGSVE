
# Function:  --------------------------------------------------------------------# 
# fit VE model and extract summary
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
  
  return(ve_summary)
}