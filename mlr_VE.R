library(ggplot2)
library(ggpubr)
library(patchwork)
library(broom)
library(dplyr)
library(car)
library(nnet)
library(purrr)
library(purrr)
source("functions.R")
###############################################################################
# VE comparison by Gxpx, protection types 
###############################################################################
# Define genotypes and protection types of interest
vaccine_groups <- c("Overall", "RV5-specific", "RV1-specific")
outcome_types <- c("GxPx", "Protection", "Distance")
genotypes <- c("G1P[8]", "G12P[8]", "G2P[4]", "G9P[8]")
genotypes.comp<- c("G1P[8]", "G2P[8]", "G3P[8]", "G9P[8]", "G12P[8]", "G1P[4]", "G2P[4]", "G3P[6]")
protection_types <- c("homotypic", "partially heterotypic", "fully heterotypic")


# data cleaning ---------------------------------------------------------------#
ve.all.data <- ve.all.data %>%
  mutate(protection = case_when(
    is.na(GxPx) ~ NA_character_,
    GxPx == "G1P[8]" ~ "homotypic",
    grepl("G1", GxPx) | grepl("P\\[8\\]", GxPx) ~ "partially heterotypic",
    TRUE ~ "fully heterotypic"
  )) %>%
  mutate(protection = factor(protection, levels = protection_types),
         GxPx = factor(GxPx, levels = genotypes.comp))



# Subset:  --------------------------------------------------------------------#
# Filter RV1 only -- RV1-specific analysis
ve.all.data <- ve.all.data %>% 
  mutate(strainSpecific = case_when(strainSpecific == "StrainS" ~ "StrainS",
                                    strainSpecific == "StrainI" ~ "StrainI",
                                    strainSpecific == "control" ~ "Control"))

rv1.only <- ve.all.data %>% filter(VaccineType != "RV5-vaccinated")
rv5.only <- ve.all.data %>% filter(VaccineType!= "RV1-vaccinated")
overall.anyvax <- ve.all.data %>% filter(VaccineType %in% c("RV1-vaccinated", "RV5-vaccinated", "Unvaccinated"))



# Output:  --------------------------------------------------------------------# 
# Overall VE to any infection (not genotype/protection/distance stratified)
ve_overall_any <- bind_rows(
  get_ve_summary(rv1.only, "RV1-specific", adjusted = FALSE),
  get_ve_summary(rv1.only, "RV1-specific", adjusted = TRUE),
  get_ve_summary(rv5.only, "RV5-specific", adjusted = FALSE),
  get_ve_summary(rv5.only, "RV5-specific", adjusted = TRUE)
) %>% mutate(Outcome = "Any Infection", OutcomeType = "Overall", Comparison = label)
ve_overall_any$Comparison <- c(rep("RV1-specific",2),rep("RV5-specific",2))

ve_oa_multi <- get_ve_summary(rv1.only, "RV1-specific", adjusted = TRUE)%>% mutate(Comparison = Group, OutcomeType= "Overall", Outcome = "Any infection") %>% select(-c(Model, Group))
ve_dist_multi <- run_multinom_ve_analysis(data = rv1.only, comparison_label = "RV1-specific",  outcome_type = "strainSpecific")
ve_protection_multi <- run_multinom_ve_analysis(data = rv1.only, comparison_label = "RV1-specific", outcome_type = "protection")
ve_gxpx_multi <- run_multinom_ve_analysis(data = rv1.only, comparison_label = "RV1-specific", outcome_type = "GxPx")

ve_oa_multi5 <- get_ve_summary(rv5.only, "RV5-specific", adjusted = TRUE) %>% mutate(Comparison = Group, OutcomeType= "Overall", Outcome = "Any infection") %>% select(-c(Model, Group))
ve_dist_multi5 <- run_multinom_ve_analysis(data = rv5.only, comparison_label = "RV5-specific",  outcome_type = "clusterAssignment.rv5")
ve_protection_multi5 <- run_multinom_ve_analysis(data = rv5.only, comparison_label = "RV5-specific", outcome_type = "protection")
ve_gxpx_multi5 <- run_multinom_ve_analysis(data = rv5.only, comparison_label = "RV5-specific", outcome_type = "GxPx")

ve_combined.multi <- bind_rows(
  ve_oa_multi,
  ve_dist_multi,
  ve_dist_multi5,
  ve_gxpx_multi,
  
  ve_oa_multi5,
  ve_gxpx_multi5,
  ve_protection_multi,
  ve_protection_multi5
)

ve_combined.multi <- ve_combined.multi %>%
  mutate(OutcomeType = factor(OutcomeType,
                              levels = c("Overall",  "strainSpecific", "clusterAssignment.rv5", "GxPx", "protection"),
                              labels = c("Overall",  "Distance_RV1", "Distance_RV5", "Genotype", "RV1-typic")),
         VE_display = sprintf("%.1f (%.1f, %.1f)", VE, CI_lower, CI_upper))
ve_combined.multi <- ve_combined.multi %>% filter(Outcome != "other") 

# Figures:  --------------------------------------------------------------------# 
ve_rv1_gxpx <- ve_combined.multi %>%
  filter(Comparison == "RV1-specific",
         # Model == "Adjusted",
         OutcomeType != "Distance_RV5")%>%
  mutate(Outcome = factor(Outcome, levels = rev(unique(Outcome))))
gg_rv1 <- ggplot(ve_rv1_gxpx, aes(x = Outcome, y = VE, ymin = CI_lower, ymax = CI_upper)) +
  geom_hline(yintercept = c(0, 25, 50, 75, 100), color = "grey70", linetype = "dotted") +
  geom_pointrange(color = "#54278f", size = 0.8, alpha = 0.4) +
  scale_y_continuous(
    limits = c(-30, 100),
    breaks = seq(0, 100, 25),
    expand = expansion(mult = c(0, 0.05))
  )+
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "",
    y = "Vaccine Effectiveness",
  ) +
  theme_pubr(base_size = 14)+
  facet_grid(OutcomeType ~ Comparison, scales = "free_y", space = "free_y", switch = "y") +
  coord_flip()+
  theme(panel.spacing = unit(1.5, "lines") )



ve_rv5_gxpx <- ve_combined.multi %>%
  filter(Comparison == "RV5-specific",
         # Model == "Adjusted",
         OutcomeType != "Distance_RV1")%>%
  mutate(Outcome = factor(Outcome, levels = rev(unique(Outcome))))
gg_rv5 <- ggplot(ve_rv5_gxpx, aes(x = Outcome, y = VE, ymin = CI_lower, ymax = CI_upper)) +
  geom_hline(yintercept = c(0, 25, 50, 75, 100), color = "grey70", linetype = "dotted") +
  geom_pointrange(color = "maroon", size = 0.8, alpha = 0.4) +
  scale_y_continuous(
    limits = c(-30, 100),
    breaks = seq(0, 100, 25),
    expand = expansion(mult = c(0, 0.05))
  )+
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "",
    y = "Vaccine Effectiveness",
  ) +
  theme_pubr(base_size = 14)+
  facet_grid(OutcomeType ~ Comparison, scales = "free_y", space = "free_y", switch = "y") +
  coord_flip()+
  theme(panel.spacing = unit(1.5, "lines") )


ve_table <- ve_combined.multi %>%
  filter(Comparison == "RV1-specific", OutcomeType != "Distance_RV5") %>%
  mutate(
    VE_display = sprintf("%.0f (%.0f, %.0f)", VE, CI_lower, CI_upper),
    Outcome = factor(Outcome, levels = rev(unique(Outcome)))  # match plot order
  ) %>%
  select(OutcomeType, Comparison, Outcome, VE_display)
ve_table <- ve_combined.multi %>%
  filter(Comparison == "RV5-specific", OutcomeType != "Distance_RV1") %>%
  mutate(
    VE_display = sprintf("%.0f (%.0f, %.0f)", VE, CI_lower, CI_upper),
    Outcome = factor(Outcome, levels = rev(unique(Outcome)))  # match plot order
  ) %>%
  select(OutcomeType, Comparison, Outcome, VE_display)



gg_rv1 + gg_rv5
ggsave(path = "~/OneDrive - Yale University/Research/PitzerLab_RVAWGS/Final_approvedVersions/2503_revision1/",
       filename = "ve_summary_vaxx_specific.tiff",
       width =11, height = 9, device='tiff', dpi=200)

# post hoc tests  -------------------------------------------------------------# 

