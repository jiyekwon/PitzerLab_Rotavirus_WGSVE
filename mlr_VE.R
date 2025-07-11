library(ggplot2)
library(ggpubr)
library(patchwork)
library(broom)
library(dplyr)
library(car)
library(nnet)
library(purrr)
library(ggridges)
library(forcats)
library(questionr)
source("functions.R")


#  Sieve analysis   ------------------------------------
select.all$colYear2<- factor(select.all$colYear, 
                                 levels = c("2014",  "2012", "2013", "2015", "2016"),
                                 labels = c("2014", "2012", "2013", "2015", "2016"))

## Run RV1-Specific Analysis -----------------------------
result_rv1 <- run_sieve_analysis(
  df = select.all,
  perc_col = "percTotal.rv1",
  vax_to_exclude = "Vaccinated_rotateq",
  breaks_vec = c(0, 9.6, 30),
  label = "RV1"
)


## Run RV5-Specific Analysis -----------------------------
result_rv5 <- run_sieve_analysis(
  df = select.all,
  perc_col = "percTotal.rv5",
  vax_to_exclude = "Vaccinated_rotarix",
  breaks_vec = c(0, 16.7, 30),
  label = "RV5"
)

## Summary  -----------------------------
cat("\n### AIC Comparison ###\n")
print(data.frame(
  Model = names(result_rv1$aic),
  AIC_RV1 = round(result_rv1$aic, 2),
  AIC_RV5 = round(result_rv5$aic, 2)
))

cat("\n### Odds Ratios (RV1) ###\n")
print(result_rv1$OR_unadjusted[2,])
print(result_rv1$OR_adj_colYear2[2,])
print(result_rv1$OR_adj_full[2,])

cat("\n### Odds Ratios (RV5) ###\n")
print(result_rv5$OR_unadjusted[2,])
print(result_rv5$OR_adj_colYear2[2,])
print(result_rv5$OR_adj_full[2,])


#  VE comparison by Gxpx, protection types  ------------------------------------
##  Define genotypes and protection types of interest --------------------------
vaccine_groups <- c("Overall", "RV5-specific", "RV1-specific")
outcome_types <- c("GxPx", "Protection", "Distance")
genotypes <- c("G1P[8]", "G12P[8]", "G2P[4]", "G9P[8]")
genotypes.comp<- c("G1P[8]", "G2P[8]", "G3P[8]", "G9P[8]", "G12P[8]", "G1P[4]", "G2P[4]", "G3P[6]")
protection_types <- c("homotypic", "partially heterotypic", "fully heterotypic")


##  Data cleaning -------------------------------------------------------------
ve.all.data <- ve.all.data %>%
  mutate(protection = case_when(
    is.na(GxPx) ~ NA_character_,
    GxPx == "G1P[8]" ~ "homotypic",
    grepl("G1", GxPx) | grepl("P\\[8\\]", GxPx) ~ "partially heterotypic",
    TRUE ~ "fully heterotypic"
  )) %>%
  mutate(protection = factor(protection, levels = protection_types),
         GxPx = factor(GxPx, levels = genotypes.comp))



## Subset:  --------------------------------------------------------------------
# Filter RV1 only -- RV1-specific analysis
ve.all.data <- ve.all.data %>% 
  mutate(strainSpecific = case_when(strainSpecific == "StrainS" ~ "StrainS",
                                    strainSpecific == "StrainI" ~ "StrainI",
                                    strainSpecific == "control" ~ "Control"))

rv1.only <- ve.all.data %>% filter(VaccineType != "RV5-vaccinated")
rv5.only <- ve.all.data %>% filter(VaccineType!= "RV1-vaccinated")
overall.anyvax <- ve.all.data %>% filter(VaccineType %in% c("RV1-vaccinated", "RV5-vaccinated", "Unvaccinated"))



##  Output:  --------------------------------------------------------------------
# Overall VE to any infection (not genotype/protection/distance stratified)
ve_overall_any <- bind_rows(
  get_ve_summary(rv1.only, "RV1-specific", adjusted = FALSE),
  get_ve_summary(rv1.only, "RV1-specific", adjusted = TRUE),
  get_ve_summary(rv5.only, "RV5-specific", adjusted = FALSE),
  get_ve_summary(rv5.only, "RV5-specific", adjusted = TRUE)
) %>% mutate(Outcome = "Any Infection", OutcomeType = "Overall", Comparison = label)
ve_overall_any$Comparison <- c(rep("RV1-specific",2),rep("RV5-specific",2))

ve_oa_multi <- get_ve_summary(rv1.only, "RV1-specific", adjusted = TRUE)%>% mutate(Comparison = Group, OutcomeType= "Overall", Outcome = "Any infection") %>% select(-c(Model, Group))
ve_dist_multi <- run_multinom_ve_analysis(data = rv1.only, comparison_label = "RV1-specific",  outcome_type = "strainSpecific")[[1]]
ve_protection_multi <- run_multinom_ve_analysis(data = rv1.only, comparison_label = "RV1-specific", outcome_type = "protection")[[1]]
ve_gxpx_multi <- run_multinom_ve_analysis(data = rv1.only, comparison_label = "RV1-specific", outcome_type = "GxPx")[[1]]

ve_oa_multi5 <- get_ve_summary(rv5.only, "RV5-specific", adjusted = TRUE) %>% mutate(Comparison = Group, OutcomeType= "Overall", Outcome = "Any infection") %>% select(-c(Model, Group))
ve_dist_multi5 <- run_multinom_ve_analysis(data = rv5.only, comparison_label = "RV5-specific",  outcome_type = "clusterAssignment.rv5")[[1]]
ve_protection_multi5 <- run_multinom_ve_analysis(data = rv5.only, comparison_label = "RV5-specific", outcome_type = "protection")[[1]]
ve_gxpx_multi5 <- run_multinom_ve_analysis(data = rv5.only, comparison_label = "RV5-specific", outcome_type = "GxPx")[[1]]

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

## Figure 3  -------------------------------------------------------------------
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

## post hoc test: wald test, model comparison  ------------------------------------------------------------- 
# H0: Test if vaccine effect j0 = j1
mod_rv1.dist <- run_multinom_ve_analysis(data = rv1.only, comparison_label = "RV1-specific",  outcome_type = "strainSpecific")[[2]]
mod_rv1.gxpx <- run_multinom_ve_analysis(data = rv1.only, comparison_label = "RV1-specific", outcome_type = "GxPx")[[2]]
AIC(mod_rv1.dist, mod_rv1.gxpx)

mod_rv5.dist<- run_multinom_ve_analysis(data = rv5.only, comparison_label = "RV5-specific",  outcome_type = "clusterAssignment.rv5")[[2]]
mod_rv5.gxpx <- run_multinom_ve_analysis(data = rv5.only, comparison_label = "RV5-specific", outcome_type = "GxPx")[[2]]
AIC(mod_multi.dist, mod_multi.gxpx)



b <- as.vector(t(coef(mod_multi)))  # coefficients for j0, then j1
Sigma <- vcov(mod_multi)
coef_names <- colnames(Sigma)
j0_idx <- which(coef_names == "j=0:VaccineBinary")
j1_idx <- which(coef_names == "j=1:VaccineBinary")

L <- rep(0, length(b))
L[j0_idx] <- 1
L[j1_idx] <- -1
wald.test(b = b, Sigma = Sigma, L = matrix(L, nrow = 1))


#  GD distribution  by site   --------------------------------------------------

## data cleaning ------
select.all$vaxx_binary<- ifelse(select.all$VaxxStatEver_4 == "Unvaccinated", "No", "Yes")

## Figure A: individual points by study site ----------------------------------
figure_state.A <- select.all %>% 
  ggplot(aes( x = percTotal.rv1, y = fct_reorder(studySite, RV1.vaxx),color = VaxxStatEver_4 )) +
  geom_jitter(height = 0.1, width = 1, size = 3, alpha = 0.7) +
  # scale_color_manual(values = c("#A8ACADFF",  "maroon"), "Vaccination status")+
  scale_color_manual(values = c("#A8ACADFF", "#1B608F", "#D6CBB5FF"), "Vaccination status")+
  labs( x = "Genetic Distance to RV1", y = "Study Site (ordered by RV1 coverage)") +
  theme_minimal() +
  theme( panel.grid.minor = element_blank(),
         panel.grid.major.y = element_blank())

## ---- Figure B: Summary per study site ----------------------------------------
city_summary_rv1 <- get_city_summary(select.all, percTotal.rv1, RV1.vaxx)

figure_state.B<- ggplot(city_summary_rv1, aes(x = coverage*100, y = mean_distance, color = coverage)) +
  geom_errorbar(aes(ymin = mean_distance - se_distance, ymax = mean_distance + se_distance), 
                width = 1, alpha = 0.6)+
  geom_point(aes(size = n_samples)) +
  geom_text_repel(aes(label = studySite), nudge_x = 2,color = "grey30", size = 3.5, max.overlaps = 10) +
  scale_size_continuous(name = "Sample Size") +
  labs( x = "RV1 Use (%)", y = "Mean Genetic Distance to RV1", title = "") +
  theme_minimal()+
  scico::scale_color_scico(palette = "navia",direction = -1, end = 0.82, 
                           begin = 0.2,
                           alpha = 0.9,
                           name = "RV1 vaccine use")

combined_plot <- figure_state.A + plot_spacer()+ figure_state.B + 
  plot_layout(ncol = 3, widths = c(1, 0.2, 1)) + 
  plot_annotation(tag_levels = "A") +  
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") 

ggsave(path = "~/OneDrive - Yale University/Research/PitzerLab_RVAWGS/Final_approvedVersions/2503_revision1/",
       filename = "location_summary_vaxx_3group.tiff",
       width =14, height = 7, device='tiff', dpi=200)

## ---- RV5 plots  --------------------------------------------------------------
figure_state.A.5 <- select.all %>% 
  ggplot(aes( x = percTotal.rv5, y = fct_reorder(studySite, RV5.vaxx),color = VaxxStatEver_4 )) +
  geom_jitter(height = 0.1, width = 1, size = 3, alpha = 0.7) +
  # scale_color_manual(values = c("#A8ACADFF",  "maroon"), "Vaccination status")+
  scale_color_manual(values = c("#A8ACADFF",  "#D6CBB5FF", "maroon" ), "Vaccination status")+
  labs( x = "Genetic Distance to RV5", y = "Study Site (ordered by RV5 coverage)") +
  theme_minimal() +
  theme( panel.grid.minor = element_blank(),
         panel.grid.major.y = element_blank())


city_summary_rv5 <- get_city_summary(select.all, percTotal.rv5, RV5.vaxx)
figure_state.B.5 <- ggplot(city_summary_rv5, 
                           aes(x = coverage * 100, 
                               y = mean_distance, color = coverage)) +
  #geom_jitter(inherit.aes = FALSE, data = select.all, aes(x = RV5.vaxx*100, y = percTotal.rv5),  color = "grey80", alpha = 0.6, size = 1)+
  geom_errorbar(aes(ymin = mean_distance - se_distance, ymax = mean_distance + se_distance), 
                width = 1, alpha = 0.6) +
  geom_point(aes(size = n_samples)) +
  geom_text_repel(aes(label = studySite), 
                  nudge_x = 2, color = "grey30", 
                  size = 3.5, max.overlaps = 10) +
  scale_size_continuous(name = "Sample Size") +
  scale_color_viridis_c(
    option = "rocket",
    name = "RV5 vaccine use",
    begin = 0.2, end = 0.9, 
    direction = -1, alpha = 0.9,
    breaks = c(0, 0.25, 0.50, 0.75),
    limits = c(0, 0.75) ) +
  labs( x = "RV5 Use (%)", y = "Mean Genetic Distance to RV5" ) +
  ylim(14, 19)+
  theme_minimal()

combined_plot.5 <- figure_state.A.5 + plot_spacer()+ figure_state.B.5 + 
  plot_layout(ncol = 3, widths = c(1, 0.2, 1)) + 
  plot_annotation(tag_levels = "A") +  
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") 

combined_plot.5
ggsave(path = "~/OneDrive - Yale University/Research/PitzerLab_RVAWGS/Final_approvedVersions/2503_revision1/",
       filename = "location_summary_vaxx_3group_rv5.tiff",
       width =14, height = 7, device='tiff', dpi=200)


combined_plot/ combined_plot.5 + plot_annotation(tag_levels = "A") 
ggsave(path = "~/OneDrive - Yale University/Research/PitzerLab_RVAWGS/Final_approvedVersions/2503_revision1/",
       filename = "location_combined_2.tiff",
       width =14, height = 12, device='tiff', dpi=200)