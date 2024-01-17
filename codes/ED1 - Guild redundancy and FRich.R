

sort_dataset <- function(data) {
  
  analysis_data <- data %>% 
    #filter(Species_richness > 4) %>%
    distinct(SSBS, .keep_all = TRUE) 
  
  analysis_data$LandUse[analysis_data$LandUse == "Primary minimal"] <- "Pristine primary"
  
  # Read in additional covariates and combine
  diversity <- readRDS("data/diversity_combined_sites_Aug2022.rds") %>%
    dplyr::select(SSBS, Predominant_habitat)
  
  
  ## Split the Urbans as they give divergent results
  analysis_data$LandUse[which(analysis_data$LandUseIntensity == "Urban-Minimal use")] <- "Urban - Minimal"
  analysis_data$LandUse[which(analysis_data$LandUse == "Urban")] <- "Urban - Intensive"
  
  # Change the Land-use with Secondary veg into the original land-use type
  # subset the secondary veg
  second_veg <- analysis_data %>% filter(LandUse == "Secondary vegetation")
  
  # remove the subset from the analysis_data
  analysis_data <- analysis_data %>% filter(LandUse != "Secondary vegetation")
  
  # restore the original LandUse type from the data
  second_veg$LandUse <- diversity$Predominant_habitat[match(second_veg$SSBS, diversity$SSBS)]
  
  # bind it back to the analysis data
  analysis_data <- rbind(analysis_data, second_veg)
  
  # abbreviate the names in the LU column 
  analysis_data$LandUse[which(grepl("Young secondary vegetation", analysis_data$LandUse))] <- "YSV"
  analysis_data$LandUse[which(grepl("Intermediate secondary vegetation", analysis_data$LandUse))] <- "ISV"
  analysis_data$LandUse[which(grepl("Mature secondary vegetation", analysis_data$LandUse))] <- "MSV"
  
  
  # Remove undeeded datasets
  analysis_data  <- analysis_data %>%
    filter(!SS %in% c("LK1_2009__Hayward 1", # not entire community ##Patrick
                      "GN1_2010__Hvenegaard 1"), # pseudo replicate
           !LandUse == "Secondary vegetation (indeterminate age)") 
  
  # relevel the LandUse column
  analysis_data$LandUse[grep("Urban", analysis_data$LandUse)] <- "Urban"
  analysis_data$LandUse <- factor(analysis_data$LandUse, levels = c("Pristine primary",
                                                                    "Disturbed primary",
                                                                    "MSV",
                                                                    "ISV",
                                                                    "YSV",
                                                                    "Plantation forest",
                                                                    "Pasture",
                                                                    "Cropland",
                                                                    "Urban"))
  
  
  # Scale and transform columns
  # SQRT transform
  analysis_data$Resistance_spline <- scale(analysis_data$Resistance_spline)
  analysis_data$Resistance_step <- scale(analysis_data$Resistance_step)
  
  analysis_data$cov_pearsons_rand_red <- as.numeric(scale(analysis_data$cov_pearsons_random_sensitivity_score_redraw))
  analysis_data$cov_spearman_rand_red <- as.numeric(scale(analysis_data$cov_spearman_random_sensitivity_score_redraw))
  
  analysis_data$cov_pearsons_sense_red <- as.numeric(scale(analysis_data$cov_pearsons_sense_redraw))
  analysis_data$cov_spearman_sense_red <- as.numeric(scale(analysis_data$cov_spearman_sense_redraw))
  
  analysis_data$cov_pearsons_abundance_red <- as.numeric(scale(-analysis_data$cov_pearsons_abundance_redraw))
  analysis_data$cov_spearman_abundance_red <- as.numeric(scale(-analysis_data$cov_spearman_abundance_redraw))
  
  
  analysis_data$Climatic_zone <- "Tropics"
  analysis_data$Climatic_zone[which(abs(analysis_data$Latitude) > 20)] <- "Temperate"
  
  
  return(analysis_data)
  
}

## function to add 0s to NAs provided that study has at least 1 non-NA site
## if study all NAs it will be excluded from models
## But if some study sites have assemblage we assume LU causes lack of species
fix_NAs <- function(guild_results) {
  colnames(guild_results)[2:4]  <- c("Resistance_guild", "Red_guild", "FRich_guild")

all_results <- left_join(analysis_data_tw %>% dplyr::select(SS, SSB, SSBS, LandUse), guild_results)

Studies_to_model <- all_results %>% 
  group_by(SS) %>%
  summarise(n_NAs = length(which(is.na(Red_guild))),
            n_sites = n()) %>%
  mutate(prop_NAs = n_NAs/n_sites) %>%
  ungroup() %>%
  filter(prop_NAs < 1) %>%
  pull(SS) 

guild_results <- all_results %>% filter(SS %in% Studies_to_model)
guild_results$FRich_guild[is.na(guild_results$FRich_guild)] <- 0
guild_results$FRich_guild <-  as.numeric(scale(sqrt(guild_results$FRich_guild)))

guild_results$Red_guild[is.na(guild_results$Red_guild)] <- 0
guild_results$Red_guild[guild_results$Red_guild < 0] <- 0
guild_results$Red_guild <- as.numeric(scale(sqrt(guild_results$Red_guild)))


return(guild_results)

}

## Get the trait weighted dataset which is ready to model
analysis_data_tw <- read.csv("data/analysis_data_trait_0123.csv") %>% sort_dataset()


frug_results <- read.csv("data/frug_resistance_trait.csv") %>%
  dplyr::select(SSBS, Resistance_spline, Red, FRich) %>% fix_NAs()

inv_results <- read.csv("data/inv_resistance_trait.csv") %>%
  dplyr::select(SSBS, Resistance_spline, Red, FRich) %>% fix_NAs()

omni_results <- read.csv("data/omniv_resistance_trait.csv") %>%
  dplyr::select(SSBS, Resistance_spline, Red, FRich) %>% fix_NAs()

gran_results <- read.csv("data/graniv_resistance_trait.csv") %>%
  dplyr::select(SSBS, Resistance_spline, Red, FRich) %>% fix_NAs()


########## FRich models ###########
frug_mod_FRich <- lmer(FRich_guild ~ LandUse + (1|SSB) + (1|SS), data = frug_results)
inv_mod_FRich <- lmer(FRich_guild ~ LandUse + (1|SSB) + (1|SS), data = inv_results)
gran_mod_FRich <- lmer(FRich_guild ~ LandUse + (1|SSB) + (1|SS), data = gran_results)
omni_mod_FRich <- lmer(FRich_guild ~ LandUse + (1|SSB) + (1|SS), data = omni_results)

frug_FRich <- as.data.frame(summary(frug_mod_FRich)$coef)
inv_FRich <- as.data.frame(summary(inv_mod_FRich)$coef)
gran_FRich <- as.data.frame(summary(gran_mod_FRich)$coef)
omni_FRich <- as.data.frame(summary(omni_mod_FRich)$coef)

all_FRich <- do.call(rbind, list(inv_FRich, frug_FRich, gran_FRich, omni_FRich))

all_FRich$Sample <- c(rep("Invertivore", 9), rep("Frugivore", 9), rep("Granivore", 9), rep("Generalist", 9))
all_FRich$LandUse <- rep(c("Pristine primary", 
                              "Disturbed primary", 
                              "Mature secondary", 
                              "Intermediate secondary", 
                              "Young secondary", 
                              "Plantation", 
                              "Pasture", 
                              "Cropland", 
                              "Urban"))  


colnames(all_FRich)[2] <- "se"
all_FRich$Estimate[all_FRich$LandUse == "Pristine primary"] <- 0
all_FRich$LandUse <- factor(all_FRich$LandUse, levels = c("Pristine primary", 
                                                                    "Disturbed primary", 
                                                                    "Mature secondary", 
                                                                    "Intermediate secondary", 
                                                                    "Young secondary", 
                                                                    "Plantation", 
                                                                    "Pasture", 
                                                                    "Cropland", 
                                                                    "Urban"))

########### Red Models ###################
frug_mod_Red <- lmer(Red_guild ~ LandUse + (1|SSB) + (1|SS), data = frug_results)
inv_mod_Red <- lmer(Red_guild ~ LandUse + (1|SSB) + (1|SS), data = inv_results)
gran_mod_Red <- lmer(Red_guild ~ LandUse + (1|SSB) + (1|SS), data = gran_results)
omni_mod_Red <- lmer(Red_guild ~ LandUse + (1|SSB) + (1|SS), data = omni_results)

frug_Red <- as.data.frame(summary(frug_mod_Red)$coef)
inv_Red <- as.data.frame(summary(inv_mod_Red)$coef)
gran_Red <- as.data.frame(summary(gran_mod_Red)$coef)
omni_Red <- as.data.frame(summary(omni_mod_Red)$coef)

all_Red <- do.call(rbind, list(inv_Red, frug_Red, gran_Red, omni_Red))

all_Red$Sample <- c(rep("Invertivore", 9), rep("Frugivore", 9), rep("Granivore", 9), rep("Generalist", 9))
all_Red$LandUse <- rep(c("Pristine primary", 
                           "Disturbed primary", 
                           "Mature secondary", 
                           "Intermediate secondary", 
                           "Young secondary", 
                           "Plantation", 
                           "Pasture", 
                           "Cropland", 
                           "Urban"))  


colnames(all_Red)[2] <- "se"
all_Red$Estimate[all_Red$LandUse == "Pristine primary"] <- 0
all_Red$LandUse <- factor(all_Red$LandUse, levels = c("Pristine primary", 
                                                          "Disturbed primary", 
                                                          "Mature secondary", 
                                                          "Intermediate secondary", 
                                                          "Young secondary", 
                                                          "Plantation", 
                                                          "Pasture", 
                                                          "Cropland", 
                                                          "Urban"))







all_Red$Metric <- "Functional redundancy"
all_FRich$Metric <- "Functional diversity"

all_results <- rbind(all_Red, all_FRich)


all_results$Sample <- factor(all_results$Sample, levels = c("Invertivore", "Frugivore", "Granivore", "Generalist"))

RichRed_fig <- ggplot(all_results) +
  geom_rect(xmin = 0.4, xmax = 1.5, ymin = -Inf, ymax = Inf,fill = "#ededed", color = "#ededed", alpha = .05) +
  geom_rect(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf,fill = "#ededed", color = "#ededed", alpha = .05) +
  geom_rect(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf,fill = "#ededed", color = "#ededed", alpha = .05) +
  geom_rect(xmin = 6.5, xmax = 7.5, ymin = -Inf, ymax = Inf,fill = "#ededed", color = "#ededed", alpha = .05) +
  geom_rect(xmin = 8.5, xmax = 9.5, ymin = -Inf, ymax = Inf,fill = "#ededed", color = "#ededed", alpha = .05) +
  geom_point(aes(y = Estimate , x = LandUse, col = LandUse, shape = Sample), position = position_dodge(0.6), size = 6) +
  geom_errorbar(aes(x = LandUse, ymin = (Estimate - 1.96*se), ymax = (Estimate + 1.96*se), col = LandUse, group = Sample), width=.2, size = 1.25, position = position_dodge(0.6)) +
  geom_hline(yintercept = 0, col = "red", cex = 0.5, linetype = "dashed") +
  #ggtitle("b") +
  xlab("") +
  ylab("Coefficient estimate") +
  ylim(-1.1, 0.72)  +
  theme_bw() + 
  # geom_text(inherit.aes = FALSE,
  #           data    = dat_text,
  #           mapping = aes(x = -Inf, y = -Inf, label = label),
  #           hjust   = -1,
  #           vjust   = -17,
  #           size = 10,
  #           col = "black",
  #           fontface = "bold") +
  facet_wrap(~Metric, ncol = 1) +
  scale_color_manual(values = c("Pristine primary" = "red", 
                                "Disturbed primary" = "firebrick", 
                                "Mature secondary" = "#4EADC7", 
                                "Intermediate secondary" = "#7EB1D8", 
                                "Young secondary" = "#A4C2FA", 
                                "Plantation" = "#EBF500", 
                                "Cropland" = "gold2", 
                                "Pasture" = "gold", 
                                "Urban" = "darkgray"))  +
  guides(color=FALSE) +
  scale_shape_manual(values = c("Invertivore" = 19, "Frugivore" = 17, "Granivore" = 15, "Generalist" = 7)) +
  geom_hline(yintercept = 0,color = "red", cex = 0.5, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 60, size = 27, hjust = 1.1),
        axis.text.y = element_text(size = 27),
        axis.title.y = element_text(size = 35, vjust = 2.2) ,
        legend.position = "top",
        strip.text = element_text(size = 30),
        legend.text = element_text(size = 28),
        legend.title =  element_blank(),
        panel.grid = element_blank())
RichRed_fig


png("figures/ED1 - Between guilds.png", width = 1190/1.5, height = 1446/1.5)
RichRed_fig
dev.off()
