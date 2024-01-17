library(lme4)
library(lmerTest)
library(dplyr)
library(ggpubr)
library(glmmTMB)
library(ggdist)
library(ggridges)
library(DHARMa)
library(tidyr)
library(plotrix)
library(ggdist)

transfomation <- "sqrt"
logtype <- "1"
# Read in the Random simulations dat  a
analysis_data<- read.csv("data/analysis_data_trait_0123.csv") %>% 
  distinct(SSBS, .keep_all = TRUE) #%>%
#filter(Species_richness > 4)

analysis_data$LandUse[analysis_data$LandUse == "Primary minimal"] <- "Pristine primary"

# Read in additional covariates and combine
diversity <- readRDS("data/diversity_combined_sites_Aug2022.rds") %>%
  dplyr::select(SSBS, Predominant_habitat)

## Split the Urbans as they give divergent results
analysis_data$LandUse[which(analysis_data$LandUseIntensity == "Urban-Minimal use")] <- "Urban - Minimal"
analysis_data$LandUse[which(analysis_data$LandUse == "Urban")] <- "Urban - Intensive"
analysis_data$LandUse[which(analysis_data$LandUse == "Plantation forest")] <- "Plantation"


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
analysis_data$LandUse[which(grepl("Young secondary vegetation", analysis_data$LandUse))] <- "Young secondary"
analysis_data$LandUse[which(grepl("Intermediate secondary vegetation", analysis_data$LandUse))] <- "Intermediate secondary"
analysis_data$LandUse[which(grepl("Mature secondary vegetation", analysis_data$LandUse))] <- "Mature secondary"


# Remove undeeded datasets
analysis_data  <- analysis_data%>%
  filter(!SS %in% c("LK1_2009__Hayward 1", # not entire community ##Patrick
                    "GN1_2010__Hvenegaard 1"), # pseudo replicate
         !LandUse == "Secondary vegetation (indeterminate age)") 

# relevel the LandUse column
analysis_data$LandUse <- factor(analysis_data$LandUse, levels = c("Pristine primary",
                                                                  "Disturbed primary",
                                                                  "Mature secondary",
                                                                  "Intermediate secondary",
                                                                  "Young secondary",
                                                                  "Plantation",
                                                                  "Pasture",
                                                                  "Cropland",
                                                                  "Urban - Minimal",
                                                                  "Urban - Intensive"))




# Scale and transform columns
# SQRT transform
analysis_data$FRich <-  as.numeric(scale(sqrt(analysis_data$FRich)))
analysis_data$Red[analysis_data$Red < 0] <- 0
analysis_data$Red <- as.numeric(scale(sqrt(analysis_data$Red)))
analysis_data$Species_richness_scaled <- as.numeric(scale(sqrt(analysis_data$Species_richness)))

################################################################################
## Plot all ##
################################################################################
######################################## EXTENDED DATA FIGURE WITH ABUNDANCE ###########################
analysis_data$cov_pearsons_abundance_redraw <- as.numeric(scale(analysis_data$cov_pearsons_abundance_redraw))

LandUses_cov_study <- lmer(cov_pearsons_abundance_redraw ~  LandUse  + (1|SS) + (1|SSB), data = analysis_data %>% filter(cov_pearsons_abundance_redraw > -10))
summary(LandUses_cov_study)
data_all_cov <- as.data.frame(summary(LandUses_cov_study)$coefficient)
data_all_cov$LandUse <- c("Pristine primary", "Disturbed primary", "Mature secondary", "Intermediate secondary", "Young secondary",
                          "Plantation", "Pasture","Cropland", "Urban - Minimal", "Urban - Intensive")

data_all_cov$LandUse <- factor(data_all_cov$LandUse, levels = c("Pristine primary", "Disturbed primary", "Mature secondary", "Intermediate secondary", "Young secondary",
                                                                "Plantation", "Pasture","Cropland", "Urban - Minimal", "Urban - Intensive"))

colnames(data_all_cov)[2] <- "se"
data_all_cov$Estimate[1] <- 0


ED_2 <- ggplot(data_all_cov) +
  geom_point(aes(y = Estimate, x = LandUse, col = LandUse), position = position_dodge(0.5), size = 6) +
  geom_errorbar(aes(x = LandUse, ymin = Estimate - (1.96*se), ymax = Estimate + (1.96*se), col = LandUse), width=.3, size = 1, position = position_dodge(0.5)) +
  geom_hline(yintercept = 0, col = "red", cex = 0.5, linetype = "dashed") +
  #xlim(-1.0, 0.6) +
  ylab("Coefficient estimate") +
  xlab("LandUse") +
  theme_classic() +
  scale_color_manual(values = c("Pristine primary" = "red", 
                                "Disturbed primary" = "firebrick", 
                                "Mature secondary" = "#4EADC7", 
                                "Intermediate secondary" = "#7EB1D8", 
                                "Young secondary" = "#A4C2FA", 
                                "Plantation" = "#EBF500", 
                                "Cropland" = "gold", 
                                "Pasture" = "#FFF700", 
                                "Urban - Minimal" = "lightgray", 
                                "Urban - Intensive" = "darkgray"))  +
  #facet_wrap(~Metric) +
  #scale_x_discrete(labels = c("-1.0" = " ", "-0.5" = " ", "0.0" = " ", "0.5" = " ")) +
  theme(axis.text.y = element_text(size = 17),
        axis.text.x = element_text(angle = 60, size = 17, hjust = 1.1),
        axis.title.x = element_blank(),
        #axis.line.x = element_blank(),
        #axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 20),
        axis.title.y = element_text(size = 25) ,
        legend.position = "none",
        plot.margin = margin(1,0.2,1,0, "cm"),
        plot.title = element_text(size = 25, face = "bold", hjust = -0.25))

dev.size("px")
#[1] 3095 2312
png("figures/ED 2 cov_ab_2023.png", width = 795/1.5, height = 733/1.5)
ED_2
dev.off()
