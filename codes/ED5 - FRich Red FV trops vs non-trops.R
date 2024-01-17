library(lme4)
library(lmerTest)
library(tidyverse)
library(ggpubr)
library(glmmTMB)

transformation <- "sqrt"
logtype <- 0.1
###############  INVERTIVORES  ###############################
# Read in the Random simulations 
analysis_data<- read.csv("data/Tropics_richness_redundancy_cov.csv")

# Read in additional covariates and combine
diversity <- readRDS("data/diversity_combined_sites_tropics_Aug2022.rds") %>%
  dplyr::select(SSBS, SS, SSB, Predominant_habitat, LandUse, LandUseIntensity)

diversity$LandUse <- as.character(diversity$LandUse)
diversity$LandUseIntensity <- as.character(diversity$LandUseIntensity)
diversity$SS <- as.character(diversity$SS)
diversity$SSB <- as.character(diversity$SSB)
diversity$SSBS <- as.character(diversity$SSBS)

analysis_data$SSB <- diversity$SSB[match(analysis_data$SSBS, diversity$SSBS)]
analysis_data$SS <- diversity$SS[match(analysis_data$SSBS, diversity$SSBS)]
analysis_data$LandUse<-diversity$LandUse[match(analysis_data$SSBS, diversity$SSBS)]
analysis_data$LandUseIntensity<-diversity$LandUseIntensity[match(analysis_data$SSBS, diversity$SSBS)]

analysis_data$LandUse <- as.character(analysis_data$LandUse)

## Split the Urbans as they give divergent results
analysis_data$LandUse[which(analysis_data$LandUse == "Primary minimal")] <- "Pristine primary"
analysis_data$LandUse[which(analysis_data$LandUseIntensity == "Urban-Minimal use")] <- "Urban"
analysis_data$LandUse[which(analysis_data$LandUse == "Urban")] <- "Urban"
analysis_data$LandUse[which(analysis_data$LandUse == "Plantation forest")] <- "Plantation"

# Change the Land-use with Secondary veg into the original land-use type
# subset the secondary veg
second_veg <- analysis_data %>% dplyr::filter(LandUse == "Secondary vegetation")

# remove the subset from the analysis_data
analysis_data <- analysis_data %>% dplyr::filter(LandUse != "Secondary vegetation")

# restore the original LandUse type from the data
second_veg$LandUse <- diversity$Predominant_habitat[match(second_veg$SSBS, diversity$SSBS)]

# bind it back to the analysis data
analysis_data <- rbind(analysis_data, second_veg)

# abbreviate the names in the LU column 
# abbreviate the names in the LU column 
analysis_data$LandUse[which(grepl("Young secondary vegetation", analysis_data$LandUse))] <- "Young secondary"
analysis_data$LandUse[which(grepl("Intermediate secondary vegetation", analysis_data$LandUse))] <- "Intermediate secondary"
analysis_data$LandUse[which(grepl("Mature secondary vegetation", analysis_data$LandUse))] <- "Mature secondary"

# Remove undeeded datasets
analysis_data  <- analysis_data%>%
  dplyr::filter(!SS %in% c("LK1_2009__Hayward 1", # not entire community ##Patrick
                    "GN1_2010__Hvenegaard 1"), # pseudo replicate
         !LandUse == "Secondary vegetation") 

# relevel the LandUse column
analysis_data$LandUse <- factor(analysis_data$LandUse, levels = c("Pristine primary",
                                                                  "Disturbed primary",
                                                                  "Secondary vegetation",
                                                                  "Mature secondary",
                                                                  "Intermediate secondary",
                                                                  "Young secondary",
                                                                  "Plantation",
                                                                  "Pasture",
                                                                  "Cropland",
                                                                  "Urban"))


# Scale and transform columns
analysis_data$cov_pearsons_sense_redraw <-  as.numeric(scale(analysis_data$cov_pearsons_sense_redraw))
analysis_data$FRich <-  as.numeric(scale(sqrt(analysis_data$FRich)))
analysis_data$Red[analysis_data$Red < 0] <- 0
analysis_data$Red <- as.numeric(scale(sqrt(analysis_data$Red)))
#analysis_data$Species_richness_scaled <- as.numeric(scale(sqrt(analysis_data$Species_richness)))

  # save an analysis_data without the RAS dataset (currently slightly erroneous)
analysis_data_nRAS <- analysis_data %>%
  dplyr::filter(SS != "2020_Lees_RAS")

analysis_data_nUrbMin <- analysis_data %>%
  dplyr::filter(LandUseIntensity != "Urban-Minimal use")


# Run a standard Richness model
LandUses_rich_study <-  lmer(FRich ~ LandUse + (1|SS) + (1|SSB), data = analysis_data)
sample_size_tropics <- length(unique(analysis_data$SSBS))

# save the results into a dataframe
model_results_rich_Tropics <- as.data.frame(summary(LandUses_rich_study)$coefficients)#[-c(2),]

# rename the LandUses
model_results_rich_Tropics$LandUse <- c("Pristine primary", "Disturbed primary", "Mature secondary",
                                        "Intermediate secondary",
                                        "Young secondary",
                                        "Plantation","Pasture", "Cropland","Urban")

# relevel it (the extra levels are currently meaningless but may bring them back in)
model_results_rich_Tropics$LandUse <- factor(model_results_rich_Tropics$LandUse, levels = c("Pristine primary", "Disturbed primary", "Mature secondary",
                                                                                            "Intermediate secondary",
                                                                                            "Young secondary",
                                                                                            "Plantation", "Pasture", "Cropland","Urban"))

# Center the first variable
model_results_rich_Tropics$Estimate[1] <- 0

# rename the same second column to the standard error
colnames(model_results_rich_Tropics)[2] <- "se"



# Run a standard Richness model
LandUses_red_study <-  lmer(Red ~ LandUse + (1|SS) + (1|SSB), data = analysis_data)

# save the results into a dataframe
model_results_red_Tropics <- as.data.frame(summary(LandUses_red_study)$coefficients)#[-c(2),]

# rename the LandUses
model_results_red_Tropics$LandUse <- c("Pristine primary", "Disturbed primary", "Mature secondary",
                                        "Intermediate secondary",
                                        "Young secondary",
                                        "Plantation", "Pasture","Cropland", "Urban")

# relevel it (the extra levels are currently meaningless but may bring them back in)
model_results_red_Tropics$LandUse <- factor(model_results_red_Tropics$LandUse, levels = c("Pristine primary", "Disturbed primary", "Mature secondary",
                                                                                            "Intermediate secondary",
                                                                                            "Young secondary",
                                                                                            "Plantation","Pasture","Cropland", "Urban"))

# Center the first variable
model_results_red_Tropics$Estimate[1] <- 0

# rename the same second column to the standard error
colnames(model_results_red_Tropics)[2] <- "se"



# Run a standard Red model
LandUses_cov_study <-  lmer(cov_pearsons_sense_redraw*-1 ~ LandUse + (1|SS) + (1|SSB), data = analysis_data)
#LandUses_cov_study_occ <-  glmer(FRich ~ LandUse + (1|SS) + (1|SSB), family = "binomial", data = analysis_data_occ)

# save the results into a dataframe
model_results_cov_Tropics <- as.data.frame(summary(LandUses_cov_study)$coefficients)#[-c(2),]

# rename the LandUses
model_results_cov_Tropics$LandUse <- c("Pristine primary", "Disturbed primary", "Mature secondary",
                                    "Intermediate secondary",
                                    "Young secondary",
                                    "Plantation","Pasture","Cropland", "Urban")

# relevel it (the extra levels are currently meaningless but may bring them back in)
model_results_cov_Tropics$LandUse <- factor(model_results_cov_Tropics$LandUse, levels = c("Pristine primary", "Disturbed primary", "Mature secondary",
                                                                                    "Intermediate secondary",
                                                                                    "Young secondary",
                                                                                    "Plantation","Pasture","Cropland", "Urban"))

# Center the first variable
model_results_cov_Tropics$Estimate[1] <- 0

# rename the same second column to the standard error
colnames(model_results_cov_Tropics)[2] <- "se"



################################################################################







#######################  Temperate  ###################################
# Read in the Random simulations 
# Read in the Random simulations 
analysis_data<- read.csv("data/Temperate_richness_redundancy_cov.csv")

#colnames(analysis_data)[1] <- "SSBS"

# Read in additional covariates and combine
diversity <- readRDS("data/diversity_combined_sites_temperate_Aug2022.rds") %>%
  dplyr::select(SSBS, SS, SSB, Predominant_habitat, LandUse, LandUseIntensity)

analysis_data$SSB <- diversity$SSB[match(analysis_data$SSBS, diversity$SSBS)]
analysis_data$SS <- diversity$SS[match(analysis_data$SSBS, diversity$SSBS)]
analysis_data$LandUse<-diversity$LandUse[match(analysis_data$SSBS, diversity$SSBS)]
analysis_data$LandUseIntensity<-diversity$LandUseIntensity[match(analysis_data$SSBS, diversity$SSBS)]

analysis_data$LandUse <- as.character(analysis_data$LandUse)

## Split the Urbans as they give divergent results
analysis_data$LandUse[which(analysis_data$LandUse == "Primary minimal")] <- "Pristine primary"
analysis_data$LandUse[which(analysis_data$LandUseIntensity == "Urban-Minimal use")] <- "Urban"
analysis_data$LandUse[which(analysis_data$LandUse == "Urban")] <- "Urban"
analysis_data$LandUse[which(analysis_data$LandUse == "Plantation forest")] <- "Plantation"


# Change the Land-use with Secondary veg into the original land-use type
# subset the secondary veg
second_veg <- analysis_data %>% dplyr::filter(LandUse == "Secondary vegetation")

# remove the subset from the analysis_data
analysis_data <- analysis_data %>% dplyr::filter(LandUse != "Secondary vegetation")

# restore the original LandUse type from the data
second_veg$LandUse <- diversity$Predominant_habitat[match(second_veg$SSBS, diversity$SSBS)]

# bind it back to the analysis data
analysis_data <- rbind(analysis_data, second_veg)

# abbreviate the names in the LU column 
# abbreviate the names in the LU column 
analysis_data$LandUse[which(grepl("Young secondary vegetation", analysis_data$LandUse))] <- "Young secondary"
analysis_data$LandUse[which(grepl("Intermediate secondary vegetation", analysis_data$LandUse))] <- "Intermediate secondary"
analysis_data$LandUse[which(grepl("Mature secondary vegetation", analysis_data$LandUse))] <- "Mature secondary"
# Remove undeeded datasets
analysis_data  <- analysis_data%>%
  dplyr::filter(!SS %in% c("LK1_2009__Hayward 1", # not entire community ##Patrick
                    "GN1_2010__Hvenegaard 1"), # pseudo replicate
         !LandUse == "Secondary vegetation") 

# relevel the LandUse column
analysis_data$LandUse <- factor(analysis_data$LandUse, levels = c("Pristine primary",
                                                                  "Disturbed primary",
                                                                  "Secondary vegetation",
                                                                  "Mature secondary",
                                                                  "Intermediate secondary",
                                                                  "Young secondary",
                                                                  "Plantation",
                                                                  "Pasture",
                                                                  "Cropland",
                                                                  "Urban"))

# Scale and transform columns
# SQRT transform
analysis_data$cov_pearsons_sense_redraw <-  as.numeric(scale(analysis_data$cov_pearsons_sense_redraw))
analysis_data$FRich <-  as.numeric(scale(sqrt(analysis_data$FRich)))
analysis_data$Red[analysis_data$Red < 0] <- 0
analysis_data$Red <- as.numeric(scale(sqrt(analysis_data$Red)))

# save an analysis_data without the RAS dataset (currently slightly erroneous)
analysis_data_nRAS <- analysis_data %>%
  dplyr::filter(SS != "2020_Lees_RAS")

analysis_data_nUrbMin <- analysis_data %>%
  dplyr::filter(LandUseIntensity != "Urban-Minimal use")


# Run a standard Richness model
LandUses_rich_study <-  lmer(FRich ~ LandUse + (1|SS) + (1|SSB), data = analysis_data)
sample_size_temp <- length(unique(analysis_data$SSBS))



# save the results into a dataframe
model_results_rich_Temp <- as.data.frame(summary(LandUses_rich_study)$coefficients)#[-c(2),]

# rename the LandUses
model_results_rich_Temp$LandUse <- c("Pristine primary", "Disturbed primary", "Mature secondary",
                                        "Intermediate secondary",
                                        "Young secondary",
                                        "Plantation", "Pasture","Cropland", "Urban")

# relevel it (the extra levels are currently meaningless but may bring them back in)
model_results_rich_Temp$LandUse <- factor(model_results_rich_Temp$LandUse, levels = c("Pristine primary", "Disturbed primary", "Mature secondary",
                                                                                            "Intermediate secondary",
                                                                                            "Young secondary",
                                                                                            "Plantation","Pasture","Cropland", "Urban"))

# Center the first variable
model_results_rich_Temp$Estimate[1] <- 0

# rename the same second column to the standard error
colnames(model_results_rich_Temp)[2] <- "se"




# Run a standard Red model
LandUses_red_study <-  lmer(Red ~ LandUse + (1|SS) + (1|SSB), data = analysis_data)
#LandUses_red_study_occ <-  glmer(FRich ~ LandUse + (1|SS) + (1|SSB), family = "binomial", data = analysis_data_occ)

# save the results into a dataframe
model_results_red_Temp <- as.data.frame(summary(LandUses_red_study)$coefficients)#[-c(2),]

# rename the LandUses
model_results_red_Temp$LandUse <- c("Pristine primary", "Disturbed primary", "Mature secondary",
                                    "Intermediate secondary",
                                    "Young secondary",
                                    "Plantation","Pasture","Cropland", "Urban")

# relevel it (the extra levels are currently meaningless but may bring them back in)
model_results_red_Temp$LandUse <- factor(model_results_red_Temp$LandUse, levels = c("Pristine primary", "Disturbed primary", "Mature secondary",
                                                                                    "Intermediate secondary",
                                                                                    "Young secondary",
                                                                                    "Plantation","Pasture", "Cropland","Urban"))

# Center the first variable
model_results_red_Temp$Estimate[1] <- 0

# rename the same second column to the standard error
colnames(model_results_red_Temp)[2] <- "se"


################################################################################

# Run a standard Red model
LandUses_cov_study <-  lmer(cov_pearsons_sense_redraw*-1 ~ LandUse + (1|SS) + (1|SSB), data = analysis_data)
#LandUses_cov_study_occ <-  glmer(FRich ~ LandUse + (1|SS) + (1|SSB), family = "binomial", data = analysis_data_occ)

# save the results into a dataframe
model_results_cov_Temp <- as.data.frame(summary(LandUses_cov_study)$coefficients)#[-c(2),]

# rename the LandUses
model_results_cov_Temp$LandUse <- c("Pristine primary", "Disturbed primary", "Mature secondary",
                                    "Intermediate secondary",
                                    "Young secondary",
                                    "Plantation","Pasture", "Cropland","Urban")

# relevel it (the extra levels are currently meaningless but may bring them back in)
model_results_cov_Temp$LandUse <- factor(model_results_cov_Temp$LandUse, levels = c("Pristine primary", "Disturbed primary", "Mature secondary",
                                                                                    "Intermediate secondary",
                                                                                    "Young secondary",
                                                                                    "Plantation", "Pasture", "Cropland", "Urban"))

# Center the first variable
model_results_cov_Temp$Estimate[1] <- 0

# rename the same second column to the standard error
colnames(model_results_cov_Temp)[2] <- "se"

 

################################################################################

model_results_red_Temp$Sample <- "Non-tropics"
model_results_red_Tropics$Sample <- "Tropics"

model_results_red_all <- rbind(model_results_red_Temp,
                               model_results_red_Tropics)

model_results_red_all$Sample <- factor(model_results_red_all$Sample, levels = c("Non-tropics", "Tropics"))


model_results_rich_Temp$Sample <- "Non-tropics"
model_results_rich_Tropics$Sample <- "Tropics"

model_results_rich_all <- rbind(model_results_rich_Temp,
                               model_results_rich_Tropics)

model_results_rich_all$Sample <- factor(model_results_rich_all$Sample, levels = c("Tropics", "Non-tropics"))


model_results_cov_Temp$Sample <- "Non-tropics"
model_results_cov_Tropics$Sample <- "Tropics"

model_results_cov_all <- rbind(model_results_cov_Temp,
                               model_results_cov_Tropics)

model_results_cov_all$Sample <- factor(model_results_cov_all$Sample, levels = c("Non-tropics", "Tropics"))


model_results_rich_Temp$Sample <- "Non-tropics"
model_results_rich_Tropics$Sample <- "Tropics"

model_results_rich_all <- rbind(model_results_rich_Temp,
                                model_results_rich_Tropics)

model_results_rich_all$Sample <- factor(model_results_rich_all$Sample, levels = c("Tropics", "Non-tropics"))

model_results_cov_all$Metric <- "FV (Functional vulnerability)"
model_results_rich_all$Metric <- "FD (Richness)"
model_results_red_all$Metric <- "FR (Redundancy)"

model_all <- rbind(model_results_rich_all, model_results_red_all, model_results_cov_all)
model_all$Metric <- factor(model_all$Metric, levels = c("FD (Richness)", "FR (Redundancy)", "FV (Functional vulnerability)"))
model_all$Sample <- factor(model_all$Sample, levels = c("Tropics", "Non-tropics"))

dat_text2 <- data.frame(Sample = c("Tropics", "Non-tropics", "Tropics", "Non-tropics", "Tropics", "Non-tropics"), 
                        Metric = c("FD (Richness)", "FD (Richness)", "FR (Redundancy)", "FR (Redundancy)", "FV (Functional vulnerability)", "FV (Functional vulnerability)"),
                        label = c("a","b", "c", "d", "e", "f"))

dat_text2$Metric <- factor(dat_text2$Metric, levels = c("FR (Redundancy)", "FD (Richness)", "FV (Functional vulnerability)"))
dat_text2$Sample <- factor(dat_text2$Sample, levels = c("Non-tropics", "Tropics"))
#dat_text2$Sample <- 

# plot the Landscape richnesses
pd <- position_dodge(0.5)
richness<- ggplot(model_all, aes(x=LandUse, y=Estimate, fill = LandUse, color = LandUse)) + 
  geom_errorbar(aes(ymin=Estimate-  (1.96*se), ymax=Estimate+ (1.96*se), fill = LandUse, colour=LandUse), width=.2, size = 1, position=pd) +
  geom_point(position=pd, size=6) +
  geom_hline(yintercept = 0,color = "red", cex = 0.5, linetype = "dashed") +
  geom_segment(aes(x = 0.5, xend = 9.5, y = -1.7, yend = -1.7), col = "black") +
  geom_segment(aes(x = 0.5, xend = 0.5, y = 0.75, yend = -1.7), col = "black") +
  theme_classic() +
  scale_color_manual(values = c("Pristine primary" = "red", 
                                "Disturbed primary" = "firebrick", 
                                "Mature secondary" = "#4EADC7", 
                                "Intermediate secondary" = "#7EB1D8", 
                                "Young secondary" = "#A4C2FA", 
                                "Plantation" = "#EBF500", 
                                "Cropland" = "gold", 
                                "Pasture" = "#FFF700", 
                                "Urban" = "gray"))  +
  ggtitle("") +
  scale_fill_manual(values = c("Pristine primary" = "red", 
                                "Disturbed primary" = "firebrick", 
                                "Mature secondary" = "#4EADC7", 
                                "Intermediate secondary" = "#7EB1D8", 
                                "Young secondary" = "#A4C2FA", 
                                "Plantation" = "#EBF500", 
                                "Cropland" = "gold", 
                                "Pasture" = "#FFF700",  
                                "Urban" = "gray"))  +
  xlab("") +
  scale_y_continuous(breaks = c(-1.2, -0.6, 0, 0.6),
                     limits = c(-1.7, 1.2)) +
  ylab("Coefficient estimate") +
  #ylim(-1.7, 0.75)  +
  theme_bw() + 
  facet_grid(Metric ~ Sample) +
  coord_cartesian(clip = "off") +
  theme(axis.text.x = element_blank(),#text(angle = 60, size = 20, hjust = 1.1),
        axis.text.y = element_text(size = 20),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 25, vjust = 2.2) ,
        #axis.ticks.y = element_line(margin = margin(t = 10)),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        plot.title = element_text(size = 17, face = "bold"),
        strip.background = element_rect(fill = "grey"),
        strip.text.y = element_text(size = 17, face = "bold", vjust = 3),
        strip.text.x = element_text(size = 20, face = "bold"),
        strip.background.y = element_blank(),
        strip.background.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border =  element_blank(),
        #strip.placement = "outside",
        panel.background = element_blank()) +
  geom_text(inherit.aes = FALSE,
            data    = dat_text2,
            mapping = aes(x = -Inf, y = -Inf, label = label),
            hjust   = -1,
            vjust   = -12,
            size = 9,
            fontface = "bold",
            col = "black") +
  guides(col = guide_legend(ncol = 3))


png("Figures/Publication/Main/ED_2 - Tropics vs Temperate.png",
    width = 798/1.25, height = 1214/1.25)
richness  
dev.off()


