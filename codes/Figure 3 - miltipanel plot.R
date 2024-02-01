library(dplyr)
library(tidyr)
library(cowplot)
library(patchwork)

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
  
  analysis_data$LandUse[analysis_data$LandUse %in% c("Pristine primary",
                                                     "Disturbed primary")] <- "Primary vegetation"
  
  analysis_data$LandUse[analysis_data$LandUse %in% c("MSV","ISV", "YSV")] <- "Secondary vegetation"
                        
  analysis_data$LandUse[analysis_data$LandUse %in% c("Plantation forest",
                                                     "Pasture",
                                                     "Cropland")] <- "Agriculture"
  
  analysis_data$LandUse[grep("Urban", analysis_data$LandUse)] <- "Urban"
  analysis_data$LandUse <- factor(analysis_data$LandUse, levels = c("Primary vegetation",
                                                                    "Secondary vegetation",
                                                                    "Agriculture",
                                                                    "Urban"))
  
  
  # Scale and transform columns
  # SQRT transform
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
    summarise(n_NAs = length(which(is.na(Resistance_guild))),
              n_sites = n()) %>%
    mutate(prop_NAs = n_NAs/n_sites) %>%
    ungroup() %>%
    filter(prop_NAs < 1) %>%
    pull(SS) 
  
  guild_results <- all_results %>% filter(SS %in% Studies_to_model)
  guild_results$FRich_guild[is.na(guild_results$FRich_guild)] <- 0
  guild_results$FRich_guild <-  as.numeric(scale(sqrt(guild_results$FRich_guild)))
  
  guild_results$Resistance_guild[is.na(guild_results$Resistance_guild)] <- 50
  guild_results$Resistance_guild <- scale(guild_results$Resistance_guild)
  
  
  
  guild_results$Red_guild[is.na(guild_results$Red_guild)] <- 0
  guild_results$Red_guild[guild_results$Red_guild < 0] <- 0
  guild_results$Red_guild <- as.numeric(scale(sqrt(guild_results$Red_guild)))
  
  
  return(guild_results)
  
}


#################### TRAIT WEIGHTED #############################################
## Get the trait weighted dataset which is ready to model
analysis_data_tw <- read.csv("data/analysis_data_trait_full.csv") %>% sort_dataset()
analysis_data_tw$Resistance_spline <- scale(analysis_data_tw$Resistance_spline)

frug_results <- read.csv("data/frug_all_results_trait.csv") %>%
  dplyr::select(SSBS, Resistance_spline, Red, FRich) %>% fix_NAs()

inv_results <- read.csv("data/inv_all_results_trait.csv") %>%
  dplyr::select(SSBS, Resistance_spline, Red, FRich) %>% fix_NAs()

omni_results <- read.csv("data/omniv_all_results_trait.csv") %>%
  dplyr::select(SSBS, Resistance_spline, Red, FRich) %>% fix_NAs()

gran_results <- read.csv("data/graniv_all_results_trait.csv") %>%
  dplyr::select(SSBS, Resistance_spline, Red, FRich) %>% fix_NAs()



########### Models ###################
all_mod_T_Resist <- lmer(Resistance_spline ~ LandUse + (1|SSB) + (1|SS), data = analysis_data_tw)
frug_mod_T_Resist <- lmer(Resistance_guild ~ LandUse + (1|SSB) + (1|SS), data = frug_results_tw)
inv_mod_T_Resist <- lmer(Resistance_guild ~ LandUse + (1|SSB) + (1|SS), data = inv_results_tw)
gran_mod_T_Resist <- lmer(Resistance_guild ~ LandUse + (1|SSB) + (1|SS), data = gran_results_tw)
omni_mod_T_Resist <- lmer(Resistance_guild ~ LandUse + (1|SSB) + (1|SS), data = omni_results_tw)

all_mod_T_Resist <- as.data.frame(summary(all_mod_T_Resist)$coef)
frug_T_Resist <- as.data.frame(summary(frug_mod_T_Resist)$coef)
inv_T_Resist <- as.data.frame(summary(inv_mod_T_Resist)$coef)
gran_T_Resist <- as.data.frame(summary(gran_mod_T_Resist)$coef)
omni_T_Resist <- as.data.frame(summary(omni_mod_T_Resist)$coef)

all_T_Resist <- do.call(rbind, list(all_mod_T_Resist, inv_T_Resist, frug_T_Resist, gran_T_Resist, omni_T_Resist))




all_T_Resist$Sample <- c(rep("All species", 4), rep("Invertivore", 4), rep("Frugivore", 4), rep("Granivore", 4), rep("Generalist", 4))
all_T_Resist$LandUse <- rep(c("Primary vegetation",
                            "Secondary vegetation",
                            "Agriculture",
                            "Urban"), 5)  


colnames(all_T_Resist)[2] <- "se"
all_T_Resist$Estimate[all_T_Resist$LandUse == "Primary vegetation"] <- 0
all_T_Resist$LandUse <- factor(all_T_Resist$LandUse, levels = c("Primary vegetation",
                                                              "Secondary vegetation",
                                                              "Agriculture",
                                                              "Urban"))

all_T_Resist$Sample <- factor(all_T_Resist$Sample, levels = c("All species", "Generalist", "Granivore", "Frugivore", "Invertivore"))
all_T_Resist$Scenario <- "Trait-based"


#################################################################################################################################
######################### RARITY BASED #####################################
analysis_data_ab <- read.csv("data/analysis_data_abundance_0123.csv") %>% sort_dataset()
analysis_data_ab$Resistance_spline <- scale(analysis_data_ab$Resistance_spline)

frug_results_ab <- read.csv("data/frug_resistance_abundance.csv") %>%
  dplyr::select(SSBS, Resistance_spline, Red, FRich) %>% fix_NAs()

inv_results_ab <- read.csv("data/invertivore_resistance_abundance.csv") %>%
  dplyr::select(SSBS, Resistance_spline, Red, FRich) %>% fix_NAs()

omni_results_ab <- read.csv("data/omniv_resistance_abundance.csv") %>%
  dplyr::select(SSBS, Resistance_spline, Red, FRich) %>% fix_NAs()

gran_results_ab <- read.csv("data/graniv_resistance_abundance.csv") %>%
  dplyr::select(SSBS, Resistance_spline, Red, FRich) %>% fix_NAs()



########### Models ###################
all_mod_R_Resist <- lmer(Resistance_spline ~ LandUse + (1|SSB) + (1|SS), data = analysis_data_ab)
frug_mod_R_Resist <- lmer(Resistance_guild ~ LandUse + (1|SSB) + (1|SS), data = frug_results_ab)
inv_mod_R_Resist <- lmer(Resistance_guild ~ LandUse + (1|SSB) + (1|SS), data = inv_results_ab)
gran_mod_R_Resist <- lmer(Resistance_guild ~ LandUse + (1|SSB) + (1|SS), data = gran_results_ab)
omni_mod_R_Resist <- lmer(Resistance_guild ~ LandUse + (1|SSB) + (1|SS), data = omni_results_ab)

all_mod_R_Resist <- as.data.frame(summary(all_mod_R_Resist)$coef)
frug_R_Resist <- as.data.frame(summary(frug_mod_R_Resist)$coef)
inv_R_Resist <- as.data.frame(summary(inv_mod_R_Resist)$coef)
gran_R_Resist <- as.data.frame(summary(gran_mod_R_Resist)$coef)
omni_R_Resist <- as.data.frame(summary(omni_mod_R_Resist)$coef)

all_R_Resist <- do.call(rbind, list(all_mod_R_Resist, inv_R_Resist, frug_R_Resist, gran_R_Resist, omni_R_Resist))




all_R_Resist$Sample <- c(rep("All species", 4), rep("Invertivore", 4), rep("Frugivore", 4), rep("Granivore", 4), rep("Generalist", 4))
all_R_Resist$LandUse <- rep(c("Primary vegetation",
                              "Secondary vegetation",
                              "Agriculture",
                              "Urban"), 5)  


colnames(all_R_Resist)[2] <- "se"
all_R_Resist$Estimate[all_R_Resist$LandUse == "Primary vegetation"] <- 0
all_R_Resist$LandUse <- factor(all_R_Resist$LandUse, levels = c("Primary vegetation",
                                                                "Secondary vegetation",
                                                                "Agriculture",
                                                                "Urban"))

all_R_Resist$Sample <- factor(all_R_Resist$Sample, levels = c("All species", "Generalist", "Granivore", "Frugivore", "Invertivore"))
all_R_Resist$Scenario <- "Rarity-based"
#################################################################################################################################

all_Resist <- rbind(all_T_Resist, all_R_Resist)
all_Resist$Scenario <- factor(all_Resist$Scenario, levels = c("Trait-based", "Rarity-based"))


dat_text <- data.frame(
  label = c("All species",
              "Invertivore",
              "Frugivore",
              "Granivore",
              "Generalist"),
  Sample  = c("All species",
               "Invertivore",
               "Frugivore",
               "Granivore",
               "Generalist" 
))

dat_text$Sample <- factor(dat_text$Sample, levels = c("All species",
                                                      "Invertivore",
                                                      "Frugivore",
                                                      "Granivore",
                                                      "Generalist"))

Resistance <- ggplot(all_Resist %>% filter(Scenario == "Trait-based")) +
  geom_point(aes(y = Estimate , x = LandUse, col = LandUse, shape = Scenario), position = position_dodge(0.6), size = 13) +
  geom_errorbar(aes(x = LandUse, ymin = (Estimate - 1.96*se), ymax = (Estimate + 1.96*se), col = LandUse, group = Scenario), width=.4, size = 2, position = position_dodge(0.6)) +
  geom_hline(yintercept = 0, col = "black", cex = 0.5, linetype = "dashed") +
  ggtitle("c") +
  xlab("") +
  ylab("Coefficient estimate") +
  ylim(-1.0, 0.45)  +
  theme_bw() + 
  geom_text(inherit.aes = FALSE,
            data    = dat_text,
            mapping = aes(x = -Inf, y = -Inf, label = label),
            hjust   = -0.1,
            vjust   = -15,
            size = 11.25,
            col = "black",
            fontface = "bold") +
 facet_wrap(~Sample, nrow = 1) +
  scale_color_manual(values = c("Primary vegetation" = "red", 
                                "Secondary vegetation" = "#7EB1D8", 
                                "Young secondary" = "#A4C2FA", 
                                "Agriculture" = "gold", 
                                "Urban" = "darkgray"))  +
  scale_shape_manual(values = c("Trait-based" = 19, "Rarity-based" = 17)) +
  guides(shape=FALSE) +
  geom_hline(yintercept = 0,color = "red", cex = 0.5, linetype = "dashed") +
  theme_minimal() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  theme(axis.text.x = element_blank(),#(angle = 60, size = 27, hjust = 1.1),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 40, vjust = 3.4) ,
        axis.ticks = element_line(size = 1),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 60, face = "bold", vjust = -1, hjust = -0.0625),
        strip.text = element_text(size = 30, hjust = -5))
        

Resistance


Resistance_shape_leg <- ggplot(all_Resist) +
  geom_point(aes(y = Estimate , x = LandUse, col = LandUse, shape = Scenario), position = position_dodge(0.6), size = 7) +
  scale_shape_manual(values = c("Trait-based" = 19, "Rarity-based" = 17)) +
  guides(color=FALSE) +
  theme(legend.text = element_text(size = 30),
        legend.title =  element_blank())



shape_legend <- get_legend(Resistance_shape_leg + 
                             theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 25),
        legend.key.height = unit(1.05, "cm")))

ggdraw(Resistance) +
  draw_plot(shape_legend, .5, .7, .875, -0.95)
################################################################################
### Panel a and b for declines ##

################################################################################
library(lme4)
library(lmerTest)
library(tidyverse)
library(ggpubr)
library(glmmTMB)
library(ggdist)
library(ggridges)
library(DHARMa)
library(tidyverse)
library(dplyr)
library(tidyr)
library(magrittr)
library(plotrix)
library(ggdist)

results_list <- readRDS("data/simulations_trait.rds")
analysis_data<- read.csv("data/analysis_data_trait_0123.csv") %>% 
  
  filter(!SS %in% c("LK1_2009__Hayward 1", # not entire community ##Patrick
                    "GN1_2010__Hvenegaard 1")) %>%# pseudo replicate=
  #!LandUse == "Secondary vegetation") %>%
  
  distinct(SSBS, .keep_all = TRUE) #%>%
#filter(SS != "2020_Lees_RAS") #%>%
#filter(Species_richness > 2)


analysis_data$LandUse[analysis_data$LandUse == "Primary minimal"] <- "Pristine primary"
diversity <- readRDS("data/diversity_combined_sites_Aug2022.rds")

pm_sites <- analysis_data %>% filter(LandUse == "Pristine primary") %>% pull(SSBS)
dp_sites <- analysis_data %>% filter(LandUse == "Disturbed primary") %>% pull(SSBS)
sv_sites <- analysis_data %>% filter(LandUse == "Secondary vegetation") %>% pull(SSBS)
plf_sites <- analysis_data %>% filter(LandUse == "Plantation forest") %>% pull(SSBS)
pst_sites <- analysis_data %>% filter(LandUse == "Pasture") %>% pull(SSBS)
crp_sites <- analysis_data %>% filter(LandUse == "Cropland") %>% pull(SSBS)
urb_sites <- analysis_data %>% filter(LandUse == "Urban") %>% pull(SSBS)


model_list <- list()
site_av_list <- list()
study_list <- list()

for (i in c(1:length(results_list))) {
  #print(paste("i:", i))
  study_list <- results_list[i][[1]]
  study_name <- names(results_list)[i]
  pm <- 0
  dp <- 0
  sv <- 0
  plf <- 0
  pst <- 0
  crp <- 0
  urb  <- 0
  modelled_data_pm <- NULL
  modelled_data_dp <- NULL
  modelled_data_sv <- NULL
  modelled_data_plf <- NULL
  modelled_data_pst <- NULL
  modelled_data_crp <- NULL
  modelled_data_urb <- NULL
  
  for (j in c(1:length(study_list))) try({
    #print(j)
    
    site_name <- names(study_list)[j]
    
    LUType <- analysis_data$LandUse[analysis_data$SSBS == site_name]
    
    site_i <- study_list[[j]]
    site_i$Rich_TPD <- as.data.frame(site_i$Rich_TPD)
    site_i$Rich_TPD$Rich_TPD_j
    
    site_i$Red_TPD <- as.data.frame(site_i$Red_TPD)
    site_i$Rich_TPD[nrow(site_i$Rich_TPD) + 1,] <- 0
    site_i$Red_TPD[nrow(site_i$Red_TPD) + 1,] <- 0
    
    site_av <- data.frame(FD = as.numeric(rowMeans(site_i$Rich_TPD)), 
                          Red = as.numeric(rowMeans(site_i$Red_TPD)), 
                          iteration = seq(100, 0, length.out = nrow(site_i$Rich_TPD)))
    
    site_av$FD_s <- site_av$FD/max(site_av$FD)
    site_av$Red_s <- site_av$Red/max(site_av$Red)
    site_av$Species_remaining <- (nrow(site_av)-1):0
    
    ssbs_mod <- smooth.spline(site_av$FD_s, site_av$Species_remaining)
    
    newdata <- data.frame(FD_s = seq(1, 0, length.out = 26))
    
    predicted <- predict(ssbs_mod, newdata$FD_s)
    
    modelled_data <- data.frame(iteration = predicted$y, predicted_FD = predicted$x)
    
    
    
    if(LUType %in% c("Pristine primary", "Disturbed primary")) try({
      pm <- pm + 1
      ifelse(pm == 1,
             modelled_data_pm <- modelled_data,
             modelled_data_pm <- cbind(modelled_data_pm, modelled_data$predicted_FD))
    })
    # 
    # if(LUType == "Disturbed primary") try({
    #   dp <- dp + 1
    #   ifelse(dp == 1,
    #          modelled_data_dp <- modelled_data,
    #          modelled_data_dp <- cbind(modelled_data_dp, modelled_data$predicted_FD))
    # })
    
    if(LUType == "Secondary vegetation") try({
      sv <- sv + 1
      ifelse(sv == 1,
             modelled_data_sv <- modelled_data,
             modelled_data_sv <- cbind(modelled_data_sv, modelled_data$predicted_FD))
    })
    
    # if(LUType == "Plantation forest") try({
    #   plf <- plf + 1
    #   ifelse(plf == 1,
    #          modelled_data_plf <- modelled_data,
    #          modelled_data_plf <- cbind(modelled_data_plf, modelled_data$predicted_FD))
    # })
    
    if(LUType %in% c("Cropland", "Pasture", "Plantation forest")) try({
      pst <- pst + 1
      ifelse(pst == 1,
             modelled_data_pst <- modelled_data,
             modelled_data_pst <- cbind(modelled_data_pst, modelled_data$predicted_FD))
    })
    
    # if(LUType == "Cropland") try({
    #   crp <- crp + 1
    #   ifelse(crp == 1,
    #          modelled_data_crp <- modelled_data,
    #          modelled_data_crp <- cbind(modelled_data_crp, modelled_data$predicted_FD))
    # })
    if(LUType == "Urban") try({
      urb <- urb + 1
      ifelse(urb == 1,
             modelled_data_urb <- modelled_data,
             modelled_data_urb <- cbind(modelled_data_urb, modelled_data$predicted_FD))
    })
    
    
    
    ## At the end of each study put into a new thing the average means of the modelled data
    if(j == length(study_list)) {
      
      ifelse(ncol(modelled_data_pm) > 2,
             modelled_data_pm_av <-  data.frame(predicted_FD = modelled_data_pm[,2], species_removed = rowMeans(modelled_data_pm[,-2])),
             modelled_data_pm_av <- data.frame(predicted_FD = modelled_data_pm[,2], species_removed = modelled_data_pm[,1]))
      # 
      # ifelse(ncol(modelled_data_dp) > 2,
      #        modelled_data_dp_av <-  data.frame(predicted_FD = modelled_data_dp[,2], species_removed = rowMeans(modelled_data_dp[,-2])),
      #        modelled_data_dp_av <- data.frame(predicted_FD = modelled_data_dp[,2], species_removed = modelled_data_dp[,1]))
      # 
      
      ifelse(ncol(modelled_data_sv) > 2,
             modelled_data_sv_av <-  data.frame(predicted_FD = modelled_data_sv[,2], species_removed = rowMeans(modelled_data_sv[,-2])),
             modelled_data_sv_av <- data.frame(predicted_FD = modelled_data_sv[,2], species_removed = modelled_data_sv[,1]))
      
      # ifelse(ncol(modelled_data_plf) > 2,
      #        modelled_data_plf_av <-  data.frame(predicted_FD = modelled_data_plf[,2], species_removed = rowMeans(modelled_data_plf[,-2])),
      #        modelled_data_plf_av <- data.frame(predicted_FD = modelled_data_plf[,2], species_removed = modelled_data_plf[,1]))
      # 
      ifelse(ncol(modelled_data_pst) > 2,
             modelled_data_pst_av <-  data.frame(predicted_FD = modelled_data_pst[,2], species_removed = rowMeans(modelled_data_pst[,-2])),
             modelled_data_pst_av <- data.frame(predicted_FD = modelled_data_pst[,2], species_removed = modelled_data_pst[,1]))
      
      # ifelse(ncol(modelled_data_crp) > 2,
      #        modelled_data_crp_av <-  data.frame(predicted_FD = modelled_data_crp[,2], species_removed = rowMeans(modelled_data_crp[,-2])),
      #        modelled_data_crp_av <- data.frame(predicted_FD = modelled_data_crp[,2], species_removed = modelled_data_crp[,1]))
      # 
      ifelse(ncol(modelled_data_urb) > 2,
             modelled_data_urb_av <- data.frame(predicted_FD = modelled_data_urb[,2], species_removed = rowMeans(modelled_data_urb[,-2])),
             modelled_data_urb_av <- data.frame(predicted_FD = modelled_data_urb[,2], species_removed = modelled_data_urb[,1]))
      
      
      ## Put the average means into a data.frame
      site_av_list[[paste(study_name, "- Pristine primary")]] <- modelled_data_pm_av
      #site_av_list[[paste(study_name, "- Disturbed primary")]] <- modelled_data_dp_av
      site_av_list[[paste(study_name, "- Secondary vegetation")]] <- modelled_data_sv_av
      #site_av_list[[paste(study_name, "- Plantation forest")]] <- modelled_data_plf_av
      site_av_list[[paste(study_name, "- Pasture")]] <- modelled_data_pst_av
      # site_av_list[[paste(study_name, "- Cropland")]] <- modelled_data_crp_av
      site_av_list[[paste(study_name, "- Urban")]] <- modelled_data_urb_av
    }
    
  }, silent = TRUE)
  #print("done") 
}

pm <- 0
dp <- 0
sv <- 0
plf <- 0
pst <- 0
crp <- 0
urb  <- 0

## bind the predicted ones into dfs of landuse
for (i in 1:length(site_av_list)) {
  #print(i)
  if(grepl("Pristine primary",names(site_av_list)[i])) {
    print(paste("i:", i))
    pm <- pm + 1
    print(paste("pm:", pm)) 
    ifelse(pm == 1,
           pm_df <-  site_av_list[[i]],
           pm_df <-  cbind(pm_df, site_av_list[[i]]$species_removed))
  }
  # if(grepl("Disturbed primary",names(site_av_list)[i])) {
  #   dp <- dp + 1
  #   ifelse(dp == 1,
  #          dp_df <-  site_av_list[[i]],
  #          dp_df <-  cbind(dp_df, site_av_list[[i]]$species_removed))
  # }
  if(grepl("Secondary vegetation",names(site_av_list)[i])) {
    sv <- sv + 1
    ifelse(sv == 1,
           sv_df <-  site_av_list[[i]],
           sv_df <-  cbind(sv_df, site_av_list[[i]]$species_removed))
  }
  # if(grepl("Plantation forest",names(site_av_list)[i])) {
  #   plf <- plf + 1
  #   ifelse(plf == 1,
  #          plf_df <-  site_av_list[[i]],
  #          plf_df <-  cbind(plf_df, site_av_list[[i]]$species_removed))
  # }
  if(grepl("Pasture",names(site_av_list)[i])) {
    pst <- pst + 1
    ifelse(pst == 1,
           pst_df <-  site_av_list[[i]],
           pst_df <-  cbind(pst_df, site_av_list[[i]]$species_removed))
  }
  # if(grepl("Cropland",names(site_av_list)[i])) {
  #   crp <- crp + 1
  #   ifelse(crp == 1,
  #          crp_df <-  site_av_list[[i]],
  #          crp_df <-  cbind(crp_df, site_av_list[[i]]$species_removed))
  # }
  if(grepl("Urban",names(site_av_list)[i])) {
    urb <- urb + 1
    ifelse(urb == 1,
           urb_df <-  site_av_list[[i]],
           urb_df <-  cbind(urb_df, site_av_list[[i]]$species_removed))
  }
}

pm_df_mean <- data.frame(FD = pm_df$predicted_FD,
                         SR = as.numeric(rowMeans(pm_df[,-c(1,68)])),
                         se = as.numeric(apply(pm_df[,-c(1,68)], 1, std.error)),
                         LandUse = "Pristine primary vegetation") %>%
  dplyr::mutate(uppr = SR + (se),
                lwr = SR - (se))


# dp_df_mean <- data.frame(FD = dp_df$predicted_FD,
#                          SR = as.numeric(rowMeans(dp_df[,-1])),
#                          se = as.numeric(apply(dp_df[,-1], 1, std.error)),
#                          LandUse = "Disturbed primary vegetation") %>%
#   dplyr::mutate(uppr = SR + (se),
#                 lwr = SR - (se))

sv_df_mean <- data.frame(FD = sv_df$predicted_FD,
                         SR = as.numeric(rowMeans(sv_df[,-1])),
                         se = as.numeric(apply(sv_df[,-1], 1, std.error)),
                         LandUse = "Secondary vegetation") %>%
  dplyr::mutate(uppr = SR + (se),
                lwr = SR - (se))

# plf_df_mean <- data.frame(FD = plf_df$predicted_FD,
#                           SR = as.numeric(rowMeans(plf_df[,-1])),
#                           se = as.numeric(apply(plf_df[,-1], 1, std.error)),
#                           LandUse = "Plantation forest") %>%
#   dplyr::mutate(uppr = SR + (se),
#                 lwr = SR - (se))
# 
pst_df_mean <- data.frame(FD = pst_df$predicted_FD,
                          SR = as.numeric(rowMeans(pst_df[,-1])),
                          se = as.numeric(apply(pst_df[,-1], 1, std.error)),
                          LandUse = "Pasture") %>%
  dplyr::mutate(uppr = SR + (se),
                lwr = SR - (se))

# crp_df_mean <- data.frame(FD = crp_df$predicted_FD,
#                           SR = as.numeric(rowMeans(crp_df[,-1])),
#                           se = as.numeric(apply(crp_df[,-1], 1, std.error)),
#                           LandUse = "Cropland") %>%
#   dplyr::mutate(uppr = SR + (se),
#                 lwr = SR - (se))
# 
urb_df_mean <- data.frame(FD = urb_df$predicted_FD,
                          SR = as.numeric(rowMeans(urb_df[,-1])),
                          se = as.numeric(apply(urb_df[,-1], 1, std.error)),
                          LandUse = "Urban") %>%
  dplyr::mutate(uppr = SR + (se),
                lwr = SR - (se))



df_all <- rbind(pm_df_mean,
                # dp_df_mean,
                sv_df_mean,
                # plf_df_mean,
                pst_df_mean,
                #crp_df_mean,
                urb_df_mean)

df_all$LandUse[which(df_all$LandUse == "Pristine primary vegetation")] <- "Primary vegetation"
df_all$LandUse[which(df_all$LandUse == "Pasture")] <- "Agriculture"

df_all$LandUse <- factor(df_all$LandUse, levels = c("Primary vegetation",
                                                    #"Disturbed primary vegetation",
                                                    "Secondary vegetation",
                                                    #"Plantation forest",
                                                    "Agriculture",
                                                    #"Cropland",
                                                    "Urban"))    

df_all1 <- as.data.frame(df_all %>% group_by(LandUse) %>% 
                           mutate(SR = abs(SR + (max(df_all$SR) - max(SR))-max(df_all$SR))) %>%
                           mutate(lwr = abs(lwr + (max(df_all$SR) - max(lwr))-max(df_all$SR))) %>%
                           mutate(uppr = abs(uppr + (max(df_all$SR) - max(uppr))-max(df_all$SR))))



df_all1$FD <- df_all$FD*100

declines_tw <- ggplot(df_all1) +
  geom_line(aes(x = FD, y = SR, colour = LandUse), size = 3) +
  geom_ribbon(aes(x=FD, y=NULL, ymin=lwr, ymax=uppr, fill = LandUse), show.legend=FALSE,alpha=0.1) +
  scale_color_manual(labels = c(
    "Primary vegetation" = bquote("Primary vegetation (" ~ italic(n) ~ " = 458)"),
    "Secondary vegetation" = bquote("Secondary vegetation (" ~ italic(n) ~ " = 252)"),
    "Agriculture" = bquote("Agriculture (" ~ italic(n) ~ " = 509)"),
    "Urban" = bquote("Urban (" ~ italic(n) ~ " = 107)")
  ),
                     name = "Land-use: ", 
                     values = c("red", "#4EADC7", "gold", "grey")) +
  scale_fill_manual(values = c("red", "#4EADC7", "gold", "grey")) +
  guides(color = guide_legend(# title.hjust = 1, # adjust title if needed
    label.position = "left")) +
  ggtitle("b") +
  ylab("Number of species removed (Traits)") +
  xlab("FD (Richness) (%)") +
  xlim(0,100) +
  #ylim(0,140) +
  theme_classic() +
  coord_flip() +
  #scale_colour_discrete("Land use type") +
  theme(legend.position = c(0.7,0.825), 
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_rect(fill = "transparent", color = "transparent"),#alpha("white", 0.0)),
        legend.key.width = unit(17, "mm"),
        legend.text = element_text(size = 35),
        legend.title = element_blank(),
        axis.text = element_text(size = 35),
        axis.title.y = element_text(size = 40, vjust = 3.4),
        axis.title.x = element_text(size = 40, vjust = -1),
        plot.margin = margin(0,1,0,2, "cm"),
        plot.title = element_text(size = 60, face = "bold", vjust = 1, hjust = -0.0625)
  ) #+ guides(color=guide_legend(nrow=1,byrow=FALSE))


################################################################################
library(lme4)
library(lmerTest)
library(tidyverse)
library(ggpubr)
library(glmmTMB)
library(ggdist)
library(ggridges)
library(DHARMa)
library(tidyverse)
library(dplyr)
library(tidyr)
library(magrittr)
library(plotrix)
library(ggdist)

results_list <- readRDS("data/simulations_abundance.rds")
analysis_data<- read.csv("data/analysis_data_abundance_0123.csv") %>% 
  
  filter(!SS %in% c("LK1_2009__Hayward 1", # not entire community ##Patrick
                    "GN1_2010__Hvenegaard 1")) %>%# pseudo replicate=
  #!LandUse == "Secondary vegetation") %>%
  
  distinct(SSBS, .keep_all = TRUE) #%>%
#filter(SS != "2020_Lees_RAS") #%>%
#filter(Species_richness > 2)


analysis_data$LandUse[analysis_data$LandUse == "Primary minimal"] <- "Pristine primary"
diversity <- readRDS("data/diversity_combined_sites_Aug2022.rds")

pm_sites <- analysis_data %>% filter(LandUse == "Pristine primary") %>% pull(SSBS)
dp_sites <- analysis_data %>% filter(LandUse == "Disturbed primary") %>% pull(SSBS)
sv_sites <- analysis_data %>% filter(LandUse == "Secondary vegetation") %>% pull(SSBS)
plf_sites <- analysis_data %>% filter(LandUse == "Plantation forest") %>% pull(SSBS)
pst_sites <- analysis_data %>% filter(LandUse == "Pasture") %>% pull(SSBS)
crp_sites <- analysis_data %>% filter(LandUse == "Cropland") %>% pull(SSBS)
urb_sites <- analysis_data %>% filter(LandUse == "Urban") %>% pull(SSBS)


model_list <- list()
site_av_list <- list()
study_list <- list()

for (i in c(1:length(results_list))) {
  #print(paste("i:", i))
  study_list <- results_list[i][[1]]
  study_name <- names(results_list)[i]
  pm <- 0
  dp <- 0
  sv <- 0
  plf <- 0
  pst <- 0
  crp <- 0
  urb  <- 0
  modelled_data_pm <- NULL
  modelled_data_dp <- NULL
  modelled_data_sv <- NULL
  modelled_data_plf <- NULL
  modelled_data_pst <- NULL
  modelled_data_crp <- NULL
  modelled_data_urb <- NULL
  
  for (j in c(1:length(study_list))) try({
    #print(j)
    
    site_name <- names(study_list)[j]
    
    LUType <- analysis_data$LandUse[analysis_data$SSBS == site_name]
    
    site_i <- study_list[[j]]
    site_i$Rich_TPD <- as.data.frame(site_i$Rich_TPD)
    site_i$Rich_TPD$Rich_TPD_j
    
    site_i$Red_TPD <- as.data.frame(site_i$Red_TPD)
    site_i$Rich_TPD[nrow(site_i$Rich_TPD) + 1,] <- 0
    site_i$Red_TPD[nrow(site_i$Red_TPD) + 1,] <- 0
    
    site_av <- data.frame(FD = as.numeric(rowMeans(site_i$Rich_TPD)), 
                          Red = as.numeric(rowMeans(site_i$Red_TPD)), 
                          iteration = seq(100, 0, length.out = nrow(site_i$Rich_TPD)))
    
    site_av$FD_s <- site_av$FD/max(site_av$FD)
    site_av$Red_s <- site_av$Red/max(site_av$Red)
    site_av$Species_remaining <- (nrow(site_av)-1):0
    
    ssbs_mod <- smooth.spline(site_av$FD_s, site_av$Species_remaining)
    
    newdata <- data.frame(FD_s = seq(1, 0, length.out = 26))
    
    predicted <- predict(ssbs_mod, newdata$FD_s)
    
    modelled_data <- data.frame(iteration = predicted$y, predicted_FD = predicted$x)
    
    
    
    if(LUType %in% c("Pristine primary", "Disturbed primary")) try({
      pm <- pm + 1
      ifelse(pm == 1,
             modelled_data_pm <- modelled_data,
             modelled_data_pm <- cbind(modelled_data_pm, modelled_data$predicted_FD))
    })
    # 
    # if(LUType == "Disturbed primary") try({
    #   dp <- dp + 1
    #   ifelse(dp == 1,
    #          modelled_data_dp <- modelled_data,
    #          modelled_data_dp <- cbind(modelled_data_dp, modelled_data$predicted_FD))
    # })
    
    if(LUType == "Secondary vegetation") try({
      sv <- sv + 1
      ifelse(sv == 1,
             modelled_data_sv <- modelled_data,
             modelled_data_sv <- cbind(modelled_data_sv, modelled_data$predicted_FD))
    })
    
    # if(LUType == "Plantation forest") try({
    #   plf <- plf + 1
    #   ifelse(plf == 1,
    #          modelled_data_plf <- modelled_data,
    #          modelled_data_plf <- cbind(modelled_data_plf, modelled_data$predicted_FD))
    # })
    
    if(LUType %in% c("Cropland", "Pasture", "Plantation forest")) try({
      pst <- pst + 1
      ifelse(pst == 1,
             modelled_data_pst <- modelled_data,
             modelled_data_pst <- cbind(modelled_data_pst, modelled_data$predicted_FD))
    })
    
    # if(LUType == "Cropland") try({
    #   crp <- crp + 1
    #   ifelse(crp == 1,
    #          modelled_data_crp <- modelled_data,
    #          modelled_data_crp <- cbind(modelled_data_crp, modelled_data$predicted_FD))
    # })
    if(LUType == "Urban") try({
      urb <- urb + 1
      ifelse(urb == 1,
             modelled_data_urb <- modelled_data,
             modelled_data_urb <- cbind(modelled_data_urb, modelled_data$predicted_FD))
    })
    
    
    
    ## At the end of each study put into a new thing the average means of the modelled data
    if(j == length(study_list)) {
      
      ifelse(ncol(modelled_data_pm) > 2,
             modelled_data_pm_av <-  data.frame(predicted_FD = modelled_data_pm[,2], species_removed = rowMeans(modelled_data_pm[,-2])),
             modelled_data_pm_av <- data.frame(predicted_FD = modelled_data_pm[,2], species_removed = modelled_data_pm[,1]))
      # 
      # ifelse(ncol(modelled_data_dp) > 2,
      #        modelled_data_dp_av <-  data.frame(predicted_FD = modelled_data_dp[,2], species_removed = rowMeans(modelled_data_dp[,-2])),
      #        modelled_data_dp_av <- data.frame(predicted_FD = modelled_data_dp[,2], species_removed = modelled_data_dp[,1]))
      # 
      
      ifelse(ncol(modelled_data_sv) > 2,
             modelled_data_sv_av <-  data.frame(predicted_FD = modelled_data_sv[,2], species_removed = rowMeans(modelled_data_sv[,-2])),
             modelled_data_sv_av <- data.frame(predicted_FD = modelled_data_sv[,2], species_removed = modelled_data_sv[,1]))
      
      # ifelse(ncol(modelled_data_plf) > 2,
      #        modelled_data_plf_av <-  data.frame(predicted_FD = modelled_data_plf[,2], species_removed = rowMeans(modelled_data_plf[,-2])),
      #        modelled_data_plf_av <- data.frame(predicted_FD = modelled_data_plf[,2], species_removed = modelled_data_plf[,1]))
      # 
      ifelse(ncol(modelled_data_pst) > 2,
             modelled_data_pst_av <-  data.frame(predicted_FD = modelled_data_pst[,2], species_removed = rowMeans(modelled_data_pst[,-2])),
             modelled_data_pst_av <- data.frame(predicted_FD = modelled_data_pst[,2], species_removed = modelled_data_pst[,1]))
      
      # ifelse(ncol(modelled_data_crp) > 2,
      #        modelled_data_crp_av <-  data.frame(predicted_FD = modelled_data_crp[,2], species_removed = rowMeans(modelled_data_crp[,-2])),
      #        modelled_data_crp_av <- data.frame(predicted_FD = modelled_data_crp[,2], species_removed = modelled_data_crp[,1]))
      # 
      ifelse(ncol(modelled_data_urb) > 2,
             modelled_data_urb_av <- data.frame(predicted_FD = modelled_data_urb[,2], species_removed = rowMeans(modelled_data_urb[,-2])),
             modelled_data_urb_av <- data.frame(predicted_FD = modelled_data_urb[,2], species_removed = modelled_data_urb[,1]))
      
      
      ## Put the average means into a data.frame
      site_av_list[[paste(study_name, "- Pristine primary")]] <- modelled_data_pm_av
      #site_av_list[[paste(study_name, "- Disturbed primary")]] <- modelled_data_dp_av
      site_av_list[[paste(study_name, "- Secondary vegetation")]] <- modelled_data_sv_av
      #site_av_list[[paste(study_name, "- Plantation forest")]] <- modelled_data_plf_av
      site_av_list[[paste(study_name, "- Pasture")]] <- modelled_data_pst_av
      # site_av_list[[paste(study_name, "- Cropland")]] <- modelled_data_crp_av
      site_av_list[[paste(study_name, "- Urban")]] <- modelled_data_urb_av
    }
    
  }, silent = TRUE)
  #print("done") 
}

pm <- 0
dp <- 0
sv <- 0
plf <- 0
pst <- 0
crp <- 0
urb  <- 0

## bind the predicted ones into dfs of landuse
for (i in 1:length(site_av_list)) {
  print(i)
  if(grepl("Pristine primary",names(site_av_list)[i])) {
    print(paste("i:", i))
    pm <- pm + 1
    print(paste("pm:", pm)) 
    ifelse(pm == 1,
           pm_df <-  site_av_list[[i]],
           pm_df <-  cbind(pm_df, site_av_list[[i]]$species_removed))
  }
  # if(grepl("Disturbed primary",names(site_av_list)[i])) {
  #   dp <- dp + 1
  #   ifelse(dp == 1,
  #          dp_df <-  site_av_list[[i]],
  #          dp_df <-  cbind(dp_df, site_av_list[[i]]$species_removed))
  # }
  if(grepl("Secondary vegetation",names(site_av_list)[i])) {
    sv <- sv + 1
    ifelse(sv == 1,
           sv_df <-  site_av_list[[i]],
           sv_df <-  cbind(sv_df, site_av_list[[i]]$species_removed))
  }
  # if(grepl("Plantation forest",names(site_av_list)[i])) {
  #   plf <- plf + 1
  #   ifelse(plf == 1,
  #          plf_df <-  site_av_list[[i]],
  #          plf_df <-  cbind(plf_df, site_av_list[[i]]$species_removed))
  # }
  if(grepl("Pasture",names(site_av_list)[i])) {
    pst <- pst + 1
    ifelse(pst == 1,
           pst_df <-  site_av_list[[i]],
           pst_df <-  cbind(pst_df, site_av_list[[i]]$species_removed))
  }
  # if(grepl("Cropland",names(site_av_list)[i])) {
  #   crp <- crp + 1
  #   ifelse(crp == 1,
  #          crp_df <-  site_av_list[[i]],
  #          crp_df <-  cbind(crp_df, site_av_list[[i]]$species_removed))
  # }
  if(grepl("Urban",names(site_av_list)[i])) {
    urb <- urb + 1
    ifelse(urb == 1,
           urb_df <-  site_av_list[[i]],
           urb_df <-  cbind(urb_df, site_av_list[[i]]$species_removed))
  }
}

pm_df_mean <- data.frame(FD = pm_df$predicted_FD,
                         SR = as.numeric(rowMeans(pm_df[,-c(1,68)])),
                         se = as.numeric(apply(pm_df[,-c(1,68)], 1, std.error)),
                         LandUse = "Pristine primary vegetation") %>%
  dplyr::mutate(uppr = SR + (se),
                lwr = SR - (se))


# dp_df_mean <- data.frame(FD = dp_df$predicted_FD,
#                          SR = as.numeric(rowMeans(dp_df[,-1])),
#                          se = as.numeric(apply(dp_df[,-1], 1, std.error)),
#                          LandUse = "Disturbed primary vegetation") %>%
#   dplyr::mutate(uppr = SR + (se),
#                 lwr = SR - (se))

sv_df_mean <- data.frame(FD = sv_df$predicted_FD,
                         SR = as.numeric(rowMeans(sv_df[,-1])),
                         se = as.numeric(apply(sv_df[,-1], 1, std.error)),
                         LandUse = "Secondary vegetation") %>%
  dplyr::mutate(uppr = SR + (se),
                lwr = SR - (se))

# plf_df_mean <- data.frame(FD = plf_df$predicted_FD,
#                           SR = as.numeric(rowMeans(plf_df[,-1])),
#                           se = as.numeric(apply(plf_df[,-1], 1, std.error)),
#                           LandUse = "Plantation forest") %>%
#   dplyr::mutate(uppr = SR + (se),
#                 lwr = SR - (se))
# 
pst_df_mean <- data.frame(FD = pst_df$predicted_FD,
                          SR = as.numeric(rowMeans(pst_df[,-1])),
                          se = as.numeric(apply(pst_df[,-1], 1, std.error)),
                          LandUse = "Pasture") %>%
  dplyr::mutate(uppr = SR + (se),
                lwr = SR - (se))

# crp_df_mean <- data.frame(FD = crp_df$predicted_FD,
#                           SR = as.numeric(rowMeans(crp_df[,-1])),
#                           se = as.numeric(apply(crp_df[,-1], 1, std.error)),
#                           LandUse = "Cropland") %>%
#   dplyr::mutate(uppr = SR + (se),
#                 lwr = SR - (se))
# 
urb_df_mean <- data.frame(FD = urb_df$predicted_FD,
                          SR = as.numeric(rowMeans(urb_df[,-1])),
                          se = as.numeric(apply(urb_df[,-1], 1, std.error)),
                          LandUse = "Urban") %>%
  dplyr::mutate(uppr = SR + (se),
                lwr = SR - (se))



df_all <- rbind(pm_df_mean,
                # dp_df_mean,
                sv_df_mean,
                # plf_df_mean,
                pst_df_mean,
                #crp_df_mean,
                urb_df_mean)

df_all$LandUse[which(df_all$LandUse == "Pristine primary vegetation")] <- "Primary vegetation"
df_all$LandUse[which(df_all$LandUse == "Pasture")] <- "Agriculture"

df_all$LandUse <- factor(df_all$LandUse, levels = c("Primary vegetation",
                                                    #"Disturbed primary vegetation",
                                                    "Secondary vegetation",
                                                    #"Plantation forest",
                                                    "Agriculture",
                                                    #"Cropland",
                                                    "Urban"))    

df_all2 <- as.data.frame(df_all %>% group_by(LandUse) %>% 
                           mutate(SR = abs(SR + (max(df_all$SR) - max(SR))-max(df_all$SR))) %>%
                           mutate(lwr = abs(lwr + (max(df_all$SR) - max(lwr))-max(df_all$SR))) %>%
                           mutate(uppr = abs(uppr + (max(df_all$SR) - max(uppr))-max(df_all$SR))))



df_all2$FD <- df_all$FD*100

declines_ab <- ggplot(df_all2) +
  geom_line(aes(x = FD, y = SR, colour = LandUse), size = 2) +
  geom_ribbon(aes(x=FD, y=NULL, ymin=lwr, ymax=uppr, fill = LandUse), show.legend=FALSE,alpha=0.1) +
  scale_color_manual(labels = c(
    "Primary vegetation" = bquote("Primary vegetation (" ~ italic(n) ~ " = 458)"),
    "Secondary vegetation" = bquote("Secondary vegetation (" ~ italic(n) ~ " = 252)"),
    "Agriculture" = bquote("Agriculture (" ~ italic(n) ~ " = 509)"),
    "Urban" = bquote("Urban (" ~ italic(n) ~ " = 107)")
  ),
  name = "Land-use: ", 
  values = c("red", "#4EADC7", "gold", "grey")) +
  scale_fill_manual(values = c("red", "#4EADC7", "gold", "grey")) +
  guides(color = guide_legend(# title.hjust = 1, # adjust title if needed
    label.position = "left")) +
  ggtitle("b") +
  ylab("Number of species removed (Rarity)") +
  xlab("FD (Richness) (%)") +
  xlim(0,100) +
  #ylim(0,140) +
  theme_classic() +
  coord_flip() +
  #scale_colour_discrete("Land use type") +
  theme(legend.position = c(0.5,0.95), 
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_rect(fill = "transparent", color = "transparent"),#alpha("white", 0.0)),
        legend.key.width = unit(12, "mm"),
        legend.text = element_text(size = 30),
        legend.title = element_blank(),
        axis.text = element_text(size = 27),
        axis.title.y = element_text(size = 35, vjust = 3.4),
        axis.title.x = element_text(size = 35, vjust = -1),
        plot.margin = margin(0,1,0,2, "cm"),
        plot.title = element_text(size = 50, face = "bold", vjust = 1, hjust = -0.15)
  ) #+ guides(color=guide_legend(nrow=1,byrow=FALSE))


library(patchwork)
ggarrange(declines_tw, declines_ab)
a <- (plot_spacer() | declines_ab | declines_tw) + plot_layout(widths = c(1, 1000,2000)) 

 a / plot_spacer() / (plot_spacer() + ggdraw(Resistance) +
                                 draw_plot(shape_legend, .5, .7, .875, -1.1) + plot_layout(widths = c(1,45))) + plot_layout(heights = c(2.5,0.2,2.5))

dev.size("px")
png("figures/Figure 4 - LU x Resistance multipanel.png",
    width = 1800/1.25, height = 1316/1.25)
declines_tw / Resistance
dev.off()
 
a /  (plot_spacer() + Resistance + plot_layout(widths = c(1,39))) + plot_layout(heights = c(10000,8500))

 
 
