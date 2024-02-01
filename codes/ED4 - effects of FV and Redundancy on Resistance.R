library(lme4)
library(lmerTest)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(glmmTMB)
library(DHARMa)
library(lme4)
library(lmerTest)
library(tidyverse)
library(ggpubr)
library(glmmTMB)
library(DHARMa)
library(cowplot)

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
  analysis_data$LandUse <- factor(analysis_data$LandUse, levels = c("Pristine primary",
                                                                    "Disturbed primary",
                                                                    "MSV",
                                                                    "ISV",
                                                                    "YSV",
                                                                    "Plantation forest",
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
  
  analysis_data$Resistance_spline <- scale(analysis_data$Resistance_spline)
  analysis_data$Resistance_step <- scale(analysis_data$Resistance_step)
  
  analysis_data$cov_pearsons_rand_red <- as.numeric(scale(analysis_data$cov_pearsons_random_sensitivity_score_redraw))
  analysis_data$cov_spearman_rand_red <- as.numeric(scale(analysis_data$cov_spearman_random_sensitivity_score_redraw))
  
  analysis_data$cov_pearsons_sense_red <- as.numeric(scale(analysis_data$cov_pearsons_sense_redraw))
  analysis_data$cov_spearman_sense_red <- as.numeric(scale(analysis_data$cov_spearman_sense_redraw))
  
  analysis_data$cov_pearsons_abundance_red <- as.numeric(scale(-analysis_data$cov_pearsons_abundance_redraw))
  analysis_data$cov_spearman_abundance_red <- as.numeric(scale(-analysis_data$cov_spearman_abundance_redraw))
  
  
  # save an analysis_data without the RAS dataset (currently slightly erroneous)
  analysis_data_nRAS <- analysis_data %>%
    filter(SS != "2020_Lees_RAS")
  
  analysis_data_nUrbMin <- analysis_data %>%
    filter(LandUseIntensity != "Urban-Minimal use")
  
  analysis_data$RAS <- "no"
  analysis_data$RAS[analysis_data$SS == "2020_Lees_RAS"] <- "yes"
  
  
  analysis_data$Climatic_zone <- "Tropics"
  analysis_data$Climatic_zone[which(abs(analysis_data$Latitude) > 20)] <- "Temperate"
  
  return(analysis_data)
  
}

analysis_data_tw <- read.csv("data/analysis_data_trait_full.csv")  %>% sort_dataset()
analysis_data_ab <- read.csv("data/analysis_data_abundance_full.csv")  %>% sort_dataset()


analysis_data_tw$cov_pearsons_sense_red <- analysis_data_tw$cov_pearsons_sense_red*-1
analysis_data_ab$cov_pearsons_abundance_red <- analysis_data_ab$cov_pearsons_abundance_red*-1


#################### Creating rand vs trati weighted graphs #############################
# run a full model accounting for Species Richness, and FRich (trait-weighted version)
resistance_model_tw_multi <- lm(Resistance_spline ~ Red + cov_pearsons_sense_red, data = analysis_data_tw)
resistance_model_ab_multi <- lm(Resistance_spline ~ Red + cov_pearsons_abundance_red, data = analysis_data_ab)

summary(resistance_model_tw_multi)
summary(resistance_model_ab_multi)

plot_model <- function(model, data, title, grob = c("plot", "legend", "red_hist")) {
  mod <- summary(model)
  
  # Extract the effect of Redundancy
  res_red <- effects::effect(term= "Red", mod= model)
  # Convert to a data frame
  res_red_df <- as.data.frame(res_red)
  
  
  analysis_data_avs <- data %>% group_by(SS) %>% summarise(mean_red = median(Red),
                                                           mean_FRich = median(FRich),
                                                           mean_res = median(Resistance_spline),
                                                           mean_SR = mean(Species_richness_scaled),
                                                           sites = n(),
                                                           Climatic_zone = Climatic_zone[1]) 
  
  #analysis_data_avs$Climatic_zone[analysis_data_avs$Climatic_zone == "Subtropics"] <- "Temperate"
  
  
  ifelse(title == "b", 
         lp <- "bottom",
         lp <- "none")
  
  ifelse(title == "b",
         red_title <- "Redundancy",
         red_title <- "")
  
  if(title == "") {
    res_title <- "Stochastic resistance"
  }
  if(title == "a") {
    res_title <- "Trait-based resistance"
  }
  if(title == "b") {
    res_title <- "Rarity-based resistance"
  }
  
  
    
  var <- round(mod$r.squared,3)
  
  ## new plot plots both trait weighted and random
  res_red_plot <- ggplot(data=res_red_df) +
    
    geom_point(inherit.aes = FALSE, data = analysis_data_avs, aes(x = mean_red, y = mean_res, size = sites), color = "mediumorchid4", alpha = 0.3) +
    
    # line of the effect of Redundancy grouped by simulation type
    geom_line(aes(x=Red, y=fit), col = "black", lwd = 1) +
    
    #geom_abline(slope=1, intercept=0, col = "red", lwd = 2, linetype = "dashed") +
    
    # Add SE ribbons
    geom_ribbon(aes(x=Red, ymin=lower, ymax=upper), alpha= 0.3, lwd = 0, fill = "black") +
    
    # Add real data (but not including RAS - okay in the stats but not the visual)
    
    scale_size_continuous(name = "Sites: ", range = c(5, 20), breaks = c(10,20,40)) +
    
    guides(color = guide_legend(override.aes = list(size = 5))) +
    
    ggtitle(title) +
    
    # Label the axis
    labs(x=red_title, y=res_title) +
    ylim(-3,5) +
    
    guides(color = guide_legend(order = 1, override.aes = list(size = 10)),
           size = guide_legend(order = 2)) +
    
    #guides(size = FALSE
    theme_classic() +
    
    theme(plot.title = element_text(size = 32, face = "bold", vjust = -4, hjust = 0.05), 
          legend.position = "none",
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.direction = "horizontal",
          legend.box = "vertical",
          legend.key.height = unit(1.5, "cm"),
          legend.spacing.x = unit(0.1, 'cm'),
          legend.spacing.y = unit(-0.5, 'cm'),
          legend.background = element_rect(fill = "transparent", colour = NA),
          axis.title.y = element_text(size = 25, vjust = 2),
          axis.title.x = element_text(size = 25),
          axis.text = element_text(size = 20),
          axis.line = element_line(size = 1.5),
          plot.margin = margin(-0.1,0.1,-0.1, 0.1, unit = "npc")) +
          
    annotate("text", x=3.28, y= -2, label =substitute(paste("R"^2, " = ", b), list(b=var)), size = 7) +
    annotate("text", x=2.85, y= -2.75, label=substitute( paste("P-value < 0.001")), size = 7)
  
  ifelse(title == "b",
         res_red_plot <- res_red_plot,
         res_red_plot <- res_red_plot + scale_x_continuous(labels=c(" "), breaks = c(0)))

  label <- paste("R^2 == ", round(mod$r.squared,3))
     
  
  p2 <- res_red_plot
  
  
  
  if(grob == "plot") {
    return(p2)
  }
  
  
  if(grob == "legend") {
     return(get_legend(res_red_plot + theme(legend.position="bottom")))
  }
  
  if(grob == "red_hist") {
    return(xbox)
  }
   
      
}


inset_plot <- function(model) {

############# INSET PLOT ############################
model_df <- as.data.frame(summary(model)$coefficient)[-1,]
model_df$Coefficient <- c("Red",
                          "CovSense")

model_df$Coefficient <- factor(model_df$Coefficient, levels = c("CovSense",
                                                                "Red"))

colnames(model_df)[2] <- "se"

model_df$Estimate[2] <- (model_df$Estimate[2])


#group.colors <- c("Abundance weighted" = "#333BFF", "Trait weighted" = "#CC6600", "Stochastic" ="#9633FF", D = "#E2FF33", E = "#E3DB71")


inset <- ggplot(model_df) +
  geom_bar(aes(x=Coefficient, y=Estimate, fill = Coefficient), stat="identity", alpha=0.7, width = 0.8, position=position_dodge(.9)) +
  geom_errorbar(aes(x=Coefficient, ymin=Estimate - (1.96*se), ymax = Estimate + (1.96*se), width = 0.5), size = 1, position=position_dodge(.8)) +
  geom_hline(yintercept = 0, col = "black", linetype = "dashed") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5), limits = c(-0.4, 1), position = "right") +#,  breaks = c(0, 0.5, 1)) +
  scale_x_discrete(labels = c("FV", "FR")) +
  scale_fill_manual(values = c("Red" = "grey60", "CovSense" = "red")) + 
  #ggtitle("d")  +
  #ylim(-0.05, 1) +
  ylab("Coefficient estimate") +
  coord_flip() +
  theme_classic() + 
  theme(#axis.text.y = element_text(size = 16),
        axis.text = element_text(colour = "black", size = 16),#, angle = 60, hjust = 1.1),
        axis.title.y = element_blank(),#size = 15, vjust = 2.2) ,
        axis.title.x = element_text(size = 20),
        legend.position = "none",
        legend.title = element_blank(),#text(size = 20),
        plot.title = element_text(size = 28, face = "bold", vjust = -8, hjust = 0.05),
        legend.text = element_text(size = 20),
        #plot.margin = margin(3,8,1,0, unit = "cm"),
        rect = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.direction = "horizontal")

return(inset)

}

inset_plot_tw <- inset_plot(resistance_model_tw_multi)
inset_plot_ab <- inset_plot(resistance_model_ab_multi)

tw_plot <- plot_model(resistance_model_tw_multi, analysis_data_tw, "a", grob = "plot")
ab_plot <- plot_model(resistance_model_ab_multi, analysis_data_ab, "b", grob = "plot")
legend <- plot_model(resistance_model_tw_multi, analysis_data_tw, "a", grob = "legend")


tw_plot_wInset <- ggdraw(tw_plot) +
  draw_plot(inset_plot_tw, .2, .7, .45, .3)

ab_plot_wInset <- ggdraw(ab_plot) +
  draw_plot(inset_plot_ab, .2, .7, .4, .3)


png("figure/ED4_resistanance_driven_by_redundancy.png", width = 781/1.35, height = 1000/1.25)
p <- plot_grid(tw_plot_wInset, ab_plot_wInset, NA, legend, NA, ncol = 1, rel_heights = c(1,1,.18,.2, .05), rel_widths = c(.2,1,1,1,1,1))
p
dev.off()  

