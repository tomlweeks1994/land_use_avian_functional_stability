##Plot of bird distribution with variable size of point describing number of birds from that site
##load correct libraries
library(ggplot2)
library(dplyr)
library(maps)
library(tidyr)
library(gridExtra)
library(gtable)
library(grid)
library(ggpubr)
library(cowplot)

#############################################################################
                            ## Panel a ##
#############################################################################
# Read in the Random simulations dat  a
analysis_data<- read.csv("data/analysis_data_trait_0123.csv") %>% 
  distinct(SSBS, .keep_all = TRUE) 


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
         !LandUse == "Secondary vegetation") 


#diversity <- read.csv("Analysis_data/analysis_data_TPD_rand_081221.csv")

corrected <- analysis_data %>%
  distinct(SSBS, .keep_all = TRUE) %>%
  group_by(SS) %>%
  mutate(Sites = n(),
         mean_lat = mean(Latitude),
         mean_long = mean(Longitude)) %>%
  ungroup() %>%
  distinct(SS, .keep_all = TRUE)


#corrected$coords <- paste(corrected$Longitude,corrected$Latitude)

#Pasture land


##download map data for the world

world <- map_data("world")

# virids package for the color palette
#library(viridis)
library(raster)
# create plot with all PREDICTS coordinates with size and colour depicting the amount of info from each site

biomes <- raster("data/biomes_raster.tif")
r <- raster(nrow = 3600, ncol = 8640, ymn = -60)
biomes_lite <- resample(biomes, r, method = "ngb")
biomes_df <- as.data.frame(biomes_lite, xy = TRUE) %>% drop_na()
biomes_df$biomes_raster[biomes_df$biomes_raster > 20] <- NA
biomes_df$biomes_raster <- as.character(biomes_df$biomes_raster)

biomes_df$biomes_raster[biomes_df$biomes_raster == "1"] <- "Tropical & Subtropical Moist Broadleaf Forests"
biomes_df$biomes_raster[biomes_df$biomes_raster == "2"] <- "Tropical & Subtropical Dry Broadleaf Forests"
biomes_df$biomes_raster[biomes_df$biomes_raster == "3"] <- "Tropical & Subtropical Coniferous Forests"
biomes_df$biomes_raster[biomes_df$biomes_raster == "4"] <- "Temperate Broadleaf & Mixed Forests"
biomes_df$biomes_raster[biomes_df$biomes_raster == "5"] <- "Temperate Conifer Forests"
biomes_df$biomes_raster[biomes_df$biomes_raster == "6"] <- "Boreal Forests/Taiga"
biomes_df$biomes_raster[biomes_df$biomes_raster == "7"] <- "Tropical & Subtropical Grasslands, Savannas & Shrublands"
biomes_df$biomes_raster[biomes_df$biomes_raster == "8"] <- "Temperate Grasslands, Savannas & Shrublands"
biomes_df$biomes_raster[biomes_df$biomes_raster == "9"] <- "Flooded Grasslands & Savannas"
biomes_df$biomes_raster[biomes_df$biomes_raster == "10"] <- "Montane Grasslands & Shrublands"
biomes_df$biomes_raster[biomes_df$biomes_raster == "11"] <- "Tundra"
biomes_df$biomes_raster[biomes_df$biomes_raster == "12"] <- "Mediterranean Forests, Woodlands & Scrub"
biomes_df$biomes_raster[biomes_df$biomes_raster == "13"] <- "Deserts & Xeric Shrublands"
biomes_df$biomes_raster[biomes_df$biomes_raster == "14"] <- "Mangroves"

biomes_df$biomes_raster <- factor(biomes_df$biomes_raster, levels = c("Tropical & Subtropical Moist Broadleaf Forests",
                                                                      "Tropical & Subtropical Dry Broadleaf Forests",
                                                                      "Tropical & Subtropical Coniferous Forests",
                                                                      "Temperate Broadleaf & Mixed Forests",
                                                                      "Temperate Conifer Forests",
                                                                      "Boreal Forests/Taiga",
                                                                      "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                                                      "Temperate Grasslands, Savannas & Shrublands",
                                                                      "Flooded Grasslands & Savannas",
                                                                      "Montane Grasslands & Shrublands",
                                                                      "Tundra",
                                                                      "Mediterranean Forests, Woodlands & Scrub",
                                                                      "Deserts & Xeric Shrublands",
                                                                      "Mangroves"))


myColors <- c("forestgreen", "darkolivegreen3", "darkseagreen1", 
              "lightgreen", "lightblue",  "slategray3",
              "khaki", "moccasin", "lightskyblue1", 
              "paleturquoise2", "lightcyan2", "rosybrown1",  
              "cornsilk1", "blue")

names(myColors) <- levels(biomes_df$biomes_raster)
colScale <- scale_fill_manual(name = paste0("Biome"), values = myColors, na.value = "grey80")

main <- ggplot() +
  geom_tile(aes(x = x, y = y, fill = biomes_raster), data = biomes_df) +
  colScale +
  geom_point(data=corrected, aes(x=mean_long, y=mean_lat, size=Sites), color = "black", alpha = 0.6) +
  scale_size_continuous(range=c(2,18), guide = "none") +
  ylim(-60,90) +
  theme_classic() +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("a") +
  guides(fill = guide_legend(nrow = 5, ncol = 3, byrow = TRUE),
         size = "none") +
  theme(axis.line =  element_blank(),
        plot.title = element_text(size = 30, face = "bold", hjust = 0.073, vjust = -12),
        legend.position = "bottom",
        legend.spacing.x = unit(0.5, 'cm'),
        legend.margin=margin(-30,0,0,0),
        legend.box.margin=margin(-30,-10, 0,-10),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),#text(size = 14),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

## create only for shape legend
shape_only <- ggplot() +
  geom_point(data=corrected, aes(x=mean_long, y=mean_lat, size=Sites), color = "black", alpha = 0.8) +
  scale_size_continuous(range=c(1,15)) +
  ylim(-60,90) +
  theme_classic() +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(size = guide_legend(title = "Number of \nAssemblages")) +
  theme(axis.line =  element_blank(),
        #legend.position = "none",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

#shape_legend <- gtable_filter(ggplot_gtable(ggplot_build(shape_only)), "guide-box") 

shape_legend <- get_legend(shape_only)

final_plot_a <- ggdraw(main) +
  draw_plot(shape_legend, x = -0.38, y = 0)

################################################################################
                              ### Panel b ###
################################################################################
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


# Run a standard Richness model
LandUses_rich_study <-  lmer(FRich ~ LandUse + (1|SS) + (1|SSB), data = analysis_data)

data_all_rich <- as.data.frame(summary(LandUses_rich_study)$coefficient)
data_all_rich$LandUse <- c("Pristine primary", "Disturbed primary", "Mature secondary", "Intermediate secondary", "Young secondary",
                           "Plantation", "Pasture","Cropland", "Urban - Minimal", "Urban - Intensive")

data_all_rich$LandUse <- factor(data_all_rich$LandUse, levels = c("Pristine primary", "Disturbed primary", "Mature secondary", "Intermediate secondary", "Young secondary",
                                                                  "Plantation", "Pasture","Cropland", "Urban - Minimal", "Urban - Intensive"))

colnames(data_all_rich)[2] <- "se"
data_all_rich$Estimate[1] <- 0

LandUses_red_study <-  lmer(Red ~ LandUse + (1|SS) + (1|SSB), data = analysis_data)


data_all_red <- as.data.frame(summary(LandUses_red_study)$coefficient)
data_all_red$LandUse <- c("Pristine primary", "Disturbed primary", "Mature secondary", "Intermediate secondary", "Young secondary",
                          "Plantation", "Pasture","Cropland", "Urban - Minimal", "Urban - Intensive")

data_all_red$LandUse <- factor(data_all_red$LandUse, levels = c("Pristine primary", "Disturbed primary", "Mature secondary", "Intermediate secondary", "Young secondary",
                                                                "Plantation", "Pasture","Cropland", "Urban - Minimal", "Urban - Intensive"))

colnames(data_all_red)[2] <- "se"
data_all_red$Estimate[1] <- 0

####################################################

#analysis_data <- analysis_data %>% filter(Species_richness > 4)

analysis_data$Resistance_spline <- as.numeric(scale(analysis_data$Resistance_spline))
analysis_data$Resistance_step <- as.numeric(scale(analysis_data$Resistance_step))
analysis_data$cov_pearsons_sense_redraw <- as.numeric(scale(analysis_data$cov_pearsons_sense_redraw))

LandUses_cov_study <- lmer(cov_pearsons_sense_redraw*-1 ~  LandUse  + (1|SS) + (1|SSB), data = analysis_data)
anova(LandUses_cov_study)
data_all_cov <- as.data.frame(summary(LandUses_cov_study)$coefficient)
data_all_cov$LandUse <- c("Pristine primary", "Disturbed primary", "Mature secondary", "Intermediate secondary", "Young secondary",
                          "Plantation", "Pasture","Cropland", "Urban - Minimal", "Urban - Intensive")

data_all_cov$LandUse <- factor(data_all_cov$LandUse, levels = c("Pristine primary", "Disturbed primary", "Mature secondary", "Intermediate secondary", "Young secondary",
                                                                "Plantation", "Pasture","Cropland", "Urban - Minimal", "Urban - Intensive"))

colnames(data_all_cov)[2] <- "se"
data_all_cov$Estimate[1] <- 0

data_all_rich$Metric <- "Richness"
data_all_red$Metric <- "Redundancy"
data_all_cov$Metric <- "Unique species vulnerability"

data_all <- rbind(data_all_rich, data_all_red, data_all_cov)
data_all$Metric <- factor(data_all$Metric, levels = c("Richness", "Redundancy", "Unique species vulnerability"))

dat_text <- data.frame(
  label = c("b - Functional diversity", "c - Functional redundancy", "d - Functional vulnerability"),
  Metric  = c("Richness", "Redundancy", "Unique species vulnerability")
)

dat_text$Metric <- factor(dat_text$Metric, levels = c("Richness", "Redundancy", "Unique species vulnerability"))


panel_b <- ggplot(data_all) +
  geom_point(aes(y = Estimate, x = LandUse, col = LandUse), position = position_dodge(0.5), size = 6) +
  geom_errorbar(aes(x = LandUse, ymin = Estimate - (1.96*se), ymax = Estimate + (1.96*se), col = LandUse), width=.3, size = 1, position = position_dodge(0.5)) +
  geom_hline(yintercept = 0, col = "red", cex = 0.5, linetype = "dashed") +
  ylim(-0.9, 0.4) +
  ylab("Coefficient estimate") +
  xlab("LandUse") +
  theme_bw() +
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
  geom_text(inherit.aes = FALSE,
            data    = dat_text,
            mapping = aes(x = -Inf, y = -Inf, label = label),
            hjust   = -0.05,
            vjust   = -16.5,
            size = 9,
            col = "black",
            fontface = "bold") +
  facet_wrap(~Metric) +
  theme_classic() +
  #scale_x_discrete(labels = c("-1.0" = " ", "-0.5" = " ", "0.0" = " ", "0.5" = " ")) +
  theme(axis.text.y = element_text(size = 16.5),
        axis.text.x = element_blank(),#angle = 60, size = 17, hjust = 1.1),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        #axis.line.x = element_blank(),
        #axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 25) ,
        legend.position = "bottom",
        legend.key.width = unit(0.02, "npc"),
        legend.text = element_text(size = 20),
        legend.title = element_blank(), 
        plot.margin = margin(2,1,1,1, "cm"),
        plot.title = element_text(size = 25, face = "bold", hjust = -0.25))

library(patchwork)
# 2 is 1200
# 3 is 1150
# 4 is 1100 
# etc etc
png("figures/Figure 2 - black_dots_0108.png", 
    width = 1400, height = 1250)

ggarrange(final_plot_a, panel_b, heights = c(0.6, 0.4), nrow = 2)

#final_plot_a/panel_b + plot_layout(heights = c(0.7, 0.3))
dev.off()

