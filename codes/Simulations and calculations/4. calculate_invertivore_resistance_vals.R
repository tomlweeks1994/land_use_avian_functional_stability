library(tidyverse)
library(MESS)
library(Bolstad2)
library(broom)
library(minpack.lm)
library(ggpubr)
library(segmented)

#results_list <- readRDS("Analysis_data/results_list_random_simulation_2308.rds")
results_list <- readRDS("Datasets/Aug2022/Invertivores/results_list_TPDs_trait_weighted_invertivores.rds") 

diversity <- readRDS("Datasets/Aug2022/Invertivores/diversity_combined_sites_Aug2022_inv.rds")

analysis_data <- diversity %>% dplyr::distinct(SSBS, .keep_all = TRUE) %>%
  dplyr::select("SS", "SSB", "SSBS", "LandUse", "Use_intensity", "Use_intensity_2", 
                "LandUseIntensity", "LandUseIntensity_2", "Species_richness", 
                "Latitude", "Longitude")

for (i in c(1:length(results_list))) {
  print(paste("i:", i))
  study_list <- results_list[i][[1]]
  j_its <- length(study_list)
  # ifelse(i == 97,
  #        j_its <- 158,
  #        j_its <- length(study_list))
  # 
  
  for (j in c(1:j_its)) {
    #print(j)
    #if(i == 99) {
    #names(study_list[[j]]) <- c("Rich_TPD", "Red_TPD")
    #site_i <- study_list[[j]]
    #SSBS = site_i$Rich_TPD$site[1]
    #site_i$Rich_TPD <- as.data.frame(site_i$Rich_TPD[,-1])
    #site_i$Red_TPD <- as.data.frame(site_i$Red_TPD[,-1])
    #} else {
    
    site_i <- study_list[[j]]
    site_i$Rich_TPD <- as.data.frame(site_i$Rich_TPD)
    site_i$Rao_TPD <- as.data.frame(site_i$Rao_TPD)
    site_i$Red_TPD <- as.data.frame(site_i$Red_TPD)
    site_i$RedRel_TPD <- as.data.frame(site_i$RedRel_TPD)
    
    SSBS <- names(study_list[j])
    #}
    
    Red <- site_i$Red_TPD[[1]][1]
    FRich <- site_i$Rich_TPD[[1]][1]
    
    
    metrics_i <- data.frame(SSBS = SSBS,
                            Red = Red,
                            FRich = FRich,
                            Resistance_spline = NA,
                            Resistance_step = NA,
                            loss_50_FRich = NA
    )
    
    
    #print(paste("i:",i,"j:",j,"- Hyp"))
    
    #site_i$Rich_TPD[nrow(site_i$Rich_TPD) + 1,] <- 0
    
    Resistance_step <- vector()
    Resistance_spline <- vector()
    
    for(k in c(1:length(site_i$Rich_TPD))) {
      if(length(site_i$Rich_TPD[[k]][is.na(site_i$Rich_TPD[[k]])]) > 0){
        print(paste("i:", i, "j :", j,"NAs in iteration:", length(site_i$Rich_TPD[[k]][is.na(site_i$Rich_TPD[[k]])])))
      }
      site_i$Rich_TPD[[k]][is.na(site_i$Rich_TPD[[k]])] <- 0
      
      simulation_i_FRich <- site_i$Rich_TPD[[k]]
      simulation_i_FRich <- simulation_i_FRich/max(simulation_i_FRich) # Scaling it to a proportion of starting FD
      simulation_i_FRich[length(simulation_i_FRich) + 1]  <- 0
      simulation_i_FRich <- data.frame(FD = simulation_i_FRich, iteration = seq(100, 0, length.out = length(simulation_i_FRich))) #Also Scaled
      
      Resistance_step[k] <- bayestestR::area_under_curve(simulation_i_FRich$iteration, y = simulation_i_FRich$FD, method = "step")
      
      Resistance_spline[k] <- sintegral(x = simulation_i_FRich$iteration, fx = simulation_i_FRich$FD)$int
      
    }
    
    metrics_i$Resistance_spline <- mean(Resistance_spline)
    metrics_i$Resistance_step <- mean(Resistance_step)
    
    site_av <- data.frame(FD = c(as.numeric(rowMeans(site_i$Rich_TPD)),0), Red = c(as.numeric(rowMeans(site_i$Red_TPD)),0), iteration = seq(100, 0, length.out = nrow(site_i$Rich_TPD) + 1))
    metrics_i$loss_50_FRich = which.min(abs(site_av$FD - max(site_av$FD)/2)) / length(site_av$FD)
    
    ifelse(i == 1 & j == 1,
           metrics <- metrics_i,
           metrics <- rbind(metrics, metrics_i))
  }
}


metrics <- metrics %>% distinct(SSBS, .keep_all = TRUE)

analysis_data <- left_join(analysis_data, metrics, by = "SSBS")

write.csv(analysis_data, "Datasets/Aug2022/Invertivores/invertivore_resistance_trait.csv")

