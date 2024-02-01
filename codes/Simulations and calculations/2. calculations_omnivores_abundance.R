library(tidyverse)
library(MESS)
library(Bolstad2)
library(broom)
library(minpack.lm)
library(ggpubr)
library(segmented)

results_list <- readRDS("data/Omnivores/simulations_abundance_omnivores.rds") 

diversity <- readRDS("data/Omnivores/diversity_combined_sites_Aug2022_omniv.rds")

analysis_data <- diversity %>% dplyr::distinct(SSBS, .keep_all = TRUE) %>%
  dplyr::select("SS", "SSB", "SSBS", "LandUse", "Use_intensity", "Use_intensity_2", 
                "LandUseIntensity", "LandUseIntensity_2", "Species_richness", 
                "Latitude", "Longitude")

for (i in c(1:length(results_list))) {

  ## print_out
  print(paste("i:", i))
  
  ## take study i
  study_list <- results_list[i][[1]]
  
  ## length study list is the number of sites in the study
  j_its <- length(study_list)

  ## new loop
  for (j in c(1:j_its)) {

    ## take the first site
    site_i <- study_list[[j]]
    
    ## split into seperate values for frich and redundancy
    site_i$Rich_TPD <- as.data.frame(site_i$Rich_TPD)
    site_i$Red_TPD <- as.data.frame(site_i$Red_TPD)
    
    ## save the name for later
    SSBS <- names(study_list[j])
    
    ## first line Red and Rich is the value for the whole community
    ## save this as our Rich and Red values
    Red <- site_i$Red_TPD[[1]][1]
    FRich <- site_i$Rich_TPD[[1]][1]
    
    ## store in a df with a few more lines to fill in
    metrics_i <- data.frame(SSBS = SSBS,
                            Red = Red,
                            FRich = FRich,
                            Resistance_spline = NA,
                            Resistance_step = NA,
                            loss_50_FRich = NA
    )
    
    ## create vectors to store the extinction curve resistance values
    ## we will take the av of these later as our resistance for each site
    Resistance_step <- vector()
    Resistance_spline <- vector()
    
    ## set loop
    for(k in c(1:length(site_i$Rich_TPD))) {
      
      ## if there are NAs set them to 0 (these occur when there are insufficient species)
      site_i$Rich_TPD[[k]][is.na(site_i$Rich_TPD[[k]])] <- 0
      
      ## take column k as a vector
      simulation_i_FRich <- site_i$Rich_TPD[[k]]
      
      # Scale it to a proportion of starting FD
      simulation_i_FRich <- simulation_i_FRich/max(simulation_i_FRich) 
      
      # add a row at the end and set to 0 for No species present
      simulation_i_FRich[length(simulation_i_FRich) + 1]  <- 0
      
      # add an iteration going from 100 to 0 (this is proportion of species remaining)
      simulation_i_FRich <- data.frame(FD = simulation_i_FRich, iteration = seq(100, 0, length.out = length(simulation_i_FRich)))
      
      ## calculate the AUC (two approaches)
      Resistance_step[k] <- bayestestR::area_under_curve(simulation_i_FRich$iteration, y = simulation_i_FRich$FD, method = "step")
      Resistance_spline[k] <- sintegral(x = simulation_i_FRich$iteration, fx = simulation_i_FRich$FD)$int
      
    }
    
    ## take the means of the resistance vectors
    metrics_i$Resistance_spline <- mean(Resistance_spline)
    metrics_i$Resistance_step <- mean(Resistance_step)
    
    ## also take the rowmeans to pull out the average halflife
    site_av <- data.frame(FD = c(as.numeric(rowMeans(site_i$Rich_TPD)),0), Red = c(as.numeric(rowMeans(site_i$Red_TPD)),0), iteration = seq(100, 0, length.out = nrow(site_i$Rich_TPD) + 1))
    metrics_i$loss_50_FRich = which.min(abs(site_av$FD - max(site_av$FD)/2)) / length(site_av$FD)
    
    ## sequentially bind
    ifelse(i == 1 & j == 1,
           metrics <- metrics_i,
           metrics <- rbind(metrics, metrics_i))
  }
}

## join my analysis data to the metrics
analysis_data <- left_join(analysis_data, metrics, by = "SSBS")

write.csv(analysis_data, "data/Omnivores/omniv_all_results_abundance.csv")

