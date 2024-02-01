library(tidyverse)   # FD metricx
library(doParallel) # Parallel loops
library(foreach)    # Parallel loops
library(doSNOW)
library(gawdis)
library(rgl)
library(TPD)
library(ks)


cbind.fill <- function(...) {
  transpoted <- lapply(list(...),t)
  transpoted_dataframe <- lapply(transpoted, as.data.frame)
  return (data.frame(t(plyr::rbind.fill(transpoted_dataframe))))
}

## read in the response trait matrix
response_trait_matrix <- read.csv("data/Main/response_trait_matrix_Aug2022.csv")


## Read in the neccessary datas
diversity <- readRDS("data/Main/diversity_combined_sites_Aug2022.rds")
com_matrix_all <- read.csv("data/Main/com_matrix_all_combined_studies_Aug2022.csv")

## set row names to com_matrix and remove the SSBS column
row.names(com_matrix_all) <- com_matrix_all[,1]
colnames(com_matrix_all) <- chartr(".", " ", colnames(com_matrix_all))
com_matrix_all <- com_matrix_all[,-1]

## read in the effect trait values as we need to recalculate FDs
trait_matrix <- read.csv("data/Main/PCA_and_diet_data_Aug2022.csv")
row.names(trait_matrix) <- trait_matrix$Birdlife_Name

trait_matrix <- trait_matrix %>%
  dplyr::select(-c("X", "Birdlife_Name"))

## add to 1 instead of 100
trait_matrix[5:13] <- trait_matrix[5:13]/100

#create a new dataset of stuff we need for the analysis
analysis_data <- diversity %>%
  dplyr::select(SS,
                SSB,
                SSBS,
                LandUse,
                LandUseIntensity,
                Latitude,
                Longitude) %>%
  distinct(SSBS, .keep_all = TRUE)

source("codes/custom_functions/gaw_rao.R")
source("codes/custom_functions/TPD_RedRich.R")

studies <- unique(analysis_data$SS)

 results_list <- list()                # for loop only       
                            
  for(i in 1:length(studies)){ # for loop only
    stud <- studies[i]         # for loop only
    print(i)
    sites <- analysis_data %>%
      filter(SS == stud) %>%
      distinct(SSBS)
    
    com_matrix_study <- com_matrix_all[which(row.names(com_matrix_all) %in% sites$SSBS),]
    com_matrix_study <- com_matrix_study[,which(colSums(com_matrix_study) != 0)]
    trait_matrix_ss <- trait_matrix[which(row.names(trait_matrix) %in% colnames(com_matrix_study)),]
    trait_matrix_ss <- trait_matrix_ss[,which(colSums(trait_matrix_ss) != 0)]
                             
    ##create the distance matrix
    set.seed(0)
    dist.functional <- gawdis::gawdis(x = trait_matrix_ss, w.type = "optimized", groups = c(1,2,3,4, rep(5, ncol(trait_matrix_ss) - 4)), fuzzy = c(5), opti.maxiter = 1000)
                          
    dm <- as.matrix(cailliez(dist.functional))
    
    ## reconvert into 3D coordinates
    coordinates <- cmdscale(dm, 3)
    coordinates_scaled <- as.data.frame(scale(coordinates))
    
    ## get the estimated intraspecific variation
    kernel_densities <- sqrt(diag(Hpi.diag(coordinates_scaled)))
    coordinates_sds_scaled <- data.frame(sd_1 <- rep(kernel_densities[1], nrow(coordinates_scaled)),
                                         sd_2 <- rep(kernel_densities[2], nrow(coordinates_scaled)),
                                         sd_3 <- rep(kernel_densities[3], nrow(coordinates_scaled)))
    
    ## create study level TPD                         
    TPD_mean <- TPDsMean(species = row.names(coordinates), means = coordinates_scaled, sds = coordinates_sds_scaled)#, n_divisions = 50)
    
    ## create an empty list
    results_i <- list()
    
    ## for each site in the study
    for (j in 1:(nrow(com_matrix_study))) try({
      
    
      ## subset the site
      com_matrix_site <- com_matrix_study[j,]
      
      ## pull the colnumbers where there is a species present
      potential_numbers <- which(com_matrix_site > 0)
        
      ## remove one species each time and recreate a df with 
      for(k in potential_numbers) {
        
        ## take the full
        com_matrix_site_live <- com_matrix_site
        
        ## remove a species
        com_matrix_site_live[, k] <- 0
        
        ##bind these all so com_matrix_site_rem will be multi row df 
        ## first row is the full community
        ## preceeding rows are with a single species removed 
        ifelse(k == potential_numbers[1],
               com_matrix_site_rem <- rbind(com_matrix_site, com_matrix_site_live),
               com_matrix_site_rem <- rbind(com_matrix_site_rem, com_matrix_site_live))
        }
      
      ## add the name of the removed speciwes to the rows
      row.names(com_matrix_site_rem) <- c(paste(row.names(com_matrix_study)[j], "-"), paste(row.names(com_matrix_study)[j], "-", colnames(com_matrix_site)[potential_numbers]))
    
      # subset to create TPDs for all iterations of the site
      TPDs_study <- TPDc(TPDs = TPD_mean, sampUnit = com_matrix_site_rem)
      
      # calculate the redundancy for each of the rows
      TPD_output <- redundancy(TPDs_study)
      site_redundancies <- TPD_output$redundancy
      
      ## we can store the change in redundancy and appropriate that much redundancy to the species 
      SSBS_species_redundancies <- data.frame(species = trimws(gsub(paste(row.names(com_matrix_study)[j], "-"), "", names(site_redundancies))),
                                              redundancy_raw = -(as.numeric(site_redundancies) - as.numeric(site_redundancies[1])),
                                              redundancy_prop = (-( as.numeric(site_redundancies) - as.numeric(site_redundancies[1])) / as.numeric(site_redundancies[j]) * 100))

      
      
      ## then allign the sensitivity scores                         
      response_trait_matrix_SSBS <- response_trait_matrix %>% filter(Birdlife_Name %in% SSBS_species_redundancies$species)
      
      response_trait_matrix_SSBS$specialism <- scale(response_trait_matrix_SSBS$specialism)
      response_trait_matrix_SSBS$body_axis <- scale(response_trait_matrix_SSBS$body_axis)
      response_trait_matrix_SSBS$Range.Size <- scale(-log(response_trait_matrix_SSBS$Range.Size))
      response_trait_matrix_SSBS$dispersal_axis <- scale(-(response_trait_matrix_SSBS$dispersal_axis))
      
      response_trait_matrix_SSBS$sensitivity <- rowMeans(response_trait_matrix_SSBS[,3:6])
      
      SSBS_species_redundancies$study_sensitivity <- response_trait_matrix_SS$sensitivity[match(SSBS_species_redundancies$species, response_trait_matrix_SS$Birdlife_Name)]
      SSBS_species_redundancies$site_sensitivity <- response_trait_matrix_SSBS$sensitivity[match(SSBS_species_redundancies$species, response_trait_matrix_SSBS$Birdlife_Name)]
      SSBS_species_redundancies <- SSBS_species_redundancies[-1,]
      
      ## align the abundances as well                           
      SSBS_species_redundancies$abundance <- as.numeric(com_matrix_site[which(com_matrix_site > 0)])
      
      ## metrics are then the covariances between species abundance and species redundancy or sensitivity and redundancy
      metric_j <- data.frame(SSBS = row.names(com_matrix_site),
                             cov_pearsons_sense_redraw = cov(SSBS_species_redundancies$site_sensitivity, SSBS_species_redundancies$redundancy_raw),
                             cov_pearsons_abundance_redraw = cov(SSBS_species_redundancies$abundance, SSBS_species_redundancies$redundancy_raw))
      
      
      ## bind sequentially for each study                     
      
        ifelse(j ==1,
               metric_i <-  metric_j,
               metric_i <- rbind(metric_i, metric_j))
        
        })
    ## add studies sequentially too
    ifelse(i ==1,
           metrics_all <-  metric_i,
           metrics_all <- rbind(metrics_all, metric_i))
    
  }

## write file with FV results   
write.csv(metrics_all, "data/covariance_metrics_fullsample.csv")

metrics_abundance <- metrics_all %>% dplyr::select(!cov_pearsons_sense_redraw)
metrics_trait <- metrics_all %>% dplyr::select(!cov_pearsons_abundance_redraw)

## combine with red rich and stability results
red_rich_stab_t <- read.csv("data/fullsample_all_results_trait.csv")
red_rich_stab_vuln_t <- left_join(red_rich_stab_t, metrics_trait, by = "SSBS")

## link to the analysis data file to use for statistical analysis
analysis_data_fullset_t <- left_join(analysis_data, red_rich_stab_vuln_t, by = "SSBS")

write.csv(analysis_data_fullset_t, "data/analysis_data_trait_full.csv")


### ALSO SAVE THESE INTO THE ABUNDANCE DATASET
## combine with red rich and stability results
red_rich_stab_a <- read.csv("data/fullsample_all_results_abundance.csv")
red_rich_stab_vuln_a <- left_join(red_rich_stab_a, metrics_abundance, by = "SSBS")

## link to the analysis data file to use for statistical analysis
analysis_data_fullset_a <- left_join(analysis_data, red_rich_stab_vuln_a, by = "SSBS")

write.csv(analysis_data_fullset_a, "data/analysis_data_abundance_full.csv")