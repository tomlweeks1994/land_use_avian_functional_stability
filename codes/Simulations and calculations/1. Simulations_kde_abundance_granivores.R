library(tidyverse)   # FD metricx
library(doParallel) # Parallel loops
library(foreach)    # Parallel loops
library(doSNOW)   
library(gawdis)
library(rgl)
library(TPD)

cbind.fill <- function(...) {
  transpoted <- lapply(list(...),t)
  transpoted_dataframe <- lapply(transpoted, as.data.frame)
  return (data.frame(t(plyr::rbind.fill(transpoted_dataframe))))
}

## Read in the neccessary datas
diversity <- readRDS("data/Granivores/diversity_combined_sites_Aug2022_graniv.rds")
com_matrix_all <- read.csv("data/Granivores/com_matrix_all_combined_studies_granivores_Aug2022.csv")

row.names(com_matrix_all) <- com_matrix_all[,1]
colnames(com_matrix_all) <- chartr(".", " ", colnames(com_matrix_all))
com_matrix_all <- com_matrix_all[,-1]


effect_trait_matrix <- read.csv("data/Granivores/PCA_and_diet_data_Aug2022_granivores.csv")
row.names(effect_trait_matrix) <- effect_trait_matrix$Birdlife_Name
effect_trait_matrix <- effect_trait_matrix %>% 
  dplyr::select(-c("X", "Birdlife_Name"))
effect_trait_matrix[5:13] <- effect_trait_matrix[5:13]/100

response_trait_matrix <- read.csv("data/Granivores/response_trait_matrix_granivores_Aug2022.csv")
row.names(response_trait_matrix) <- response_trait_matrix$Birdlife_Name
response_trait_matrix <- response_trait_matrix %>% 
  dplyr::select(-c("X", "Birdlife_Name"))


response_trait_matrix$dispersal_axis <- (scale(-response_trait_matrix$dispersal_axis))
response_trait_matrix$Range.Size <- (scale(-(log(response_trait_matrix$Range.Size))))
response_trait_matrix$specialism <- scale(response_trait_matrix$specialism)
response_trait_matrix$body_axis <- scale(response_trait_matrix$body_axis)

response_trait_matrix$sensitivity <- rowMeans(response_trait_matrix)
response_trait_matrix$sensitivity_shifted <- response_trait_matrix$sensitivity + abs(min(response_trait_matrix$sensitivity)) + 0.1


#create a new dataset of stuff we need for the analysis
analysis_data <- diversity %>%
  dplyr::select(SS,
                SSB,
                SSBS,
                LandUse,
                Use_intensity,
                Use_intensity_2,
                LandUseIntensity,
                LandUseIntensity_2,
                Species_richness) %>%
  distinct(SSBS, .keep_all = TRUE)


studies <- as.character(unique(analysis_data$SS))

##Average weights are then given to gawdis in order to ensure similarity 
source("codes/custom_functions/gaw_rao.R")
source("codes/custom_functions/TPD_RedRich.R")

study_list_all <- list()

for(i in c(1:length(studies))) {
  
  ## set iterator
  stud <- studies[i] 
  print(i)
  
  ## select all sites within study i
  sites <- analysis_data %>%
    filter(SS == stud) %>%
    distinct(SSBS)
  
  ## reduce the community matrix to only have only appropriate sites and species present in study
  com_matrix_study <- com_matrix_all[which(row.names(com_matrix_all) %in% sites$SSBS),]
  com_matrix_study <- com_matrix_study[,which(colSums(com_matrix_study) != 0)]
  
  ## reduce the effect matrix to do the same
  effect_trait_matrix_ss <- effect_trait_matrix[which(row.names(effect_trait_matrix) %in% colnames(com_matrix_study)),]
  effect_trait_matrix_ss <- effect_trait_matrix_ss[,which(colSums(effect_trait_matrix_ss) != 0)]
  
  # using gawdis calculate the functional disimilarity matrix                        
  set.seed(0)
  dist.functional <- gawdis::gawdis(x = effect_trait_matrix_ss[,1:4])#, w.type = "optimized", groups = c(1,2,3,4), opti.maxiter = 1000)
  
  # perform cailliez correction                       
  dm <- as.matrix(cailliez(dist.functional))
  
  ## turn dm into 3-dimensional coordinates
  coordinates <- cmdscale(dm, 3)
  
  ## scale the coordinates
  coordinates_scaled <- as.data.frame(scale(coordinates))
  
  ## get the kernel densities estimates for intraspecific variation
  kernel_densities <- diag(Hpi.diag(coordinates_scaled))
  
  ## create a df of the coordinates
  coordinates_sds_scaled <- data.frame(sd_1 <- rep(kernel_densities[1], nrow(coordinates_scaled)),
                                       sd_2 <- rep(kernel_densities[2], nrow(coordinates_scaled)),
                                       sd_3 <- rep(kernel_densities[3], nrow(coordinates_scaled)))
  
  ## create the trait probability density for the entire study
  ## this parametises the entire study TPD
  ## The code following this will split the study TPDs into sites 
  ## for each site species will be removed 1 by 1 according to sensitivity scores
  ## this will create a large df 
  ## the df will firstly be all the sites in the study
  ## secondly the df will have multiple rows for each site with one species removed each time until no species remain in the site
  ## We calculate the FD at every row - this is our raw extinction curve for each site
  ## We will repeat 100 times to have 100 FD calculations for each row
  
  TPD_mean <- TPDsMean(species = row.names(coordinates), means = coordinates_scaled, sds = coordinates_sds_scaled)#, n_divisions = 50)
  
  ## set a 100 iteration loop                        
  for (j in 1:100) {
    
    ## print outs to track
    print(j)
    print(paste0("simulation: ", j,"-", stud)) 
    
    ## pull the community marix into the correct study
    com_matrix_study <- com_matrix_all[which(row.names(com_matrix_all) %in% sites$SSBS),]
    
    ## remove all species that arent present in study
    com_matrix_study <- com_matrix_study[,-which(colSums(com_matrix_study) == 0)]
    
    ## set an iterations vector as we use as a shortcut to format final output
    iterations <- vector()
    
    ## set another loop that goes through each site in the study matrix
    ## set another loop that goes through each site in the study matrix
    for (k in 1:(nrow(com_matrix_study))) {
      
      ## print out to track
      print(k)
      
      ## get the col numbers of all the species present in the study
      potential_numbers <- which(com_matrix_study[k,] > 0)
      
      ## add a column ABUNDANCE to the sp_probability 
      ## format seems strange as do not need the whole sp_probability df 
      ## however it is made to looks similar to TRAIT weighted codes)
      sp_probability_SSBS$abundance <- as.numeric(com_matrix_study[k, potential_numbers])
      
      # add randomization as sometimes same value of abundance need to randomize which is removed first
      sp_probability_SSBS <- sp_probability_SSBS[sample(nrow(sp_probability_SSBS)),]
      
      # then order by abundance
      sp_probability_SSBS <- sp_probability_SSBS[order(sp_probability_SSBS$abundance),]
      ordered <- match(row.names(sp_probability_SSBS), colnames(com_matrix_study))
      
      ## fill iterations vector as the number of species in the site
      iterations[k] <- length(potential_numbers)
      
      ## com_matrix_study_k is the com_matrix for the site
      com_matrix_study_k <- com_matrix_study[k,]
      
      ## now we run the extinctions simulations for this site
      for (l in 1:length(ordered)) {
        ## print out to track
        print(l)
        
        ## get the colnumber for the most sensitive species
        number <- ordered[l]
        
        ## set their abundance to 0
        com_matrix_study[k, number] <- 0
        
        ## add line underneath the original to be l+1 species removed
        com_matrix_study_k[1 + l,] <- com_matrix_study[k,]
        
      }
      
      ## com_matrix_study_k is now site k within study j with sequential species removed per line
      
      ## if this is the first site of the study set it to com_matrix_study_j
      ## otherwise bind each site until all sites in the study are completed
      ifelse(k == 1,
             com_matrix_study_j <- com_matrix_study_k,
             com_matrix_study_j <- rbind(com_matrix_study_j, com_matrix_study_k))
    }
    
    ## subset the study level TPD into a TPD for all sites
    TPDs_study <- TPDc(TPDs = TPD_mean, sampUnit = com_matrix_study_j)
    
    ## calculate the richness and redundancy at each step
    REND_results <-  TPDRich(TPDc = TPDs_study)
    Rich_TPD <- REND_results$communities$FRichness
    TPD_output <- redundancy(TPDs_study)
    Red_TPD <- TPD_output$redundancy
    
    ## store as a df
    Rich_TPD_j <- data.frame(site = names(Rich_TPD),
                             Rich_TPD_j = as.numeric(Rich_TPD))
    
    Red_TPD_j <- data.frame(site = names(Red_TPD),
                            Red_TPD_j = as.numeric(Red_TPD))
    ## bind
    results_j <- left_join(Rich_TPD_j,Red_TPD_j)
    
    ## if this is the first iteration then set as resuts_rich_TPD_j
    ## but next we will have to bind all 100 iterations
    ifelse(j == 1,
           results_Rich_TPD_j <- results_j %>% dplyr::select(Rich_TPD_j),
           results_Rich_TPD_j <- cbind(results_Rich_TPD_j, results_j %>% dplyr::select(Rich_TPD_j)))
    
    rownames(results_Rich_TPD_j) <- results_j$site
    
    ifelse(j == 1,
           results_Red_TPD_j <- results_j %>% dplyr::select(Red_TPD_j),
           results_Red_TPD_j <- cbind(results_Red_TPD_j, results_j %>% dplyr::select(Red_TPD_j)))
    
    rownames(results_Red_TPD_j) <- results_j$site
    
  }
  
  ## cumsum will allow us to know the break points in the iterations to split into TPDs
  ## this process feels longwinded but is faster
  iterations_cum <- cumsum(iterations)
  
  ## set lists 
  study_list <- list()
  site_datas <- list()
  
  ## lenth iterations is number of sites in study
  for(m in 1:length(iterations)) {
    print(m)
    ## split the FRich and Redundancy values at the break points of the sites
    if(m == 1) {
      site_Rich_TPD_i <- results_Rich_TPD_j[1:iterations_cum[m],]
      site_Red_TPD_i <- results_Red_TPD_j[1:iterations_cum[m],]
    } else {
      site_Rich_TPD_i <- results_Rich_TPD_j[(iterations_cum[m-1]+1):(iterations_cum[m]),]
      site_Red_TPD_i <- results_Red_TPD_j[(iterations_cum[m-1]+1):(iterations_cum[m]),]
      
    }
    
    ## save them as lists 
    site_datas[["Rich_TPD"]] <- site_Rich_TPD_i
    site_datas[["Red_TPD"]] <- site_Red_TPD_i
    
    ## set names of the study_list and put the site data inside 
    ifelse(m == 1,
           study_list[[row.names(results_Rich_TPD_j)[1]]] <- site_datas,
           study_list[[row.names(results_Rich_TPD_j)[iterations_cum[m-1]+1]]] <- site_datas
    )
    
  }
  
  ## sequentially add each study value into the full list as we iterate through i/stud
  
  study_list_all[[stud]] <- study_list
  
}


## save this as the results list
saveRDS(study_list_all, "data/Granivores/rresults_list_TPDs_abundance_weighted_granivores.rds")
