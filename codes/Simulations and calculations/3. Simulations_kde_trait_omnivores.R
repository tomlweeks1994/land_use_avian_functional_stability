library(tidyverse)   # FD metricx
library(doParallel) # Parallel loops
library(foreach)    # Parallel loops
library(doSNOW)   
library(gawdis)
library(rgl)
library(TPD)

setwd("/home/tomlweeks/PREDICTS")

cbind.fill <- function(...) {
  transpoted <- lapply(list(...),t)
  transpoted_dataframe <- lapply(transpoted, as.data.frame)
  return (data.frame(t(plyr::rbind.fill(transpoted_dataframe))))
}

## Read in the neccessary datas
diversity <- readRDS("Datasets/Aug2022/Omnivores/diversity_combined_sites_Aug2022_omniv.rds")
com_matrix_all <- read.csv("Datasets/Aug2022/Omnivores/com_matrix_all_combined_studies_omnivores_Aug2022.csv")

row.names(com_matrix_all) <- com_matrix_all[,1]
colnames(com_matrix_all) <- chartr(".", " ", colnames(com_matrix_all))
com_matrix_all <- com_matrix_all[,-1]


effect_trait_matrix <- read.csv("Datasets/Aug2022/Omnivores/PCA_and_diet_data_Aug2022_omnivores.csv")
row.names(effect_trait_matrix) <- effect_trait_matrix$Birdlife_Name
effect_trait_matrix <- effect_trait_matrix %>% 
  dplyr::select(-c("X", "Birdlife_Name"))
effect_trait_matrix[5:13] <- effect_trait_matrix[5:13]/100

response_trait_matrix <- read.csv("Datasets/Aug2022/Omnivores/response_trait_matrix_omnivores_Aug2022.csv")
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

# analysis_data <- analysis_data %>%
#   filter(SS == "2020_Lees_RAS")

no_birds <- analysis_data %>%
  filter(Species_richness == 0)

analysis_data <- analysis_data %>%
  filter(Species_richness != 0)

# write.csv(analysis_data, "Analysis_data/template_for_study_metrics.csv")
# write.csv(no_birds, "Analysis_data/template_for_study_metrics_nb.csv")

studies <- as.character(unique(analysis_data$SS))

##Average weights are then given to gawdis in order to ensure similarity 
source("Codes/Functions/gaw_rao.R")
source("Codes/Functions/TPD_RedRich.R")

# for(i in c(1:length(studies))) {
#   print(paste0("study: ", i))
#   stud <- studies[i]

#studies_temp <- studies[1:6]

#studies <- studies[c(1,5)]

# print("Starting parallel")
# registerDoParallel(cores = 50)                                                     # foreach only
# results_list <- foreach(stud = studies,                                           # foreach only
#                         .packages = c("dplyr", "ade4","SYNCSA", "gawdis", "TPD")) %dopar% {   # foreach only
study_list_all <- list()
for(i in c(1:length(studies))) try({
  
  stud <- studies[i] 
  print(i)
                          #results_list <- list()                # for loop only       
                          #for(i in 1:length(studies)){ # for loop only
                          #stud <- studies[i]         # for loop only
                          
                          sites <- analysis_data %>%
                            filter(SS == stud) %>%
                            distinct(SSBS)
                          
                          com_matrix_study <- com_matrix_all[which(row.names(com_matrix_all) %in% sites$SSBS),]
                          com_matrix_study <- com_matrix_study[,which(colSums(com_matrix_study) != 0)]
                          effect_trait_matrix_ss <- effect_trait_matrix[which(row.names(effect_trait_matrix) %in% colnames(com_matrix_study)),]
                          effect_trait_matrix_ss <- effect_trait_matrix_ss[,which(colSums(effect_trait_matrix_ss) != 0)]
                          ##Get the weights
                          set.seed(0)
                          dist.functional <- gawdis::gawdis(x = effect_trait_matrix_ss[,1:4])#, w.type = "optimized", groups = c(1,2,3,4), opti.maxiter = 1000)
                          
                          dm <- as.matrix(cailliez(dist.functional))
                          
                          coordinates <- cmdscale(dm, 3)
                          coordinates_scaled <- as.data.frame(scale(coordinates))
                          kernel_densities <- diag(Hpi.diag(coordinates_scaled))
                          
                          
                          coordinates_sds_scaled <- data.frame(sd_1 <- rep(kernel_densities[1], nrow(coordinates_scaled)),
                                                               sd_2 <- rep(kernel_densities[2], nrow(coordinates_scaled)),
                                                               sd_3 <- rep(kernel_densities[3], nrow(coordinates_scaled)))
                          
                          # coordinates_sds <- data.frame(sd_1 = rep(sd(coordinates[,1])/2, nrow(coordinates)), 
                          #                              sd_2 = rep(sd(coordinates[,2])/2, nrow(coordinates)),
                          #                              sd_3 = rep(sd(coordinates[,2])/2, nrow(coordinates)))
                          
                          TPD_mean <- TPDsMean(species = row.names(coordinates), means = coordinates_scaled, sds = coordinates_sds_scaled)#, n_divisions = 50)
                          #TPD_diss <- TPD::dissim(TPD_mean)
                          
                          for (j in 1:1) {
                            print(j)
                            print(paste0("simulation: ", j,"-", stud)) # for loop only
                            
                            com_matrix_study <- com_matrix_all[which(row.names(com_matrix_all) %in% sites$SSBS),]
                            
                            com_matrix_study <- com_matrix_study[,-which(colSums(com_matrix_study) == 0)]
                            
                            iterations <- vector()
                            
                            for (k in 1:(nrow(com_matrix_study))) {
                              print(k)
                              #print(paste("k", k))
                              
                              potential_numbers <- which(com_matrix_study[k,] > 0)
                              sp_probability_SSBS <- response_trait_matrix  %>% slice(which(row.names(response_trait_matrix) %in% colnames(com_matrix_study)[potential_numbers]))
                              sp_probability_SSBS$likelihood <- max(sp_probability_SSBS$sensitivity_shifted)/sp_probability_SSBS$sensitivity_shifted
                              sp_probability_SSBS$probability <- sp_probability_SSBS$likelihood/sum(sp_probability_SSBS$likelihood)
                              
                              ordered <- potential_numbers[sample(1:nrow(sp_probability_SSBS), size=length(potential_numbers), replace=FALSE, prob=c(sp_probability_SSBS$probability))]
                              
                              iterations[k] <- length(potential_numbers)
                              
                              com_matrix_study_k <- com_matrix_study[k,]
                              
                              for (l in 1:length(ordered)) {
                                print(l)
                                number <- ordered[l]
                                #print(paste("l", l))
                                
                                #potential_numbers[sample(1:nrow(sp_probability_SSBS), size=1, replace=FALSE, prob=c(sp_probability_SSBS$probability))]
                                
                                #random_num <- potential_random_numbers[round(runif(1, 0.5, length(potential_random_numbers) + 0.499999999 ))]
                                
                                
                                com_matrix_study[k, number] <- 0
                                
                                com_matrix_study_k[1 + l,] <- com_matrix_study[k,]
                                
                              }
                              
                              ifelse(k == 1,
                                     com_matrix_study_j <- com_matrix_study_k,
                                     com_matrix_study_j <- rbind(com_matrix_study_j, com_matrix_study_k))
                              
                            }
                            
                            TPDs_study <- TPDc(TPDs = TPD_mean, sampUnit = com_matrix_study_j)
                            
                            REND_results <-  TPDRich(TPDc = TPDs_study)
                            Rich_TPD <- REND_results$communities$FRichness
                            #Q_TPD <- Rao(diss = TPD_diss, TPDc = TPDs_study)
                            TPD_output <- redundancy(TPDs_study)
                            
                            #Rao_TPD <- Q_TPD$alpha_rao
                            Red_TPD <- TPD_output$redundancy
                            RedRel_TPD <- TPD_output$redundancyRelative
                            
                            
                            
                            #Rao_TPD_j <- data.frame(site = names(Rao_TPD),
                            #                        Rao_TPD_j = as.numeric(Rao_TPD))
                            
                            Rich_TPD_j <- data.frame(site = names(Rich_TPD),
                                                     Rich_TPD_j = as.numeric(Rich_TPD))
                            
                            Red_TPD_j <- data.frame(site = names(Red_TPD),
                                                    Red_TPD_j = as.numeric(Red_TPD))
                            
                            RedRel_TPD_j <- data.frame(site = names(RedRel_TPD),
                                                       RedRel_TPD_j = as.numeric(RedRel_TPD))
                            
                            
                            #results_j <- left_join(Rao_TPD_j, Rich_TPD_j)
                            #results_j <- left_join(results_j, Red_TPD_j)
                            #results_j <- left_join(results_j, RedRel_TPD_j)
                            results_j <- left_join(Rich_TPD_j,Red_TPD_j)
                            results_j <- left_join(results_j, RedRel_TPD_j)
                            
                            
                            
                            #ifelse(j == 1,
                            #       results_Rao_TPD_j <- results_j %>% dplyr::select(Rao_TPD_j),
                            #       results_Rao_TPD_j <- cbind(results_Rao_TPD_j, results_j %>% dplyr::select(Rao_TPD_j)))
                            #rownames(results_Rao_TPD_j) <- results_j$site
                            
                            ifelse(j == 1,
                                   results_Rich_TPD_j <- results_j %>% dplyr::select(Rich_TPD_j),
                                   results_Rich_TPD_j <- cbind(results_Rich_TPD_j, results_j %>% dplyr::select(Rich_TPD_j)))
                            rownames(results_Rich_TPD_j) <- results_j$site
                            
                            ifelse(j == 1,
                                   results_Red_TPD_j <- results_j %>% dplyr::select(Red_TPD_j),
                                   results_Red_TPD_j <- cbind(results_Red_TPD_j, results_j %>% dplyr::select(Red_TPD_j)))
                            rownames(results_Red_TPD_j) <- results_j$site
                            
                            ifelse(j == 1,
                                   results_RedRel_TPD_j <- results_j %>% dplyr::select(RedRel_TPD_j),
                                   results_RedRel_TPD_j <- cbind(results_RedRel_TPD_j, results_j %>% dplyr::select(RedRel_TPD_j)))
                            rownames(results_RedRel_TPD_j) <- results_j$site
                            
                          }
                          
                          iterations_cum <- cumsum(iterations)
                          
                          study_list <- list()
                          site_datas <- list()
                          
                          for(m in 1:length(iterations)) {
                            print(m)
                            if(m == 1) {
                              #site_Rao_TPD_i <- results_Rao_TPD_j[1:iterations_cum[m],]
                              site_Rich_TPD_i <- results_Rich_TPD_j[1:iterations_cum[m],]
                              site_Red_TPD_i <- results_Red_TPD_j[1:iterations_cum[m],]
                              site_RedRel_TPD_i <- results_RedRel_TPD_j[1:iterations_cum[m],]
                            } else {
                              #site_Rao_TPD_i <- results_Rao_TPD_j[(iterations_cum[m-1]+1):(iterations_cum[m]),]
                              site_Rich_TPD_i <- results_Rich_TPD_j[(iterations_cum[m-1]+1):(iterations_cum[m]),]
                              site_Red_TPD_i <- results_Red_TPD_j[(iterations_cum[m-1]+1):(iterations_cum[m]),]
                              site_RedRel_TPD_i <- results_RedRel_TPD_j[(iterations_cum[m-1]+1):(iterations_cum[m]),]
                              
                            }
                            
                            #site_datas[["Rao_TPD"]] <- site_Rao_TPD_i
                            site_datas[["Rich_TPD"]] <- site_Rich_TPD_i
                            site_datas[["Red_TPD"]] <- site_Red_TPD_i
                            site_datas[["RedRel_TPD"]] <- site_RedRel_TPD_i
                            
                            ifelse(m == 1,
                               study_list[[row.names(results_Rich_TPD_j)[1]]] <- site_datas,
                               study_list[[row.names(results_Rich_TPD_j)[iterations_cum[m-1]+1]]] <- site_datas
                            )
                               
                               
                          }
                          
                          study_list_all[[stud]] <- study_list
                        })



saveRDS(study_list_all, "Datasets/Aug2022/Omnivores/results_list_TPDs_trait_weighted_omnivores.rds")
