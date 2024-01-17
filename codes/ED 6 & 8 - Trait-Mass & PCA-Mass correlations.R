#====================================================================================================================
# 2 step-PCA and take trophic, loco and size create one for dispersal then do a pcoA to get one for diet choice
#====================================================================================================================
library(tidyverse)
library(sf)
library(rgeos)
library(patchwork)
library(sp)
library(ggpubr)
library(rgdal)
library(lwgeom)

diversity <- readRDS("data/diversity_combined_sites_Aug2022.rds")


trait_matrix_all <- read.csv("data/trait_matrix_all.csv")

trait_matrix_troph <-trait_matrix_all %>% 
  dplyr::select(Species1, 
                Beak.Length_Culmen, 
                Beak.Length_Nares, 
                Beak.Width, 
                Beak.Depth,
                Mass)

trait_matrix_loco <-trait_matrix_all %>% 
  dplyr::select(Species1, 
                Wing.Length,
                Tail.Length,
                Tarsus.Length,
                Mass)

row.names(trait_matrix_troph) <- trait_matrix_troph$Species1
trait_matrix_troph[-1] <- apply(trait_matrix_troph[-1], 2, as.numeric)
trait_matrix_troph[-1] <- apply(trait_matrix_troph[-1], 2, log)
trait_matrix_troph[-1] <- apply(trait_matrix_troph[-1], 2, scale)
trait_matrix_troph <- trait_matrix_troph %>% dplyr::select(-"Species1")

row.names(trait_matrix_loco) <- trait_matrix_loco$Species1
trait_matrix_loco[-1] <- apply(trait_matrix_loco[-1], 2, as.numeric)
trait_matrix_loco[-1] <- apply(trait_matrix_loco[-1], 2, log)
trait_matrix_loco[-1] <- apply(trait_matrix_loco[-1], 2, scale)
trait_matrix_loco <- trait_matrix_loco %>% dplyr::select(-"Species1")

row.names(trait_matrix_all) <- trait_matrix_all$Species1
trait_matrix_all[-1] <- apply(trait_matrix_all[-1], 2, as.numeric)
trait_matrix_all[-1] <- apply(trait_matrix_all[-1], 2, log)
trait_matrix_all[-1] <- apply(trait_matrix_all[-1], 2, scale)
trait_matrix_all <- trait_matrix_all %>% dplyr::select(-"Species1")


long_version <- trait_matrix_all %>% 
  gather(var, value, -Mass)

long_version$var[long_version$var == "Beak.Length_Culmen"] <- "Beak length (Culmen), mm"
long_version$var[long_version$var == "Beak.Length_Nares"] <- "Beak length (Nares), mm"
long_version$var[long_version$var == "Beak.Width"] <- "Beak width, mm"
long_version$var[long_version$var == "Beak.Depth"] <- "Beak depth, mm"
long_version$var[long_version$var == "Wing.Length"] <- "Wing length, mm"
long_version$var[long_version$var == "Tail.Length"] <- "Tail length, mm"
long_version$var[long_version$var == "Tarsus.Length"] <- "Tarsus length, mm"
long_version$var[long_version$var == "Hand.Wing.Index"] <- "Hand-wing index"

long_version$var <- factor(long_version$var, levels = c("Beak length (Culmen), mm",
                                                        "Beak length (Nares), mm",
                                                        "Beak width, mm",
                                                        "Beak depth, mm",
                                                        "Wing length, mm",
                                                        "Tail length, mm",
                                                        "Tarsus length, mm",
                                                        "Hand-wing index"))

dat_text_1 <- data.frame(
  label = c("a - Beak length (Culmen)", "b - Beak length (Nares)", "c - Beak width", "d - Beak depth"),
  var  = c("Beak length (Culmen), mm",
           "Beak length (Nares), mm",
           "Beak width, mm",
           "Beak depth, mm")
)

dat_text_1$var <- factor(dat_text_1$var, levels =c("Beak length (Culmen), mm",
                                               "Beak length (Nares), mm",
                                               "Beak width, mm",
                                               "Beak depth, mm"))


dat_text_2 <- data.frame(
  label = c("e - Wing length", "f - Tail length", "g - Tail length", "h - Hand-wing index"),
  var  = c("Wing length, mm",
           "Tail length, mm",
           "Tarsus length, mm",
           "Hand-wing index")
)

dat_text_2$var <- factor(dat_text_2$var, levels =c("Wing length, mm",
                                               "Tail length, mm",
                                               "Tarsus length, mm",
                                               "Hand-wing index"))
dev.size("px")
#[1] 1605  760

long_version$d <- densCols(x = long_version$Mass, 
                       y = long_version$value, 
                       colramp = colorRampPalette(viridis::viridis(10)))

long_version_1 <- long_version %>% filter(var %in% c("Beak length (Culmen), mm",
                                                     "Beak length (Nares), mm",
                                                     "Beak width, mm",
                                                     "Beak depth, mm"))

long_version_1$var <-  factor(long_version_1$var, levels =c("Beak length (Culmen), mm",
                                                      "Beak length (Nares), mm",
                                                      "Beak width, mm",
                                                      "Beak depth, mm"))


long_version_2 <- long_version %>% filter(var %in% c("Wing length, mm",
                                                     "Tail length, mm",
                                                     "Tarsus length, mm",
                                                     "Hand-wing index"))

long_version_2$var <-  factor(long_version_2$var, levels =c("Wing length, mm",
                                                            "Tail length, mm",
                                                            "Tarsus length, mm",
                                                            "Hand-wing index"))



plot_high <- ggplot(long_version_1) +
  geom_point(aes(Mass, value), col = "grey80")  +
  coord_cartesian(clip = "off") +
  #scale_color_identity(aes(fill = d)) +
  ylab("Trait size") +
  ylim(-3.5, 5.5) +
  xlab("Body mass") +
  geom_smooth(aes(Mass, value), data = long_version_1 %>% filter(var != "Hand-wing index"), method = "lm", se = TRUE, col = "red") +
  theme_classic() +
  stat_cor(aes(Mass, value, label = ..r.label..),
               label.x = 2.9, label.y = -2.75, size = 8) +
  facet_wrap(~var, nrow = 1) +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(), #_text(size = 20),
        legend.position = "none",
        panel.grid = element_blank()) +
    geom_text(inherit.aes = FALSE,
              data    = dat_text_1,
              mapping = aes(x = -Inf, y = -Inf, label = label),
              hjust   = -0.05,
              vjust   = -15,
              size = 8,
              fontface = "bold",
              col = "black") 


plot_low <- ggplot(long_version_2) +
  geom_point(aes(Mass, value), col = "grey80")  +
  coord_cartesian(clip = "off") +
  #scale_color_identity(aes(fill = d)) +
  ylab("Trait size") +
  ylim(-3.5, 5.5) +
  xlab("Body mass") +
  geom_smooth(aes(Mass, value), data = long_version_2 %>% filter(var != "Hand-wing index"), method = "lm", se = TRUE, col = "red") +
  theme_classic() +
  stat_cor(aes(Mass, value, label = ..r.label..),
           label.x = 2.9, label.y = -2.75, size = 8) +
  facet_wrap(~var, nrow = 1) +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(), #(size = 20),
        legend.position = "none",
        panel.grid = element_blank()) +
  geom_text(inherit.aes = FALSE,
            data    = dat_text_2,
            mapping = aes(x = -Inf, y = -Inf, label = label),
            hjust   = -0.075,
            vjust   = -15,
            size = 8,
            fontface = "bold",
            col = "black") 


png("figures/ED 6 - BM x Traits.png",
    width = 1385, height = 895/1.25)
plot(plot_high / plot_low)
dev.off()

## Just get the columns with the morphometric trait data
#### Log transform - then standardise and centre for a mean of zero and a standard deviation of 1


#Perform the PCA
Full.pca <- princomp(trait_matrix_troph, scores = TRUE)

#plot to see how correlated to size it was
#plot(Scaled_BM, Full.pca$scores[,1])
body_PC1 <- Full.pca$scores[,1]
troph_axis <- Full.pca$scores[,2]

#Perform the PCA
Full.pca2 <- princomp(trait_matrix_loco, scores = TRUE)

#plot to see how correlated to size it was
#plot(Scaled_BM, Full.pca$scores[,1])
body_PC2 <- Full.pca2$scores[,1]
loco_axis <- Full.pca2$scores[,2]

beak_bod <- data.frame(PC = body_PC1, body_size = trait_matrix_troph$Mass, axis = "Beak PCA", component = "PC1")
beak_bod2 <- data.frame(PC = troph_axis, body_size = trait_matrix_troph$Mass, axis = "Beak PCA", component = "PC2")
loco_bod <-  data.frame(PC = body_PC2,  body_size = trait_matrix_loco$Mass, axis = "Locomotory PCA", component = "PC1")
loco_bod2 <-  data.frame(PC = loco_axis, body_size = trait_matrix_loco$Mass, axis = "Locomotory PCA", component = "PC2")

PC1s_bod <- rbind(beak_bod, beak_bod2, loco_bod, loco_bod2)

Body_axis <- princomp(cbind(body_PC2, body_PC2), scores = TRUE)

bod_bod <- data.frame(PC = Body_axis$scores[,1], body_size = trait_matrix_loco$Mass, axis = "Size PCA", component = "PC1")



dat_text2 <- data.frame(axis = c("Beak PCA", "Locomotory PCA", "Beak PCA", "Locomotory PCA"), 
                        component = c("PC1", "PC2", "PC2", "PC1"),
                        label = c("a","d", "c", "b"))

dev.size("px")
#[1] 909 780



PC1s_bod$d <- densCols(x = PC1s_bod$body_size, 
                           y = PC1s_bod$PC, 
                           colramp = colorRampPalette(viridis::viridis(10)))

ED_8 <- ggplot(PC1s_bod, aes(x = body_size, y = PC)) +
  geom_point(col = "grey80") + 
  scale_color_identity(aes(fill = d)) +
  geom_smooth(aes(body_size, PC), data = PC1s_bod %>% filter(component == "PC1"), method = "lm", se = TRUE, col = "red") +
  facet_grid(component ~ axis) +
  ylab("Eigenvalue") +
  xlab("Body mass") +
  stat_cor(aes(label = ..r.label..),
           r.accuracy = 0.001,
           label.x = 2, label.y = -4.2, size = 8) +
  theme_bw() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        strip.background = element_blank(), #rect(fill = "grey"),
        strip.text = element_text(size = 25, face = "bold"),
        strip.text.y = element_text(vjust = 2),
        panel.grid = element_blank()) +
  geom_text(inherit.aes = FALSE,
            data    = dat_text2,
            mapping = aes(x = -Inf, y = -Inf, label = label),
            hjust   = -1.5,
            vjust   = -12,
            size = 9,
            fontface = "bold",
            col = "black") 


png("figures/ED 8 - PCA x Mass.png",
    width = 1025/1.25, height = 895/1.25)

ED_8

dev.off()





ggplot(bod_bod, aes(x = body_size, y = PC)) + geom_point() + 
  #facet_grid(component ~ axis) +
  ylab("Size Axis") +
  xlab("log(Body Mass)") +
  stat_cor(aes(label = ..r.label..),
           label.x = 2.75, label.y = -4.2, size = 8) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        strip.background = element_rect(fill = "grey"),
        strip.text = element_text(size = 17))

