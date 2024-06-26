  #'  =============================================================
  #'  Code from Bassing et al. (2024) Global Ecology & Conservation
  #'  
  #'  Co-Occurrence plots 
  #'  Washington Predator-Prey Project
  #'  =============================================================
  #'  Predict and plot conditional occupancy probabilities for interacting species
  #'  based on significant co-occurrence results from multi-species occupancy models.
  #'  
  #'  NOTE: Need to create a folder within working directory called "Figures" to 
  #'  save co-occurrence plots.
  #'  =============================================================
  
  #'  Load libraries
  library(unmarked)
  library(ggplot2)
  library(tidyverse)
  library(khroma)
  library(patchwork)
  
  #'  Read in model outputs
  load(file = "./MultiSpp_CoOcc_Model_Outputs.RData")
  
  scale_cov <- function(cov) {
    #'  Identify range of covariate of interest
    r <- range(cov)
    #'  Create sequence of values starting and ending with range values
    x <- seq(r[1], r[2], length.out = 100)
    #'  Scale values based on covariate mean and standard deviation
    x_scaled <- (x - mean(cov)) / (sd(cov))
    #'  Make single data frame of scaled and unscaled covariate
    x_cov <- as.data.frame(cbind(x, x_scaled))
    
    return(x_cov)
  }
  scaled_elev <- scale_cov(stations$Elev)
  scaled_forest <- scale_cov(stations$PercForest)
  scaled_graze <- scale_cov(stations$GrazingActivity)
  scaled_hunt <- scale_cov(stations$HuntingActivity)
  
  
  #'  -----------------------------------------------------
  ####  Predict effect of cattle on species interactions  ####
  #'  -----------------------------------------------------
  #'  Function to predict species interactions in response to covariate of interest
  #'  Note: nsims = number of bootstraps
  spp_interactions_g <- function(mod, elev, act, forest, area, spp1, spp2, cov) {  
    
    #'  Create data frame using the scaled covariate of interest while holding
    #'  all others at their mean (0 when scaled) or desired category (0 or 1)
    cov_df <- data.frame(Elev = elev, GrazingActivity = act, PercForest = forest,
                         Study_Area = area)  
    #'  Create characters for each species that include a "-", necessary for cond 
    #'  argument in predict when species is not present
    no_spp2 <- paste0("-",spp2); no_spp1 <- paste0("-",spp1)
    
    #'  Predict conditional occupancy when spp2 is absent 
    spp2_absent <- predict(mod, type = "state", species = spp1, cond = no_spp2,
                           newdata = cov_df, se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp1,
             Interaction = paste0(spp2, " absent"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    #'  Predict conditional occupancy when spp2 is present 
    spp2_present <- predict(mod, type = "state", species = spp1, cond = spp2,
                            newdata = cov_df, se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp1,
             Interaction = paste0(spp2, " present"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    
    #'  Predict conditional occupancy when spp1 is absent 
    spp1_absent <- predict(mod, type = "state", species = spp2, cond = no_spp1,
                           newdata = cov_df, se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp2,
             Interaction = paste0(spp1, " absent"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    #'  Predict conditional occupancy when spp1 is present 
    spp1_present <- predict(mod, type = "state", species = spp2, cond = spp1,
                            newdata = cov_df, se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp2,
             Interaction = paste0(spp1, " present"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    
    #'  Create one large data frame with all marginal probabilities
    sppX_df <- rbind(spp2_absent, spp2_present, spp1_absent, spp1_present)
    
    return(sppX_df)
  }
  #'  Predict cond. occupancy for each pairwise interaction over range of cattle activity
  #'  NOTE: setting all other continuous variables to 0 (their mean value) and 
  #'  setting categorical variables as either 0 or 1
  #'  For elk, study area = 0 (NE) because there were so few detections in OK
  #'  For md, wtd, and moose, study area = 1 (OK) b/c that is where the bulk of
  #'  cattle activity occurred and there were sufficient detections of these species
  sppX_coug_wtd_cattle <- spp_interactions_g(gs_cougwtd_graze2, elev = 0, act = scaled_graze[,2], 
                                             forest = 0, area = 1, spp1 = "cougar",  
                                             spp2 = "wtd", cov = scaled_graze[,1])
  sppX_bob_md_cattle <- spp_interactions_g(gs_bobmd_graze2, elev = 0, act = scaled_graze[,2], 
                                           forest = 0, area = 1, spp1 = "bobcat",  
                                           spp2 = "mule_deer", cov = scaled_graze[,1])
  sppX_bob_wtd_cattle <- spp_interactions_g(gs_bobwtd_graze2, elev = 0, act = scaled_graze[,2], 
                                            forest = 0, area = 1, spp1 = "bobcat", 
                                            spp2 = "wtd", cov = scaled_graze[,1]) 
  sppX_coy_md_cattle <- spp_interactions_g(gs_coymd_graze2, elev = 0, act = scaled_graze[,2], 
                                           forest = 0, area = 1, spp1 = "coyote", 
                                           spp2 = "mule_deer", cov = scaled_graze[,1]) 
  sppX_coy_wtd_cattle <- spp_interactions_g(gs_coywtd_graze2, elev = 0, act = scaled_graze[,2], 
                                            forest = 0, area = 1, spp1 = "coyote", 
                                            spp2 = "wtd", cov = scaled_graze[,1]) 
  
  #'  Save these since they take so long to generate!
  sppX_cattle_list <- list(sppX_coug_wtd_cattle, sppX_bob_md_cattle, sppX_bob_wtd_cattle, 
                           sppX_coy_md_cattle, sppX_coy_wtd_cattle)
  names(sppX_cattle_list) <- c("sppX_coug_wtd_cattle", "sppX_bob_md_cattle", "sppX_bob_wtd_cattle", 
                               "sppX_coy_md_cattle", "sppX_coy_wtd_cattle")
  save(sppX_cattle_list, file = "./sppX_cattle_for_visualizing.RData")
  
  
  #'  ------------------------------------------------------
  ####  Predict effect of hunters on species interactions  ####
  #'  ------------------------------------------------------
  #'  Function to predict species interactions in response to covariate of interest
  #'  Note: nsims = number of bootstraps
  spp_interactions_h <- function(mod, elev, act, forest, area, spp1, spp2, cov) {  
    
    #'  Create data frame using the scaled covariate of interest while holding
    #'  all others at their mean (0 when scaled) or desired category (0 or 1)
    cov_df <- data.frame(Elev = elev, HuntingActivity = act, PercForest = forest,
                         Study_Area = area) 
    #Public = factor(pub, levels = c(0, 1)), Study_Area = factor(area, levels = c(0, 1)))
    #'  Create characters for each species that include a "-", necessary for cond 
    #'  argument in predict when species is not present
    no_spp2 <- paste0("-",spp2); no_spp1 <- paste0("-",spp1)
    
    #'  Predict conditional occupancy when spp2 is absent
    spp2_absent <- predict(mod, type = "state", species = spp1, cond = no_spp2,
                           newdata = cov_df, se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp1,
             Interaction = paste0(spp2, " absent"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    #'  Predict conditional occupancy when spp2 is present
    spp2_present <- predict(mod, type = "state", species = spp1, cond = spp2,
                            newdata = cov_df, se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp1,
             Interaction = paste0(spp2, " present"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    
    #'  Predict conditional occupancy when spp1 is absent
    spp1_absent <- predict(mod, type = "state", species = spp2, cond = no_spp1,
                           newdata = cov_df, se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp2,
             Interaction = paste0(spp1, " absent"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    #'  Predict conditional occupancy when spp1 is present
    spp1_present <- predict(mod, type = "state", species = spp2, cond = spp1,
                            newdata = cov_df, se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp2,
             Interaction = paste0(spp1, " present"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    
    #'  Create one large data frame with all marginal probabilities
    sppX_df <- rbind(spp2_absent, spp2_present, spp1_absent, spp1_present)
    
    return(sppX_df)
  }
  #'  Predict cond. occupancy for each pairwise interaction over range of hunter activity
  #'  NOTE: setting all other continuous variables to 0 (their mean value) and 
  #'  setting categorical variables as either 0 or 1
  #'  For elk, moose, & wtd study area = 0 (NE) vs md study area = 1 (OK) because 
  #'  higher species-specific densities in the respective study areas
  sppX_wolf_moose_hunter <- spp_interactions_h(hs_wolfmoose_hunt2, elev = 0, act = scaled_hunt[,2], 
                                               forest = 0, area = 0, spp1 = "wolf",   
                                               spp2 = "moose", cov = scaled_hunt[,1])
  
  #'  Save these since they take so long to generate!
  save(sppX_wolf_moose_hunter, file = "./sppX_hunter_for_visualizing.RData")
  
  
    #'  -------------------------------------
  #####  Visualize conditional occupancy  #####
  #'  -------------------------------------
  #'  Load predicted conditional occupancy
  load("./sppX_cattle_for_visualizing.RData")
  load("./sppX_hunter_for_visualizing.RData")
  load("./sppX_landownership_for_visualizing.RData")
  
  ######  Cattle activity plots  ######
  #'  ----------------------------
  #'  Cougar-wtd co-occurrence with cattle activity
  coug_wtd_data <- sppX_cattle_list$sppX_coug_wtd_cattle %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("wtd absent", "wtd present", "cougar absent", "cougar present")))
  newlabs <- c("cougar" = "Cougar", "wtd" = "White-tailed deer") 
  coug_wtd_graze_facet <- ggplot(coug_wtd_data, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("cougar absent" = "Cougar absent", "cougar present" = "Cougar present", 
                                   "wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("cougar absent" = "Cougar absent", "cougar present" = "Cougar present",
                                 "wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom",
          text = element_text(size = 14)) +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Conditional occupancy") +
    ggtitle("Cougar - white-tailed deer co-occurrence, grazing season")
  coug_wtd_graze_facet
  
  # ggsave("./Figures/OccX_coug_wtd_graze.tiff", coug_wtd_graze_facet, 
  #        units = "in", width = 8, height = 6, dpi = 600, device = 'tiff', compression = 'lzw') 
  
  #'  Coyote-md co-occurrence with cattle activity
  coy_md_data <- sppX_cattle_list$sppX_coy_md_cattle %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("mule_deer absent", "mule_deer present", "coyote absent", "coyote present")),
           Interaction = gsub("_", " ", Interaction),
           InterXSpp = gsub("_", " ", InterXSpp),
           Species = gsub("_", " ", Species))
  newlabs <- c("coyote" = "Coyote", "mule deer" = "Mule deer") 
  coy_md_graze_facet <- ggplot(coy_md_data, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("coyote absent" = "Coyote \nabsent", "coyote present" = "Coyote \npresent", 
                                   "mule deer absent" = "Mule deer \nabsent", "mule deer present" = "Mule deer \npresent"), name = "Species interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("coyote absent" = "Coyote \nabsent", "coyote present" = "Coyote \npresent",
                                 "mule deer absent" = "Mule deer \nabsent", "mule deer present" = "Mule deer \npresent"), name = "Species interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom",
          text = element_text(size = 14)) +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Conditional occupancy") +
    ggtitle("Coyote - mule deer co-occurrence")
  coy_md_graze_facet
  
  # ggsave("./Figures/OccX_coy_md_graze.tiff", coy_md_graze_facet, 
  #        units = "in", width = 8, height = 6, dpi = 600, device = 'tiff', compression = 'lzw')
  
  #'  Coyote-wtd co-occurrence with cattle activity
  coy_wtd_data <- sppX_cattle_list$sppX_coy_wtd_cattle %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("wtd absent", "wtd present", "coyote absent", "coyote present")))
  newlabs <- c("coyote" = "Coyote", "wtd" = "White-tailed deer") 
  coy_wtd_graze_facet <- ggplot(coy_wtd_data, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("coyote absent" = "Coyote \nabsent", "coyote present" = "Coyote \npresent", 
                                   "wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("coyote absent" = "Coyote \nabsent", "coyote present" = "Coyote \npresent",
                                 "wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom",
          text = element_text(size = 14)) +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Conditional occupancy") +
    ggtitle("Coyote - white-tailed deer co-occurrence")
  coy_wtd_graze_facet
  
  # ggsave("./Figures/OccX_coy_wtd_graze.tiff", coy_wtd_graze_facet, 
  #        units = "in", width = 8, height = 6, dpi = 600, device = 'tiff', compression = 'lzw')
  
  #'  Bobcat-md co-occurrence with cattle activity
  bob_md_data <- sppX_cattle_list$sppX_bob_md_cattle %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("mule_deer absent", "mule_deer present", "bobcat absent", "bobcat present")),
           Interaction = gsub("_", " ", Interaction),
           InterXSpp = gsub("_", " ", InterXSpp),
           Species = gsub("_", " ", Species))
  newlabs <- c("bobcat" = "Bobcat", "mule deer" = "Mule deer") 
  bob_md_graze_facet <- ggplot(bob_md_data, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("bobcat absent" = "Bobcat \nabsent", "bobcat present" = "Bobcat \npresent", 
                                   "mule deer absent" = "Mule deer \nabsent", "mule deer present" = "Mule deer \npresent"), name = "Species interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("bobcat absent" = "Bobcat \nabsent", "bobcat present" = "Bobcat \npresent",
                                 "mule deer absent" = "Mule deer \nabsent", "mule deer present" = "Mule deer \npresent"), name = "Species interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom",
          text = element_text(size = 14),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    # xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Conditional occupancy") +
    ggtitle("Bobcat - mule deer co-occurrence")
  bob_md_graze_facet
  
  # ggsave("./Figures/OccX_bob_md_graze.tiff", bob_md_graze_facet, 
  #        units = "in", width = 8, height = 6, dpi = 600, device = 'tiff', compression = 'lzw')
  
  #'  Bobcat-wtd co-occurrence with cattle activity
  bob_wtd_data <- sppX_cattle_list$sppX_bob_wtd_cattle %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("wtd absent", "wtd present", "bobcat absent", "bobcat present")))
  newlabs <- c("bobcat" = "Bobcat", "wtd" = "White-tailed deer") 
  bob_wtd_graze_facet <- ggplot(bob_wtd_data, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("bobcat absent" = "Bobcat \nabsent", "bobcat present" = "Bobcat \npresent", 
                                   "wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("bobcat absent" = "Bobcat \nabsent", "bobcat present" = "Bobcat \npresent",
                                 "wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom",
          text = element_text(size = 14),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    # xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Conditional occupancy") +
    ggtitle("Bobcat - white-tailed deer co-occurrence")
  bob_wtd_graze_facet
  
  # ggsave("./Figures/OccX_bob_wtd_graze.tiff", bob_wtd_graze_facet, 
  #        units = "in", width = 8, height = 6, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  ###### Hunter activity plots  ######
  #'  ---------------------------
  #'  Wolf-moose co-occurrence with hunter activity
  wolf_moose_data <-  sppX_wolf_moose_hunter %>%  
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("moose absent", "moose present", "wolf absent", "wolf present")),
           Species = factor(Species, levels = c("wolf", "moose")))
  newlabs <- c("wolf" = "Wolf", "moose" = "Moose") 
  wolf_moose_hunt_facet <- ggplot(wolf_moose_data, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                   "moose absent" = "Moose absent", "moose present" = "Moose present"), name = "Species interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                 "moose absent" = "Moose absent", "moose present" = "Moose present"), name = "Species interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom",
          text = element_text(size = 14)) +
    xlab("Hunter activity (hunter detections/day)") + 
    ylab("Conditional occupancy") +
    ggtitle("Wolf - moose co-occurrence, hunting season")
  wolf_moose_hunt_facet
  
  # ggsave("./Figures/OccX_wolf_moose_hunt.tiff", wolf_moose_hunt_facet, 
  #        units = "in", width = 8, height = 6, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  ###### Patchwork plots together  ######
  #'  ------------------------------
  coocc_patwork <- coug_wtd_graze_facet / wolf_moose_hunt_facet + 
    plot_annotation(tag_levels = 'a') & theme(text = element_text(size = 22))
  coocc_patwork
  coocc_patwork_meso_wtd <- bob_wtd_graze_facet / coy_wtd_graze_facet + 
    plot_annotation(tag_levels = 'a') & theme(text = element_text(size = 22))
  coocc_patwork_meso_wtd
  coocc_patwork_meso_md <- bob_md_graze_facet / coy_md_graze_facet + 
    plot_annotation(tag_levels = 'a') & theme(text = element_text(size = 22))
  coocc_patwork_meso_md
  
  # ggsave("./Figures/OccX_pred-prey_cattl&hunter.tiff", coocc_patwork, 
  #        units = "in", width = 12, height = 16, dpi = 600, device = 'tiff', compression = 'lzw')
  # ggsave("./Figures/OccX_meso-wtd_cattle.tiff", coocc_patwork_meso_wtd, 
  #        units = "in", width = 12, height = 16, dpi = 600, device = 'tiff', compression = 'lzw')
  # ggsave("./Figures/OccX_meso-md_cattle.tiff", coocc_patwork_meso_md, 
  #        units = "in", width = 12, height = 16, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  
  #'  ------------------------------------------------------------
  ####  Predict effect of landownership on species interactions  ####
  #'  ------------------------------------------------------------
  #'  Function to predict species interactions in response to covariate of interest
  #'  Note: nsims = number of bootstraps
  spp_interactions_pub <- function(mod, elev, act, forest, pub, area, spp1, spp2, cov) {  
    
    #'  Create data frame using the scaled covariate of interest while holding
    #'  all others at their mean (0 when scaled) or desired category (0 or 1)
    cov_df <- data.frame(Elev = elev, HuntingActivity = act, PercForest = forest,
                         Public = factor(pub, levels = c(0, 1)), Study_Area = area)
    #'  Create characters for each species that include a "-", necessary for cond 
    #'  argument in predict when species is not present
    no_spp2 <- paste0("-",spp2); no_spp1 <- paste0("-",spp1)
    
    #'  Predict conditional occupancy when spp2 is absent
    spp2_absent <- predict(mod, type = "state", species = spp1, cond = no_spp2,
                           newdata = cov_df, se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp1,
             Interaction = paste0(spp2, " absent"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    #'  Predict conditional occupancy when spp2 is present
    spp2_present <- predict(mod, type = "state", species = spp1, cond = spp2,
                            newdata = cov_df, se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp1,
             Interaction = paste0(spp2, " present"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    
    #'  Predict conditional occupancy when spp1 is absent
    spp1_absent <- predict(mod, type = "state", species = spp2, cond = no_spp1,
                           newdata = cov_df, se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp2,
             Interaction = paste0(spp1, " absent"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    #'  Predict conditional occupancy when spp1 is present
    spp1_present <- predict(mod, type = "state", species = spp2, cond = spp1,
                            newdata = cov_df, se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp2,
             Interaction = paste0(spp1, " present"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    
    #'  Create one large data frame with all marginal probabilities
    sppX_df <- rbind(spp2_absent, spp2_present, spp1_absent, spp1_present)
    
    return(sppX_df)
  }
  
  #'  Predict cond. occupancy for each pairwise interaction on public vs private land
  #'  NOTE: setting all continuous variables to 0 (their mean value) and setting 
  #'  other categorical variables as either 0 or 1
  #'  For elk, moose, & wtd study area = 0 (NE) vs md study area = 1 (OK) because 
  #'  higher species-specific densities in the respective study areas
  sppX_coug_md_pub <- spp_interactions_pub(hs_cougmd_pub2, elev = 0, act = 0,
                                           forest = 0, pub = 1, area = 1, spp1 = "cougar",
                                           spp2 = "muledeer", cov = "Public")
  sppX_coug_md_priv <- spp_interactions_pub(hs_cougmd_pub2, elev = 0, act = 0,
                                            forest = 0, pub = 0, area = 1, spp1 = "cougar",
                                            spp2 = "muledeer", cov = "Private")
  sppX_coug_wtd_pub <- spp_interactions_pub(hs_cougwtd_pub2, elev = 0, act = 0,
                                            forest = 0, pub = 1, area = 0, spp1 = "cougar",
                                            spp2 = "wtd", cov = "Public")
  sppX_coug_wtd_priv <- spp_interactions_pub(hs_cougwtd_pub2, elev = 0, act = 0,
                                             forest = 0, pub = 0, area = 0, spp1 = "cougar",
                                             spp2 = "wtd", cov = "Private")
  sppX_bob_wtd_pub <- spp_interactions_pub(hs_bobwtd_pub2, elev = 0, act = 0,
                                           forest = 0, pub = 1, area = 0, spp1 = "bobcat",
                                           spp2 = "wtd", cov = "Public")
  sppX_bob_wtd_priv <- spp_interactions_pub(hs_bobwtd_pub2, elev = 0, act = 0,
                                            forest = 0, pub = 0, area = 0, spp1 = "bobcat",
                                            spp2 = "wtd", cov = "Private")
  sppX_coy_wtd_pub <- spp_interactions_pub(hs_coywtd_pub2, elev = 0, act = 0,
                                           forest = 0, pub = 1, area = 0, spp1 = "coyote",
                                           spp2 = "wtd", cov = "Public")
  sppX_coy_wtd_priv <- spp_interactions_pub(hs_coywtd_pub2, elev = 0, act = 0,
                                            forest = 0, pub = 0, area = 0, spp1 = "coyote",
                                            spp2 = "wtd", cov = "Private")
  
  #'  Save these since they take so long to generate!
  sppX_coug_md_landownership <- rbind(sppX_coug_md_pub, sppX_coug_md_priv)
  sppX_coug_wtd_landownership <- rbind(sppX_coug_wtd_pub, sppX_coug_wtd_priv)
  sppX_bob_wtd_landownership <- rbind(sppX_bob_wtd_pub, sppX_bob_wtd_priv)
  sppX_coy_wtd_landownership <- rbind(sppX_coy_wtd_pub, sppX_coy_wtd_priv)
  
  sppX_landownership_list <- list(sppX_coug_md_landownership, sppX_coug_wtd_landownership, sppX_bob_wtd_landownership, sppX_coy_wtd_landownership)
  names(sppX_landownership_list) <- c("sppX_coug_md_landownership", "sppX_coug_wtd_landowner", "sppX_bob_wtd_landownership", "sppX_coy_wtd_landowner")
  save(sppX_landownership_list, file = "./sppX_landownership_for_visualizing.RData")
  
  
  ##### Land ownership plots  #####
  #' --------------------------
  #'  Cougar-md co-occurrence on private vs public land
  coug_md_land_data <- sppX_landownership_list$sppX_coug_md_landownership %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("muledeer absent", "muledeer present", "cougar absent", "cougar present")))
  newlabs <- c("cougar" = "Cougar", "muledeer" = "Mule deer") 
  coug_md_land_facet <- ggplot(coug_md_land_data, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_errorbar(aes(ymin=lower, ymax = upper, color = Interaction), width = 0, position = position_dodge(width = 0.4)) +
    scale_fill_bright(labels = c("cougar absent" = "Cougar absent", "cougar present" = "Cougar present",
                                 "muledeer absent" = "Mule deer absent", "muledeer present" = "Mule deer present"), name = "Species Interaction") +
    geom_point(stat = "identity", aes(col = Interaction), size = 2.5, position = position_dodge(width = 0.4)) +
    scale_colour_bright(labels = c("cougar absent" = "Cougar absent", "cougar present" = "Cougar present", 
                                   "muledeer absent" = "Mule deer absent", "muledeer present" = "Mule deer present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Land ownership") + 
    ylab("Conditional occupancy") +
    ggtitle("Cougar - mule deer co-occurrence, hunting season")
  coug_md_land_facet
  
  # ggsave("./Figures/OccX_coug_md_landownership.tiff", coug_md_land_facet, 
  #        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw') 
  
  #'  Cougar-wtd co-occurrence on private vs public land
  coug_wtd_land_data <- sppX_landownership_list$sppX_coug_wtd_landownership %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("wtd absent", "wtd present", "cougar absent", "cougar present"))) 
  newlabs <- c("cougar" = "Cougar", "wtd" = "White-tailed deer") 
  coug_wtd_land_facet <- ggplot(coug_wtd_land_data, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_errorbar(aes(ymin=lower, ymax = upper, color = Interaction), width = 0, position = position_dodge(width = 0.4)) +
    scale_fill_bright(labels = c("cougar absent" = "Cougar absent", "cougar present" = "Cougar present",
                                 "wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species Interaction") +
    geom_point(stat = "identity", aes(col = Interaction), size = 2.5, position = position_dodge(width = 0.4)) +
    scale_colour_bright(labels = c("cougar absent" = "Cougar absent", "cougar present" = "Cougar present", 
                                   "wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Land ownership") + 
    ylab("Conditional occupancy") +
    ggtitle("Cougar - white-tailed deer co-occurrence, hunting season")
  coug_wtd_land_facet
  
  # ggsave("./Figures/OccX_coug_wtd_landownership.tiff", coug_wtd_land_facet, 
  #        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw') 
  
  #'  Bobcat-wtd co-occurrence on private vs public land
  bob_wtd_land_data <- sppX_landownership_list$sppX_bob_wtd_landownership %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("wtd absent", "wtd present", "bobcat absent", "bobcat present")))
  newlabs <- c("bobcat" = "Bobcat", "wtd" = "White-tailed deer") 
  bob_wtd_land_facet <- ggplot(bob_wtd_land_data, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_errorbar(aes(ymin=lower, ymax = upper, color = Interaction), width = 0, position = position_dodge(width = 0.4)) +
    scale_fill_bright(labels = c("bobcat absent" = "Bobcat absent", "bobcat present" = "Bobcat present",
                                 "wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species Interaction") +
    geom_point(stat = "identity", aes(col = Interaction), size = 2.5, position = position_dodge(width = 0.4)) +
    scale_colour_bright(labels = c("bobcat absent" = "Bobcat absent", "bobcat present" = "Bobcat present", 
                                   "wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Land ownership") + 
    ylab("Conditional occupancy") +
    ggtitle("Bobcat - white-tailed deer co-occurrence, hunting season")
  bob_wtd_land_facet
  
  # ggsave("./Figures/OccX_bob_wtd_landownership.tiff", bob_wtd_land_facet, 
  #        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw') 
  
  #'  Coyote-wtd co-occurrence on private vs public land
  coy_wtd_land_data <- sppX_landownership_list$sppX_coy_wtd_landownership %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("wtd absent", "wtd present", "coyote absent", "coyote present")))
  newlabs <- c("coyote" = "Coyote", "wtd" = "White-tailed deer") 
  coy_wtd_land_facet <- ggplot(coy_wtd_land_data, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_errorbar(aes(ymin=lower, ymax = upper, color = Interaction), width = 0, position = position_dodge(width = 0.4)) +
    scale_fill_bright(labels = c("coyote absent" = "Coyote absent", "coyote present" = "Coyote present",
                                 "wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species Interaction") +
    geom_point(stat = "identity", aes(col = Interaction), size = 2.5, position = position_dodge(width = 0.4)) +
    scale_colour_bright(labels = c("coyote absent" = "Coyote absent", "coyote present" = "Coyote present", 
                                   "wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Land ownership") + 
    ylab("Conditional occupancy") +
    ggtitle("Coyote - white-tailed deer co-occurrence, hunting season")
  coy_wtd_land_facet
  
  # ggsave("./Figures/OccX_coy_wtd_landownership.tiff", coy_wtd_land_facet, 
  #        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw') 
  
    
  
