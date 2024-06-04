  #'  =============================================================
  #'  Code from Bassing et al. (2024) Global Ecology & Conservation
  #'  
  #'  Multi-Species Co-Occurance Models
  #'  Washington Predator-Prey Project
  #'  =============================================================
  #'  Script to run multi-species, single-season occupancy models for deer, elk, 
  #'  moose, black bears, cougars, wolves, coyotes, and bobcats during the grazing 
  #'  season 2018-2020 (July - Sept) and hunting season 2018-2020 (Oct - Nov), 
  #'  respectively. Grazing season co-occurrence models include 13 7-day sampling 
  #'  occasions comprising the peak of livestock activity detected on camera. 
  #'  Hunting season co-occurrence models include 8 7-day sampling occasions 
  #'  comprising the two general rifle hunting seasons in eastern Washington. 
  #'  Co-occurrence models test whether predator-prey co-occurrence is not 
  #'  independent and whether their occurrence, co-occurrence, and detection are 
  #'  influenced by livestock and hunting activity at camera sites.
  #'  
  #'  Encounter histories are generated with the Detection_histories_for_unmarked.R
  #'  and Data_formatting_cattle_hunter_activity.R scripts. Covariate data formatted 
  #'  with the Data_formatting_for_co-occ_model.R script.
  #'  
  #'  NOTE: Need to create a folder within working directory called "Tables" to 
  #'  save occupancy model result tables.
  #'  =============================================================
  
  #'  Clean workspace & load libraries
  rm(list = ls())
  
  library(unmarked)
  library(tidyverse)
  
  #'  Source scripts that 
  #'  1. Generate detection histories for each species (predator & prey)
  source("./Scripts/Sourcing_scripts/1_Detection_histories_for_unmarked.R")
  
  #'  2. Generate detection histories representing grazing and hunting activity 
  #'  Note: 5 min intervals between independent detections of cattle or humans
  source("./Scripts/Sourcing_scripts/2_Data_formatting_cattle_hunter_activity.R")
  
  #'  3. Format covariate data and detection histories for multi-species occupancy 
  #'  models in unmarked with two-species interactions
  source("./Scripts/Sourcing_scripts/3_Data_formatting_for_co-occ_model.R")
  
  
  ####  Multi-Species Occupancy models  ####
  #'  ==================================
  #'  Multi-species occupancy model --> occuMulti (pg 83) in unmarked manual
  #'  Occupancy formulas should match number/order of columns in fDesign matrix
  #'  i.e., one formula for each species/interaction of interest
  #'  Detection formulas should match number/order of species in list of DH
  #'  Use ~1 to estimate intercept without covariates
  #'  Use ~0 to fix a natural parameter to 0
  #'  E.g., occFormulas <- c("~1", "~1", "~0") estimates intercept for 1st order 
  #'  natural parameters (2 spp) but fixes 2nd order natural parameter to 0.
  #'  Covariates: Can use different covariates on different natural parameters, 
  #'  E.g., covs on 1st order parameters to explain single-spp occurrence 
  #'  regardless of other spp, covs on 2nd order parameters to explain co-occ
  #'  
  #'  Testing hypothesis that co-occurrence is non-independent and that cattle/
  #'  hunter activity impacts occurrence and/or co-occurrence patterns between
  #'  predators and prey.
  #'  
  #'  Include a consistent set of additional covariates to account for habitat
  #'  variation and other factors we know influence occurrence and detection.
  #'  =============================
  
  ####  GRAZING SEASON MODELS  ####
  #'  -------------------------
  #'  Detection formulas
  detFormulas_graze <- c("~Trail + WeeklyGrazing", "~Trail + WeeklyGrazing")
  detFormulas_allot <- c("~Trail + PublicGrazing", "~Trail + PublicGrazing")
  
  #'  Occupancy formulas (decreasing in complexity from effect of cattle activity 
  #'  on co-occurrence (graze2) to no effect but co-occurrence (graze1) to 
  #'  independent occurrence (graze0))
  occFormulas_graze2 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + GrazingActivity", 
                          "~Study_Area + Elev + I(Elev^2) + PercForest + GrazingActivity", 
                          "~GrazingActivity")
  occFormulas_graze1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + GrazingActivity", 
                          "~Study_Area + Elev + I(Elev^2) + PercForest + GrazingActivity", 
                          "~1")
  occFormulas_graze0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + GrazingActivity",
                          "~Study_Area + Elev + I(Elev^2) + PercForest + GrazingActivity",
                          "~0")
  #'  Occupancy formulas for additional exploration of grazing allotment effects
  occFormulas_allot2 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + PublicGrazing", 
                          "~Study_Area + Elev + I(Elev^2) + PercForest + PublicGrazing", 
                          "~PublicGrazing")
  occFormulas_allot1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + PublicGrazing", 
                          "~Study_Area + Elev + I(Elev^2) + PercForest + PublicGrazing", 
                          "~1")
  occFormulas_allot0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + PublicGrazing",
                          "~Study_Area + Elev + I(Elev^2) + PercForest + PublicGrazing",
                          "~0")
  
  #####  Cougar-Mule Deer Grazing Season  #####
  (gs_cougmd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coug_md_grazing_UMF, silent = TRUE))
  (gs_cougmd_allot2 <- occuMulti(detFormulas_allot, occFormulas_allot2, coug_md_grazing_UMF, silent = TRUE))
   
  #####  Cougar-ELK Grazing Season  #####
  (gs_cougelk_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coug_elk_grazing_UMF, silent = TRUE))
  # (gs_cougelk_allot2 <- occuMulti(detFormulas_allot, occFormulas_allot2, coug_elk_grazing_UMF, silent = TRUE)) # fails
  # (gs_cougelk_allot1 <- occuMulti(detFormulas_allot, occFormulas_allot1, coug_elk_grazing_UMF, silent = TRUE)) # fails
  (gs_cougelk_allot0 <- occuMulti(detFormulas_allot, occFormulas_allot0, coug_elk_grazing_UMF, silent = TRUE))  
    
  #####  Cougar-White-tailed Deer Grazing Season  #####
  (gs_cougwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coug_wtd_grazing_UMF, silent = TRUE))
  (gs_cougwtd_allot2 <- occuMulti(detFormulas_allot, occFormulas_allot2, coug_wtd_grazing_UMF, silent = TRUE))
    
  #####  Cougar-Moose Grazing Season  #####
  (gs_cougmoose_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coug_moose_grazing_UMF, silent = TRUE))
  (gs_cougmoose_allot2 <- occuMulti(detFormulas_allot, occFormulas_allot2, coug_moose_grazing_UMF, silent = TRUE))
    
  #####  Wolf-Mule Deer Grazing Season  #####
  (gs_wolfmd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, wolf_md_grazing_UMF, silent = TRUE))
  (gs_wolfmd_allot2 <- occuMulti(detFormulas_allot, occFormulas_allot2, wolf_md_grazing_UMF, silent = TRUE))
    
  #####  Wolf-Elk Grazing Season  #####
  # (gs_wolfelk_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, wolf_elk_grazing_UMF, silent = TRUE)) # fails
  (gs_wolfelk_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, wolf_elk_grazing_UMF, silent = TRUE))
  # (gs_wolfelk_allot2 <- occuMulti(detFormulas_allot, occFormulas_allot2, wolf_elk_grazing_UMF, silent = TRUE)) # fails
  (gs_wolfelk_allot1 <- occuMulti(detFormulas_allot, occFormulas_allot1, wolf_elk_grazing_UMF, silent = TRUE))
    
  #####  Wolf-White-tailed Deer Grazing Season  #####
  (gs_wolfwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, wolf_wtd_grazing_UMF, silent = TRUE))
  (gs_wolfwtd_allot2 <- occuMulti(detFormulas_allot, occFormulas_allot2, wolf_wtd_grazing_UMF, silent = TRUE))
  
  #####  Wolf-Moose Grazing Season  #####
  (gs_wolfmoose_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, wolf_moose_grazing_UMF, silent = TRUE))
  (gs_wolfmoose_allot2 <- occuMulti(detFormulas_allot, occFormulas_allot2, wolf_moose_grazing_UMF, silent = TRUE))
    
  #####  Bear-Mule Deer Grazing Season  #####
  (gs_bearmd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bear_md_grazing_UMF, silent = TRUE))
  (gs_bearmd_allot2 <- occuMulti(detFormulas_allot, occFormulas_allot2, bear_md_grazing_UMF, silent = TRUE))
    
  #####  Bear-Elk Grazing Season  #####
  # (gs_bearelk_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bear_elk_grazing_UMF, silent = TRUE)) # GrazingActivity on f12 not converging well
  (gs_bearelk_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, bear_elk_grazing_UMF, silent = TRUE)) 
  (gs_bearelk_allot2 <- occuMulti(detFormulas_allot, occFormulas_allot2, bear_elk_grazing_UMF, silent = TRUE)) 
      
  #####  Bear-White-tailed Deer Grazing Season  #####
  (gs_bearwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bear_wtd_grazing_UMF, silent = TRUE))
  (gs_bearwtd_allot2 <- occuMulti(detFormulas_allot, occFormulas_allot2, bear_wtd_grazing_UMF, silent = TRUE))
    
  #####  Bear-Moose Grazing Season  #####
  # (gs_bearmoose_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bear_moose_grazing_UMF, silent = TRUE)) # GrazingActivity on f2 & f12 not converging well
  (gs_bearmoose_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, bear_moose_grazing_UMF, silent = TRUE))
  (gs_bearmoose_allot2 <- occuMulti(detFormulas_allot, occFormulas_allot2, bear_moose_grazing_UMF, silent = TRUE)) 
    
  #####  Coyote-Mule Deer Grazing Season  #####
  (gs_coymd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coy_md_grazing_UMF, silent = TRUE))
  (gs_coymd_allot2 <- occuMulti(detFormulas_allot, occFormulas_allot2, coy_md_grazing_UMF, silent = TRUE))
  
  #####  Coyote-White-tailed Deer Grazing Season  #####
  (gs_coywtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coy_wtd_grazing_UMF, silent = TRUE))
  (gs_coywtd_allot2 <- occuMulti(detFormulas_allot, occFormulas_allot2, coy_wtd_grazing_UMF, silent = TRUE))
    
  #####  Bobcat-Mule Deer Grazing Season  #####
  (gs_bobmd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bob_md_grazing_UMF, silent = TRUE))
  (gs_bobmd_allot2 <- occuMulti(detFormulas_allot, occFormulas_allot2, bob_md_grazing_UMF, silent = TRUE))
    
  #####  Bobcat-White-tailed Deer Grazing Season  #####
  (gs_bobwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bob_wtd_grazing_UMF, silent = TRUE))
  (gs_bobwtd_allot2 <- occuMulti(detFormulas_allot, occFormulas_allot2, bob_wtd_grazing_UMF, silent = TRUE))
    
  
  ####  HUNTING SEASON MODELS  ####
  #'  -------------------------
  #'  Detection formulas
  detFormulas_pub <- c("~Trail + Public", "~Trail + Public")
  detFormulas_hunt <- c("~Trail + WeeklyHunting", "~Trail + WeeklyHunting") 
  #'  Remove Public vs Private covariate from predator sub-model
  detFormulas_pubish <- c("~Trail", "~Trail + Public")
  
  #'  Occupancy formulas (decreasing in complexity from effect of hunter activity 
  #'  on co-occurrence (hunt2) to no effect but co-occurrence (hunt1) to 
  #'  independent occurrence (hunt0))
  occFormulas_hunt2 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                         "~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                         "~HuntingActivity")
  occFormulas_hunt1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                         "~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                         "~1")
  occFormulas_hunt0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                         "~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                         "~0")
  #'  Occupancy formulas for additional exploration of landownership effects
  occFormulas_pub2 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                        "~Public")
  occFormulas_pub1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                        "~1")
  occFormulas_pub0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                        "~0")
  #'  Remove Public vs Private covariate from predator sub-model
  occFormulas_pubish2 <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                           "~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                           "~Public")
  occFormulas_pubish1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                           "~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                           "~1")
  occFormulas_pubish0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                           "~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                           "~0")
   
  
  #####  Cougar-Mule Deer Hunting Season  #####
  (hs_cougmd_hunt2 <- occuMulti(detFormulas_hunt, occFormulas_hunt2, coug_md_hunting_UMF, silent = TRUE)) 
  (hs_cougmd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pubish2, coug_md_hunting_UMF, silent = TRUE))
  
  #####  Cougar-Elk Hunting Season  #####
  # (hs_cougelk_hunt2 <- occuMulti(detFormulas_hunt, occFormulas_hunt2, coug_elk_hunting_UMF, silent = TRUE)) # fails
  # (hs_cougelk_hunt1 <- occuMulti(detFormulas_hunt, occFormulas_hunt1, coug_elk_hunting_UMF, silent = TRUE)) # fails
  (hs_cougelk_hunt0 <- occuMulti(detFormulas_hunt, occFormulas_hunt0, coug_elk_hunting_UMF, silent = TRUE)) 
  # (hs_cougelk_pub2 <- occuMulti(detFormulas_pub, occFormulas_pubish2, coug_elk_hunting_UMF, silent = TRUE)) # fails
  # (hs_cougelk_pub1 <- occuMulti(detFormulas_pub, occFormulas_pubish1, coug_elk_hunting_UMF, silent = TRUE)) # fails
  (hs_cougelk_pub0 <- occuMulti(detFormulas_pub, occFormulas_pubish0, coug_elk_hunting_UMF, silent = TRUE))
  
  #####  Cougar-White-tailed Deer Hunting Season  #####
  # (hs_cougwtd_hunt2 <- occuMulti(detFormulas_hunt, occFormulas_hunt2, coug_wtd_hunting_UMF, silent = TRUE)) # fails
  (hs_cougwtd_hunt1 <- occuMulti(detFormulas_hunt, occFormulas_hunt1, coug_wtd_hunting_UMF, silent = TRUE)) 
  (hs_cougwtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pubish2, coug_wtd_hunting_UMF, silent = TRUE))
  
  #####  Cougar-Moose Hunting Season  ###E#
  # (hs_cougmoose_hunt2 <- occuMulti(detFormulas_hunt, occFormulas_hunt2, coug_moose_hunting_UMF, silent = TRUE)) # fails
  # (hs_cougmoose_hunt1 <- occuMulti(detFormulas_hunt, occFormulas_hunt1, coug_moose_hunting_UMF, silent = TRUE)) # fails
  (hs_cougmoose_hunt0 <- occuMulti(detFormulas_hunt, occFormulas_hunt0, coug_moose_hunting_UMF, silent = TRUE)) 
  # (hs_cougmoose_pub2 <- occuMulti(detFormulas_pub, occFormulas_pubish2, coug_moose_hunting_UMF, silent = TRUE)) # fails
  # (hs_cougmoose_pub1 <- occuMulti(detFormulas_pub, occFormulas_pubish1, coug_moose_hunting_UMF, silent = TRUE)) # f2 & f12 intercepts not converging well
  (hs_cougmoose_pub0 <- occuMulti(detFormulas_pub, occFormulas_pubish0, coug_moose_hunting_UMF, silent = TRUE))
  
  #####  Wolf-Mule Deer Hunting Season  #####
  (hs_wolfmd_hunt2 <- occuMulti(detFormulas_hunt, occFormulas_hunt2, wolf_md_hunting_UMF, silent = TRUE))
  # (hs_wolfmd_pub2 <- occuMulti(detFormulas_pubish, occFormulas_pubish2, wolf_md_hunting_UMF, silent = TRUE)) # fails
  (hs_wolfmd_pub1 <- occuMulti(detFormulas_pubish, occFormulas_pubish1, wolf_md_hunting_UMF, silent = TRUE))
  
  #####  Wolf-Elk Hunting Season  #####
  (hs_wolfelk_hunt2 <- occuMulti(detFormulas_hunt, occFormulas_hunt2, wolf_elk_hunting_UMF, silent = TRUE))
  # (hs_wolfelk_pub2 <- occuMulti(detFormulas_pubish, occFormulas_pubish2, wolf_elk_hunting_UMF, silent = TRUE)) # fails
  (hs_wolfelk_pub1 <- occuMulti(detFormulas_pubish, occFormulas_pubish1, wolf_elk_hunting_UMF, silent = TRUE))
  
  #####  Wolf-White-tailed Deer Hunting Season  #####
  (hs_wolfwtd_hunt2 <- occuMulti(detFormulas_hunt, occFormulas_hunt2, wolf_wtd_hunting_UMF, silent = TRUE))  
  # (hs_wolfwtd_pub2 <- occuMulti(detFormulas_pubish, occFormulas_pubish2, wolf_wtd_hunting_UMF, silent = TRUE)) # fails
  (hs_wolfwtd_pub1 <- occuMulti(detFormulas_pubish, occFormulas_pubish1, wolf_wtd_hunting_UMF, silent = TRUE))
  
  #####  Wolf-Moose Hunting Season  #####
  (hs_wolfmoose_hunt2 <- occuMulti(detFormulas_hunt, occFormulas_hunt2, wolf_moose_hunting_UMF, silent = TRUE))
  # (hs_wolfmoose_pub2 <- occuMulti(detFormulas_pubish, occFormulas_pubish2, wolf_moose_hunting_UMF, silent = TRUE)) # fails
  (hs_wolfmoose_pub1 <- occuMulti(detFormulas_pubish, occFormulas_pubish1, wolf_moose_hunting_UMF, silent = TRUE))
   
  #####  Bear-Mule Deer Hunting Season  #####
  (hs_bearmd_hunt2 <- occuMulti(detFormulas_hunt, occFormulas_hunt2, bear_md_hunting_UMF, silent = TRUE))
  (hs_bearmd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bear_md_hunting_UMF, silent = TRUE))
  
  #####  Bear-Elk Hunting Season  #####
  (hs_bearelk_hunt2 <- occuMulti(detFormulas_hunt, occFormulas_hunt2, bear_elk_hunting_UMF, silent = TRUE))
  # (hs_bearelk_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bear_elk_hunting_UMF, silent = TRUE)) # fails
  (hs_bearelk_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bear_elk_hunting_UMF, silent = TRUE))
  
  #####  Bear-White-tailed Deer Hunting Season  #####
  (hs_bearwtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bear_wtd_hunting_UMF, silent = TRUE))
  (hs_bearwtd_hunt2 <- occuMulti(detFormulas_hunt, occFormulas_hunt2, bear_wtd_hunting_UMF, silent = TRUE)) 
  
  #####  Bear-Moose Hunting Season  #####
  (hs_bearmoose_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bear_moose_hunting_UMF, silent = TRUE)) # fails
  (hs_bearmoose_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bear_moose_hunting_UMF, silent = TRUE))
  (hs_bearmoose_hunt2 <- occuMulti(detFormulas_hunt, occFormulas_hunt2, bear_moose_hunting_UMF, silent = TRUE)) 
  
  #####  Coyote-Mule Deer Hunting Season  #####
  (hs_coymd_hunt2 <- occuMulti(detFormulas_hunt, occFormulas_hunt2, coy_md_hunting_UMF, silent = TRUE)) 
  (hs_coymd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, coy_md_hunting_UMF, silent = TRUE))
  
  #####  Coyote-White-tailed Deer Hunting Season  #####
  # (hs_coywtd_hunt2 <- occuMulti(detFormulas_hunt, occFormulas_hunt2, coy_wtd_hunting_UMF, silent = TRUE)) # fails
  (hs_coywtd_hunt1 <- occuMulti(detFormulas_hunt, occFormulas_hunt1, coy_wtd_hunting_UMF, silent = TRUE))
  (hs_coywtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, coy_wtd_hunting_UMF, silent = TRUE))
  
  #####  Bobcat-Mule Deer Hunting Season  #####
  (hs_bobmd_hunt2 <- occuMulti(detFormulas_hunt, occFormulas_hunt2, bob_md_hunting_UMF, silent = TRUE)) 
  (hs_bobmd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bob_md_hunting_UMF, silent = TRUE))
  
  #####  Bobcat-White-tailed Deer Hunting Season  #####
  # (hs_bobwtd_hunt2 <- occuMulti(detFormulas_hunt, occFormulas_hunt2, bob_wtd_hunting_UMF, silent = TRUE)) # HunterActivity not converging well
  (hs_bobwtd_hunt1 <- occuMulti(detFormulas_hunt, occFormulas_hunt1, bob_wtd_hunting_UMF, silent = TRUE))
  (hs_bobwtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bob_wtd_hunting_UMF, silent = TRUE))
  
  #' Save model outputs in one giant R image (helpful for creating figures later)
  save.image(file = "./MultiSpp_CoOcc_Model_Outputs.RData")
  
  
  ####  Summary tables  ####
  #'  ------------------
  #'  Save model outputs in table format 
  #'  Function to save occupancy results and append species/season info
  rounddig <- 2
  occ_out <- function(mod, spp1, spp2, season) {
    out <- summary(mod@estimates)$state %>%
      mutate(
        Parameter = row.names(summary(mod@estimates)$state),
        Species1 = rep(spp1, nrow(.)),
        Species2 = rep(spp2, nrow(.)),
        Season = rep(season, nrow(.)),
        Estimate = round(Estimate, rounddig),
        SE = round(SE, rounddig),
        z = round(z, rounddig),
        Pval = round(`P(>|z|)`, rounddig)
      ) %>%
      dplyr::select(-`P(>|z|)`) %>%
      relocate(Parameter, .before = Estimate) %>%
      relocate(Species1, .before = Parameter) %>%
      relocate(Species2, .before = Parameter) %>%
      relocate(Season, .before = Parameter) 
    
    return(out)
  }
  #####  Grazing season models - psi  #####
  occ_cougmd_grazing_allot <- occ_out(gs_cougmd_allot2, "Cougar", "Mule Deer", "Grazing")
  occ_cougmd_grazing_cattle <- occ_out(gs_cougmd_graze2, "Cougar", "Mule Deer", "Grazing")
  occ_cougelk_grazing_allot <- occ_out(gs_cougelk_allot0, "Cougar", "Elk", "Grazing")
  occ_cougelk_grazing_cattle <- occ_out(gs_cougelk_graze2, "Cougar", "Elk", "Grazing")
  occ_cougwtd_grazing_allot <- occ_out(gs_cougwtd_allot2, "Cougar", "White-tailed Deer", "Grazing")
  occ_cougwtd_grazing_cattle <- occ_out(gs_cougwtd_graze2, "Cougar", "White-tailed Deer", "Grazing")
  occ_cougmoose_grazing_allot <- occ_out(gs_cougmoose_allot2, "Cougar", "Moose", "Grazing")
  occ_cougmoose_grazing_cattle <- occ_out(gs_cougmoose_graze2, "Cougar", "Moose", "Grazing")
  occ_wolfmd_grazing_allot <- occ_out(gs_wolfmd_allot2, "Wolf", "Mule Deer", "Grazing")
  occ_wolfmd_grazing_cattle <- occ_out(gs_wolfmd_graze2, "Wolf", "Mule Deer", "Grazing")
  occ_wolfelk_grazing_allot <- occ_out(gs_wolfelk_allot1, "Wolf", "Elk", "Grazing")
  occ_wolfelk_grazing_cattle <- occ_out(gs_wolfelk_graze1, "Wolf", "Elk", "Grazing")
  occ_wolfwtd_grazing_allot <- occ_out(gs_wolfwtd_allot2, "Wolf", "White-tailed Deer", "Grazing")
  occ_wolfwtd_grazing_cattle <- occ_out(gs_wolfwtd_graze2, "Wolf", "White-tailed Deer", "Grazing")
  occ_wolfmoose_grazing_allot <- occ_out(gs_wolfmoose_allot2, "Wolf", "Moose", "Grazing")
  occ_wolfmoose_grazing_cattle <- occ_out(gs_wolfmoose_graze2, "Wolf", "Moose", "Grazing")
  occ_bearmd_grazing_allot <- occ_out(gs_bearmd_allot2, "Black Bear", "Mule Deer", "Grazing")
  occ_bearmd_grazing_cattle <- occ_out(gs_bearmd_graze2, "Black Bear", "Mule Deer", "Grazing")
  occ_bearelk_grazing_allot <- occ_out(gs_bearelk_allot2, "Black Bear", "Elk", "Grazing")
  occ_bearelk_grazing_cattle <- occ_out(gs_bearelk_graze1, "Black Bear", "Elk", "Grazing")
  occ_bearwtd_grazing_allot <- occ_out(gs_bearwtd_allot2, "Black Bear", "White-tailed Deer", "Grazing")
  occ_bearwtd_grazing_cattle <- occ_out(gs_bearwtd_graze2, "Black Bear", "White-tailed Deer", "Grazing")
  occ_bearmoose_grazing_allot <- occ_out(gs_bearmoose_allot2, "Black Bear", "Moose", "Grazing")
  occ_bearmoose_grazing_cattle <- occ_out(gs_bearmoose_graze1, "Black Bear", "Moose", "Grazing")
  occ_coymd_grazing_allot <- occ_out(gs_coymd_allot2, "Coyote", "Mule Deer", "Grazing")
  occ_coymd_grazing_cattle <- occ_out(gs_coymd_graze2, "Coyote", "Mule Deer", "Grazing")
  occ_coywtd_grazing_allot <- occ_out(gs_coywtd_allot2, "Coyote", "White-tailed Deer", "Grazing")
  occ_coywtd_grazing_cattle <- occ_out(gs_coywtd_graze2, "Coyote", "White-tailed Deer", "Grazing")
  occ_bobmd_grazing_allot <- occ_out(gs_bobmd_allot2, "Bobcat", "Mule Deer Deer", "Grazing")
  occ_bobmd_grazing_cattle <- occ_out(gs_bobmd_graze2, "Bobcat", "Mule Deer Deer", "Grazing")
  occ_bobwtd_grazing_allot <- occ_out(gs_bobwtd_allot2, "Bobcat", "White-tailed Deer", "Grazing")
  occ_bobwtd_grazing_cattle <- occ_out(gs_bobwtd_graze2, "Bobcat", "White-tailed Deer", "Grazing")
  #####  Hunting season models - psi  #####
  occ_cougmd_hunting_hunt <- occ_out(hs_cougmd_hunt2, "Cougar", "Mule Deer", "Hunting")
  occ_cougmd_hunting_pub <- occ_out(hs_cougmd_pub2, "Cougar", "Mule Deer", "Hunting")
  occ_cougelk_hunting_hunt <- occ_out(hs_cougelk_hunt0, "Cougar", "Elk", "Hunting")
  occ_cougelk_hunting_pub <- occ_out(hs_cougelk_pub0, "Cougar", "Elk", "Hunting")
  occ_cougwtd_hunting_hunt <- occ_out(hs_cougwtd_hunt1, "Cougar", "White-tailed Deer", "Hunting")
  occ_cougwtd_hunting_pub <- occ_out(hs_cougwtd_pub2, "Cougar", "White-tailed Deer", "Hunting")
  occ_cougmoose_hunting_hunt <- occ_out(hs_cougmoose_hunt0, "Cougar", "Moose", "Hunting")
  occ_cougmoose_hunting_pub <- occ_out(hs_cougmoose_pub0, "Cougar", "Moose", "Hunting")
  occ_wolfmd_hunting_hunt <- occ_out(hs_wolfmd_hunt2, "Wolf", "Mule Deer", "Hunting")
  occ_wolfmd_hunting_pub <- occ_out(hs_wolfmd_pub1, "Wolf", "Mule Deer", "Hunting")
  occ_wolfelk_hunting_hunt <- occ_out(hs_wolfelk_hunt2, "Wolf", "Elk", "Hunting")
  occ_wolfelk_hunting_pub <- occ_out(hs_wolfelk_pub1, "Wolf", "Elk", "Hunting")
  occ_wolfwtd_hunting_hunt <- occ_out(hs_wolfwtd_hunt2, "Wolf", "White-tailed Deer", "Hunting")
  occ_wolfwtd_hunting_pub <- occ_out(hs_wolfwtd_pub1, "Wolf", "White-tailed Deer", "Hunting")
  occ_wolfmoose_hunting_hunt <- occ_out(hs_wolfmoose_hunt2, "Wolf", "Moose", "Hunting")
  occ_wolfmoose_hunting_pub <- occ_out(hs_wolfmoose_pub1, "Wolf", "Moose", "Hunting")
  occ_bearmd_hunting_hunt <- occ_out(hs_bearmd_hunt2, "Black Bear", "Mule Deer", "Hunting")
  occ_bearmd_hunting_pub <- occ_out(hs_bearmd_pub2, "Black Bear", "Mule Deer", "Hunting")
  occ_bearelk_hunting_hunt <- occ_out(hs_bearelk_hunt2, "Black Bear", "Elk", "Hunting")
  occ_bearelk_hunting_pub <- occ_out(hs_bearelk_pub1, "Black Bear", "Elk", "Hunting")
  occ_bearwtd_hunting_hunt <- occ_out(hs_bearwtd_hunt2, "Black Bear", "White-tailed Deer", "Hunting")
  occ_bearwtd_hunting_pub <- occ_out(hs_bearwtd_pub2, "Black Bear", "White-tailed Deer", "Hunting")
  occ_bearmoose_hunting_hunt <- occ_out(hs_bearmoose_hunt2, "Black Bear", "Moose", "Hunting")
  occ_bearmoose_hunting_pub <- occ_out(hs_bearmoose_pub1, "Black Bear", "Moose", "Hunting")
  occ_coymd_hunting_hunt <- occ_out(hs_coymd_hunt2, "Coyote", "Mule Deer", "Hunting")
  occ_coymd_hunting_pub <- occ_out(hs_coymd_pub2, "Coyote", "Mule Deer", "Hunting")
  occ_coywtd_hunting_hunt <- occ_out(hs_coywtd_hunt1, "Coyote", "White-tailed Deer", "Hunting")
  occ_coywtd_hunting_pub <- occ_out(hs_coywtd_pub2, "Coyote", "White-tailed Deer", "Hunting")
  occ_bobmd_hunting_hunt <- occ_out(hs_bobmd_hunt2, "Bobcat", "Mule Deer Deer", "Hunting")
  occ_bobmd_hunting_pub <- occ_out(hs_bobmd_pub2, "Bobcat", "Mule Deer Deer", "Hunting")
  occ_bobwtd_hunting_hunt <- occ_out(hs_bobwtd_hunt1, "Bobcat", "White-tailed Deer", "Hunting")
  occ_bobwtd_hunting_pub <- occ_out(hs_bobwtd_pub2, "Bobcat", "White-tailed Deer", "Hunting")
  
  #'  Merge into larger data frames for easy comparison
  #'  Full models
  graze_allot_results <- rbind(occ_cougmd_grazing_allot, occ_cougelk_grazing_allot, occ_cougwtd_grazing_allot, occ_cougmoose_grazing_allot,
                               occ_wolfmd_grazing_allot, occ_wolfelk_grazing_allot, occ_wolfwtd_grazing_allot, occ_wolfmoose_grazing_allot,
                               occ_bearmd_grazing_allot, occ_bearelk_grazing_allot, occ_bearwtd_grazing_allot, occ_bearmoose_grazing_allot,
                               occ_coymd_grazing_allot, occ_coywtd_grazing_allot, occ_bobmd_grazing_allot, occ_bobwtd_grazing_allot)
  graze_cattle_results <- rbind(occ_cougmd_grazing_cattle, occ_cougelk_grazing_cattle, occ_cougwtd_grazing_cattle, occ_cougmoose_grazing_cattle,
                                occ_wolfmd_grazing_cattle, occ_wolfelk_grazing_cattle, occ_wolfwtd_grazing_cattle, occ_wolfmoose_grazing_cattle,
                                occ_bearmd_grazing_cattle, occ_bearelk_grazing_cattle, occ_bearwtd_grazing_cattle, occ_bearmoose_grazing_cattle,
                                occ_coymd_grazing_cattle, occ_coywtd_grazing_cattle, occ_bobmd_grazing_cattle, occ_bobwtd_grazing_cattle)
  hunt_hunt_results <- rbind(occ_cougmd_hunting_hunt, occ_cougelk_hunting_hunt, occ_cougwtd_hunting_hunt, occ_cougmoose_hunting_hunt,
                             occ_wolfmd_hunting_hunt, occ_wolfelk_hunting_hunt, occ_wolfwtd_hunting_hunt, occ_wolfmoose_hunting_hunt,
                             occ_bearmd_hunting_hunt, occ_bearelk_hunting_hunt, occ_bearwtd_hunting_hunt, occ_bearmoose_hunting_hunt,
                             occ_coymd_hunting_hunt, occ_coywtd_hunting_hunt, occ_bobmd_hunting_hunt, occ_bobwtd_hunting_hunt)
  hunt_pub_results <- rbind(occ_cougmd_hunting_pub, occ_cougelk_hunting_pub, occ_cougwtd_hunting_pub, occ_cougmoose_hunting_pub,
                            occ_wolfmd_hunting_pub, occ_wolfelk_hunting_pub, occ_wolfwtd_hunting_pub, occ_wolfmoose_hunting_pub,
                            occ_bearmd_hunting_pub, occ_bearelk_hunting_pub, occ_bearwtd_hunting_pub, occ_bearmoose_hunting_pub,
                            occ_coymd_hunting_pub, occ_coywtd_hunting_pub, occ_bobmd_hunting_pub, occ_bobwtd_hunting_pub)
  
  #'  Spread this out so the coefficient effects are easier to compare across species
  format_se <- function(out) {
    reformat_results <- out %>%  
      dplyr::select(-z) %>%
      mutate(SE = round(SE, 2),
             SE = paste0("(", SE, ")")) 
    return(reformat_results)
  }
  results_graze_allot <- format_se(graze_allot_results)
  results_graze_cattle <- format_se(graze_cattle_results)
  results_hunt_hunt <- format_se(hunt_hunt_results)
  results_hunt_pub <- format_se(hunt_pub_results)
  
  #'  Convert results to wide format
  format_wide <- function(out) {
    wide_results <- out %>%
      unite(Est_SE, Estimate, SE, sep = " ") %>%
      unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
      #'  Change species names to general classes
      mutate(
        Parameter = gsub("cougar", "Species 1", Parameter),
        Parameter = gsub("wolf", "Species 1", Parameter),
        Parameter = gsub("blackbear", "Species 1", Parameter),
        Parameter = gsub("coyote", "Species 1", Parameter),
        Parameter = gsub("bobcat", "Species 1", Parameter),
        Parameter = gsub("muledeer", "Species 2", Parameter),
        Parameter = gsub("mule_deer", "Species 2", Parameter),
        Parameter = gsub("elk", "Species 2", Parameter),
        Parameter = gsub("wtd", "Species 2", Parameter),
        Parameter = gsub("moose", "Species 2", Parameter)
      ) %>%
      spread(Parameter, Est_SE_Pval)
    return(wide_results)
  }
  wideresults_graze_allot <- format_wide(results_graze_allot)
  wideresults_graze_cattle <- format_wide(results_graze_cattle)
  wideresults_hunt_hunt <- format_wide(results_hunt_hunt)
  wideresults_hunt_pub <- format_wide(results_hunt_pub)
  
  #'  Grazing results tables
  results_graze_allot_wide <- wideresults_graze_allot %>%
    relocate("[Species 1:Species 2] (Intercept)", .after = "[Species 2] Study_AreaOK") %>%
    relocate("[Species 1:Species 2] PublicGrazing1", .after = "[Species 1:Species 2] (Intercept)") %>%
    relocate("[Species 1] I(Elev^2)", .after = "[Species 1] Elev") %>%
    relocate("[Species 2] I(Elev^2)", .after = "[Species 2] Elev") %>%
    relocate("[Species 1] PublicGrazing1", .after = "[Species 1] I(Elev^2)") %>%
    relocate("[Species 2] PublicGrazing1", .after = "[Species 2] I(Elev^2)") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] Elev", c("[Species 1] Elevation (SE)", "[Species 1] Elevation Pval"), sep = "_") %>%
    separate("[Species 2] Elev", c("[Species 2] Elevation (SE)", "[Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1] I(Elev^2)", c("[Species 1] Elevation^2 (SE)", "[Species 1] Elevation^2 Pval"), sep = "_") %>%
    separate("[Species 2] I(Elev^2)", c("[Species 2] Elevation^2 (SE)", "[Species 2] Elevation^2 Pval"), sep = "_") %>%
    separate("[Species 1] PercForest", c("[Species 1] PercentForest (SE)", "[Species 1] PercentForest Pval"), sep = "_") %>%
    separate("[Species 2] PercForest", c("[Species 2] PercentForest (SE)", "[Species 2] PercentForest Pval"), sep = "_") %>%
    separate("[Species 1] Study_AreaOK", c("[Species 1] Study_AreaOK (SE)", "[Species 1] Study_AreaOK Pval"), sep = "_") %>%
    separate("[Species 2] Study_AreaOK", c("[Species 2] Study_AreaOK (SE)", "[Species 2] Study_AreaOK Pval"), sep = "_") %>%
    separate("[Species 1] PublicGrazing1", c("[Species 1] GrazingAllotment (SE)", "[Species 1] GrazingAllotment Pval"), sep = "_") %>%
    separate("[Species 2] PublicGrazing1", c("[Species 2] GrazingAllotment (SE)", "[Species 2] GrazingAllotment Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] (Intercept)", c("[Species 1:Species 2] Intercept (SE)", "[Species 1:Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] PublicGrazing1", c("[Species 1:Species 2] GrazingAllotment (SE)", "[Species 1:Species 2] GrazingAllotment Pval"), sep = "_") %>%
    arrange(match(Species1, c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf"))) 
  
  results_graze_cattle_wide <- wideresults_graze_cattle %>%
    relocate("[Species 1:Species 2] (Intercept)", .after = "[Species 2] Study_AreaOK") %>%
    relocate("[Species 1:Species 2] GrazingActivity", .after = "[Species 1:Species 2] (Intercept)") %>%
    relocate("[Species 1] I(Elev^2)", .after = "[Species 1] Elev") %>%
    relocate("[Species 2] I(Elev^2)", .after = "[Species 2] Elev") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] Elev", c("[Species 1] Elevation (SE)", "[Species 1] Elevation Pval"), sep = "_") %>%
    separate("[Species 2] Elev", c("[Species 2] Elevation (SE)", "[Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1] I(Elev^2)", c("[Species 1] Elevation^2 (SE)", "[Species 1] Elevation^2 Pval"), sep = "_") %>%
    separate("[Species 2] I(Elev^2)", c("[Species 2] Elevation^2 (SE)", "[Species 2] Elevation^2 Pval"), sep = "_") %>%
    separate("[Species 1] PercForest", c("[Species 1] PercentForest (SE)", "[Species 1] PercentForest Pval"), sep = "_") %>%
    separate("[Species 2] PercForest", c("[Species 2] PercentForest (SE)", "[Species 2] PercentForest Pval"), sep = "_") %>%
    separate("[Species 1] Study_AreaOK", c("[Species 1] Study_AreaOK (SE)", "[Species 1] Study_AreaOK Pval"), sep = "_") %>%
    separate("[Species 2] Study_AreaOK", c("[Species 2] Study_AreaOK (SE)", "[Species 2] Study_AreaOK Pval"), sep = "_") %>%
    separate("[Species 1] GrazingActivity", c("[Species 1] GrazingActivity (SE)", "[Species 1] GrazingActivity Pval"), sep = "_") %>%
    separate("[Species 2] GrazingActivity", c("[Species 2] GrazingActivity (SE)", "[Species 2] GrazingActivity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] (Intercept)", c("[Species 1:Species 2] Intercept (SE)", "[Species 1:Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] GrazingActivity", c("[Species 1:Species 2] GrazingActivity (SE)", "[Species 1:Species 2] GrazingActivity Pval"), sep = "_") %>%
    arrange(match(Species1, c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf"))) 
  
  #'  Hunting results tables
  results_hunt_hunt_wide <- wideresults_hunt_hunt %>%
    relocate("[Species 1:Species 2] (Intercept)", .after = "[Species 2] Study_AreaOK") %>%
    relocate("[Species 1:Species 2] HuntingActivity", .after = "[Species 1:Species 2] (Intercept)") %>%
    relocate("[Species 1] I(Elev^2)", .after = "[Species 1] Elev") %>%
    relocate("[Species 2] I(Elev^2)", .after = "[Species 2] Elev") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] Elev", c("[Species 1] Elevation (SE)", "[Species 1] Elevation Pval"), sep = "_") %>%
    separate("[Species 2] Elev", c("[Species 2] Elevation (SE)", "[Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1] I(Elev^2)", c("[Species 1] Elevation^2 (SE)", "[Species 1] Elevation^2 Pval"), sep = "_") %>%
    separate("[Species 2] I(Elev^2)", c("[Species 2] Elevation^2 (SE)", "[Species 2] Elevation^2 Pval"), sep = "_") %>%
    separate("[Species 1] PercForest", c("[Species 1] PercentForest (SE)", "[Species 1] PercentForest Pval"), sep = "_") %>%
    separate("[Species 2] PercForest", c("[Species 2] PercentForest (SE)", "[Species 2] PercentForest Pval"), sep = "_") %>%
    separate("[Species 1] Study_AreaOK", c("[Species 1] Study_AreaOK (SE)", "[Species 1] Study_AreaOK Pval"), sep = "_") %>%
    separate("[Species 2] Study_AreaOK", c("[Species 2] Study_AreaOK (SE)", "[Species 2] Study_AreaOK Pval"), sep = "_") %>%
    separate("[Species 1] HuntingActivity", c("[Species 1] HuntingActivity (SE)", "[Species 1] HuntingActivity Pval"), sep = "_") %>%
    separate("[Species 2] HuntingActivity", c("[Species 2] HuntingActivity (SE)", "[Species 2] HuntingActivity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] (Intercept)", c("[Species 1:Species 2] Intercept (SE)", "[Species 1:Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] HuntingActivity", c("[Species 1:Species 2] HuntingActivity (SE)", "[Species 1:Species 2] HuntingActivity Pval"), sep = "_") %>%
    arrange(match(Species1, c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf"))) 
  
  results_hunt_pub_wide <- wideresults_hunt_pub %>%
    relocate("[Species 1:Species 2] (Intercept)", .after = "[Species 2] Study_AreaOK") %>%
    relocate("[Species 1:Species 2] Public1", .after = "[Species 1:Species 2] (Intercept)") %>%
    relocate("[Species 1] I(Elev^2)", .after = "[Species 1] Elev") %>%
    relocate("[Species 2] I(Elev^2)", .after = "[Species 2] Elev") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] Elev", c("[Species 1] Elevation (SE)", "[Species 1] Elevation Pval"), sep = "_") %>%
    separate("[Species 2] Elev", c("[Species 2] Elevation (SE)", "[Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1] I(Elev^2)", c("[Species 1] Elevation^2 (SE)", "[Species 1] Elevation^2 Pval"), sep = "_") %>%
    separate("[Species 2] I(Elev^2)", c("[Species 2] Elevation^2 (SE)", "[Species 2] Elevation^2 Pval"), sep = "_") %>%
    separate("[Species 1] PercForest", c("[Species 1] PercentForest (SE)", "[Species 1] PercentForest Pval"), sep = "_") %>%
    separate("[Species 2] PercForest", c("[Species 2] PercentForest (SE)", "[Species 2] PercentForest Pval"), sep = "_") %>%
    separate("[Species 1] Study_AreaOK", c("[Species 1] Study_AreaOK (SE)", "[Species 1] Study_AreaOK Pval"), sep = "_") %>%
    separate("[Species 2] Study_AreaOK", c("[Species 2] Study_AreaOK (SE)", "[Species 2] Study_AreaOK Pval"), sep = "_") %>%
    separate("[Species 1] Public1", c("[Species 1] PublicLand (SE)", "[Species 1] PublicLand Pval"), sep = "_") %>%
    separate("[Species 2] Public1", c("[Species 2] PublicLand (SE)", "[Species 2] PublicLand Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] (Intercept)", c("[Species 1:Species 2] Intercept (SE)", "[Species 1:Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Public1", c("[Species 1:Species 2] PublicLand (SE)", "[Species 1:Species 2] PublicLand Pval"), sep = "_") %>%
    arrange(match(Species1, c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf"))) 
  
  #'  Save!
  write.csv(results_graze_allot, paste0("./Tables/CoOcc_OccProb_GrazingResults_Allotment_", Sys.Date(), ".csv"))
  write.csv(results_graze_allot_wide, paste0("./Tables/CoOcc_OccProb_GrazingResults_Allotment_wide_", Sys.Date(), ".csv"))
  write.csv(results_graze_cattle, paste0("./Tables/CoOcc_OccProb_GrazingResults_CattleAct_", Sys.Date(), ".csv"))
  write.csv(results_graze_cattle_wide, paste0("./Tables/CoOcc_OccProb_GrazingResults_CattleAct_wide_", Sys.Date(), ".csv"))
  write.csv(results_hunt_hunt, paste0("./Tables/CoOcc_OccProb_HuntingResults_HunterAct_", Sys.Date(), ".csv"))
  write.csv(results_hunt_hunt_wide, paste0("./Tables/CoOcc_OccProb_HuntingResults_HunterAct_wide_", Sys.Date(), ".csv"))
  write.csv(results_hunt_pub, paste0("./Tables/CoOcc_OccProb_HuntingResults_Public_", Sys.Date(), ".csv"))
  write.csv(results_hunt_pub_wide, paste0("./Tables/CoOcc_OccProb_HuntingResults_Public_wide_", Sys.Date(), ".csv"))
  
  
  #'  Function to save detection results
  det_out <- function(mod, spp1, spp2, season) {
    out <- summary(mod@estimates)$det %>%
      mutate(
        Parameter = row.names(summary(mod@estimates)$det),
        Species1 = rep(spp1, nrow(.)),
        Species2 = rep(spp2, nrow(.)),
        Season = rep(season, nrow(.)),
        Estimate = round(Estimate, 2),
        SE = round(SE, 2),
        z = round(z, 2),
        Pval = round(`P(>|z|)`, 2)
      ) %>%
      dplyr::select(-`P(>|z|)`) %>%
      relocate(Parameter, .before = Estimate) %>%
      relocate(Species1, .before = Parameter) %>%
      relocate(Species2, .before = Parameter) %>%
      relocate(Season, .before = Parameter) 
    return(out)
  }
  #####  Grazing season models - p  #####
  det_cougmd_grazing_allot <- det_out(gs_cougmd_allot2, "Cougar", "Mule Deer", "Grazing")
  det_cougmd_grazing_cattle <- det_out(gs_cougmd_graze2, "Cougar", "Mule Deer", "Grazing")
  det_cougelk_grazing_allot <- det_out(gs_cougelk_allot0, "Cougar", "Elk", "Grazing")
  det_cougelk_grazing_cattle <- det_out(gs_cougelk_graze2, "Cougar", "Elk", "Grazing")
  det_cougwtd_grazing_allot <- det_out(gs_cougwtd_allot2, "Cougar", "White-tailed Deer", "Grazing")
  det_cougwtd_grazing_cattle <- det_out(gs_cougwtd_graze2, "Cougar", "White-tailed Deer", "Grazing")
  det_cougmoose_grazing_allot <- det_out(gs_cougmoose_allot2, "Cougar", "Moose", "Grazing")
  det_cougmoose_grazing_cattle <- det_out(gs_cougmoose_graze2, "Cougar", "Moose", "Grazing")
  det_wolfmd_grazing_allot <- det_out(gs_wolfmd_allot2, "Wolf", "Mule Deer", "Grazing")
  det_wolfmd_grazing_cattle <- det_out(gs_wolfmd_graze2, "Wolf", "Mule Deer", "Grazing")
  det_wolfelk_grazing_allot <- det_out(gs_wolfelk_allot1, "Wolf", "Elk", "Grazing")
  det_wolfelk_grazing_cattle <- det_out(gs_wolfelk_graze1, "Wolf", "Elk", "Grazing")
  det_wolfwtd_grazing_allot <- det_out(gs_wolfwtd_allot2, "Wolf", "White-tailed Deer", "Grazing")
  det_wolfwtd_grazing_cattle <- det_out(gs_wolfwtd_graze2, "Wolf", "White-tailed Deer", "Grazing")
  det_wolfmoose_grazing_allot <- det_out(gs_wolfmoose_allot2, "Wolf", "Moose", "Grazing")
  det_wolfmoose_grazing_cattle <- det_out(gs_wolfmoose_graze2, "Wolf", "Moose", "Grazing")
  det_bearmd_grazing_allot <- det_out(gs_bearmd_allot2, "Black Bear", "Mule Deer", "Grazing")
  det_bearmd_grazing_cattle <- det_out(gs_bearmd_graze2, "Black Bear", "Mule Deer", "Grazing")
  det_bearelk_grazing_allot <- det_out(gs_bearelk_allot2, "Black Bear", "Elk", "Grazing")
  det_bearelk_grazing_cattle <- det_out(gs_bearelk_graze1, "Black Bear", "Elk", "Grazing")
  det_bearwtd_grazing_allot <- det_out(gs_bearwtd_allot2, "Black Bear", "White-tailed Deer", "Grazing")
  det_bearwtd_grazing_cattle <- det_out(gs_bearwtd_graze2, "Black Bear", "White-tailed Deer", "Grazing")
  det_bearmoose_grazing_allot <- det_out(gs_bearmoose_allot2, "Black Bear", "Moose", "Grazing")
  det_bearmoose_grazing_cattle <- det_out(gs_bearmoose_graze1, "Black Bear", "Moose", "Grazing")
  det_coymd_grazing_allot <- det_out(gs_coymd_allot2, "Coyote", "Mule Deer", "Grazing")
  det_coymd_grazing_cattle <- det_out(gs_coymd_graze2, "Coyote", "Mule Deer", "Grazing")
  det_coywtd_grazing_allot <- det_out(gs_coywtd_allot2, "Coyote", "White-tailed Deer", "Grazing")
  det_coywtd_grazing_cattle <- det_out(gs_coywtd_graze2, "Coyote", "White-tailed Deer", "Grazing")
  det_bobmd_grazing_allot <- det_out(gs_bobmd_allot2, "Bobcat", "Mule Deer Deer", "Grazing")
  det_bobmd_grazing_cattle <- det_out(gs_bobmd_graze2, "Bobcat", "Mule Deer Deer", "Grazing")
  det_bobwtd_grazing_allot <- det_out(gs_bobwtd_allot2, "Bobcat", "White-tailed Deer", "Grazing")
  det_bobwtd_grazing_cattle <- det_out(gs_bobwtd_graze2, "Bobcat", "White-tailed Deer", "Grazing")
  #####  Hunting season models - p  #####
  det_cougmd_hunting_hunt <- det_out(hs_cougmd_hunt2, "Cougar", "Mule Deer", "Hunting")
  det_cougmd_hunting_pub <- det_out(hs_cougmd_pub2, "Cougar", "Mule Deer", "Hunting")
  det_cougelk_hunting_hunt <- det_out(hs_cougelk_hunt0, "Cougar", "Elk", "Hunting")
  det_cougelk_hunting_pub <- det_out(hs_cougelk_pub0, "Cougar", "Elk", "Hunting")
  det_cougwtd_hunting_hunt <- det_out(hs_cougwtd_hunt1, "Cougar", "White-tailed Deer", "Hunting")
  det_cougwtd_hunting_pub <- det_out(hs_cougwtd_pub2, "Cougar", "White-tailed Deer", "Hunting")
  det_cougmoose_hunting_hunt <- det_out(hs_cougmoose_hunt0, "Cougar", "Moose", "Hunting")
  det_cougmoose_hunting_pub <- det_out(hs_cougmoose_pub0, "Cougar", "Moose", "Hunting")
  det_wolfmd_hunting_hunt <- det_out(hs_wolfmd_hunt2, "Wolf", "Mule Deer", "Hunting")
  det_wolfmd_hunting_pub <- det_out(hs_wolfmd_pub1, "Wolf", "Mule Deer", "Hunting")
  det_wolfelk_hunting_hunt <- det_out(hs_wolfelk_hunt2, "Wolf", "Elk", "Hunting")
  det_wolfelk_hunting_pub <- det_out(hs_wolfelk_pub1, "Wolf", "Elk", "Hunting")
  det_wolfwtd_hunting_hunt <- det_out(hs_wolfwtd_hunt2, "Wolf", "White-tailed Deer", "Hunting")
  det_wolfwtd_hunting_pub <- det_out(hs_wolfwtd_pub1, "Wolf", "White-tailed Deer", "Hunting")
  det_wolfmoose_hunting_hunt <- det_out(hs_wolfmoose_hunt2, "Wolf", "Moose", "Hunting")
  det_wolfmoose_hunting_pub <- det_out(hs_wolfmoose_pub1, "Wolf", "Moose", "Hunting")
  det_bearmd_hunting_hunt <- det_out(hs_bearmd_hunt2, "Black Bear", "Mule Deer", "Hunting")
  det_bearmd_hunting_pub <- det_out(hs_bearmd_pub2, "Black Bear", "Mule Deer", "Hunting")
  det_bearelk_hunting_hunt <- det_out(hs_bearelk_hunt2, "Black Bear", "Elk", "Hunting")
  det_bearelk_hunting_pub <- det_out(hs_bearelk_pub1, "Black Bear", "Elk", "Hunting")
  det_bearwtd_hunting_hunt <- det_out(hs_bearwtd_hunt2, "Black Bear", "White-tailed Deer", "Hunting")
  det_bearwtd_hunting_pub <- det_out(hs_bearwtd_pub2, "Black Bear", "White-tailed Deer", "Hunting")
  det_bearmoose_hunting_hunt <- det_out(hs_bearmoose_hunt2, "Black Bear", "Moose", "Hunting")
  det_bearmoose_hunting_pub <- det_out(hs_bearmoose_pub1, "Black Bear", "Moose", "Hunting")
  det_coymd_hunting_hunt <- det_out(hs_coymd_hunt2, "Coyote", "Mule Deer", "Hunting")
  det_coymd_hunting_pub <- det_out(hs_coymd_pub2, "Coyote", "Mule Deer", "Hunting")
  det_coywtd_hunting_hunt <- det_out(hs_coywtd_hunt1, "Coyote", "White-tailed Deer", "Hunting")
  det_coywtd_hunting_pub <- det_out(hs_coywtd_pub2, "Coyote", "White-tailed Deer", "Hunting")
  det_bobmd_hunting_hunt <- det_out(hs_bobmd_hunt2, "Bobcat", "Mule Deer Deer", "Hunting")
  det_bobmd_hunting_pub <- det_out(hs_bobmd_pub2, "Bobcat", "Mule Deer Deer", "Hunting")
  det_bobwtd_hunting_hunt <- det_out(hs_bobwtd_hunt1, "Bobcat", "White-tailed Deer", "Hunting")
  det_bobwtd_hunting_pub <- det_out(hs_bobwtd_pub2, "Bobcat", "White-tailed Deer", "Hunting")
  
  #'  Merge into larger data frames for easy comparison
  #'  Full models
  graze_det_allot_results <- rbind(det_cougmd_grazing_allot, det_cougelk_grazing_allot, det_cougwtd_grazing_allot, det_cougmoose_grazing_allot,
                                   det_wolfmd_grazing_allot, det_wolfelk_grazing_allot, det_wolfwtd_grazing_allot, det_wolfmoose_grazing_allot,
                                   det_bearmd_grazing_allot, det_bearelk_grazing_allot, det_bearwtd_grazing_allot, det_bearmoose_grazing_allot,
                                   det_coymd_grazing_allot, det_coywtd_grazing_allot, det_bobmd_grazing_allot, det_bobwtd_grazing_allot)
  graze_det_cattle_results <- rbind(det_cougmd_grazing_cattle, det_cougelk_grazing_cattle, det_cougwtd_grazing_cattle, det_cougmoose_grazing_cattle,
                                    det_wolfmd_grazing_cattle, det_wolfelk_grazing_cattle, det_wolfwtd_grazing_cattle, det_wolfmoose_grazing_cattle,
                                    det_bearmd_grazing_cattle, det_bearelk_grazing_cattle, det_bearwtd_grazing_cattle, det_bearmoose_grazing_cattle,
                                    det_coymd_grazing_cattle, det_coywtd_grazing_cattle, det_bobmd_grazing_cattle, det_bobwtd_grazing_cattle)
  hunt_det_hunt_results <- rbind(det_cougmd_hunting_hunt, det_cougelk_hunting_hunt, det_cougwtd_hunting_hunt, det_cougmoose_hunting_hunt,
                                 det_wolfmd_hunting_hunt, det_wolfelk_hunting_hunt, det_wolfwtd_hunting_hunt, det_wolfmoose_hunting_hunt,
                                 det_bearmd_hunting_hunt, det_bearelk_hunting_hunt, det_bearwtd_hunting_hunt, det_bearmoose_hunting_hunt,
                                 det_coymd_hunting_hunt, det_coywtd_hunting_hunt, det_bobmd_hunting_hunt, det_bobwtd_hunting_hunt)
  hunt_det_pub_results <- rbind(det_cougmd_hunting_pub, det_cougelk_hunting_pub, det_cougwtd_hunting_pub, det_cougmoose_hunting_pub,
                                det_wolfmd_hunting_pub, det_wolfelk_hunting_pub, det_wolfwtd_hunting_pub, det_wolfmoose_hunting_pub,
                                det_bearmd_hunting_pub, det_bearelk_hunting_pub, det_bearwtd_hunting_pub, det_bearmoose_hunting_pub,
                                det_coymd_hunting_pub, det_coywtd_hunting_pub, det_bobmd_hunting_pub, det_bobwtd_hunting_pub)
  
  results_graze_det_allot <- format_se(graze_det_allot_results)
  results_graze_det_cattle <- format_se(graze_det_cattle_results)
  results_hunt_det_hunt <- format_se(hunt_det_hunt_results)
  results_hunt_det_pub <- format_se(hunt_det_pub_results)
  
  wideresults_graze_det_allot <- format_wide(results_graze_det_allot)
  wideresults_graze_det_cattle <- format_wide(results_graze_det_cattle)
  wideresults_hunt_det_hunt <- format_wide(results_hunt_det_hunt)
  wideresults_hunt_det_pub <- format_wide(results_hunt_det_pub)
  
  results_det_graze_allot_wide <- wideresults_graze_det_allot %>%
    relocate("[Species 1] PublicGrazing1", .after = "[Species 1] TrailDirt road") %>%
    relocate("[Species 2] PublicGrazing1", .after = "[Species 2] TrailDirt road") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] TrailDirt road", c("[Species 1] Dirt road (SE)", "[Species 1] Dirt road Pval"), sep = "_") %>%
    separate("[Species 2] TrailDirt road", c("[Species 2] Dirt road (SE)", "[Species 2] Dirt road Pval"), sep = "_") %>%
    separate("[Species 1] PublicGrazing1", c("[Species 1] GrazingAllotment (SE)", "[Species 1] GrazingAllotment Pval"), sep = "_") %>%
    separate("[Species 2] PublicGrazing1", c("[Species 2] GrazingAllotment (SE)", "[Species 2] GrazingAllotment Pval"), sep = "_") %>%
    arrange(match(Species1, c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf"))) 
  results_det_graze_cattle_wide <- wideresults_graze_det_cattle %>%
    relocate("[Species 1] WeeklyGrazing", .after = "[Species 1] TrailDirt road") %>%
    relocate("[Species 2] WeeklyGrazing", .after = "[Species 2] TrailDirt road") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] TrailDirt road", c("[Species 1] Dirt road (SE)", "[Species 1] Dirt road Pval"), sep = "_") %>%
    separate("[Species 2] TrailDirt road", c("[Species 2] Dirt road (SE)", "[Species 2] Dirt road Pval"), sep = "_") %>%
    separate("[Species 1] WeeklyGrazing", c("[Species 1] WeeklyGrazing (SE)", "[Species 1] WeeklyGrazing Pval"), sep = "_") %>%
    separate("[Species 2] WeeklyGrazing", c("[Species 2] WeeklyGrazing (SE)", "[Species 2] WeeklyGrazing Pval"), sep = "_") %>%
    arrange(match(Species1, c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf"))) 
  results_det_hunt_hunt_wide <- wideresults_hunt_det_hunt %>%
    relocate("[Species 1] WeeklyHunting", .after = "[Species 1] TrailDirt road") %>%
    relocate("[Species 2] WeeklyHunting", .after = "[Species 2] TrailDirt road") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] TrailDirt road", c("[Species 1] Dirt road (SE)", "[Species 1] Dirt road Pval"), sep = "_") %>%
    separate("[Species 2] TrailDirt road", c("[Species 2] Dirt road (SE)", "[Species 2] Dirt road Pval"), sep = "_") %>%
    separate("[Species 1] WeeklyHunting", c("[Species 1] WeeklyHunting (SE)", "[Species 1] WeeklyHunting Pval"), sep = "_") %>%
    separate("[Species 2] WeeklyHunting", c("[Species 2] WeeklyHunting (SE)", "[Species 2] WeeklyHunting Pval"), sep = "_") %>%
    arrange(match(Species1, c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf"))) 
  results_det_hunt_pub_wide <- wideresults_hunt_det_pub %>%
    relocate("[Species 1] Public1", .after = "[Species 1] TrailDirt road") %>%
    relocate("[Species 2] Public1", .after = "[Species 2] TrailDirt road") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] TrailDirt road", c("[Species 1] Dirt road (SE)", "[Species 1] Dirt road Pval"), sep = "_") %>%
    separate("[Species 2] TrailDirt road", c("[Species 2] Dirt road (SE)", "[Species 2] Dirt road Pval"), sep = "_") %>%
    separate("[Species 1] Public1", c("[Species 1] PublicLand (SE)", "[Species 1] PublicLand Pval"), sep = "_") %>%
    separate("[Species 2] Public1", c("[Species 2] PublicLand (SE)", "[Species 2] PublicLand Pval"), sep = "_") %>%
    arrange(match(Species1, c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf"))) 
  
  #'  Save!
  write.csv(results_graze_det_allot, paste0("./Tables/CoOcc_DetProb_GrazingResults_Allotment_", Sys.Date(), ".csv"))
  write.csv(results_det_graze_allot_wide, paste0("./Tables/CoOcc_DetProb_GrazingResults_Allotment_wide", Sys.Date(), ".csv"))
  write.csv(results_graze_det_cattle, paste0("./Tables/CoOcc_DetProb_GrazingResults_CattleAct_", Sys.Date(), ".csv"))
  write.csv(results_det_graze_cattle_wide, paste0("./Tables/CoOcc_DetProb_GrazingResults_CattleAct_wide", Sys.Date(), ".csv"))
  write.csv(results_hunt_det_hunt, paste0("./Tables/CoOcc_DetProb_HuntingResults_HunterAct_", Sys.Date(), ".csv"))
  write.csv(results_det_hunt_hunt_wide, paste0("./Tables/CoOcc_DetProb_HuntingResults_HunterAct_wide", Sys.Date(), ".csv"))
  write.csv(results_hunt_det_pub, paste0("./Tables/CoOcc_DetProb_HuntingResults_Public_", Sys.Date(), ".csv"))
  write.csv(results_det_hunt_pub_wide, paste0("./Tables/CoOcc_DetProb_HuntingResults_Public_wide", Sys.Date(), ".csv"))
  
  
  
  
