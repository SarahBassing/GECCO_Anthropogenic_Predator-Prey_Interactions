  #'  =============================================================
  #'  Code from Bassing et al. (2024) Global Ecology & Conservation
  #'  
  #'  Data formatting for Two-Species Co-Occurrence Models 
  #'  Washington Predator-Prey Project
  #'  =============================================================
  #'  Script to create unmarked data frames and format covariate data for multi-
  #'  species occupancy models including deer, elk, moose, black bears, cougars, 
  #'  wolves, coyotes, and bobcats plus cattle during the grazing season 2018-2020 
  #'  (July - Sept) and hunters during the hunting season 2018-2020 (Oct - Nov), 
  #'  respectively, with two-way interaction between predator and prey species. 
  #'  Grazing season co-occurrence models include 13 7-day sampling occasions 
  #'  comprising the peak of livestock activity detected on camera. 
  #'  Hunting season co-occurrence models include 8 7-day sampling occasions 
  #'  comprising the two general rifle hunting seasons in eastern Washington, USA. 
  #'  
  #'  Covariate data included in occupancy models were collected at each camera 
  #'  site or extracted from remotely sensed data.
  #'  
  #'  Must run detection history and cattle/hunter activity scripts first for 
  #'  this script to create unmarked data frames
  #'  --> Detection_histories_for_unmarked.R
  #'  --> Cattle_hunter_activity.R
  #'  =============================================================
  
  library(unmarked)
  library(tidyverse)
  
  #'  Read in covariate data
  #'  Land mgnt, habitat, roads, sum anthro activity, etc. => site-level occupancy covs
  #'  Dist. to focal pt, height, monitoring, weekly anthro activity = > survey/site-level detection covs
  stations <- read.csv("./Data/Camera_Location_Covariates.csv") %>% 
    #'  Get rid of mysterious space after one of the NEs (ugh)
    mutate(
      Study_Area = ifelse(Study_Area == "NE ", "NE", as.character(Study_Area)),
    ) %>%
    #'  Consolidate trail types
    mutate(
      Monitoring = ifelse(Monitoring == "Closed road", "Dirt road", as.character(Monitoring)),
      Monitoring = ifelse(Monitoring == "Decommissioned road", "Dirt road", as.character(Monitoring)),
      Monitoring = ifelse(Monitoring == "Game trail", "Trail", as.character(Monitoring)),
    ) %>%
    #'  Consolidate Land_Owner to binary private (0) vs public (1)
    #'  Note this allows private timberlands to be considered "public" for access
    #'  Use Land_Mgnt if want to keep private timberlands private
    mutate(
      Public = ifelse(Land_Owner == "Private", 0, 1)
    ) %>%
    #'  Identify grazing allotments/leases/permits on USFS, WA DNR, & WDFW lands
    mutate(USFS_grazing = ifelse(Allot_Active == "ACTIVE" & Cattle_Allot == "YES", 1, 0),
           USFS_grazing = ifelse(is.na(USFS_grazing), 0, USFS_grazing),
           DNR_grazing = ifelse(Study_Area == "OK" & DNR_parcel == 1 & Land_Owner == "WA DNR", 1, 0),
           #'  Had to do these by hand based on mapped grazing areas on DNR & WDFW websites
           DNR_grazing = ifelse(CameraLocation == "NE6685_85" | CameraLocation == "NE5884_39" | 
                                  CameraLocation == "NE4774_34", 1, DNR_grazing),
           WDFW_grazing = ifelse(CameraLocation == "OK4945_93" |CameraLocation == "OK4561_86" | 
                                   CameraLocation == "OK6316_62" | CameraLocation == "OK6317_73" | 
                                   CameraLocation == "OK4306_78" | CameraLocation == "OK4499_E" | 
                                   CameraLocation == "OK3733_14", 1, 0),
           PublicGrazing = ifelse(USFS_grazing == 1 | DNR_grazing == 1 | WDFW_grazing == 1, 1, 0)) %>%
    arrange(Year, CameraLocation) #NECESSARY TO MATCH DH's CAMERALOCATION ORDER
  
  #'  Remove rows when camera was completely inoperable so missing detection data
  #'  for ALL sampling occasions at that site
  #'  Grazing season missing all data at sites: 10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336
  #'  Hunting season missing all data at sites: 16, 23, 25, 27, 29, 38, 85, 111, 119, 128, 129, 144, 146, 156, 162, 216, 274, 282, 283
  stations_graze <- stations[-c(10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336),]
  stations_hunt <- stations[-c(16, 23, 25, 27, 29, 38, 85, 111, 119, 128, 129, 144, 146, 156, 162, 216, 274, 282, 283),]
  
  
  #'  Split site-level covariates up by study area
  stations_grazeNE <- filter(stations_graze, Study_Area == "NE")
  stations_grazeOK <- filter(stations_graze, Study_Area == "OK")
  stations_huntNE <- filter(stations_hunt, Study_Area == "NE")
  stations_huntOK <- filter(stations_hunt, Study_Area == "OK")
  
  #'  Scale and format site-level covariates
  format_covs <- function(cams) {
    #'  Rename, format, and scale as needed
    formatted <- transmute(cams,
                           CameraLocation = as.factor(CameraLocation),
                           Year = as.factor(Year),
                           Study_Area = as.factor(Study_Area),
                           Distance = scale(Distance_Focal_Point),
                           Height = scale(Height_frm_grnd),
                           Trail = as.factor(Monitoring),
                           Land_Mgnt = as.factor(Land_Mgnt),
                           Land_Owner = as.factor(Land_Owner),
                           Public = as.factor(Public),
                           PublicGrazing = as.factor(PublicGrazing),
                           GrazingActivity = scale(GrazingActivity), # Summed over 13 week period
                           HuntingActivity = scale(HuntingActivity), # Summed over 8 week period
                           PercForest = scale(PercForest),
                           Elev = scale(Elev)
    ) %>%
      arrange(Year, CameraLocation) #NECESSARY TO MATCH DH's CAMERALOCATION ORDER
    
    #'  Adjust reference category for Trail factors
    order_trail <- c("Trail", "Dirt road") 
    formatted <- formatted %>%
      mutate(
        Trail = fct_relevel(Trail, order_trail)
      )
    
    #'  Identify which sites have missing data for distance & height covariates
    noDist <- which(is.na(formatted$Distance)); print(noDist)
    noHgt <- which(is.na(formatted$Height)); print(noHgt)
    
    #'  Replace missing values with mean once covariates are z-transformed
    formatted$Distance[is.na(formatted$Distance),] <- 0
    formatted$Height[is.na(formatted$Height),] <- 0
    
    return(formatted)
  }
  stations_graze <- format_covs(stations_graze)
  stations_hunt <- format_covs(stations_hunt)
  stations_grazeNE <- format_covs(stations_grazeNE)
  stations_grazeOK <- format_covs(stations_grazeOK)
  stations_huntNE <- format_covs(stations_huntNE)
  stations_huntOK <- format_covs(stations_huntOK)
  
  #'  Make sure order of rows for DHs match order of rows for covariates
  head(DH_wtd_graze1820)
  head(DH_cattle_graze1820)
  head(stations_graze)
  
  DH_bear_hunt1820[115:120,] # should be end of Year1, start of Year2 cams
  DH_all_hunt1820[115:120,]
  stations_hunt[115:120,]
  
  head(DH_md_graze1820_OK)
  head(DH_cattle_graze1820_OK)
  head(stations_grazeOK)
  
  DH_moose_hunt1820_NE[52:57,] # should be end of Year1, start of Year2 NE cams
  DH_all_hunt1820_NE[52:57,]
  stations_huntNE[52:57,]
  
  
  #'  Check for correlation among covariates
  #'  Watch out for NAs (use = "complete.obs")
  site_covs <- stations_graze[ , c("Elev", "PercForest", "GrazingActivity", "HuntingActivity", "Distance", "Height")]
  (cov_matrix <- cor(site_covs, use = "complete.obs")) 
  site_covs_NE <- stations_huntNE[ , c("Elev", "PercForest", "GrazingActivity", "HuntingActivity", "Distance", "Height")]
  (cov_matrix_NE <- cor(site_covs_NE, use = "complete.obs"))
  site_covs_OK <- stations_huntOK[ , c("Elev", "PercForest","GrazingActivity", "HuntingActivity", "Distance", "Height")]
  (cov_matrix_OK <- cor(site_covs_OK, use = "complete.obs"))
  
  #'  Scale survey-level covariates
  scale_srvy_cov <- function(DH) {
    #'  Find mean & standard deviation of covariates across all sites & occasions
    mu <- mean(as.matrix(DH), na.rm = TRUE)
    sd <- sd(as.matrix(DH), na.rm = TRUE)
    
    #'  Z-transform (center observations around mean & scale by 1 SD)
    scaled <- ((DH - mu) / sd)
    
    return(scaled)
  }
  #'  Scale each survey-level covariate
  cattle_week_scaled <- scale_srvy_cov(DH_cattle_graze1820)
  cattle_week_scaled_NE <- scale_srvy_cov(DH_cattle_graze1820_NE)
  cattle_week_scaled_OK <- scale_srvy_cov(DH_cattle_graze1820_OK)
  
  hunter_week_scaled <- scale_srvy_cov(DH_all_hunt1820)
  hunter_week_scaled_NE <- scale_srvy_cov(DH_all_hunt1820_NE)
  hunter_week_scaled_OK <- scale_srvy_cov(DH_all_hunt1820_OK)
  
  Effort_graze_scaled <- scale_srvy_cov(Effort_graze1820)
  Effort_graze_scaled_NE <- scale_srvy_cov(Effort_graze1820_NE)
  Effort_graze_scaled_OK <- scale_srvy_cov(Effort_graze1820_OK)
  
  Effort_hunt_scaled <- scale_srvy_cov(Effort_hunt1820)
  Effort_hunt_scaled_NE <- scale_srvy_cov(Effort_hunt1820_NE)
  Effort_hunt_scaled_OK <- scale_srvy_cov(Effort_hunt1820_OK)
  
  #'  Create list of survey-level covariates for unmarked
  srvy_grazing_covs <- list(WeeklyGrazing = cattle_week_scaled,
                            SamplingEffort = Effort_graze_scaled)
  srvy_grazing_covs_NE <- list(WeeklyGrazing = cattle_week_scaled_NE,
                               SamplingEffort = Effort_graze_scaled_NE)
  srvy_grazing_covs_OK <- list(WeeklyGrazing = cattle_week_scaled_OK,
                               SamplingEffort = Effort_graze_scaled_OK)
  
  srvy_hunting_covs <- list(WeeklyHunting = hunter_week_scaled,
                            SamplingEffort = Effort_hunt_scaled)
  srvy_hunting_covs_NE <- list(WeeklyHunting = hunter_week_scaled_NE,
                               SamplingEffort = Effort_hunt_scaled_NE)
  srvy_hunting_covs_OK <- list(WeeklyHunting = hunter_week_scaled_OK,
                               SamplingEffort = Effort_hunt_scaled_OK)
  
  #'  Double check each list looks correct
  head(srvy_grazing_covs[[1]])
  head(srvy_grazing_covs[[2]])
  head(srvy_hunting_covs[[1]])
  head(srvy_hunting_covs[[2]])
  
  
  ####  Setup data for unmarked  ####
  #'  ---------------------------
  #'  Multi-species unmarkedDF --> unmarkedFrameOccuMulti (pg 151 of unmarked manual)
  #'  List relevant detection histories and create unmarked data frame for each 
  #'  multi-species occupancy model.
  #'  Define maximum interaction order with maxOrder. Defaults to all possible
  #'  interactions if not defined. 
  
  ####  COUGAR-MULE DEER UMFs  ####
  #'  Grazing Season
  coug_md_graze_DH <- list(cougar = DH_coug_graze1820, muledeer = DH_md_graze1820)
  coug_md_grazing_UMF <- unmarkedFrameOccuMulti(y = coug_md_graze_DH,
                                                siteCovs = stations_graze,
                                                obsCovs = srvy_grazing_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coug_md_grazing_UMF)
  #'  Review covariates
  summary(coug_md_grazing_UMF)
  #'  Look at natural parameters f design matrix
  coug_md_grazing_UMF@fDesign
  
  #'  Hunting Season
  coug_md_hunt_DH <- list(cougar = DH_coug_hunt1820, muledeer = DH_md_hunt1820)
  coug_md_hunting_UMF <- unmarkedFrameOccuMulti(y = coug_md_hunt_DH,
                                                siteCovs = stations_hunt,
                                                obsCovs = srvy_hunting_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coug_md_hunting_UMF)
  #'  Review covariates
  summary(coug_md_hunting_UMF)
  #'  Look at natural parameters f design matrix
  coug_md_hunting_UMF@fDesign
  
  ####  COUGAR-ELK UMFs  ####
  #'  Grazing Season
  coug_elk_graze_DH <- list(cougar = DH_coug_graze1820, elk = DH_elk_graze1820)
  coug_elk_grazing_UMF <- unmarkedFrameOccuMulti(y = coug_elk_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coug_elk_grazing_UMF)
  
  #'  Hunting Season
  coug_elk_hunt_DH <- list(cougar = DH_coug_hunt1820, elk = DH_elk_hunt1820)
  coug_elk_hunting_UMF <- unmarkedFrameOccuMulti(y = coug_elk_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coug_elk_hunting_UMF)
  
  ####  COUGAR-WHITE-TAILED DEER UMFs  ####
  #'  Grazing Season
  coug_wtd_graze_DH <- list(cougar = DH_coug_graze1820, wtd = DH_wtd_graze1820)
  coug_wtd_grazing_UMF <- unmarkedFrameOccuMulti(y = coug_wtd_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coug_wtd_grazing_UMF)
  
  #'  Hunting Season
  coug_wtd_hunt_DH <- list(cougar = DH_coug_hunt1820, wtd = DH_wtd_hunt1820)
  coug_wtd_hunting_UMF <- unmarkedFrameOccuMulti(y = coug_wtd_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coug_wtd_hunting_UMF)
  
  ####  COUGAR-MOOSE UMFs  ####
  #'  Grazing Season
  coug_moose_graze_DH <- list(cougar = DH_coug_graze1820, moose = DH_moose_graze1820)
  coug_moose_grazing_UMF <- unmarkedFrameOccuMulti(y = coug_moose_graze_DH,
                                                   siteCovs = stations_graze,
                                                   obsCovs = srvy_grazing_covs,
                                                   maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coug_moose_grazing_UMF)
  
  #'  Hunting Season
  coug_moose_hunt_DH <- list(cougar = DH_coug_hunt1820, moose = DH_moose_hunt1820)
  coug_moose_hunting_UMF <- unmarkedFrameOccuMulti(y = coug_moose_hunt_DH,
                                                   siteCovs = stations_hunt,
                                                   obsCovs = srvy_hunting_covs,
                                                   maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coug_moose_hunting_UMF)
  
  
  ####  WOLF-MULE DEER UMFs  ####
  #'  Grazing Season
  wolf_md_graze_DH <- list(wolf = DH_wolf_graze1820, muledeer = DH_md_graze1820)
  wolf_md_grazing_UMF <- unmarkedFrameOccuMulti(y = wolf_md_graze_DH,
                                                siteCovs = stations_graze,
                                                obsCovs = srvy_grazing_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(wolf_md_grazing_UMF)
  
  #'  Hunting Season
  wolf_md_hunt_DH <- list(wolf = DH_wolf_hunt1820, muledeer = DH_md_hunt1820)
  wolf_md_hunting_UMF <- unmarkedFrameOccuMulti(y = wolf_md_hunt_DH,
                                                siteCovs = stations_hunt,
                                                obsCovs = srvy_hunting_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(wolf_md_hunting_UMF)
  
  ####  WOLF-ELK UMFs  ####
  #'  Grazing Season
  wolf_elk_graze_DH <- list(wolf = DH_wolf_graze1820, elk = DH_elk_graze1820)
  wolf_elk_grazing_UMF <- unmarkedFrameOccuMulti(y = wolf_elk_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(wolf_elk_grazing_UMF)
  
  #'  Hunting Season
  wolf_elk_hunt_DH <- list(wolf = DH_wolf_hunt1820, elk = DH_elk_hunt1820)
  wolf_elk_hunting_UMF <- unmarkedFrameOccuMulti(y = wolf_elk_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(wolf_elk_hunting_UMF)
  
  ####  WOLF-WHITE-TAILED DEER UMFs  ####
  #'  Grazing Season
  wolf_wtd_graze_DH <- list(wolf = DH_wolf_graze1820, wtd = DH_wtd_graze1820)
  wolf_wtd_grazing_UMF <- unmarkedFrameOccuMulti(y = wolf_wtd_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(wolf_wtd_grazing_UMF)
  
  #'  Hunting Season
  wolf_wtd_hunt_DH <- list(wolf = DH_wolf_hunt1820, wtd = DH_wtd_hunt1820)
  wolf_wtd_hunting_UMF <- unmarkedFrameOccuMulti(y = wolf_wtd_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(wolf_wtd_hunting_UMF)
  
  ####  WOLF-MOOSE UMFs  ####
  #'  Grazing Season
  wolf_moose_graze_DH <- list(wolf = DH_wolf_graze1820, moose = DH_moose_graze1820)
  wolf_moose_grazing_UMF <- unmarkedFrameOccuMulti(y = wolf_moose_graze_DH,
                                                   siteCovs = stations_graze,
                                                   obsCovs = srvy_grazing_covs,
                                                   maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(wolf_moose_grazing_UMF)
  
  #'  Hunting Season
  wolf_moose_hunt_DH <- list(wolf = DH_wolf_hunt1820, moose = DH_moose_hunt1820)
  wolf_moose_hunting_UMF <- unmarkedFrameOccuMulti(y = wolf_moose_hunt_DH,
                                                   siteCovs = stations_hunt,
                                                   obsCovs = srvy_hunting_covs,
                                                   maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(wolf_moose_hunting_UMF)
  
  
  
  ####  BLACK Bear-MULE DEER UMFs  ####
  #'  Grazing Season
  bear_md_graze_DH <- list(blackbear = DH_bear_graze1820, muledeer = DH_md_graze1820)
  bear_md_grazing_UMF <- unmarkedFrameOccuMulti(y = bear_md_graze_DH,
                                                siteCovs = stations_graze,
                                                obsCovs = srvy_grazing_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bear_md_grazing_UMF)
  
  #'  Hunting Season
  bear_md_hunt_DH <- list(blackbear = DH_bear_hunt1820, muledeer = DH_md_hunt1820)
  bear_md_hunting_UMF <- unmarkedFrameOccuMulti(y = bear_md_hunt_DH,
                                                siteCovs = stations_hunt,
                                                obsCovs = srvy_hunting_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bear_md_hunting_UMF)
  
  ####  BLACK BEAR-ELK UMFs  ####
  #'  Grazing Season
  bear_elk_graze_DH <- list(blackbear = DH_bear_graze1820, elk = DH_elk_graze1820)
  bear_elk_grazing_UMF <- unmarkedFrameOccuMulti(y = bear_elk_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bear_elk_grazing_UMF)
  
  #'  Hunting Season
  bear_elk_hunt_DH <- list(blackbear = DH_bear_hunt1820, elk = DH_elk_hunt1820)
  bear_elk_hunting_UMF <- unmarkedFrameOccuMulti(y = bear_elk_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bear_elk_hunting_UMF)
  
  ####  BLACK BEAR-WHITE-TAILED DEER UMFs  ####
  #'  Grazing Season
  bear_wtd_graze_DH <- list(blackbear = DH_bear_graze1820, wtd = DH_wtd_graze1820)
  bear_wtd_grazing_UMF <- unmarkedFrameOccuMulti(y = bear_wtd_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bear_wtd_grazing_UMF)
  
  #'  Hunting Season
  bear_wtd_hunt_DH <- list(blackbear = DH_bear_hunt1820, wtd = DH_wtd_hunt1820)
  bear_wtd_hunting_UMF <- unmarkedFrameOccuMulti(y = bear_wtd_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bear_wtd_hunting_UMF)
  
  ####  BLACK BEAR-MOOSE UMFs  ####
  #'  Grazing Season
  bear_moose_graze_DH <- list(blackbear = DH_bear_graze1820, moose = DH_moose_graze1820)
  bear_moose_grazing_UMF <- unmarkedFrameOccuMulti(y = bear_moose_graze_DH,
                                                   siteCovs = stations_graze,
                                                   obsCovs = srvy_grazing_covs,
                                                   maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bear_moose_grazing_UMF)
  
  #'  Hunting Season
  bear_moose_hunt_DH <- list(blackbear = DH_bear_hunt1820, moose = DH_moose_hunt1820)
  bear_moose_hunting_UMF <- unmarkedFrameOccuMulti(y = bear_moose_hunt_DH,
                                                   siteCovs = stations_hunt,
                                                   obsCovs = srvy_hunting_covs,
                                                   maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bear_moose_hunting_UMF)
  
  
  
  ####  COYOTE-MULE DEER UMFs  ####
  #'  Grazing Season
  coy_md_graze_DH <- list(coyote = DH_coy_graze1820, mule_deer = DH_md_graze1820)
  coy_md_grazing_UMF <- unmarkedFrameOccuMulti(y = coy_md_graze_DH,
                                               siteCovs = stations_graze,
                                               obsCovs = srvy_grazing_covs,
                                               maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coy_md_grazing_UMF)
  
  #'  Hunting Season
  coy_md_hunt_DH <- list(coyote = DH_coy_hunt1820, mule_deer = DH_md_hunt1820)
  coy_md_hunting_UMF <- unmarkedFrameOccuMulti(y = coy_md_hunt_DH,
                                               siteCovs = stations_hunt,
                                               obsCovs = srvy_hunting_covs,
                                               maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coy_md_hunting_UMF)
  
  ####  COYOTE-WHITE-TAILED DEER UMFs  ####
  #'  Grazing Season
  coy_wtd_graze_DH <- list(coyote = DH_coy_graze1820, wtd = DH_wtd_graze1820)
  coy_wtd_grazing_UMF <- unmarkedFrameOccuMulti(y = coy_wtd_graze_DH,
                                                siteCovs = stations_graze,
                                                obsCovs = srvy_grazing_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coy_wtd_grazing_UMF)
  
  #'  Hunting Season
  coy_wtd_hunt_DH <- list(coyote = DH_coy_hunt1820, wtd = DH_wtd_hunt1820)
  coy_wtd_hunting_UMF <- unmarkedFrameOccuMulti(y = coy_wtd_hunt_DH,
                                                siteCovs = stations_hunt,
                                                obsCovs = srvy_hunting_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coy_wtd_hunting_UMF)
  
  
  
  ####  BOBCAT-MULE DEER UMFs  ####
  #'  Grazing Season
  bob_md_graze_DH <- list(bobcat = DH_bob_graze1820, mule_deer = DH_md_graze1820)
  bob_md_grazing_UMF <- unmarkedFrameOccuMulti(y = bob_md_graze_DH,
                                               siteCovs = stations_graze,
                                               obsCovs = srvy_grazing_covs,
                                               maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bob_md_grazing_UMF)
  
  #'  Hunting Season
  bob_md_hunt_DH <- list(bobcat = DH_bob_hunt1820, mule_deer = DH_md_hunt1820)
  bob_md_hunting_UMF <- unmarkedFrameOccuMulti(y = bob_md_hunt_DH,
                                               siteCovs = stations_hunt,
                                               obsCovs = srvy_hunting_covs,
                                               maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bob_md_hunting_UMF)
  
  ####  BOBCAT-WHITE-TAILED DEER UMFs  ####
  #'  Grazing Season
  bob_wtd_graze_DH <- list(bobcat = DH_bob_graze1820, wtd = DH_wtd_graze1820)
  bob_wtd_grazing_UMF <- unmarkedFrameOccuMulti(y = bob_wtd_graze_DH,
                                                siteCovs = stations_graze,
                                                obsCovs = srvy_grazing_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bob_wtd_grazing_UMF)
  
  #'  Hunting Season
  bob_wtd_hunt_DH <- list(bobcat = DH_bob_hunt1820, wtd = DH_wtd_hunt1820)
  bob_wtd_hunting_UMF <- unmarkedFrameOccuMulti(y = bob_wtd_hunt_DH,
                                                siteCovs = stations_hunt,
                                                obsCovs = srvy_hunting_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bob_wtd_hunting_UMF)
  
  
  #'  ====================================================
  #'  Use UMFs in co-occurrence models in Occupancy_Models.R script.
  #'  This runs multi-species occupancy models with 2-way species interactions.
  #'
  #'  END
  #'  ====================================================
