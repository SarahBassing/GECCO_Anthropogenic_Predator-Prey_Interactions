  #'  =============================================================
  #'  Code from Bassing et al. (2024) Global Ecology & Conservation
  #'  
  #'  Cattle & Hunter Detection Histories
  #'  Washington Predator-Prey Project
  #'  =============================================================
  #'  Script to create detection histories for cattle and hunters to be used as
  #'  survey and site-level covariates in occupancy models. Key difference is using
  #'  COUNTS of detected cattle and hunters instead of binary detection/non-detection
  #'  data to represent grazing and hunting activity at each camera site over time.
  #'  
  #'  Combines: 
  #'  "full_camdata_DATE.csv" from Detections_by_Camera_Station.R in 
  #'      WPPP_CameraTrapping.Rproj
  #'     -Contains ALL detections of animals, humans, & vehicles (no empties),
  #'      camera coordinates and deployment covariates (cam height, etc.)
  #'  "All_Camera_Stations_18-21_updated_DATE.csv"
  #'     -Contains camera locations (including updated names/locations when a
  #'     camera was moved), deployment & pull dates, and problem dates
  #'  =============================================================
  
  #'  Load libraries
  library(camtrapR)
  library(chron)
  library(lubridate)
  library(tidyverse)
  
  #'  Read in camera station operation dates
  cam_stations <- read_csv("./Data/All_Camera_Stations.csv")
  
  #'  Read in detection data and ensure date/tine are formatted correctly
  megadata <- read_csv("./Data/Raw_detections.csv") %>%
    mutate(
      DateTime = as.POSIXct(DateTime,
                            format="%Y-%m-%d %H:%M:%S",tz="America/Los_Angeles"),
      Date = as.Date(Date, format = "%Y-%m-%d"),
      Time = chron(times = Time)
    )
  
  #'  Split out into animal (wildlife & cattle) and human images
  animals <- filter(megadata, !is.na(Species)) %>%
    #'  Add column indicating kind of anthropogenic activity that is of interest
    mutate(
      Activity = ifelse(Species == "Cattle", "Grazing", "Other")) 
  
  humans <- filter(megadata, !is.na(HumanActivity)) %>%
    #'  Add column indicating kind of anthropogenic activity that is of interest
    mutate(
      Activity = ifelse(HumanActivity == "Hunter Bow" | HumanActivity == "Hunter Rifle", "Bow_Rifle_Hunting", "Other"))
  
  ####  Extract independent detections for cattle #### 
  #'  --------------------------------------------
  #'  Use 5 minute interval to identify detection events
  cow_dat <- filter(animals, Species == "Cattle") %>%
    arrange(CameraLocation, DateTime)
  
  caps <- c()
  caps[1] <- 1
  for (i in 2:nrow(cow_dat)){
    if (cow_dat$CameraLocation[i-1] != cow_dat$CameraLocation[i]) caps[i] = i
    else (if (cow_dat$Species[i-1] != cow_dat$Species[i]) caps[i] = i
          else (if (difftime(cow_dat$DateTime[i], cow_dat$DateTime[i-1], units = c("mins")) > 5) caps[i] = i
                else caps[i] = caps[i-1]))
  }
  
  caps <- as.factor(caps)
  
  #'  Add new column to larger data set
  capdata <- cbind(as.data.frame(cow_dat), caps)
  
  #'  Retain only the first image from each unique detection event
  cattle_detections <- capdata %>%
    group_by(caps) %>%
    slice(1L) %>%
    ungroup()
  
  ####  Extract independent detections for human activity  #### 
  #'  -----------------------------------------------------
  #'  Use 5 minute interval to identify detection events
  human_dat <- filter(humans, !is.na(HumanActivity)) %>%
    arrange(CameraLocation, DateTime)
  
  caps <- c()
  caps[1] <- 1
  for (i in 2:nrow(human_dat)){
    if (human_dat$CameraLocation[i-1] != human_dat$CameraLocation[i]) caps[i] = i
    else (if (human_dat$Species[i-1] != human_dat$Species[i]) caps[i] = i
          else (if (difftime(human_dat$DateTime[i], human_dat$DateTime[i-1], units = c("mins")) > 5) caps[i] = i
                else caps[i] = caps[i-1]))
  }
  
  caps <- as.factor(caps)
  
  #'  Add new column to larger data set
  capdata <- cbind(as.data.frame(human_dat), caps)
  
  #'  Retain only the first image from each unique detection event
  human_detections <- capdata %>%
    group_by(caps) %>%
    slice(1L) %>%
    ungroup()
  
  ####  Filter dates to specific range  ####
  #'  Cattle grazing activity 
  #'  Grazing Season 2018: 07/01/2018 - 09/29/2018 (thirteen 7-day sampling periods)
  cattle2018 <- cattle_detections %>%
    filter(Date > "2018-06-30") %>%
    filter(Date < "2018-09-30") %>%
    dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "Species", "Activity")  
  #'  Subset by study area
  NE_cattle18 <- filter(cattle2018, grepl("NE", CameraLocation))
  OK_cattle18 <- filter(cattle2018, grepl("OK", CameraLocation))
  #'  Grazing Season 2019: 07/01/2019 - 09/29/2019 (thirteen 7-day sampling periods)
  cattle2019 <- cattle_detections %>%
    filter(Date > "2019-06-30") %>%
    filter(Date < "2019-09-30") %>%
    dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "Species", "Activity")  
  #'  Subset by study area
  NE_cattle19 <- filter(cattle2019, grepl("NE", CameraLocation))
  OK_cattle19 <- filter(cattle2019, grepl("OK", CameraLocation))
  #'  Grazing Season 2020: 07/01/2020 - 09/29/2020 (thirteen 7-day sampling periods)
  cattle2020 <- cattle_detections %>%
    filter(Date > "2020-06-30") %>%
    filter(Date < "2020-09-30") %>%
    dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "Species", "Activity")  
  #'  Subset by study area
  NE_cattle20 <- filter(cattle2020, grepl("NE", CameraLocation))
  OK_cattle20 <- filter(cattle2020, grepl("OK", CameraLocation))
  
  #'  Hunter Activity
  #'  Hunting Season 2018: 10/1/2018 - 11/25/2018 (eight 7-day sampling periods)
  hunters2018 <- human_detections %>%
    filter(Date > "2018-09-30") %>%
    filter(Date < "2018-11-26") %>%
    dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "Species", "HumanActivity", "Activity") 
  #'  Subset by study area
  NE_hunters18 <- filter(hunters2018, grepl("NE", CameraLocation))
  OK_hunters18 <- filter(hunters2018, grepl("OK", CameraLocation))
  #'  Hunting Season 2019: 10/1/2019 - 11/25/2019 (eight 7-day sampling periods)
  hunters2019 <- human_detections %>%
    filter(Date > "2019-09-30") %>%
    filter(Date < "2019-11-26") %>%
    dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "Species", "HumanActivity", "Activity")  
  #'  Subset by study area
  NE_hunters19 <- filter(hunters2019, grepl("NE", CameraLocation))
  OK_hunters19 <- filter(hunters2019, grepl("OK", CameraLocation))
  #'  Hunting Season 2020: 10/1/2020 - 11/25/2020 (eight 7-day sampling periods)
  hunters2020 <- human_detections %>%
    filter(Date > "2020-09-30") %>%
    filter(Date < "2020-11-26") %>%
    dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "Species", "HumanActivity", "Activity")  
  #'  Subset by study area
  NE_hunters20 <- filter(hunters2020, grepl("NE", CameraLocation))
  OK_hunters20 <- filter(hunters2020, grepl("OK", CameraLocation))
  
  
  #'  Save for data summary and making data available for publication
  hunter_dets <- as.data.frame(rbind(hunters2018, hunters2019, hunters2020)) %>%
    dplyr::select(-c(HumanActivity))
  CowHunter_Detections <- as.data.frame(rbind(cattle2018, cattle2019, 
                                              cattle2020, hunter_dets)) %>%
    mutate(
      Season = ifelse(Date < "2018-09-30", "Grazing18", "Grazing19"),
      Season = ifelse(Date > "2020-01-01" & Date < "2020-09-30", "Grazing20", Season),
      Season = ifelse(Date > "2018-09-30" & Date < "2018-11-30", "Hunt18", Season),
      Season = ifelse(Date > "2019-09-30" & Date < "2019-11-30", "Hunt19", Season),
      Season = ifelse(Date > "2020-09-30", "Hunt20", Season)
    ) %>%
    filter(Activity == "Bow_Rifle_Hunting" | Activity == "Grazing")
  
  
  ####  Camera Operation Table  ####
  #'  ------------------------------ 
  #'  Creates a matrix with each camera & dates it deployed
  #'  Define date format at this step, not before!
  #'  1 = operating; 0 = not operating but deployed; NA = not deployed
  camop_problem <- cameraOperation(CTtable = cam_stations,
                                   stationCol = "CameraLocation",
                                   setupCol = "Setup_date",
                                   retrievalCol = "Retrieval_date",
                                   hasProblems = TRUE,
                                   dateFormat = "%m/%d/%Y", # Define date format here!
                                   writecsv = FALSE) 
  
  probs <- as.data.frame(camop_problem)
  
  #'  Filter data by study area
  NEcams <- filter(cam_stations, grepl("NE", CameraLocation))
  camop_problem_NE <- cameraOperation(CTtable = NEcams,
                                      stationCol = "CameraLocation",
                                      setupCol = "Setup_date",
                                      retrievalCol = "Retrieval_date",
                                      hasProblems = TRUE,
                                      dateFormat = "%m/%d/%Y", 
                                      writecsv = FALSE) 
  
  OKcams <- filter(cam_stations, grepl("OK", CameraLocation))
  camop_problem_OK <- cameraOperation(CTtable = OKcams,
                                      stationCol = "CameraLocation",
                                      setupCol = "Setup_date",
                                      retrievalCol = "Retrieval_date",
                                      hasProblems = TRUE,
                                      dateFormat = "%m/%d/%Y", 
                                      writecsv = FALSE) 
  
  ####  Detection Histories  ####
  #'  ---------------------------
  #'  Function to create season-specific detection histories for cattle/hunters
  #'  for each season of interest.
  #'  
  #'  Notes about arguments:
  #'  -occasionLength: 7 day sampling occasions
  #'  -day1: sampling occasion begins on a specific date (e.g., "2018-07-01")
  #'   FYI this defines start date but NOT end date so DH goes until camera pulls
  #'  -includeEffort: compute # active trap days/station/occasion
  #'  -scaleEffort: center & scale effort (FALSE)
  #'  -output: return binary detections or counts of detections; using "count"
  #'   b/c images are thinned to independent detection events
  #'   
  #'  FYI, cannot have any NAs in the Species column or this doesn't work
  #'  Need to remove columns that extend beyond date range of interest!
  #'  Grazing: July 1 - Sept 29 = 13 weeks
  #'  Hunting: Oct 1 - Nov 26 = 8 weeks 
  DH_counts <- function(images, spp, start_date) {
    det_hist <- detectionHistory(recordTable = images,
                                 camOp = camop_problem,
                                 stationCol = "CameraLocation",
                                 speciesCol = "Activity",
                                 recordDateTimeCol = "DateTime",
                                 recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                 species = spp,
                                 occasionLength = 7,
                                 day1 = start_date, 
                                 timeZone = "America/Los_Angeles",
                                 output = "count",
                                 includeEffort = TRUE,
                                 scaleEffort = FALSE,
                                 outDir = "./Data/Detection_Histories")
    
    return(det_hist)
  }
  #'  Generate detection histories for cattle and hunter activity
  #'  Note: grazing DH for grazing season vs hunter DHs for hunting season
  
  ####  CATTLE  ####
  #'  Survey-level covariate for grazing season only
  cattle_graze18 <- DH_counts(cattle2018, "Grazing", "2018-07-01")
  DH_cattle_graze18 <- cattle_graze18[[1]][1:125,1:13]   
  cattle_graze19 <- DH_counts(cattle2019, "Grazing", "2019-07-01")
  DH_cattle_graze19 <- cattle_graze19[[1]][126:242,1:13] 
  cattle_graze20 <- DH_counts(cattle2020, "Grazing", "2020-07-01")
  DH_cattle_graze20 <- cattle_graze20[[1]][243:361,1:13] 
  
  DH_cattle_graze1820 <- rbind(DH_cattle_graze18, DH_cattle_graze19, DH_cattle_graze20)
  #'  Remove rows missing detection data for ALL occasions (12 cameras completely inoperable)
  DH_cattle_graze1820 <- DH_cattle_graze1820[-c(10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336),]
  DH_cattle_graze1820_NE <- DH_cattle_graze1820[grepl("NE", row.names(DH_cattle_graze1820)),]
  DH_cattle_graze1820_OK <- DH_cattle_graze1820[grepl("OK", row.names(DH_cattle_graze1820)),]
  
  
  ####  ALL HUNTERS ON-FOOT  ####
  #'  Survey-level covariate for hunting season only
  all_hunt18 <- DH_counts(hunters2018, "Bow_Rifle_Hunting", "2018-10-01")
  DH_all_hunt18 <- all_hunt18[[1]][1:125,1:8]  
  all_hunt19 <- DH_counts(hunters2019, "Bow_Rifle_Hunting", "2019-10-01")
  DH_all_hunt19 <- all_hunt19[[1]][126:242,1:8] 
  all_hunt20 <- DH_counts(hunters2020, "Bow_Rifle_Hunting", "2020-10-01")
  DH_all_hunt20 <- all_hunt20[[1]][243:361,1:8]
  
  DH_all_hunt1820 <- rbind(DH_all_hunt18, DH_all_hunt19, DH_all_hunt20)
  #'  Remove rows missing detection data for ALL occasions (19 cameras completely inoperable)
  DH_all_hunt1820 <- DH_all_hunt1820[-c(16, 23, 25, 27, 29, 38, 85, 111, 119, 128, 129, 144, 146, 156, 162, 216, 274, 282, 283),]
  DH_all_hunt1820_NE <- DH_all_hunt1820[grepl("NE", row.names(DH_all_hunt1820)),]
  DH_all_hunt1820_OK <- DH_all_hunt1820[grepl("OK", row.names(DH_all_hunt1820)),]
  
 
  ####  SAMPLING EFFORT  ####
  #'  -------------------
  #'  Pull out sampling effort for annual grazing season and hunting season
  SEffort_graze18 <- cattle_graze18[[2]][1:125,1:13]
  SEffort_graze19 <- cattle_graze19[[2]][126:242,1:13]
  SEffort_graze20 <- cattle_graze20[[2]][243:361,1:13]
  SEffort_hunt18 <- all_hunt18[[2]][1:125,1:8]
  SEffort_hunt19 <- all_hunt19[[2]][126:242,1:8]
  SEffort_hunt20 <- all_hunt20[[2]][243:361,1:8]
  
  #'  Sum number of trap nights per camera
  SEffort_graze1820 <- rbind(SEffort_graze18, SEffort_graze19, SEffort_graze20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  SEffort_graze1820 <- SEffort_graze1820[-c(10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336),]
  TrpNgts_graze1820 <- as.data.frame(rowSums(SEffort_graze1820,  na.rm = TRUE))
  colnames(TrpNgts_graze1820) <- "Trap_Nights"
  TrpNgts_graze1820_NE <- TrpNgts_graze1820[grepl("NE", row.names(TrpNgts_graze1820)),]
  TrpNgts_graze1820_OK <- TrpNgts_graze1820[grepl("OK", row.names(TrpNgts_graze1820)),]
  
  SEffort_hunt1820 <- rbind(SEffort_hunt18, SEffort_hunt19, SEffort_hunt20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  SEffort_hunt1820 <- SEffort_hunt1820[-c(16, 23, 25, 27, 29, 38, 85, 111, 119, 128, 129, 144, 146, 156, 162, 216, 274, 282, 283),]
  TrpNgts_hunt1820 <- as.data.frame(rowSums(SEffort_hunt1820,  na.rm = TRUE))
  colnames(TrpNgts_hunt1820) <- "Trap_Nights"
  TrpNgts_hunt1820_NE <- TrpNgts_hunt1820[grepl("NE", row.names(TrpNgts_hunt1820)),]
  TrpNgts_hunt1820_OK <- TrpNgts_hunt1820[grepl("OK", row.names(TrpNgts_hunt1820)),]
  
  
  ####  Site-Level Anthropogenic Covariates  ####
  #'  ---------------------------------------
  #'  Generate site-level measure of anthropogenic activity for each camera,
  #'  standardized by number of trap nights per site 
  #'  DON'T remove rows with all missing observations yet- save for full covariate 
  #'  data frame in Occupancy_Models.R script
  
  #'  Total cattle detection events across entire study period
  sum_graze18 <- rowSums(DH_cattle_graze18, na.rm = TRUE)
  sum_graze19 <- rowSums(DH_cattle_graze19, na.rm = TRUE)
  sum_graze20 <- rowSums(DH_cattle_graze20, na.rm = TRUE)
  #'  Sum number of trap nights per camera (including inoperable cams)
  SEffort_graze <- rbind(SEffort_graze18, SEffort_graze19, SEffort_graze20)
  TrpNgts_graze <- as.data.frame(rowSums(SEffort_graze,  na.rm = TRUE))
  #'  Create data frame with counts of detections and sampling effort per site
  GrazingActivity <- as.data.frame(c(sum_graze18, sum_graze19, sum_graze20)) %>%
    cbind(row.names(.)) %>%
    cbind(TrpNgts_graze) 
  colnames(GrazingActivity) <- c("CattleCounts", "CameraLocation", "TrapNights")
  #'  Calculate mean detections/trap night so standardized by effort
  GrazingActivity <- GrazingActivity %>%
    transmute(CameraLocation = CameraLocation,
              GrazingActivity = CattleCounts/TrapNights) %>%
    #'  0 detections/0 trap nights causes problems
    mutate(GrazingActivity = ifelse(GrazingActivity == "NaN", 0, GrazingActivity),
           GrazingActivity = round(GrazingActivity, 2))
  
  #'  Total hunter detection events across entire study period
  sum_hunt18 <- rowSums(DH_all_hunt18, na.rm = TRUE)
  sum_hunt19 <- rowSums(DH_all_hunt19, na.rm = TRUE)
  sum_hunt20 <- rowSums(DH_all_hunt20, na.rm = TRUE)
  #'  Sum number of trap nights per camera (including inoperable cams)
  SEffort_hunt <- rbind(SEffort_hunt18, SEffort_hunt19, SEffort_hunt20)
  TrpNgts_hunt <- as.data.frame(rowSums(SEffort_hunt,  na.rm = TRUE))
  #'  Create data frame with counts of detections and sampling effort per site
  HuntingActivity <- as.data.frame(c(sum_hunt18, sum_hunt19, sum_hunt20)) %>%
    cbind(row.names(.)) %>%
    cbind(TrpNgts_hunt)
  colnames(HuntingActivity) <- c("HunterCounts", "CameraLocation", "TrapNights")
  #'  Calculate mean detections/trap night so standardized by effort
  HuntingActivity <- HuntingActivity %>%
    transmute(CameraLocation = CameraLocation,
              HuntingActivity = HunterCounts/TrapNights) %>%
    #'  0 detections/0 trap nights causes problems
    mutate(HuntingActivity = ifelse(HuntingActivity == "NaN", 0, HuntingActivity),
           HuntingActivity = round(HuntingActivity, 2))
  
  #'  Create single site-level df of anthropogenic activity to add to occupancy covs
  anthro_covs <- full_join(GrazingActivity, HuntingActivity, by = "CameraLocation") %>%
    relocate(CameraLocation, .before = "GrazingActivity")

  
  ####  Summary Stats  ####
  #'  -----------------
  #'  Number of independent detections per season from camera traps
  ndet <- CowHunter_Detections %>%
    group_by(Season, Activity) %>%
    summarise(n = n()) %>%
    ungroup()
  summary_dets <- group_by(ndet, Activity) %>% 
    summarize(mu_locs = mean(n), sd = sd(n), se_locs = sd(n)/sqrt(n())) %>% 
    ungroup()
  
  #'  Percent of cameras where a species was detected
  #'  Based on number of active cameras in each study area and season
  ncams_graze <- nrow(eff_graze1820)
  NEcams_graze <- nrow(eff_graze1820[grepl("NE", eff_graze1820$CameraLocation),])
  OKcams_graze <- nrow(eff_graze1820[grepl("OK", eff_graze1820$CameraLocation),])
  ncams_hunt <- nrow(eff_hunt1820)
  NEcams_hunt <- nrow(eff_hunt1820[grepl("NE", eff_hunt1820$CameraLocation),])
  OKcams_hunt <- nrow(eff_hunt1820[grepl("OK", eff_hunt1820$CameraLocation),])
  
  perc_cams <- CowHunter_Detections %>%
    dplyr::select(CameraLocation, Season, Activity) %>%
    group_by(Activity, CameraLocation) %>%
    filter(row_number(CameraLocation) == 1) %>%
    ungroup() %>%
    group_by(Activity) %>%
    summarise(ncams = n()) %>%
    ungroup() %>%
    mutate(
      propcams = ifelse(Activity == "Grazing", ncams/ncams_graze, ncams/ncams_hunt),
      propcams = round(propcams, 2)) %>%
    dplyr::select(-ncams)
  colnames(perc_cams) <- c("Activity", "Percent of cameras")
  
  graze_cams <- CowHunter_Detections %>%
    filter(Activity == "Grazing") %>%
    dplyr::select(CameraLocation, Season, Activity) %>%
    group_by(Activity, CameraLocation) %>%
    filter(row_number(CameraLocation) == 1) %>%
    ungroup() 
  nOKgraze_cams <- filter(graze_cams, grepl("OK", CameraLocation))
  (percOKgraze_cams <- nrow(nOKgraze_cams)/nrow(graze_cams))
  
  hunter_cams <- CowHunter_Detections %>%
    filter(Activity != "Grazing") %>%
    dplyr::select(CameraLocation, Season, Activity) %>%
    group_by(Activity, CameraLocation) %>%
    filter(row_number(CameraLocation) == 1) %>%
    ungroup() 
  nOKhunter_cams <- filter(hunter_cams, grepl("OK", CameraLocation))
  (percOKhunter_cams <- nrow(nOKhunter_cams)/nrow(hunter_cams))
  
  
  
  
