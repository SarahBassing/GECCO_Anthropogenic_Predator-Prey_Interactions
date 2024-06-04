  #'  =============================================================
  #'  Code from Bassing et al. (2024) Global Ecology & Conservation
  #'  
  #'  Temporal overlap in response to cattle/hunters
  #'  Washington Predator-Prey Project
  #'  =============================================================
  #'  Script to estimate diel activity patterns and temporal overlap between
  #'  predators and prey during the cattle grazing and hunting seasons, as 
  #'  well as temporal overlap for each individual species at sites with and
  #'  without cattle/hunter activity. Then generates result tables and various
  #'  figures comparing coefficient of overlap and activity curves.
  #'  
  #'  NOTE: Need to create a folder within working directory called "Figures" to 
  #'  save temporal overlap plots.
  #'  =============================================================
  
  #'  Load libraries
  library(lubridate)
  library(chron)
  library(overlap)
  library(circular)
  library(ggplot2)
  library(khroma)
  library(patchwork)
  library(sp)
  library(tidyverse)

  #'  Read in formatted detection and covariate data with sun time added in 
  #'  Sun Time: accounts for seasonal changes in sunlight by converting clock time 
  #'  to "sun time" by mapping sunrise to π/2 and sunset to 3π/2. Sunrise & sunset 
  #'  times are determined based on the date and camera locations 
  megadata <- read_csv("./Data/Detection_data_sunTime.csv")
  
  ####  Extract independent detections for wildlife species  ####
  #'  -------------------------------------------------------
  dat <- arrange(megadata, CameraLocation, DateTime)
  caps <- c()
  caps[1] <- 1
  for (i in 2:nrow(dat)){
    if (dat$CameraLocation[i-1] != dat$CameraLocation[i]) caps[i] = i
    else (if (dat$Species[i-1] != dat$Species[i]) caps[i] = i
          else (if (difftime(dat$DateTime[i], dat$DateTime[i-1], units = c("mins")) > 30) caps[i] = i
                else caps[i] = caps[i-1]))
  }
  
  caps <- as.factor(caps)
  
  #'  Add new column to larger data set
  capdata <- cbind(as.data.frame(dat), caps)
  
  #'  Retain only the first image from each unique detection event
  detections <- capdata %>%
    group_by(caps) %>%
    slice(1L) %>%
    ungroup()
  detections_OK <- detections %>%
    filter(grepl("OK", CameraLocation))
  
  #'  Retain only the last image from each unique detection event
  last_det <- capdata %>%
    group_by(caps) %>%
    slice_tail(n = 1) %>%
    ungroup()
  last_det_OK <- last_det %>%
    filter(grepl("OK", CameraLocation))
  
  #'  Filter data to desired date ranges
  grazing_filter <- function(dets) {
    cattle2018 <- dets %>%
      filter(Date > "2018-06-30") %>%
      filter(Date < "2018-09-30") %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", 
                    "Species", "HumanActivity", "PublicGrazing")  
    cattle2019 <- dets %>%
      filter(Date > "2019-06-30") %>%
      filter(Date < "2019-09-30") %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", 
                    "Species", "HumanActivity", "PublicGrazing")  
    cattle2020 <- dets %>%
      filter(Date > "2020-06-30") %>%
      filter(Date < "2020-09-30") %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", 
                    "Species", "HumanActivity", "PublicGrazing")  
    detections <- rbind(cattle2018, cattle2019, cattle2020)
    return(detections)
  }
  grazing_first <- grazing_filter(detections)
  grazing_last <- grazing_filter(last_det)
  grazing_first_OK <- grazing_filter(detections_OK)
  grazing_last_OK <- grazing_filter(last_det_OK)
  
  #'  Hunting season
  hunting_filter <- function(dets) {
    hunters2018 <- dets %>%
      filter(Date > "2018-09-30") %>%
      filter(Date < "2018-11-26") %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", 
                    "Species", "HumanActivity", "Public1") 
    hunters2019 <- dets %>%
      filter(Date > "2019-09-30") %>%
      filter(Date < "2019-11-26") %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", 
                    "Species", "HumanActivity", "Public1")  
    hunters2020 <- dets %>%
      filter(Date > "2020-09-30") %>%
      filter(Date < "2020-11-26") %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", 
                    "Species", "HumanActivity", "Public1")  
    last_dets <- rbind(hunters2018, hunters2019, hunters2020)
    return(last_dets)
  }
  hunting_first <- hunting_filter(detections)
  hunting_last <- hunting_filter(last_det)
  
  
  #'  ------------------------------------------
  ####  Predator-prey temporal overlap analysis ####
  #'  ------------------------------------------
  #'  Function to estimate temporal overlap between predators (spp1) and prey (spp2)
  #'  at camera sites where cattle/humans (spp3) are present (detected) vs absent (not detected).
  pred_prey_overlap <- function(spp1, spp2, spp3, name1, name2, name3, nboot, dhat) { 
    
    #'  Create logical vectors (T/F) indicating whether spp1 & spp2 were detected 
    #'  at the same site and reduce detection events to just those cameras --> These 
    #'  species need to spatially overlap for any temporal overlap to be meaningful
    both.present <- spp1$CameraLocation %in% spp2$CameraLocation
    spp1_dat <- cbind(spp1, both.present) 
    spp1_dat <- spp1_dat[spp1_dat$both.present == T,]
    both.present <- spp2$CameraLocation %in% spp1$CameraLocation
    spp2_dat <- cbind(spp2, both.present)
    spp2_dat <- spp2_dat[spp2_dat$both.present == T,]
    #'  Double check the same number of cameras are included in each
    length(unique(spp1_dat$CameraLocation)); length(unique(spp2_dat$CameraLocation))
    
    #'  Create a logical vector (T/F) indicating whether spp3 was detected at each
    #'  site where spp1 was detected (based on cameras where spp1 & spp2 were detected)
    spp3.present <- spp1_dat$CameraLocation %in% spp3$CameraLocation
    spp1_dat <- cbind(spp1_dat, spp3.present)
    #'  Split out data into camera locations where both spp1 & spp3 are present vs spp3 absent
    spp1_spp3.present <- spp1_dat[spp1_dat$spp3.present == T,]
    spp1_spp3.absent <- spp1_dat[spp1_dat$spp3.present == F,]
    
    #'  Create a logical vector (T/F) indicating whether spp3 was detected at each
    #'  site where spp2 was detected (based on cameras where spp1 & spp2 were detected)
    spp3.present <- spp2_dat$CameraLocation %in% spp3$CameraLocation
    spp2_dat <- cbind(spp2_dat, spp3.present)
    #'  Split out data into camera locations where both spp1 & spp3 are present vs spp3 absent
    spp2_spp3.present <- spp2_dat[spp2_dat$spp3.present == T,]
    spp2_spp3.absent <- spp2_dat[spp2_dat$spp3.present == F,]
    
    #'  Review sample size per species- smaller sample will determine which coefficient
    #'  of overlap estimator to use (Meredith & Ridout 2017 suggest using delta1 
    #'  for small samples [<50 detection events] and delta4 for larger samples 
    #'  [>50 detection events]).
    ndet_spp1_spp3.present <- nrow(spp1_spp3.present)
    ndet_spp2_spp3.present <- nrow(spp2_spp3.present)
    ndet_spp1_spp3.absent <- nrow(spp1_spp3.absent)
    ndet_spp2_spp3.absent <- nrow(spp2_spp3.absent)
    print(ndet_spp1_spp3.present); print(ndet_spp2_spp3.present)
    print(ndet_spp1_spp3.absent); print(ndet_spp2_spp3.absent)
    
    #'  Visualize general temporal activity with density plots
    densityPlot(spp1$sunTime, rug = T, col = "blue", main = paste0("Density Plot of ", name1, " Daily Activity"))
    densityPlot(spp2$sunTime, rug = T, col = "blue", main = paste0("Density Plot of ", name2, " Daily Activity"))
    densityPlot(spp3$sunTime, rug = T, col = "blue", main = paste0("Density Plot of ", name3, " Daily Activity"))
    
    #'  Visualize temporal overlap
    #'  Overlap when anthropogenic activity is present
    saveOverlap_present <- overlapPlot(spp1_spp3.present$sunTime, spp2_spp3.present$sunTime, rug = T, 
                                       xscale = NA, xcenter = "noon", linet = c(1, 1), linec = c("red", "blue"), 
                                       linew = c(2, 2), main = paste0("Overlap Plots of ", name1, " and ", name2, " \ndiel activity when ", name3, " are present")) 
    #'  Density plot of anthropogenic activity (meaningless if activity is allotments or public land)
    saveDensity <- densityPlot(spp3$sunTime, add = T, xscale = NA, linec = "black", lwd = 2, lty = 2, extend = NULL)
    legend("topleft", c("Predator", "Prey", "Anthro activity"), lty=c(1, 1, 2), col=c("red", "blue", "black"), bg = "white", bty = "n")
    
    #'  Wrangle density data from wide to long format
    DensityA_p <- saveOverlap_present[,1:2] %>%
      mutate(Species = "Predator",
             Anthro_Activity = "Present")
    DensityB_p <- saveOverlap_present[,c(1,3)] %>%
      mutate(Species = "Prey",
             Anthro_Activity = "Present")
    overlap_present <- full_join(DensityA_p, DensityB_p, by = c("x", "Anthro_Activity")) %>%
      full_join(saveDensity, by ="x") %>%
      mutate(Species.z = name3)
    
    #'  Overlap when anthropogenic activity is absent (meaningless when allotment or public land)
    saveOverlap_absent <- overlapPlot(spp1_spp3.absent$sunTime, spp2_spp3.absent$sunTime, rug = T, 
                                      xscale = NA, xcenter = "noon", linet = c(1, 1), linec = c("red", "blue"), 
                                      linew = c(2, 2), main = paste0("Overlap Plots of ", name1, " and ", name2, " \ndiel activity when ", name3, " are absent")) 
    #'  Replot anthropogenic activity density data just for comparison even though it's absent at these sites
    saveDensity <- densityPlot(spp3$sunTime, add = T, xscale = NA, linec = "black", lwd = 2, lty = 2, extend = NULL)
    legend("topleft", c("Predator", "Prey", "Anthro activity"), lty=c(1, 1, 2), col=c("red", "blue", "black"), bg = "white", bty = "n")
    
    #'  Wrangle from wide to long format
    DensityA_a <- saveOverlap_absent[,1:2] %>%
      mutate(Species = "Predator",
             Anthro_Activity = "Absent")
    DensityB_a <- saveOverlap_absent[,c(1,3)] %>%
      mutate(Species = "Prey",
             Anthro_Activity = "Absent")
    overlap_absent <- full_join(DensityA_a, DensityB_a, by = c("x", "Anthro_Activity")) %>%
      full_join(saveDensity, by ="x") %>%
      mutate(Species.z = name3)
    #'  Bind into single long data set of density estimates to make custome overlap plots
    plotdata <- rbind(overlap_present, overlap_absent)
    
    #'  Calculate coefficient of overlap
    dhats_spp1.spp2.spp3 <- overlapEst(A = spp1_spp3.present$sunTime,
                                       B = spp2_spp3.present$sunTime, type = dhat)
    dhats_spp1.spp2.NOspp3 <- overlapEst(A = spp1_spp3.absent$sunTime,
                                         B = spp2_spp3.absent$sunTime, type = dhat)
    
    #'  Bootstrap to estimate standard errors
    #'  FYI: smooth = TRUE is default and allows bootstrap to randomly sample from
    #'  a distribution of times that have a wider range than the original sample
    #'  (see pg. 5 in Overlap package vignette for details)
    spp12.spp3.boot <- bootstrap(spp1_spp3.present$sunTime, spp2_spp3.present$sunTime,
                                 nboot, smooth = TRUE, type = dhat)
    spp12.NOspp3.boot <- bootstrap(spp1_spp3.absent$sunTime, spp2_spp3.absent$sunTime,
                                   nboot, smooth = TRUE, type = dhat)
    #'  Bootstrap mean will be a little different then detla coefficient due to
    #'  bootstrap bias (BSmean - delta) that needs to be accounted for in 95% CIs
    BSmean.present <- mean(spp12.spp3.boot)
    BSmean.absent <- mean(spp12.NOspp3.boot)
    
    #'  Bootstrap 95% Confidence Intervals
    #'  norm0 uses the standard deviation of bootstrap results to calculate CI (delta +/- 1.96*SDboot)
    #'  basic0 takes the 2.5% and 97.5% percentiles and adjusts based on BS bias (percentile - BSbias)
    #'  If sampling distribution is normal, norm0 and basic0 should be similar;
    #'  if sampling distribution is skewed (i.e., if delta is close to 0 or 1) then
    #'  basic0 is the better estimator - using basic0 b/c should be good in either situation
    #'  Using bootCIlogit instead of bootCI so that bias corrections are done on
    #'  the logit scale, then backtransformed. Without this, 95% CIs can fall
    #'  outside (0, 1) interval. See Overlap vignette for more details.
    spp3.present_CI <- bootCIlogit(dhats_spp1.spp2.spp3, spp12.spp3.boot) #[i]
    spp3.absent_CI <- bootCIlogit(dhats_spp1.spp2.NOspp3, spp12.NOspp3.boot) #[i]
    
    #'  Print results
    #'  Effect of spp3 being present
    print("Overlap coefficients when spp3 is present"); print(dhats_spp1.spp2.spp3)
    print("Bootstrap mean"); print(BSmean.present)
    print("Bootstrap 95% CI"); print(spp3.present_CI)
    
    #'  Effect of spp3 being absent
    print("Overlap coefficients when spp3 is present"); print(dhats_spp1.spp2.NOspp3)
    print("Bootstrap mean"); print(BSmean.absent)
    print("Bootstrap 95% CI"); print(spp3.absent_CI)
    
    #'  Save as a giant list
    overlap_list <- list(dhats_spp1.spp2.spp3, dhats_spp1.spp2.NOspp3,
                         spp12.spp3.boot, spp12.NOspp3.boot,
                         spp3.present_CI, spp3.absent_CI, ndet_spp1_spp3.present,
                         ndet_spp2_spp3.present, ndet_spp1_spp3.absent, ndet_spp2_spp3.absent,
                         plotdata)
    names(overlap_list) <- c("dhat_spp3.present", "dhat_spp3.absent", "dhat_spp3.present_boot",
                             "dhat_spp3.absent_boot", "spp3.present_CI", "spp3.absent_CI",
                             "ndet_spp1_spp3.present", "ndet_spp2_spp3.present",
                             "ndet_spp1_spp3.absent", "ndet_spp2_spp3.absent",
                             "overlap.plot.data")
    
    return(overlap_list)
  }
  #####  Predator-Prey Overlap Grazing Season  #####
  #'  -----------------------------------------
  #'  Estimate temporal overlap between predators and prey when cattle are/aren't detected
  #'  NOTES:
  #'   -Focusing on only OK study area since big difference in number of cameras 
  #'    with cattle in NE vs OK, pooling across study areas can bias results
  #'   -Some species pairs commented out owing to too few detections of one or both spp
  #'   -Ran some overlap estimates twice but with different dhat estimators because 
  #'    different sample sizes for detections at sites with/without cattle 
  
  #'  Define number of bootstraps
  nboot <- 10000
  
  #'  Call function for each predator-prey pairing
  coug_md_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Cougar"), 
                                          spp2 = filter(grazing_first_OK, Species == "Mule Deer"), 
                                          spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                          name1 = "Cougar", name2 = "Mule Deer", 
                                          name3 = "Cattle", nboot = nboot, dhat = "Dhat1") 
  # coug_elk_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Cougar"), 
  #                                         spp2 = filter(grazing_first_OK, Species == "Elk"), 
  #                                         spp3 = filter(grazing_first_OK, Species == "Cattle"), 
  #                                         name1 = "Cougar", name2 = "Elk", 
  #                                         name3 = "Cattle", nboot = nboot, dhat = "Dhat1")  
  coug_wtd_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Cougar"), 
                                           spp2 = filter(grazing_first_OK, Species == "White-tailed Deer"), 
                                           spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                           name1 = "Cougar", name2 = "White-tailed Deer", 
                                           name3 = "Cattle", nboot = nboot, dhat = "Dhat1")   
  coug_moose_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Cougar"), 
                                             spp2 = filter(grazing_first_OK, Species == "Moose"), 
                                             spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                             name1 = "Cougar", name2 = "Moose", 
                                             name3 = "Cattle", nboot = nboot, dhat = "Dhat1")  
  # wolf_md_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Wolf"), 
  #                                         spp2 = filter(grazing_first_OK, Species == "Mule Deer"), 
  #                                         spp3 = filter(grazing_first_OK, Species == "Cattle"), 
  #                                         name1 = "Wolf", name2 = "Mule Deer", 
  #                                         name3 = "Cattle", nboot = nboot, dhat = "Dhat1") 
  # wolf_elk_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Wolf"), 
  #                                          spp2 = filter(grazing_first_OK, Species == "Elk"), 
  #                                          spp3 = filter(grazing_first_OK, Species == "Cattle"), 
  #                                          name1 = "Wolf", name2 = "Elk", 
  #                                          name3 = "Cattle", nboot = nboot, dhat = "Dhat1") 
  # wolf_wtd_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Wolf"), 
  #                                          spp2 = filter(grazing_first_OK, Species == "White-tailed Deer"), 
  #                                          spp3 = filter(grazing_first_OK, Species == "Cattle"), 
  #                                          name1 = "Wolf", name2 = "White-tailed Deer", 
  #                                          name3 = "Cattle", nboot = nboot, dhat = "Dhat1") 
  # wolf_moose_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Wolf"), 
  #                                            spp2 = filter(grazing_first_OK, Species == "Moose"), 
  #                                            spp3 = filter(grazing_first_OK, Species == "Cattle"), 
  #                                            name1 = "Wolf", name2 = "Moose", 
  #                                            name3 = "Cattle", nboot = nboot, dhat = "Dhat1")  
  bear_md_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"), 
                                          spp2 = filter(grazing_first_OK, Species == "Mule Deer"), 
                                          spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                          name1 = "Black bear", name2 = "Mule Deer", 
                                          name3 = "Cattle", nboot = nboot, dhat = "Dhat4") 
  # bear_elk_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"), 
  #                                          spp2 = filter(grazing_first_OK, Species == "Elk"), 
  #                                          spp3 = filter(grazing_first_OK, Species == "Cattle"), 
  #                                          name1 = "Black bear", name2 = "Elk", 
  #                                          name3 = "Cattle", nboot = nboot, dhat = "Dhat1")    
  bear_wtd_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"), 
                                           spp2 = filter(grazing_first_OK, Species == "White-tailed Deer"), 
                                           spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                           name1 = "Black bear", name2 = "White-tailed Deer", 
                                           name3 = "Cattle", nboot = nboot, dhat = "Dhat4")   
  bear_moose_graze_over_dhat1 <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"), 
                                                   spp2 = filter(grazing_first_OK, Species == "Moose"), 
                                                   spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                                   name1 = "Black bear", name2 = "Moose", 
                                                   name3 = "Cattle", nboot = nboot, dhat = "Dhat1")   # bear-moose: cattle present all >50, cattle absent moose <50
  bear_moose_graze_over_dhat4 <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"), 
                                                   spp2 = filter(grazing_first_OK, Species == "Moose"), 
                                                   spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                                   name1 = "Black bear", name2 = "Moose", 
                                                   name3 = "Cattle", nboot = nboot, dhat = "Dhat4")   # bear-moose: cattle present all >50, cattle absent moose <50
  bob_md_graze_over_dhat1 <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Bobcat"), 
                                               spp2 = filter(grazing_first_OK, Species == "Mule Deer"), 
                                               spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                               name1 = "Bobcat", name2 = "Mule Deer", 
                                               name3 = "Cattle", nboot = nboot, dhat = "Dhat1")   # bob-md: cattle present bobcat >> 50, cattle absent bobcat <50
  bob_md_graze_over_dhat4 <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Bobcat"), 
                                               spp2 = filter(grazing_first_OK, Species == "Mule Deer"), 
                                               spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                               name1 = "Bobcat", name2 = "Mule Deer", 
                                               name3 = "Cattle", nboot = nboot, dhat = "Dhat4")   # bob-md: cattle present bobcat >> 50, cattle absent bobcat <50
  # bob_wtd_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Bobcat"), 
  #                                          spp2 = filter(grazing_first_OK, Species == "White-tailed Deer"), 
  #                                          spp3 = filter(grazing_first_OK, Species == "Cattle"), 
  #                                          name1 = "Bobcat", name2 = "White-tailed Deer", 
  #                                          name3 = "Cattle", nboot = 10000, dhat = "Dhat1") 
  coy_md_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Coyote"), 
                                         spp2 = filter(grazing_first_OK, Species == "Mule Deer"), 
                                         spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                         name1 = "Coyote", name2 = "Mule Deer", 
                                         name3 = "Cattle", nboot = nboot, dhat = "Dhat4") 
  coy_wtd_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Coyote"), 
                                          spp2 = filter(grazing_first_OK, Species == "White-tailed Deer"), 
                                          spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                          name1 = "Coyote", name2 = "White-tailed Deer", 
                                          name3 = "Cattle", nboot = nboot, dhat = "Dhat4") 
  
  #'  Save all overlap results
  #'  Keep in mind bear-moose and bobcat-md repeat with different dhat estimates
  #'  due to differences in sample size when cattle are/are not detected
  #'  Some pairs excluded owing to small sample sizes:
  #'  coug_elk_graze_over, wolf_md_graze_over, wolf_wtd_graze_over, wolf_moose_graze_over, 
  #'  wolf_elk_graze_over, bear_elk_graze_over, bob_wtd_graze_over
  pred_prey_graze_overlap <- list(coug_md_graze_over, coug_wtd_graze_over, coug_moose_graze_over, 
                                  bear_md_graze_over, bear_wtd_graze_over, 
                                  bear_moose_graze_over_dhat1, bear_moose_graze_over_dhat4,  
                                  bob_md_graze_over_dhat1, bob_md_graze_over_dhat4, 
                                  coy_md_graze_over, coy_wtd_graze_over) 
  # save(pred_prey_graze_overlap, file = "./pred_prey_graze_overlap_OK.RData")
  
  
  
  #####  Predator-Prey Overlap Hunting Season #####
  #'  ----------------------------------------
  #'  Estimate temporal overlap between predators and prey when hunters are/aren't detected
  #'  NOTES:
  #'   -Some species pairs commented out owing to too few detections
  #'   -Ran some overlap estimates twice but with different dhat estimators because 
  #'    different sample sizes for detections at sites with/without hunters
  
  #'  Define number of bootstraps
  nboot <- 10000
  
  #'  Call function for each predator-prey pairing
  coug_md_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Cougar"), 
                                         spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                         spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                         name1 = "Cougar", name2 = "Mule Deer", 
                                         name3 = "Hunters", nboot = nboot, dhat = "Dhat1") 
  coug_elk_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Cougar"), 
                                          spp2 = filter(hunting_first, Species == "Elk"), 
                                          spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                          name1 = "Cougar", name2 = "Elk", 
                                          name3 = "Hunters", nboot = nboot, dhat = "Dhat1")
  coug_wtd_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Cougar"), 
                                          spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                          spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                          name1 = "Cougar", name2 = "White-tailed Deer", 
                                          name3 = "Hunters", nboot = nboot, dhat = "Dhat1")
  coug_moose_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Cougar"), 
                                            spp2 = filter(hunting_first, Species == "Moose"), 
                                            spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                            name1 = "Cougar", name2 = "Moose", 
                                            name3 = "Hunters", nboot = nboot, dhat = "Dhat1")
  wolf_md_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Wolf"),
                                         spp2 = filter(hunting_first, Species == "Mule Deer"),
                                         spp3 = filter(hunting_first, HumanActivity == "Hunter"),
                                         name1 = "Wolf", name2 = "Mule Deer",
                                         name3 = "Hunters", nboot = nboot, dhat = "Dhat1")
  # wolf_elk_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Wolf"),
  #                                   spp2 = filter(hunting_first, Species == "Elk"),
  #                                   spp3 = filter(hunting_first, HumanActivity == "Hunter"),
  #                                   name1 = "Wolf", name2 = "Elk",
  #                                   name3 = "Hunters", nboot = nboot, dhat = "Dhat1")
  wolf_wtd_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Wolf"),
                                          spp2 = filter(hunting_first, Species == "White-tailed Deer"),
                                          spp3 = filter(hunting_first, HumanActivity == "Hunter"),
                                          name1 = "Wolf", name2 = "White-tailed Deer",
                                          name3 = "Hunters", nboot = nboot, dhat = "Dhat1")
  wolf_moose_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Wolf"), 
                                            spp2 = filter(hunting_first, Species == "Moose"), 
                                            spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                            name1 = "Wolf", name2 = "Moose", 
                                            name3 = "Hunters", nboot = nboot, dhat = "Dhat1")
  bear_md_hunt_over_dhat1 <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                               spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                               spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                               name1 = "Black Bear", name2 = "Mule Deer", 
                                               name3 = "Hunters", nboot = nboot, dhat = "Dhat1") # bear-md: Hunter present bear <50, hunter absent all >50
  bear_md_hunt_over_dhat4 <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                               spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                               spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                               name1 = "Black Bear", name2 = "Mule Deer", 
                                               name3 = "Hunters", nboot = nboot, dhat = "Dhat4") # bear-md: Hunter present bear <50, hunter absent all >50
  bear_elk_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                          spp2 = filter(hunting_first, Species == "Elk"), 
                                          spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                          name1 = "Black Bear", name2 = "Elk", 
                                          name3 = "Hunters", nboot = nboot, dhat = "Dhat1")
  bear_wtd_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                          spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                          spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                          name1 = "Black Bear", name2 = "White-tailed Deer", 
                                          name3 = "Hunters", nboot = nboot, dhat = "Dhat1")
  bear_moose_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                            spp2 = filter(hunting_first, Species == "Moose"), 
                                            spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                            name1 = "Black Bear", name2 = "Moose", 
                                            name3 = "Hunters", nboot = nboot, dhat = "Dhat1")
  bob_md_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Bobcat"), 
                                        spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                        spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                        name1 = "Bobcat", name2 = "Mule Deer", 
                                        name3 = "Hunters", nboot = nboot, dhat = "Dhat1")
  bob_wtd_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Bobcat"), 
                                         spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                         spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                         name1 = "Bobcat", name2 = "White-tailed Deer", 
                                         name3 = "Hunters", nboot = nboot, dhat = "Dhat1")
  coy_md_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Coyote"), 
                                        spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                        spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                        name1 = "Coyote", name2 = "Mule Deer", 
                                        name3 = "Hunters", nboot = nboot, dhat = "Dhat4")
  coy_wtd_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Coyote"), 
                                         spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                         spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                         name1 = "Coyote", name2 = "White-tailed Deer", 
                                         name3 = "Hunters", nboot = nboot, dhat = "Dhat4")
  
  #'  Bundle and save!
  pred_prey_hunt_overlap <- list(coug_md_hunt_over, coug_elk_hunt_over, coug_wtd_hunt_over, coug_moose_hunt_over,
                                 wolf_md_hunt_over, wolf_moose_hunt_over, wolf_wtd_hunt_over, 
                                 bear_md_hunt_over_dhat1, bear_md_hunt_over_dhat4, 
                                 bear_elk_hunt_over, bear_wtd_hunt_over, bear_moose_hunt_over,
                                 bob_md_hunt_over, bob_wtd_hunt_over, coy_md_hunt_over, coy_wtd_hunt_over) #wolf_elk_hunt_over, 
  # save(pred_prey_hunt_overlap, file = "./pred_prey_hunt_overlap.RData")
  
  
  #####  Format and visualized predator-prey result  #####
  #'  ----------------------------------------------
  # load("./pred_prey_graze_overlap_OK.RData")
  # load("./pred_prey_hunt_overlap.RData")
  
  ######  Results table for predator-prey overlap  ######
  #'  ----------------------------------------------
  #'  Create results tables from overlap estimates
  #'  Indexing [2,i] gives norm0 CIs, [4,i] gives basic0 CIs
  results_table <- function(overlap_out, spp1, spp2, spp3) {
    dhat_spp3.present <- round(overlap_out[[1]], 2)
    spp3.present_lci <- round(overlap_out[[5]][4,1], 2)
    spp3.present_uci <- round(overlap_out[[5]][4,2], 2)
    dhat_spp3.absent <- round(overlap_out[[2]], 2)
    spp3.absent_lci <- round(overlap_out[[6]][4,1], 2)
    spp3.absent_uci <- round(overlap_out[[6]][4,2], 2)
    pair <- paste0(spp1, "-", spp2)
    spp <- c(pair, pair)
    predator <- c(spp1, spp1)
    prey <- c(spp2, spp2)
    activity <- c("Detected", "Not detected")
    Dhat <- c(dhat_spp3.present, dhat_spp3.absent)
    l95 <- c(spp3.present_lci, spp3.absent_lci)
    u95 <- c(spp3.present_uci, spp3.absent_uci)
    ndet_predator <- c(overlap_out[[7]], overlap_out[[9]])
    ndet_prey <- c(overlap_out[[8]], overlap_out[[10]])
    df <- as.data.frame(cbind(spp, predator, prey, activity, Dhat, l95, u95, 
                              ndet_predator, ndet_prey))
    rownames(df) <- NULL
    names(df)[names(df) == "spp"] <- "Species.pair"
    names(df)[names(df) == "activity"] <- paste0(spp3, ".activity")
    df <- mutate(df, Dhat = as.numeric(Dhat),
                 l95 = as.numeric(l95),
                 u95 = as.numeric(u95))
    return(df)
  }
  #'  Grazing results: Okanogan study area only
  #'  When using two different overlap estimators for same species pairing (e.g.,
  #'  when presence vs absence sample sizes requires different overlap estimators), 
  #'  need to filter to the appropriate estimator given sample size- 
  #'  dhat1 for when cattle are NOT DETECTED, dhat4 for when cattle are DETECTED
  #'  Note: overlap estimate when cattle are detected is always first row,
  #'  overlap estimate when cattle are not detected is always second row 
  coug_md_graze_out <- results_table(pred_prey_graze_overlap[[1]], spp1 = "Cougar", spp2 = "Mule Deer", spp3 = "Grazing")
  coug_wtd_graze_out <- results_table(pred_prey_graze_overlap[[2]], spp1 = "Cougar", spp2 = "White-tailed Deer", spp3 = "Grazing")
  coug_moose_graze_out <- results_table(pred_prey_graze_overlap[[3]], spp1 = "Cougar", spp2 = "Moose", spp3 = "Grazing")
  bear_md_graze_out <- results_table(pred_prey_graze_overlap[[4]], spp1 = "Black bear", spp2 = "Mule Deer", spp3 = "Grazing")
  bear_wtd_graze_out <- results_table(pred_prey_graze_overlap[[5]], spp1 = "Black bear", spp2 = "White-tailed Deer", spp3 = "Grazing")
  #'  bear-moose: cattle present both n>50, cattle absent moose n<50
  bear_moose_graze_out1dhat <- results_table(pred_prey_graze_overlap[[6]], spp1 = "Black bear", spp2 = "Moose", spp3 = "Grazing")
  bear_moose_graze_out1 <- bear_moose_graze_out1dhat[2,]
  bear_moose_graze_out4dhat <- results_table(pred_prey_graze_overlap[[7]], spp1 = "Black bear", spp2 = "Moose", spp3 = "Grazing")
  bear_moose_graze_out4 <- bear_moose_graze_out4dhat[1,]
  #'  bob-md: cattle present both n>50, cattle absent bobcat n<50
  bob_md_graze_out1dhat <- results_table(pred_prey_graze_overlap[[8]], spp1 = "Bobcat", spp2 = "Mule Deer", spp3 = "Grazing")
  bob_md_graze_out1 <- bob_md_graze_out1dhat[2,]
  bob_md_graze_out4dhat <- results_table(pred_prey_graze_overlap[[9]], spp1 = "Bobcat", spp2 = "Mule Deer", spp3 = "Grazing")
  bob_md_graze_out4 <- bob_md_graze_out4dhat[1,]
  coy_md_graze_out <- results_table(pred_prey_graze_overlap[[10]], spp1 = "Coyote", spp2 = "Mule Deer", spp3 = "Grazing")
  coy_wtd_graze_out <- results_table(pred_prey_graze_overlap[[11]], spp1 = "Coyote", spp2 = "White-tailed Deer", spp3 = "Grazing")
  
  cattle_overlap_tbl <- rbind(coug_md_graze_out, coug_wtd_graze_out, coug_moose_graze_out,
                              bear_md_graze_out, bear_wtd_graze_out, bear_moose_graze_out4, 
                              bear_moose_graze_out1, bob_md_graze_out4, bob_md_graze_out1, 
                              coy_md_graze_out, coy_wtd_graze_out)
  # write.csv(cattle_overlap_tbl, file = "./pred-prey_graze_overlap.csv")
  
  #'  Hunting results: NE & OK study areas
  coug_md_hunt_out <- results_table(pred_prey_hunt_overlap[[1]], spp1 = "Cougar", spp2 = "Mule Deer", spp3 = "Hunter")
  coug_elk_hunt_out <- results_table(pred_prey_hunt_overlap[[2]], spp1 = "Cougar", spp2 = "Elk", spp3 = "Hunter")
  coug_wtd_hunt_out <- results_table(pred_prey_hunt_overlap[[3]], spp1 = "Cougar", spp2 = "White-tailed Deer", spp3 = "Hunter")
  coug_moose_hunt_out <- results_table(pred_prey_hunt_overlap[[4]], spp1 = "Cougar", spp2 = "Moose", spp3 = "Hunter")
  wolf_md_hunt_out <- results_table(pred_prey_hunt_overlap[[5]], spp1 = "Wolf", spp2 = "Mule Deer", spp3 = "Hunter")
  wolf_moose_hunt_out <- results_table(pred_prey_hunt_overlap[[6]], spp1 = "Wolf", spp2 = "Moose", spp3 = "Hunter")
  wolf_wtd_hunt_out <- results_table(pred_prey_hunt_overlap[[7]], spp1 = "Wolf", spp2 = "White-tailed Deer", spp3 = "Hunter")
  # bear-md: hunter present bear n<50, hunter absent both n>50
  bear_md_hunt_out1dhat <- results_table(pred_prey_hunt_overlap[[8]], spp1 = "Black bear", spp2 = "Mule Deer", spp3 = "Hunter")
  bear_md_hunt_out1 <- bear_md_hunt_out1dhat[1,]
  bear_md_hunt_out4dhat <- results_table(pred_prey_hunt_overlap[[9]], spp1 = "Black bear", spp2 = "Mule Deer", spp3 = "Hunter")
  bear_md_hunt_out4 <- bear_md_hunt_out4dhat[2,]
  bear_elk_hunt_out <- results_table(pred_prey_hunt_overlap[[10]], spp1 = "Black bear", spp2 = "Elk", spp3 = "Hunter")
  bear_wtd_hunt_out <- results_table(pred_prey_hunt_overlap[[11]], spp1 = "Black bear", spp2 = "White-tailed Deer", spp3 = "Hunter")
  bear_moose_hunt_out <- results_table(pred_prey_hunt_overlap[[12]], spp1 = "Black bear", spp2 = "Moose", spp3 = "Hunter")
  bob_md_hunt_out <- results_table(pred_prey_hunt_overlap[[13]], spp1 = "Bobcat", spp2 = "Mule Deer", spp3 = "Hunter")
  bob_wtd_hunt_out <- results_table(pred_prey_hunt_overlap[[14]], spp1 = "Bobcat", spp2 = "White-tailed Deer", spp3 = "Hunter")
  coy_md_hunt_out <- results_table(pred_prey_hunt_overlap[[15]], spp1 = "Coyote", spp2 = "Mule Deer", spp3 = "Hunter")
  coy_wtd_hunt_out <- results_table(pred_prey_hunt_overlap[[16]], spp1 = "Coyote", spp2 = "White-tailed Deer", spp3 = "Hunter")
  
  hunter_overlap_tbl <- rbind(coug_md_hunt_out, coug_elk_hunt_out, coug_wtd_hunt_out, coug_moose_hunt_out,  
                              wolf_moose_hunt_out, bear_md_hunt_out4, bear_md_hunt_out1, 
                              bear_elk_hunt_out, bear_wtd_hunt_out, bear_moose_hunt_out, 
                              bob_md_hunt_out, bob_wtd_hunt_out, coy_md_hunt_out, coy_wtd_hunt_out)
  # write.csv(hunter_overlap_tbl, file = "./pred-prey_hunt_overlap.csv")

  
  ######  Plot predator-prey coefficients of overlap  ######
  #'  ------------------------------------------------
  #'  Make one single facet_grid plot by grouped by ungulate species
  cattle_overlap_tbl$prey <- factor(cattle_overlap_tbl$prey, levels = c("Moose", "Mule Deer", "White-tailed Deer", "Elk"))
  overlap_grazing_effect <- ggplot(cattle_overlap_tbl, aes(x = predator, y = Dhat, group = Grazing.activity)) +   
    geom_errorbar(aes(ymin = l95, ymax = u95, col = prey), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = 'identity', aes(col = prey, shape = Grazing.activity), size = 3.25, position = position_dodge(width = 0.4)) + 
    scale_colour_bright() +
    ylim(0,1) + theme_bw() +
    theme(text = element_text(size = 22)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(color = "none", shape = guide_legend(title = "Cattle activity")) + 
    ggtitle("Effect of cattle activity on temporal overlap") +
    xlab(NULL) + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~prey, scales = "free", space = "free") 
  overlap_grazing_effect
  
  #'  Make one single facet_grid plot by grouped by ungulate species
  hunter_overlap_tbl$prey <- factor(hunter_overlap_tbl$prey, levels = c("Moose", "Mule Deer", "White-tailed Deer", "Elk"))
  overlap_hunting_effect <- ggplot(hunter_overlap_tbl, aes(x = predator, y = Dhat, group = Hunter.activity)) +   
    geom_errorbar(aes(ymin = l95, ymax = u95, col = prey), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = 'identity', aes(col = prey, shape = Hunter.activity), size = 3.25, position = position_dodge(width = 0.4)) + 
    scale_colour_bright() +
    ylim(0,1) + theme_bw() +
    theme(text = element_text(size = 22)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(color = "none", shape = guide_legend(title = "Hunter activity")) +
    ggtitle("Effect of hunter activity on temporal overlap") +
    xlab("Predator species") + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~prey, scales = "free", space = "free") 
  overlap_hunting_effect
  
  #'  Patchwork figure
  predprey_overlap <- overlap_grazing_effect / overlap_hunting_effect +
    plot_annotation(tag_levels = 'a') +
    plot_annotation(theme = theme(plot.title = element_text(size = 24)))
  predprey_overlap
  # ggsave("./Figures/Overlap_All_Covs_Pred-Prey.tiff", predprey_overlap, 
  #        units = "in", width = 12, height = 15, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  ######  Overlap plots  ######
  #'  --------------------
  #'  Plot temporal overlap when cattle and hunter activity are present vs absent
  overlap_anthro_activity_plots <- function(dat, name1, name2, name3, dhat, y_up, overlapcolor, curvecolor1, curvecolor2) {
    #'  Sample sizes for predators[1] and prey[2] when cattle/hunters are [p]resent or [a]bsent
    n1p <- dat[[7]]; n1a <- dat[[9]]
    n2p <- dat[[8]]; n2a <- dat[[10]]
    spp1p <- paste0(name1, " (n = ", n1p, ")"); spp1a <- paste0(name1, " (n = ", n1a, ")")
    spp2p <- paste0(name2, " (n = ", n2p, ")"); spp2a <- paste0(name2, " (n = ", n2a, ")")
    #'  Temporal overlap between predators and prey when cattle/hunters are present[1] or absent[2]
    dhatp <- dhat[1,5]; dhatpl <- dhat[1,6]; dhatpu<- dhat[1,7]
    dhata <- dhat[2,5]; dhatal <- dhat[2,6]; dhatau<- dhat[2,7]
    #'  Density data for overlap plots
    overdensity <- dat[[11]]
    #'  Separate data sets based on whether cattle/hunter activity is present
    pres <- overdensity[overdensity$Anthro_Activity == "Present",]
    abs <- overdensity[overdensity$Anthro_Activity == "Absent",]
    
    overlap_p <- ggplot(pres, aes(x, densityA, colour = Species.x)) +
      geom_line(lwd = 0.75) + 
      geom_line(aes(x, densityB, colour = Species.y), lwd = 0.75) +  
      geom_area(aes(y = pmin(densityA, densityB)),
                alpha = 0.3, color = NA, fill = overlapcolor) +
      geom_line(aes(x, y, colour =  Species.z), linetype = "dashed", lwd = 0.75) +  
      scale_x_continuous(breaks = c(0, 1.57, 3.0, 4.71, 6.0),
                         labels = c('Midnight', 'Dawn', 'Noon', 'Dusk', 'Midnight')) +
      geom_vline(xintercept = pi/2, linetype="dotted") +
      geom_vline(xintercept = (3*pi)/2, linetype="dotted") +
      theme_bw() +
      theme(text = element_text(size = 20)) +
      theme(legend.background = element_rect(fill = "transparent"),
            legend.key = element_rect(colour = NA, fill = NA)) +
      ylim(0, y_up) +
      labs(x = "Time of day", y = "Density", color = paste0("\u0394 = ", dhatp, " (", dhatpl, " - ", dhatpu, ")"), title = paste0(name3, " present")) + 
      scale_color_manual(labels = c(name3, spp1p, spp2p), values = curvecolor1) 
    plot(overlap_p)
    
    overlap_a <- ggplot(abs, aes(x, densityA, colour = Species.x)) +
      geom_line(lwd = 0.75) + 
      geom_line(aes(x, densityB, colour = Species.y), lwd = 0.75) +  
      geom_area(aes(y = pmin(densityA, densityB)),
                alpha = 0.3, color = NA, fill = overlapcolor) +
      scale_x_continuous(breaks = c(0, 1.57, 3.0, 4.71, 6.0),
                         labels = c('Midnight', 'Dawn', 'Noon', 'Dusk', 'Midnight')) +
      geom_vline(xintercept = pi/2, linetype="dotted") +
      geom_vline(xintercept = (3*pi)/2, linetype="dotted") +
      theme_bw() +
      theme(text = element_text(size = 20)) +
      theme(legend.background = element_rect(fill = "transparent"),
            legend.key = element_rect(colour = NA, fill = NA)) +
      ylim(0, y_up) +
      labs(x = "Time of day", y = "Density", color = paste0("\u0394 = ", dhata, " (", dhatal, " - ", dhatau, ")"), title = paste0(name3, " absent")) + 
      scale_color_manual(labels = c(spp1a, spp2a), values = curvecolor2)  
    plot(overlap_a)
    
    plots <- list(overlap_p, overlap_a)
    return(plots)
  }
  ######  Cattle Activity Overlap Plots  ######
  #'  Keep track of list positions when dhat1 and dhat4 are being combined
  #'  Dhat1 for sample sizes <50, Dhat4 for sample sizes >50, fig [[1]] = present, fig [[2]] = absent
  coug_md_overPlot_g <- overlap_anthro_activity_plots(pred_prey_graze_overlap[[1]], name1 = "Cougar", name2 = "Mule deer", name3 = "Cattle", 
                                                      dhat = coug_md_graze_out, y_up = 0.6, overlapcolor = "hotpink", 
                                                      curvecolor1 = c("black", "maroon", "pink1"), curvecolor2 = c("maroon", "pink1"))
  (coug_md_graze_overlap_plot <- coug_md_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.88)) + coug_md_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.86)))
  # ggsave(coug_md_graze_overlap_plot, filename = "./Figures/Overlap_Plot_coug_md_cattle.tiff", width = 10, height = 6, dpi = 600, units = "in", device='tiff')
  coug_wtd_overPlot_g <- overlap_anthro_activity_plots(pred_prey_graze_overlap[[2]], name1 = "Cougar", name2 = "White-tailed deer", name3 = "Cattle", 
                                                       dhat = coug_wtd_graze_out, y_up = 0.7, overlapcolor = "forestgreen", 
                                                       curvecolor1 = c("black", "darkgreen", "chartreuse3"), curvecolor2 = c("darkgreen", "chartreuse3"))
  (coug_wtd_graze_overlap_plot<- coug_wtd_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.88)) + coug_wtd_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.86)))
  # ggsave(coug_wtd_graze_overlap_plot, filename = "./Figures/Overlap_Plot_coug_wtd_cattle.tiff", width = 10, height = 6, dpi = 600, units = "in", device='tiff')
  coug_moose_overPlot_g <- overlap_anthro_activity_plots(pred_prey_graze_overlap[[3]], name1 = "Cougar", name2 = "Moose", name3 = "Cattle", 
                                                         dhat = coug_moose_graze_out, y_up = 0.6, overlapcolor = "skyblue3", 
                                                         curvecolor1 = c("black", "royalblue4", "steelblue1"), curvecolor2 = c("royalblue4", "steelblue1"))
  (coug_moose_graze_overlap_plot<- coug_moose_overPlot_g[[2]] + theme(legend.position = c(0.32, 0.88)) + coug_moose_overPlot_g[[1]] + theme(legend.position = c(0.32, 0.86)))
  # ggsave(coug_moose_graze_overlap_plot, filename = "./Figures/Overlap_Plot_coug_moose_cattle.tiff", width = 10, height = 6, dpi = 600, units = "in", device='tiff')
  bear_md_overPlot_g <- overlap_anthro_activity_plots(pred_prey_graze_overlap[[4]], name1 = "Black bear", name2 = "Mule deer", name3 = "Cattle", 
                                                      dhat = bear_md_graze_out, y_up = 0.6, overlapcolor = "hotpink", 
                                                      curvecolor1 = c("black", "maroon", "pink1"), curvecolor2 = c("maroon", "pink1"))
  (bear_md_graze_overlap_plot <- bear_md_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.88)) + bear_md_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.86)))
  # ggsave(bear_md_graze_overlap_plot, filename = "./Figures/Overlap_Plot_bear_md_cattle.tiff", width = 10, height = 6, dpi = 600, units = "in", device='tiff')
  bear_wtd_overPlot_g <- overlap_anthro_activity_plots(pred_prey_graze_overlap[[5]], name1 = "Black bear", name2 = "White-tailed deer", name3 = "Cattle", 
                                                       dhat = bear_wtd_graze_out, y_up = 0.6, overlapcolor = "forestgreen", 
                                                       curvecolor1 = c("black", "darkgreen", "chartreuse3"), curvecolor2 = c("darkgreen", "chartreuse3"))
  (bear_wtd_graze_overlap_plot <- bear_wtd_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.88)) + bear_wtd_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.86)))
  # ggsave(bear_wtd_graze_overlap_plot, filename = "./Figures/Overlap_Plot_bear_wtd_cattle.tiff", width = 10, height = 6, dpi = 600, units = "in", device='tiff')
  #'  bear-moose: cattle present both n>50, cattle absent moose n<50
  bear_moose_overPlot_g1 <- overlap_anthro_activity_plots(pred_prey_graze_overlap[[6]], name1 = "Black bear", name2 = "Moose", name3 = "Cattle", 
                                                          dhat = bear_moose_graze_out1dhat, y_up = 0.6, overlapcolor = "skyblue3", 
                                                          curvecolor1 = c("black", "royalblue4", "steelblue1"), curvecolor2 = c("royalblue4", "steelblue1"))
  bear_moose_overPlot_g4 <- overlap_anthro_activity_plots(pred_prey_graze_overlap[[7]], name1 = "Black bear", name2 = "Moose", name3 = "Cattle", 
                                                          dhat = bear_moose_graze_out4dhat, y_up = 0.6, overlapcolor = "skyblue3", 
                                                          curvecolor1 = c("black", "royalblue4", "steelblue1"), curvecolor2 = c("royalblue4", "steelblue1"))
  (bear_moose_graze_overlap_plot <- bear_moose_overPlot_g1[[2]] + theme(legend.position = c(0.30, 0.88)) + bear_moose_overPlot_g4[[1]] + theme(legend.position = c(0.30, 0.86)))
  # ggsave(bear_moose_graze_overlap_plot, filename = "./Figures/Overlap_Plot_bear_moose_cattle.tiff", width = 10, height = 6, dpi = 600, units = "in", device='tiff')
  #'  bob-md: cattle present both n>50, cattle absent bobcat n<50
  bob_md_overPlot_g1 <- overlap_anthro_activity_plots(pred_prey_graze_overlap[[8]], name1 = "Bobcat", name2 = "Mule deer", name3 = "Cattle", 
                                                      dhat = bob_md_graze_out1dhat, y_up = 0.6, overlapcolor = "hotpink", 
                                                      curvecolor1 = c("black", "maroon", "pink1"), curvecolor2 = c("maroon", "pink1"))
  bob_md_overPlot_g4 <- overlap_anthro_activity_plots(pred_prey_graze_overlap[[9]], name1 = "Bobcat", name2 = "Mule deer", name3 = "Cattle", 
                                                      dhat = bob_md_graze_out4dhat, y_up = 0.6, overlapcolor = "hotpink", 
                                                      curvecolor1 = c("black", "maroon", "pink1"), curvecolor2 = c("maroon", "pink1"))
  (bob_md_graze_overlap_plot <- bob_md_overPlot_g1[[2]] + theme(legend.position = c(0.30, 0.88)) + bob_md_overPlot_g4[[1]] + theme(legend.position = c(0.30, 0.86)))
  # ggsave(bob_md_graze_overlap_plot, filename = "./Figures/Overlap_Plot_bob_md_cattle.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coy_md_overPlot_g <- overlap_anthro_activity_plots(pred_prey_graze_overlap[[10]], name1 = "Coyote", name2 = "Mule deer", name3 = "Cattle", 
                                                     dhat = coy_md_graze_out, y_up = 0.6, overlapcolor = "hotpink", 
                                                     curvecolor1 = c("black", "maroon", "pink1"), curvecolor2 = c("maroon", "pink1"))
  (coy_md_graze_overlap_plot <- coy_md_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.88)) + coy_md_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.86)))
  # ggsave(coy_md_graze_overlap_plot, filename = "./Figures/Overlap_Plot_coy_md_cattle.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coy_wtd_overPlot_g <- overlap_anthro_activity_plots(pred_prey_graze_overlap[[11]], name1 = "Coyote", name2 = "White-tailed deer", name3 = "Cattle", 
                                                      dhat = coy_wtd_graze_out, y_up = 0.6, overlapcolor = "forestgreen", 
                                                      curvecolor1 = c("black", "darkgreen", "chartreuse3"), curvecolor2 = c("darkgreen", "chartreuse3"))
  (coy_wtd_graze_overlap_plot <- coy_wtd_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.88)) + coy_wtd_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.86)))
  # ggsave(coy_wtd_graze_overlap_plot, filename = "./Figures/Overlap_Plot_coy_wtd_cattle.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  
  ######  Hunter Activity Overlap Plots  ######
  coug_md_overPlot_h <- overlap_anthro_activity_plots(pred_prey_hunt_overlap[[1]], name1 = "Cougar", name2 = "Mule deer", name3 = "Hunters", 
                                                      dhat = coug_md_hunt_out, y_up = 0.6, overlapcolor = "hotpink", 
                                                      curvecolor1 = c("black", "maroon", "pink1"), curvecolor2 = c("maroon", "pink1"))
  (coug_md_hunt_overlap_plot <- coug_md_overPlot_h[[2]] + theme(legend.position = c(0.35, 0.88)) + coug_md_overPlot_h[[1]] + theme(legend.position = c(0.35, 0.86)))
  # ggsave(coug_md_hunt_overlap_plot, filename = "./Figures/Overlap_Plot_coug_md_hunter.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coug_elk_overPlot_h <- overlap_anthro_activity_plots(pred_prey_hunt_overlap[[2]], name1 = "Cougar", name2 = "Elk", name3 = "Hunters", 
                                                       dhat = coug_elk_hunt_out, y_up = 0.6, overlapcolor = "darkgoldenrod1", 
                                                       curvecolor1 = c("black", "darkgoldenrod2", "gold"), curvecolor2 = c("darkgoldenrod2", "gold"))
  (coug_elk_hunt_overlap_plot <- coug_elk_overPlot_h[[2]] + theme(legend.position = c(0.32, 0.88)) + coug_elk_overPlot_h[[1]] + theme(legend.position = c(0.32, 0.86)))
  # ggsave(coug_elk_hunt_overlap_plot, filename = "./Figures/Overlap_Plot_coug_elk_hunter.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coug_wtd_overPlot_h <- overlap_anthro_activity_plots(pred_prey_hunt_overlap[[3]], name1 = "Cougar", name2 = "White-tailed deer", name3 = "Hunters", 
                                                       dhat = coug_wtd_hunt_out, y_up = 0.6, overlapcolor = "forestgreen", 
                                                       curvecolor1 = c("black", "darkgreen", "chartreuse3"), curvecolor2 = c("darkgreen", "chartreuse3"))
  (coug_wtd_hunt_overlap_plot <- coug_wtd_overPlot_h[[2]] + theme(legend.position = c(0.40, 0.88)) + coug_wtd_overPlot_h[[1]] + theme(legend.position = c(0.40, 0.86)))
  # ggsave(coug_wtd_hunt_overlap_plot, filename = "./Figures/Overlap_Plot_coug_wtd_hunter.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coug_moose_overPlot_h <- overlap_anthro_activity_plots(pred_prey_hunt_overlap[[4]], name1 = "Cougar", name2 = "Moose", name3 = "Hunters", 
                                                         dhat = coug_moose_hunt_out, y_up = 0.6, overlapcolor = "skyblue3", 
                                                         curvecolor1 = c("black", "royalblue4", "steelblue1"), curvecolor2 = c("royalblue4", "steelblue1"))
  (coug_moose_hunt_overlap_plot <- coug_moose_overPlot_h[[2]] + theme(legend.position = c(0.30, 0.88)) + coug_moose_overPlot_h[[1]] + theme(legend.position = c(0.30, 0.86)))
  # ggsave(coug_moose_hunt_overlap_plot, filename = "./Figures/Overlap_Plot_coug_moose_hunter.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  wolf_md_overPlot_h <- overlap_anthro_activity_plots(pred_prey_hunt_overlap[[5]], name1 = "Wolf", name2 = "Mule deer", name3 = "Hunters", 
                                                      dhat = wolf_md_hunt_out, y_up = 0.6, overlapcolor = "hotpink", 
                                                      curvecolor1 = c("black", "maroon", "pink1"), curvecolor2 = c("maroon", "pink1"))
  (wolf_md_hunt_overlap_plot <- wolf_md_overPlot_h[[2]] + theme(legend.position = c(0.30, 0.88)) + wolf_md_overPlot_h[[1]] + theme(legend.position = c(0.32, 0.86)))
  # ggsave(wolf_md_hunt_overlap_plot, filename = "./Figures/Overlap_Plot_wolf_md_hunter.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  wolf_wtd_overPlot_h <- overlap_anthro_activity_plots(pred_prey_hunt_overlap[[6]], name1 = "Wolf", name2 = "White-tailed deer", name3 = "Hunters", 
                                                       dhat = wolf_wtd_hunt_out, y_up = 0.7, overlapcolor = "forestgreen", 
                                                       curvecolor1 = c("black", "darkgreen", "chartreuse3"), curvecolor2 = c("darkgreen", "chartreuse3"))
  (wolf_wtd_hunt_overlap_plot <- wolf_wtd_overPlot_h[[2]] + theme(legend.position = c(0.40, 0.88)) + wolf_wtd_overPlot_h[[1]] + theme(legend.position = c(0.40, 0.86)))
  # ggsave(wolf_wtd_hunt_overlap_plot, filename = "./Figures/Overlap_Plot_wolf_wtd_hunter.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  wolf_moose_overPlot_h <- overlap_anthro_activity_plots(pred_prey_hunt_overlap[[7]], name1 = "Wolf", name2 = "Moose", name3 = "Hunters", 
                                                         dhat = wolf_moose_hunt_out, y_up = 0.7, overlapcolor = "skyblue3", 
                                                         curvecolor1 = c("black", "royalblue4", "steelblue1"), curvecolor2 = c("royalblue4", "steelblue1"))
  (wolf_moose_hunt_overlap_plot <- wolf_moose_overPlot_h[[2]] + theme(legend.position = c(0.35, 0.88)) + wolf_moose_overPlot_h[[1]] + theme(legend.position = c(0.35, 0.86)))
  # ggsave(wolf_moose_hunt_overlap_plot, filename = "./Figures/Overlap_Plot_wolf_moose_hunter.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # bear-md: hunter present bear n<50, hunter absent both n>50
  bear_md_overPlot_h1 <- overlap_anthro_activity_plots(pred_prey_hunt_overlap[[8]], name1 = "Black bear", name2 = "Mule deer", name3 = "Hunters", 
                                                       dhat = bear_md_hunt_out1dhat, y_up = 0.6, overlapcolor = "hotpink", 
                                                       curvecolor1 = c("black", "maroon", "pink1"), curvecolor2 = c("maroon", "pink1"))
  bear_md_overPlot_h4 <- overlap_anthro_activity_plots(pred_prey_hunt_overlap[[9]], name1 = "Black bear", name2 = "Mule deer", name3 = "Hunters", 
                                                       dhat = bear_md_hunt_out4dhat, y_up = 0.6, overlapcolor = "hotpink", 
                                                       curvecolor1 = c("black", "maroon", "pink1"), curvecolor2 = c("maroon", "pink1"))
  (bear_md_hunt_overlap_plot <- bear_md_overPlot_h4[[2]] + theme(legend.position = c(0.25, 0.88)) + bear_md_overPlot_h1[[1]] + theme(legend.position = c(0.25, 0.86)))
  # ggsave(bear_md_hunt_overlap_plot, filename = "./Figures/Overlap_Plot_bear_md_hunter.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  bear_elk_overPlot_h <- overlap_anthro_activity_plots(pred_prey_hunt_overlap[[10]], name1 = "Black bear", name2 = "Elk", name3 = "Hunters", 
                                                       dhat = bear_elk_hunt_out, y_up = 0.6, overlapcolor = "darkgoldenrod1", 
                                                       curvecolor1 = c("black", "darkgoldenrod2", "gold"), curvecolor2 = c("darkgoldenrod2", "gold"))
  (bear_elk_hunt_overlap_plot <- bear_elk_overPlot_h[[2]] + theme(legend.position = c(0.32, 0.88)) + bear_elk_overPlot_h[[1]] + theme(legend.position = c(0.32, 0.86)))
  # ggsave(bear_elk_hunt_overlap_plot, filename = "./Figures/Overlap_Plot_bear_elk_hunter.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  bear_wtd_overPlot_h <- overlap_anthro_activity_plots(pred_prey_hunt_overlap[[11]], name1 = "Black bear", name2 = "White-tailed deer", name3 = "Hunters", 
                                                       dhat = bear_wtd_hunt_out, y_up = 0.6, overlapcolor = "forestgreen", 
                                                       curvecolor1 = c("black", "darkgreen", "chartreuse3"), curvecolor2 = c("darkgreen", "chartreuse3"))
  (bear_wtd_hunt_overlap_plot <- bear_wtd_overPlot_h[[2]] + theme(legend.position = c(0.40, 0.88)) + bear_wtd_overPlot_h[[1]] + theme(legend.position = c(0.40, 0.86)))
  # ggsave(bear_wtd_hunt_overlap_plot, filename = "./Figures/Overlap_Plot_bear_wtd_hunter.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  bear_moose_overPlot_h <- overlap_anthro_activity_plots(pred_prey_hunt_overlap[[12]], name1 = "Black bear", name2 = "Moose", name3 = "Hunters", 
                                                         dhat = bear_moose_hunt_out, y_up = 0.6, overlapcolor = "skyblue3", 
                                                         curvecolor1 = c("black", "royalblue4", "steelblue1"), curvecolor2 = c("royalblue4", "steelblue1"))
  (bear_moose_hunt_overlap_plot <- bear_moose_overPlot_h[[2]] + theme(legend.position = c(0.30, 0.88)) + bear_moose_overPlot_h[[1]] + theme(legend.position = c(0.30, 0.86)))
  # ggsave(bear_moose_hunt_overlap_plot, filename = "./Figures/Overlap_Plot_bear_moose_hunter.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  bob_md_overPlot_h <- overlap_anthro_activity_plots(pred_prey_hunt_overlap[[13]], name1 = "Bobcat", name2 = "Mule deer", name3 = "Hunters", 
                                                     dhat = bob_md_hunt_out, y_up = 0.6, overlapcolor = "hotpink", 
                                                     curvecolor1 = c("black", "maroon", "pink1"), curvecolor2 = c("maroon", "pink1"))
  (bob_md_hunt_overlap_plot <- bob_md_overPlot_h[[2]] + theme(legend.position = c(0.35, 0.88)) + bob_md_overPlot_h[[1]] + theme(legend.position = c(0.35, 0.86)))
  # ggsave(bob_md_hunt_overlap_plot, filename = "./Figures/Overlap_Plot_bob_md_hunter.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  bob_wtd_overPlot_h <- overlap_anthro_activity_plots(pred_prey_hunt_overlap[[14]], name1 = "Bobcat", name2 = "White-tailed deer", name3 = "Hunters", 
                                                      dhat = bob_wtd_hunt_out, y_up = 0.6, overlapcolor = "forestgreen", 
                                                      curvecolor1 = c("black", "darkgreen", "chartreuse3"), curvecolor2 = c("darkgreen", "chartreuse3"))
  (bob_wtd_hunt_overlap_plot <- bob_wtd_overPlot_h[[2]] + theme(legend.position = c(0.40, 0.88)) + bob_wtd_overPlot_h[[1]] + theme(legend.position = c(0.40, 0.86)))
  # ggsave(bob_wtd_hunt_overlap_plot, filename = "./Figures/Overlap_Plot_bob_wtd_hunter.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coy_md_overPlot_h <- overlap_anthro_activity_plots(pred_prey_hunt_overlap[[15]], name1 = "Coyote", name2 = "Mule deer", name3 = "Hunters", 
                                                     dhat = coy_md_hunt_out, y_up = 0.6, overlapcolor = "hotpink", 
                                                     curvecolor1 = c("black", "maroon", "pink1"), curvecolor2 = c("maroon", "pink1"))
  (coy_md_hunt_overlap_plot <- coy_md_overPlot_h[[2]] + theme(legend.position = c(0.35, 0.88)) + coy_md_overPlot_h[[1]] + theme(legend.position = c(0.35, 0.86)))
  # ggsave(coy_md_hunt_overlap_plot, filename = "./Figures/Overlap_Plot_coy_md_hunter.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coy_wtd_overPlot_h <- overlap_anthro_activity_plots(pred_prey_hunt_overlap[[16]], name1 = "Coyote", name2 = "White-tailed deer", name3 = "Hunters", 
                                                      dhat = coy_wtd_hunt_out, y_up = 0.6, overlapcolor = "forestgreen", 
                                                      curvecolor1 = c("black", "darkgreen", "chartreuse3"), curvecolor2 = c("darkgreen", "chartreuse3"))
  (coy_wtd_hunt_overlap_plot <- coy_wtd_overPlot_h[[2]] + theme(legend.position = c(0.40, 0.88)) + coy_wtd_overPlot_h[[1]] + theme(legend.position = c(0.40, 0.86)))
  # ggsave(coy_wtd_hunt_overlap_plot, filename = "./Figures/Overlap_Plot_coy_wtd_hunter.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  
  
  ######  Patchwork plots for publication  ######
  #'  Figure for main text
  pred_prey_overlap_example <- coug_md_graze_overlap_plot / bear_md_hunt_overlap_plot +
    plot_annotation(tag_levels = 'a') + 
    plot_annotation(title = "Effects of cattle and hunter activity on predator - mule deer overlap",
                    theme = theme(plot.title = element_text(size = 24)))
  pred_prey_overlap_example
  # ggsave(pred_prey_overlap_example, filename = "./Figures/Overlap_pred_prey_cattle_hunter.tiff", width = 12, height = 12, dpi = 600, units = "in", device='tiff')
  
  #'  Cattle effects (supplemental figures)
  pred_moose_overlap_cattle <- bear_moose_graze_overlap_plot / coug_moose_graze_overlap_plot +
    plot_annotation(tag_levels = 'i') + 
    plot_annotation(title = "Effects of cattle activity on predator - moose overlap",
                    theme = theme(plot.title = element_text(size = 24)))
  pred_moose_overlap_cattle
  # ggsave(pred_moose_overlap_cattle, filename = "./Figures/Overlap_pred_moose_cattle.tiff", width = 12, height = 12, dpi = 600, units = "in", device='tiff')
  
  pred_md_overlap_cattle_a <- bob_md_graze_overlap_plot / coug_md_graze_overlap_plot +
    plot_annotation(tag_levels = 'i') +  
    plot_annotation(title = "Effects of cattle activity on predator - mule deer overlap",
                    theme = theme(plot.title = element_text(size = 24)))
  pred_md_overlap_cattle_a
  # ggsave(pred_md_overlap_cattle_a, filename = "./Figures/Overlap_pred_md_cattle_1.tiff", width = 12, height = 12, dpi = 600, units = "in", device='tiff')
  
  pred_md_overlap_cattle_b <- coy_md_graze_overlap_plot +
    plot_annotation(tag_levels = 'i') + plot_layout(nrow = 2, ncol = 2) +
    plot_annotation(title = "Effects of cattle activity on predator - mule deer overlap",
                    theme = theme(plot.title = element_text(size = 24)))
  pred_md_overlap_cattle_b
  # ggsave(pred_md_overlap_cattle_b, filename = "./Figures/Overlap_pred_md_cattle_2.tiff", width = 12, height = 12, dpi = 600, units = "in", device='tiff')
  
  pred_wtd_overlap_cattle_a <- bear_wtd_graze_overlap_plot / coug_wtd_graze_overlap_plot +
    plot_annotation(tag_levels = 'i') +  
    plot_annotation(title = "Effects of hunter activity on predator - white-tailed deer overlap",
                    theme = theme(plot.title = element_text(size = 24)))
  pred_wtd_overlap_cattle_a
  # ggsave(pred_wtd_overlap_cattle_a, filename = "./Figures/Overlap_pred_wtd_cattle_1.tiff", width = 12, height = 12, dpi = 600, units = "in", device='tiff')
  
  pred_wtd_overlap_cattle_b <- coy_wtd_graze_overlap_plot +
    plot_annotation(tag_levels = 'i') +  plot_layout(nrow = 2, ncol = 2) +
    plot_annotation(title = "Effects of hunter activity on predator - white-tailed deer overlap",
                    theme = theme(plot.title = element_text(size = 24)))
  pred_wtd_overlap_cattle_b
  # ggsave(pred_wtd_overlap_cattle_b, filename = "./Figures/Overlap_pred_wtd_cattle_2.tiff", width = 12, height = 12, dpi = 600, units = "in", device='tiff')
  
  
  #'  Hunter effects (supplemental figures)
  pred_elk_overlap_hunt <- bear_elk_hunt_overlap_plot / coug_elk_hunt_overlap_plot +
    plot_annotation(tag_levels = 'i') + 
    plot_annotation(title = "Effects of hunter activity on predator - elk overlap",
                    theme = theme(plot.title = element_text(size = 24)))
  pred_elk_overlap_hunt
  # ggsave(pred_elk_overlap_hunt, filename = "./Figures/Overlap_pred_elk_hunter.tiff", width = 12, height = 12, dpi = 600, units = "in", device='tiff')
  
  pred_moose_overlap_hunt <- coug_moose_hunt_overlap_plot / bear_moose_hunt_overlap_plot +
    plot_annotation(tag_levels = 'i') + 
    plot_annotation(title = "Effects of hunter activity on predator - moose overlap",
                    theme = theme(plot.title = element_text(size = 24)))
  pred_moose_overlap_hunt
  # ggsave(pred_moose_overlap_hunt, filename = "./Figures/Overlap_pred_moose_hunter.tiff", width = 12, height = 12, dpi = 600, units = "in", device='tiff')
  
  pred_md_overlap_hunt_a <- bob_md_hunt_overlap_plot / coug_md_hunt_overlap_plot +
    plot_annotation(tag_levels = 'i') +  #plot_layout(widths = c(1, 1, 1, 1)) + #plot_layout(nrow = 2, ncol = 2) +
    plot_annotation(title = "Effects of hunter activity on predator - mule deer overlap",
                    theme = theme(plot.title = element_text(size = 24)))
  pred_md_overlap_hunt_a
  # ggsave(pred_md_overlap_hunt_a, filename = "./Figures/Overlap_pred_md_hunter_1.tiff", width = 12, height = 12, dpi = 600, units = "in", device='tiff')
  
  pred_md_overlap_hunt_b <- coy_md_hunt_overlap_plot / wolf_md_hunt_overlap_plot +
    plot_annotation(tag_levels = 'i') +  
    plot_annotation(title = "Effects of hunter activity on predator - mule deer overlap",
                    theme = theme(plot.title = element_text(size = 24)))
  pred_md_overlap_hunt_b
  # ggsave(pred_md_overlap_hunt_b, filename = "./Figures/Overlap_pred_md_hunter_2.tiff", width = 12, height = 12, dpi = 600, units = "in", device='tiff')
  
  pred_wtd_overlap_hunt_a <- bear_wtd_hunt_overlap_plot / bob_wtd_hunt_overlap_plot +
    plot_annotation(tag_levels = 'i') +  
    plot_annotation(title = "Effects of hunter activity on predator - mule deer overlap",
                    theme = theme(plot.title = element_text(size = 24)))
  pred_wtd_overlap_hunt_a
  # ggsave(pred_wtd_overlap_hunt_a, filename = "./Figures/Overlap_pred_wtd_hunter_1.tiff", width = 12, height = 12, dpi = 600, units = "in", device='tiff')
  
  pred_wtd_overlap_hunt_b <- coug_wtd_hunt_overlap_plot / coy_wtd_hunt_overlap_plot +
    plot_annotation(tag_levels = 'i') +  
    plot_annotation(title = "Effects of hunter activity on predator - mule deer overlap",
                    theme = theme(plot.title = element_text(size = 24)))
  pred_wtd_overlap_hunt_b
  # ggsave(pred_wtd_overlap_hunt_b, filename = "./Figures/Overlap_pred_wtd_hunter_2.tiff", width = 12, height = 12, dpi = 600, units = "in", device='tiff')
  
  
  #'  --------------------------------------------
  ####  Single species temporal overlap analysis  ####
  #'  --------------------------------------------
  #'  Function to estimate differences in temporal activity for a species at camera  
  #'  sites where cattle/humans are present (detected) vs absent (not detected).
  spp_overlap <- function(spp1, spp2, name1, name2, nboot, dhat) { 
    #'  Create logical vectors (T/F) indicating which cameras spp1 was detected 
    #'  at with and without spp2
    both.present <- spp1$CameraLocation %in% spp2$CameraLocation
    spp1_dat <- cbind(spp1, both.present) 
    #'  Split out data into camera locations where both are present vs spp2 absent
    spp1_spp2.present <- spp1_dat[spp1_dat$both.present == T,]
    spp1_spp2.absent <- spp1_dat[spp1_dat$both.present == F,]
    
    #'  Review sample size per species- smaller sample will determine which coefficient
    #'  of overlap estimator to use (delta1 for small samples [<50 detection events], 
    #'  delta4 for larger samples [>75 detection events])
    ndet_spp2.present <- nrow(spp1_spp2.present)
    ndet_spp2.absent <- nrow(spp1_spp2.absent)
    print(ndet_spp2.present); print(ndet_spp2.absent)
    
    #' #'  Visualize general temporal activity with density plots
    densityPlot(spp1$sunTime, rug = T, col = "blue", main = paste0("Density Plot of ", name1, " Daily Activity"))
    densityPlot(spp2$sunTime, rug = T, col = "blue", main = paste0("Density Plot of ", name2, " Daily Activity"))
    
    #'  Visualize temporal overlap
    saveOverlap <- overlapPlot(spp1_spp2.present$sunTime, spp1_spp2.absent$sunTime, rug = T, 
                               xscale = NA, xcenter = "noon", linet = c(1, 1), linec = c("red", "blue"), 
                               linew = c(2, 2), main = paste0("Overlap Plots of ", name1, " diel activity \nwhen ", name2, " are Present and Absente")) 
    saveDensity <- densityPlot(spp2$sunTime, add = T, xscale = NA, linec = "black", lwd = 2, lty = 2, extend = NULL)
    legend("topleft", c("Anthro activity present", "Anthro activity  absent", "Anthro activity"), lty=c(1, 1, 2), col=c("red", "blue", "black"), bg = "white", bty = "n")
    
    #'  Wrangle density data from wide to long format
    DensityA <- saveOverlap[,1:2] %>%
      mutate(Anthro_Activity = "Present")
    DensityB <- saveOverlap[,c(1,3)] %>%
      mutate(Anthro_Activity = "Absent")
    plotdata <- full_join(DensityA, DensityB, by = "x") %>%
      full_join(saveDensity, by ="x") %>%
      mutate(Species.z = name2)
    
    #'  Calculate coefficient of overlap
    dhats_spp1.spp2 <- overlapEst(A = spp1_spp2.present$sunTime, 
                                  B = spp1_spp2.absent$sunTime, type = dhat) 
    
    #'  Bootstrap to estimate standard errors
    #'  FYI: smooth = TRUE is default and allows bootstrap to randomly sample from 
    #'  a distribution of times that have a wider range than the original sample
    #'  (see pg. 5 in Overlap package vignette for details) 
    spp1.spp2.boot <- bootstrap(spp1_spp2.present$sunTime, spp1_spp2.absent$sunTime, 
                                nboot, smooth = TRUE, type = dhat)  
    #'  Bootstrap mean will be a little different then detla coefficient due to
    #'  bootstrap bias (BSmean - delta) that needs to be accounted for in 95% CIs
    BSmean <- mean(spp1.spp2.boot)
    
    #'  Bootstrap 95% Confidence Intervals
    #'  norm0 uses the standard deviation of bootstrap results to calculate CI (delta +/- 1.96*SDboot)
    #'  basic0 takes the 2.5% and 97.5% percentiles and adjusts based on BS bias (percentile - BSbias)
    #'  If sampling distribution is normal, norm0 and basic0 should be similar;
    #'  if sampling distribution is skewed (i.e., if delta is close to 0 or 1) then
    #'  basic0 is the better estimator
    #'  Using bootCIlogit instead of bootCI so that bias corrections are done on
    #'  the logit scale, then backtransformed. Without this, 95% CIs can fall
    #'  outside (0, 1) interval. See Overlap vignette for more details.
    CI <- bootCIlogit(dhats_spp1.spp2, spp1.spp2.boot) #dhats_spp1.spp2[i]
    
    #'  Print results
    #'  Effect of spp2 being present
    print("Overlap coefficients when spp2 is present"); print(dhats_spp1.spp2)
    print("Bootstrap mean"); print(BSmean)
    print("Bootstrap 95% CI"); print(CI)
    
    #'  Save as a giant list
    overlap_list <- list(dhats_spp1.spp2, spp1.spp2.boot, CI, ndet_spp2.present, ndet_spp2.absent, plotdata)
    names(overlap_list) <- c("dhats_spp1.spp2", "spp1.spp2.boot", "CI", "ndet_spp2.present", "ndet_spp2.absent", "overlap.plot.data")
    
    return(overlap_list)
  }
  #'  Estimate temporal overlap for a species when cattle are/aren't detected
  #'  Focusing on only OK study area since big difference in number of cameras 
  #'  with cattle in NE vs OK, pooling across study areas could be confounding
  #'  any apparent temporal patterns
  
  #'  Define number of bootstraps
  nboot <- 10000
  
  #####  Single-Species Overlap Grazing Season  #####
  #'  ------------------------------------------
  coug_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Cougar"),
                                 spp2 = filter(grazing_first_OK, Species == "Cattle"), 
                                 name1 = "Cougar", name2 = "Cattle", nboot = nboot, dhat = "Dhat1") #i = 1
  # wolf_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Wolf"), 
  #                                spp2 = filter(grazing_first_OK, Species == "Cattle"), 
  #                                name1 = "Wolf", name2 = "Cattle", nboot = nboot, dhat = "Dhat1") # only 2 observations when cattle are present
  bear_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"), 
                                 spp2 = filter(grazing_first_OK, Species == "Cattle"), 
                                 name1 = "Black bear", name2 = "Cattle", nboot = nboot, dhat = "Dhat4")
  bob_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Bobcat"), 
                                spp2 = filter(grazing_first_OK, Species == "Cattle"), 
                                name1 = "Bobcat", name2 = "Cattle", nboot = nboot, dhat = "Dhat4")
  coy_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Coyote"), 
                                spp2 = filter(grazing_first_OK, Species == "Cattle"), 
                                name1 = "Coyote", name2 = "Cattle", nboot = nboot, dhat = "Dhat4")
  md_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Mule Deer"), 
                               spp2 = filter(grazing_first_OK, Species == "Cattle"), 
                               name1 = "Mule deer", name2 = "Cattle", nboot = nboot, dhat = "Dhat4")
  # elk_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Elk"), 
  #                               spp2 = filter(grazing_first_OK, Species == "Cattle"), 
  #                               name1 = "Elk", name2 = "Cattle", nboot = nboot, dhat = "Dhat1")  # only 4 observations where cattle absent - OK study area!
  wtd_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "White-tailed Deer"), 
                                spp2 = filter(grazing_first_OK, Species == "Cattle"), 
                                name1 = "White-tailed deer", name2 = "Cattle", nboot = nboot, dhat = "Dhat4")
  moose_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Moose"), 
                                  spp2 = filter(grazing_first_OK, Species == "Cattle"), 
                                  name1 = "Moose", name2 = "Cattle", nboot = nboot, dhat = "Dhat4")
  
  #####  Single-Species Overlap Hunting Season  #####
  #'  ------------------------------------------
  #'  Estimate temporal overlap for a species when hunters are/are not present
  coug_hunt_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Cougar"),
                                spp2 = filter(hunting_first, HumanActivity == "Hunter"), 
                                name1 = "Cougar", name2 = "Hunters", nboot = nboot, dhat = "Dhat4") #i = 1
  wolf_hunt_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Wolf"), 
                                spp2 = filter(hunting_first, HumanActivity == "Hunter"), 
                                name1 = "Wolf", name2 = "Hunters", nboot = nboot, dhat = "Dhat1")
  bear_hunt_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                spp2 = filter(hunting_first, HumanActivity == "Hunter"), 
                                name1 = "Black bear", name2 = "Hunters", nboot = nboot, dhat = "Dhat4")
  bob_hunt_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Bobcat"), 
                               spp2 = filter(hunting_first, HumanActivity == "Hunter"), 
                               name1 = "Bobcat", name2 = "Hunters", nboot = nboot, dhat = "Dhat4")
  coy_hunt_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Coyote"), 
                               spp2 = filter(hunting_first, HumanActivity == "Hunter"), 
                               name1 = "Coyote", name2 = "Hunters", nboot = nboot, dhat = "Dhat4")
  md_hunt_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Mule Deer"), 
                              spp2 = filter(hunting_first, HumanActivity == "Hunter"), 
                              name1 = "Mule deer", name2 = "Hunters", nboot = nboot, dhat = "Dhat4")
  elk_hunt_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Elk"), 
                               spp2 = filter(hunting_first, HumanActivity == "Hunter"), 
                               name1 = "Elk", name2 = "Hunters", nboot = nboot, dhat = "Dhat1")
  wtd_hunt_over <- spp_overlap(spp1 = filter(hunting_first, Species == "White-tailed Deer"), 
                               spp2 = filter(hunting_first, HumanActivity == "Hunter"), 
                               name1 = "White-tailed deer", name2 = "Hunters", nboot = nboot, dhat = "Dhat4")
  moose_hunt_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Moose"), 
                                 spp2 = filter(hunting_first,HumanActivity == "Hunter"), 
                                 name1 = "Moose", name2 = "Hunters", nboot = nboot, dhat = "Dhat4")
  
  #'  List results
  graze_overlap <- list(coug_graze_over, bear_graze_over, bob_graze_over, #wolf_graze_over, elk_graze_over, 
                        coy_graze_over, md_graze_over, wtd_graze_over, moose_graze_over)
  # save(graze_overlap, file = "./grazing_effect_overlap_OK.RData")
  
  hunt_overlap <- list(coug_hunt_over, wolf_hunt_over, bear_hunt_over, bob_hunt_over, 
                       coy_hunt_over, md_hunt_over, elk_hunt_over, wtd_hunt_over, moose_hunt_over)
  # save(hunt_overlap, file = "./hunter_effect_overlap.RData")
  
  
  #####  Format and visualized predator-prey result  #####
  #'  -----------------------------------------------
  # load("./grazing_effect_overlap_OK.RData")
  # load("./hunter_effect_overlap.RData")
  
  ######  Results table for single-species overlap  ######
  #'  ----------------------------------------------
  #'  Create results tables from overlap estimates
  results_table <- function(overlap_out, spp1) {
    Dhat <- round(overlap_out[[1]], 2)
    l95 <- round(overlap_out[[3]][2,1], 2)
    u95 <- round(overlap_out[[3]][2,2], 2)
    Species <- spp1
    ndet_present <- overlap_out[[4]]
    ndet_absent <- overlap_out[[5]]
    df <- as.data.frame(cbind(Species, Dhat, l95, u95, ndet_present, ndet_absent))
    rownames(df) <- NULL
    df <- mutate(df, Dhat = as.numeric(Dhat),
                 l95 = as.numeric(l95),
                 u95 = as.numeric(u95))
    return(df)
  }
  #'  Coefficient of overlap for each species when grazing is and is not detected
  coug_graze_out <- results_table(graze_overlap[[1]], spp1 = "Cougar")
  # wolf_graze_out <- results_table(graze_overlap[[2]], spp1 = "Wolf")
  bear_graze_out <- results_table(graze_overlap[[2]], spp1 = "Black bear")
  bob_graze_out <- results_table(graze_overlap[[3]], spp1 = "Bobcat")
  coy_graze_out <- results_table(graze_overlap[[4]], spp1 = "Coyote")
  md_graze_out <- results_table(graze_overlap[[5]], spp1 = "Mule deer")
  # elk_graze_out <- results_table(graze_overlap[[7]], spp1 = "Elk")
  wtd_graze_out <- results_table(graze_overlap[[6]], spp1 = "White-tailed deer")
  moose_graze_out <- results_table(graze_overlap[[7]], spp1 = "Moose")
  
  #'  Coefficient of overlap for each species when hunters are and are not detected
  coug_hunt_out <- results_table(hunt_overlap[[1]], spp1 = "Cougar")
  wolf_hunt_out <- results_table(hunt_overlap[[2]], spp1 = "Wolf")
  bear_hunt_out <- results_table(hunt_overlap[[3]], spp1 = "Black bear")
  bob_hunt_out <- results_table(hunt_overlap[[4]], spp1 = "Bobcat")
  coy_hunt_out <- results_table(hunt_overlap[[5]], spp1 = "Coyote")
  md_hunt_out <- results_table(hunt_overlap[[6]], spp1 = "Mule deer")
  elk_hunt_out <- results_table(hunt_overlap[[7]], spp1 = "Elk")
  wtd_hunt_out <- results_table(hunt_overlap[[8]], spp1 = "White-tailed deer")
  moose_hunt_out <- results_table(hunt_overlap[[9]], spp1 = "Moose")
  
  
  grazing_effects <- rbind(coug_graze_out, bear_graze_out, bob_graze_out, #wolf_graze_out, elk_graze_out, 
                           coy_graze_out, md_graze_out, wtd_graze_out, moose_graze_out) %>%
    mutate(`Anthropogenic \ndisturbance` = "Cattle activity")
  # write.csv(grazing_effects, file = "./graze_effect_overlap_tbl_OK.csv")
  hunter_effects <- rbind(coug_hunt_out, wolf_hunt_out, bear_hunt_out, bob_hunt_out,
                          coy_hunt_out, md_hunt_out, elk_hunt_out, wtd_hunt_out, moose_hunt_out) %>%
    mutate(`Anthropogenic \ndisturbance` = "Hunter activity")
  # write.csv(hunter_effects, file = "./hunter_effect_overlap_tbl.csv")
  
  
  ######  Plot coefficient of overlap estimates for each species  ######
  #'  ------------------------------------------------------------
  #'  Plot coefficient of overlap estimates for each species
  #'  Make one single facet_grid plot 
  spp_overlap <- rbind(grazing_effects, hunter_effects) #allot_effects, public_effects
  
  spp_overlap_facet <- ggplot(spp_overlap, aes(x = Species, y = Dhat)) +   
    geom_errorbar(aes(ymin = l95, ymax = u95, col = `Anthropogenic \ndisturbance`), width = 0.4) +
    geom_point(stat = 'identity', aes(col = `Anthropogenic \ndisturbance`), size = 3) + 
    scale_colour_vibrant() +
    ylim(0,1) + theme_bw() +
    theme(text = element_text(size = 22),
          axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(color = "none") + 
    ggtitle("Species-specific differences in diel activity patterns") +
    xlab("Species") + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~`Anthropogenic \ndisturbance`, scales = "free", space = "free")
  spp_overlap_facet
  # ggsave(spp_overlap_facet, filename = "./Figures/SpeciesSpecific_Overlap_Plot.tiff", width = 10, height = 8, dpi = 600, units = "in", device='tiff')
  
  
  ######  Overlap plots  ######
  #'  -------------------
  #'  Plot temporal overlap when cattle and hunter activity are present vs absent
  overlap_singlespp_plots <- function(dat, name1, name3, dhat, anthro, lp, y_up, overlapcolor, curvecolor) {
    #'  Sample sizes for predators[1] and prey[2] when cattle/hunters are [p]resent or [a]bsent
    n1p <- dat[[4]]; n1a <- dat[[5]]
    #'  Create labels
    spp1p <- paste0(name3, " present (n = ", n1p, ")"); spp1a <- paste0(name3, " absent (n = ", n1a, ")")
    spp3 <- paste0(anthro, " activity")
    #'  Temporal overlap between predators and prey when cattle/hunters are present[1] or absent[2]
    Dhat <- dhat[2]; Dhatl <- dhat[3]; Dhatu<- dhat[4]
    #'  Density data for overlap plots
    overdensity <- dat[[6]]
    #'  Separate data sets based on whether cattle/hunter activity is present
    pres <- overdensity[overdensity$Anthro_Activity == "Present",]
    abs <- overdensity[overdensity$Anthro_Activity == "Absent",]
    
    overlap <- ggplot(overdensity, aes(x, densityA, colour = Anthro_Activity.x)) +
      geom_line(lwd = 1) + 
      geom_line(aes(x, densityB, colour = Anthro_Activity.y), lwd = 1) +  
      geom_area(aes(y = pmin(densityA, densityB)),
                alpha = 0.3, color = NA, fill = overlapcolor) + 
      geom_line(aes(x, y, colour =  Species.z), linetype = "dashed", lwd = 1) +  
      scale_x_continuous(breaks = c(0, 1.57, 3.0, 4.71, 6.0),
                         labels = c('Midnight', 'Dawn', 'Noon', 'Dusk', 'Midnight')) +
      geom_vline(xintercept = pi/2, linetype="dotted") +
      geom_vline(xintercept = (3*pi)/2, linetype="dotted") +
      theme_bw() +
      theme(legend.background = element_rect(fill = "transparent"),
            legend.key = element_rect(colour = NA, fill = NA),
            legend.position = c(lp, 0.85)) +
      theme(text = element_text(size = 20)) +
      ylim(0, y_up) +
      labs(x = "Time of day", y = "Density", color = paste0("\u0394 = ", Dhat, " (", Dhatl, " - ", Dhatu, ")"), title = paste0(name1, " activity")) +
      scale_colour_manual(breaks = c("Absent", "Present", name3), labels = c(spp1a, spp1p, spp3), values = curvecolor)
    plot(overlap)
    
    return(overlap)
  }
  ######  Cattle Activity Overlap Plots  ######
  coug_overPlot_g <- overlap_singlespp_plots(graze_overlap[[1]], name1 = "Cougar", name3 = "Cattle", dhat = coug_graze_out, 
                                             anthro = "Cattle", lp = 0.30, y_up = 0.85, overlapcolor = "darkorange3", curvecolor = c("orangered3", "sienna1", "black"))
  # ggsave(coug_overPlot_g, filename = "./Figures/Overlap_coug_graze.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  bear_overPlot_g <- overlap_singlespp_plots(graze_overlap[[2]], name1 = "Black bear", name3 = "Cattle", dhat = bear_graze_out, 
                                             anthro = "Cattle", lp = 0.30, y_up = 0.85, overlapcolor = "darkorange3", curvecolor = c("orangered3", "sienna1", "black"))
  # ggsave(bear_overPlot_g, filename = "./Figures/Overlap_bear_graze.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  bob_overPlot_g <- overlap_singlespp_plots(graze_overlap[[3]], name1 = "Bobcat", name3 = "Cattle", dhat = bob_graze_out, 
                                            anthro = "Cattle", lp = 0.30, y_up = 0.85, overlapcolor = "darkorange3", curvecolor = c("orangered3", "sienna1", "black"))
  # ggsave(bob_overPlot_g, filename = "./Figures/Overlap_bob_graze.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  coy_overPlot_g <- overlap_singlespp_plots(graze_overlap[[4]], name1 = "Coyote", name3 = "Cattle", dhat = coy_graze_out, 
                                            anthro = "Cattle", lp = 0.30, y_up = 0.85, overlapcolor = "darkorange3", curvecolor = c("orangered3", "sienna1", "black"))
  # ggsave(coy_overPlot_g, filename = "./Figures/Overlap_coy_graze.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_overPlot_g <- overlap_singlespp_plots(graze_overlap[[5]], name1 = "Mule deer", name3 = "Cattle", dhat = md_graze_out, 
                                           anthro = "Cattle", lp = 0.30, y_up = 0.85, overlapcolor = "darkorange3", curvecolor = c("orangered3", "sienna1", "black"))
  # ggsave(md_overPlot_g, filename = "./Figures/Overlap_md_graze.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_overPlot_g <- overlap_singlespp_plots(graze_overlap[[6]], name1 = "White-tailed deer", name3 = "Cattle", dhat = wtd_graze_out, 
                                            anthro = "Cattle", lp = 0.30, y_up = 0.85, overlapcolor = "darkorange3", curvecolor = c("orangered3", "sienna1", "black"))
  # ggsave(wtd_overPlot_g, filename = "./Figures/Overlap_wtd_graze.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  moose_overPlot_g <- overlap_singlespp_plots(graze_overlap[[7]], name1 = "Moose", name3 = "Cattle", dhat = moose_graze_out, 
                                              anthro = "Cattle", lp = 0.30, y_up = 0.85, overlapcolor = "darkorange3", curvecolor = c("orangered3", "sienna1", "black"))
  # ggsave(moose_overPlot_g, filename = "./Figures/Overlap_moose_graze.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  ######  Hunter Activity Overlap Plots  ######
  coug_overPlot_h <- overlap_singlespp_plots(hunt_overlap[[1]], name1 = "Cougar", name3 = "Hunters", dhat = coug_hunt_out, 
                                             anthro = "Hunter", lp = 0.30, y_up = 0.85, overlapcolor = "steelblue", curvecolor = c("navy", "deepskyblue", "black"))
  # ggsave(coug_overPlot_h, filename = "./Figures/Overlap_coug_hunt.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wolf_overPlot_h <- overlap_singlespp_plots(hunt_overlap[[2]], name1 = "Wolf", name3 = "Hunters", dhat = wolf_hunt_out, 
                                             anthro = "Hunter", lp = 0.30, y_up = 0.85, overlapcolor = "steelblue", curvecolor = c("navy", "deepskyblue", "black"))
  # ggsave(wolf_overPlot_h, filename = "./Figures/Overlap_wolf_hunt.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  bear_overPlot_h <- overlap_singlespp_plots(hunt_overlap[[3]], name1 = "Black bear", name3 = "Hunters", dhat = bear_hunt_out, 
                                             anthro = "Hunter", lp = 0.30, y_up = 0.85, overlapcolor = "steelblue", curvecolor = c("navy", "deepskyblue", "black"))
  # ggsave(bear_overPlot_h, filename = "./Figures/Overlap_bear_hunt.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  bob_overPlot_h <- overlap_singlespp_plots(hunt_overlap[[4]], name1 = "Bobcat", name3 = "Hunters", dhat = bob_hunt_out, 
                                            anthro = "Hunter", lp = 0.30, y_up = 0.85, overlapcolor = "steelblue", curvecolor = c("navy", "deepskyblue", "black"))
  # ggsave(bob_overPlot_h, filename = "./Figures/Overlap_bob_hunt.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  coy_overPlot_h <- overlap_singlespp_plots(hunt_overlap[[5]], name1 = "Coyote", name3 = "Hunters", dhat = coy_hunt_out, 
                                            anthro = "Hunter", lp = 0.30, y_up = 0.85, overlapcolor = "steelblue", curvecolor = c("navy", "deepskyblue", "black"))
  # ggsave(coy_overPlot_h, filename = "./Figures/Overlap_coy_hunt.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_overPlot_h <- overlap_singlespp_plots(hunt_overlap[[6]], name1 = "Mule deer", name3 = "Hunters", dhat = md_hunt_out, 
                                           anthro = "Hunter", lp = 0.30, y_up = 0.85, overlapcolor = "steelblue", curvecolor = c("navy", "deepskyblue", "black"))
  # ggsave(md_overPlot_h, filename = "./Figures/Overlap_md_hunt.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  elk_overPlot_h <- overlap_singlespp_plots(hunt_overlap[[7]], name1 = "Elk", name3 = "Hunters", dhat = elk_hunt_out, 
                                            anthro = "Hunter", lp = 0.30, y_up = 0.85, overlapcolor = "steelblue", curvecolor = c("navy", "deepskyblue", "black"))
  # ggsave(elk_overPlot_h, filename = "./Figures/Overlap_elk_hunt.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_overPlot_h <- overlap_singlespp_plots(hunt_overlap[[8]], name1 = "White-tailed deer", name3 = "Hunters", dhat = wtd_hunt_out, 
                                            anthro = "Hunter", lp = 0.30, y_up = 0.85, overlapcolor = "steelblue", curvecolor = c("navy", "deepskyblue", "black"))
  # ggsave(wtd_overPlot_h, filename = "./Figures/Overlap_wtd_hunt.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  moose_overPlot_h <- overlap_singlespp_plots(hunt_overlap[[9]], name1 = "Moose", name3 = "Hunters", dhat = moose_hunt_out, 
                                              anthro = "Hunter", lp = 0.30, y_up = 0.85, overlapcolor = "steelblue", curvecolor = c("navy", "deepskyblue", "black"))
  # ggsave(moose_overPlot_h, filename = "./Figures/Overlap_moose_hunt.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  
  ######  Patchwork plots for publication  ######
  #'  --------------------------------------
  individual_spp_overlap <- (coug_overPlot_g + wtd_overPlot_g) / (wolf_overPlot_h + moose_overPlot_h) +
    plot_annotation(tag_levels = 'a') + 
    plot_annotation(title = "Effects of cattle and hunter activity on species-specific activity curves",
                    theme = theme(plot.title = element_text(size = 24)))
  individual_spp_overlap
  # ggsave(individual_spp_overlap, filename = "./Figures/Overlap_spp_specific_example.tiff", width = 12, height = 12, dpi = 600, units = "in", device='tiff')
  
  #'  Patchwork figures for supplemental materials
  cattle_overlap_spp_specific <- bear_overPlot_g + bob_overPlot_g + coy_overPlot_g + 
    md_overPlot_g + moose_overPlot_g + plot_annotation(tag_levels = 'a') + plot_layout(nrow = 3, ncol = 2) +
    plot_annotation(title = "Effects of cattle activity on species-specific activity curves",
                    theme = theme(plot.title = element_text(size = 24)))
  cattle_overlap_spp_specific
  # ggsave(cattle_overlap_spp_specific, filename = "./Figures/Overlap_spp_specific_cattle.tiff", width = 14, height = 18, dpi = 600, units = "in", device='tiff')
  
  hunter_overlap_spp_specific <- bear_overPlot_h + bob_overPlot_h + coug_overPlot_h + coy_overPlot_h + 
    elk_overPlot_h + md_overPlot_h + wtd_overPlot_h + plot_annotation(tag_levels = 'a') + plot_layout(nrow = 4, ncol = 2) +
    plot_annotation(title = "Effects of hunter activity on species-specific activity curves",
                    theme = theme(plot.title = element_text(size = 24)))
  hunter_overlap_spp_specific
  # ggsave(hunter_overlap_spp_specific, filename = "./Figures/Overlap_spp_specific_hunter.tiff", width = 13, height = 21, dpi = 600, units = "in", device='tiff')
  
  

  
  
